library(R6)
#' @importFrom stats approx approxfun splinefun model.response model.frame model.matrix predict pnorm qnorm runif
#' @importFrom survival survreg survreg.control coxph basehaz
#' @importFrom randomForestSRC rfsrc
#' @importFrom grf survival_forest
#' @importFrom xgboost xgb.train xgb.DMatrix setinfo
#' @importFrom R6 R6Class
NULL

#' Initialize a Survival Model
#'
#' @param model_type One of "grf", "grf2", "survreg", "rsf", "cox", "xgb".
#' @param use_covariates Optional character vector of covariate names to use.
#' @return An object of class inheriting from \code{SurvivalModelWrapper}.
#' @export
init_surv_model <- function(model_type, use_covariates = NULL) {
  switch(
    model_type,
    "grf"     = GRF_SurvivalForestWrapper$new(use_covariates = use_covariates),
    "grf2"    = GRF_SurvivalForestWrapper$new(
                  use_covariates = use_covariates,
                  params = list(min.node.size = 5, num.trees = 20, honesty = FALSE)
                ),
    "survreg" = SurvregModelWrapper$new(use_covariates = use_covariates, dist = "lognormal"),
    "rsf"     = randomForestSRC_SurvivalWrapper$new(use_covariates = use_covariates),
    "cox"     = CoxphModelWrapper$new(use_covariates = use_covariates),
    "xgb"     = XgbAftModelWrapper$new(use_covariates = use_covariates),
    stop("Unknown model type: ", model_type)
  )
}

ensure_matrix <- function(new_data) {
  if (is.vector(new_data)) return(matrix(new_data, nrow = 1))
  if (is.matrix(new_data)) return(new_data)
  stop("Input must be either a vector or a matrix")
}

safe_clip <- function(x, lower = 0, upper = 1) {
  clipped <- pmax(lower, pmin(x, upper))
  if (!is.null(dim(x))) dim(clipped) <- dim(x)
  clipped
}

# Monotone interpolation of rfsrc survival curves at arbitrary time points.
get_survival_prob_at_time <- function(pred, time_points) {
  times <- pred$time.interest
  survival_probs <- pred$survival
  out <- matrix(NA, nrow = nrow(survival_probs), ncol = length(time_points))
  for (i in seq_len(nrow(survival_probs))) {
    f <- splinefun(x = times, y = survival_probs[i, ], method = "monoH.FC")
    out[i, ] <- pmin(pmax(f(time_points), 0), 1)
  }
  out
}

# ---- Abstract wrapper --------------------------------------------------------

SurvivalModelWrapper <- R6Class("SurvivalModelWrapper",
  public = list(
    model = NULL,
    formula = NULL,
    time.points = NULL,
    use_covariates = NULL,
    bootstrap_models = NULL,
    data_train = NULL,
    params = NULL,

    initialize = function(use_covariates = NULL, params = NULL) {
      self$use_covariates <- use_covariates
      self$params <- params
    },

    fit = function(formula, data) {
      stop("This method should be implemented in the subclass.")
    },

    # Restricted mean survival time up to tau, from the predicted survival
    # curve on a common time grid. A horizon-free, model-agnostic risk score
    # (larger = lower risk).
    predict_rmst = function(new_data, tau, n_grid = 100) {
      grid <- seq(0, tau, length.out = n_grid)
      S <- self$predict(new_data, grid)$predictions
      S <- pmin(pmax(S, 0), 1)
      dt <- diff(grid)
      as.numeric((S[, -1, drop = FALSE] + S[, -n_grid, drop = FALSE]) %*% dt / 2)
    },

    # Model-agnostic permutation importance: decrease in the mean time-dependent
    # AUC (IPCW, cumulative/dynamic, via timeROC) over the horizons `times` when
    # each covariate is permuted in held-out data. The marker is the RMST up to
    # tau, so the ranking does not depend on a single horizon. Factors are
    # permuted as whole columns. Returns a named numeric vector over covariates.
    #
    # With no arguments, a fraction `valid_frac` of the stored training data is
    # held out: a copy of this model is refitted on the remainder and evaluated
    # on the held-out part. Alternatively pass `data_valid` explicitly, in which
    # case the fitted model itself is used.
    importance = function(data_valid = NULL, times = c(3, 6, 9), tau = NULL,
                          n_perm = 3, valid_frac = 0.2, seed = 1) {
      if (is.null(tau)) tau <- max(times) * 2

      if (is.null(data_valid)) {
        if (is.null(self$data_train)) stop("No training data stored; fit the model first.")
        set.seed(seed)
        n <- nrow(self$data_train)
        idx_valid <- sample.int(n, size = floor(valid_frac * n))
        model_eval <- self$clone(deep = TRUE)
        model_eval$fit(self$formula, self$data_train[-idx_valid, , drop = FALSE])
        data_valid <- self$data_train[idx_valid, , drop = FALSE]
      } else {
        model_eval <- self
      }

      data_valid <- model_eval$select_columns(data_valid)
      vars <- setdiff(colnames(data_valid), c("time", "status"))

      mean_auc <- function(d) {
        marker <- -model_eval$predict_rmst(d, tau)     # higher = higher risk
        roc <- timeROC::timeROC(T = d$time, delta = d$status, marker = marker,
                                cause = 1, times = times, iid = FALSE)
        mean(as.numeric(roc$AUC[-1]))                  # first entry is t = 0
      }
      base <- mean_auc(data_valid)

      set.seed(seed)
      sapply(setNames(vars, vars), function(v) {
        mean(replicate(n_perm, {
          d <- data_valid
          d[[v]] <- sample(d[[v]])
          base - mean_auc(d)
        }))
      })
    },

    # Default: survival curve obtained by inverting predict_quantiles().
    predict = function(new_data, time.points = NULL) {
      if (is.null(time.points)) time.points <- self$time.points
      if (is.null(time.points)) stop("time.points not provided and no default set.")

      probs <- seq(0.01, 0.99, by = 0.01)
      survival_times <- self$predict_quantiles(new_data, probs = probs)

      survival_probs <- matrix(NA, nrow = nrow(survival_times), ncol = length(time.points))
      for (i in seq_len(nrow(survival_times))) {
        if (!all(is.na(survival_times[i, ]))) {
          f <- approxfun(rev(as.numeric(survival_times[i, ])), rev(probs),
                         rule = 2, method = "linear", ties = "ordered")
          survival_probs[i, ] <- f(time.points)
        } else {
          survival_probs[i, ] <- 1
        }
      }
      list(predictions = survival_probs, time.points = time.points)
    },

    predict_confidence_bound = function(new_data, time.points = NULL, level = 0.95,
                                        side = c("upper", "lower"), B = 100, seed = NULL) {
      side <- match.arg(side)
      new_data <- self$select_columns(new_data)
      if (is.null(time.points)) time.points <- self$time.points

      if (is.null(self$bootstrap_models)) {
        if (is.null(self$data_train)) stop("Training data was not saved; cannot bootstrap.")
        if (!is.null(seed)) set.seed(seed)
        formula <- self$formula
        data <- self$data_train
        self$bootstrap_models <- lapply(seq_len(B), function(b) {
          idx <- sample(seq_len(nrow(data)), replace = TRUE)
          model_b <- self$clone(deep = TRUE)
          model_b$fit(formula, data[idx, , drop = FALSE])
          model_b
        })
      }

      n <- nrow(new_data); m <- length(time.points)
      preds <- array(NA, dim = c(n, m, B))
      for (b in seq_len(B)) {
        preds[, , b] <- self$bootstrap_models[[b]]$predict(new_data, time.points)$predictions
      }
      point_est <- self$predict(new_data, time.points)$predictions
      se <- apply(preds, c(1, 2), sd)
      z <- qnorm(1 - (1 - level) / 2)
      bound <- safe_clip(if (side == "upper") point_est + z * se else point_est - z * se)
      list(bound = bound, time.points = time.points, side = side)
    },

    predict_quantiles = function(new_data, probs = c(0.25, 0.5, 0.75)) {
      predictions <- self$predict(new_data)
      curves <- cbind(1, predictions$predictions, 0)
      tp <- c(0, predictions$time.points, max(predictions$time.points) + 1)

      find_quantile <- function(surv, times, p) {
        approx(x = rev(surv), y = rev(times), xout = 1 - p, rule = 2, ties = "ordered")$y
      }
      q <- t(sapply(seq_len(nrow(curves)), function(i) {
        sapply(probs, function(p) find_quantile(curves[i, ], tp, p))
      }))
      q <- ensure_matrix(q)
      colnames(q) <- paste0("Q", probs * 100, "%")
      rownames(q) <- paste0("Individual_", seq_len(nrow(q)))
      as.data.frame(q)
    },

    predict_interpolate = function(survival_probs, original_failure_times, time.points) {
      out <- matrix(NA, nrow = nrow(survival_probs), ncol = length(time.points))
      for (i in seq_len(nrow(survival_probs))) {
        f <- approxfun(original_failure_times, survival_probs[i, ], rule = 2, ties = "ordered")
        out[i, ] <- f(time.points)
      }
      out
    },

    parse_formula = function(formula, data) {
      self$formula <- formula
      response <- model.response(model.frame(formula, data))
      if (!inherits(response, "Surv")) {
        stop("The left-hand side of the formula must be a survival object (Surv(time, status)).")
      }
      list(time = response[, 1], status = response[, 2],
           covariates = model.matrix(formula, data)[, -1, drop = FALSE])
    },

    select_columns = function(new_data) {
      if (!is.null(self$use_covariates)) {
        new_data %>% select(time, status, all_of(self$use_covariates))
      } else {
        new_data
      }
    }
  )
)

# ---- GRF ---------------------------------------------------------------------

GRF_SurvivalForestWrapper <- R6Class("GRF_SurvivalForestWrapper",
  inherit = SurvivalModelWrapper,
  public = list(

    fit = function(formula, data) {
      self$data_train <- data
      data <- self$select_columns(data)
      parsed <- self$parse_formula(formula, data)

      params <- if (is.null(self$params)) list() else as.list(self$params)
      if (is.null(params$num.trees))     params$num.trees <- 100
      if (is.null(params$min.node.size)) params$min.node.size <- 15
      if (is.null(params$honesty))       params$honesty <- TRUE

      self$model <- grf::survival_forest(parsed$covariates, Y = parsed$time, D = parsed$status,
                                         num.trees = params$num.trees,
                                         min.node.size = params$min.node.size,
                                         honesty = params$honesty)
      self$time.points <- self$model$failure.times
      if (length(self$time.points) == 0) self$time.points <- unique(sort(parsed$time))
    },

    predict = function(new_data, time.points = NULL) {
      new_data <- self$select_columns(new_data)
      X_new <- model.matrix(self$formula, new_data)[, -1, drop = FALSE]
      survival_probs <- predict(self$model, newdata = X_new)$predictions
      if (length(survival_probs) == 0) {
        survival_probs <- matrix(1, nrow = nrow(new_data), ncol = length(self$time.points))
      }
      if (!is.null(time.points)) {
        return(list(predictions = self$predict_interpolate(survival_probs, self$time.points, time.points),
                    time.points = time.points))
      }
      list(predictions = survival_probs, time.points = self$time.points)
    }
  )
)

# ---- Random survival forest (randomForestSRC) --------------------------------

randomForestSRC_SurvivalWrapper <- R6Class("randomForestSRC_SurvivalWrapper",
  inherit = SurvivalModelWrapper,
  public = list(

    subsampled_model = NULL,

    fit = function(formula, data, ntree = 100, ...) {
      self$data_train <- data
      # rfsrc requires a base data.frame (tibbles can crash the C code)
      data <- as.data.frame(self$select_columns(data))
      parsed <- self$parse_formula(formula, data)

      self$model <- randomForestSRC::rfsrc(formula, data = data, ntree = ntree, ...)
      self$time.points <- self$model$time.interest
      if (length(self$time.points) == 0) self$time.points <- unique(sort(parsed$time))
    },

    predict = function(new_data, time.points = NULL) {
      new_data <- as.data.frame(self$select_columns(new_data))
      if (is.null(time.points)) time.points <- self$time.points

      pred <- predict(self$model, newdata = new_data)
      survival_probs <- get_survival_prob_at_time(pred, time.points)
      if (length(survival_probs) == 0) {
        survival_probs <- matrix(1, nrow = nrow(new_data), ncol = length(time.points))
      }
      list(predictions = survival_probs, time.points = time.points)
    }
    # predict_confidence_bound() falls back to the bootstrap implementation
    # of the parent class.
  )
)

# ---- Parametric AFT (survreg) ------------------------------------------------

SurvregModelWrapper <- R6Class("SurvregModelWrapper",
  inherit = SurvivalModelWrapper,
  public = list(
    dist = NULL,

    initialize = function(use_covariates = NULL, dist = "weibull") {
      self$use_covariates <- use_covariates
      self$dist <- dist
    },

    fit = function(formula, data) {
      self$data_train <- data
      data <- self$select_columns(data)
      parsed <- self$parse_formula(formula, data)
      self$model <- survival::survreg(formula, data = data, dist = self$dist,
                                      control = survival::survreg.control(maxiter = 1000))
      self$time.points <- unique(sort(parsed$time))
    },

    surv_function = function() {
      switch(self$dist,
        "weibull"   = function(t, mu, scale) 1 - pweibull(t, shape = 1 / scale, scale = exp(mu)),
        "lognormal" = function(t, mu, scale) 1 - plnorm(t, meanlog = mu, sdlog = scale),
        stop("Unsupported distribution: ", self$dist))
    },

    predict_quantiles = function(new_data, probs = c(0.25, 0.5, 0.75)) {
      q <- sapply(probs, function(p) predict(self$model, newdata = new_data, type = "quantile", p = p))
      q <- ensure_matrix(q)
      colnames(q) <- paste0("Q", probs * 100, "%")
      rownames(q) <- paste0("Individual_", seq_len(nrow(q)))
      as.data.frame(q)
    },

    predict = function(new_data, time.points = NULL) {
      new_data <- self$select_columns(new_data)
      if (is.null(time.points)) time.points <- self$time.points
      lp <- predict(self$model, newdata = new_data, type = "lp")
      S <- self$surv_function()
      out <- t(sapply(lp, function(mu) S(time.points, mu, self$model$scale)))
      if (length(time.points) == 1) out <- matrix(out, ncol = 1)
      list(predictions = out, time.points = time.points)
    },

    predict_confidence_bound = function(new_data, time.points = NULL, level = 0.95,
                                        side = c("upper", "lower")) {
      side <- match.arg(side)
      new_data <- self$select_columns(new_data)
      if (is.null(time.points)) time.points <- self$time.points
      z <- qnorm(1 - (1 - level) / 2)
      lp <- predict(self$model, newdata = new_data, type = "lp", se.fit = TRUE)
      S <- self$surv_function()
      mu_ci <- if (side == "upper") lp$fit + z * lp$se.fit else lp$fit - z * lp$se.fit
      out <- t(sapply(mu_ci, function(mu) S(time.points, mu, self$model$scale)))
      if (length(time.points) == 1) out <- matrix(out, ncol = 1)
      list(bound = out, time.points = time.points, side = side)
    }
  )
)

# ---- Cox ---------------------------------------------------------------------

CoxphModelWrapper <- R6Class("CoxphModelWrapper",
  inherit = SurvivalModelWrapper,
  public = list(

    fit = function(formula, data) {
      self$data_train <- data
      data <- self$select_columns(data)
      parsed <- self$parse_formula(formula, data)
      self$model <- survival::coxph(formula, data = data, x = TRUE, y = TRUE)
      self$time.points <- unique(sort(parsed$time))
    },

    predict = function(new_data, time.points = NULL) {
      new_data <- self$select_columns(new_data)
      if (is.null(time.points)) time.points <- self$time.points
      sf <- survival::survfit(self$model, newdata = new_data, se.fit = FALSE)
      out <- matrix(NA_real_, nrow = nrow(new_data), ncol = length(time.points))
      for (i in seq_len(nrow(new_data))) {
        out[i, ] <- summary(sf[i], times = time.points, extend = TRUE)$surv
      }
      list(predictions = out, time.points = time.points)
    },

    predict_confidence_bound = function(new_data, time.points = NULL, level = 0.95,
                                        side = c("upper", "lower"),
                                        conf.type = c("plain", "log", "log-log", "logit", "arcsin")) {
      side <- match.arg(side); conf.type <- match.arg(conf.type)
      new_data <- self$select_columns(new_data)
      if (is.null(time.points)) time.points <- self$time.points

      sf <- survival::survfit(self$model, newdata = new_data, conf.int = level,
                              conf.type = conf.type, se.fit = (conf.type == "plain"))
      out <- matrix(NA_real_, nrow = nrow(new_data), ncol = length(time.points))
      z <- qnorm(1 - (1 - level) / 2)
      for (i in seq_len(nrow(new_data))) {
        s <- summary(sf[i], times = time.points, extend = TRUE)
        out[i, ] <- if (conf.type == "plain") {
          if (side == "upper") pmin(1, s$surv + z * s$std.err) else pmax(0, s$surv - z * s$std.err)
        } else {
          if (side == "upper") s$upper else s$lower
        }
      }
      list(bound = out, time.points = time.points, side = side)
    }
  )
)

# ---- XGBoost AFT -------------------------------------------------------------

XgbAftModelWrapper <- R6::R6Class(
  "XgbAftModelWrapper",
  inherit = SurvivalModelWrapper,
  public = list(

    # Fit an AFT model. The number of boosting rounds is chosen by early
    # stopping on a held-out fraction of the training data; the final model
    # is then refitted on all training data with that number of rounds.
    fit = function(formula, data,
                   params = list(
                     objective = "survival:aft",
                     eval_metric = "aft-nloglik",
                     aft_loss_distribution = "logistic",
                     aft_loss_distribution_scale = 1.5,
                     tree_method = "hist",
                     eta = 0.05,
                     max_depth = 6,
                     subsample = 0.8,
                     colsample_bytree = 0.8
                   ),
                   nrounds = 2000,
                   early_stopping_rounds = 50,
                   valid_frac = 0.2,
                   verbose = 0) {

      self$data_train <- data
      data <- self$select_columns(data)
      parsed <- self$parse_formula(formula, data)
      time <- parsed$time; status <- parsed$status
      self$time.points <- unique(sort(time))

      rhs_vars <- attr(stats::terms(formula, data = data), "term.labels")
      if (length(rhs_vars) == 0L) stop("Formula must include predictors for xgboost AFT.")

      X <- mk_x_matrix_xgb_aft(data, rhs_vars)
      b <- make_aft_bounds_xgb_aft(time, status)
      make_dmat <- function(idx) {
        d <- xgboost::xgb.DMatrix(X[idx, , drop = FALSE])
        xgboost::setinfo(d, "label_lower_bound", b$lb[idx])
        xgboost::setinfo(d, "label_upper_bound", b$ub[idx])
        d
      }

      n <- nrow(X)
      if (!is.null(early_stopping_rounds) && valid_frac > 0 && n >= 50) {
        idx_valid <- sample.int(n, size = max(10L, floor(valid_frac * n)))
        idx_train <- setdiff(seq_len(n), idx_valid)
        bst_es <- xgboost::xgb.train(params = params, data = make_dmat(idx_train),
                                     nrounds = nrounds,
                                     watchlist = list(valid = make_dmat(idx_valid)),
                                     early_stopping_rounds = early_stopping_rounds,
                                     verbose = verbose)
        best_nrounds <- max(1L, bst_es$best_iteration)
      } else {
        best_nrounds <- nrounds
      }

      bst <- xgboost::xgb.train(params = params, data = make_dmat(seq_len(n)),
                                nrounds = best_nrounds, verbose = verbose)

      self$model <- list(bst = bst, features = rhs_vars, x_colnames = colnames(X),
                         dist = params$aft_loss_distribution,
                         sigma = params$aft_loss_distribution_scale,
                         nrounds = best_nrounds)
      invisible(self)
    },

    predict = function(new_data, time.points = NULL) {
      if (is.null(self$model)) stop("Model has not been fitted yet. Call $fit() first.")
      new_data <- self$select_columns(new_data)
      if (is.null(time.points)) time.points <- self$time.points

      times <- suppressWarnings(as.numeric(time.points))
      times <- times[is.finite(times)]
      if (!length(times)) stop("No valid time.points supplied.")
      times[times < 0] <- 0

      X <- mk_x_matrix_xgb_aft(new_data, self$model$features, self$model$x_colnames)
      mu <- as.numeric(predict(self$model$bst, xgboost::xgb.DMatrix(X), outputmargin = TRUE))
      Z <- outer(mu, log(times), function(m, lt) (lt - m) / self$model$sigma)
      S <- .Sbar_from_z_xgb_aft(Z, dist = self$model$dist)
      if (is.null(dim(S))) S <- matrix(S, ncol = length(times))
      colnames(S) <- paste0("t=", times)
      list(predictions = S, time.points = times)
    },

    predict_quantile = function(new_data, probs) {
      if (is.null(self$model)) stop("Model has not been fitted yet. Call $fit() first.")
      new_data <- self$select_columns(new_data)
      probs <- suppressWarnings(as.numeric(probs)); probs <- probs[is.finite(probs)]
      if (!length(probs)) stop("No valid probs supplied.")

      X <- mk_x_matrix_xgb_aft(new_data, self$model$features, self$model$x_colnames)
      mu <- as.numeric(predict(self$model$bst, xgboost::xgb.DMatrix(X), outputmargin = TRUE))
      z <- .z_from_S_xgb_aft(probs, dist = self$model$dist)
      Tmat <- exp(outer(mu, z, function(m, zz) m + self$model$sigma * zz))
      colnames(Tmat) <- paste0("p=", probs)
      list(quantiles = Tmat, probs = probs)
    }
  )
)

# ---- XGBoost helpers ---------------------------------------------------------

mk_x_matrix_xgb_aft <- function(df, features, ref_colnames = NULL) {
  X <- stats::model.matrix(~ . - 1, data = df[, features, drop = FALSE])
  if (!is.null(ref_colnames)) {
    miss <- setdiff(ref_colnames, colnames(X))
    if (length(miss)) {
      X <- cbind(X, matrix(0, nrow = nrow(X), ncol = length(miss), dimnames = list(NULL, miss)))
    }
    X <- X[, ref_colnames, drop = FALSE]
  }
  X
}

make_aft_bounds_xgb_aft <- function(time, status) {
  list(lb = time, ub = ifelse(status == 1, time, Inf))
}

# Survival S = 1 - F(z) for the AFT error distribution
.Sbar_from_z_xgb_aft <- function(z, dist) {
  switch(dist,
    logistic = stats::plogis(-z),
    normal   = stats::pnorm(-z),
    extreme  = exp(-exp(z)),
    stop("Unsupported aft_loss_distribution: ", dist))
}

# Inverse: z such that S(z) = p
.z_from_S_xgb_aft <- function(p, dist) {
  eps <- 1e-8
  p <- pmin(pmax(p, eps), 1 - eps)
  switch(dist,
    logistic = log((1 - p) / p),
    normal   = -stats::qnorm(p),
    extreme  = log(-log(p)),
    stop("Unsupported aft_loss_distribution: ", dist))
}
