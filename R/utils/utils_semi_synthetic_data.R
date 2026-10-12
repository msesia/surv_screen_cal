library(R6)
library(dplyr)

RealDataGenerator <- R6Class("RealDataGenerator",
  public = list(
    data_raw = NULL,

    initialize = function(data) {
      stopifnot("time" %in% colnames(data), "status" %in% colnames(data))
      self$data_raw <- data
    },

    sample = function(shuffle = TRUE, return_oracle = FALSE) {
      n <- nrow(self$data_raw)
      df <- if (shuffle) self$data_raw[sample(n), ] else self$data_raw
      df.oracle <- df
      df.oracle$event_time <- NA
      df.oracle$censoring_time <- NA
      if (return_oracle) {
        return(list(observed = df, oracle = df.oracle))
      } else {
        return(df)
      }
    }
  )
)

SemiSyntheticDataGenerator <- R6Class("SemiSyntheticDataGenerator",
  public = list(
    surv_model_type = NULL,
    cens_model_type = NULL,
    surv_model = NULL,
    cens_model = NULL,
    data_raw = NULL,
    surv_gen = NULL,
    cens_gen = NULL,

    initialize = function(data, surv_model_type, cens_model_type) {
      stopifnot("time" %in% colnames(data), "status" %in% colnames(data))
      self$surv_model_type <- surv_model_type
      self$cens_model_type <- cens_model_type
      self$data_raw <- data

      # Initialize and fit survival model
      self$surv_model <- init_surv_model(surv_model_type)
      self$surv_model$fit(Surv(time, status) ~ ., data = data)

      # Initialize and fit censoring model
      cens_model_base <- init_censoring_model(cens_model_type)
      self$cens_model <- CensoringModel$new(model = cens_model_base)
      self$cens_model$fit(data = data)

      # Wrap in imputation models for sampling
      self$surv_gen <- ImputationModel$new(self$surv_model)
      self$cens_gen <- ImputationModel$new(self$cens_model)
    },

    sample = function(shuffle = TRUE, return_oracle = FALSE) {
      n <- nrow(self$data_raw)
      T_new <- self$surv_gen$sample(self$data_raw)
      C_new <- self$cens_gen$sample(self$data_raw)

      df <- self$data_raw %>%
        select(-time, -status) %>%
        mutate(
          event_time = T_new,
          censoring_time = C_new,
          time = pmin(event_time, censoring_time),
          status = as.integer(event_time <= censoring_time)
        )

      if (shuffle) {
        df <- df[sample(n), ]
      }

      data_observed <- df %>% select(time, status, everything(), -censoring_time, -event_time)
      data_oracle <- df %>% select(time, status, event_time, censoring_time, everything())

      if (return_oracle) {
        return(list(observed = data_observed, oracle = data_oracle))
      } else {
        return(data_observed)
      }
    }
  )
)

library(survival)
library(purrr)

SemiSyntheticDataGeneratorSimple <- R6Class("SemiSyntheticDataGeneratorSimple",
  public = list(
    data_raw = NULL,
    surv_model = NULL,
    cens_model = NULL,

    initialize = function(data) {
      stopifnot("time" %in% colnames(data), "status" %in% colnames(data))
      self$data_raw <- data

      # Fit survival model (event = 1)
      self$surv_model <- coxph(Surv(time, status) ~ ., data = data, x = TRUE)

      # Fit censoring model (event = 0)
      data_cens <- data
      data_cens$status <- 1 - data$status
      self$cens_model <- coxph(Surv(time, status) ~ ., data = data_cens, x = TRUE)
    },

    sample_survival_time = function(model, newdata) {
      base_haz <- basehaz(model, centered = FALSE)
      linpred <- predict(model, newdata, type = "lp")
      n <- nrow(newdata)

      map_dbl(1:n, function(i) {
        u <- runif(1)
        S_target <- u
        haz_target <- -log(S_target)
        cum_haz_target <- haz_target / exp(linpred[i])
        idx <- which(base_haz$hazard >= cum_haz_target)[1]
        if (!is.na(idx)) {
          return(base_haz$time[idx])
        } else {
          return(max(base_haz$time))
        }
      })
    },

    sample = function(shuffle = TRUE, return_oracle = FALSE) {
      data <- self$data_raw

      T_new <- self$sample_survival_time(self$surv_model, data)
      C_new <- self$sample_survival_time(self$cens_model, data)

      df <- data %>%
        select(-time, -status) %>%
        mutate(
          event_time = T_new,
          censoring_time = C_new,
          time = pmin(event_time, censoring_time),
          status = as.integer(event_time <= censoring_time)
        )

      if (shuffle) {
        df <- df[sample(nrow(df)), ]
      }

      data_observed <- df %>% select(time, status, everything(), -event_time, -censoring_time)
      data_oracle <- df %>% select(time, status, event_time, censoring_time, everything())

      if (return_oracle) {
        return(list(observed = data_observed, oracle = data_oracle))
      } else {
        return(data_observed)
      }
    }
  )
)

# =============================================================================
# Tilted semi-synthetic data generator.
#
# Extends SemiSyntheticDataGenerator: the survival and censoring models are
# fitted on the raw data as before, but event and censoring times are drawn
# from TILTED versions of the fitted conditional curves,
#
#     S_tilt(t | x) = S_hat(t | x)^{exp(f(x))},
#
# i.e. a proportional-hazards shift with covariate-dependent log hazard ratio
# f(x), or more generally a time-varying log hazard ratio f(x, t) applied to
# the fitted cumulative-hazard increments (see sample_tilted_times). Each tilt
# is given as list(f = ~ <expression in the covariates, optionally t and tau>),
# e.g.
#   tilt_surv = list(f = ~ -1 + 4 * as.numeric(X32 > q.z))      # proportional
#   tilt_surv = list(f = ~ gamma * z * (1 - 2 * t / tau))        # crossing hazards
# where free symbols (q.z, gamma, z) are looked up in the formula's
# environment. A NULL tilt, or f = 0, reproduces the untilted generator.
#
# Two generators sharing tilt_surv but with different tilt_cens describe the
# same survival law under two censoring mechanisms, which is how the
# experiments produce a deliberately misspecified censoring model.
# =============================================================================

# Evaluate a one-sided tilt formula on a data frame -> log hazard ratio.
# The formula may use the covariates and, optionally, the symbol `t` (time)
# and `tau` (the end of the time grid), e.g.
#   ~ gamma * z * (1 - 2 * t / tau)          # time-varying (non-proportional)
# Returns an n x L matrix f(x_i, t_l) when `times` is given; a tilt without
# `t` is constant across columns (proportional hazards).
eval_tilt <- function(tilt, data, times = NULL) {
  n <- nrow(data)
  if (is.null(tilt) || is.null(tilt$f)) {
    return(if (is.null(times)) rep(0, n) else matrix(0, n, length(times)))
  }
  f <- tilt$f
  stopifnot(inherits(f, "formula"), length(f) == 2L)
  uses_t <- "t" %in% all.vars(f)
  env_vars <- as.list(data)
  if (!is.null(times)) env_vars$tau <- max(times)

  eval_at <- function(t_val) {
    out <- eval(f[[2L]], envir = c(env_vars, list(t = t_val)), enclos = environment(f))
    out <- as.numeric(out)
    if (length(out) == 1L) out <- rep(out, n)
    stopifnot(length(out) == n, all(is.finite(out)))
    out
  }
  if (is.null(times)) return(eval_at(NA_real_))
  if (!uses_t) return(matrix(eval_at(NA_real_), n, length(times)))
  sapply(times, eval_at)                               # n x L
}

# Draw one time per row from tilted survival curves by inverse transform.
# `curves` is an n x L matrix of S_hat(t_l | x_i) on the increasing grid
# `times`; `log_hr` is a length-n vector (proportional tilt) or an n x L
# matrix (time-varying tilt) of log hazard ratios. The tilted cumulative
# hazard is sum_l exp(log_hr[i, l]) * dLambda_hat[i, l]; for a constant
# log_hr this reduces to S_hat^exp(log_hr). Times are linearly interpolated
# on the grid; if a curve never drops below U the time is set to the last
# grid point.
sample_tilted_times <- function(curves, times, log_hr) {
  n <- nrow(curves); L <- length(times)
  stopifnot(ncol(curves) == L)
  if (!is.matrix(log_hr)) { stopifnot(length(log_hr) == n); log_hr <- matrix(log_hr, n, L) }
  stopifnot(all(dim(log_hr) == c(n, L)))

  S0 <- pmin(pmax(curves, 1e-12), 1)
  S0 <- t(apply(S0, 1, cummin))                        # enforce monotonicity
  H0 <- -log(S0)                                       # cumulative hazard on the grid
  dH <- cbind(H0[, 1], H0[, -1, drop = FALSE] - H0[, -L, drop = FALSE])   # increments
  H  <- t(apply(exp(log_hr) * dH, 1, cumsum))          # tilted cumulative hazard
  S  <- exp(-H)
  U <- runif(n)
  out <- numeric(n)
  for (i in seq_len(n)) {
    s <- c(1, S[i, ]); tt <- c(0, times)
    j <- which(s <= U[i])[1]                           # first grid point where S <= U
    if (is.na(j)) { out[i] <- times[L]; next }
    if (j == 1L) { out[i] <- 0; next }
    ds <- s[j - 1] - s[j]
    frac <- if (ds > 0) (s[j - 1] - U[i]) / ds else 1
    out[i] <- tt[j - 1] + frac * (tt[j] - tt[j - 1])
  }
  out
}

TiltedSemiSyntheticDataGenerator <- R6Class("TiltedSemiSyntheticDataGenerator",
  inherit = SemiSyntheticDataGenerator,
  public = list(
    tilt_surv = NULL,
    tilt_cens = NULL,
    time_grid = NULL,

    initialize = function(data, surv_model_type, cens_model_type,
                          tilt_surv = NULL, tilt_cens = NULL, n_grid = 200) {
      super$initialize(data, surv_model_type, cens_model_type)
      self$tilt_surv <- tilt_surv
      self$tilt_cens <- tilt_cens
      t_max <- max(data$time)
      self$time_grid <- seq(0, t_max, length.out = n_grid)[-1]
    },

    # Fitted conditional curves on the grid: survival S(t|x) and censoring
    # G(t|x) = P(C > t | x). Assumes the censoring wrapper exposes predict()
    # in the same layout as the survival wrappers.
    curves_surv = function(data) {
      self$surv_model$predict(data, self$time_grid)$predictions
    },
    curves_cens = function(data) {
      self$cens_model$predict(data, self$time_grid)$predictions
    },

    sample = function(shuffle = TRUE, return_oracle = FALSE) {
      data <- self$data_raw
      n <- nrow(data)

      T_new <- sample_tilted_times(self$curves_surv(data), self$time_grid,
                                   eval_tilt(self$tilt_surv, data, self$time_grid))
      C_new <- sample_tilted_times(self$curves_cens(data), self$time_grid,
                                   eval_tilt(self$tilt_cens, data, self$time_grid))

      df <- data %>%
        select(-time, -status) %>%
        mutate(
          event_time     = T_new,
          censoring_time = C_new,
          time   = pmin(event_time, censoring_time),
          status = as.integer(event_time <= censoring_time)
        )

      if (shuffle) df <- df[sample(n), ]

      data_observed <- df %>% select(time, status, everything(), -censoring_time, -event_time)
      data_oracle   <- df %>% select(time, status, event_time, censoring_time, everything())

      if (return_oracle) list(observed = data_observed, oracle = data_oracle) else data_observed
    }
  )
)

# Quick sanity check of two generators: marginal event / censoring rates and
# survival at a horizon, from one draw each.
compare_generators <- function(gen_a, gen_b, t0, labels = c("A", "B"), seed = 1) {
  set.seed(seed)
  summ <- function(g, lab) {
    d <- g$sample(shuffle = FALSE, return_oracle = TRUE)$oracle
    tibble(generator = lab,
           event_rate     = mean(d$status),
           surv_t0        = mean(d$event_time > t0),
           cens_before_t0 = mean(d$censoring_time <= t0),
           median_T       = median(d$event_time),
           median_C       = median(d$censoring_time))
  }
  out <- bind_rows(summ(gen_a, labels[1]), summ(gen_b, labels[2]))
  print(out)
  invisible(out)
}
