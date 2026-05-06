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
#' Instantiates a survival model wrapper object based on the selected model type.
#'
#' @param model_type A character string specifying the model type. One of:
#'   \itemize{
#'     \item \code{"grf"} — Generalized Random Forest via the \code{grf} package.
#'     \item \code{"grf2"} — Variant of GRF with different default parameters.
#'     \item \code{"survreg"} — Parametric survival model via \code{survreg}.
#'     \item \code{"rf"} — Random survival forest via \code{randomForestSRC}.
#'     \item \code{"cox"} — Cox proportional hazards model.
#'     \item \code{"xgb_aft"} — XGBoost Accelerated Failure Time survival model.
#'   }
#' @param use_covariates Optional character vector of covariate names to use.
#'
#' @return An object of class inheriting from \code{SurvivalModelWrapper}.
#' @export
init_surv_model <- function(model_type, use_covariates = NULL) {
  surv_model <- switch(
    model_type,
    "grf"     = GRF_SurvivalForestWrapper$new(use_covariates = use_covariates),
    "grf2"    = GRF_SurvivalForestWrapper$new(
                  use_covariates = use_covariates,
                  params = list(min.node.size = 5, num.trees = 20, honesty = FALSE)
                ),
    "survreg" = SurvregModelWrapper$new(use_covariates = use_covariates, dist = "lognormal"),
    "rf"      = randomForestSRC_SurvivalWrapper$new(use_covariates = use_covariates),
    "cox"     = CoxphModelWrapper$new(use_covariates = use_covariates),
    "xgb"     = XgbAftModelWrapper$new(use_covariates = use_covariates),
    stop("Unknown model type: ", model_type)
  )
  return(surv_model)
}

#' Ensure Input is Matrix
#'
#' Utility function to safely convert vectors to matrices for consistent input handling.
#'
#' @param new_data A vector or matrix.
#' @return A matrix with one row (if input was a vector) or the original matrix.
#' @keywords internal
ensure_matrix <- function(new_data) {
  ## Check if the input is a vector
  if (is.vector(new_data)) {
    ## Convert the vector to a matrix with one row
    new_data <- matrix(new_data, nrow = 1)
  ## Check if the input is already a matrix
  } else if (is.matrix(new_data)) {
    ## If it's already a matrix, return it unchanged
    return(new_data)
  ## Raise an error if the input is neither a vector nor a matrix
  } else {
    stop("Input must be either a vector or a matrix")
  }
  ## Return the resulting matrix
  return(new_data)
}

safe_clip <- function(x, lower = 0, upper = 1) {
  clipped <- pmax(lower, pmin(x, upper))
  # Restore dimensions if x is a matrix or array
  if (!is.null(dim(x))) {
    dim(clipped) <- dim(x)
  }
  return(clipped)
}

#' Interpolate Survival Probabilities at Given Time Points
#'
#' Applies monotonic interpolation (PCHIP-style) to estimate survival probabilities
#' at arbitrary time points from a survival object produced by \code{randomForestSRC::rfsrc()}.
#'
#' @param pred A prediction object containing survival curves and time points.
#' @param time_points A numeric vector of time points at which to interpolate survival probabilities.
#'
#' @return A numeric matrix of interpolated survival probabilities.
#' @keywords internal
get_survival_prob_at_time <- function(pred, time_points) {
  ## Extract the survival times and survival probabilities from the prediction object
  times <- pred$time.interest
  survival_probs <- pred$survival

  ## Initialize a matrix to hold the interpolated probabilities
  num_individuals <- nrow(survival_probs)
  num_time_points <- length(time_points)
  interpolated_probs <- matrix(NA, nrow = num_individuals, ncol = num_time_points)

  ## Perform monotone interpolation using PCHIP for each individual and each time point
  for (i in 1:num_individuals) {
    ## Create the monotone interpolation function for the current individual
    interp_function <- splinefun(x = times, y = survival_probs[i, ], method = "monoH.FC")

    ## Interpolate at each time point and ensure the probabilities are within [0, 1]
    interpolated_probs[i, ] <- pmin(pmax(interp_function(time_points), 0), 1)
  }

  return(interpolated_probs)
}

#' @title Abstract Survival Model Wrapper
#'
#' @description
#' This R6 class provides a common interface for fitting and predicting survival models.
#' Subclasses should implement the \code{fit()} method, and optionally override \code{predict()} and \code{predict_quantiles()}.
#'
#' @field model Fitted model object.
#' @field formula Formula used for model fitting.
#' @field time.points Default failure times used for interpolation.
#' @field use_covariates Optional vector of covariate names to include.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{fit(formula, data)}}{Abstract method to fit a survival model.}
#'   \item{\code{predict(new_data, time.points = NULL)}}{Predict survival curves.}
#'   \item{\code{predict_quantiles(new_data, probs)}}{Return survival quantiles for each individual.}
#'   \item{\code{predict_interpolate(survival_probs, original_failure_times, time.points)}}{Linear interpolation utility.}
#'   \item{\code{parse_formula(formula, data)}}{Extracts time, status, and covariates from formula.}
#'   \item{\code{select_columns(data)}}{Restricts columns to specified covariates.}
#' }
#'
#' @export
SurvivalModelWrapper <- R6Class("SurvivalModelWrapper",
  public = list(
    model = NULL,                ## Holds the trained model object.
    formula = NULL,              ## Stores the formula used to fit the model.
    time.points = NULL,          ## A sequence of time points for which survival probabilities are calculated.
    use_covariates = NULL,       ## List of relevant covariates (e.g., c("X1", "X3")) to use when fitting the censoring model
    bootstrap_models = NULL,     # cache list of bootstrapped models
    data_train = NULL,           # store training data to reuse for boostrapping
    params = NULL,               # store list of parameters     

    #' @description Constructor for SurvivalModelWrapper.
    #' @param use_covariates Optional character vector of covariates to be used.
    initialize = function(use_covariates = NULL, params=NULL) {
      self$use_covariates <- use_covariates
      self$params <- params
    },

    #' @description Abstract method for fitting a survival model. Should be implemented in subclass.
    #' @param formula A survival formula (e.g., Surv(time, status) ~ predictors).
    #' @param data A data.frame with columns for time, status, and covariates.
    fit = function(formula, data) {
      stop("This method should be implemented in the subclass.")
    },


    ## ' @description Predicts survival probabilities for individuals in new data.
    #' @param new_data A data.frame of new individuals with covariates.
    #' @param time.points Optional vector of time points. Defaults to model's time.points.
    #' @return A list with `predictions` (matrix of survival probabilities) and `time.points`.
    predict = function(new_data, time.points = NULL) {
      ## If time.points is not provided, use the default values from the model
      if (is.null(time.points)) {
        time.points <- self$time.points
      }

      ## Ensure time.points is provided either by the user or from the model
      if (is.null(time.points)) {
        stop("Error: time.points are not provided, and there is no default value set.")
      }

      ## Step 1: Predict quantiles for each probability in `probs`
      probs <- seq(0.01, 0.99, by = 0.01)
      survival_times <- self$predict_quantiles(new_data, probs = probs)

      ## Step 2: Initialize a matrix to hold survival probabilities
      survival_probs <- matrix(NA, nrow = nrow(survival_times), ncol = length(time.points))

      ## Step 3: Loop through each individual and interpolate the survival probabilities
      for (i in 1:nrow(survival_times)) {
        if(!all(is.na(survival_times[i, ]))) {
          ## Create an interpolation function for each individual
          interp_fun <- approxfun(rev(as.numeric(survival_times[i, ])), rev(probs), rule = 2, method = "linear", ties = "ordered")
          ## Apply the interpolation function to the provided time.points
          survival_probs[i, ] <- interp_fun(time.points)
        } else {
            survival_probs[i, ] <- 1
        }
      }

      ## Step 4: Return the predicted survival curves and failure times
      list(predictions = survival_probs, time.points = time.points)
    },

    #' @description Predict confidence intervals via bootstrap.
    #' @param new_data New data.frame for prediction.
    #' @param time.points Optional vector of time points.
    #' @param level Confidence level (e.g., 0.95).
    #' @param side "upper" or "lower" confidence bound.
    #' @param B Number of bootstrap samples.
    #' @param seed Random seed for reproducibility.
    predict_confidence_bound = function(new_data, time.points = NULL, level = 0.95, side = c("upper", "lower"), B = 100, seed = NULL) {
        side <- match.arg(side)
        new_data <- self$select_columns(new_data)
        if (is.null(time.points)) {
            time.points <- self$time.points
        }

        if (is.null(self$bootstrap_models)) {
            if (is.null(self$data_train)) {
                stop("Training data was not saved; cannot perform bootstrap. Please set self$data_train during fit.")
            }

            message("Bootstrapping ", B, " models...")

            if(!is.null(seed)) set.seed(seed)
            formula <- self$formula
            data <- self$data_train

            self$bootstrap_models <- lapply(1:B, function(b) {
                idx <- sample(1:nrow(data), replace = TRUE)
                model_b <- self$clone(deep = TRUE)  # new instance of this wrapper
                model_b$fit(formula, data[idx, , drop = FALSE])
                model_b
            })
        }

        ## Predict from each bootstrap model
        n <- nrow(new_data)
        m <- length(time.points)
        bootstrap_preds <- array(NA, dim = c(n, m, B))

        for (b in seq_len(B)) {
            pred_b <- self$bootstrap_models[[b]]$predict(new_data, time.points)
            bootstrap_preds[, , b] <- pred_b$predictions
        }

        ## Point estimate from original model
        point_est <- self$predict(new_data, time.points)$predictions
        se <- apply(bootstrap_preds, c(1, 2), sd)
        z <- qnorm(1 - (1 - level) / 2)

        if(side=="upper") {
            bound_raw <- point_est + z * se
        } else {
            bound_raw <- point_est - z * se
        }
        bound <- safe_clip(bound_raw)

        list(bound = bound, time.points = time.points, side = side)
    },

    #' @description Predicts survival quantiles (e.g. median, quartiles).
    #' @param new_data New data with covariates.
    #' @param probs Numeric vector of quantile probabilities. Defaults to c(0.25, 0.5, 0.75).
    #' @return A data.frame of predicted quantiles for each individual.
    predict_quantiles = function(new_data, probs = c(0.25, 0.5, 0.75)) {
        ## Predict survival curves
        predictions <- self$predict(new_data)
        survival_curves <- predictions$predictions
        time_points <- predictions$time.points  ## Time points associated with the survival curves

        ## Add padding to ensure interpolation will work
        survival_curves <- cbind(1,survival_curves,0)
        time_points <- c(0,time_points,max(time_points)+1)

        ## Function to find the survival time corresponding to a given survival percentile using built-in interpolation
        find_quantile <- function(survival_probs, time_points, percentile) {
            target_prob <- 1 - percentile  ## Convert percentile to survival probability threshold
            ## Use linear interpolation
            interpolated_time <- approx(x = rev(survival_probs), y = rev(time_points), xout = target_prob, rule = 2, ties="ordered")$y
            return(interpolated_time)
        }

        ## Initialize a list to store quantiles for each individual
        quantiles_list <- list()

        ## Loop over each individual
        for (i in 1:nrow(survival_curves)) {
            ## For each individual, find the survival times at the specified percentiles
            quantiles <- sapply(probs, function(p) find_quantile(survival_curves[i, ], time_points, p))
            quantiles_list[[i]] <- quantiles
        }

        ## Convert the list of quantiles to a data frame
        quantiles_df <- do.call(rbind, quantiles_list)
        colnames(quantiles_df) <- paste0("Q", probs * 100, "%")
        rownames(quantiles_df) <- paste0("Individual_", 1:nrow(quantiles_df))

        return(as.data.frame(quantiles_df))
    },

    #' @description Linearly interpolates survival probabilities at new time points.
    #' @param survival_probs Matrix of survival probabilities.
    #' @param original_failure_times Original time points corresponding to survival_probs.
    #' @param time.points New time points to interpolate to.
    #' @return Matrix of interpolated probabilities.
    predict_interpolate = function(survival_probs, original_failure_times, time.points) {
      ## Initialize a matrix to store interpolated survival probabilities
      survival_probs_interp <- matrix(NA, nrow = nrow(survival_probs), ncol = length(time.points))

      ## Loop through each individual and interpolate the survival probabilities at the custom times
      for (i in 1:nrow(survival_probs)) {
        ## Create an interpolation function for each individual
        interp_fun <- approxfun(original_failure_times, survival_probs[i, ], rule = 2, ties = "ordered")

        ## Apply the interpolation function to the custom failure times
        survival_probs_interp[i, ] <- interp_fun(time.points)
      }

      ## Return the interpolated survival probabilities
      return(survival_probs_interp)
    },

    #' @description Parses a survival formula and extracts components.
    #' @param formula A survival formula, e.g., Surv(time, status) ~ X1 + X2.
    #' @param data A data.frame containing survival and covariate columns.
    #' @return A list with `time`, `status`, and `covariates`.
    parse_formula = function(formula, data) {
      self$formula <- formula
      response <- model.response(model.frame(formula, data))
      time <- response[, 1]
      status <- response[, 2]
      covariates <- model.matrix(formula, data)[, -1, drop = FALSE]  ## Remove intercept
      list(time = time, status = status, covariates = covariates)
    },

    #' @description Subsets the data to only include selected covariates.
    #' @param new_data A data.frame to filter.
    #' @return A filtered data.frame with only the relevant covariates.
    select_columns = function(new_data) {
      if(!is.null(self$use_covariates)) {
          new_data_sel <- new_data %>% select(time, status, self$use_covariates)
      } else {
          new_data_sel <- new_data
      }
      return(new_data_sel)
  }

  ),

)

#' @title GRF Survival Forest Model
#'
#' @description
#' Implements a wrapper for generalized random forest survival models using the \pkg{grf} package.
#'
#' @export
GRF_SurvivalForestWrapper <- R6Class("GRF_SurvivalForestWrapper",
  inherit = SurvivalModelWrapper,
  public = list(

    #' @description
    #' Fit a generalized random forest survival model using the \pkg{grf} package.
    #'
    #' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
    #' @param data A data.frame containing columns for time, status, and covariates.
    fit = function(formula, data) {
      self$data_train = data
      data <- self$select_columns(data)
      parsed_data <- self$parse_formula(formula, data)

      ## Parameters
      params <- if (is.null(self$params)) list() else as.list(self$params)
      if (is.null(params$num.trees))     params$num.trees <- 100
      if (is.null(params$min.node.size)) params$min.node.size <- 15
      if (is.null(params$honesty))       params$honesty <- TRUE
      ## Fit the survival forest model
      self$model <- grf::survival_forest(parsed_data$covariates,
                                         Y = parsed_data$time,
                                         D = parsed_data$status,
                                         num.trees = params$num.trees,
                                         min.node.size = params$min.node.size,
                                         honesty = params$honesty
                                         )
      self$time.points <- self$model$failure.times  ## Extract the default failure times
      if(length(self$time.points)==0) {
          self$time.points <- unique(sort(parsed_data$time))
      }
    },

    #' @description
    #' Predicts survival probabilities at specified time points for new data.
    #'
    #' @param new_data A data.frame of new observations with the same structure as the training data.
    #' @param time.points Optional numeric vector of time points at which to compute predictions.
    #'
    #' @return A list with:
    #' \describe{
    #'   \item{\code{predictions}}{Matrix of survival probabilities (rows = individuals, cols = time points).}
    #'   \item{\code{time.points}}{Vector of time points used in prediction.}
    #' }
    predict = function(new_data, time.points = NULL) {
      new_data <- self$select_columns(new_data)

      ## Generate the design matrix from the new data using the stored formula
      covariates_new <- model.matrix(self$formula, new_data)[, -1, drop = FALSE]

      ## Predict survival curves
      predictions <- predict(self$model, newdata = covariates_new)

      ## Use default failure times if custom ones are not provided
      original_failure_times <- self$time.points
      survival_probs <- predictions$predictions

      ## If predictions fail, replace survival_probs with ones
      if(length(survival_probs)==0){
          survival_probs = matrix(1, nrow = nrow(new_data), ncol = length(self$time.points))
      }

      ## If custom failure times are provided, use the default interpolation method
      if (!is.null(time.points)) {
        survival_probs_interp <- self$predict_interpolate(survival_probs, original_failure_times, time.points)
        return(list(predictions = survival_probs_interp, time.points = time.points))
      }

      ## If no custom failure times are provided, return the original predictions
      return(list(predictions = survival_probs, time.points = original_failure_times))
    }

  )
)

#' @title RandomForestSRC Survival Wrapper
#'
#' @description
#' Wrapper for the \code{randomForestSRC::rfsrc()} survival model.
#'
#' @export
randomForestSRC_SurvivalWrapper <- R6Class("randomForestSRC_SurvivalWrapper",
  inherit = SurvivalModelWrapper,
  public = list(

    #' @description
    #' Fits a survival forest model using the \pkg{randomForestSRC} package.
    #' This method parses the survival formula and fits a random survival forest to the data.
    #'
    #' @param formula A survival formula, typically \code{Surv(time, status) ~ predictors}.
    #' @param data A data.frame containing time, status, and covariates.
    #' @param ntree Number of trees to grow in the forest (default = 100).
    #' @param ... Additional arguments passed to \code{randomForestSRC::rfsrc()}.
    fit = function(formula, data, ntree = 100, ...) {
      self$data_train = data
      data <- self$select_columns(data)

      ## Extract the model frame based on the formula and data
      mf <- model.frame(formula, data)

      ## Ensure the response is a survival object (Surv)
      response <- model.response(mf)
      if (!inherits(response, "Surv")) {
        stop("The left-hand side of the formula must be a survival object (Surv(time, status)).")
      }
      ## Fit the survival forest model using randomForestSRC
      self$model <- randomForestSRC::rfsrc(formula, data = data, ntree = ntree, save.memory=TRUE, ...)

      ## Store the failure times from the model
      self$time.points <- self$model$time.interest

      if(length(self$time.points)==0) {
          self$time.points <- unique(sort(parsed_data$time))
      }
    },

    #' @description
    #' Predicts survival probabilities at specified time points for new data.
    #'
    #' @param new_data A data.frame of new observations with the same structure as the training data.
    #' @param time.points Optional numeric vector of time points at which to compute predictions.
    #'
    #' @return A list with:
    #' \describe{
    #'   \item{\code{predictions}}{Matrix of survival probabilities (rows = individuals, cols = time points).}
    #'   \item{\code{time.points}}{Vector of time points used in prediction.}
    #' }
    predict = function(new_data, time.points = NULL) {
      new_data <- self$select_columns(new_data)
      ## Ensure that new_data is correctly formatted
      if (!is.data.frame(new_data)) {
        stop("new_data must be a data frame.")
      }

      ## Use the model's failure times if custom failure times are not provided
      if (is.null(time.points)) {
        time.points <- self$time.points
      }

      ## Predict survival curves using the trained random forest model
      predictions <- predict(self$model, newdata = new_data)
      ## Extract survival probabilities for each individual at each failure time
      survival_probs <- get_survival_prob_at_time(predictions, time.points)

      ## If predictions fail, replace survival_probs with ones
      if(length(survival_probs)==0){
          survival_probs = matrix(1, nrow = nrow(new_data), ncol = length(self$time.points))
      }

      ## Return predictions and failure times
      list(predictions = survival_probs, time.points = time.points)
    },
    #' @description Predict confidence bounds for survival probability.
    #' @param new_data New covariate data.
    #' @param time.points Time grid to evaluate survival curves.
    #' @param level Confidence level (e.g., 0.95).
    #' @param side Side of confidence bound (either upper or lower)
    #' @return Matrix of upper bounds for survival probability.
    predict_confidence_bound = function(new_data, time.points = NULL, level = 0.95, side = c("upper", "lower")) {
        side <- match.arg(side)
        new_data <- self$select_columns(new_data)

        if (is.null(time.points)) {
            time.points <- self$time.points
        }

        ## Perform subsampling if not already done
        if (is.null(self$subsampled_model)) {
            self$subsampled_model <- randomForestSRC::subsample(self$model, B = 10, importance = "permute")
        }

        ## Predict with subsampled model
        pred <- predict(self$subsampled_model, newdata = new_data)

        ## Extract survival predictions and variance
        survival <- get_survival_prob_at_time(pred, time.points)
        survival_var <- get_survival_prob_at_time(pred, time.points, type = "var")

        ## Confidence interval
        z <- qnorm(1 - (1 - level) / 2)
        se <- sqrt(survival_var)

        bound <- switch(
            side,
            "upper" = pmin(1, survival + z * se),
            "lower" = pmax(0, survival - z * se)
        )

        return(list(bound = bound, time.points = time.points, side = side))
    }


  )
)

#' @title Parametric Survival Model Wrapper
#'
#' @description
#' Wrapper for parametric survival models via \code{survival::survreg()}.
#'
#' @field dist A character string indicating the distribution used by the model (e.g., "weibull", "lognormal").
#'
#' @export
SurvregModelWrapper <- R6Class("SurvregModelWrapper",
  inherit = SurvivalModelWrapper,
  public = list(
    dist = NULL,  ## Distribution parameter

    #' @description Constructor for SurvregModelWrapper.
    #' @param use_covariates Optional character vector of covariate names to use.
    #' @param dist Distribution to use in the survreg model (e.g., \"weibull\", \"lognormal\").
    initialize = function(use_covariates = NULL, dist = "weibull") {
      self$use_covariates <- use_covariates
      self$dist <- dist
    },

    #' @description Fit a parametric survival model using survreg.
    #' @param formula A formula of the form \code{Surv(time, status) ~ predictors}.
    #' @param data A data.frame with survival outcome and covariates.
    fit = function(formula, data) {
      self$data_train = data
      data <- self$select_columns(data)
      parsed_data <- self$parse_formula(formula, data)
      self$model <- survival::survreg(formula, data = data, dist = self$dist, control = survival::survreg.control(maxiter = 1000))
      self$time.points <- unique(sort(parsed_data$time))
    },

    #' @description Predict survival time quantiles.
    #' @param new_data A data.frame of covariates.
    #' @param probs A numeric vector of quantile probabilities.
    #' @return A data.frame of quantile predictions.
    predict_quantiles = function(new_data, probs = c(0.25, 0.5, 0.75)) {

      ## Predict quantiles for each probability in probs
      quantiles_matrix <- sapply(probs, function(p) {
        predict(self$model, newdata = new_data, type = "quantile", p = p)
      })

      ## Ensure quantiles_matrix is a matrix
      quantiles_matrix <- ensure_matrix(quantiles_matrix)

      ## Set column names to represent the quantiles
      colnames(quantiles_matrix) <- paste0("Q", probs * 100, "%")
      rownames(quantiles_matrix) <- paste0("Individual_", 1:nrow(quantiles_matrix))

      ## Return the quantile estimates as a data frame
      as.data.frame(quantiles_matrix)
    },
    #' @description Predict survival probabilities at specified time points.
    #' @param new_data A data.frame of covariates.
    #' @param time.points Optional vector of time points. Defaults to fitted model's grid.
    #' @return A list with `predictions` matrix and `time.points`.
    predict = function(new_data, time.points = NULL) {
      new_data <- self$select_columns(new_data)

      if (is.null(time.points)) {
        time.points <- self$time.points
      }

      lp <- predict(self$model, newdata = new_data, type = "lp")
      scale <- self$model$scale
      dist <- self$dist
      n <- nrow(new_data)

      # Define survival function depending on the distribution
      surv_func <- switch(dist,
        "weibull" = function(t, mu, scale) {
          1 - pweibull(t, shape = 1 / scale, scale = exp(mu))
        },
        "lognormal" = function(t, mu, scale) {
          1 - plnorm(t, meanlog = mu, sdlog = scale)
        },
        stop("Unsupported distribution: ", dist)
      )

      # Compute survival probabilities
      surv_matrix <- matrix(NA, nrow = n, ncol = length(time.points))
      for (i in seq_len(n)) {
        surv_matrix[i, ] <- surv_func(time.points, mu = lp[i], scale = scale)
      }

      return(list(predictions = surv_matrix, time.points = time.points))
    },
    #' @description Predict confidence bounds for survival probability.
    #' @param new_data New covariate data.
    #' @param time.points Time grid to evaluate survival curves.
    #' @param level Confidence level (e.g., 0.95).
    #' @param side Side of confidence bound (either upper or lower)
    #' @return Matrix of upper bounds for survival probability.
    predict_confidence_bound = function(new_data, time.points = NULL, level = 0.95, side = c("upper", "lower")) {
        side <- match.arg(side)
        new_data <- self$select_columns(new_data)

        if (is.null(time.points)) {
            time.points <- self$time.points
        }

        z <- qnorm(1 - (1 - level) / 2)
        lp <- predict(self$model, newdata = new_data, type = "lp", se.fit = TRUE)

        scale <- self$model$scale
        dist <- self$dist
        n <- nrow(new_data)
        bound_matrix <- matrix(NA, nrow = n, ncol = length(time.points))

                                        # Select survival function depending on distribution
        surv_func <- switch(dist,
                            "weibull" = function(t, mu, scale) {
                                1 - pweibull(t, shape = 1 / scale, scale = exp(mu))
                            },
                            "lognormal" = function(t, mu, scale) {
                                1 - plnorm(t, meanlog = mu, sdlog = scale)
                            },
                            stop("Unsupported distribution: ", dist)
                            )

        for (i in seq_len(n)) {
            mu <- lp$fit[i]
            se <- lp$se.fit[i]

            mu_ci <- switch(
                side,
                "upper" = mu + z * se,
                "lower" = mu - z * se
            )

            bound_matrix[i, ] <- surv_func(time.points, mu_ci, scale)
        }

        return(list(bound = bound_matrix, time.points = time.points, side = side))
    }

  )
  )

#' @title Cox Proportional Hazards Wrapper
#'
#' @description
#' Wrapper for Cox proportional hazards models via \code{survival::coxph()}.
#'
#' @field model The fitted Cox model object from \code{survival::coxph()}.
#' @field formula_env Environment where the formula is evaluated (optional).
#'
#' @export
CoxphModelWrapper <- R6Class("CoxphModelWrapper",
  inherit = SurvivalModelWrapper,
  public = list(

    model = NULL,             ## To store the fitted coxph model
    formula_env = NULL,        ## To store the formula environment

    #' @description Fit a Cox proportional hazards model using coxph().
    #' @param formula A formula like \code{Surv(time, status) ~ predictors}.
    #' @param data A data.frame with survival times, censoring indicator, and predictors.
    fit = function(formula, data) {
        self$data_train = data
        data <- self$select_columns(data)
        ## Fit the coxph model and store it
        parsed_data <- self$parse_formula(formula, data)
        self$model <- survival::coxph(formula, data = data, x = TRUE, y = TRUE)
        self$time.points <- unique(sort(parsed_data$time))

        ## Capture and store the formula's environment
                                        #self$formula_env <- environment(formula)

        ## Ensure the environment is set correctly for terms evaluation
                                        #environment(self$model$terms) <- self$formula_env
    },

    #' @description Predict survival probabilities for new individuals using Cox model.
    #' @param new_data A data.frame with covariates.
    #' @param time.points Optional numeric vector of time points to predict at.
    #' @return A list with `predictions` matrix and `time.points`.
    predict = function(new_data, time.points = NULL) {
        new_data <- self$select_columns(new_data)
        if (is.null(time.points)) time.points <- self$time.points

                                        # survfit handles baseline + centering consistently with the fitted model
        sf <- survival::survfit(self$model, newdata = new_data, se.fit = FALSE)
                                        # extract survival at requested times, extending flat at ends
        pred_mat <- matrix(NA_real_, nrow = nrow(new_data), ncol = length(time.points))
        for (i in seq_len(nrow(new_data))) {
            s <- summary(sf[i], times = time.points, extend = TRUE)
            pred_mat[i, ] <- s$surv
        }
        list(predictions = pred_mat, time.points = time.points)
    },

    predict_confidence_bound = function(new_data, time.points = NULL,
                                        level = 0.95, side = c("upper","lower"),
                                        conf.type = c("plain","log","log-log","logit","arcsin")) {
        new_data <- self$select_columns(new_data)
        side <- match.arg(side)
        conf.type <- match.arg(conf.type)
        if (is.null(time.points)) time.points <- self$time.points

                                        # Ask survfit for CIs in the same scale you want to report
        sf <- survival::survfit(self$model, newdata = new_data,
                                conf.int = level, conf.type = conf.type, se.fit = (conf.type=="plain"))
        out <- matrix(NA_real_, nrow = nrow(new_data), ncol = length(time.points))

        for (i in seq_len(nrow(new_data))) {
            s <- summary(sf[i], times = time.points, extend = TRUE)
            if (conf.type == "plain") {
                z <- qnorm(1 - (1 - level) / 2)
                                        # std.err from summary() is on S-scale; match your previous math
                b <- switch(side,
                            upper = pmin(1, s$surv + z * s$std.err),
                            lower = pmax(0, s$surv - z * s$std.err))
            } else {
                                        # use survfit’s transformed CIs, then pick a side
                b <- switch(side,
                            upper = s$upper,  # already transformed back to S-scale
                            lower = s$lower)
            }
            out[i, ] <- b
        }
        list(bound = out, time.points = time.points, side = side)
    }

  )
  )

#' @title XGBoost AFT Survival Wrapper
#'
#' @description
#' Wrapper for Accelerated Failure Time (AFT) survival models via
#' \code{xgboost::xgb.train(objective = "survival:aft")}.
#'
#' @details
#' This class fits an AFT model with xgboost and provides methods to
#' predict survival probabilities at given time points and survival
#' quantiles (inverse of the survival curve) for new individuals.
#'
#' @field model A list containing the fitted xgboost booster and
#'   meta-information (features, design-matrix column names, distribution,
#'   and scale parameter).
#' @field time.points Default time grid used for prediction.
#'
#' @export
XgbAftModelWrapper <- R6::R6Class(
  "XgbAftModelWrapper",
  inherit = SurvivalModelWrapper,
  public = list(

    model = NULL,        # list: bst, features, x_colnames, dist, sigma
    time.points = NULL,  # default time grid (from training data)

    #' @description Fit an xgboost AFT survival model.
    #'
    #' @param formula A formula like \code{Surv(time, status) ~ predictors}.
    #' @param data A data.frame with survival times, status (1=event, 0=censor),
    #'   and predictors.
    #' @param params List of xgboost parameters. Must include
    #'   \code{objective = "survival:aft"} and the AFT settings.
    #' @param nrounds Number of boosting iterations.
    #' @param verbose Verbosity level passed to \code{xgboost::xgb.train()}.
    fit = function(formula,
                   data,
                   params = list(
                     objective = "survival:aft",
                     eval_metric = "aft-nloglik",
                     aft_loss_distribution = "logistic",  # "logistic" | "normal" | "extreme"
                     aft_loss_distribution_scale = 1.5,
                     tree_method = "hist",
                     eta = 0.05,
                     max_depth = 6,
                     subsample = 0.8,
                     colsample_bytree = 0.8
                   ),
                   nrounds = 10000,
                   early_stopping_rounds = 50,
                   verbose = 0) {

      self$data_train <- data
      data <- self$select_columns(data)

      # Parse Surv(...) and extract times/status (relies on SurvivalModelWrapper)
      parsed_data <- self$parse_formula(formula, data)
      time   <- parsed_data$time
      status <- parsed_data$status
      self$time.points <- unique(sort(time))

      # Extract predictor variable names from RHS of formula
      terms_obj <- stats::terms(formula, data = data)
      rhs_vars  <- attr(terms_obj, "term.labels")
      if (length(rhs_vars) == 0L) {
        stop("Formula must include predictors on the right-hand side for xgboost AFT.")
      }

      # Design matrix from predictors
      X <- mk_x_matrix_xgb_aft(data, rhs_vars)

      # AFT bounds and DMatrix
      b <- make_aft_bounds_xgb_aft(time, status)
      dtr <- xgboost::xgb.DMatrix(X)
      xgboost::setinfo(dtr, "label_lower_bound", b$lb)
      xgboost::setinfo(dtr, "label_upper_bound", b$ub)

      bst <- xgboost::xgb.train(
        params  = params,
        data    = dtr,
        nrounds = nrounds,
        verbose = verbose
      )

      # Store fitted model and meta-info
      self$model <- list(
        bst        = bst,
        features   = rhs_vars,
        x_colnames = colnames(X),
        dist       = params$aft_loss_distribution,
        sigma      = params$aft_loss_distribution_scale
      )

      invisible(self)
    },

    #' @description Predict survival probabilities for new individuals.
    #'
    #' @param new_data A data.frame with covariates.
    #' @param time.points Optional numeric vector of time points at which to
    #'   predict survival. Defaults to the unique event times from training.
    #' @return A list with `predictions` matrix (nrow = nrow(new_data),
    #'   ncol = length(time.points)) and `time.points`.
    predict = function(new_data, time.points = NULL) {
      if (is.null(self$model)) {
        stop("Model has not been fitted yet. Call $fit() first.")
      }

      new_data <- self$select_columns(new_data)
      if (is.null(time.points)) time.points <- self$time.points

      # Coerce and sanitize time points
      times <- suppressWarnings(as.numeric(time.points))
      times <- times[is.finite(times)]
      if (!length(times)) stop("No valid time.points supplied.")
      times[times < 0] <- 0

      # Build design matrix using training-time feature map
      X <- mk_x_matrix_xgb_aft(
        df           = new_data,
        features     = self$model$features,
        ref_colnames = self$model$x_colnames
      )

      # μ (log-time) from xgboost; outputmargin = TRUE gives raw margin
      mu <- as.numeric(
        predict(self$model$bst, xgboost::xgb.DMatrix(X), outputmargin = TRUE)
      )

      # z_ik = (log t_k - mu_i) / sigma
      lt <- log(times)
      Z  <- outer(mu, lt, function(m, ltk) (ltk - m) / self$model$sigma)

      S <- .Sbar_from_z_xgb_aft(Z, dist = self$model$dist)
      if (is.null(dim(S))) S <- matrix(S, ncol = length(times))
      colnames(S) <- paste0("t=", times)

      list(predictions = S, time.points = times)
    },

    #' @description Predict survival quantiles (inverse survival curve).
    #'
    #' @details
    #' For each individual and each survival probability \eqn{p}, this
    #' method returns the corresponding time \eqn{t_p} such that
    #' \eqn{S(t_p \mid x) = p}.
    #'
    #' @param new_data A data.frame with covariates.
    #' @param probs Numeric vector of survival probabilities
    #'   (e.g., \code{c(0.1, 0.2, 0.5, 0.8, 0.9)}).
    #' @return A list with `quantiles` matrix (nrow = nrow(new_data),
    #'   ncol = length(probs)) and `probs`.
    predict_quantile = function(new_data, probs) {
      if (is.null(self$model)) {
        stop("Model has not been fitted yet. Call $fit() first.")
      }

      new_data <- self$select_columns(new_data)

      # Sanitize probabilities
      probs <- suppressWarnings(as.numeric(probs))
      probs <- probs[is.finite(probs)]
      if (!length(probs)) stop("No valid probs supplied.")

      # Build design matrix using training-time feature map
      X <- mk_x_matrix_xgb_aft(
        df           = new_data,
        features     = self$model$features,
        ref_colnames = self$model$x_colnames
      )

      # μ (log-time) from xgboost
      mu <- as.numeric(
        predict(self$model$bst, xgboost::xgb.DMatrix(X), outputmargin = TRUE)
      )

      # Standardized z_p for each probability p (same for all individuals)
      z_vec <- .z_from_S_xgb_aft(probs, dist = self$model$dist)  # length = length(probs)

      # log t_{i,p} = mu_i + sigma * z_p
      logT <- outer(mu, z_vec, function(m, z) m + self$model$sigma * z)
      Tmat <- exp(logT)

      colnames(Tmat) <- paste0("p=", probs)
      list(quantiles = Tmat, probs = probs)
    }

  )
)

# -------------------------------------------------------------------------
# Helper functions for XgbAftModelWrapper
# -------------------------------------------------------------------------

#' @keywords internal
mk_x_matrix_xgb_aft <- function(df, features, ref_colnames = NULL) {
  X <- stats::model.matrix(~ . - 1, data = df[, features, drop = FALSE])
  if (!is.null(ref_colnames)) {
    miss <- setdiff(ref_colnames, colnames(X))
    if (length(miss)) {
      X <- cbind(
        X,
        matrix(
          0,
          nrow  = nrow(X),
          ncol  = length(miss),
          dimnames = list(NULL, miss)
        )
      )
    }
    X <- X[, ref_colnames, drop = FALSE]
  }
  X
}

#' @keywords internal
make_aft_bounds_xgb_aft <- function(time, status) {
  lb <- time
  ub <- ifelse(status == 1, time, Inf)
  list(lb = lb, ub = ub)
}

#' @keywords internal
.Sbar_from_z_xgb_aft <- function(z, dist) {
  switch(
    dist,
    logistic = stats::plogis(-z),   # 1 - plogis(z) = plogis(-z)
    normal   = stats::pnorm(-z),    # 1 - pnorm(z)  = pnorm(-z)
    extreme  = exp(-exp(z)),        # 1 - (1 - exp(-exp(z))) = exp(-exp(z))
    stop("Unsupported aft_loss_distribution: ", dist)
  )
}

#' @keywords internal
.z_from_S_xgb_aft <- function(p, dist) {
  # p = S(t) in (0, 1) -> standardized z
  eps <- 1e-8
  p <- pmin(pmax(p, eps), 1 - eps)

  switch(
    dist,
    # S(z) = plogis(-z) = 1 / (1 + exp(z))
    # => p = 1 / (1 + exp(z))
    # => exp(z) = (1 - p) / p
    # => z = log((1 - p) / p)
    logistic = log((1 - p) / p),

    # S(z) = pnorm(-z)
    # => p = pnorm(-z)
    # => -z = qnorm(p)
    # => z = -qnorm(p)
    normal   = -stats::qnorm(p),

    # S(z) = exp(-exp(z))
    # => p = exp(-exp(z))
    # => log p = -exp(z)
    # => -log p = exp(z)
    # => z = log(-log p)
    extreme  = log(-log(p)),

    stop("Unsupported aft_loss_distribution: ", dist)
  )
}
