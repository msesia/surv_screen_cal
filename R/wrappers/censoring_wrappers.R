library(R6)

#' #' Initialize a Censoring Model
#'
#' Returns a survival model wrapper suitable for modeling censoring distributions.
#' The model inverts the event indicator and supports optional covariate selection.
#'
#' @param model_type A string specifying the model type. One of `"grf"`, `"survreg"`, `"rf"`, or `"cox"`.
#' @param use_covariates Optional character vector of covariate names to include (default: `NULL`, uses all).
#'
#' @return An instance of a \code{SurvivalModelWrapper} subclass for censoring modeling.
#' @export
init_censoring_model <- function(model_type, use_covariates = NULL) {
  cens_model <- switch(model_type,
    "grf" = GRF_SurvivalForestWrapper$new(use_covariates = use_covariates),
    "survreg" = SurvregModelWrapper$new(dist = "lognormal", use_covariates = use_covariates),
    "rf" = randomForestSRC_SurvivalWrapper$new(use_covariates = use_covariates),
    "cox" = CoxphModelWrapper$new(use_covariates = use_covariates),
    "xgb" = XgbAftModelWrapper$new(use_covariates = use_covariates),
    stop("Unknown censoring model type!")
  )
  return(cens_model)
}

#' CensoringModel Class
#'
#' Wraps a `SurvivalModelWrapper` to model censoring distributions by flipping event indicators.
#'
#' @description
#' This class inherits from `SurvivalModelWrapper` and wraps any compatible survival model.
#' It fits a censoring distribution (i.e., the survival function of the censoring time)
#' by inverting the event indicator: `status := 1 - status`.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(model)}}{Creates a new `CensoringModel` with a specified survival model wrapper.}
#'   \item{\code{fit(data)}}{Fits the censoring distribution using the flipped event indicators.}
#'   \item{\code{predict(new_data, time.points = NULL)}}{Predicts the censoring survival curve.}
#' }
#'
#' @export
CensoringModel <- R6::R6Class("CensoringModel",
  inherit = SurvivalModelWrapper,
  public = list(

    #' @description Constructor for CensoringModel.
    #' @param model A survival model wrapper (e.g., output of `init_censoring_model()`).
    initialize = function(model) {
      self$model <- model
    },

    #' @description Fit the censoring model using flipped event indicators.
    #' @param data A data frame containing time, status, and covariates.
    fit = function(data) {
      data_cens <- data
      data_cens$status <- 1 - data_cens$status
      self$model$fit(Surv(time, status) ~ ., data = data_cens)
    },

    #' @description Predict the censoring survival function.
    #' @param new_data New covariate data for prediction.
    #' @param time.points Optional time grid for survival curve evaluation.
    #' @return A list with predicted survival probabilities and time points.
    predict = function(new_data, time.points = NULL) {
      self$model$predict(new_data, time.points = time.points)
    }
  )
)
