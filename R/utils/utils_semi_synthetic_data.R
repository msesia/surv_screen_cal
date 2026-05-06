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
      if (shuffle) {
        df <- self$data_raw[sample(n), ]
      }
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

library(R6)
library(survival)
library(dplyr)
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

## # Augment dataset with permuted "noise" copies of original variables
## # - Do this *before* any training.
## # - Excludes outcome columns (default: time, status).
## augment_noise_vars <- function(data,
##                                copies = 1,
##                                exclude = c("time", "status"),
##                                joint_permute = TRUE,   # one permutation for all vars (like your code)
##                                suffix = "c",
##                                seed = NULL) {
##   stopifnot(is.data.frame(data), copies >= 1)
##   if (!is.null(seed)) set.seed(as.integer(seed))

##   n <- nrow(data)
##   out <- data

##   # variables to copy (all predictors by default)
##   var_names <- setdiff(colnames(data), exclude)
##   if (length(var_names) == 0L) return(out)

##   for (rep in seq_len(copies)) {
##     if (joint_permute) {
##       # one permutation applied to all variables -> preserves their joint structure
##       idx <- sample.int(n)
##       new_df <- data[idx, var_names, drop = FALSE]
##     } else {
##       # independent permutation per variable
##       new_df <- as.data.frame(lapply(data[var_names], function(col) col[sample.int(n)]),
##                               stringsAsFactors = FALSE)
##     }

##     colnames(new_df) <- paste(var_names, suffix, rep, sep = ".")
##     out <- cbind(out, new_df, stringsAsFactors = FALSE)
##   }

##   out
## }
