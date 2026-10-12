options(width = 250)

library(tidyverse)

parse_input <- TRUE

source("../R/load_methods.R")

# =============================================================================
# Variable-reduction analysis (Section 4.8 and Supplement S2.2).
# For each model in the set: fit on all variables, rank variables by
# model-specific importance, and refit on the top-k variables. Time-dependent
# AUC on the full test set is computed for every k in num_vars_auc; the full
# pipeline (simultaneous bands + selected thresholds per target yield) only
# for k in num_vars_bands (all variables and 5 variables).
# =============================================================================

# ---- Input parameters --------------------------------------------------------

if (parse_input) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 8) stop("Expected 8 arguments.")
  setup             <- args[1]
  cens_model_type   <- args[2]
  surv_model_type   <- args[3]              # a model-set name, see MODEL_SETS
  num_samples_train <- as.integer(args[4])
  num_samples_cal   <- as.integer(args[5])
  cohort_size       <- as.integer(args[6])
  screening_time    <- as.numeric(args[7])
  batch             <- as.integer(args[8])  # random seed
} else {
  setup <- "real_varsel_1"
  cens_model_type <- "grf"
  surv_model_type <- "grf"
  num_samples_train <- 5000
  num_samples_cal <- 10000
  cohort_size <- 100
  screening_time <- 6
  batch <- 1
}

MODEL_SETS <- list(
  cox  = "cox",
  grf  = "grf",
  rsf  = "rsf",
  xgb  = "xgb",
  mix2 = c("grf", "cox"),
  mix  = c("cox", "grf", "rsf", "xgb")
)
if (!surv_model_type %in% names(MODEL_SETS)) stop("Unknown model set: ", surv_model_type)
model_names <- MODEL_SETS[[surv_model_type]]
K <- length(model_names)

real_data        <- 1
num_samples_test <- cohort_size
alpha            <- 0.05
B_boot           <- 500
lambda_seq       <- seq(0, 1, by = 0.01)
n_min            <- 10       # paper: n_min
pi_min           <- 0.9      # paper: pi_min

target_yield_list <- c(0.02, 0.05, 0.2, 0.3, 0.4, 0.5)
num_vars_bands    <- c(0, 5)                              # full pipeline (bands + selection); 0 = all
num_vars_auc      <- c(0, 10, 7, 5, 3, 2, 1)              # AUC only

# ---- Output ------------------------------------------------------------------

header <- tibble(real_data = real_data, surv_model_type = surv_model_type,
                 cens_model_type = cens_model_type, n_models = K,
                 n_train = num_samples_train, n_cal = num_samples_cal,
                 m = cohort_size, n_test = num_samples_test, alpha = alpha,
                 screening_time = screening_time, batch = batch)

output_file <- sprintf("results/%s/surv_%s_surv_%s_train%d_cal%d_m%d_time%s_batch%d.txt",
                       setup, cens_model_type, surv_model_type, num_samples_train, num_samples_cal,
                       cohort_size, screening_time, batch)
cat("Output file name:", output_file, "\n")

# ---- Data and models ---------------------------------------------------------

load_raw_data <- function() {
  data.path.1 <- "/media/msesia/ExternalSSD/data/fhrd/train_test.rds"
  data.path.2 <- "/project2/sesia_1123/data/exploration/train_test.rds"
  data.path <- if (file.exists(data.path.1)) data.path.1 else data.path.2
  combined <- readRDS(data.path)
  combined <- as_tibble(do.call("rbind", combined)) %>%
    mutate(time = month, status = event) %>%
    select(-month, -event) %>%
    select(time, status, everything())
  # Keep the original variable names for reporting
  var_names <<- tibble(code = paste0("X", seq_len(ncol(combined) - 2)),
                       variable = setdiff(colnames(combined), c("time", "status")))
  colnames(combined) <- c("time", "status", var_names$code)
  combined
}
data.raw <- load_raw_data()

data.generator <- RealDataGenerator$new(data = data.raw)

cens_model <- CensoringModel$new(model = init_censoring_model(cens_model_type))

# ---- One replication ---------------------------------------------------------

analyze_data <- function(data.train, data.cal, data.test) {

  cens_model$fit(data = data.train)

  ## Calibration outcomes and IPCW inputs (do not depend on the survival model)
  prob_O.cal <- compute_censoring_survival(data.cal, cens_model, data.cal$time, data.cal$status, method = "et")
  O.cal <- as.integer(data.cal$status)
  E.cal <- as.integer(data.cal$time <= screening_time)
  E.cal[O.cal == 0L] <- NA_integer_

  ## ---- Full fit and model-specific variable ranking ------------------------
  full_models <- setNames(lapply(model_names, init_surv_model), model_names)
  for (k in model_names) full_models[[k]]$fit(Surv(time, status) ~ ., data = data.train)

  ranking <- map_dfr(model_names, function(k) {
    imp <- sort(full_models[[k]]$importance(), decreasing = TRUE)
    tibble(model = k, rank = seq_along(imp), code = names(imp), importance = as.numeric(imp)) %>%
      left_join(var_names, by = "code")
  })

  ## ---- Bands, selection and AUC for a given set of fitted models -----------
  sim_bands <- function(scores) {
    bootstrap_simultaneous(scores = scores, E = E.cal, O = O.cal, prob_O = prob_O.cal,
                           B = B_boot, m = cohort_size, alpha = alpha, seed = 1,
                           fs_smoothing = TRUE, ipcw_method = "ht",
                           lambda_seq = lambda_seq, n_min = n_min, pi_min = pi_min)
  }

  sim_to_long <- function(estim, band_method, models) {
    map_dfr(seq_along(models), function(j) {
      estim_j <- if (length(models) == 1) estim else unpack_model(estim, j)
      est <- tibble(lambda = estim_j$estimate$lambda,
                    yield  = estim_j$estimate$pS,
                    ppv    = 1 - estim_j$estimate$pE)
      convert_band_to_long_targets(estim = estim_j, band_method = band_method, m = cohort_size) %>%
        left_join(est, by = "lambda") %>%
        mutate(model = models[j], valid_models = paste(models, collapse = "|"))
    })
  }

  valid_set <- paste(model_names, collapse = "|")
  selection <- function(bands, target_yield) {                 # rule (7)
    bands %>%
      filter(band_type=="prediction", yield_lo >= target_yield, target_m == cohort_size, valid_models == valid_set) %>%
      slice_max(lambda, n = 1, with_ties = FALSE) %>%
      mutate(target_yield = target_yield)
  }

  ## Bands and threshold selection (full pipeline)
  run_bands <- function(models, num_vars) {
    scores.cal <- lapply(models, function(mod) as.numeric(mod$predict(data.cal, screening_time)$predictions))
    bands <- sim_to_long(sim_bands(scores.cal), "Bootstrap-sim", model_names) %>%
      mutate(t0 = screening_time)
    map_dfr(target_yield_list, ~ selection(bands, .x)) %>% mutate(num_vars = num_vars)
  }

  ## Time-dependent AUC (IPCW, cumulative/dynamic) on the full test set.
  run_auc <- function(models, num_vars) {
    scores.test <- lapply(models, function(mod) as.numeric(mod$predict(data.test, screening_time)$predictions))
    map_dfr(model_names, function(k) {
      roc <- timeROC::timeROC(T = data.test$time, delta = data.test$status,
                              marker = -scores.test[[k]],   # marker increases with risk
                              cause = 1, times = screening_time, iid = FALSE)
      tibble(model = k, auc = as.numeric(roc$AUC[2]))
    }) %>% mutate(num_vars = num_vars)
  }

  ## ---- Loop over predictor-set sizes ---------------------------------------
  ## AUC for every size in num_vars_auc; bands + selection only for num_vars_bands.
  fit_top_k <- function(nv) {
    if (nv == 0) return(full_models)
    setNames(lapply(model_names, function(k) {
      top_k <- ranking %>% filter(model == k, rank <= nv) %>% pull(code)
      mod <- init_surv_model(k, use_covariates = top_k)
      mod$fit(Surv(time, status) ~ ., data = data.train)
      mod
    }), model_names)
  }

  selections <- list(); aucs <- list()
  for (nv in sort(unique(c(num_vars_bands, num_vars_auc)), decreasing = TRUE)) {
    cat(sprintf("num_vars = %d\n", nv))
    models <- fit_top_k(nv)
    if (nv %in% num_vars_auc)   aucs[[length(aucs) + 1]]             <- run_auc(models, nv)
    if (nv %in% num_vars_bands) selections[[length(selections) + 1]] <- run_bands(models, nv)
  }

  ## ---- Write -----------------------------------------------------------------
  write.csv(bind_cols(header, bind_rows(selections)), output_file, row.names = FALSE)
  write.csv(bind_cols(header, bind_rows(aucs)),
            str_replace(output_file, "\\.txt$", "_auc.txt"), row.names = FALSE)
  write.csv(bind_cols(header, ranking),
            str_replace(output_file, "\\.txt$", "_importance.txt"), row.names = FALSE)
  cat(sprintf("Written results on %s (+ _auc, _importance).\n", output_file))
}

run_experiment <- function(random.state) {
  set.seed(random.state)
  oracle <- data.generator$sample(shuffle = TRUE, return_oracle = TRUE)$oracle
  splits <- split_data_n(oracle, n_train = num_samples_train, n_cal = num_samples_cal,
                         n_test = num_samples_test)
  drop_oracle <- function(d) d %>% select(-event_time, -censoring_time)
  analyze_data(drop_oracle(splits$data.train), drop_oracle(splits$data.cal), drop_oracle(splits$data.test))
}

## ---- Run ---------------------------------------------------------------------
run_experiment(batch)
