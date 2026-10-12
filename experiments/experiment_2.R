library(tidyverse)

parse_input <- TRUE

if (parse_input) {
  source("../R/load_methods.R")
} else {
  source("../R/utils/utils_splitting.R")
  source("../R/utils/utils_weights_scores.R")
  source("../R/utils/utils_misc.R")
  source("../R/utils/utils_plotting.R")
  source("../R/utils/utils_imputation.R")
  source("../R/utils/utils_semi_synthetic_data.R")
  source("../R/utils/utils_evaluation.R")
  source("../R/wrappers/censoring_wrappers.R")
  source("../R/wrappers/survival_wrappers.R")
  source("../R/methods/calibration_bands.R")
  source("../R/methods/pointwise_bootstrap.R")
  source("../R/methods/simultaneous_bootstrap.R")
}

# ---- Input parameters --------------------------------------------------------

if (parse_input) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 9) stop("Expected 9 to 11 arguments.")
  setup             <- args[1]
  real_data         <- as.integer(args[2])
  gen_model_type    <- args[3]
  surv_model_type   <- args[4]              # a model-set name, see MODEL_SETS
  num_samples_train <- as.integer(args[5])
  num_samples_cal   <- as.integer(args[6])
  num_samples_test  <- as.integer(args[7])
  screening_time    <- as.numeric(args[8])
  batch             <- as.integer(args[9])
  cens_tilt         <- if (length(args) >= 10) as.numeric(args[10]) else 0   # see setup "v0t"
  cens_shift        <- if (length(args) >= 11) as.numeric(args[11]) else 0   # see setup "v0t"
} else {
  setup <- "v0t"
  real_data <- 0
  gen_model_type <- "grf"
  surv_model_type <- "grf"
  num_samples_train <- 5000
  num_samples_cal <- 100
  num_samples_test <- 100
  screening_time <- 3
  batch <- 1
  cens_tilt <- 1
  cens_shift <- log(5)
}

# Model sets: surv_model_type names a set; K = number of models in the set.
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

alpha           <- 0.05
cens_model_type <- "grf"
B_boot          <- 500
batch_size      <- 20
lambda_seq      <- seq(0, 1, by = 0.01)
n_min           <- 10       # paper: n_min
pi_min          <- 0.9      # paper: pi_min

# ---- Output ------------------------------------------------------------------

header <- tibble(real_data = real_data, gen_model_type = gen_model_type,
                 surv_model_type = surv_model_type, cens_model_type = cens_model_type,
                 n_models = K, cens_tilt = cens_tilt, cens_shift = cens_shift,
                 cens_hidden = as.integer(setup == "v0h"),
                 n_train = num_samples_train, n_cal = num_samples_cal,
                 n_test = num_samples_test, alpha = alpha,
                 screening_time = screening_time, batch = batch)

# Under "v0t" every file is tagged with its tilt strength (including 0, the
# control point); other setups keep their names. Must match the launcher.
mix_tag <- sprintf("_tilt%g_shift%g", cens_tilt, cens_shift)

output_file <- sprintf("results/%s/real_%d_gen_%s_surv_%s_train%d%s_cal%d_test%d_time%s_batch%d.txt",
                       setup, real_data, gen_model_type, surv_model_type,
                       num_samples_train, mix_tag, num_samples_cal, num_samples_test,
                       screening_time, batch)
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
cat("Output file name:", output_file, "\n")

# ---- Data and models ---------------------------------------------------------

load_raw_data <- function() {
    data.path.1 <- "/media/msesia/ExternalSSD/data/fhrd/train_test.rds"
    data.path.2 <- "/project2/sesia_1123/data/exploration/train_test.rds"
    if(file.exists(data.path.1)) {
        data.path <- data.path.1
    } else {
        data.path <- data.path.2
    }
    combined <- readRDS(data.path)
    combined <- as_tibble(do.call("rbind", combined)) %>%
        mutate(time = month, status = event) %>%
        select(-month, -event) %>%
        select(time, status, everything())
    colnames(combined) <- c("time", "status", paste0("X", seq_len(ncol(combined) - 2)))
    combined
}
data.raw <- load_raw_data()

# ---- Setups "v0t" / "v0h": censoring mechanism of controllable difficulty ---
# The analysis data keep the fitted semi-synthetic SURVIVAL law (no tilt) but
# CENSORING is tilted by a proportional-hazards factor
#     exp(cens_shift + cens_tilt * z),   z = standardized `tilt_var`.
# cens_shift raises the overall censoring rate so that censoring before t0 is
# common enough for the IPCW correction to matter; cens_tilt makes the
# censoring mechanism depend more strongly on `tilt_var` (hazard ratios of
# ~1.6, 2.7, 4.5 per SD for cens_tilt = 0.5, 1, 1.5).
#   v0t: tilt on X32 (log bilirubin); the censoring model sees all covariates,
#        so it can recover the dependence (censoring harder, but estimable).
#   v0h: tilt on X10 (serum albumin, the strongest prognostic variable); the
#        censoring model is fitted WITHOUT X10, i.e. informative censoring
#        driven by a variable the censoring model does not observe, so G_hat
#        is increasingly misspecified as cens_tilt grows.
# cens_tilt = cens_shift = 0 reproduces the untilted generator in both cases.
tilt_var    <- if (setup == "v0h") "X10" else "X32"
cens_hidden <- setup == "v0h"

z_mean <- mean(data.raw[[tilt_var]]); z_sd <- sd(data.raw[[tilt_var]])
tilt_cens_analysis <- list(
  f = as.formula(sprintf("~ cens_shift + cens_tilt * (%s - z_mean) / z_sd", tilt_var))
)   # cens_shift, cens_tilt, z_mean, z_sd are resolved in the global environment

if (real_data) {
  data.generator <- RealDataGenerator$new(data = data.raw)
} else if (setup %in% c("v0t", "v0h")) {
  data.generator <- TiltedSemiSyntheticDataGenerator$new(
    data = data.raw, surv_model_type = gen_model_type, cens_model_type = gen_model_type,
    tilt_surv = NULL, tilt_cens = tilt_cens_analysis)
  cat(sprintf("Setup %s: censoring log-hazard shift %.2f, tilt %.2f per SD of %s (%s)\n",
              setup, cens_shift, cens_tilt, tilt_var,
              if (cens_hidden) "hidden from the censoring model" else "visible to the censoring model"))
} else {
  data.generator <- SemiSyntheticDataGenerator$new(data = data.raw, surv_model_type = gen_model_type,
                                                   cens_model_type = gen_model_type)
}

surv_models <- setNames(lapply(model_names, init_surv_model), model_names)

# Censoring model: all covariates, except that under "v0h" the tilt variable
# is withheld
all_vars  <- setdiff(colnames(data.raw), c("time", "status"))
cens_vars <- if (cens_hidden) setdiff(all_vars, tilt_var) else all_vars
cens_model <- CensoringModel$new(model = init_censoring_model(cens_model_type, use_covariates = cens_vars))

# ---- One replication ---------------------------------------------------------

analyze_data <- function(data.train, data.cal, data.test, data.test.oracle) {

  # Fit K survival models and the censoring model on training data
  for (k in model_names) surv_models[[k]]$fit(Surv(time, status) ~ ., data = data.train)
  cens_model$fit(data = data.train)

  # One score vector per model
  scores.cal  <- lapply(surv_models, function(mod) as.numeric(mod$predict(data.cal,  screening_time)$predictions))
  scores.test <- lapply(surv_models, function(mod) as.numeric(mod$predict(data.test, screening_time)$predictions))

  # Calibration outcomes and IPCW inputs
  prob_O.cal <- compute_censoring_survival(data.cal, cens_model, data.cal$time, data.cal$status, method = "et")
  O.cal <- as.integer(data.cal$status)
  E.cal <- as.integer(data.cal$time <= screening_time)
  E.cal[O.cal == 0L] <- NA_integer_

  # ---- Oracle truth on the test cohort, per (model, lambda) ----
  df_eval <- map_dfr(model_names, function(k) {
    map_dfr(lambda_seq, function(lambda) {
      sel  <- which(scores.test[[k]] >= lambda)
      eval <- evaluate_yield_and_event_rate(data.test.oracle, selected_idx = sel, t0 = screening_time)
      tibble(model = k, lambda = lambda, num_sel = eval$num_selected,
             yield = eval$yield, ppv = 1 - eval$event_rate)
    })
  })

  # ---- Band methods ----
  # Each row records `model` (the model defining the cell) and `valid_models`
  # (the set of models over which the band is simultaneously valid).
  sim_bands <- function(scores, fs_smoothing = TRUE) {
    bootstrap_simultaneous(scores = scores, E = E.cal, O = O.cal, prob_O = prob_O.cal,
                           B = B_boot, m = num_samples_test, alpha = alpha, seed = 1,
                           fs_smoothing = fs_smoothing, ipcw_method = "ht",
                           lambda_seq = lambda_seq, n_min = n_min, pi_min = pi_min)
  }
  # Long format with bounds AND point estimates (yield_hat, ppv_hat)
  sim_to_long <- function(estim, band_method, models) {
    map_dfr(seq_along(models), function(j) {
      estim_j <- if (length(models) == 1) estim else unpack_model(estim, j)
      est <- tibble(lambda    = estim_j$estimate$lambda,
                    yield_hat = estim_j$estimate$pS,
                    ppv_hat   = 1 - estim_j$estimate$pE)
      convert_band_to_long_targets(estim = estim_j, band_method = band_method, m = num_samples_test) %>%
        left_join(est, by = "lambda") %>%
        mutate(model = models[j], valid_models = paste(models, collapse = "|"))
    })
  }

  # (a) Simultaneous bands, jointly valid over all K models and thresholds
  bands_joint      <- sim_to_long(sim_bands(scores.cal, TRUE),  "Bootstrap-sim",          model_names)
  bands_joint_nofs <- sim_to_long(sim_bands(scores.cal, FALSE), "Bootstrap-sim-baseline", model_names)

  # (b) Simultaneous bands per model (valid over thresholds only); skipped when K = 1
  bands_single <- if (K > 1) {
    map_dfr(model_names, function(k) sim_to_long(sim_bands(scores.cal[[k]], TRUE), "Bootstrap-sim", k))
  } else tibble()

  # (c) Pointwise bootstrap baseline, per cell
  bands_point <- map_dfr(model_names, function(k) {
    band_df <- map_dfr(lambda_seq, function(lambda) {
      S <- as.integer(scores.cal[[k]] >= lambda)
      estim <- bootstrap_pointwise(S = S, E = E.cal, O = O.cal, prob_O = prob_O.cal,
                                   B = B_boot, m = num_samples_test, alpha = alpha, seed = 1,
                                   ipcw_method = "ht", n_min = n_min, pi_min = pi_min)
      make_lambda_row(lambda, S, estim, E.cal, O.cal, prob_O.cal, num_samples_test)
    })
    convert_pointwise_to_long_targets(band_df = band_df, band_method = "Bootstrap-point",
                                      m = num_samples_test) %>%
      mutate(model = k, valid_models = k)
  })

  df_bands <- bind_rows(bands_joint, bands_joint_nofs, bands_single, bands_point)
  left_join(df_eval, df_bands, by = c("model", "lambda"))
}

run_experiment <- function(random.state) {
  set.seed(random.state)
  oracle <- data.generator$sample(shuffle = TRUE, return_oracle = TRUE)$oracle
  splits <- split_data_n(oracle, n_train = num_samples_train, n_cal = num_samples_cal,
                         n_test = num_samples_test)
  drop_oracle <- function(d) d %>% select(-event_time, -censoring_time)

  # Censoring summary of the calibration data (recorded with the results):
  # overall censored fraction, and fraction censored before t0
  cens_summary <- tibble(
    frac_censored    = mean(splits$data.cal$status == 0),
    frac_censored_t0 = mean(splits$data.cal$status == 0 & splits$data.cal$time <= screening_time)
  )

  res <- analyze_data(drop_oracle(splits$data.train), drop_oracle(splits$data.cal),
                      drop_oracle(splits$data.test), splits$data.test)
  bind_cols(cens_summary, res)
}

# ---- Run ---------------------------------------------------------------------

results_df <- tibble()
pb <- txtProgressBar(min = 0, max = batch_size, style = 3)
for (i in seq_len(batch_size)) {
  random.state <- batch * 1000 + i
  res <- run_experiment(random.state)
  results_df <- bind_rows(results_df, bind_cols(tibble(Seed = random.state), header, res))
  write.csv(results_df, output_file, row.names = FALSE)
  setTxtProgressBar(pb, i)
}
close(pb)
