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
  if (length(args) < 9) stop("Expected 9 arguments.")
  setup             <- args[1]
  real_data         <- as.integer(args[2])
  gen_model_type    <- args[3]
  surv_model_type   <- args[4]              # a model-set name, see MODEL_SETS
  num_samples_train <- as.integer(args[5])
  num_samples_cal   <- as.integer(args[6])
  num_samples_test  <- as.integer(args[7])
  screening_time    <- as.numeric(args[8])
  batch             <- as.integer(args[9])
  surv_nph          <- if (length(args) >= 10) as.numeric(args[10]) else 0   # see setup "v0n"
} else {
  setup <- "v0"
  real_data <- 0
  gen_model_type <- "grf"
  surv_model_type <- "grf"
  num_samples_train <- 5000
  num_samples_cal <- 100
  num_samples_test <- 100
  screening_time <- 3
  batch <- 1
  surv_nph <- 0
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
                 n_models = K,
                 n_train = num_samples_train, n_cal = num_samples_cal,
                 n_test = num_samples_test, alpha = alpha,
                 screening_time = screening_time, batch = batch,
                 surv_nph = surv_nph)

nph_tag <- if (surv_nph != 0) sprintf("_nph%g", surv_nph) else ""
output_file <- sprintf("results/%s/real_%d_gen_%s_surv_%s_train%d%s_cal%d_test%d_time%s_batch%d.txt",
                       setup, real_data, gen_model_type, surv_model_type,
                       num_samples_train, nph_tag, num_samples_cal, num_samples_test,
                       screening_time, batch)

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

# ---- Setup "v0n": non-proportional survival hazards of controllable strength -
# The fitted survival hazard is multiplied by exp(surv_nph * z * (1 - 2 t / tau)),
# z = standardized serum albumin (X10): albumin raises the hazard early and
# lowers it late (crossing hazards). A Cox model, which assumes a constant
# log hazard ratio, becomes increasingly miscalibrated at t0 as surv_nph
# grows; surv_nph = 0 reproduces the fitted generator. Censoring is as fitted.
z_mean <- mean(data.raw$X10)
z_sd <- sd(data.raw$X10)
t_flip <- screening_time
tilt_surv_nph <- list(f = ~ surv_nph * ((X10 - z_mean) / z_sd) * (1 - 2 * pmin(t, t_flip) / t_flip))

data.generator <- if (real_data) {
  RealDataGenerator$new(data = data.raw)
} else if (setup == "v0n") {
  TiltedSemiSyntheticDataGenerator$new(data = data.raw, surv_model_type = gen_model_type,
                                       cens_model_type = gen_model_type,
                                       tilt_surv = tilt_surv_nph, tilt_cens = NULL)
} else {
  SemiSyntheticDataGenerator$new(data = data.raw, surv_model_type = gen_model_type,
                                 cens_model_type = gen_model_type)
}


surv_models <- setNames(lapply(model_names, init_surv_model), model_names)
cens_model  <- CensoringModel$new(model = init_censoring_model(cens_model_type))

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
  sim_to_long <- function(estim, band_method, models) {
    map_dfr(seq_along(models), function(j) {
      estim_j <- if (length(models) == 1) estim else unpack_model(estim, j)
      convert_band_to_long_targets(estim = estim_j, band_method = band_method, m = num_samples_test) %>%
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
  analyze_data(drop_oracle(splits$data.train), drop_oracle(splits$data.cal),
               drop_oracle(splits$data.test), splits$data.test)
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
