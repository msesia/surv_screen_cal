library(tidyverse)

parse_input <- TRUE

source("../R/load_methods.R")

# =============================================================================
# Computational cost of the simultaneous bootstrap (Algorithms 1-3) as a
# function of the calibration size n and the cohort size m, on real FHRD
# data. The survival and censoring models are fitted once on a training
# split and their cost is NOT included; the timed part is the bootstrap
# itself (with and without the finite-sample adjustment) on calibration
# data of size n, for a cohort of size m.
# =============================================================================

# ---- Input parameters --------------------------------------------------------

if (parse_input) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 5) stop("Expected 5 arguments: setup surv_model n_cal m batch")
  setup           <- args[1]
  surv_model_type <- args[2]              # a model-set name, see MODEL_SETS
  num_samples_cal <- as.integer(args[3])
  cohort_size     <- as.integer(args[4])
  batch           <- as.integer(args[5])
} else {
  setup <- "timing"
  surv_model_type <- "grf"
  num_samples_cal <- 1000
  cohort_size <- 100
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

num_samples_train <- 5000
screening_time    <- 3
alpha             <- 0.05
cens_model_type   <- "grf"
B_boot            <- 500
batch_size        <- 5          # timed repetitions (fresh calibration subsample each)
lambda_seq        <- seq(0, 1, by = 0.01)
n_min             <- 10
pi_min            <- 0.9

# ---- Output ------------------------------------------------------------------

header <- tibble(surv_model_type = surv_model_type, n_models = K,
                 cens_model_type = cens_model_type, n_train = num_samples_train,
                 n_cal = num_samples_cal, m = cohort_size, B = B_boot,
                 n_lambda = length(lambda_seq), screening_time = screening_time,
                 batch = batch)

output_file <- sprintf("results/%s/surv_%s_cal%d_m%d_batch%d.txt",
                       setup, surv_model_type, num_samples_cal, cohort_size, batch)
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
cat("Output file name:", output_file, "\n")

# ---- Data and models (fitted once, not timed) --------------------------------

load_raw_data <- function() {
  data.path.1 <- "/media/msesia/ExternalSSD/data/fhrd/train_test.rds"
  data.path.2 <- "/project2/sesia_1123/data/exploration/train_test.rds"
  data.path <- if (file.exists(data.path.1)) data.path.1 else data.path.2
  combined <- readRDS(data.path)
  combined <- as_tibble(do.call("rbind", combined)) %>%
    mutate(time = month, status = event) %>%
    select(-month, -event) %>%
    select(time, status, everything())
  colnames(combined) <- c("time", "status", paste0("X", seq_len(ncol(combined) - 2)))
  combined
}
data.raw <- load_raw_data()

set.seed(batch)
data.all   <- data.raw[sample(nrow(data.raw)), ]
data.train <- data.all[seq_len(num_samples_train), ]
data.pool  <- data.all[-seq_len(num_samples_train), ]      # calibration data drawn from here
stopifnot(nrow(data.pool) >= num_samples_cal)

surv_models <- setNames(lapply(model_names, init_surv_model), model_names)
for (k in model_names) surv_models[[k]]$fit(Surv(time, status) ~ ., data = data.train)
cens_model <- CensoringModel$new(model = init_censoring_model(cens_model_type))
cens_model$fit(data = data.train)

# ---- Timed replication -------------------------------------------------------

time_one <- function(random.state) {
  set.seed(random.state)
  data.cal <- data.pool[sample(nrow(data.pool), num_samples_cal), ]

  # Inputs to the bootstrap (scores, IPCW quantities): computed once, not timed
  scores.cal <- lapply(surv_models, function(mod) as.numeric(mod$predict(data.cal, screening_time)$predictions))
  prob_O.cal <- compute_censoring_survival(data.cal, cens_model, data.cal$time, data.cal$status, method = "et")
  O.cal <- as.integer(data.cal$status)
  E.cal <- as.integer(data.cal$time <= screening_time)
  E.cal[O.cal == 0L] <- NA_integer_

  run <- function(fs_smoothing) {
    bootstrap_simultaneous(scores = scores.cal, E = E.cal, O = O.cal, prob_O = prob_O.cal,
                           B = B_boot, m = cohort_size, alpha = alpha, seed = 1,
                           fs_smoothing = fs_smoothing, ipcw_method = "ht",
                           lambda_seq = lambda_seq, n_min = n_min, pi_min = pi_min)
  }

  t_nofs <- system.time(run(FALSE))[["elapsed"]]   # Algorithms 1-2 only
  t_fs   <- system.time(run(TRUE))[["elapsed"]]    # with Algorithm 3 (three internal runs)

  tibble(seed = random.state, time_sim_nofs = t_nofs, time_sim_fs = t_fs)
}

# ---- Run ---------------------------------------------------------------------

results_df <- tibble()
for (i in seq_len(batch_size)) {
  random.state <- batch * 1000 + i
  results_df <- bind_rows(results_df, bind_cols(header, time_one(random.state)))
  write.csv(results_df, output_file, row.names = FALSE)
  cat(sprintf("n = %d, m = %d, rep %d: %.1f s (no FS), %.1f s (FS)\n",
              num_samples_cal, cohort_size, i,
              results_df$time_sim_nofs[i], results_df$time_sim_fs[i]))
}
