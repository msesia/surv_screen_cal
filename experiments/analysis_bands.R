library(tidyverse)

parse_input <- TRUE

source("../R/load_methods.R")

# ---- Input parameters --------------------------------------------------------

if (parse_input) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 7) stop("Expected 6 arguments.")
  setup             <- args[1]
  cens_model_type   <- args[2]
  surv_model_type   <- args[3]              # a model-set name, see MODEL_SETS
  num_samples_train <- as.integer(args[4])
  num_samples_cal   <- as.integer(args[5])
  cohort_size       <- as.integer(args[6])
  screening_time    <- as.numeric(args[7])
} else {
  setup <- "v0"  
  cens_model_type <- "grf"
  surv_model_type <- "grf"
  num_samples_train <- 5000
  num_samples_cal <- 10000
  cohort_size <- 100
  screening_time <- 6
}


## LOOP OVER surv_model_type in "grf" "mix2"

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

real_data <- 1
num_samples_test <- cohort_size
alpha            <- 0.05
B_boot           <- 500
lambda_seq       <- seq(0, 1, by = 0.01)
n_min            <- 10       # paper: n_min
pi_min           <- 0.9      # paper: pi_min

# ---- Output ------------------------------------------------------------------

header <- tibble(real_data = real_data, surv_model_type = surv_model_type,
                 cens_model_type = cens_model_type, n_models = K,
                 n_train = num_samples_train, n_cal = num_samples_cal,
                 m = cohort_size, n_test = num_samples_test, alpha = alpha,
                 screening_time = screening_time)


output_file <- sprintf("results/%s/cens_%s_surv_%s_train%d_cal%d_m%d_time%s.txt",
                       setup, cens_model_type, surv_model_type, num_samples_train, num_samples_cal, cohort_size, screening_time)
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

data.generator <- if (real_data) {
  RealDataGenerator$new(data = data.raw)
} else {
  SemiSyntheticDataGenerator$new(data = data.raw, surv_model_type = gen_model_type,
                                 cens_model_type = gen_model_type)
}

surv_models <- setNames(lapply(model_names, init_surv_model), model_names)
cens_model  <- CensoringModel$new(model = init_censoring_model(cens_model_type))

# ---- One replication ---------------------------------------------------------

analyze_data <- function(data.train, data.cal, data.test, data.test.oracle) {

    ## Fit K survival models and the censoring model on training data
    for (k in model_names) surv_models[[k]]$fit(Surv(time, status) ~ ., data = data.train)
    cens_model$fit(data = data.train)

    ## One score vector per model
    scores.cal  <- lapply(surv_models, function(mod) as.numeric(mod$predict(data.cal,  screening_time)$predictions))
    scores.test <- lapply(surv_models, function(mod) as.numeric(mod$predict(data.test, screening_time)$predictions))

    ## Calibration outcomes and IPCW inputs
    prob_O.cal <- compute_censoring_survival(data.cal, cens_model, data.cal$time, data.cal$status, method = "et")
    O.cal <- as.integer(data.cal$status)
    E.cal <- as.integer(data.cal$time <= screening_time)
    E.cal[O.cal == 0L] <- NA_integer_

    ## ---- Band methods ----
    sim_bands <- function(scores, m, fs_smoothing = TRUE) {
        bootstrap_simultaneous(scores = scores, E = E.cal, O = O.cal, prob_O = prob_O.cal,
                               B = B_boot, m = m, alpha = alpha, seed = 1,
                               fs_smoothing = fs_smoothing, ipcw_method = "ht",
                               lambda_seq = lambda_seq, n_min = n_min, pi_min = pi_min)
    }

                                        # Long format with bounds AND point estimates; target_m is the cohort
                                        # size the prediction bands were computed for.
    sim_to_long <- function(estim, band_method, models, m) {
        map_dfr(seq_along(models), function(j) {
            estim_j <- if (length(models) == 1) estim else unpack_model(estim, j)
            est <- tibble(lambda = estim_j$estimate$lambda,
                          yield  = estim_j$estimate$pS,
                          ppv    = 1 - estim_j$estimate$pE)
            convert_band_to_long_targets(estim = estim_j, band_method = band_method, m = m) %>%
                left_join(est, by = "lambda") %>%
                mutate(model = models[j], valid_models = paste(models, collapse = "|"))
        })
    }

    ## (a) Simultaneous bands for m = 100 and m = 1000; the confidence band
    ##     does not depend on m, so it is kept from the first call only.
    bands <- sim_to_long(sim_bands(scores.cal, cohort_size,  TRUE), "Bootstrap-sim", model_names, m = cohort_size) %>%
        mutate(t0 = screening_time)

    write.csv(bands, output_file, row.names = FALSE)

    cat(sprintf("Written results on %s.\n", output_file))        
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

## ---- Run ---------------------------------------------------------------------
run_experiment(1)
