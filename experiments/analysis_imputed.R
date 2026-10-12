library(tidyverse)

parse_input <- TRUE

source("../R/load_methods.R")

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
  imputed           <- as.logical(as.numeric(args[8]))
} else {
  setup <- "v0"
  cens_model_type <- "grf"
  surv_model_type <- "grf"
  num_samples_train <- 5000
  num_samples_cal <- 10000
  cohort_size <- 100
  screening_time <- 6
  imputed <- TRUE
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
                 screening_time = screening_time, imputed = imputed)

output_file <- sprintf("results/%s/cens_%s_surv_%s_train%d_cal%d_m%d_time%s_imp%d.txt",
                       setup, cens_model_type, surv_model_type, num_samples_train,
                       num_samples_cal, cohort_size, screening_time, as.integer(imputed))
cat("Output file name:", output_file, "\n")

# ---- Data and models ---------------------------------------------------------
load_raw_data <- function(imputed = TRUE) {
    fname <- if (imputed) "train_test.rds" else "train_test_miss.rds"
    data.path.1 <- file.path("/media/msesia/ExternalSSD/data/fhrd", fname)
    data.path.2 <- file.path("/project2/sesia_1123/data/exploration", fname)
    data.path <- if (file.exists(data.path.1)) data.path.1 else data.path.2
    if (!file.exists(data.path)) stop("Data file not found: ", fname)
    combined <- readRDS(data.path)
    combined <- as_tibble(do.call("rbind", combined)) %>%
        mutate(time = month, status = event) %>%
        select(-month, -event) %>%
        select(time, status, everything())
    # Mapping from generic names (X1, X2, ...) to the original column names
    orig    <- colnames(combined)[-(1:2)]
    x_names <- paste0("X", seq_along(orig))
    colnames(combined) <- c("time", "status", x_names)
    attr(combined, "var_map") <- tibble(X = x_names, original = orig)
    combined
}

# ---- Top-10 variables per survival model (original column names) -----------
# Check each name against print(var_map, n = Inf); the stopifnot below
# fails if any of them is not a column of the data.
TOP_VARS <- list(
  cox = c("ecog_value", "alb_ctns", "chloride", "NLR", "log_heartrate",
          "smoking_code", "age", "sex_code", "tumor_stage", "hgb_ctns"),
  rsf = c("alb_ctns", "ecog_value", "age", "NLR", "chloride",
          "log_heartrate", "hgb_ctns", "oxg", "log_alp_ctns", "log_neut_ctns"),
  grf = c("alb_ctns", "ecog_value", "NLR", "chloride", "log_heartrate",
          "age", "hgb_ctns", "log_alp_ctns", "log_neut_ctns", "log_rbc_ctns"),
  xgb = c("alb_ctns", "ecog_value", "age", "chloride", "NLR",
          "smoking_code", "log_heartrate", "hgb_ctns", "tumor_stage", "log_alp_ctns")
)

# Restrict a data frame to the top variables of one model
restrict_to_top <- function(data, model, var_map, top_vars = TOP_VARS) {
  wanted  <- top_vars[[model]]
  missing <- setdiff(wanted, var_map$original)
  if (length(missing) > 0) stop("Not found in data: ", paste(missing, collapse = ", "))
  x_cols <- var_map$X[match(wanted, var_map$original)]
  data %>% select(time, status, all_of(x_cols))
}

# ---- Use -------------------------------------------------------------------
data.raw <- load_raw_data(imputed = imputed)
var_map  <- attr(data.raw, "var_map")

data.top <- restrict_to_top(data.raw, surv_model_type, var_map)


data.generator <- if (real_data) {
  RealDataGenerator$new(data = data.top)
} else {
  SemiSyntheticDataGenerator$new(data = data.top, surv_model_type = gen_model_type,
                                 cens_model_type = gen_model_type)
}

surv_models <- setNames(lapply(model_names, init_surv_model), model_names)
cens_model  <- CensoringModel$new(model = init_censoring_model(cens_model_type))

# ---- One replication ---------------------------------------------------------

analyze_data <- function(data.train, data.cal, data.test, data.test.oracle) {

    ## Keep only complete cases for model fitting and calibration
    n_tr0 <- nrow(data.train)
    n_cal0 <- nrow(data.cal)
    data.train <- drop_na(data.train)
    data.cal   <- drop_na(data.cal)
    data.test  <- drop_na(data.test)
    cat(sprintf("Complete cases: training %d/%d, calibration %d/%d\n",
                nrow(data.train), n_tr0, nrow(data.cal), n_cal0))
    
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
    bands <- sim_to_long(sim_bands(scores.cal, cohort_size, TRUE), "Bootstrap-sim", model_names, m = cohort_size) %>%
        mutate(t0 = screening_time,
               n_train_cc = nrow(data.train),   # complete cases actually used
               n_cal_cc   = nrow(data.cal)) %>%
        cbind(header)

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
