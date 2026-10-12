options(width=250)

library(tidyverse)

parse_input <- TRUE

source("../R/load_methods.R")

# =============================================================================
# Same as the selection experiment, but each survival model is first fitted on
# all variables to rank them by model-specific importance, then refitted on
# its top-`num_vars` variables; bands, threshold selection and test-set
# evaluation use the reduced models. The importance ranking is saved.
# =============================================================================

# ---- Input parameters --------------------------------------------------------

if (parse_input) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 7) stop("Expected 7 arguments.")
  setup             <- args[1]
  cens_model_type   <- args[2]
  surv_model_type   <- args[3]              # a model-set name, see MODEL_SETS
  num_samples_train <- as.integer(args[4])
  num_samples_cal   <- as.integer(args[5])
  cohort_size       <- as.integer(args[6])
  screening_time    <- as.numeric(args[7])
} else {
  setup <- "real_selection_1_top5"
  cens_model_type <- "grf"
  surv_model_type <- "grf"
  num_samples_train <- 5000
  num_samples_cal <- 10000
  cohort_size <- 100
  screening_time <- 6
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

real_data <- 1
num_samples_test <- cohort_size
alpha            <- 0.05
B_boot           <- 500
lambda_seq       <- seq(0, 1, by = 0.01)
n_min            <- 10       # paper: n_min
pi_min           <- 0.9      # paper: pi_min
num_vars         <- 5        # number of top-ranked variables used to refit each model

target_yield_list <- c(0.02,0.05,0.2,0.3,0.4,0.5)

# ---- Output ------------------------------------------------------------------

header <- tibble(real_data = real_data, surv_model_type = surv_model_type,
                 cens_model_type = cens_model_type, n_models = K,
                 n_train = num_samples_train, n_cal = num_samples_cal,
                 m = cohort_size, n_test = num_samples_test, alpha = alpha,
                 screening_time = screening_time, num_vars = num_vars)

output_file <- sprintf("results/%s/cens_%s_surv_%s_train%d_cal%d_m%d_time%s_vars%d.txt",
                       setup, cens_model_type, surv_model_type, num_samples_train, num_samples_cal,
                       cohort_size, screening_time, num_vars)
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

    ## ---- Fit on all variables, rank them, refit on the top-num_vars ----------
    full_models <- setNames(lapply(model_names, init_surv_model), model_names)
    for (k in model_names) full_models[[k]]$fit(Surv(time, status) ~ ., data = data.train)

    ranking <- map_dfr(model_names, function(k) {
        imp <- sort(full_models[[k]]$importance(), decreasing = TRUE)
        tibble(model = k, rank = seq_along(imp), code = names(imp), importance = as.numeric(imp)) %>%
            left_join(var_names, by = "code")
    })
    write.csv(bind_cols(header, ranking),
              str_replace(output_file, "\\.txt$", "_importance.txt"), row.names = FALSE)

    surv_models <- setNames(lapply(model_names, function(k) {
        top_k <- ranking %>% filter(model == k, rank <= num_vars) %>% pull(code)
        mod <- init_surv_model(k, use_covariates = top_k)
        mod$fit(Surv(time, status) ~ ., data = data.train)
        mod
    }), model_names)

    ## One score vector per (reduced) model
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

    bands <- sim_to_long(sim_bands(scores.cal, cohort_size, TRUE), "Bootstrap-sim", model_names, m = cohort_size) %>%
        mutate(t0 = screening_time)

    valid_set <- paste(model_names, collapse = "|")
    if (grepl("^real_selection_1", setup)) {
        selection <- function(bands, target_yield) {                 # rule (7)
            bands %>%
                filter(band_type=="prediction", yield_lo >= target_yield, target_m == cohort_size, valid_models == valid_set) %>%
                slice_max(lambda, n = 1, with_ties = FALSE) %>%
                mutate(target_yield = target_yield)
        }
    } else {
        selection <- function(bands, target_yield) {                 # joint (model, lambda) rule
            bands %>%
                filter(band_type=="prediction", yield_lo > target_yield, target_m == cohort_size, valid_models == valid_set) %>%
                arrange(desc(ppv_lo), desc(lambda)) %>%
                slice_head(n = 1) %>%
                mutate(target_yield = target_yield)
        }
    }

    selections <- map_dfr(target_yield_list, ~ selection(bands, .x))

    write.csv(bind_cols(header, selections), output_file, row.names = FALSE)
    cat(sprintf("Written selections on %s.\n", output_file))

    ## ---- Test-set evaluation of the selected rules ----------------------------
    n_rep <- 100

    prob_O.test <- compute_censoring_survival(data.test, cens_model, data.test$time, data.test$status, method = "et")
    O.test <- as.integer(data.test$status)
    E.test <- as.integer(data.test$time <= screening_time)
    E.test[O.test == 0L] <- NA_integer_
    w.test <- ifelse(O.test == 1L, 1 / pmax(prob_O.test, 1e-6), 0)   # IPCW weights (HT)

    ipcw_ppv <- function(selected) {
        if (length(selected) == 0) return(NA_real_)
        r <- mean(ifelse(is.na(E.test[selected]), 0, E.test[selected]) * w.test[selected])
        1 - min(1, max(0, r))
    }

    evaluate_rule <- function(row) {
        sel_test <- scores.test[[row$model]] >= row$lambda
        cb <- bands %>%
            filter(model == row$model, lambda == row$lambda, band_type == "confidence")

        map_dfr(seq_len(n_rep), function(r) {
            set.seed(1000 + r)
            idx <- sample.int(num_samples_test, cohort_size, replace = FALSE)
            s   <- idx[sel_test[idx]]
            tibble(
                target_yield  = row$target_yield,
                model         = row$model,
                lambda        = row$lambda,
                rep           = r,
                yield_real    = length(s) / cohort_size,
                ppv_real      = ipcw_ppv(s),
                yield_pred_lo = row$yield_lo, yield_pred_hi = row$yield_hi,
                ppv_pred_lo   = row$ppv_lo,   ppv_pred_hi   = row$ppv_hi,
                yield_conf_lo = cb$yield_lo,  yield_conf_hi = cb$yield_hi,
                ppv_conf_lo   = cb$ppv_lo,    ppv_conf_hi   = cb$ppv_hi
            )
        })
    }

    evaluation <- map_dfr(seq_len(nrow(selections)), function(i) evaluate_rule(selections[i, ]))

    output_file_eval <- str_replace(output_file, "\\.txt$", "_eval.txt")
    write.csv(bind_cols(header, evaluation), output_file_eval, row.names = FALSE)
    cat(sprintf("Written evaluation on %s.\n", output_file_eval))
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
run_experiment(1)
