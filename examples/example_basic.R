###############################################################################
# Quickstart: calibration bands for yield and PPV of survival screening rules
#
# A screening rule selects patients whose predicted probability of surviving
# past a horizon t0 is at least lambda:  A(x; k, lambda) = 1{ z_k(x) >= lambda },
# where z_k is the score of survival model k. This script shows how to:
#   1) Load data and split it into training / calibration / test sets
#   2) Fit one or more survival models and a censoring model
#   3) Compute the calibration-set ingredients (scores, O/E labels, IPCW inputs)
#   4) Construct simultaneous confidence and prediction bands for yield and
#      PPV over a grid of thresholds, either for a single model or jointly
#      over several candidate models (valid under post hoc selection of both
#      the model and the threshold)
#   5) Compare with the pointwise bootstrap baseline
#   6) Select a model and threshold with an explicit rule
#   7) Plot the bands
#
# Core methods live in ../R/ (loaded via load_methods.R).
#
# Key objects (calibration set only):
#   - scores : list of numeric scores, one per model (z_k = P(T > t0 | X))
#   - O      : event-observed indicator, O = 1{T <= C}
#   - E      : observed-by-t0 indicator, E = 1{tildeT <= t0}; NA if O = 0
#   - prob_O : estimated censoring survival used for IPCW
#
# Notes:
#   - All inference uses ONLY the calibration set (data.cal); the test set is
#     not needed and is kept here only for optional downstream evaluation.
#   - "m" is the size of the future cohort targeted by the prediction bands.
###############################################################################

# ----- Libraries -----
library(tidyverse)
library(survival)

# ----- Load methods -----
source("../R/load_methods.R")

######################
## User parameters  ##
######################

# Candidate survival models (any subset of "cox", "grf", "rsf", "xgb").
# With a single model, the bands are simultaneous over thresholds only;
# with several, they are simultaneous over all models and thresholds.
surv_models_list <- c("cox", "grf", "rsf", "xgb")
cens_model_type  <- "grf"   # censoring model used for IPCW

# Sample sizes
n_train <- 2000
n_cal   <- 2000
n_test  <- 1000             # only for optional downstream evaluation

# Screening horizon (months)
t0 <- 3

# Inference controls
alpha  <- 0.05
B_boot <- 500
n_min  <- 10                # minimum number of selected calibration patients for a PPV band
pi_min <- 0.9               # minimum fraction of non-empty simulated cohorts for a PPV prediction band

# Future cohort size for the prediction bands
m <- 100

# Grid of thresholds
lambda_seq <- seq(0, 1, by = 0.01)

# Reproducibility
seed_split <- 2025
seed_boot  <- 1


##################
## 1) Load data ##
##################

data.full <- read_csv("../data/synthetic_survival_data.csv", show_col_types = FALSE) %>%
  mutate(across(where(is.character), as.factor))

set.seed(seed_split)
splits <- split_data_n(data.full, n_train = n_train, n_cal = n_cal, n_test = n_test)
data.train <- splits$data.train
data.cal   <- splits$data.cal
data.test  <- splits$data.test


##############################################
## 2) Fit survival models + censoring model ##
##############################################

surv_models <- setNames(lapply(surv_models_list, init_surv_model), surv_models_list)
for (k in surv_models_list) surv_models[[k]]$fit(Surv(time, status) ~ ., data = data.train)

cens_model <- CensoringModel$new(model = init_censoring_model(cens_model_type))
cens_model$fit(data = data.train)


##################################################
## 3) Calibration-set ingredients               ##
##################################################

# One score vector per model: P(T > t0 | X)
scores <- lapply(surv_models, function(mod) as.numeric(mod$predict(data.cal, t0)$predictions))

# IPCW input: estimated censoring survival at the observed times
prob_O <- compute_censoring_survival(data.cal, cens_model, data.cal$time, data.cal$status, method = "et")

# Labels
O <- as.integer(data.cal$status)
E <- as.integer(data.cal$time <= t0)
E[O == 0L] <- NA_integer_


##################################################
## 4) Simultaneous bands (our method)           ##
##################################################

# Runs the simultaneous bootstrap (with the finite-sample adjustment) for a
# list of score vectors, and returns the bands in long format together with
# the point estimates. `valid_models` records the set of models over which
# the bands are jointly valid.
simultaneous_bands <- function(scores, m, fs_smoothing = TRUE) {
  models <- names(scores)
  estim  <- bootstrap_simultaneous(
    scores = if (length(scores) == 1) scores[[1]] else scores,
    E = E, O = O, prob_O = prob_O,
    B = B_boot, m = m, alpha = alpha, seed = seed_boot,
    fs_smoothing = fs_smoothing, ipcw_method = "ht",
    lambda_seq = lambda_seq, n_min = n_min, pi_min = pi_min
  )
  map_dfr(seq_along(models), function(j) {
    estim_j <- if (length(models) == 1) estim else unpack_model(estim, j)
    est <- tibble(lambda    = estim_j$estimate$lambda,
                  yield_hat = estim_j$estimate$pS,
                  ppv_hat   = 1 - estim_j$estimate$pE)
    convert_band_to_long_targets(estim = estim_j, band_method = "Simultaneous", m = m) %>%
      left_join(est, by = "lambda") %>%
      mutate(model = models[j], valid_models = paste(models, collapse = "|"))
  })
}

# (a) Bands jointly valid over all candidate models and thresholds
bands_joint <- simultaneous_bands(scores, m = m)

# (b) For comparison: bands valid over thresholds only, one model at a time
bands_single <- map_dfr(surv_models_list, function(k) simultaneous_bands(scores[k], m = m))


##################################################
## 5) Pointwise bootstrap baseline              ##
##################################################

# Intervals calibrated threshold by threshold: valid at a FIXED, pre-specified
# lambda, but not under data-driven selection of the threshold or the model.
bands_pointwise <- map_dfr(surv_models_list, function(k) {
  band_df <- map_dfr(lambda_seq, function(lambda) {
    S <- as.integer(scores[[k]] >= lambda)
    estim <- bootstrap_pointwise(S = S, E = E, O = O, prob_O = prob_O,
                                 B = B_boot, m = m, alpha = alpha, seed = seed_boot,
                                 ipcw_method = "ht", n_min = n_min, pi_min = pi_min)
    make_lambda_row(lambda, S, estim, E, O, prob_O, m)
  })
  convert_pointwise_to_long_targets(band_df = band_df, band_method = "Pointwise", m = m) %>%
    mutate(model = k, valid_models = k)
})


##################################################
## 6) Selecting a model and threshold           ##
##################################################

# Because the joint bands are simultaneous over models and thresholds, any
# data-driven choice keeps a valid uncertainty statement. Example rule:
# among all (model, lambda) whose predictive lower bound on yield is at least
# mu0, pick the one maximizing the predictive lower bound on PPV.
mu0 <- 0.2

selection <- bands_joint %>%
  filter(band_type == "prediction", yield_lo >= mu0) %>%
  slice_max(ppv_lo, n = 1, with_ties = FALSE) %>%
  select(model, lambda, yield_hat, yield_lo, yield_hi, ppv_hat, ppv_lo, ppv_hi)

print(selection)


########################################
## 7) Plot the bands                  ##
########################################

model_labels <- c(cox = "Cox", grf = "GRF", rsf = "RSF", xgb = "XGB-AFT")

plot_df <- bind_rows(
  bands_joint  %>% mutate(validity = "All models"),
  bands_single %>% mutate(validity = "Single model")
) %>%
  filter(band_type == "prediction") %>%
  transmute(model = factor(model_labels[model], levels = model_labels[surv_models_list]),
            lambda, validity,
            Yield_lo = yield_lo, Yield_hi = yield_hi, Yield_est = yield_hat,
            PPV_lo = ppv_lo, PPV_hi = ppv_hi, PPV_est = ppv_hat) %>%
  pivot_longer(-c(model, lambda, validity),
               names_to = c("quantity", ".value"), names_sep = "_") %>%
  mutate(quantity = factor(quantity, levels = c("Yield", "PPV")),
         validity = factor(validity, levels = c("All models", "Single model")))

p <- ggplot(plot_df, aes(x = lambda)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = validity)) +
  geom_line(data = filter(plot_df, validity == "Single model"),
            aes(y = est), linewidth = 0.8) +
  facet_grid(quantity ~ model) +
  scale_fill_manual(values = c("All models" = "#DF8E89", "Single model" = "#FDDBC7")) +
  coord_cartesian(xlim = c(0.2, 1)) +
  labs(x = expression("Screening threshold " * lambda), y = NULL,
       fill = sprintf("Prediction band (m = %d)\nvalid over", m),
       title = sprintf("Simultaneous 95%% prediction bands (n = %d, t0 = %g months)", n_cal, t0)) +
  theme_bw()

print(p)

dir.create("figures", showWarnings = FALSE)
ggsave(sprintf("figures/quickstart_bands_m%d_t%g.png", m, t0),
       plot = p, width = 9, height = 4.5, units = "in", bg = "white")
