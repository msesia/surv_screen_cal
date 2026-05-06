###############################################################################
# Quickstart: Calibration bands for yield and PPV
#
# This script shows how to:
#   1) Load FHRD sample data
#   2) Split into train / validation (cal) / test sizes
#   3) Fit survival + censoring models (GRF by default)
#   4) Compute calibration-time ingredients (scores, O/E labels, IPCW weights)
#   5) Run three inference methods:
#        - bootstrap_pointwise()       : pointwise CI/PI at a fixed lambda
#        - bootstrap_simultaneous()    : uniform (simultaneous) bands over lambda
#        - finite_sample_pointwise()   : finite-sample pointwise bounds (conservative)
#   6) Produce bands over a grid of lambdas and save plots
#
# Where to find the implementations:
#   - Core methods live in ../R/ (loaded via load_methods.R)
#
# Key objects (calibration/validation set only):
#   - scores : numeric risk score used for screening (here: P(T > t0 | X))
#   - S      : selection indicator at a threshold lambda (S = 1{scores >= lambda})
#   - O      : observed event-time indicator (O = 1{T <= C})
#   - E      : observed-by-t0 indicator (E = 1{tildeT <= t0}; set to NA if O=0)
#   - prob_O : estimated P(O=1 | X) or related censoring survival used for IPCW
#
# Notes:
#   - All estimation/inference uses ONLY the validation set (data.cal).
#   - "m" is the *future* cohort size used for prediction intervals; no test data
#     outcomes are required for inference.
###############################################################################

# ----- Libraries -----
library(tidyverse)
library(survival)

# ----- Load methods -----
source("../R/load_methods.R")

######################
## User parameters  ##
######################

# Models used to build the score and IPCW weights
surv_model_type <- "grf"   # survival model used to compute scores
cens_model_type <- "grf"   # censoring model used to compute prob_O

# Sample sizes (drawn from the FHRD cohort)
n_train <- 2000
n_val   <- 2000   # "calibration/validation" size used for inference
n_test  <- 1000   # only used for downstream evaluation (not needed for inference)

# Screening time horizon (months)
t0 <- 3

# Inference controls
alpha  <- 0.05
B_boot <- 500

# Future cohort sizes for prediction intervals
m_small <- 100
m_large <- 1000

# Grid of thresholds for bands
lambda_seq <- seq(0, 1, by = 0.01)

# Reproducibility
seed_split <- 2025
seed_boot  <- 1


##################
## 1) Load data ##
##################

data.full <- read_csv("../data/synthetic_survival_data.csv")

##########################
## 2) Split the dataset  ##
##########################

set.seed(seed_split)
splits <- split_data_n(
  data.full,
  n_train = n_train,
  n_cal   = n_val,
  n_test  = n_test
)

data.train <- splits$data.train
data.val   <- splits$data.cal   # naming: this is the validation/calibration set
data.test  <- splits$data.test


#############################################
## 3) Fit survival model + censoring model  ##
#############################################

surv_model <- init_surv_model(surv_model_type)
surv_model$fit(Surv(time, status) ~ ., data = data.train)

cens_base_model <- init_censoring_model(cens_model_type)
cens_model <- CensoringModel$new(model = cens_base_model)
cens_model$fit(data = data.train)


###########################################################
## 4) Build validation-set ingredients: scores, O/E, prob_O
###########################################################

# IPCW ingredient: estimated censoring survival / observation probability
# (used internally by ET-IPCW; see compute_censoring_survival)
prob_O <- compute_censoring_survival(
  data.val, cens_model,
  data.val$time, data.val$status,
  method = "et"
)

# Score used for screening: P(T > t0 | X) from the survival model
scores <- as.numeric(surv_model$predict(data.val, t0)$predictions)

# Labels for IPCW-style estimation at horizon t0
# O = 1 if event time is observed (uncensored); O=0 if censored before event
O <- as.integer(data.val$status)

# E = 1 if observed time is before t0; undefined if censored (O=0) for event status
E <- as.integer(data.val$time <= t0)
E[O == 0L] <- NA_integer_


#############################################################
## 5) Example: inference at a single threshold lambda
#############################################################

lambda <- 0.5
S <- as.integer(scores >= lambda)

# Pointwise bootstrap: returns confidence + prediction intervals (if m provided)
estim_boot_m100 <- bootstrap_pointwise(
  S = S, E = E, O = O, prob_O = prob_O,
  B = B_boot, m = m_small, alpha = alpha, seed = seed_boot,
  ipcw_method = "hajek"
)

estim_boot_m1000 <- bootstrap_pointwise(
  S = S, E = E, O = O, prob_O = prob_O,
  B = B_boot, m = m_large, alpha = alpha, seed = seed_boot,
  ipcw_method = "hajek"
)

# Finite-sample (pointwise): conservative bounds (confidence + prediction if m provided)
estim_fs_m100 <- finite_sample_pointwise(
  S = S, E = E, O = O, prob_O = prob_O,
  m = m_small, alpha = alpha
)

estim_fs_m1000 <- finite_sample_pointwise(
  S = S, E = E, O = O, prob_O = prob_O,
  m = m_large, alpha = alpha
)

# Print summaries
estim_boot_m100
estim_fs_m100


####################################################################
## 6) Bands over lambda: pointwise bootstrap, simultaneous bootstrap,
##    and finite-sample pointwise
####################################################################

# Helper: build per-lambda rows in the same format as plotting expects
# (make_lambda_row is provided in ../R/)
make_band_rows_pointwise <- function(lambda_seq, scores, E, O, prob_O, B_boot, m, alpha, seed, method_fun, ipcw_method = "ht") {
  purrr::map_dfr(lambda_seq, function(lambda) {
    S <- as.integer(scores >= lambda)
    estim <- method_fun(
      S = S, E = E, O = O, prob_O = prob_O,
      B = B_boot, m = m, alpha = alpha, seed = seed,
      ipcw_method = ipcw_method
    )
    make_lambda_row(lambda, S, estim, E, O, prob_O, m)
  })
}

# --- (a) Pointwise bootstrap bands ---
band_boot_m100 <- make_band_rows_pointwise(
  lambda_seq = lambda_seq, scores = scores, E = E, O = O, prob_O = prob_O,
  B_boot = B_boot, m = m_small, alpha = alpha, seed = seed_boot,
  method_fun = bootstrap_pointwise, ipcw_method = "ht"
)

# --- (b) Simultaneous bootstrap bands (uniform over lambda) ---
# This returns a full data frame already indexed by lambda
band_boot_sim_m100 <- as_tibble(
  bootstrap_simultaneous(
    scores = scores, E = E, O = O, prob_O = prob_O,
    B = B_boot, m = m_small, alpha = alpha, seed = seed_boot,
    ipcw_method = "hajek", lambda_seq = lambda_seq
  )
)

########################################
## 7) Plot and save calibration bands  ##
########################################

p_m100 <- plot_calibration_bands(
  band_list = list(
    "Bootstrap (pointwise)"       = band_boot_m100,
    "Bootstrap (simultaneous)"    = band_boot_sim_m100
  ),
  test_cohort_size = m_small,
  plot_title    = sprintf("Calibration bands (n=%d, \u03B1=%.2f)", n_val, alpha),
  plot_subtitle = sprintf("Selecting patients with P(T>%.1f | X) \u2265 \u03BB", t0)
)

print(p_m100)

dir.create("figures", showWarnings = FALSE)
ggsave(
  filename = sprintf("figures/plot_m%d_t%.1f.png", m_small, t0),
  plot = p_m100, width = 7, height = 6, units = "in", bg = "white"
)


###############################################################################
# Next steps / common modifications
#
# 1) Change t0 (screening horizon):
#      t0 <- 9
#      scores <- as.numeric(surv_model$predict(data.val, t0)$predictions)
#      E <- as.integer(data.val$time <= t0); E[O==0] <- NA
#
# 2) Change m (future cohort size) to see prediction bands widen/narrow:
#      m_small <- 50; m_large <- 2000
#
# 3) If you only need pointwise inference at a single lambda, skip the grid loop
#    and call bootstrap_pointwise() directly.
###############################################################################                                        
