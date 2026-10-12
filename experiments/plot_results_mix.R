library(tidyverse)
library(ggh4x)
library(showtext)
font_add_google("Noto Sans", "noto")
showtext_auto()

# =============================================================================
# Sensitivity to the censoring mechanism (setup "v0t"): coverage, width and
# point estimate of the prediction bands (test sizes side by side) as a
# function of one censoring parameter at a FIXED calibration size. The
# data-generating censoring hazard is multiplied by
# exp(cens_shift + cens_tilt * z), z = standardized X32: `cens_shift` raises
# the overall censoring rate, `cens_tilt` its dependence on covariates. The
# censoring model is always fitted on data from the same law. Same layout as
# the main validation figures; each row of PLOTS selects which parameter is
# on the x-axis, the value of the other one, and the design shown.
# =============================================================================

dir.save   <- "figures"
cov_limits <- c(0.5, 1)              # y-axis range of the coverage row
dir.create(dir.save, showWarnings = FALSE, recursive = TRUE)

# Figures to produce, one row each. `setup` is the results directory
# (results_hpc/<setup>), `x_var` the censoring parameter on the x-axis
# ("cens_shift" or "cens_tilt") and `fixed` the value of the other one; the
# remaining columns identify the design, calibration size, horizon and band
# method. Each figure shows yield and PPV side by side.
PLOTS <- tribble(
  ~setup, ~x_var,       ~fixed, ~surv, ~cens,  ~gen,   ~n_train, ~n_cal, ~t0, ~method,
  "v0t",  "cens_shift", 0,      "grf", "grf",  "grf",  2000,     1000,   3,   "Bootstrap-sim",
  "v0h",  "cens_tilt",  0,      "grf", "grf",  "grf",  2000,     1000,   3,   "Bootstrap-sim"#,
#  "v0t",  "cens_tilt",  1.6,    "cox", "grf",  "grf",  2000,     1000,   3,   "Bootstrap-sim"
)

# For x_var = "cens_shift", the x-axis shows the resulting fraction of
# calibration patients censored before t0 (averaged over replicates) rather
# than the shift itself.
x_labels <- c(cens_shift = expression("Fraction censored before " * t[0]),
              cens_tilt  = "Covariate-dependent censoring (tilt strength)")

# ---- Constants (as in plot_results_with_ablation.R) --------------------------

method.values <- c("Bootstrap-sim", "Bootstrap-sim-baseline", "Bootstrap-point")
method.labels <- c("Simultaneous bootstrap", "Simultaneous bootstrap (no FS adjustment)",
                   "Pointwise bootstrap")
method.slugs  <- c("Bootstrap-sim" = "sim", "Bootstrap-sim-baseline" = "sim-nofs",
                   "Bootstrap-point" = "pt")
main.method.values <- c("Bootstrap-sim", "Bootstrap-point")

band.values <- c("prediction", "confidence")
band.labels <- c("Prediction band", "Confidence band")

valid_lambda <- seq(0.5, 1, by = 0.1)
valid_labels <- c("Worst-case", paste0("λ=", formatC(valid_lambda, format = "f", digits = 2)))
lambda_colors <- c("Worst-case" = "#000000", "λ=0.50" = "#E69F00", "λ=0.60" = "#56B4E9",
                   "λ=0.70" = "#009E73", "λ=0.80" = "#F0E442", "λ=0.90" = "#0072B2", "λ=1.00" = "#D55E00")
lambda_shapes <- c("Worst-case" = 17, "λ=0.50" = 16, "λ=0.60" = 15, "λ=0.70" = 18,
                   "λ=0.80" = 8, "λ=0.90" = 7, "λ=1.00" = 4)

MODEL_SETS <- list(cox = "cox", grf = "grf", rsf = "rsf", xgb = "xgb",
                   mix = c("cox", "grf", "rsf", "xgb"), mix2 = c("grf", "cox"))
valid_set_label <- function(set_name) {
  models <- if (set_name %in% names(MODEL_SETS)) MODEL_SETS[[set_name]] else set_name
  paste(models, collapse = "|")
}

# Columns identifying a design; both censoring parameters are included so the
# aggregation never mixes them (the x-axis one is varied within a figure)
DESIGN_COLS <- c("real_data", "gen_model_type", "surv_model_type", "cens_model_type", "n_train",
                 "cens_shift", "cens_tilt")

# ---- Loading and pre-processing ---------------------------------------------

load_results <- function(setup) {
  idir <- sprintf("results_hpc/%s", setup)
  map_dfr(list.files(idir, full.names = TRUE),
          ~ read_delim(.x, delim = ",", col_types = cols(), guess_max = 2))
}

covers <- function(x, lo, hi) {
  case_when(
    is.na(x) & lo == 0 & hi == 1               ~ TRUE,
    !is.na(x) & !is.na(lo) & x >= lo & x <= hi ~ TRUE,
    !is.na(x) & !is.na(lo)                     ~ FALSE,
    TRUE                                       ~ NA
  )
}

prepare_results <- function(results) {
  # Point estimates: written by newer runs; otherwise fall back to the band
  # midpoint (which is biased where the band is clipped at 0 or 1)
  if (!"ppv_hat"   %in% names(results)) results$ppv_hat   <- (results$ppv_lo   + results$ppv_hi)   / 2
  if (!"yield_hat" %in% names(results)) results$yield_hat <- (results$yield_lo + results$yield_hi) / 2

  results_cov <- results %>%
    mutate(cover_yield = covers(yield, yield_lo, yield_hi),
           cover_ppv   = covers(ppv,   ppv_lo,   ppv_hi),
           width_yield = yield_hi - yield_lo,
           width_ppv   = ppv_hi   - ppv_lo,
           lambda_type = "pointwise",
           band_method_raw = band_method,
           band_method = factor(band_method, method.values, method.labels),
           band_type   = factor(band_type,   band.values,   band.labels))

  group_cols <- c("Seed", DESIGN_COLS, "n_models", "n_cal", "n_test",
                  "screening_time", "alpha", "batch",
                  "band_method", "band_method_raw", "band_type", "target_m", "valid_models")
  in_valid_set <- function(model, valid_models) {
    mapply(function(m, v) m %in% strsplit(v, "\\|")[[1]], model, valid_models)
  }

  # Worst case over all cells in the band's validity set
  sim_cov_seed <- results_cov %>%
    filter(in_valid_set(model, valid_models)) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(N = n(), model = NA_character_, lambda = NA_real_,
              cover_yield = all(cover_yield), cover_ppv = all(cover_ppv),
              width_yield = max(width_yield, na.rm = TRUE),
              width_ppv   = max(width_ppv,   na.rm = TRUE),
              yield = NA_real_, ppv = NA_real_, yield_hat = NA_real_, ppv_hat = NA_real_,
              yield_lo = NA_real_, yield_hi = NA_real_, ppv_lo = NA_real_, ppv_hi = NA_real_,
              .groups = "drop") %>%
    mutate(lambda_type = "worst-case")

  bind_rows(results_cov, sim_cov_seed)
}

summarise_results <- function(results_all) {
  results_all %>%
    group_by(across(all_of(c(DESIGN_COLS, "n_models", "n_cal", "n_test", "alpha",
                             "screening_time", "lambda_type", "lambda", "model", "valid_models",
                             "band_method", "band_method_raw", "band_type", "target_m")))) %>%
    summarise(R                = n(),
              frac_censored_t0 = mean(frac_censored_t0, na.rm = TRUE),
              yield_coverage   = mean(cover_yield, na.rm = TRUE),
              se_yield_cov     = sd(cover_yield, na.rm = TRUE) / sqrt(R),
              ppv_coverage     = mean(cover_ppv,   na.rm = TRUE),
              se_ppv_cov       = sd(cover_ppv,     na.rm = TRUE) / sqrt(R),
              mean_width_yield = mean(width_yield, na.rm = TRUE),
              se_width_yield   = sd(width_yield,   na.rm = TRUE) / sqrt(R),
              mean_width_ppv   = mean(width_ppv,   na.rm = TRUE),
              se_width_ppv     = sd(width_ppv,     na.rm = TRUE) / sqrt(R),
              mean_est_yield   = mean(yield_hat,   na.rm = TRUE),
              se_est_yield     = sd(yield_hat,     na.rm = TRUE) / sqrt(R),
              mean_est_ppv     = mean(ppv_hat,     na.rm = TRUE),
              se_est_ppv       = sd(ppv_hat,       na.rm = TRUE) / sqrt(R),
              mean_true_yield  = mean(yield,       na.rm = TRUE),
              se_true_yield    = sd(yield,         na.rm = TRUE) / sqrt(R),
              mean_true_ppv    = mean(ppv,         na.rm = TRUE),
              se_true_ppv      = sd(ppv,           na.rm = TRUE) / sqrt(R),
              .groups = "drop") %>%
    mutate(lambda_label = ifelse(is.na(lambda), "Worst-case",
                                 paste0("λ=", formatC(lambda, format = "f", digits = 2))),
           n_test_label = paste0("Test size = ", n_test))
}

# ---- Plot: coverage, width and estimate vs censoring parameter ---------------

plot_coverage_width_mix <- function(df, method_label, x_var) {
  stopifnot(length(method_label) == 1)
  d <- df %>% filter(band_method == method_label, band_type == "Prediction band")

  # x positions on a linear scale: for cens_shift, the fraction censored
  # before t0 at each shift level (one value per level, averaged over all
  # replicates at that level); for cens_tilt, the tilt itself
  if (x_var == "cens_shift") {
    x_map <- d %>% group_by(cens_shift) %>%
      summarise(x = weighted.mean(frac_censored_t0, R, na.rm = TRUE), .groups = "drop")
    d <- d %>% left_join(x_map, by = "cens_shift")
  } else {
    d <- d %>% mutate(x = .data[[x_var]])
  }

  metric_levels <- c("Coverage", "Width", "Estimate", "Realized")
  to_long <- function(qty, cov_val, cov_se, wid_val, wid_se, est_val, est_se, true_val, true_se) {
    bind_rows(
      d %>% transmute(x, n_test_label, lambda_label, quantity = qty, metric = "Coverage",
                      value = .data[[cov_val]], se = .data[[cov_se]]),
      d %>% transmute(x, n_test_label, lambda_label, quantity = qty, metric = "Width",
                      value = .data[[wid_val]], se = .data[[wid_se]]),
      d %>% transmute(x, n_test_label, lambda_label, quantity = qty, metric = "Estimate",
                      value = .data[[est_val]], se = .data[[est_se]]),
      d %>% transmute(x, n_test_label, lambda_label, quantity = qty, metric = "Realized",
                      value = .data[[true_val]], se = .data[[true_se]])
    )
  }
  long <- bind_rows(
    to_long("Yield", "yield_coverage", "se_yield_cov", "mean_width_yield", "se_width_yield",
            "mean_est_yield", "se_est_yield", "mean_true_yield", "se_true_yield"),
    to_long("PPV",   "ppv_coverage",   "se_ppv_cov",   "mean_width_ppv",   "se_width_ppv",
            "mean_est_ppv",   "se_est_ppv",   "mean_true_ppv",   "se_true_ppv")
  ) %>%
    filter(!(metric %in% c("Width", "Estimate", "Realized") & lambda_label == "Worst-case"), is.finite(value),
           lambda_label %in% valid_labels) %>%
    mutate(quantity = factor(quantity, levels = c("Yield", "PPV")),
           metric   = factor(metric, levels = metric_levels),
           ymin = ifelse(metric == "Coverage", pmax(0, value - 2 * se), value - 2 * se),
           ymax = ifelse(metric == "Coverage", pmin(1, value + 2 * se), value + 2 * se))

  hline_df <- tibble(metric = factor("Coverage", levels = metric_levels), yintercept = 0.95)

  # x ticks: round breaks spanning the data, with the two extremes always
  # included (and any pretty break too close to an extreme dropped)
  x_rng <- range(long$x)
  x_breaks <- scales::breaks_pretty(n = 4)(x_rng)
  x_breaks <- x_breaks[x_breaks > x_rng[1] & x_breaks < x_rng[2]]
  x_breaks <- x_breaks[abs(x_breaks - x_rng[1]) > 0.08 * diff(x_rng) &
                       abs(x_breaks - x_rng[2]) > 0.08 * diff(x_rng)]
  x_breaks <- sort(c(x_rng[1], x_breaks, x_rng[2]))

  ggplot(long, aes(x = x, y = value, color = lambda_label, shape = lambda_label)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.01, alpha = 0.4) +
    geom_hline(data = hline_df, aes(yintercept = yintercept),
               linetype = "dashed", inherit.aes = FALSE) +
    facet_nested(metric ~ quantity + n_test_label, scales = "free_y",
                 nest_line = element_line(linetype = 1), strip = strip_nested(bleed = FALSE)) +
    facetted_pos_scales(y = list(
      metric == "Coverage" ~ scale_y_continuous(limits = cov_limits, oob = scales::oob_keep),
      metric == "Width"    ~ scale_y_continuous(),
      metric == "Estimate" ~ scale_y_continuous(),
      metric == "Realized" ~ scale_y_continuous()
    )) +
    labs(x = x_labels[[x_var]], y = NULL,
         color = "Screening threshold", shape = "Screening threshold") +
    scale_color_manual(values = lambda_colors) +
    scale_shape_manual(values = lambda_shapes) +
    scale_x_continuous(breaks = x_breaks, labels = function(b) formatC(b, format = "f", digits = 2),
                       expand = expansion(mult = 0.04)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---- Run ---------------------------------------------------------------------

# Load and summarise each results directory once
summaries <- map(set_names(unique(PLOTS$setup)),
                 ~ summarise_results(prepare_results(load_results(.x))))

for (i in seq_len(nrow(PLOTS))) {
  spec      <- PLOTS[i, ]
  other_var <- setdiff(c("cens_shift", "cens_tilt"), spec$x_var)

  df_i <- summaries[[spec$setup]] %>%
    filter(surv_model_type == spec$surv, cens_model_type == spec$cens,
           gen_model_type == spec$gen, n_train == spec$n_train,
           n_cal == spec$n_cal, screening_time == spec$t0,
           .data[[other_var]] == spec$fixed,
           valid_models == valid_set_label(spec$surv))

  fname <- sprintf("%s/exp_%s_%s_%s_surv-%s_cens-%s_n%d_t%g_%s-%g.pdf",
                   dir.save, spec$setup, sub("cens_", "cens", spec$x_var),
                   method.slugs[[spec$method]], spec$surv, spec$cens, spec$n_cal, spec$t0,
                   sub("cens_", "", other_var), spec$fixed)
  if (nrow(df_i) == 0) { cat("Skipping", fname, ": no results\n"); next }
  cat(sprintf("%s: %s levels = %s\n", basename(fname), spec$x_var,
              paste(sort(unique(df_i[[spec$x_var]])), collapse = ", ")))

  m_label <- method.labels[match(spec$method, method.values)]
  p <- plot_coverage_width_mix(df_i, method_label = m_label, x_var = spec$x_var)
  ggsave(fname, p, width = 8, height = 5.4, units = "in", bg = "white")
  cat("Written", fname, "\n")
}
