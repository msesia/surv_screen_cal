library(tidyverse)
library(ggh4x)

# =============================================================================
# Empirical coverage of the prediction bands as a function of the
# realized number of selected patients in the test cohort, at a fixed
# calibration size n. The number selected is binned; within each simulation
# replicate and bin, the band is deemed to cover if it covers at EVERY
# threshold lambda whose realized selection count falls in the bin (worst
# case within the bin). Coverage is then averaged over replicates.
#
# Input: the simulation results used by plot_results_with_ablation.R
# (results_hpc/<setup>/*.txt), with columns Seed, n_cal, n_test, lambda,
# model, num_sel, yield, ppv, band_method, band_type, target_m, *_lo, *_hi,
# valid_models and the design columns.
# =============================================================================

# Figures to produce, one row each. `setup` is the results directory,
# `surv`/`cens`/`gen` identify the design (survival model, censoring model,
# data-generating model), `n_cal` the fixed calibration size and `t0` the
# horizon. Each figure shows yield and PPV side by side, prediction bands only.
PLOTS <- tribble(
  ~setup, ~surv,  ~cens,  ~gen,   ~n_train, ~n_cal, ~t0,
  "v0",   "cox",  "grf",  "grf",  2000,     1000,   3,
  "v0",   "cox",  "grf",  "grf",  2000,     1000,   9
)

sel_breaks <- c(0, 5, 10, 20, 50, 100, 200, 500, 1000)   # bins for number selected (a, b]
dir.save   <- "figures"
dir.create(dir.save, showWarnings = FALSE)

method.values <- c("Bootstrap-sim", "Bootstrap-point")
method.labels <- c("Simultaneous bootstrap", "Pointwise bootstrap")
method.colors <- c("Simultaneous bootstrap" = "#0072B2", "Pointwise bootstrap" = "#D55E00")
method.shapes <- c("Simultaneous bootstrap" = 16,        "Pointwise bootstrap" = 17)
band.values   <- c("prediction", "confidence")
band.labels   <- c("Prediction band", "Confidence band")

DESIGN_COLS <- c("real_data", "gen_model_type", "surv_model_type", "cens_model_type", "n_train")

# ---- Loading -----------------------------------------------------------------

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

bin_label <- function(x, breaks) {
  cut(x, breaks = breaks, include.lowest = FALSE, right = TRUE,
      labels = paste0(head(breaks, -1) + 1, "-", tail(breaks, -1)))
}

# ---- Coverage by number selected -------------------------------------------

summarise_by_num_selected <- function(results, n_cal_fix, breaks) {
  results %>%
    filter(n_cal == n_cal_fix, band_method %in% method.values, num_sel > 0) %>%
    mutate(
      cover_yield = covers(yield, yield_lo, yield_hi),
      cover_ppv   = covers(ppv,   ppv_lo,   ppv_hi),
      sel_bin     = bin_label(num_sel, breaks),
      band_method = factor(band_method, method.values, method.labels),
      band_type   = factor(band_type,   band.values,   band.labels)
    ) %>%
    filter(!is.na(sel_bin)) %>%
    # Stage 1: worst case within (replicate, bin): covered at every lambda in the bin
    group_by(across(all_of(c(DESIGN_COLS, "Seed", "n_cal", "n_test", "screening_time",
                             "band_method", "band_type", "target_m", "model", "sel_bin")))) %>%
    summarise(n_lambda    = n(),
              cover_yield = all(cover_yield, na.rm = TRUE),
              cover_ppv   = all(cover_ppv,   na.rm = TRUE),
              .groups = "drop") %>%
    # Stage 2: average over replicates
    group_by(across(all_of(c(DESIGN_COLS, "n_cal", "n_test", "screening_time",
                             "band_method", "band_type", "target_m", "model", "sel_bin")))) %>%
    summarise(R              = n(),
              yield_coverage = mean(cover_yield), se_yield_cov = sd(cover_yield) / sqrt(R),
              ppv_coverage   = mean(cover_ppv),   se_ppv_cov   = sd(cover_ppv)   / sqrt(R),
              .groups = "drop") %>%
    mutate(n_test_label = paste0("Test size = ", n_test))
}

# ---- Plot --------------------------------------------------------------------

plot_coverage_by_num_selected <- function(summ, min_R = 20) {

  d <- summ %>%
    filter(band_type == "Prediction band", R >= min_R) %>%   # bins seen in too few replicates are dropped
    select(sel_bin, n_test_label, band_method, yield_coverage, se_yield_cov, ppv_coverage, se_ppv_cov) %>%
    pivot_longer(c(yield_coverage, ppv_coverage), names_to = "quantity", values_to = "value") %>%
    mutate(se = if_else(quantity == "yield_coverage", se_yield_cov, se_ppv_cov),
           quantity = factor(if_else(quantity == "yield_coverage", "Yield", "PPV"),
                             levels = c("Yield", "PPV")),
           ymin = pmax(0, value - 2 * se), ymax = pmin(1, value + 2 * se))

  ggplot(d, aes(x = sel_bin, y = value, color = band_method, shape = band_method,
                group = band_method)) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, alpha = 0.5) +
    facet_nested(. ~ quantity + n_test_label, scales = "free_x",
                 nest_line = element_line(linetype = 1),
                 strip = strip_nested(bleed = FALSE)) +
    scale_color_manual(values = method.colors) +
    scale_shape_manual(values = method.shapes) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Number of selected patients in the test cohort", y = "Coverage",
         color = "Method", shape = "Method") +
    theme_bw() +
    theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---- Run ---------------------------------------------------------------------

results <- map(set_names(unique(PLOTS$setup)), load_results)

for (i in seq_len(nrow(PLOTS))) {
  spec <- PLOTS[i, ]

  res <- results[[spec$setup]] %>%
    filter(surv_model_type == spec$surv, cens_model_type == spec$cens,
           gen_model_type == spec$gen, n_train == spec$n_train,
           screening_time == spec$t0)
  summ <- summarise_by_num_selected(res, spec$n_cal, sel_breaks)

  fname <- sprintf("%s/cov_by_numsel_surv-%s_cens-%s_n%d_t%g.pdf",
                   dir.save, spec$surv, spec$cens, spec$n_cal, spec$t0)
  if (nrow(summ) == 0) { cat("Skipping", fname, ": no results\n"); next }

  p <- plot_coverage_by_num_selected(summ)
  ggsave(fname, p, width = 8, height = 2.5, units = "in", bg = "white")
  cat("Written", fname, "\n")
}
