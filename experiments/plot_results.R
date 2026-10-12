library(tidyverse)
library(ggh4x)
library(kableExtra)
library(showtext)
font_add_google("Noto Sans", "noto")
showtext_auto()

# ---- Loading -----------------------------------------------------------------

load_results <- function(setup) {
  idir <- sprintf("results_hpc/%s", setup)
  ifile.list <- list.files(idir)
  do.call("rbind", lapply(ifile.list, function(ifile) {
    read_delim(sprintf("%s/%s", idir, ifile),
               delim = ",", col_types = cols(), guess_max = 2)
  }))
}

# ---- Constants ---------------------------------------------------------------

# NEW: third method = simultaneous bootstrap WITHOUT the finite-sample
# adjustment (Algorithm 3), used only for the ablation figure.
method.values <- c("Bootstrap-sim", "Bootstrap-sim-baseline", "Bootstrap-point")
method.labels <- c("Simultaneous bootstrap",
                   "Simultaneous bootstrap (no FS adjustment)",
                   "Pointwise bootstrap")
method.slugs  <- c("Bootstrap-sim"          = "sim",
                   "Bootstrap-sim-baseline" = "sim-nofs",
                   "Bootstrap-point"        = "pt")

# Methods for the main per-method figures (unchanged output)
main.method.values <- c("Bootstrap-sim", "Bootstrap-point")

# NEW: methods compared in the ablation figure
ablation.method.values <- c("Bootstrap-sim", "Bootstrap-sim-baseline")
ablation.method.labels <- c("Included", "Excluded")
ablation_width_lambda <- "λ=0.80"   # threshold at which width is shown
ablation_cov_limits   <- c(0.6, 1)  # y-axis range for the coverage row
ablation_colors    <- c("Included" = "#0072B2",
                        "Excluded" = "#D55E00")
ablation_shapes    <- c("Included" = 16,
                        "Excluded" = 17)

# x-axis (calibration size), log scale
n_cal_breaks <- c(10, 100, 1000)

# Model sets (must match the experiment script); used to identify the bands
# that are simultaneously valid over the whole set named by surv_model_type.
MODEL_SETS <- list(
  cox  = "cox",
  grf  = "grf",
  rsf  = "rsf",
  xgb  = "xgb",
  mix  = c("cox", "grf", "rsf", "xgb"),
  mix2 = c("grf", "cox")
)
valid_set_label <- function(set_name) {
  models <- if (set_name %in% names(MODEL_SETS)) MODEL_SETS[[set_name]] else set_name
  paste(models, collapse = "|")
}

band.values <- c("prediction", "confidence")
band.labels <- c("Prediction band", "Confidence band")

# Thresholds shown in the per-threshold figures (the worst case is always
# computed over the full grid)
valid_lambda <- seq(0.5, 1, by = 0.1)
valid_labels <- c("Worst-case", paste0("λ=", formatC(valid_lambda, format = "f", digits = 2)))

lambda_colors <- c(
  "Worst-case" = "#000000",
  "λ=0.50"     = "#E69F00",
  "λ=0.60"     = "#56B4E9",
  "λ=0.70"     = "#009E73",
  "λ=0.80"     = "#F0E442",
  "λ=0.90"     = "#0072B2",
  "λ=1.00"     = "#D55E00"
)

lambda_shapes <- c(
  "Worst-case" = 17,
  "λ=0.50"     = 16,
  "λ=0.60"     = 15,
  "λ=0.70"     = 18,
  "λ=0.80"     = 8,
  "λ=0.90"     = 7,
  "λ=1.00"     = 4
)

DESIGN_COLS <- c("real_data", "gen_model_type",
                 "surv_model_type", "cens_model_type",
                 "n_train", "surv_nph")

DESIGN_ABBR <- c(real_data       = "rd",
                 gen_model_type  = "gen",
                 surv_model_type = "surv",
                 cens_model_type = "cens",
                 n_train         = "ntr",
                 surv_nph        = "nph")

# ---- Pre-processing ----------------------------------------------------------

covers <- function(x, lo, hi) {
  case_when(
    is.na(x) & lo == 0 & hi == 1               ~ TRUE,
    !is.na(x) & !is.na(lo) & x >= lo & x <= hi ~ TRUE,
    !is.na(x) & !is.na(lo)                     ~ FALSE,
    TRUE                                       ~ NA
  )
}

prepare_results <- function(results) {

  # Columns written only by newer runs: survival-tilt strength (0 = fitted
  # generator) and point estimates (fallback: band midpoint)
  if (!"surv_nph"  %in% names(results)) results$surv_nph  <- 0
  if (!"ppv_hat"   %in% names(results)) results$ppv_hat   <- (results$ppv_lo   + results$ppv_hi)   / 2
  if (!"yield_hat" %in% names(results)) results$yield_hat <- (results$yield_lo + results$yield_hi) / 2

  results_cov <- results %>%
    mutate(
      cover_yield = covers(yield, yield_lo, yield_hi),
      cover_ppv   = covers(ppv,   ppv_lo,   ppv_hi),
      width_yield = yield_hi - yield_lo,
      width_ppv   = ppv_hi   - ppv_lo,
      lambda_type = "pointwise",
      band_method_raw = band_method,
      band_method = factor(band_method, method.values, method.labels),
      band_type   = factor(band_type,   band.values,   band.labels)
    )

  # Worst case over all cells (model, lambda) whose model belongs to the set
  # `valid_models` over which the band is simultaneously valid.
  group_cols <- c(
    "Seed", DESIGN_COLS, "n_models",
    "n_cal", "n_test", "screening_time", "alpha", "batch",
    "band_method", "band_method_raw", "band_type", "target_m",
    "valid_models"
  )

  in_valid_set <- function(model, valid_models) {
    mapply(function(m, v) m %in% strsplit(v, "\\|")[[1]], model, valid_models)
  }

  worst_case <- function(d, cols) {
    d %>%
      group_by(across(all_of(cols))) %>%
      summarise(
        N           = n(),
        model       = NA_character_,
        lambda      = NA_real_,
        cover_yield = all(cover_yield),
        cover_ppv   = all(cover_ppv),
        width_yield = max(width_yield, na.rm = TRUE),
        width_ppv   = max(width_ppv,   na.rm = TRUE),
        yield = NA_real_, ppv = NA_real_, yield_hat = NA_real_, ppv_hat = NA_real_,
        yield_lo = NA_real_, yield_hi = NA_real_,
        ppv_lo   = NA_real_, ppv_hi   = NA_real_,
        .groups = "drop"
      ) %>%
      mutate(lambda_type = "worst-case")
  }

  sim_cov_seed <- results_cov %>%
    filter(in_valid_set(model, valid_models)) %>%
    worst_case(group_cols)

  # Multi-model sets: the pointwise baseline has no joint band, so its worst
  # case over ALL models and thresholds of the set is computed here, and
  # labelled with the set's valid_models so it sits next to the joint
  # simultaneous band
  point_cov_set <- results_cov %>%
    filter(band_method_raw == "Bootstrap-point", n_models > 1) %>%
    worst_case(setdiff(group_cols, "valid_models")) %>%
    mutate(valid_models = map_chr(surv_model_type, valid_set_label))

  bind_rows(results_cov, sim_cov_seed, point_cov_set)
}

summarise_results <- function(results_all) {
  results_all %>%
    group_by(across(all_of(c(DESIGN_COLS, "n_models",
                             "n_cal", "n_test", "alpha", "screening_time",
                             "lambda_type", "lambda", "model", "valid_models",
                             "band_method", "band_method_raw",
                             "band_type", "target_m")))) %>%
    summarise(
      R                = n(),
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
      .groups = "drop"
    ) %>%
    mutate(
      lambda_label = ifelse(
        is.na(lambda),
        "Worst-case",
        paste0("λ=", formatC(lambda, format = "f", digits = 2))
      ),
      n_test_label = paste0("Test size = ", n_test)
    )
}

design_slug <- function(row) {
  parts <- vapply(DESIGN_COLS, function(col) {
    sprintf("%s-%s", DESIGN_ABBR[[col]], as.character(row[[col]]))
  }, character(1))
  paste(parts, collapse = "_")
}

# ---- Plotting ---------------------------------------------------------------

plot_coverage_width_one <- function(df, method_label, quantity,
                                    split_test = TRUE) {

  stopifnot(length(method_label) == 1, quantity %in% c("Yield", "PPV"))

  d <- df %>% filter(band_method == method_label)

  if (quantity == "Yield") {
    cov_val <- "yield_coverage";   cov_se <- "se_yield_cov"
    wid_val <- "mean_width_yield"; wid_se <- "se_width_yield"
  } else {
    cov_val <- "ppv_coverage";     cov_se <- "se_ppv_cov"
    wid_val <- "mean_width_ppv";   wid_se <- "se_width_ppv"
  }

  long <- bind_rows(
    d %>% transmute(n_cal, n_test, n_test_label, band_type, lambda_label,
                    metric = "Coverage",
                    value = .data[[cov_val]], se = .data[[cov_se]]),
    d %>% transmute(n_cal, n_test, n_test_label, band_type, lambda_label,
                    metric = "Width",
                    value = .data[[wid_val]], se = .data[[wid_se]])
  ) %>%
    filter(!(metric == "Width" & lambda_label == "Worst-case")) %>%
    filter(is.finite(value)) %>%
    mutate(
      metric = factor(metric, levels = c("Coverage", "Width")),
      ymin = ifelse(metric == "Coverage",
                    pmax(0, value - 1.96 * se),
                    value - 1.96 * se),
      ymax = ifelse(metric == "Coverage",
                    pmin(1, value + 1.96 * se),
                    value + 1.96 * se)
    )

  hline_df <- tibble(
    metric     = factor("Coverage", levels = c("Coverage", "Width")),
    yintercept = 0.95
  )

  facet <- if (split_test) {
    facet_nested(metric ~ band_type + n_test_label,
                 scales = "free_y",
                 nest_line = element_line(linetype = 1),
                 strip = strip_nested(bleed = FALSE))
  } else {
    facet_grid(metric ~ band_type, scales = "free_y")
  }

  long %>%
    filter(lambda_label %in% valid_labels) %>%
    ggplot(aes(x = n_cal, y = value,
               color = lambda_label, shape = lambda_label)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  width = 0.02, alpha = 0.4) +
    geom_hline(data = hline_df, aes(yintercept = yintercept),
               linetype = "dashed", inherit.aes = FALSE) +
    facet +
    labs(x = "Calibration size (n)", y = NULL,
         color = "Screening threshold", shape = "Screening threshold") +
    scale_color_manual(values = lambda_colors) +
    scale_shape_manual(values = lambda_shapes) +
    scale_x_continuous(trans = "log10", breaks = n_cal_breaks) +
    theme_bw()
}

# NEW ------------------------------------------------------------------------
#' Ablation figure: simultaneous bootstrap with vs. without the finite-sample
#' adjustment, prediction bands only. Top row: worst-case coverage. Bottom
#' row: width at a single threshold (ablation_width_lambda). Colour/shape
#' encode whether the adjustment is included.
plot_ablation_coverage_width <- function(df, quantity = "PPV",
                                         split_test = TRUE,
                                         cov_limits = ablation_cov_limits) {

  stopifnot(quantity %in% c("Yield", "PPV"))

  d <- df %>%
    filter(band_method_raw %in% ablation.method.values,
           lambda_label %in% c("Worst-case", ablation_width_lambda),
           band_type == "Prediction band") %>%
    mutate(method_ab = factor(band_method_raw,
                              ablation.method.values,
                              ablation.method.labels))

  if (nrow(d) == 0) return(NULL)

  if (quantity == "Yield") {
    cov_val <- "yield_coverage";   cov_se <- "se_yield_cov"
    wid_val <- "mean_width_yield"; wid_se <- "se_width_yield"
  } else {
    cov_val <- "ppv_coverage";     cov_se <- "se_ppv_cov"
    wid_val <- "mean_width_ppv";   wid_se <- "se_width_ppv"
  }

  width_label <- paste0("Width (", ablation_width_lambda, ")")

  long <- bind_rows(
    d %>% filter(lambda_label == "Worst-case") %>%
      transmute(n_cal, n_test, n_test_label, method_ab,
                metric = "Coverage (WC)",
                value = .data[[cov_val]], se = .data[[cov_se]]),
    d %>% filter(lambda_label == ablation_width_lambda) %>%
      transmute(n_cal, n_test, n_test_label, method_ab,
                metric = width_label,
                value = .data[[wid_val]], se = .data[[wid_se]])
  ) %>%
    filter(is.finite(value)) %>%
    mutate(
      metric = factor(metric, levels = c("Coverage (WC)", width_label)),
      ymin = ifelse(metric == "Coverage (WC)",
                    pmax(0, value - 1.96 * se), value - 1.96 * se),
      ymax = ifelse(metric == "Coverage (WC)",
                    pmin(1, value + 1.96 * se), value + 1.96 * se)
    )

  hline_df <- tibble(
    metric     = factor("Coverage (WC)",
                        levels = c("Coverage (WC)", width_label)),
    yintercept = 0.95
  )

  facet <- if (split_test) {
    facet_grid(metric ~ n_test_label, scales = "free_y")
  } else {
    facet_grid(metric ~ ., scales = "free_y")
  }

  # Per-row y limits: fixed range for coverage, free for width.
  # oob_keep avoids dropping points that fall outside cov_limits.
  y_scales <- list(
    metric == "Coverage (WC)" ~ scale_y_continuous(limits = cov_limits,
                                                   oob = scales::oob_keep),
    metric == width_label     ~ scale_y_continuous()
  )

  ggplot(long, aes(x = n_cal, y = value,
                   color = method_ab, shape = method_ab)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.02, alpha = 0.4) +
    geom_hline(data = hline_df, aes(yintercept = yintercept),
               linetype = "dashed", inherit.aes = FALSE) +
    facet +
    facetted_pos_scales(y = y_scales) +
    labs(x = "Calibration size (n)", y = NULL,
         color = "Finite-sample\nadjustment",
         shape = "Finite-sample\nadjustment") +
    scale_color_manual(values = ablation_colors) +
    scale_shape_manual(values = ablation_shapes) +
    scale_x_continuous(trans = "log10", breaks = n_cal_breaks) +
    theme_bw() +
    theme(legend.position = "right")
}


# NEW ------------------------------------------------------------------------
#' Prediction bands as a function of a design parameter (e.g. the strength
#' of the non-proportional survival tilt, `surv_nph`) at a FIXED calibration
#' size: rows = coverage (worst case and fixed thresholds), width, point
#' estimate, realized value; columns = quantity (yield / PPV) x test size.
plot_vs_param <- function(df, method_label, x_var, x_label,
                          cov_limits = c(0.5, 1)) {
  stopifnot(length(method_label) == 1)
  d <- df %>%
    filter(band_method == method_label, band_type == "Prediction band") %>%
    mutate(x = .data[[x_var]])

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
    filter(!(metric != "Coverage" & lambda_label == "Worst-case"), is.finite(value),
           lambda_label %in% valid_labels) %>%
    mutate(quantity = factor(quantity, levels = c("Yield", "PPV")),
           metric   = factor(metric, levels = metric_levels),
           ymin = ifelse(metric == "Coverage", pmax(0, value - 2 * se), value - 2 * se),
           ymax = ifelse(metric == "Coverage", pmin(1, value + 2 * se), value + 2 * se))

  hline_df <- tibble(metric = factor("Coverage", levels = metric_levels), yintercept = 0.95)

  ggplot(long, aes(x = x, y = value, color = lambda_label, shape = lambda_label)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.02 * diff(range(long$x)), alpha = 0.4) +
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
    labs(x = x_label, y = NULL, color = "Screening threshold", shape = "Screening threshold") +
    scale_color_manual(values = lambda_colors) +
    scale_shape_manual(values = lambda_shapes) +
    scale_x_continuous(breaks = sort(unique(long$x))) +
    theme_bw()
}

# NEW ------------------------------------------------------------------------
#' Multiple candidate models (e.g. setup "v0K", surv_model_type = "mix"):
#' coverage and width vs calibration size. The black line is the worst-case
#' coverage over ALL models and thresholds of the set; the coloured lines show
#' each model separately at a single fixed threshold `fixed_lambda`.
#' For the simultaneous bootstrap these are the joint bands (valid over the
#' whole set); for the pointwise baseline, the per-model intervals.
model_labels <- c(cox = "Cox", grf = "GRF", rsf = "RSF", xgb = "XGB-AFT")
wc_label     <- "Worst case"
model_colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CC79A7")
model_shapes <- c(17, 16, 15, 18, 8)

plot_coverage_width_models <- function(df, method_raw, quantity, set_name,
                                       fixed_lambda = 0.8, split_test = TRUE) {

  stopifnot(quantity %in% c("Yield", "PPV"))
  set_label <- valid_set_label(set_name)

  d <- df %>%
    filter(band_method_raw == method_raw) %>%
    filter(
      # worst case over the whole set
      (lambda_label == "Worst-case" & valid_models == set_label) |
      # per-model cells at the fixed threshold: joint band for the simultaneous
      # method, per-model interval for the pointwise baseline
      (lambda_label != "Worst-case" & abs(lambda - fixed_lambda) < 1e-8 &
         (valid_models == set_label | (method_raw == "Bootstrap-point" & valid_models == model)))
    ) %>%
    mutate(series = ifelse(lambda_label == "Worst-case", wc_label,
                           sprintf("%s (λ = %.2f)", model_labels[model], fixed_lambda)))
  series_levels <- c(wc_label, sprintf("%s (λ = %.2f)", model_labels, fixed_lambda))
  d <- d %>% mutate(series = factor(series, levels = series_levels))
  if (nrow(d) == 0) return(NULL)

  if (quantity == "Yield") {
    cov_val <- "yield_coverage";   cov_se <- "se_yield_cov"
    wid_val <- "mean_width_yield"; wid_se <- "se_width_yield"
  } else {
    cov_val <- "ppv_coverage";     cov_se <- "se_ppv_cov"
    wid_val <- "mean_width_ppv";   wid_se <- "se_width_ppv"
  }

  long <- bind_rows(
    d %>% transmute(n_cal, n_test_label, band_type, series, metric = "Coverage",
                    value = .data[[cov_val]], se = .data[[cov_se]]),
    d %>% transmute(n_cal, n_test_label, band_type, series, metric = "Width",
                    value = .data[[wid_val]], se = .data[[wid_se]])
  ) %>%
    filter(!(metric == "Width" & series == wc_label), is.finite(value)) %>%
    mutate(metric = factor(metric, levels = c("Coverage", "Width")),
           ymin = ifelse(metric == "Coverage", pmax(0, value - 1.96 * se), value - 1.96 * se),
           ymax = ifelse(metric == "Coverage", pmin(1, value + 1.96 * se), value + 1.96 * se))

  hline_df <- tibble(metric = factor("Coverage", levels = c("Coverage", "Width")),
                     yintercept = 0.95)

  facet <- if (split_test) {
    facet_nested(metric ~ band_type + n_test_label, scales = "free_y",
                 nest_line = element_line(linetype = 1), strip = strip_nested(bleed = FALSE))
  } else {
    facet_grid(metric ~ band_type, scales = "free_y")
  }

  ggplot(long, aes(x = n_cal, y = value, color = series, shape = series, group = series)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.02, alpha = 0.4) +
    geom_hline(data = hline_df, aes(yintercept = yintercept),
               linetype = "dashed", inherit.aes = FALSE) +
    facet +
    labs(x = "Calibration size (n)", y = NULL,
         color = "Screening rule", shape = "Screening rule") +
    scale_color_manual(values = setNames(model_colors, series_levels), drop = TRUE) +
    scale_shape_manual(values = setNames(model_shapes, series_levels), drop = TRUE) +
    scale_x_continuous(trans = "log10", breaks = n_cal_breaks) +
    theme_bw()
}

# NEW ------------------------------------------------------------------------
#' Table of worst-case (simultaneous) coverage of the prediction bands over
#' the threshold grid, for yield and PPV, as intervals (mean +/- 2 SE across
#' simulation replicates). Rows: calibration size n. Column groups: method
#' (simultaneous / pointwise) x test size m.
make_coverage_tabular <- function(summ, methods = main.method.values) {

  # Interval mean +/- 2 SE; in red (requires \usepackage{xcolor}) when the
  # upper bound falls below the nominal level
  fmt_int <- function(x, se, nominal = 0.95) {
    lo <- pmax(0, x - 2 * se); hi <- pmin(1, x + 2 * se)   # coverage is a proportion
    txt <- sprintf("[%.3f, %.3f]", lo, hi)
    ifelse(!is.na(hi) & hi < nominal, paste0("\\textcolor{red}{", txt, "}"), txt)
  }

  d <- summ %>%
    filter(band_method_raw %in% methods, band_type == "Prediction band",
           lambda_label == "Worst-case") %>%
    transmute(method = factor(band_method_raw, methods, method.labels[match(methods, method.values)]),
              n_cal, n_test,
              cov_yield = fmt_int(yield_coverage, se_yield_cov),
              cov_ppv   = fmt_int(ppv_coverage,   se_ppv_cov)) %>%
    arrange(method, n_test, n_cal)
  if (nrow(d) == 0) return(NULL)

  groups <- d %>% distinct(method, n_test) %>% arrange(method, n_test)
  key    <- function(x) paste(x$method, x$n_test, sep = "|")

  wide <- d %>%
    mutate(g = match(key(.), key(groups))) %>%
    select(n_cal, g, cov_yield, cov_ppv) %>%
    pivot_wider(names_from = g, values_from = c(cov_yield, cov_ppv), names_glue = "{.value}_{g}") %>%
    arrange(n_cal)
  wide <- wide[, c("n_cal", unlist(map(seq_len(nrow(groups)), ~ paste0(c("cov_yield_", "cov_ppv_"), .x))))]

  col_names <- c("$n$", rep(c("Yield", "PPV"), nrow(groups)))   # repeats are fine in col.names

  span <- function(labels, width) {                # run-length spans of `width` columns each
    r <- rle(labels)
    setNames(c(1, width * r$lengths), c(" ", r$values))   # duplicate names are allowed
  }
  header_m    <- span(sprintf("$m = %d$", groups$n_test), 2)
  header_meth <- span(as.character(groups$method), 2)

  kbl(wide, format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",
      col.names = col_names, align = c("r", rep("c", ncol(wide) - 1))) %>%
    add_header_above(header_m,    escape = FALSE) %>%
    add_header_above(header_meth, escape = FALSE) %>%
    as.character()
}

# ---- Run ---------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("tables",  showWarnings = FALSE, recursive = TRUE)

# Output groups to produce in this run (subset of the tables below). Only the
# results directories needed by these groups are loaded.
RUN <- c("MULTI", "TABLES_MULTI")  # any of "MAIN", "ABLATION", "TABLES", "MULTI", "TABLES_MULTI", "NPH"

# Outputs to produce, one row each. Every row identifies a results directory
# (`setup`) and a design (survival model set, censoring model, generating
# model, n_train, horizon). Only the listed outputs are produced.
#   PLOTS_MAIN     coverage/width vs calibration size, per quantity and method
#   PLOTS_ABLATION finite-sample adjustment ablation, per quantity
#   TABLES_WC      worst-case coverage tables (both methods side by side)
#   PLOTS_MULTI    multiple candidate models: worst case over all models and
#                  thresholds, plus each model at a fixed threshold
#   TABLES_MULTI   worst-case coverage tables for multiple candidate models
#                  (worst case over ALL models and thresholds of the set)
#   PLOTS_NPH      prediction bands vs non-proportionality strength (fixed n_cal)
PLOTS_MAIN <- tribble(
  ~setup, ~surv, ~cens, ~gen,  ~n_train, ~t0, ~quantity, ~method,
  "v0",   "cox", "grf", "grf", 2000,     3,   "PPV",     "Bootstrap-sim",
  "v0",   "cox", "grf", "grf", 2000,     3,   "Yield",   "Bootstrap-sim",
  "v0",   "cox", "grf", "grf", 2000,     3,   "PPV",     "Bootstrap-point",
  "v0",   "cox", "grf", "grf", 2000,     3,   "Yield",   "Bootstrap-point",
  "v0",   "cox", "grf", "grf", 2000,     9,   "PPV",     "Bootstrap-sim",
  "v0",   "cox", "grf", "grf", 2000,     9,   "Yield",   "Bootstrap-sim",
  "v0",   "grf", "grf", "grf", 2000,     3,   "PPV",     "Bootstrap-sim",
  "v0",   "xgb", "grf", "grf", 2000,     3,   "PPV",     "Bootstrap-sim"
)

PLOTS_ABLATION <- tribble(
  ~setup, ~surv, ~cens, ~gen,  ~n_train, ~t0, ~quantity,
  "v0",   "cox", "grf", "grf", 2000,     3,   "PPV",
  "v0",   "cox", "grf", "grf", 2000,     3,   "Yield"
)

TABLES_WC <- tribble(
  ~setup, ~surv, ~cens, ~gen,  ~n_train, ~t0,
  "v0",   "cox", "grf", "grf", 2000,     3,
  "v0",   "cox", "grf", "grf", 2000,     9
)

PLOTS_MULTI <- tribble(
  ~setup, ~surv, ~cens, ~gen,  ~n_train, ~t0, ~quantity, ~method,           ~lambda,
  "v0K",  "mix", "grf", "grf", 2000,     3,   "PPV",     "Bootstrap-sim",   0.8,
  "v0K",  "mix", "grf", "grf", 2000,     3,   "Yield",   "Bootstrap-sim",   0.8,
  "v0K",  "mix", "grf", "grf", 2000,     3,   "PPV",     "Bootstrap-point", 0.8,
  "v0K",  "mix", "grf", "grf", 2000,     3,   "Yield",   "Bootstrap-point", 0.8
)

TABLES_MULTI <- tribble(
  ~setup, ~surv, ~cens, ~gen,  ~n_train, ~t0,
  "v0K",  "mix", "grf", "grf", 2000,     3
)

PLOTS_NPH <- tribble(
  ~setup, ~surv, ~cens, ~gen,  ~n_train, ~n_cal, ~t0, ~method,
  "v0n",  "cox", "grf", "grf", 2000,     1000,   3,   "Bootstrap-sim",
  "v0n",  "grf", "grf", "grf", 2000,     1000,   3,   "Bootstrap-sim"
)

quantity_abbr <- c(Yield = "yi", PPV = "ppv")

# Lazy loading: a results directory is loaded and summarised only the first
# time a figure needs it, then kept in memory for the rest of the run. The
# summary is also cached on disk (summaries_cache/<setup>.rds) and reused as
# long as no results file in the directory is newer than the cache.
.summaries <- new.env()
dir.create("summaries_cache", showWarnings = FALSE)

get_summary <- function(setup) {
  if (!is.null(.summaries[[setup]])) return(.summaries[[setup]])

  cache  <- file.path("summaries_cache", paste0(setup, ".rds"))
  files  <- list.files(file.path("results_hpc", setup), full.names = TRUE)
  fresh  <- file.exists(cache) && length(files) > 0 &&
            file.mtime(cache) >= max(file.mtime(files))
  if (fresh) {
    cat("Loading cached summary for", setup, "\n")
    summ <- readRDS(cache)
  } else {
    cat("Loading and summarising results for", setup, "(", length(files), "files )\n")
    summ <- summarise_results(prepare_results(load_results(setup)))
    saveRDS(summ, cache)
  }
  .summaries[[setup]] <- summ
  summ
}

# Rows of the summary matching a design spec (bands valid over exactly the
# model set named by `surv`); extra filters (e.g. n_cal) applied by the caller
select_design <- function(spec, nph = 0) {
  get_summary(spec$setup) %>%
    filter(surv_model_type == spec$surv, cens_model_type == spec$cens,
           gen_model_type == spec$gen, n_train == spec$n_train,
           screening_time == spec$t0, surv_nph %in% nph,
           valid_models == valid_set_label(spec$surv))
}
design_tag <- function(spec) {
  sprintf("t%.1f_rd-0_gen-%s_surv-%s_cens-%s_ntr-%d", spec$t0, spec$gen, spec$surv, spec$cens, spec$n_train)
}

# Main figures: coverage and width vs calibration size
if ("MAIN" %in% RUN) for (i in seq_len(nrow(PLOTS_MAIN))) {
  spec <- PLOTS_MAIN[i, ]
  df_i <- select_design(spec)
  fname <- sprintf("figures/exp_%s_%s_%s.pdf", quantity_abbr[[spec$quantity]],
                   method.slugs[[spec$method]], design_tag(spec))
  if (nrow(df_i) == 0) { cat("Skipping", fname, ": no results\n"); next }
  p <- plot_coverage_width_one(df_i, method_label = method.labels[match(spec$method, method.values)],
                               quantity = spec$quantity, split_test = TRUE)
  ggsave(fname, p, width = 7, height = 3, units = "in", bg = "white")
  cat("Written", fname, "\n")
}

# Ablation figures: with vs. without the finite-sample adjustment
if ("ABLATION" %in% RUN) for (i in seq_len(nrow(PLOTS_ABLATION))) {
  spec <- PLOTS_ABLATION[i, ]
  df_i <- select_design(spec)
  fname <- sprintf("figures/ablation_fs_%s_%s.pdf", quantity_abbr[[spec$quantity]], design_tag(spec))
  p <- plot_ablation_coverage_width(df_i, quantity = spec$quantity, split_test = TRUE)
  if (is.null(p)) { cat("Skipping", fname, ": no results\n"); next }
  ggsave(fname, p, width = 5.5, height = 3.2, units = "in", bg = "white")
  cat("Written", fname, "\n")
}

# Worst-case coverage tables
if ("TABLES" %in% RUN) for (i in seq_len(nrow(TABLES_WC))) {
  spec <- TABLES_WC[i, ]
  tex <- make_coverage_tabular(select_design(spec))
  fname <- sprintf("tables/worstcase_coverage_%s.tex", design_tag(spec))
  if (is.null(tex)) { cat("Skipping", fname, ": no results\n"); next }
  writeLines(tex, fname)
  cat("Written", fname, "\n")
}

# Non-proportional hazards: prediction bands vs tilt strength at fixed n_cal
if ("NPH" %in% RUN) for (i in seq_len(nrow(PLOTS_NPH))) {
  spec <- PLOTS_NPH[i, ]
  df_i <- get_summary(spec$setup) %>%
    filter(surv_model_type == spec$surv, cens_model_type == spec$cens,
           gen_model_type == spec$gen, n_train == spec$n_train,
           n_cal == spec$n_cal, screening_time == spec$t0,
           valid_models == valid_set_label(spec$surv))
  fname <- sprintf("figures/exp_%s_nph_%s_surv-%s_cens-%s_n%d_t%g.pdf",
                   spec$setup, method.slugs[[spec$method]], spec$surv, spec$cens, spec$n_cal, spec$t0)
  if (nrow(df_i) == 0) { cat("Skipping", fname, ": no results\n"); next }
  p <- plot_vs_param(df_i, method_label = method.labels[match(spec$method, method.values)],
                     x_var = "surv_nph",
                     x_label = "Non-proportionality of survival hazards (tilt strength)")
  ggsave(fname, p, width = 8, height = 5.4, units = "in", bg = "white")
  cat("Written", fname, "\n")
}

# Multiple candidate models: worst case over all models and thresholds, and
# each model at a fixed threshold, vs calibration size
if ("MULTI" %in% RUN) for (i in seq_len(nrow(PLOTS_MULTI))) {
  spec <- PLOTS_MULTI[i, ]
  df_i <- get_summary(spec$setup) %>%
    filter(surv_model_type == spec$surv, cens_model_type == spec$cens,
           gen_model_type == spec$gen, n_train == spec$n_train,
           screening_time == spec$t0, surv_nph == 0)
  fname <- sprintf("figures/exp_%s_%s_%s_lam%.2f_%s.pdf", spec$setup, quantity_abbr[[spec$quantity]],
                   method.slugs[[spec$method]], spec$lambda, design_tag(spec))
  p <- plot_coverage_width_models(df_i, method_raw = spec$method, quantity = spec$quantity,
                                  set_name = spec$surv, fixed_lambda = spec$lambda)
  if (is.null(p)) { cat("Skipping", fname, ": no results\n"); next }
  ggsave(fname, p, width = 7, height = 3, units = "in", bg = "white")
  cat("Written", fname, "\n")
}

# Worst-case coverage tables for multiple candidate models: worst case over
# all models and thresholds of the set, simultaneous (joint bands) and
# pointwise bootstrap side by side
if ("TABLES_MULTI" %in% RUN) for (i in seq_len(nrow(TABLES_MULTI))) {
  spec <- TABLES_MULTI[i, ]
  df_i <- get_summary(spec$setup) %>%
    filter(surv_model_type == spec$surv, cens_model_type == spec$cens,
           gen_model_type == spec$gen, n_train == spec$n_train,
           screening_time == spec$t0, surv_nph == 0,
           valid_models == valid_set_label(spec$surv))
  tex <- make_coverage_tabular(df_i)
  fname <- sprintf("tables/worstcase_coverage_%s_%s.tex", spec$setup, design_tag(spec))
  if (is.null(tex)) { cat("Skipping", fname, ": no results\n"); next }
  writeLines(tex, fname)
  cat("Written", fname, "\n")
}
