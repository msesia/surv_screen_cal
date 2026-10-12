library(tidyverse)

dir.save <- "results_hpc/real_imputed"

# =============================================================================
# Simultaneous prediction bands for yield and PPV as functions of lambda, on
# the real FHRD data, comparing the analysis on imputed data with the
# complete-case analysis (records with missing covariates removed).
# Facets: quantity (rows) x survival model (columns); colour = data version.
#
# Input: .txt files written by the imputation experiment. The data version
# and nominal calibration size are parsed from the file name
# (..._cal<n>_..._imp<0|1>.txt), so older files without header columns work.
# =============================================================================

# ---- Figures to produce ------------------------------------------------------

PLOTS <- tribble(
  ~n_cal, ~horizon, ~m_val, ~models,
  10000,  6,        100,    c("cox", "grf", "rsf", "xgb"),
  10000,  6,        1000,   c("cox", "grf", "rsf", "xgb"),
  1000,   6,        100,    c("cox", "grf", "rsf", "xgb"),
  1000,   6,        1000,   c("cox", "grf", "rsf", "xgb")
)

# ---- Loading -----------------------------------------------------------------

load_imputation_bands <- function(dir) {
  files <- list.files(dir, pattern = "\\.txt$", full.names = TRUE)
  map_dfr(files, function(f) {
    fn <- basename(f)
    read_delim(f, delim = ",", col_types = cols(), show_col_types = FALSE) %>%
      # parsed from the file name; overrides any header columns of the same name
      mutate(imputed = as.integer(str_match(fn, "_imp(\\d)\\.txt$")[, 2]),
             n_cal   = as.integer(str_match(fn, "_cal(\\d+)_")[, 2]))
  }) %>%
    distinct()
}

bands_all <- load_imputation_bands(dir.save)

# ---- Plot --------------------------------------------------------------------

data_levels  <- c("Imputed", "Complete cases")
data_colors  <- c("Imputed" = "#0072B2", "Complete cases" = "#D55E00")
model_labels <- c(cox = "Cox", grf = "GRF", rsf = "RSF", xgb = "XGB-AFT")

plot_bands_by_imputation <- function(bands, n_cal_val, horizon, m_val,
                                     models      = c("cox", "grf", "rsf", "xgb"),
                                     band_method = "Bootstrap-sim",
                                     xlim        = c(0.2, 1),
                                     font_size   = 10) {

  qty_levels <- c("Yield", "PPV")

  df0 <- bands %>%
    filter(n_cal == n_cal_val, t0 == horizon, model %in% models,
           band_method == !!band_method,
           band_type == "prediction", target_m == m_val) %>%
    mutate(data_version = factor(ifelse(imputed == 1, "Imputed", "Complete cases"),
                                 levels = data_levels),
           model_f = factor(model_labels[model], levels = model_labels[models]))
  if (nrow(df0) == 0) return(NULL)

  to_long <- function(d, lo, hi, est, qty) {
    d %>% transmute(model_f, lambda, data_version, quantity = qty,
                    lo = .data[[lo]], hi = .data[[hi]], estimate = .data[[est]])
  }
  long <- bind_rows(to_long(df0, "yield_lo", "yield_hi", "yield", "Yield"),
                    to_long(df0, "ppv_lo",   "ppv_hi",   "ppv",   "PPV")) %>%
    mutate(quantity = factor(quantity, levels = qty_levels))

  p <- ggplot(long, aes(x = lambda, colour = data_version, fill = data_version)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, colour = NA) +
    geom_line(aes(y = estimate), linewidth = 0.8) +
    facet_grid(quantity ~ model_f, scales = "fixed") +
    coord_cartesian(xlim = xlim) +
    scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1)) +
    scale_colour_manual(values = data_colors, breaks = data_levels) +
    scale_fill_manual(values = data_colors, breaks = data_levels) +
    labs(x = expression("Selection threshold " * lambda), y = NULL,
         colour = sprintf("Prediction band\n(m = %d)", m_val),
         fill   = sprintf("Prediction band\n(m = %d)", m_val)) +
    theme_bw() +
    theme(legend.position  = "right",
          panel.grid.minor = element_blank(),
          axis.text        = element_text(size = font_size),
          axis.title       = element_text(size = font_size),
          strip.text       = element_text(size = font_size),
          legend.text      = element_text(size = font_size),
          legend.title     = element_text(size = font_size))

  # Effective (complete-case) calibration sizes, if recorded, for the caption
  if ("n_cal_cc" %in% names(df0)) {
    df0 %>% distinct(model, data_version, n_cal_cc) %>% arrange(model, data_version) %>% print()
  }

  list(plot = p, data = df0)
}

# ---- Run ---------------------------------------------------------------------

dir.create("figures", showWarnings = FALSE)

for (i in seq_len(nrow(PLOTS))) {
  spec <- PLOTS[i, ]
  out  <- plot_bands_by_imputation(bands_all, n_cal_val = spec$n_cal, horizon = spec$horizon,
                                   m_val = spec$m_val, models = spec$models[[1]])
  fname <- file.path("figures", sprintf("Yield_PPV_vs_lambda_imputation_n%d_t%d_m%d.pdf",
                                        spec$n_cal, spec$horizon, spec$m_val))
  if (is.null(out)) { cat("Skipping", fname, ": no results\n"); next }
  ggsave(fname, plot = out$plot, device = cairo_pdf, width = 9, height = 3.5)
  cat("Written", fname, "\n")
}


library(kableExtra)

# ---- Table: band widths at selected thresholds and overall band area --------

table_lambdas <- c(0.6, 0.7, 0.8, 0.9)   # thresholds at which widths are reported
area_range    <- c(0, 1)                  # lambda range for the area under the width curve

# Rows for one configuration and band type. Confidence bands do not depend on
# m; they are taken from the runs with the given m_val.
select_bands <- function(bands, n_cal_val, horizon, m_val, models,
                         band_type = c("prediction", "confidence"),
                         band_method = "Bootstrap-sim") {
  band_type <- match.arg(band_type)
  bands %>%
    filter(n_cal == n_cal_val, t0 == horizon, model %in% models,
           band_method == !!band_method, band_type == !!band_type,
           target_m == m_val) %>%
    mutate(data_version = factor(ifelse(imputed == 1, "Imputed", "Complete cases"),
                                 levels = data_levels),
           model_f = factor(model_labels[model], levels = model_labels[models]))
}

# Area under the width curve by the trapezoidal rule. Over [0, 1] this equals
# the average band width across thresholds; vacuous PPV bands ([0, 1])
# contribute width 1.
band_area <- function(lambda, width, range = area_range) {
  keep <- lambda >= range[1] & lambda <= range[2] & is.finite(width)
  x <- lambda[keep]; y <- width[keep]
  o <- order(x); x <- x[o]; y <- y[o]
  if (length(x) < 2) return(NA_real_)
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

make_imputation_table <- function(df0, lambdas = table_lambdas, range = area_range) {
  d <- df0 %>%
    mutate(w_yield = yield_hi - yield_lo,
           w_ppv   = ppv_hi   - ppv_lo)

  # Widths at the selected thresholds (tolerant match on the lambda grid)
  at_lambda <- d %>%
    filter(map_lgl(lambda, ~ any(abs(.x - lambdas) < 1e-8))) %>%
    mutate(lam = sprintf("%.1f", lambda)) %>%
    select(model_f, data_version, lam, w_yield, w_ppv) %>%
    pivot_wider(names_from = lam, values_from = c(w_yield, w_ppv), names_glue = "{.value}_{lam}")

  # Area under the width curve
  areas <- d %>%
    group_by(model_f, data_version) %>%
    summarise(area_yield = band_area(lambda, w_yield, range),
              area_ppv   = band_area(lambda, w_ppv,   range), .groups = "drop")

  lam_lab <- sprintf("%.1f", lambdas)
  tab <- at_lambda %>%
    left_join(areas, by = c("model_f", "data_version")) %>%
    arrange(model_f, data_version) %>%
    select(model_f, data_version,
           all_of(paste0("w_yield_", lam_lab)), area_yield,
           all_of(paste0("w_ppv_",   lam_lab)), area_ppv) %>%
    mutate(across(where(is.numeric), ~ sprintf("%.3f", .x)))

  area_lab <- if (all(range == c(0, 1))) "Area" else sprintf("Area [%.1f, %.1f]", range[1], range[2])
  col_names <- c("Model", "Data",
                 paste0("$\\lambda = ", lam_lab, "$"), area_lab,
                 paste0("$\\lambda = ", lam_lab, "$"), area_lab)
  n_l <- length(lambdas)

  kbl(tab, format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",
      col.names = col_names, align = c("l", "l", rep("c", 2 * (n_l + 1)))) %>%
    add_header_above(c(" " = 2, "Yield" = n_l + 1, "PPV" = n_l + 1)) %>%
    collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
    as.character()
}

# ---- Run tables --------------------------------------------------------------

dir.create("tables", showWarnings = FALSE)

# Prediction bands: one table per (n_cal, t0, m)
for (i in seq_len(nrow(PLOTS))) {
  spec <- PLOTS[i, ]
  df_i <- select_bands(bands_all, spec$n_cal, spec$horizon, spec$m_val,
                       spec$models[[1]], band_type = "prediction")
  fname <- file.path("tables", sprintf("imputation_widths_pred_n%d_t%d_m%d.tex",
                                       spec$n_cal, spec$horizon, spec$m_val))
  if (nrow(df_i) == 0) { cat("Skipping", fname, ": no results\n"); next }
  writeLines(make_imputation_table(df_i), fname)
  cat("Written", fname, "\n")
}

# Confidence bands: one table per (n_cal, t0), independent of m
conf_specs <- PLOTS %>% group_by(n_cal, horizon) %>% slice(1) %>% ungroup()
for (i in seq_len(nrow(conf_specs))) {
  spec <- conf_specs[i, ]
  df_i <- select_bands(bands_all, spec$n_cal, spec$horizon, spec$m_val,
                       spec$models[[1]], band_type = "confidence")
  fname <- file.path("tables", sprintf("imputation_widths_conf_n%d_t%d.tex",
                                       spec$n_cal, spec$horizon))
  if (nrow(df_i) == 0) { cat("Skipping", fname, ": no results\n"); next }
  writeLines(make_imputation_table(df_i), fname)
  cat("Written", fname, "\n")
}


library(kableExtra)

# ---- Table: band widths at selected thresholds and overall band area --------

table_lambdas <- c(0.6, 0.7, 0.8, 0.9)   # thresholds at which widths are reported
area_range    <- c(0, 1)                  # lambda range for the area under the width curve

# Rows for one configuration and band type. Prediction bands are selected by
# their target cohort size m. Confidence bands do not depend on m (their
# target_m may be NA); every run writes the same band, so they are taken
# from the runs with the given m_val when that is recorded, and otherwise
# deduplicated across runs.
select_bands <- function(bands, n_cal_val, horizon, m_val, models,
                         band_type = c("prediction", "confidence"),
                         band_method = "Bootstrap-sim") {
  band_type <- match.arg(band_type)
  d <- bands %>%
    filter(n_cal == n_cal_val, t0 == horizon, model %in% models,
           band_method == !!band_method, band_type == !!band_type)

  if (band_type == "prediction") {
    d <- d %>% filter(target_m == m_val)
  } else {
    if ("m" %in% names(d) && any(d$m == m_val, na.rm = TRUE)) d <- d %>% filter(m == m_val)
    d <- d %>% distinct(model, imputed, lambda, .keep_all = TRUE)
  }

  d %>%
    mutate(data_version = factor(ifelse(imputed == 1, "Imputed", "Complete cases"),
                                 levels = data_levels),
           model_f = factor(model_labels[model], levels = model_labels[models]))
}

# Area under the width curve by the trapezoidal rule. Over [0, 1] this equals
# the average band width across thresholds; vacuous PPV bands ([0, 1])
# contribute width 1.
band_area <- function(lambda, width, range = area_range) {
  keep <- lambda >= range[1] & lambda <= range[2] & is.finite(width)
  x <- lambda[keep]; y <- width[keep]
  o <- order(x); x <- x[o]; y <- y[o]
  if (length(x) < 2) return(NA_real_)
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

make_imputation_table <- function(df0, lambdas = table_lambdas, range = area_range) {
  d <- df0 %>%
    mutate(w_yield = yield_hi - yield_lo,
           w_ppv   = ppv_hi   - ppv_lo,
           # effective calibration size: complete cases if recorded, else nominal
           n_eff   = if ("n_cal_cc" %in% names(df0)) n_cal_cc else n_cal)

  # Effective calibration size per model and data version
  n_eff_tab <- d %>% distinct(model_f, data_version, n_eff)

  # Widths at the selected thresholds (tolerant match on the lambda grid)
  at_lambda <- d %>%
    filter(map_lgl(lambda, ~ any(abs(.x - lambdas) < 1e-8))) %>%
    mutate(lam = sprintf("%.1f", lambda)) %>%
    select(model_f, data_version, lam, w_yield, w_ppv) %>%
    pivot_wider(names_from = lam, values_from = c(w_yield, w_ppv), names_glue = "{.value}_{lam}")

  # Area under the width curve
  areas <- d %>%
    group_by(model_f, data_version) %>%
    summarise(area_yield = band_area(lambda, w_yield, range),
              area_ppv   = band_area(lambda, w_ppv,   range), .groups = "drop")

  lam_lab  <- sprintf("%.1f", lambdas)
  val_cols <- c(paste0("w_yield_", lam_lab), "area_yield",
                paste0("w_ppv_",   lam_lab), "area_ppv")

  # Format to 3 decimals and bold, within each model, the smaller of the two
  # values (imputed vs complete cases); no bold if they are equal as printed
  fmt_bold_min <- function(x) {
    r   <- round(x, 3)
    txt <- ifelse(is.na(x), "--", sprintf("%.3f", x))
    is_min <- !is.na(r) & r == min(r, na.rm = TRUE) & n_distinct(r[!is.na(r)]) > 1
    ifelse(is_min, paste0("\\textbf{", txt, "}"), txt)
  }

  tab <- at_lambda %>%
    left_join(areas,     by = c("model_f", "data_version")) %>%
    left_join(n_eff_tab, by = c("model_f", "data_version")) %>%
    arrange(model_f, data_version) %>%
    select(model_f, data_version, n_eff, all_of(val_cols)) %>%
    group_by(model_f) %>%
    mutate(across(all_of(val_cols), fmt_bold_min)) %>%
    ungroup() %>%
    mutate(n_eff = format(n_eff, big.mark = ","))

  area_lab <- if (all(range == c(0, 1))) "Area" else sprintf("Area [%.1f, %.1f]", range[1], range[2])
  col_names <- c("Model", "Data", "$n$",
                 paste0("$\\lambda = ", lam_lab, "$"), area_lab,
                 paste0("$\\lambda = ", lam_lab, "$"), area_lab)
  n_l <- length(lambdas)

  kbl(tab, format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",
      col.names = col_names, align = c("l", "l", "r", rep("c", 2 * (n_l + 1)))) %>%
    add_header_above(c(" " = 3, "Yield" = n_l + 1, "PPV" = n_l + 1)) %>%
    collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
    as.character()
}

# ---- Run tables --------------------------------------------------------------

dir.create("tables", showWarnings = FALSE)

# Prediction bands: one table per (n_cal, t0, m)
for (i in seq_len(nrow(PLOTS))) {
  spec <- PLOTS[i, ]
  df_i <- select_bands(bands_all, spec$n_cal, spec$horizon, spec$m_val,
                       spec$models[[1]], band_type = "prediction")
  fname <- file.path("tables", sprintf("imputation_widths_pred_n%d_t%d_m%d.tex",
                                       spec$n_cal, spec$horizon, spec$m_val))
  if (nrow(df_i) == 0) { cat("Skipping", fname, ": no results\n"); next }
  writeLines(make_imputation_table(df_i), fname)
  cat("Written", fname, "\n")
}

# Confidence bands: one table per (n_cal, t0), independent of m
conf_specs <- PLOTS %>% group_by(n_cal, horizon) %>% slice(1) %>% ungroup()
for (i in seq_len(nrow(conf_specs))) {
  spec <- conf_specs[i, ]
  df_i <- select_bands(bands_all, spec$n_cal, spec$horizon, spec$m_val,
                       spec$models[[1]], band_type = "confidence")
  fname <- file.path("tables", sprintf("imputation_widths_conf_n%d_t%d.tex",
                                       spec$n_cal, spec$horizon))
  if (nrow(df_i) == 0) { cat("Skipping", fname, ": no results\n"); next }
  writeLines(make_imputation_table(df_i), fname)
  cat("Written", fname, "\n")
}
