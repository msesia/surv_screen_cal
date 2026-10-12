library(tidyverse)
library(ggh4x)

# =============================================================================
# Empirical width and coverage of confidence and prediction bands as a
# function of calibration sample size, at the thresholds selected for each
# target yield (Figures "width_yield_ppv_m<m>_t<t0>", one per m and t0).
#
# Input: results/<setup>/surv_<set>_train<n>_cal<n>_m<m>_time<t0>[_batch<b>]_eval.txt,
# one row per (target_yield, repetition) with realized yield / IPCW PPV on a
# test subsample of size m and the bounds of both bands. Coverage and width
# are computed within each batch (random data split) and then averaged over
# batches; error bars show +/- 2 standard errors across batches (n = number of batches).
# =============================================================================

dir.save <- "figures"
dir.create(dir.save, showWarnings = FALSE)

# Figures to produce, one row each: results directory, survival model set
# (surv_model_type), censoring model (cens_model_type), cohort size, horizon,
# and target yields. The file name is built from these:
#   width_yield_ppv_surv-<surv>_cens-<cens>_m<m>_t<t0>[_joint].pdf
PLOTS <- tribble(
  ~setup,                 ~surv,  ~cens,  ~m,   ~t0, ~target_yields,
  "real_selection_1b_cal", "grf",  "grf",  1000, 3,   c(0.02, 0.05),
  "real_selection_1b_cal", "grf",  "grf",  100,  3,   c(0.2, 0.5),
  "real_selection_1b_cal", "grf",  "grf",  1000, 6,   c(0.02, 0.05),
  "real_selection_1b_cal", "grf",  "grf",  1000, 9,   c(0.02, 0.05)
#  "real_selection_1_cal", "grf",  "grf",  1000, 6,   c(0.2, 0.5),
#  "real_selection_1_cal", "grf",  "grf",  1000, 9,   c(0.2, 0.5),
#  "real_selection_1_cal", "grf",  "cox",  1000, 3,   c(0.2, 0.5),
#  "real_selection_2_cal", "mix",  "grf",  1000, 3,   c(0.2, 0.5)
)

# Selection rule used in each setup, appended to the file name when not rule (7)
RULE_TAG <- c(real_selection_1_cal = "", real_selection_1b_cal = "", real_selection_2_cal = "")

figure_name <- function(spec) {
  sprintf("width_yield_ppv_surv-%s_cens-%s_m%d_t%g%s",
          spec$surv, spec$cens, spec$m, spec$t0, RULE_TAG[[spec$setup]])
}

# ---- Loading and aggregation -------------------------------------------------

covers <- function(x, lo, hi) !is.na(x) & x >= lo & x <= hi

load_summary <- function(dir) {
  files <- list.files(dir, pattern = "_eval\\.txt$", full.names = TRUE)
  if (length(files) == 0) {
    warning("No *_eval.txt files found in ", dir)
    return(tibble())
  }
  raw <- map_dfr(files, ~ read_delim(.x, delim = ",", col_types = cols(), show_col_types = FALSE))
  # Independent random splits: `Seed` (one per split, written by the
  # experiment); fall back to `batch` or a single group if absent.
  raw <- raw %>% mutate(seed = if ("Seed" %in% names(raw)) Seed
                               else if ("batch" %in% names(raw)) batch else 1)
  cat(sprintf("%s: %d files, %d independent splits (%s)\n", dir, length(files), n_distinct(raw$seed),
              if ("Seed" %in% names(raw)) "Seed column" else "no Seed column: error bars unavailable"))

  keys <- c("surv_model_type", "cens_model_type", "n_train", "n_cal", "m", "screening_time", "target_yield")

  # Stage 1: within each split (seed), coverage over repetitions and band widths
  per_seed <- raw %>%
    group_by(across(all_of(c(keys, "seed")))) %>%
    summarise(
      model            = first(model),        # selected model (varies for the joint rule)
      lambda           = first(lambda),
      conf_cover_yield = mean(covers(yield_real, yield_conf_lo, yield_conf_hi)),
      conf_cover_ppv   = mean(covers(ppv_real,   ppv_conf_lo,   ppv_conf_hi), na.rm = TRUE),
      pred_cover_yield = mean(covers(yield_real, yield_pred_lo, yield_pred_hi)),
      pred_cover_ppv   = mean(covers(ppv_real,   ppv_pred_lo,   ppv_pred_hi), na.rm = TRUE),
      conf_width_yield = mean(yield_conf_hi - yield_conf_lo),
      conf_width_ppv   = mean(ppv_conf_hi   - ppv_conf_lo),
      pred_width_yield = mean(yield_pred_hi - yield_pred_lo),
      pred_width_ppv   = mean(ppv_pred_hi   - ppv_pred_lo),
      .groups = "drop"
    )

  # Stage 2: mean and standard error across splits
  metric_cols <- c("conf_cover_yield", "conf_cover_ppv", "pred_cover_yield", "pred_cover_ppv",
                   "conf_width_yield", "conf_width_ppv", "pred_width_yield", "pred_width_ppv")
  per_seed %>%
    group_by(across(all_of(keys))) %>%
    summarise(
      n_seeds = n(),
      across(all_of(metric_cols), list(mean = ~ mean(.x, na.rm = TRUE),
                                       se   = ~ if (n() > 1) sd(.x, na.rm = TRUE) / sqrt(n()) else NA_real_)),
      .groups = "drop"
    )
}

# ---- Plot --------------------------------------------------------------------

plot_bands_two_row <- function(summ, set_value, t0_value, m_value, cens_value,
                               target_yields = NULL, font_size = 10) {

  pred_label  <- paste0("Prediction (m = ", m_value, ")")
  band_levels <- c("Confidence", pred_label)
  band_colors <- c("Confidence" = "#2166AC", setNames("#DF8E89", pred_label))
  band_shapes <- c("Confidence" = 16,        setNames(17,        pred_label))

  cols_wide <- c("conf_width_yield", "conf_width_ppv", "pred_width_yield", "pred_width_ppv",
                 "conf_cover_yield", "conf_cover_ppv", "pred_cover_yield", "pred_cover_ppv")

  d <- summ %>%
    filter(surv_model_type == set_value, cens_model_type == cens_value,
           screening_time == t0_value, m == m_value)
  if (!is.null(target_yields)) d <- d %>% filter(target_yield %in% target_yields)

  plot_df <- d %>%
    mutate(facet_label = factor(paste0("Target yield = ", target_yield),
                                levels = paste0("Target yield = ", sort(unique(target_yield))))) %>%
    select(n_cal, facet_label, all_of(paste0(cols_wide, "_mean")), all_of(paste0(cols_wide, "_se"))) %>%
    pivot_longer(cols = -c(n_cal, facet_label),
                 names_to = c("name", ".value"), names_pattern = "^(.*)_(mean|se)$") %>%
    rename(value = mean) %>%
    filter(!is.na(value)) %>%
    mutate(
      metric   = if_else(str_detect(name, "width"), "Width",      "Coverage"),
      quantity = if_else(str_detect(name, "yield"), "Yield",      "PPV"),
      band     = if_else(str_detect(name, "^conf"), "Confidence", pred_label),
      metric   = factor(metric,   levels = c("Width", "Coverage")),
      quantity = factor(quantity, levels = c("Yield", "PPV")),
      band     = factor(band,     levels = band_levels)
    )

  x_range <- range(plot_df$n_cal)

  hline_df <- plot_df %>%
    filter(metric == "Coverage") %>%
    distinct(metric, quantity, facet_label) %>%
    mutate(x = x_range[1], xend = x_range[2], y = 0.95)

  # invisible points to pin the coverage panels to [0, 1]
  dummy_df <- plot_df %>%
    filter(metric == "Coverage") %>%
    distinct(metric, quantity, facet_label) %>%
    cross_join(tibble(value = c(0, 1)))

  ggplot(plot_df, aes(x = n_cal, y = value, color = band, shape = band, group = band)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = value - 2 * se, ymax = value + 2 * se),
                  width = 0.05, alpha = 0.6, na.rm = TRUE) +
    geom_point(data = dummy_df, aes(x = x_range[1], y = value),
               inherit.aes = FALSE, alpha = 0, size = 0) +
    geom_segment(data = hline_df, aes(x = x, xend = xend, y = y, yend = y),
                 linetype = "dashed", linewidth = 0.8, color = "black", inherit.aes = FALSE) +
    facet_nested(metric ~ quantity + facet_label,
                 scales = "free_y", independent = "y",
                 nest_line = element_line(color = "grey40"), labeller = label_value) +
    scale_color_manual(values = band_colors) +
    scale_shape_manual(values = band_shapes) +
    scale_x_continuous(breaks = c(500, 5000, 10000), labels = c("500", "5K", "10K")) +
    labs(x = "Calibration size", y = NULL, color = "Band type", shape = "Band type") +
    theme_bw() +
    theme(legend.position  = "right",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey85"),
          ggh4x.facet.nestline = element_line(color = "grey40"),
          axis.text        = element_text(size = font_size),
          axis.title       = element_text(size = font_size),
          strip.text       = element_text(size = font_size),
          legend.text      = element_text(size = font_size),
          legend.title     = element_text(size = font_size))
}

# ---- Run ---------------------------------------------------------------------

# Load each results directory once
summaries <- map(set_names(unique(PLOTS$setup)),
                 ~ load_summary(file.path("results_hpc", .x)))

for (i in seq_len(nrow(PLOTS))) {
  spec <- PLOTS[i, ]
  summ <- summaries[[spec$setup]]
  if (nrow(summ) == 0) { cat("Skipping", figure_name(spec), ": no results loaded\n"); next }

  avail <- summ %>%
    filter(surv_model_type == spec$surv, cens_model_type == spec$cens,
           m == spec$m, screening_time == spec$t0)
  fname <- figure_name(spec)
  if (nrow(avail) == 0) {
    cat(sprintf("Skipping %s: no results in %s\n", fname, spec$setup))
    next
  }

  p <- plot_bands_two_row(summ, spec$surv, t0_value = spec$t0, m_value = spec$m,
                          cens_value = spec$cens, target_yields = spec$target_yields[[1]])
  ggsave(file.path(dir.save, paste0(fname, ".pdf")), p,
         device = cairo_pdf, width = 10, height = 3)
  cat("Written", fname, "\n")
}
