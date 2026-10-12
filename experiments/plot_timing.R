library(tidyverse)
library(showtext)
font_add_google("Noto Sans", "noto")
showtext_auto()

# =============================================================================
# Wall-clock time of the simultaneous bootstrap (B = 500, 101 thresholds) as a
# function of the calibration size n, for cohort sizes m (facets), with and
# without the finite-sample adjustment of Algorithm 3 (colours). Points are
# means over repetitions, error bars +/- 2 SE.
#
# Input: results/timing/surv_<set>_cal<n>_m<m>_batch<b>.txt, one row per
# repetition with `time_sim_nofs` and `time_sim_fs` (seconds).
# =============================================================================

dir.save <- "figures"
dir.create(dir.save, showWarnings = FALSE)

PLOTS <- tribble(
  ~setup,   ~surv,
  "timing", "grf"
)

fs.labels <- c(time_sim_fs = "Included", time_sim_nofs = "Excluded")
fs.colors <- c("Included" = "#0072B2", "Excluded" = "#D55E00")
fs.shapes <- c("Included" = 16,        "Excluded" = 17)

# ---- Loading and aggregation -------------------------------------------------

load_timing <- function(setup) {
  list.files(file.path("results_hpc", setup), pattern = "\\.txt$", full.names = TRUE) %>%
    map_dfr(~ read_delim(.x, delim = ",", col_types = cols(), show_col_types = FALSE))
}

summarise_timing <- function(raw) {
  raw %>%
    pivot_longer(c(time_sim_fs, time_sim_nofs), names_to = "fs", values_to = "time") %>%
    group_by(surv_model_type, n_cal, m, B, n_lambda, fs) %>%
    summarise(R = n(), mean_time = mean(time), se_time = sd(time) / sqrt(R), .groups = "drop") %>%
    mutate(fs      = factor(fs.labels[fs], levels = fs.labels),
           m_label = paste0("Cohort size m = ", m))
}

# ---- Plot --------------------------------------------------------------------

plot_timing <- function(summ, font_size = 10) {
  ggplot(summ, aes(x = n_cal, y = mean_time, color = fs, shape = fs, group = fs)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    #geom_errorbar(aes(ymin = mean_time - 2 * se_time, ymax = mean_time + 2 * se_time),
    #              width = 0.05, alpha = 0.5) +
    facet_grid(. ~ m_label) +
    scale_x_continuous(trans = "log10", breaks = c(100, 1000, 10000),
                       labels = c("100", "1K", "10K")) +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(values = fs.colors) +
    scale_shape_manual(values = fs.shapes) +
    labs(x = "Calibration size (n)", y = "Computation time (s)",
         color = "Finite-sample\nadjustment", shape = "Finite-sample\nadjustment") +
    theme_bw() +
    theme(legend.position  = "right",
          panel.grid.minor = element_blank(),
          axis.text        = element_text(size = font_size),
          axis.title       = element_text(size = font_size),
          strip.text       = element_text(size = font_size),
          legend.text      = element_text(size = font_size),
          legend.title     = element_text(size = font_size))
}

# ---- Run ---------------------------------------------------------------------

for (i in seq_len(nrow(PLOTS))) {
  spec <- PLOTS[i, ]
  summ <- summarise_timing(load_timing(spec$setup)) %>% filter(surv_model_type == spec$surv)
  fname <- file.path(dir.save, sprintf("timing_%s_surv-%s.pdf", spec$setup, spec$surv))
  if (nrow(summ) == 0) { cat("Skipping", fname, ": no results\n"); next }

  ggsave(fname, plot_timing(summ), device = cairo_pdf, width = 6.5, height = 2)
  cat("Written", fname, "\n")

  # Fitted scaling exponent (time ~ n^a) per m and adjustment, for the text
  summ %>%
    group_by(m, fs) %>%
    summarise(exponent = coef(lm(log(mean_time) ~ log(n_cal)))[2], .groups = "drop") %>%
    print()
}
