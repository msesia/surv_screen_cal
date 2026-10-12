library(tidyverse)
library(colorspace)

# =============================================================================
# Yield and PPV vs. lambda, with simultaneous confidence and prediction bands.
#
# Input: long-format bands with columns
#   lambda, band_method, target_m (Inf for confidence bands, m for prediction
#   bands), band_type ("confidence" / "prediction"), yield_lo, yield_hi,
#   ppv_lo, ppv_hi, model
# and t0 (screening time), as written by the launcher-driven R script.
# Point estimates are read from `yield` / `ppv` columns if present, otherwise
# taken as the midpoint of the confidence band.
# =============================================================================

dir.save <- "results_hpc/real_bands_1"

# ---- Loading -----------------------------------------------------------------

# Load every result file written by the launcher (one per (m, t0) argument
# combination; each contains the bands for all horizons and both cohort
# sizes, so duplicates are removed). The row-level horizon is `t0`.
load_all_bands <- function(dir) {
  files <- list.files(dir, pattern = "\\.txt$", full.names = TRUE)
  map_dfr(files, ~ read_delim(.x, delim = ",", col_types = cols(), show_col_types = FALSE)) %>%
    select(-any_of(c("screening_time", "m", "batch"))) %>%   # launcher arguments, not row-level
    distinct() %>%
    rename(screening_time = t0)
}

bands_all <- load_all_bands(dir.save)

# ---- Plot --------------------------------------------------------------------

plot_yield_ppv_bands <- function(bands,
                                 model_name  = "grf",
                                 band_method = "Bootstrap-sim",
                                 m_values    = c(100, 1000),
                                 horizons    = c(3, 6, 9),
                                 xlim        = c(0.2, 1),
                                 font_size   = 10) {

  qty_levels <- c("Yield", "PPV")
  t_levels   <- paste0(horizons, " months")
  t_labels   <- setNames(paste0("t[0]*': ", horizons, " months'"), t_levels)

  df0 <- bands %>%
    filter(model == model_name, band_method == !!band_method,
           screening_time %in% horizons) %>%
    mutate(t_star_f = factor(paste0(screening_time, " months"), levels = t_levels))

  # ---- Point estimates: from yield/ppv columns if present, else CB midpoint
  conf0 <- df0 %>% filter(band_type == "confidence")
  has_est <- all(c("yield", "ppv") %in% names(df0))
  est_df <- conf0 %>%
    transmute(t_star_f, lambda,
              yield_hat = if (has_est) yield else (yield_lo + yield_hi) / 2,
              ppv_hat   = if (has_est) ppv   else (ppv_lo   + ppv_hi)   / 2)
  if (!has_est) message("No point-estimate columns found; using confidence-band midpoints.")

  line_df <- bind_rows(
    est_df %>% transmute(t_star_f, lambda, quantity = "Yield", estimate = yield_hat),
    est_df %>% transmute(t_star_f, lambda, quantity = "PPV",   estimate = ppv_hat)
  ) %>% mutate(quantity = factor(quantity, levels = qty_levels))

  # ---- Bands in long form: one row per (horizon, lambda, quantity, band)
  to_long <- function(d, band_label) {
    bind_rows(
      d %>% transmute(t_star_f, lambda, quantity = "Yield", band = band_label,
                      lo = yield_lo, hi = yield_hi),
      d %>% transmute(t_star_f, lambda, quantity = "PPV",   band = band_label,
                      lo = ppv_lo, hi = ppv_hi)
    )
  }

  conf_df <- to_long(conf0, "Confidence")
  pred_df <- map_dfr(m_values, function(m_val) {
    df0 %>%
      filter(band_type == "prediction", target_m == m_val) %>%
      to_long(paste0("Prediction (m = ", m_val, ")"))
  })

  # Legend order: confidence first; drawing order: widest (smallest m) first
  band_levels <- c("Confidence", paste0("Prediction (m = ", rev(m_values), ")"))
  band_colors <- setNames(c("#2166AC", "#DF8E89", "#FDDBC7")[seq_along(band_levels)],
                          band_levels)

  plot_df <- bind_rows(conf_df, pred_df) %>%
    mutate(quantity = factor(quantity, levels = qty_levels),
           band     = factor(band, levels = rev(band_levels)))

  p <- ggplot(line_df, aes(x = lambda, y = estimate)) +
    geom_ribbon(data = arrange(plot_df, band),
                aes(x = lambda, ymin = lo, ymax = hi, fill = band),
                alpha = 1, inherit.aes = FALSE) +
    scale_fill_manual(values = band_colors, breaks = band_levels) +
    geom_line(color = "black", linewidth = 1.0) +
    facet_grid(quantity ~ t_star_f, scales = "fixed",
               labeller = labeller(t_star_f = as_labeller(t_labels, label_parsed))) +
    coord_cartesian(xlim = xlim) +
    scale_x_continuous(breaks = seq(0.25,1, by = 0.25)) +
    labs(x = expression("Selection threshold " * lambda), y = NULL, fill = "Band type") +
    theme_bw() +
    theme(legend.position  = "right",
          panel.grid.minor = element_blank(),
          strip.placement  = "outside",
          axis.text        = element_text(size = font_size),
          axis.title       = element_text(size = font_size),
          strip.text       = element_text(size = font_size),
          legend.text      = element_text(size = font_size),
          legend.title     = element_text(size = font_size))

  p_bw <- suppressMessages(
    p + scale_fill_manual(values = setNames(colorspace::desaturate(band_colors), band_levels),
                          breaks = band_levels)
  )

  list(plot = p, plot_bw = p_bw)
}

# ---- Run ---------------------------------------------------------------------

out <- plot_yield_ppv_bands(bands_all, model_name = "grf", m_values = c(100, 1000))
out$plot

dir.create("figures", showWarnings = FALSE)
ggsave(file.path("figures", "Yield_PPV_vs_lambda_cal10K_3times_m100_m1000_3by2.pdf"),
       plot = out$plot, device = cairo_pdf, width = 9, height = 3.5)
ggsave("preview_bw.pdf", out$plot_bw, width = 9, height = 4)
