library(tidyverse)
library(colorspace)

dir.save <- "results_hpc/real_bands_2"

# =============================================================================
# Prediction bands (m = 100) at a given t0 for four survival models, comparing
# bands simultaneously valid over one model only (single-model runs) with
# bands jointly valid over all four models (the "mix" run).
# Facets: quantity (rows) x model (columns); shading = validity set.
#
# Input: all .txt files in `dir.save` written by the launcher-driven R script;
# duplicated rows (e.g. repeated confidence bands) are removed after loading.
# =============================================================================

# ---- Loading -----------------------------------------------------------------

load_all_bands <- function(dir) {
  files <- list.files(dir, pattern = "\\.txt$", full.names = TRUE)
  map_dfr(files, ~ read_delim(.x, delim = ",", col_types = cols(), show_col_types = FALSE)) %>%
    distinct()
}

bands_all <- load_all_bands(dir.save)
horizon <- 6

# ---- Plot --------------------------------------------------------------------

plot_bands_by_validity <- function(bands, horizon,
                                   models      = c("cox", "grf", "rsf", "xgb"),
                                   band_method = "Bootstrap-sim",
                                   m_val       = 100,
                                   xlim        = c(0.2, 1),
                                   font_size   = 10) {

  qty_levels   <- c("Yield", "PPV")
  model_labels <- c(cox = "Cox", grf = "GRF", rsf = "RSF", xgb = "XGB-AFT")
  K <- length(models)

  n_valid <- function(v) lengths(strsplit(v, "\\|"))

  validity_levels <- c("Single model", "All models")

  df0 <- bands %>%
    filter(t0 == horizon, model %in% models,
           band_method == !!band_method,
           band_type == "prediction", target_m == m_val) %>%
    mutate(k_valid = n_valid(valid_models)) %>%
    filter(k_valid %in% c(1, K)) %>%          # keep single-model and all-model bands only
    mutate(
      validity = factor(ifelse(k_valid == 1, validity_levels[1], validity_levels[2]),
                        levels = validity_levels),
      model_f  = factor(model_labels[model], levels = model_labels[models])
    )

  # Point estimates do not depend on the validity set; take them from the
  # single-model rows.
  line_df <- df0 %>%
    filter(k_valid == 1) %>%
    {bind_rows(
      transmute(., model_f, lambda, quantity = "Yield", estimate = yield),
      transmute(., model_f, lambda, quantity = "PPV",   estimate = ppv)
    )} %>%
    mutate(quantity = factor(quantity, levels = qty_levels))

  band_df <- bind_rows(
    df0 %>% transmute(model_f, lambda, quantity = "Yield", validity, lo = yield_lo, hi = yield_hi),
    df0 %>% transmute(model_f, lambda, quantity = "PPV",   validity, lo = ppv_lo,   hi = ppv_hi)
  ) %>%
    mutate(quantity = factor(quantity, levels = qty_levels),
           # draw the wider (all-model) band first, the single-model band on top
           validity = factor(validity, levels = rev(validity_levels)))

  band_colors <- setNames(c("#FDDBC7", "#DF8E89"), validity_levels)

  p <- ggplot(line_df, aes(x = lambda, y = estimate)) +
    geom_ribbon(data = arrange(band_df, validity),
                aes(x = lambda, ymin = lo, ymax = hi, fill = validity),
                alpha = 1, inherit.aes = FALSE) +
    scale_fill_manual(values = band_colors, breaks = validity_levels) +
    geom_line(color = "black", linewidth = 1.0) +
    facet_grid(quantity ~ model_f, scales = "fixed") +
    coord_cartesian(xlim = xlim) +
    scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1)) +
    labs(x = expression("Selection threshold " * lambda), y = NULL,
         fill = sprintf("Prediction band (m = %d)\nvalid over", m_val)) +
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
    p + scale_fill_manual(values = setNames(colorspace::desaturate(band_colors), validity_levels),
                          breaks = validity_levels)
  )

  list(plot = p, plot_bw = p_bw, data = df0)
}

# ---- Run ---------------------------------------------------------------------

out <- plot_bands_by_validity(bands_all, horizon = horizon, m_val = 100)
out$plot

dir.create("figures", showWarnings = FALSE)
ggsave(file.path("figures", sprintf("Yield_PPV_vs_lambda_validity_t%d_m100.pdf", horizon)),
       plot = out$plot, device = cairo_pdf, width = 9, height = 3.5)
