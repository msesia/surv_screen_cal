library(tidyverse)

## =============================================================================
## Predicted yield and PPV with error bars given by the simultaneous prediction
## intervals, at the thresholds selected for each target yield, as bar charts
## (rows: Yield/PPV; columns: t0; bars: one per model in a group).
##
## Each row of PLOTS specifies one figure:
##   setup    results directory (results_hpc/<setup>)
##   color    which model type the bars/colours distinguish:
##            "surv" -> different survival models, fixed censoring model `cens`
##            "cens" -> different censoring models, fixed survival model `surv`
##   surv/cens the fixed one (pass NA for the varying dimension)
##   m        cohort size
##   num_vars 0 = models fitted on all variables; k = top-k variables
##            (only present in variable-reduction results)
## Target yields are {0.2, 0.5} for m = 100 and {0.02, 0.05} for m = 1000.
## File name: fig_bars_surv-<surv|all>_cens-<cens|all>_m<m>[_vars<k>].pdf
## =============================================================================

dir.save <- "figures"
dir.create(dir.save, showWarnings = FALSE)

PLOTS <- tribble(
  ~setup,                     ~color, ~surv, ~cens, ~m,   ~num_vars,
  "real_selection_1",         "surv", NA,    "grf", 100,  0,
  "real_selection_1",         "surv", NA,    "grf", 1000, 0,
#  "real_selection_varsel_1",  "surv", NA,    "grf", 100,  5,
  "real_selection_3",         "cens", "grf", NA,    1000, 0, 
  "real_selection_3",         "cens", "xgb", NA,    1000, 0 
)

model_labels <- c(cox = "Cox", grf = "GRF", rsf = "RSF", xgb = "XGB-AFT")
target_yields_by_m <- list("100" = c(0.2, 0.5), "1000" = c(0.02, 0.05))

## ---- Loading -----------------------------------------------------------------

load_selections <- function(dir) {
  d <- list.files(dir, pattern = "\\.txt$", full.names = TRUE) %>%
    discard(~ str_detect(.x, "_(eval|auc|importance)\\.txt$")) %>%
    map_dfr(~ read_delim(.x, delim = ",", col_types = cols(), show_col_types = FALSE)) %>%
    distinct() %>%
    filter(model %in% names(model_labels))          # single-model runs only
  if (!"num_vars" %in% names(d)) d$num_vars <- 0     # full-variable runs
  if (!"batch" %in% names(d))    d$batch    <- 1
  d
}

## ---- Plot --------------------------------------------------------------------

make_bar_plot <- function(d, color_var, target_yields, font_size = 10) {

  plot_df <- bind_rows(
    d %>% transmute(t0, group = .data[[color_var]], target_yield, var_label = "Yield",
                    pred = yield, pred_lo = yield_lo, pred_hi = yield_hi),
    d %>% transmute(t0, group = .data[[color_var]], target_yield, var_label = "PPV",
                    pred = ppv,   pred_lo = ppv_lo,   pred_hi = ppv_hi)
  ) %>%
    mutate(
      group        = factor(model_labels[group], levels = model_labels),
      target_label = factor(target_yield, levels = sort(target_yields)),
      var_label    = factor(var_label, levels = c("Yield", "PPV")),
      t0_label     = paste0("t[0]*': ", t0, " months'")
    )

  legend_title <- if (color_var == "model") "Survival model" else "Censoring model"

  ggplot(plot_df, aes(x = target_label, y = pred, fill = group)) +
    geom_col(position = "dodge") +
    geom_errorbar(aes(ymin = pred_lo, ymax = pred_hi),
                  position = position_dodge(width = 0.9), width = 0.2) +
    facet_grid(var_label ~ t0_label, scales = "free_y",
               labeller = labeller(var_label = label_value, t0_label = label_parsed)) +
    labs(x = "Target yield", y = NULL, fill = legend_title) +
    theme_bw() +
    theme(legend.position  = "right",
          panel.grid.minor = element_blank(),
          strip.placement  = "outside",
          axis.text        = element_text(size = font_size),
          axis.title       = element_text(size = font_size),
          strip.text       = element_text(size = font_size),
          legend.text      = element_text(size = font_size),
          legend.title     = element_text(size = font_size))
}

## ---- Run ---------------------------------------------------------------------

selections <- map(set_names(unique(PLOTS$setup)),
                  ~ load_selections(file.path("results_hpc", .x)))

for (i in seq_len(nrow(PLOTS))) {
  spec <- PLOTS[i, ]
  target_yields <- target_yields_by_m[[as.character(spec$m)]]

  d <- selections[[spec$setup]] %>%
    filter(target_m == spec$m, num_vars == spec$num_vars, target_yield %in% target_yields)

  if (spec$color == "surv") {
    d <- d %>% filter(cens_model_type == spec$cens)
    color_var <- "model"
  } else {
    d <- d %>% filter(model == spec$surv)
    color_var <- "cens_model_type"
  }

  # One selection per (t0, group, target yield): keep the first batch if several
  d <- d %>%
    group_by(t0, .data[[color_var]], target_yield) %>%
    slice_min(batch, n = 1, with_ties = FALSE) %>%
    ungroup()

  fname <- sprintf("fig_bars_surv-%s_cens-%s_m%d%s",
                   if (spec$color == "surv") "all" else spec$surv,
                   if (spec$color == "cens") "all" else spec$cens,
                   spec$m,
                   if (spec$num_vars > 0) sprintf("_vars%d", spec$num_vars) else "")

  if (nrow(d) == 0) { cat("Skipping", fname, ": no results\n"); next }

  p <- make_bar_plot(d, color_var = color_var, target_yields = target_yields)
  ggsave(file.path(dir.save, paste0(fname, ".pdf")), p,
         device = cairo_pdf, height = 3, width = 7, units = "in")
  cat("Written", fname, "\n")
}
