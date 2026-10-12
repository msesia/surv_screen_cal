library(tidyverse)
library(kableExtra)

# =============================================================================
# Table 2 / S5: relative change (%) in predicted PPV and in the width of the
# PPV prediction interval when reducing from the full (all variables) to the
# reduced (top-k) model, at selected target yields, one table per cohort size.
#
# Input: results/real_varsel_1/surv_<model>_..._m<m>_time<t0>_batch<b>.txt
# (selection files), rows for num_vars = 0 (all) and num_vars = k.
# Changes are computed within batch and then averaged over batches.
# =============================================================================

dir.results <- "results_hpc/real_varsel_1"
dir.tables  <- "tables"
dir.create(dir.tables, showWarnings = FALSE)

model_labels  <- c(cox = "Cox", grf = "GRF", rsf = "RSF", xgb = "XGB-AFT")
num_vars_red  <- 10
target_yields_by_m <- list("100" = c(0.2, 0.5), "1000" = c(0.02, 0.05))

# ---- Loading -----------------------------------------------------------------

sel_all <- list.files(dir.results, pattern = "\\.txt$", full.names = TRUE) %>%
  discard(~ str_detect(.x, "_(eval|auc|importance)\\.txt$")) %>%
  map_dfr(~ read_delim(.x, delim = ",", col_types = cols(), show_col_types = FALSE)) %>%
  distinct() %>%
  filter(model %in% names(model_labels), num_vars %in% c(0, num_vars_red)) %>%
  mutate(width_ppv = ppv_hi - ppv_lo)

# ---- Relative changes --------------------------------------------------------

changes <- sel_all %>%
  select(model, t0, target_m, target_yield, batch, num_vars, ppv, width_ppv) %>%
  pivot_wider(names_from = num_vars, values_from = c(ppv, width_ppv),
              names_glue = "{.value}_v{num_vars}") %>%
  mutate(change_ppv   = 100 * (ppv_v10 / ppv_v0 - 1),
         change_width = 100 * (width_ppv_v10 / width_ppv_v0 - 1)) %>%
  group_by(model, t0, target_m, target_yield) %>%
  summarise(change_ppv = mean(change_ppv), change_width = mean(change_width),
            n_batches = n(), .groups = "drop")

# ---- Tabular -----------------------------------------------------------------

make_change_tabular <- function(ch, m_value, target_yields) {
  d <- ch %>%
    filter(target_m == m_value, target_yield %in% target_yields) %>%
    mutate(model = factor(model, levels = names(model_labels))) %>%
    arrange(model, t0) %>%
    select(model, t0, target_yield, change_ppv, change_width) %>%
    pivot_wider(names_from = target_yield, values_from = c(change_ppv, change_width),
                names_glue = "{.value}_{target_yield}") %>%
    transmute(Model = model_labels[as.character(model)], t0,
              across(starts_with("change_"), ~ sprintf("%.3f", .x)))

  yields <- sort(target_yields)
  col_names <- c("Model", "T0HEADER",
                 paste0("Yield ", yields), paste0("Yield ", yields, " "))
  header <- c(" " = 2)
  header["Change (\\\\%) in PPV"]   <- length(yields)
  header["Change (\\\\%) in width"] <- length(yields)

  kbl(d, format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",
      col.names = col_names, align = c("l", "c", rep("r", 2 * length(yields)))) %>%
    add_header_above(header, escape = FALSE) %>%
    collapse_rows(columns = 1, valign = "middle", latex_hline = "major") %>%  # requires \usepackage{multirow}
    as.character() %>%
    str_replace(fixed("T0HEADER"), "\\makecell{$t_0$\\\\(months)}")   # requires \usepackage{makecell}
}

# ---- Run ---------------------------------------------------------------------

for (m in sort(unique(changes$target_m))) {
  ty <- target_yields_by_m[[as.character(m)]]
  if (is.null(ty)) { cat("No target yields defined for m =", m, "\n"); next }
  outfile <- file.path(dir.tables, sprintf("ppv_diff_vars%d_m%d.tex", num_vars_red, m))
  writeLines(make_change_tabular(changes, m_value = m, target_yields = ty), outfile)
  cat("Written", outfile, "\n")
}
