library(tidyverse)
library(kableExtra)

# Threshold-selection tables (Table 1 / S4 and their joint-rule analogues).
#   real_selection_1: rule (7) on a single survival model.
#   real_selection_2: joint rule over (model, lambda) based on prediction bounds.
# One table per (setup, cohort size m); target yields {0.2,...,0.5} for
# m = 100 and {0.02,...,0.05} for m = 1000.

dir.tables <- "tables"
dir.create(dir.tables, showWarnings = FALSE)

model_labels <- c(cox = "Cox", grf = "GRF", rsf = "RSF", xgb = "XGB-AFT")
target_yields_by_m <- list("100"  = c(0.2, 0.3, 0.4, 0.5),
                           "1000" = c(0.02, 0.03, 0.04, 0.05))

# ---- Loading -----------------------------------------------------------------

load_selection <- function(setup) {
  list.files(file.path("results_hpc", setup), pattern = "\\.txt$", full.names = TRUE) %>%
    discard(~ str_detect(.x, "_(eval|auc|importance)\\.txt$")) %>%   # selection files only
    map_dfr(~ read_delim(.x, delim = ",", col_types = cols(), show_col_types = FALSE)) %>%
    distinct()
}

# ---- Tabular -----------------------------------------------------------------

make_selection_tabular <- function(sel, m, target_yields, horizons = c(3, 6, 9)) {
  m_val <- m
  d <- sel %>%
    filter(target_m == m_val, t0 %in% horizons, target_yield %in% target_yields) %>%
    arrange(t0, target_yield) %>%
    transmute(t0, target_yield,
              model_txt  = model_labels[model],
              lambda_txt = sprintf("%.2f", lambda),
              yield_txt  = sprintf("%.4f", yield),
              ppv_txt    = sprintf("%.3f", ppv),
              yield_pi   = sprintf("[%.2f, %.2f]", yield_lo, yield_hi),
              ppv_pi     = sprintf("[%.3f, %.3f]", ppv_lo, ppv_hi),
              nsel       = sprintf("[%d, %d]", round(m_val * yield_lo), round(m_val * yield_hi)))

  block_ends <- cumsum(rle(d$t0)$lengths)
  block_ends <- block_ends[-length(block_ends)]

  header <- c(" " = 4, "Estimate" = 2)
  header[sprintf("95\\\\%% Prediction Interval ($m = %d$)", m)] <- 3

  kbl(d, format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",
      col.names = c("T0HEADER", "TYHEADER", "Model", "$\\hat{\\lambda}$",
                    "Yield", "PPV", "Yield ", "PPV ", "No. selected"),
      align = rep("c", ncol(d))) %>%
    add_header_above(header, escape = FALSE) %>%
    row_spec(block_ends, hline_after = TRUE) %>%
    as.character() %>%
    str_replace_all(fixed("\\hline"), "\\midrule") %>%
    str_replace(fixed("T0HEADER"), "\\makecell{$t_0$\\\\(months)}") %>%   # requires \usepackage{makecell}
    str_replace(fixed("TYHEADER"), "\\makecell{Target\\\\yield}")
}

write_tex <- function(tex, name) {
  outfile <- file.path(dir.tables, name)
  writeLines(tex, outfile)
  cat("Written", outfile, "\n")
}

# ---- Run ---------------------------------------------------------------------

# Setup 1: rule (7), GRF model; one table per cohort size
sel_1 <- load_selection("real_selection_1") %>% filter(model == "grf", cens_model_type == "grf")
for (m in sort(unique(sel_1$target_m))) {
  write_tex(make_selection_tabular(sel_1, m = m, target_yields = target_yields_by_m[[as.character(m)]]),
            sprintf("lambda_selection_1_m%d.tex", m))
}

# Setup 2: joint (model, lambda) rule over the "mix" set
sel_2 <- load_selection("real_selection_2") %>% filter(surv_model_type == "mix")
for (m in sort(unique(sel_2$target_m))) {
  write_tex(make_selection_tabular(sel_2, m = m, target_yields = target_yields_by_m[[as.character(m)]]),
            sprintf("lambda_selection_2_m%d.tex", m))
}
