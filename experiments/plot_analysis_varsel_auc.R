library(tidyverse)
library(kableExtra)

# =============================================================================
# Variable-reduction analysis outputs:
#  (1) Median time-dependent AUC as a function of predictor-set size, by model,
#      for t0 = 3, 6, 9 months (Figure S1), one figure per cohort size m.
#  (2) Table S3: top-10 variables per model by importance averaged over batches.
#
# Input: results/real_varsel_1/surv_<model>_..._time<t0>_batch<b>_auc.txt and
#        ..._importance.txt; num_vars = 0 means all variables.
# =============================================================================

dir.results <- "results_hpc/real_varsel_1"
dir.save    <- "figures"
dir.tables  <- "tables"
dir.create(dir.save,   showWarnings = FALSE)
dir.create(dir.tables, showWarnings = FALSE)

model_labels <- c(cox = "Cox", rsf = "RSF", grf = "GRF", xgb = "XGB-AFT")
n_vars_total <- 43

# Clinical labels for the raw variable names (Table S1)
variable_labels <- c(
  alb_ctns = "Albumin", ecog_value = "ECOG", age = "Age", NLR = "NLR",
  oxg = "Oxygen", log_heartrate = "Heart rate", chloride = "Chloride",
  hgb_ctns = "Hemoglobin", log_mono_ctns = "Monocytes",
  log_alp_ctns = "Alkaline phosphatase", log_rbc_ctns = "Red blood cells",
  prot_ctns = "Protein", log_sysbp = "Systolic BP", LDH = "LDH",
  log_ca_ctns = "Calcium", sodium_ctns = "Sodium", log_neut_ctns = "Neutrophils",
  log_creat_ctns = "Creatinine", log_eos_leu_ctns = "Eosinophil/leukocyte ratio",
  log_ast_alt_ctns = "AST/ALT ratio", log_alt_ctns = "ALT", plat_ctns = "Platelets",
  log_bili_ctns = "Bilirubin", k_ctns = "Potassium", temperature = "Temperature",
  Number_of_metastatic_sites = "No. metastatic sites", log_gluc_ctns = "Glucose",
  log_baso_ctns = "Basophils", log_bmi = "BMI", liver_cag = "Liver metastases",
  SquamousCell = "Histology", race = "Race", tumor_stage = "Tumor stage",
  bone_cag = "Bone metastases", sex_code = "Sex", insure = "Insurance",
  smoking_code = "Smoking", ses_index_2015_2019 = "SES index",
  trtN30_cag = "Medications (30 d)", brain_cag = "Brain metastases",
  labN30_cag = "Lab visits (30 d)", vitalN30_cag = "Vital records (30 d)",
  CCI_cag = "CCI"
)
label_variable <- function(v) ifelse(v %in% names(variable_labels), variable_labels[v], v)

read_all <- function(pattern) {
  list.files(dir.results, pattern = pattern, full.names = TRUE) %>%
    map_dfr(~ read_delim(.x, delim = ",", col_types = cols(), show_col_types = FALSE))
}

# ---- (1) AUC vs number of variables -----------------------------------------

auc_all <- read_all("_auc\\.txt$") %>%
  mutate(num_vars = if_else(num_vars == 0, n_vars_total, num_vars))

auc_summ <- auc_all %>%
  group_by(m, screening_time, model, num_vars) %>%
  summarise(median_auc = median(auc, na.rm = TRUE),
            q25 = quantile(auc, 0.25, na.rm = TRUE),
            q75 = quantile(auc, 0.75, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  mutate(model = factor(model_labels[model], levels = model_labels),
         t0_label = factor(paste0(screening_time, "m"),
                           levels = paste0(sort(unique(screening_time)), "m")))

plot_auc_vs_vars <- function(summ, m_value, font_size = 10) {
  d <- summ %>% filter(m == m_value)
  breaks <- c(43, 10, 3)

  ggplot(d, aes(x = num_vars, y = median_auc, color = model, shape = model, group = model)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_grid(. ~ t0_label) +
    scale_x_continuous(trans = scales::compose_trans("log10", "reverse"),
                       breaks = breaks, labels = breaks) +
    labs(x = "Number of variables", y = "Median AUC (across batches)", color = "Model", shape = "Model") +
    theme_bw() +
    theme(legend.position  = "right",
          panel.grid.minor = element_blank(),
          axis.text        = element_text(size = font_size),
          axis.title       = element_text(size = font_size),
          strip.text       = element_text(size = font_size),
          legend.text      = element_text(size = font_size),
          legend.title     = element_text(size = font_size))
}

for (m_val in sort(unique(auc_summ$m))) {
  p <- plot_auc_vs_vars(auc_summ, m_value = m_val)
  ggsave(file.path(dir.save, sprintf("auc_vs_num_vars_m%d.pdf", m_val)), p,
         device = cairo_pdf, width = 9, height = 3)
}

# ---- (2) Table S3: top-10 variables by mean importance -----------------------

# The ranking does not depend on t0 or m (it comes from the full-variable fit),
# so keep one copy per (model, batch) before averaging over batches.
imp_all <- read_all("_importance\\.txt$") %>%
  distinct(model, batch, code, variable, importance)

top10 <- imp_all %>%
  group_by(model, code, variable) %>%
  summarise(median_importance = median(importance), .groups = "drop") %>%
  group_by(model) %>%
  slice_max(median_importance, n = 10, with_ties = FALSE) %>%
  mutate(rank = row_number(-median_importance)) %>%
  ungroup() %>%
  transmute(model, rank, label = label_variable(variable)) %>%
  pivot_wider(names_from = model, values_from = label) %>%
  arrange(rank)

# Column order as in the paper: Rank, Cox, RSF, GRF, XGB
cols <- intersect(names(model_labels), names(top10))
tab <- top10 %>% select(rank, all_of(cols))

tex <- kbl(tab, format = "latex", booktabs = TRUE, escape = TRUE, linesep = "",
           col.names = c("Rank", model_labels[cols]),
           align = c("c", rep("l", length(cols)))) %>%
  as.character()
writeLines(tex, file.path(dir.tables, "var_importance_top10.tex"))
cat("Written", file.path(dir.tables, "var_importance_top10.tex"), "\n")
