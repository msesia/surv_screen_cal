library(tidyverse)
library(kableExtra)

# Table S1: baseline covariates, summarised on observed (non-imputed) values

data_path <- if (file.exists("/media/msesia/ExternalSSD/data/fhrd/train_test_miss.rds")) {
  "/media/msesia/ExternalSSD/data/fhrd/train_test_miss.rds"
} else {
  "/project2/sesia_1123/data/exploration/train_test_miss.rds"
}
back_transform <- exp   # inverse of the log used for "log_" variables
max_levels     <- 10    # at most this many distinct values -> categorical

# Variables used in the analysis (code -> description), in table order
VAR_LABELS <- c(
  alb_ctns = "Albumin", ecog_value = "ECOG performance status",
  age = "Age", NLR = "Neutrophil-to-lymphocyte ratio", oxg = "Oxygen saturation",
  log_heartrate = "Heart rate", chloride = "Chloride", hgb_ctns = "Hemoglobin",
  log_mono_ctns = "Monocyte count", log_alp_ctns = "Alkaline phosphatase",
  log_rbc_ctns = "Red blood cell count", prot_ctns = "Total protein",
  log_sysbp = "Systolic blood pressure", LDH = "Lactate dehydrogenase",
  log_ca_ctns = "Calcium", sodium_ctns = "Sodium", log_neut_ctns = "Neutrophil count",
  log_creat_ctns = "Creatinine", log_eos_leu_ctns = "Eosinophil-to-leukocyte ratio",
  log_ast_alt_ctns = "AST-to-ALT ratio", log_alt_ctns = "Alanine aminotransferase",
  plat_ctns = "Platelet count", log_bili_ctns = "Bilirubin", k_ctns = "Potassium",
  temperature = "Body temperature", Number_of_metastatic_sites = "Number of metastatic sites",
  log_gluc_ctns = "Glucose", log_baso_ctns = "Basophil count", log_bmi = "Body mass index",
  liver_cag = "Liver metastases", SquamousCell = "Tumor histology", race = "Race",
  tumor_stage = "Tumor stage", bone_cag = "Bone metastases", sex_code = "Sex",
  insure = "Insurance type", smoking_code = "Smoking status",
  ses_index_2015_2019 = "Socioeconomic index", trtN30_cag = "Medications in past 30 days",
  brain_cag = "Brain metastases", labN30_cag = "Lab visits in past 30 days",
  vitalN30_cag = "Vital sign records in past 30 days", CCI_cag = "Comorbidity index"
)

# Readable level labels, overriding the raw values in the data.
# CHECK the 0/1 codings (e.g. table(raw$sex_code)) before using these.
LEVEL_LABELS <- c(
  sex_code     = "Female; Male",
  smoking_code = "No; Yes",
  SquamousCell = "Non-squamous; Squamous",
  insure       = "Commercial with Medicare/Medicaid; Commercial without Medicare/Medicaid; Medicare/Medicaid and other; Medicare/Medicaid only; Other only",
  ses_index_2015_2019 = "1 (lowest)--5 (highest)",
  CCI_cag      = "$\\le 6$; $>6$",
  vitalN30_cag = "$<2$; $\\ge 2$",
  trtN30_cag   = "$<2$; $\\ge 2$",
  labN30_cag   = "$<2$; $\\ge 2$"
)

# Integer-coded ordinal variables reported as a range (min--max)
ORDINAL_VARS <- c("ecog_value", "tumor_stage", "Number_of_metastatic_sites")

# ---- Load data: keep only the analysis variables ------------------------------

raw <- as_tibble(do.call("rbind", readRDS(data_path)))
covars <- raw %>% select(all_of(intersect(names(raw), names(VAR_LABELS))))

extra <- setdiff(names(raw), c(names(VAR_LABELS), "month", "event"))
if (length(extra) > 0) message("Dropped columns not used in the analysis: ", paste(extra, collapse = ", "))
not_found <- setdiff(names(VAR_LABELS), names(raw))
if (length(not_found) > 0) message("Not found in data: ", paste(not_found, collapse = ", "))

# ---- Summaries -----------------------------------------------------------------

fmt_num <- function(x) formatC(signif(x, 3), format = "fg", digits = 3)

# Escape LaTeX special characters in levels taken from the data
tex_escape <- function(x) {
  x %>%
    str_replace_all("\\\\", "\\\\textbackslash{}") %>%
    str_replace_all("([&%$#_{}])", "\\\\\\1") %>%
    str_replace_all("<=", "$\\\\le$") %>%
    str_replace_all(">=", "$\\\\ge$") %>%
    str_replace_all("(?<![$\\\\a-z])<", "$<$") %>%
    str_replace_all("(?<![$\\\\a-z])>", "$>$")
}

summarise_var <- function(x, code) {
  x_obs <- x[!is.na(x)]
  if (code %in% names(LEVEL_LABELS)) {
    type <- 2; summ <- LEVEL_LABELS[[code]]
  } else if (code %in% ORDINAL_VARS) {
    lv  <- unique(as.character(x_obs))
    num <- suppressWarnings(as.numeric(lv))
    lv  <- if (!anyNA(num)) lv[order(num)] else sort(lv)   # numeric order when possible
    type <- 2; summ <- sprintf("%s--%s", lv[1], lv[length(lv)])
  } else if (is.numeric(x_obs) && n_distinct(x_obs) > max_levels) {
    type <- 1
    if (str_starts(code, "log_")) x_obs <- back_transform(x_obs)
    summ <- sprintf("%s (%s)", fmt_num(median(x_obs)), fmt_num(IQR(x_obs)))
  } else {
    type <- 2; summ <- tex_escape(paste(sort(unique(as.character(x_obs))), collapse = "; "))
  }
  tibble(type = type,
         Code = paste0("\\texttt{", str_replace_all(code, "_", "\\\\_"), "}"),
         Variable = VAR_LABELS[[code]],
         Summary = summ,
         Missing = if (any(is.na(x))) sprintf("%.1f", 100 * mean(is.na(x))) else "--")
}

tab <- map2_dfr(covars, names(covars), summarise_var) %>%
  select(-type)

# ---- LaTeX tabular -------------------------------------------------------------

tex <- kbl(tab, format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",
           col.names = c("Code", "Variable", "Median (IQR) or levels", "Missing (\\%)"),
           align = c("l", "l", "l", "r")) %>%
  column_spec(2, width = "5cm") %>%     # wrap long descriptions
  column_spec(3, width = "6cm") %>%     # wrap long level lists
  as.character()

dir.create("tables", showWarnings = FALSE)
writeLines(tex, "tables/variables.tex")
cat("Written tables/variables.tex (", nrow(tab), "variables )\n")
