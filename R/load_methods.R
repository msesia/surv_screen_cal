library(survival)

## Get directory of this file
this_file <- normalizePath(sys.frame(1)$ofile)
this_dir  <- dirname(this_file)

# Helper to source relative to this file
source_rel <- function(path) {
  source(file.path(this_dir, path))
}

source_rel("utils/utils_splitting.R")
source_rel("utils/utils_weights_scores.R")
source_rel("utils/utils_misc.R")
source_rel("utils/utils_plotting.R")
source_rel("utils/utils_imputation.R")
source_rel("utils/utils_semi_synthetic_data.R")
source_rel("utils/utils_evaluation.R")

source_rel("wrappers/censoring_wrappers.R")
source_rel("wrappers/survival_wrappers.R")

source_rel("methods/calibration_bands.R")
source_rel("methods/pointwise_bootstrap.R")
source_rel("methods/simultaneous_bootstrap.R")
source_rel("methods/pointwise_finite_sample.R")
source_rel("methods/pointwise_selection_event.R")
