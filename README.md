# SurvScreenCal

*Uncertainty quantification for machine learning-derived survival screening rules in oncology*

---

## What this does

Many oncology trials require patients to have a minimum life expectancy (e.g. survive at least 3, 6, or 9 months after enrollment). In prospective settings this is assessed by clinician judgment; in retrospective EHR studies it is rarely documented and cannot be directly replicated.

Survival prediction models offer a natural solution: screen in a patient if and only if their predicted probability of surviving past a target horizon $t_0$ exceeds a threshold $\lambda$. Two quantities then matter clinically: **yield** (proportion of patients selected) and **PPV** (proportion of selected patients who actually survive past $t_0$). Raising $\lambda$ tends to increase PPV but reduces yield; the goal is to find a threshold that achieves an acceptable balance.

This is complicated in practice because survival models are often miscalibrated outside their training data (a threshold of 0.9 does not guarantee that 90% of selected patients survive), survival data are right-censored, and any future cohort has a random composition. Yield and PPV must therefore be estimated empirically from a held-out calibration dataset, with appropriate uncertainty.

This repository implements a bootstrap-based method that does exactly this. Given one or more predictive survival models and a right-censored calibration dataset from the target population, it estimates yield and PPV across all candidate thresholds, correcting for censoring by inverse probability of censoring weighting, and constructs uncertainty bands that are **simultaneously valid across all thresholds and all candidate models**. A clinician can therefore inspect the full yield–PPV trade-offs, choose both a model and a threshold after seeing the data, and still obtain a valid uncertainty statement for the chosen screening rule. Two band types are provided:

- **Confidence bands**: uncertainty about the population-level yield and PPV, arising from the finite calibration sample
- **Prediction bands**: additionally account for random variation in the composition of a specific future cohort of size $m$ (wider for smaller cohorts)

The method is model-agnostic: survival models are treated as black boxes, and no assumptions on their calibration are required. A finite-sample adjustment keeps the bands reliable when the calibration sample, or the number of selected patients, is small. A pointwise bootstrap baseline, valid only at a pre-specified threshold, is also included for comparison.

---

## Structure

```
R/           # Core method implementation
examples/    # Worked usage examples
data/        # Synthetic dataset
analysis/    # Code for the simulations and data analysis in the paper
```

---

## Quick start

```r
install.packages(c("tidyverse", "survival", "R6", "randomForestSRC", "grf", "xgboost", "survminer"))
setwd("examples")
source("example_basic.R")
```

---

## Example: `examples/example_basic.R`

The example walks through a complete analysis on the synthetic dataset, using a screening horizon of $t_0 = 3$ months:

1. **Split the data** into training (2,000 patients), calibration (2,000), and test (1,000) sets. Only the calibration set is used for inference; the test set is optional.
2. **Fit the models on the training data**: four candidate survival models (Cox, generalized random forest, random survival forest, XGBoost-AFT), which score each patient by their predicted probability of surviving past $t_0$, and a censoring model used to correct for right-censoring.
3. **Prepare the calibration inputs**: the scores of each model, the event indicators, and the estimated censoring probabilities used for inverse probability of censoring weighting.
4. **Compute simultaneous bands** for yield and PPV over a grid of thresholds $\lambda \in \{0, 0.01, \dots, 1\}$, for a future cohort of size $m = 100$:
   - jointly valid over all four models and all thresholds (`bands_joint`), and
   - for comparison, valid over thresholds for each model on its own (`bands_single`).
5. **Compute the pointwise bootstrap baseline** (`bands_pointwise`), whose intervals are valid only at a pre-specified threshold.
6. **Select a model and threshold** with an explicit rule: among all (model, threshold) pairs whose predictive lower bound on yield is at least 20%, choose the one maximizing the predictive lower bound on PPV. Because the joint bands are simultaneous, the reported intervals remain valid for this (or any other) data-driven choice.
7. **Plot the prediction bands** for each model, contrasting the joint and single-model bands, and save the figure to `examples/figures/`.

Each band table is in long format, with one row per model, threshold, and band type (`confidence` or `prediction`). Its columns include the point estimates (`yield_hat`, `ppv_hat`) and the band limits (`yield_lo`, `yield_hi`, `ppv_lo`, `ppv_hi`).

Common modifications, described at the end of the script, include using a single model, changing the horizon $t_0$ or the cohort size $m$, and disabling the finite-sample adjustment.

---

## Data

`data/synthetic_survival_data.csv`: 7,000 synthetic patients (time in months, ~72% event rate, 45 variables).

The real Flatiron Health data analyzed in the accompanying paper are not publicly available.

---

## Citation

> V. Svetnik, M. Sesia, M. Johnson, P. Tao, and M. Burcu. *Reproducible life-expectancy eligibility criteria for oncology studies: uncertainty quantification for machine learning-derived survival screening rules.* Preprint, 2026.
