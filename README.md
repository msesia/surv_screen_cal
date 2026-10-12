# SurvScreenCal

*Uncertainty-aware survival screening rules for oncology trial eligibility*

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
source("quickstart.R")
```

The quickstart fits several survival models (Cox, generalized random forest, random survival forest, XGBoost-AFT) and a censoring model on training data. It then computes simultaneous confidence and prediction bands for yield and PPV on calibration data, jointly over all models and thresholds. It compares these with single-model and pointwise bands, illustrates a data-driven choice of model and threshold, and plots the results.

---

## Data

`data/synthetic_survival_data.csv`: 7,000 synthetic patients (time in months, ~72% event rate, 45 variables).

The real Flatiron Health data analyzed in the accompanying paper are not publicly available.

---

## Citation

> V. Svetnik, M. Sesia, M. Johnson, P. Tao, and M. Burcu. *Reproducible life-expectancy eligibility criteria for oncology studies: uncertainty quantification for machine learning-derived survival screening rules.* Preprint, 2026.

