# SurvScreenCal

*Uncertainty-aware survival screening rules for oncology trial eligibility*

---

## What this does

Many oncology trials require patients to have a minimum life expectancy (e.g. survive at least 3, 6, or 9 months post-enrollment). In prospective settings this is assessed by clinician judgment; in retrospective EHR studies it is rarely documented and cannot be directly replicated.

Survival prediction models offer a natural solution: screen in a patient if and only if their predicted probability of surviving past a target horizon $t_0$ exceeds a threshold $\lambda$. Two quantities then matter clinically: **yield** (proportion of patients selected) and **PPV** (proportion of selected patients who actually survive past $t_0$). Raising $\lambda$ increases PPV but reduces yield; the goal is to find a threshold that achieves an acceptable balance.

This is complicated in practice because survival models are often miscalibrated outside their training data — a threshold of 0.9 does not guarantee 90% of selected patients survive — so yield and PPV must be estimated empirically from a held-out calibration dataset.

This repository implements a bootstrap-based method that does exactly this. Given any predictive survival model and a retrospective calibration dataset from the target population, it estimates yield and PPV across all candidate thresholds and constructs uncertainty bands that are **simultaneously valid across all thresholds**. This means a clinician can inspect the full yield–PPV trade-off curve, choose a threshold after seeing the data, and still obtain a valid uncertainty statement. Two band types are provided:

- **Confidence bands** — uncertainty from the finite calibration sample
- **Prediction bands** — additionally account for random variation in the composition of any specific future cohort (wider for smaller future cohorts)

No assumptions on model calibration are required.

---

## Structure

```
R/           # Core method implementation
examples/    # Worked usage examples  
data/        # Synthetic dataset
analysis/    # Code for data analysis
```

---

## Quick start

```r
install.packages(c("tidyverse", "survival", "R6", "survival", "randomForestSRC", "grf", "xgboost","survminer"))
source("examples/example_basic.R")
```

---

## Data

`data/synthetic_survival_data.csv` — 7,000 synthetic patients (time in months, ~72% event rate, 45 variables).

The real Flatiron data set analyzed in the accompanying paper is not publicly available.

---

## Citation

> V. Svetnik, M. Sesia, M. Johnson, P. Tao, and M. Burku. *Reproducible Life-Expectancy Eligibility Criteria for Oncology Studies: Uncertainty Quantification for Machine Learning-Derived Survival Screening Rules.* preprint
