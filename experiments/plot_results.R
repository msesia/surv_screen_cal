library(tidyverse)

load_results <- function(setup) {
    idir <- sprintf("results_hpc/%s", setup)
    ifile.list <- list.files(idir)
    results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
    }))
    return(results)
}
results <- load_results("v0")

method.values <- c("Bootstrap-sim", "Bootstrap-point")
method.labels <- c("Simultaneous bootstrap", "Pointwise bootstrap")

band.values <- c("prediction", "confidence")
band.labels <- c("Prediction band", "Confidence band")

lambda_colors <- c(
  "Worst-case" = "#000000",
  "λ=0.50"     = "#E69F00",
  "λ=0.60"     = "#56B4E9",
  "λ=0.70"     = "#009E73",
  "λ=0.80"     = "#F0E442",
  "λ=0.90"     = "#0072B2",
  "λ=1.00"     = "#D55E00"
)

lambda_shapes <- c(
  "Worst-case" = 17,
  "λ=0.50"     = 16,
  "λ=0.60"     = 15,
  "λ=0.70"     = 18,
  "λ=0.80"     = 8,
  "λ=0.90"     = 7,
  "λ=1.00"     = 4
)


t0 <- 9

## Evaluate coverage
results_cov <- results %>%
  mutate(
    cover_yield = case_when(
      is.na(yield) & yield_lo == 0 & yield_hi == 1 ~ TRUE,
      !is.na(yield) & !is.na(yield_lo) &
        yield >= yield_lo & yield <= yield_hi ~ TRUE,
      !is.na(yield) & !is.na(yield_lo) ~ FALSE,
      TRUE ~ NA
    ),
    cover_ppv = case_when(
      is.na(ppv) & ppv_lo == 0 & ppv_hi == 1 ~ TRUE,
      !is.na(ppv) & !is.na(ppv_lo) &
        ppv >= ppv_lo & ppv <= ppv_hi ~ TRUE,
      !is.na(ppv) & !is.na(ppv_lo) ~ FALSE,
      TRUE ~ NA
    ),
    width_yield = yield_hi - yield_lo,
    width_ppv    = ppv_hi - ppv_lo
  ) %>%
  mutate(lambda_type = "pointwise") %>%
    mutate(band_method = factor(band_method, method.values, method.labels),
           band_type = factor(band_type, band.values, band.labels)
           )

## Evaluate worst-case coverage (over lambda)
sim_cov_seed <- results_cov %>%
    group_by(
        Seed, real_data, gen_model_type,
        surv_model_type, cens_model_type,
        n_train, n_cal, n_test, screening_time, alpha, batch,
        band_method, band_type, target_m
    ) %>%
    summarise(
        N = n(),
        lambda = NA_real_,  ## keep numeric type
        cover_yield = all(cover_yield),
        cover_ppv    = all(cover_ppv),
        width_yield = max(width_yield, na.rm = TRUE),
        width_ppv    = max(width_ppv, na.rm = TRUE),
        yield = NA_real_,
        ppv    = NA_real_,
        yield_lo = NA_real_,
        yield_hi = NA_real_,
        ppv_lo = NA_real_,
        ppv_hi = NA_real_,
        .groups = "drop"
    ) %>%
    mutate(lambda_type = "worst-case")


results_all <- dplyr::bind_rows(results_cov, sim_cov_seed)

    
summary_all <- results_all %>%
  group_by(n_cal, n_test, alpha, screening_time, lambda_type, lambda, band_method, band_type, target_m) %>%
  summarise(
    R = n(),
    ## Coverage
    yield_coverage = mean(cover_yield, na.rm = TRUE),
    se_yield_cov   = sd(cover_yield, na.rm = TRUE) / sqrt(R),
    ppv_coverage = mean(cover_ppv, na.rm = TRUE),
    se_ppv_cov   = sd(cover_ppv, na.rm = TRUE) / sqrt(R),
    ## Width
    mean_width_yield = mean(width_yield, na.rm = TRUE),
    se_width_yield   = sd(width_yield, na.rm = TRUE) / sqrt(R),
    mean_width_ppv = mean(width_ppv, na.rm = TRUE),
    se_width_ppv   = sd(width_ppv, na.rm = TRUE) / sqrt(R),
    .groups = "drop"
  )


summary.plot <- summary_all %>%
    filter(screening_time == t0) %>%
    mutate(
        lambda_label = ifelse(
            is.na(lambda),
            "Worst-case",
            paste0("λ=", formatC(lambda, format = "f", digits = 2))
        ),
        n_test_label = paste0("Test size = ", n_test)
    )

p.yi.c <- summary.plot %>%
    ggplot(aes(x = n_cal, y = yield_coverage, color = lambda_label, shape = lambda_label)) +
    geom_point() +
    geom_line() +
    geom_errorbar(
        aes(
            ymin = pmax(0, yield_coverage - 1.96 * se_yield_cov),
            ymax = pmin(1, yield_coverage + 1.96 * se_yield_cov)
        ),
        width = 0.02,
        alpha = 0.4
    ) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    facet_grid(band_type ~ band_method + n_test_label) +
    labs(
        y = "Yield coverage",
        x = "Calibration size (n)",
        color = "Screening threshold", shape = "Screening threshold"#,
                                        #title = sprintf("Time: %.1f, Yield coverage with 95%% bands", t0)
    ) +
    scale_color_manual(values = lambda_colors) +
    scale_shape_manual(values = lambda_shapes) +    
    theme_bw()
#p.yi.cc

ggsave(sprintf("figures/experiments_yield_coverage_t%.1f.png", t0), p.yi.c, width=9, height=3.5, units="in", bg = 'white')


p.er.c <- summary.plot %>%
  ggplot(aes(x = n_cal, y = ppv_coverage, color = lambda_label, shape = lambda_label)) +
  geom_point() +
  geom_line() +
  geom_errorbar(
    aes(
      ymin = pmax(0, ppv_coverage - 1.96 * se_ppv_cov),
      ymax = pmin(1, ppv_coverage + 1.96 * se_ppv_cov)
    ),
    width = 0.02,
    alpha = 0.4
  ) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  facet_grid(band_type ~ band_method + n_test_label) +
  labs(
    y = "PPV coverage",
    x = "Calibration size (n)",
    color = "Screening threshold", shape = "Screening threshold"#,
    #title = sprintf("Time: %.1f, PPV coverage with 95%% bands", t0)
  ) +
    scale_color_manual(values = lambda_colors) +
    scale_shape_manual(values = lambda_shapes) +    
  theme_bw()
##p.er.c

ggsave(sprintf("figures/experiments_ppv_coverage_t%.1f.png", t0), p.er.c, width=9, height=3.5, units="in", bg = 'white')

p.yi.w <- summary.plot %>%
  filter(summary.plot$lambda_label !="Worst-case") %>%
  ggplot(aes(x = n_cal, y = mean_width_yield, color = lambda_label, shape = lambda_label)) +
  geom_point() +
  geom_line() +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean_width_yield - 1.96 * se_width_yield),
      ymax = pmin(1, mean_width_yield + 1.96 * se_width_yield)
    ),
    width = 0.02,
    alpha = 0.4
  ) +
  facet_grid(band_type ~ band_method + n_test_label) +
  labs(
    y = "Yield width",
    x = "Calibration size (n)",
    color = "Screening threshold", shape = "Screening threshold"#,
    #title = sprintf("Time: %.1f, Yield width with 95%% bands", t0)
  ) +
    scale_color_manual(values = lambda_colors) +
    scale_shape_manual(values = lambda_shapes) +    
  theme_bw()

ggsave(sprintf("figures/experiments_yield_width_t%.1f.png", t0), p.yi.w, width=9, height=3.5, units="in", bg = 'white')

p.er.w <- summary.plot %>%
  filter(summary.plot$lambda_label !="Worst-case") %>%
  ggplot(aes(x = n_cal, y = mean_width_ppv, color = lambda_label, shape = lambda_label)) +
  geom_point() +
  geom_line() +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean_width_ppv - 1.96 * se_width_ppv),
      ymax = pmin(1, mean_width_ppv + 1.96 * se_width_ppv)
    ),
    width = 0.02,
    alpha = 0.4
  ) +
  facet_grid(band_type ~ band_method + n_test_label) +
  labs(
    y = "PPV width",
    x = "Calibration size (n)",
    color = "Screening threshold", shape = "Screening threshold"#,
    #title = sprintf("Time: %.1f, PPV width with 95%% bands", t0)
  ) +
    scale_color_manual(values = lambda_colors) +
    scale_shape_manual(values = lambda_shapes) +    
  theme_bw()

ggsave(sprintf("figures/experiments_ppv_width_t%.1f.png", t0), p.er.w, width=9, height=3.5, units="in", bg = 'white')
