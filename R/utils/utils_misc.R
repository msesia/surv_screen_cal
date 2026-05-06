## Helper function
make_lambda_row <- function(lambda, S, estim, E, O, prob_O, test_cohort_size) {
  out <- tibble(
    lambda = lambda,
    num_cal = length(S),
    num_sel = sum(S),
    ## point estimates
    pS_hat = estim$estimate$pS,
    pE_hat = estim$estimate$pE,
    ## confidence bands
    pS_lo  = estim$confint$pS[1],
    pS_hi  = estim$confint$pS[2],
    pE_lo  = estim$confint$pE[1],
    pE_hi  = estim$confint$pE[2]
  )
  if (!is.null(test_cohort_size) && !is.null(estim$predint)) {
    out <- out |>
      mutate(
        num_test = estim$predint$m,
        pS_pred_lo = estim$predint$prop_selected[1],
        pS_pred_hi = estim$predint$prop_selected[2],
        pE_pred_lo = estim$predint$prop_event_among_selected[1],
        pE_pred_hi = estim$predint$prop_event_among_selected[2]
      )
  } else {
    out <- out |>
      mutate(
        num_test = NA_real_,
        pS_pred_lo = NA_real_,
        pS_pred_hi = NA_real_,
        pE_pred_lo = NA_real_,
        pE_pred_hi = NA_real_
      )
  }      
  out <- out %>% mutate(yield_hat = pS_hat, ppv_hat=1-pE_hat, yield_pop_lo = pS_lo, yield_pop_hi = pS_hi,
                        ppv_pop_lo = 1-pE_hi, ppv_pop_hi = 1-pE_lo,
                        selected_test_lo = floor(num_test * pS_pred_lo), selected_test_hi = ceiling(num_test * pS_pred_hi),
                        ppv_test_lo = 1-pE_pred_hi, ppv_test_hi = 1-pE_pred_lo) %>%
      select(-pS_hat, -pE_hat, -pS_lo, -pS_hi, -pE_lo, -pE_hi, -pS_pred_lo, -pS_pred_hi, -pE_pred_lo, -pE_pred_hi) %>%
      select(lambda, num_test, everything())
  out
}

convert_band_to_long_targets <- function(estim, band_method, m = NULL) {
  ## Confidence bands (population / calibration uncertainty)
  conf_df <- tibble::tibble(
    lambda = estim$estimate$lambda,
    band_method = band_method,
    target_m = Inf,
    band_type = "confidence",
    yield_lo = estim$confint_simul$pS[, 1],
    yield_hi = estim$confint_simul$pS[, 2],
    ppv_lo    = estim$confint_simul$pE[, 1],
    ppv_hi    = estim$confint_simul$pE[, 2]
  )
  ## Prediction bands (future cohort), only if present
  pred_df <- NULL
  if (!is.null(m) &&
      !is.null(estim$predint_simul) &&
      !is.null(estim$predint_simul$prop_selected)) {
    pred_df <- tibble::tibble(
      lambda = estim$predint_simul$lambda,
      band_method = band_method,
      target_m = m,
      band_type = "prediction",
      yield_lo = estim$predint_simul$prop_selected[, 1],
      yield_hi = estim$predint_simul$prop_selected[, 2],
      ppv_lo    = if (!is.null(estim$predint_simul$prop_event_among_selected))
                      estim$predint_simul$prop_event_among_selected[, 1] else NA_real_,
      ppv_hi    = if (!is.null(estim$predint_simul$prop_event_among_selected))
                      estim$predint_simul$prop_event_among_selected[, 2] else NA_real_
    )
  }
  dplyr::bind_rows(conf_df, pred_df)
}

convert_pointwise_to_long_targets <- function(band_df, band_method, m = NULL) {

  ## ---- Confidence bands ----
  conf_df <- band_df %>%
    dplyr::transmute(
      lambda,
      band_method = band_method,
      target_m = Inf,
      band_type = "confidence",
      yield_lo = yield_pop_lo,
      yield_hi = yield_pop_hi,
      ppv_lo    = ppv_pop_lo,
      ppv_hi    = ppv_pop_hi
    )

  ## ---- Prediction bands (if available) ----
  pred_df <- NULL
  if (!is.null(m) &&
      all(c("yield_pred_lo","yield_pred_hi") %in% names(band_df))) {

    pred_df <- band_df %>%
      dplyr::transmute(
        lambda,
        band_method = band_method,
        target_m = m,
        band_type = "prediction",
        yield_lo = yield_pred_lo,
        yield_hi = yield_pred_hi,
        ppv_lo    = if ("ppv_pred_lo" %in% names(band_df)) ppv_pred_lo else NA_real_,
        ppv_hi    = if ("ppv_pred_hi" %in% names(band_df)) ppv_pred_hi else NA_real_
      )
  }

  dplyr::bind_rows(conf_df, pred_df)
}

add_worstcase_prediction_bands <- function(df,
                                          alpha = 0.05,
                                          ppv_grid = 5,
                                          yield_grid = 5,
                                          mc_reps = 1000,
                                          seed = 1,
                                          ppv_when_S0 = c("NA", "01")) {
  stopifnot(is.data.frame(df))

  # Ensure required columns are present
  req <- c("num_test", "yield_pop_lo", "yield_pop_hi", "ppv_pop_lo", "ppv_pop_hi")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

  ppv_when_S0 <- match.arg(ppv_when_S0)
  a2 <- alpha / 2

  # Helper: simulate PPV_hat = T/S under (pS, pPPV)
  sim_ppv_hat <- function(m, pS, pPPV, R) {
    S <- stats::rbinom(R, size = m, prob = pS)
    # If S=0, PPV is undefined; we return NA and handle later.
    T <- stats::rbinom(R, size = S, prob = pPPV)
    ppv_hat <- ifelse(S > 0, T / S, NA_real_)
    ppv_hat
  }

  set.seed(seed)

  df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      # ------------------------------------------------------------
      # Yield prediction band (counts) from pop yield CI via binomial
      # S ~ Bin(m, pS), with pS in [yield_pop_lo, yield_pop_hi]
      # We take worst-case tails by using the extreme pS values.
      # ------------------------------------------------------------
      selected_test_lo = stats::qbinom(a2,       size = num_test, prob = yield_pop_lo),
      selected_test_hi = stats::qbinom(1 - a2,   size = num_test, prob = yield_pop_hi),

      # Optional: yield prediction band on probability scale
      yield_test_lo = selected_test_lo / num_test,
      yield_test_hi = selected_test_hi / num_test,

      # ------------------------------------------------------------
      # PPV prediction band:
      # PPV_hat = T/S where
      #   S ~ Bin(m, pS)
      #   T | S ~ Bin(S, pPPV)
      #
      # We have two pop CIs: pS in [..] and pPPV in [..]
      # We build a conservative band by taking worst-case quantiles
      # over a grid of (pS, pPPV) values spanning that rectangle.
      # ------------------------------------------------------------
      ppv_test_lo = {
        pS_vals   <- seq(yield_pop_lo, yield_pop_hi, length.out = yield_grid)
        pPPV_vals <- seq(ppv_pop_lo,   ppv_pop_hi,   length.out = ppv_grid)

        # For each parameter pair, compute lower quantile of PPV_hat
        q_los <- c()
        for (pS in pS_vals) for (pPPV in pPPV_vals) {
          sims <- sim_ppv_hat(num_test, pS, pPPV, mc_reps)
          sims <- sims[!is.na(sims)]  # drop S=0 cases
          if (length(sims) == 0) {
            # If almost always S=0, PPV is effectively undefined
            q_los <- c(q_los, NA_real_)
          } else {
            q_los <- c(q_los, stats::quantile(sims, probs = a2, names = FALSE, type = 1))
          }
        }
        # Worst-case lower bound: smallest quantile across parameter rectangle
        out <- suppressWarnings(min(q_los, na.rm = TRUE))
        if (is.infinite(out)) out <- NA_real_
        out
      },

      ppv_test_hi = {
        pS_vals   <- seq(yield_pop_lo, yield_pop_hi, length.out = yield_grid)
        pPPV_vals <- seq(ppv_pop_lo,   ppv_pop_hi,   length.out = ppv_grid)

        # For each parameter pair, compute upper quantile of PPV_hat
        q_his <- c()
        for (pS in pS_vals) for (pPPV in pPPV_vals) {
          sims <- sim_ppv_hat(num_test, pS, pPPV, mc_reps)
          sims <- sims[!is.na(sims)]
          if (length(sims) == 0) {
            q_his <- c(q_his, NA_real_)
          } else {
            q_his <- c(q_his, stats::quantile(sims, probs = 1 - a2, names = FALSE, type = 1))
          }
        }
        # Worst-case upper bound: largest quantile across parameter rectangle
        out <- suppressWarnings(max(q_his, na.rm = TRUE))
        if (is.infinite(out)) out <- NA_real_
        out
      },

      # ------------------------------------------------------------
      # Optional convention when S=0 can happen.
      # If you prefer a deterministic fallback band when PPV is undefined,
      # set ppv_when_S0="01" and we widen to [0,1] whenever yield is tiny.
      # Here we use a crude trigger: selected_test_hi == 0.
      # ------------------------------------------------------------
      ppv_test_lo = dplyr::if_else(ppv_when_S0 == "01" & selected_test_hi == 0, 0, ppv_test_lo),
      ppv_test_hi = dplyr::if_else(ppv_when_S0 == "01" & selected_test_hi == 0, 1, ppv_test_hi)
    ) |>
    dplyr::ungroup()
}
