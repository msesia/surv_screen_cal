finite_sample_pointwise <- function(
  S,
  E,
  O = NULL,
  prob_O = NULL,
  prob_O_min = 1e-6,
  m = NULL,
  alpha = 0.05,
  ratio_if_A0 = NA_real_,
  clip_pE = TRUE,
  method_pE = "empirical_bernstein"
) {
  stopifnot(length(S) == length(E))
  if (!is.null(O))      stopifnot(length(O) == length(S))
  if (!is.null(prob_O)) stopifnot(length(prob_O) == length(S))

  S <- as.integer(S)
  n <- length(S)
  
  # ---- point estimates ----
  point_est <- point_estimate_pS_pE(S, E, O, prob_O, prob_O_min, clip_pE, ipcw_method="ht")
  pS_hat <- point_est$pS
  pE_hat <- point_est$pE

  # ---- finite-sample confidence bounds ----
  confint.pS <- bound_binomial_exact(sum(S == 1L), n, alpha = alpha)
  confint.pE <- bound_pE_ipcw_fs(
    S = S, E = E, O = O, prob_O = prob_O,
    alpha = alpha, method = method_pE
  )
  if (clip_pE) confint.pE <- pmax(0, pmin(1, confint.pE))

  # ---- predictive intervals (optional) ----
  predint <- NULL
  if (!is.null(m)) {
    tmp <- predint_wc_yield_and_event_rate(
      m,
      confint.pS[1], confint.pS[2],
      confint.pE[1], confint.pE[2],
      alpha = alpha,
      ratio_if_A0 = ratio_if_A0
    )
    predint <- tmp$predint
  }

  out <- list(
    estimate = list(pS = pS_hat, pE = pE_hat),
    confint  = list(pS = confint.pS, pE = confint.pE),
    predint  = predint,
    settings = list(
      alpha = alpha,
      m = m,
      ratio_if_A0 = ratio_if_A0,
      clip_pE = clip_pE,
      prob_O_min = prob_O_min,
      ipcw_used = !is.null(O) || !is.null(prob_O),
      finite_sample = TRUE,
      method_pE = method_pE
    )
  )

  class(out) <- "selection_event_bounds"
  out
}

#' Exact (Clopper–Pearson) binomial confidence bounds
#'
#' Computes finite-sample confidence bounds for p in Binomial(n, p),
#' using the Clopper–Pearson (Beta inversion) method.
#'
#' Inputs:
#'   k : number of successes (integer, 0 <= k <= n)
#'   n : number of trials (integer, n >= 1)
#'
#' Interval tails:
#'   alpha_lower and alpha_upper control lower/upper tail probabilities.
#'   Either may be NA for one-sided bounds.
#'
#' Defaults:
#'   alpha = 0.05
#'   alpha_lower = alpha/2
#'   alpha_upper = alpha/2
#'
#' Output:
#'   Named numeric vector c(lo = ..., hi = ...)
#'
bound_binomial_exact <- function(k, n, alpha = 0.05, alpha_lower = alpha / 2, alpha_upper = alpha / 2) {
  if (n <= 0L) stop("n must be positive.")
  if (k < 0L || k > n) stop("k must satisfy 0 <= k <= n.")
  # ---- lower bound ----
  lo <- if (is.na(alpha_lower)) {
    NA_real_
  } else if (k == 0L) {
    0
  } else {
    stats::qbeta(alpha_lower, k, n - k + 1)
  }
  # ---- upper bound ----
  hi <- if (is.na(alpha_upper)) {
    NA_real_
  } else if (k == n) {
    1
  } else {
    stats::qbeta(1 - alpha_upper, k + 1, n - k)
  }
  c(lo = lo, hi = hi)
}

#' Finite-sample bounds for pE = P(E=1 | S=1) using IPCW + concentration
#'
#' Goal:
#'   Bound the event risk among selected:
#'     pE = P(E=1 | S=1)
#'   when event status at the horizon may be censored.
#'
#' Inputs (all length n):
#'   S      : {0,1} selection indicator. S[i]=1 means subject i is selected.
#'   E      : {0,1,NA} event-by-t indicator. May be NA when O[i]=0.
#'   O      : {0,1} observed-at-t indicator. O[i]=1 means E[i] is observed.
#'   prob_O : estimated P(O=1 | Xi) = P(C >= t | Xi), i.e., probability outcome is observed at t.
#'
#' IPCW contribution:
#'   Define Y[i] = 0 if O[i]=0, else Y[i] = E[i] / max(prob_O[i], prob_O_min).
#'   Then (under standard IPCW conditions) pE = E[Y | S=1] and
#'     pE = phi / mu where:
#'       mu  = E[S]
#'       phi = E[S * Y]
#'
#' Bound construction:
#'   - Bound mu = E[S] using exact binomial (Clopper–Pearson) bounds on k/n, k=sum(S).
#'   - Bound phi = E[S*Y] using either empirical Bernstein or Hoeffding for bounded Z=S*Y in [0,M].
#'   - Combine via ratio bounds:
#'       lower:  max(0, phi_low / mu_upp)
#'       upper:  min(1, phi_upp / mu_low)
#'
#' IMPORTANT about M (boundedness):
#'   Concentration bounds require Z = S*Y in [0, M] almost surely.
#'   For strict finite-sample validity, pass a deterministic M (e.g., from weight clipping).
#'   If M is NULL, this function uses max(Z) from the data as a pragmatic fallback.
#'
#' Interval tails:
#'   Controlled by alpha_lower and alpha_upper. Either may be NA for one-sided bounds.
#'
#' Returns:
#'   Named numeric vector c(lo = ..., hi = ...). NA on a side means "not requested".
bound_pE_ipcw_fs <- function(
  S, E, O, prob_O,
  alpha = 0.05,
  alpha_lower = alpha / 2,
  alpha_upper = alpha / 2,
  prob_O_min = 1e-6,
  M = NULL,
  method = c("empirical_bernstein", "hoeffding")
) {
  method <- match.arg(method)

  stopifnot(length(S) == length(E),
            length(S) == length(O),
            length(S) == length(prob_O))

  n <- length(S)
  S <- as.integer(S)
  O <- as.integer(O)
  prob_O <- as.numeric(prob_O)

  # --- input checks (match your NA logic) ---
  if (any(is.na(O))) stop("O contains NA.")
  if (any(!(O %in% c(0L, 1L)))) stop("O must be in {0,1}.")
  if (any(is.na(E) & O == 1L)) stop("E has NA where O==1; E must be observed when O==1.")
  if (any(!(is.na(E) | E %in% c(0, 1)))) stop("E must be 0/1/NA.")
  if (any(is.na(prob_O[O == 1L]))) stop("prob_O has NA where O==1.")
  if (any(prob_O[O == 1L] < 0)) stop("prob_O must be >= 0 where O==1.")

  # IPCW contribution Y: 0 if O==0; else E/prob_O
  denom <- pmax(prob_O, prob_O_min)
  Y <- numeric(n)
  idx_obs <- which(O == 1L)
  Y[idx_obs] <- as.numeric(E[idx_obs]) / denom[idx_obs]

  # Z = S * Y, so phi = E[Z]
  Z <- as.numeric(S) * Y
  phi_hat <- mean(Z)
  varZ <- stats::var(Z)

  # Need Z in [0, M]
  if (is.null(M)) {
    # Pragmatic fallback; for strict finite-sample validity, pass deterministic M.
    M <- max(Z, na.rm = TRUE)
  }
  if (!is.finite(M) || M <= 0) {
    # If M is degenerate, we can't form a meaningful concentration bound
    return(c(lo = if (is.na(alpha_lower)) NA_real_ else 0,
             hi = if (is.na(alpha_upper)) NA_real_ else 1))
  }

  # ---- bounds for mu = E[S] using exact binomial ----
  k <- sum(S == 1L)

  # allocate tail budgets:
  # - for upper bound on pE we need mu_low and phi_upp
  # - for lower bound on pE we need mu_upp and phi_low
  #
  # Use a simple split within each side: half to mu and half to phi.
  mu_low <- mu_upp <- NA_real_
  phi_low <- phi_upp <- NA_real_

  # ---- upper bound side ----
  if (!is.na(alpha_upper)) {
    mu_bounds_u <- bound_binomial_exact(k, n, alpha_lower = NA, alpha_upper = alpha_upper / 2)
    mu_low <- mu_bounds_u["lo"]  # NA by construction
    mu_low <- bound_binomial_exact(k, n, alpha_lower = alpha_upper / 2, alpha_upper = NA)["lo"]

    # phi upper bound with tail alpha_upper/2
    if (method == "empirical_bernstein") {
      phi_upp <- phi_hat +
        sqrt(2 * varZ * log(2 / (alpha_upper / 2)) / n) +
        (7 * M * log(2 / (alpha_upper / 2))) / (3 * max(1, n - 1))
    } else {
      phi_upp <- phi_hat + M * sqrt(log(2 / (alpha_upper / 2)) / (2 * n))
    }
  }

  # ---- lower bound side ----
  if (!is.na(alpha_lower)) {
    mu_upp <- bound_binomial_exact(k, n,
                                  alpha_lower = NA, alpha_upper = alpha_lower / 2)["hi"]

    # phi lower bound with tail alpha_lower/2
    if (method == "empirical_bernstein") {
      rad <- sqrt(2 * varZ * log(2 / (alpha_lower / 2)) / n) +
        (7 * M * log(2 / (alpha_lower / 2))) / (3 * max(1, n - 1))
      phi_low <- phi_hat - rad
    } else {
      phi_low <- phi_hat - M * sqrt(log(2 / (alpha_lower / 2)) / (2 * n))
    }
  }

  # Combine ratio bounds
  lo <- if (is.na(alpha_lower)) {
    NA_real_
  } else {
    max(0, phi_low / mu_upp)
  }

  hi <- if (is.na(alpha_upper)) {
    NA_real_
  } else {
    # if mu_low is tiny/nonpositive, upper bound becomes 1
    if (!is.finite(mu_low) || mu_low <= 0) 1 else min(1, phi_upp / mu_low)
  }

  # clip
  lo <- if (!is.na(lo)) max(0, min(1, lo)) else lo
  hi <- if (!is.na(hi)) max(0, min(1, hi)) else hi

  c(lo = lo, hi = hi)
}


#' Worst-case prediction intervals for a future cohort given bounds on (pS, pE)
#'
#' This is a worst-case analogue of the predictive part of
#' bootstrap_bounds_selection_and_event_rate(), but instead of using bootstrap
#' draws for (pS, pE), it uses user-supplied confidence bounds:
#'   pS in [pS_lo, pS_hi] where pS = P(S=1)
#'   pE in [pE_lo, pE_hi] where pE = P(E=1 | S=1)
#'
#' Future-cohort model (same as in bootstrap function):
#'   A ~ Binomial(m, pS)                 # number selected
#'   B | A=a ~ Binomial(a, pE)           # number of events among selected
#'
#' The returned prediction intervals are conservative ("worst-case") over the
#' parameter rectangles [pS_lo,pS_hi] and [pE_lo,pE_hi].
#'
#' Interval tails:
#'   alpha_lower and alpha_upper behave like in the bootstrap function; either
#'   can be NA for one-sided bounds.
#'
#' Output naming matches the bootstrap function:
#'   - prop_selected corresponds to A/m
#'   - prop_event_among_selected corresponds to B/A (ratio_if_A0 when A=0)
#'
predint_wc_yield <- function(m, pS_lo, pS_hi, alpha = 0.05, alpha_lower = alpha/2, alpha_upper = alpha/2) {
  # Worst-case binomial prediction interval for A over pS in [pS_lo, pS_hi]
  A_lo <- if (is.na(alpha_lower)) NA_integer_ else stats::qbinom(alpha_lower, m, pS_lo)
  A_hi <- if (is.na(alpha_upper)) NA_integer_ else stats::qbinom(1 - alpha_upper, m, pS_hi)
  c(lo = A_lo/m, hi = A_hi/m)
}


predint_wc_yield_and_event_rate <- function(
  m,
  pS_lo, pS_hi,          # bounds for pS = P(S=1)
  pE_lo, pE_hi,          # bounds for pE = P(E=1 | S=1)
  alpha = 0.05,
  alpha_lower = alpha/2, alpha_upper = alpha/2,
  split = 0.5,           # fraction of tail probability allocated to predicting A
  ratio_if_A0 = 0
) {
  # ---- allocate tail probabilities (Bonferroni split) ----
  # We want joint coverage for:
  #   (i) A in its prediction interval, and
  #   (ii) B/A in its prediction interval (given A).
  # Use a conservative split of tail error across the two steps.
  alphaA_L <- if (is.na(alpha_lower)) NA else split * alpha_lower
  alphaA_U <- if (is.na(alpha_upper)) NA else split * alpha_upper
  alphaB_L <- if (is.na(alpha_lower)) NA else (1 - split) * alpha_lower
  alphaB_U <- if (is.na(alpha_upper)) NA else (1 - split) * alpha_upper

  # ---- Step 1: predict A (selected count) under worst-case pS in [pS_lo,pS_hi] ----
  A_lo <- if (is.na(alphaA_L)) 0L else stats::qbinom(alphaA_L, m, pS_lo)
  A_hi <- if (is.na(alphaA_U)) m  else stats::qbinom(1 - alphaA_U, m, pS_hi)

  # Consider all plausible selected counts in the predicted range
  A_grid <- A_lo:A_hi

  # ---- Step 2: for each plausible A=a, predict B|A=a under worst-case pE bounds ----
  # We then map B bounds to bounds on the ratio B/A and take the union over a.
  rate_lo <- +Inf
  rate_hi <- -Inf

  for (a in A_grid) {
    if (a == 0L) {
      # Convention when no one is selected
      rate_lo <- min(rate_lo, ratio_if_A0)
      rate_hi <- max(rate_hi, ratio_if_A0)
      next
    }

    # Worst-case binomial prediction interval for B|A=a over pE in [pE_lo,pE_hi]
    B_lo <- if (is.na(alphaB_L)) 0L else stats::qbinom(alphaB_L, a, pE_lo)
    B_hi <- if (is.na(alphaB_U)) a  else stats::qbinom(1 - alphaB_U, a, pE_hi)

    # Convert to event rate among selected
    rate_lo <- min(rate_lo, B_lo / a)
    rate_hi <- max(rate_hi, B_hi / a)
  }

    ## ---- finalize bounds (replace NA with [0,1]) ----
    prop_sel_lo <- A_lo / m
    prop_sel_hi <- A_hi / m

    if (!is.finite(prop_sel_lo) || !is.finite(prop_sel_hi)) {
        prop_sel_lo <- 0
        prop_sel_hi <- 1
        
    }

    if (!is.finite(rate_lo) || !is.finite(rate_hi)) {
        rate_lo <- 0
        rate_hi <- 1
    }

    list(
        predint = list(
            m = m,
            prop_selected = c(lo = prop_sel_lo, hi = prop_sel_hi),
            prop_event_among_selected = c(lo = rate_lo, hi = rate_hi)
        ),
        settings = list(
            alpha = alpha,
            alpha_lower = alpha_lower,
            alpha_upper = alpha_upper,
            split = split
        )
    )
}

