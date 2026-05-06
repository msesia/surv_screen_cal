#' Bootstrap intervals for selection rate and event risk among selected
#'
#' Wrapper around `bootstrap_simultaneous_internal()` that optionally
#' applies "FS smoothing" by appending one pseudo-observation and forming an envelope CI:
#'   - pessimistic run uses pseudo E*=0  (lower bounds)
#'   - optimistic run uses pseudo E*=1  (upper bounds)
#'
#' Returns
#'   - confint: lower from pessimistic, upper from optimistic
#'   - predint: same envelope logic (if `m` is provided)
#'
#'
#' @param S Integer vector in {0,1}.
#'   Selection indicator. `S[i] = 1` means subject i is selected.
#'
#' @param E Numeric/integer vector in {0,1,NA}.
#'   Event indicator by time t.
#'   - 1 = event occurred
#'   - 0 = no event
#'   - NA allowed only when using IPCW (i.e., when `O[i] = 0`)
#'
#' @param O Optional integer vector in {0,1}.
#'   Observation indicator for event status at time t.
#'   `O[i] = 1` means event status is observed.
#'   Required if `prob_O` is supplied.
#'
#' @param prob_O Optional numeric vector in [0,1].
#'   Estimated probability that event status is observed
#'   (e.g., `P(C >= t | X_i)`).
#'   Required if `O` is supplied.
#'
#' @param prob_O_min Numeric scalar > 0.
#'   Lower bound applied to `prob_O` to stabilize IPCW weights.
#'
#' @param m Optional positive integer.
#'   Size of a future cohort for predictive intervals.
#'
#' @param alpha Numeric scalar in (0,1).
#'   Total error level for two-sided intervals.
#'
#' @param B Positive integer.
#'   Number of bootstrap resamples.
#'
#' @param seed Integer.
#'   Random seed for reproducibility.
#'
#' @param clip_pE Logical.
#'   If TRUE, constrains bootstrap draws of `pE` to [0,1].
#'
#' @param fs_smoothing Logical.
#'   If TRUE, applies finite-sample smoothing via a pseudo-observation and
#'   returns envelope confidence bounds.
#'
#' @param lambda_seq Numeric vector.
#'   Sequence of threshold values \eqn{\lambda} at which
#'   selection indicators \eqn{S_\lambda = 1(scores \ge \lambda)}
#'   are evaluated.
#'   Simultaneous confidence and prediction bands are constructed
#'   uniformly over this grid. If \code{NULL}, a default evenly spaced
#'   grid of length 100 over \eqn{[0,1]} is used.
#'
#' @param ipcw_method Character.
#'   Specifies which inverse-probability–of–censoring weighted (IPCW)
#'   estimator to use for the event risk among selected individuals when
#'   \code{O} and \code{prob_O} are supplied.
#'   \describe{
#'     \item{"ht"}{Horvitz–Thompson–style estimator:
#'       mean of \eqn{E_i / \max(prob_O[i], prob_O_min)} among selected,
#'       with zero contribution for unobserved outcomes.}
#'     \item{"hajek"}{Hájek (ratio / normalized-weight) estimator:
#'       \eqn{\sum w_i E_i / \sum w_i}, where \eqn{w_i = 1 / \max(prob_O[i], prob_O_min)}
#'       among selected and observed individuals.}
#'   }
#'   Ignored if IPCW inputs (\code{O}, \code{prob_O}) are not provided.
#'
#' @return A list with components:
#'
#' \describe{
#'
#'   \item{settings}{
#'     List of input settings used to compute the results
#'     (e.g., `alpha`, `B`, `seed`, `ipcw_used`, `fs_smoothing`).
#'   }
#'
#'   \item{estimate}{
#'     Point estimates computed from the observed data:
#'     \itemize{
#'       \item `pS`: estimated selection proportion.
#'       \item `pE`: estimated event risk among selected.
#'     }
#'   }
#'
#' #' \item{confint}{
#'     Percentile confidence intervals:
#'     \itemize{
#'       \item `pS`: bounds for the selection proportion.
#'       \item `pE`: bounds for the event risk among selected.
#'     }
#'     If `fs_smoothing = TRUE`, lower bounds come from the pessimistic run
#'     and upper bounds from the optimistic run.
#'   }
#'
#'   \item{predint}{
#'     (Only if `m` is provided.)
#'     Percentile predictive intervals:
#'     \itemize{
#'       \item `prop_selected`: interval for A/m.
#'       \item `prop_event_among_selected`: interval for B/A.
#'     }
#'     Envelope logic is applied when `fs_smoothing = TRUE`.
#'   }
#'
#' }
#'
bootstrap_simultaneous <- function(
  scores,
  E,
  O = NULL,
  prob_O = NULL,
  prob_O_min = 1e-6,
  m = NULL,
  alpha = 0.05,
  B = 5000,
  seed = 1,
  clip_pE = TRUE,
  fs_smoothing = TRUE,
  ipcw_method = c("ht", "hajek"),
  lambda_seq = NULL,
  K_min = 0
) {
  stopifnot(length(scores) == length(E))
  if (!is.null(O))      stopifnot(length(O) == length(scores))
  if (!is.null(prob_O)) stopifnot(length(prob_O) == length(scores))

  ## Merge Lx2 bands
  envelope <- function(x, y) {
      stopifnot(
          is.matrix(x), is.matrix(y),
          ncol(x) == 2, ncol(y) == 2,
          nrow(x) == nrow(y)
      )
      cbind(
          lo = pmin(x[, 1], y[, 1]),
          hi = pmax(x[, 2], y[, 2])
      )
  }

  run_internal <- function(scores, E, O, prob_O, prob_O_min) {
    bootstrap_simultaneous_internal(
      scores = scores, E = E, O = O, prob_O = prob_O, prob_O_min = prob_O_min,
      m = m,
      alpha = alpha,
      B = B, seed = seed,
      clip_pE = clip_pE,
      ipcw_method = ipcw_method,
      lambda_seq = lambda_seq, K_min = K_min
    )
  }

  # Baseline (no smoothing) run, used for point estimates and returned settings
  base <- run_internal(scores, E, O, prob_O, prob_O_min)

  if (!isTRUE(fs_smoothing)) {
    return(base)
  }

  # ---- FS smoothing: append pseudo observation ----
  prob_O_min2 <- min(1, prob_O_min)

  estim_pess <- run_internal(c(scores, -Inf, Inf), c(E,0,0), c(O,1,1), c(prob_O,1,1), prob_O_min2)
  estim_opt  <- run_internal(c(scores, Inf), c(E, 1), c(O,1), c(prob_O,1), prob_O_min2)

  # ---- Merge envelope bands ----
  out <- base

  out$confint_simul$pS <- envelope(estim_pess$confint_simul$pS, estim_opt$confint_simul$pS)
  out$confint_simul$pE <- envelope(estim_pess$confint_simul$pE, estim_opt$confint_simul$pE)

  # ---- Merge prediction envelope if present ----
  if (!is.null(m) &&
      !is.null(estim_pess$predint_simul) &&
      !is.null(estim_opt$predint_simul)) {

    out$predint_simul <- base$predint_simul

    out$predint_simul$prop_selected <-
      envelope(estim_pess$predint_simul$prop_selected, estim_opt$predint_simul$prop_selected)

    out$predint_simul$prop_event_among_selected <-
      envelope(estim_pess$predint_simul$prop_event_among_selected, estim_opt$predint_simul$prop_event_among_selected)

    # keep lambda in predint_simul
    out$predint_simul$lambda <- base$estimate$lambda
  }

  # ---- Record smoothing metadata ----
  out$settings$fs_smoothing <- TRUE
  out$settings$fs_smoothing_detail <- "Envelope bands: lower from pessimistic (E*=0), upper from optimistic (E*=1)"
  out$settings$seed <- seed
  out$settings$B <- B

  class(out) <- "selection_event_bands"
  out
}

#' Internal worker: bootstrap draws + percentile CIs (single run)
#'
#' This function computes bootstrap draws of:
#'   - pS = mean(S==1)
#'   - pE = mean(contribution | S==1), where contribution depends on IPCW usage
#'
#' IPCW mode:
#'   contribution[i] = E[i] / max(prob_O[i], prob_O_min)  when O[i]==1
#'   contribution[i] = 0                                  when O[i]==0
#'
#' Notes on missingness:
#'   - In IPCW mode, E may be NA only where O==0, and those rows contribute 0.
#'   - Bootstrap draws of pE are set to NA when a resample contains no selected (K==0).
#'   - Confidence intervals ignore NA draws; if all draws are NA, bounds are NA.
bootstrap_simultaneous_internal <- function(
  scores,
  E,
  O = NULL,
  prob_O = NULL,
  prob_O_min = 1e-6,
  m = NULL,
  alpha = 0.05,
  B = 5000,
  seed = 1,
  clip_pE = TRUE,
  ipcw_method = c("ht", "hajek"),
  lambda_seq = NULL,
  lambda_min_denom = 1e-8,
  min_nonNA_frac = 0.9,
  K_min = 10,
  enforce_PI_contains_CI = TRUE
) {
  stopifnot(length(scores) == length(E))
  n <- length(scores)
  ipcw_method <- match.arg(ipcw_method)

  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0,1).")
  }

  if (is.null(lambda_seq)) lambda_seq <- seq(0, 1, length.out = 10)
  L <- length(lambda_seq)

  # ---- IPCW preprocessing ----
  use_ipcw <- !is.null(O) || !is.null(prob_O)
  if (use_ipcw) {
    if (is.null(O) || is.null(prob_O)) stop("If using IPCW, both O and prob_O must be supplied.")
    stopifnot(length(O) == n, length(prob_O) == n)

    O <- as.integer(O)
    prob_O <- as.numeric(prob_O)

    if (any(is.na(O))) stop("O contains NA.")
    if (any(!(O %in% c(0L, 1L)))) stop("O must be in {0,1}.")
    if (any(is.na(E) & O == 1L)) stop("E has NA where O==1; E must be observed when O==1.")
    if (any(!(is.na(E) | E %in% c(0, 1)))) stop("E must be 0/1/NA.")
    if (any(is.na(prob_O[O == 1L]))) stop("prob_O has NA where O==1.")
    if (any(prob_O[O == 1L] < 0)) stop("prob_O must be >= 0 where O==1.")

    denom <- pmax(prob_O, prob_O_min)
    w_all <- 1 / denom

    # HT contribution for E if needed
    contrib_ht <- numeric(n)
    idx_obs <- which(O == 1L)
    contrib_ht[idx_obs] <- as.numeric(E[idx_obs]) / denom[idx_obs]

  } else {
    if (any(is.na(E))) stop("E contains NA but IPCW inputs (O, prob_O) were not provided.")
    if (any(!(E %in% c(0, 1)))) stop("E must be in {0,1}.")
    denom <- NULL
    w_all <- NULL
    O <- NULL
    contrib_ht <- NULL
  }

  # ---- Build S matrix: S_{i,ℓ} = 1(scores_i >= lambda_ℓ) ----
  S_num <- 1.0 * outer(scores, lambda_seq, FUN = ">=")   # n x L

  # ---- Point estimates ----
  pS_hat <- colMeans(S_num)

  # Build numerator/denominator contributions for ER in matrix form
  if (!use_ipcw) {
    E_num <- as.numeric(E)
    Z1_mat <- S_num * E_num   # n x L
    Z0_mat <- S_num           # n x L
  } else if (ipcw_method == "ht") {
    Z1_mat <- S_num * contrib_ht
    Z0_mat <- S_num
  } else { # hajek
    Ew <- numeric(n); Ow <- numeric(n)
    idx_obs <- which(O == 1L)
    Ew[idx_obs] <- as.numeric(E[idx_obs]) * w_all[idx_obs]
    Ow[idx_obs] <- w_all[idx_obs]
    Z1_mat <- S_num * Ew
    Z0_mat <- S_num * Ow
  }

  mu1 <- colMeans(Z1_mat)
  mu0 <- colMeans(Z0_mat)

  mu0_cut <- max(lambda_min_denom, K_min / n)
  ok_E_strict <- is.finite(mu0) & (mu0 >= mu0_cut)

  pE_hat <- rep(NA_real_, L)
  pE_hat[ok_E_strict] <- mu1[ok_E_strict] / mu0[ok_E_strict]
  if (clip_pE) pE_hat <- pmin(1, pmax(0, pE_hat))

  # ---- Influence functions ----
  psiS <- sweep(S_num, 2, pS_hat, FUN = "-")   # n x L

  psiE <- matrix(0, nrow = n, ncol = L)
  if (any(ok_E_strict)) {
    Z1c <- sweep(Z1_mat[, ok_E_strict, drop = FALSE], 2, mu1[ok_E_strict], FUN = "-")
    Z0c <- sweep(Z0_mat[, ok_E_strict, drop = FALSE], 2, mu0[ok_E_strict], FUN = "-")

    a <- 1 / mu0[ok_E_strict]
    b <- mu1[ok_E_strict] / (mu0[ok_E_strict]^2)

    psiE[, ok_E_strict] <- sweep(Z1c, 2, a, FUN = "*") - sweep(Z0c, 2, b, FUN = "*")
  }

  # ---- Standard errors for studentization ----
  sigmaS <- sqrt(colMeans(psiS^2))
  sigmaE <- rep(NA_real_, L)
  if (any(ok_E_strict)) sigmaE[ok_E_strict] <- sqrt(colMeans(psiE[, ok_E_strict, drop = FALSE]^2))

  okS  <- is.finite(sigmaS) & (sigmaS > 0)
  okE2 <- ok_E_strict & is.finite(sigmaE) & (sigmaE > 0)

  # ---- Multiplier draws (separate calibration) ----
  set.seed(seed)
  G <- matrix(stats::rnorm(B * n), nrow = B, ncol = n)

  ZS <- (G %*% psiS) / sqrt(n)  # B x L
  ZE <- matrix(NA_real_, nrow = B, ncol = L)
  if (any(okE2)) ZE[, okE2] <- (G %*% psiE[, okE2, drop = FALSE]) / sqrt(n)

  # two-sided sup stats: max over lambda of |Z|
  TS_abs <- rep(0, B)
  if (any(okS)) {
    ZS_std <- sweep(ZS[, okS, drop = FALSE], 2, sigmaS[okS], FUN = "/")
    TS_abs <- apply(abs(ZS_std), 1, max, na.rm = TRUE)
  }

  TE_abs <- rep(0, B)
  if (any(okE2)) {
    ZE_std <- sweep(ZE[, okE2, drop = FALSE], 2, sigmaE[okE2], FUN = "/")
    TE_abs <- apply(abs(ZE_std), 1, max, na.rm = TRUE)
  }

  cS <- if (any(is.finite(TS_abs))) as.numeric(stats::quantile(TS_abs[is.finite(TS_abs)], probs = 1 - alpha, names = FALSE, type = 7)) else NA_real_

  cE <- if (any(is.finite(TE_abs))) as.numeric(stats::quantile(TE_abs[is.finite(TE_abs)],
                                                               probs = 1 - alpha, names = FALSE, type = 7)) else NA_real_

  ## ---- Simultaneous CI bands (separate cS, cE) ----

  ## Start with degenerate CIs equal to point estimates
  pS_lo <- pS_hat
  pS_hi <- pS_hat

  ## Widen only where variance > 0
  if (any(okS) && is.finite(cS)) {
      pS_lo[okS] <- pS_hat[okS] - cS * sigmaS[okS] / sqrt(n)
      pS_hi[okS] <- pS_hat[okS] + cS * sigmaS[okS] / sqrt(n)
  }

  ## Clamp to [0,1]
  pS_lo <- pmax(0, pS_lo)
  pS_hi <- pmin(1, pS_hi)

  # ER CI: strict-lambda CI; non-estimable lambdas get vacuous [0,1]
  pE_lo <- rep(0, L)
  pE_hi <- rep(1, L)
  if (any(okE2) && is.finite(cE)) {
    pE_lo[okE2] <- pE_hat[okE2] - cE * sigmaE[okE2] / sqrt(n)
    pE_hi[okE2] <- pE_hat[okE2] + cE * sigmaE[okE2] / sqrt(n)
    if (clip_pE) {
      pE_lo <- pmax(0, pE_lo)
      pE_hi <- pmin(1, pE_hi)
    }
  }

  out <- list(
    settings = list(
      alpha = alpha,
      B = B,
      seed = seed,
      n = n,
      m = m,
      prob_O_min = prob_O_min,
      ipcw_used = use_ipcw,
      ipcw_method = ipcw_method,
      lambda_seq = lambda_seq,
      lambda_min_denom = lambda_min_denom,
      K_min = K_min,
      mu0_cut = mu0_cut,
      critical_value_yield = cS,
      critical_value_er = cE
    ),
    estimate = list(
      lambda = lambda_seq,
      pS = pS_hat,
      pE = pE_hat
    ),
    confint_simul = list(
      lambda = lambda_seq,
      pS = cbind(lo = pS_lo, hi = pS_hi),
      pE = cbind(lo = pE_lo, hi = pE_hi)
    )
  )

  # ---- Unconditional simultaneous prediction bands from realized predictive draws ----
  if (!is.null(m)) {

    set.seed(seed)

    # Bootstrap resampling indices for calibration sample
    idx_mat <- matrix(sample.int(n, n * B, replace = TRUE), nrow = B, ncol = n)

    # Shared uniforms for selection (dependence across lambda)
    U_sel_mat <- matrix(stats::runif(B * m), nrow = B, ncol = m)
    # Independent uniforms for events
    U_evt_mat <- matrix(stats::runif(B * m), nrow = B, ncol = m)

    # Store realized predictive draws (curves over lambda)
    A_rate_draws <- matrix(NA_real_, nrow = B, ncol = L)   # A/m
    R_draws      <- matrix(NA_real_, nrow = B, ncol = L)   # B/A (NA when A=0)

    for (b in 1:B) {

      idx <- idx_mat[b, ]

      # Bootstrap pS curve
      pS_b <- colMeans(S_num[idx, , drop = FALSE])
      pS_b <- pmin(1, pmax(0, pS_b))

      # Bootstrap pE curve (plug-in on bootstrap sample)
      S_b <- S_num[idx, , drop = FALSE]

      if (!use_ipcw) {
        E_b  <- as.numeric(E)[idx]
        mu1b <- colMeans(S_b * E_b)
        mu0b <- colMeans(S_b)
      } else if (ipcw_method == "ht") {
        contrib_b <- contrib_ht[idx]
        mu1b <- colMeans(S_b * contrib_b)
        mu0b <- colMeans(S_b)
      } else {
        O_b <- O[idx]
        w_b <- 1 / denom[idx]
        Ew_b <- ifelse(O_b == 1L, as.numeric(E)[idx] * w_b, 0)
        Ow_b <- ifelse(O_b == 1L, w_b, 0)
        mu1b <- colMeans(S_b * Ew_b)
        mu0b <- colMeans(S_b * Ow_b)
      }

      pE_b <- rep(NA_real_, L)
      ok0  <- is.finite(mu0b) & (mu0b > 0)
      pE_b[ok0] <- mu1b[ok0] / mu0b[ok0]
      if (clip_pE) pE_b <- pmin(1, pmax(0, pE_b))

      # Shared uniforms for this replicate
      U_sel <- U_sel_mat[b, ]
      U_evt <- U_evt_mat[b, ]

      # Future selection indicators over lambda: m x L
      Sel_mat <- outer(U_sel, pS_b, FUN = "<=")

      # A(lambda) and A/m
      A_cnt  <- colSums(Sel_mat)
      A_rate <- A_cnt / m

      # Future events among selected
      Evt_mat <- outer(U_evt, pE_b, FUN = "<=")
      B_cnt <- colSums(Sel_mat & Evt_mat)

      # Ratio R = B/A (undefined if A=0)
      R_hat <- rep(NA_real_, L)
      okA <- (A_cnt > 0) & is.finite(pE_b)
      R_hat[okA] <- B_cnt[okA] / A_cnt[okA]

      A_rate_draws[b, ] <- A_rate
      R_draws[b, ]      <- R_hat
    }
      
    # Center curves
    center_A <- pS_hat
    center_R <- pE_hat
    
    # Studentization scales: empirical predictive SD at each lambda
    kS_vec <- colSums(S_num)
    pS_lambda_bounds <- t(sapply(kS_vec, function(k) {bound_binomial_exact(k, n, alpha=0.05, alpha_lower=0.05, alpha_upper=0.05)}))    
    sdA_lambda_bounds <- sqrt(pS_lambda_bounds*(1-pS_lambda_bounds)/m)
    sdA_lambda <- apply(sdA_lambda_bounds, 1, max)             

    sdR_lambda <- apply(R_draws,      2, stats::sd, na.rm = TRUE)
           
    centerA_mat <- matrix(center_A, nrow = B, ncol = L, byrow = TRUE)
    centerR_mat <- matrix(center_R, nrow = B, ncol = L, byrow = TRUE)
    sdA_mat     <- matrix(sdA_lambda, nrow = B, ncol = L, byrow = TRUE)
    sdR_mat     <- matrix(sdR_lambda, nrow = B, ncol = L, byrow = TRUE)

    # two-sided |Z| sup stats (calibrated separately)

    # Yield
    ZA_abs <- apply(abs((A_rate_draws - centerA_mat) / sdA_mat), 1, max, na.rm = TRUE)
    qA <- as.numeric(stats::quantile(ZA_abs[is.finite(ZA_abs)], probs = 1 - alpha, names = FALSE, type = 7))

    A_lo <- center_A - qA * sdA_lambda
    A_hi <- center_A + qA * sdA_lambda
    A_lo <- pmax(0, A_lo)
    A_hi <- pmin(1, A_hi)
      
      ## ER
      qR <- NA_real_
      R_lo <- rep(0, L); R_hi <- rep(1, L)

      nonNA_frac_R <- colMeans(is.finite(R_draws))
      okR <- is.finite(center_R) & (nonNA_frac_R >= min_nonNA_frac)
      if (any(okR)) {
          ZR_abs <- apply(abs((R_draws[, okR, drop = FALSE] - centerR_mat[, okR, drop = FALSE]) /
                              sdR_mat[, okR, drop = FALSE]),
                          1, max, na.rm = TRUE)
          qR <- as.numeric(stats::quantile(ZR_abs[is.finite(ZR_abs)],
                                           probs = 1 - alpha, names = FALSE, type = 7))

          R_lo[okR] <- center_R[okR] - qR * sdR_lambda[okR]
          R_hi[okR] <- center_R[okR] + qR * sdR_lambda[okR]
          R_lo <- pmax(0, R_lo); R_hi <- pmin(1, R_hi)
      }

    # Optional: enforce PI contains CI pointwise
    if (enforce_PI_contains_CI) {
      CI_pS_lo <- out$confint_simul$pS[, "lo"]
      CI_pS_hi <- out$confint_simul$pS[, "hi"]
      CI_pE_lo <- out$confint_simul$pE[, "lo"]
      CI_pE_hi <- out$confint_simul$pE[, "hi"]

      okS2 <- is.finite(CI_pS_lo) & is.finite(CI_pS_hi) & is.finite(A_lo) & is.finite(A_hi)
      A_lo[okS2] <- pmin(A_lo[okS2], CI_pS_lo[okS2])
      A_hi[okS2] <- pmax(A_hi[okS2], CI_pS_hi[okS2])

      okE2 <- is.finite(CI_pE_lo) & is.finite(CI_pE_hi)
      R_lo[okE2] <- pmin(R_lo[okE2], CI_pE_lo[okE2])
      R_hi[okE2] <- pmax(R_hi[okE2], CI_pE_hi[okE2])

      A_lo <- pmax(0, A_lo); A_hi <- pmin(1, A_hi)
      R_lo <- pmax(0, R_lo); R_hi <- pmin(1, R_hi)
    }

    out$predint_simul <- list(
      lambda = lambda_seq,
      prop_selected = cbind(lo = A_lo, hi = A_hi),
      prop_event_among_selected = cbind(lo = R_lo, hi = R_hi),
      critical_value_pred_yield = qA,
      critical_value_pred_er = qR,
      center = list(prop_selected = center_A, prop_event_among_selected = center_R),
      scale = list(prop_selected_sd = sdA_lambda, prop_event_among_selected_sd = sdR_lambda),
      method = "Unconditional predictive simultaneous bands from realized draws; two-sided |Z|; calibrated separately for Yield and ER"
    )
  }

  class(out) <- "selection_event_bands"
  out
}
