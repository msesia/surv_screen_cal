# =============================================================================
# Pointwise bootstrap intervals (baseline) for yield and event rate among
# selected patients, at a single threshold.
#
# Mirrors bootstrap_simultaneous() except that intervals are percentile
# intervals calibrated at one threshold (no supremum over a grid).
# Vacuous-interval rules are the same as in the simultaneous method:
#   - n_min : minimum number of selected calibration patients for the event
#             rate to be estimable (otherwise CI and PI for pE are [0,1]).
#   - pi_min: minimum fraction of simulated cohorts with >= 1 selected patient
#             for the predictive interval of pE to be reported.
# =============================================================================


# ---- Wrapper: finite-sample smoothing via envelope of two augmented runs ----

bootstrap_pointwise <- function(S, E,
                                O = NULL, prob_O = NULL, prob_O_min = 1e-6,
                                m = NULL,
                                alpha = 0.05, B = 5000, seed = 1,
                                clip_pE = TRUE,
                                fs_smoothing = TRUE,
                                ipcw_method = c("ht", "hajek"),
                                n_min = 0, pi_min = 0.9) {

  stopifnot(length(S) == length(E))
  if (!is.null(O))      stopifnot(length(O) == length(S))
  if (!is.null(prob_O)) stopifnot(length(prob_O) == length(S))
  ipcw_method <- match.arg(ipcw_method)

  envelope <- function(x, y) c(min(x[1], y[1]), max(x[2], y[2]))

  run <- function(S, E, O, prob_O, prob_O_min) {
    bootstrap_pointwise_internal(
      S = S, E = E, O = O, prob_O = prob_O, prob_O_min = prob_O_min,
      m = m, alpha = alpha, B = B, seed = seed, clip_pE = clip_pE,
      ipcw_method = ipcw_method, n_min = n_min, pi_min = pi_min
    )
  }

  base <- run(S, E, O, prob_O, prob_O_min)
  if (!isTRUE(fs_smoothing)) return(base)

  # Pseudo-observations (uncensored, prob_O = 1):
  #   run A: one non-selected and one selected patient without event
  #   run B: one selected patient with event
  prob_O_min2 <- min(1, prob_O_min)
  run_A <- run(c(S, 0, 1), c(E, 0, 0), c(O, 1, 1), c(prob_O, 1, 1), prob_O_min2)
  run_B <- run(c(S, 1),    c(E, 1),    c(O, 1),    c(prob_O, 1),    prob_O_min2)

  out <- base   # point estimates and settings from the unsmoothed run

  out$confint$pS <- envelope(run_A$confint$pS, run_B$confint$pS)
  out$confint$pE <- envelope(run_A$confint$pE, run_B$confint$pE)

  if (!is.null(m)) {
    out$predint$prop_selected <-
      envelope(run_A$predint$prop_selected, run_B$predint$prop_selected)
    out$predint$prop_event_among_selected <-
      envelope(run_A$predint$prop_event_among_selected,
               run_B$predint$prop_event_among_selected)
  }

  out$settings$fs_smoothing <- TRUE
  class(out) <- "selection_event_bounds"
  out
}


# ---- Internal worker: one bootstrap run, no smoothing ----

bootstrap_pointwise_internal <- function(S, E,
                                         O = NULL, prob_O = NULL, prob_O_min = 1e-6,
                                         m = NULL,
                                         alpha = 0.05, B = 5000, seed = 1,
                                         clip_pE = TRUE,
                                         ipcw_method = c("ht", "hajek"),
                                         n_min = 0, pi_min = 0.9) {

  stopifnot(length(S) == length(E))
  stopifnot(is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1)
  S <- as.integer(S)
  n <- length(S)
  ipcw_method <- match.arg(ipcw_method)

  clip01 <- function(x) pmin(1, pmax(0, x))

  # Two-sided percentile interval over finite values; [0,1] if none.
  percentile_interval <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(c(0, 1))
    clip01(as.numeric(stats::quantile(x, probs = c(alpha / 2, 1 - alpha / 2),
                                      names = FALSE, type = 7)))
  }

  # ---- IPCW contributions ------------------------------------------------
  use_ipcw <- !is.null(O) || !is.null(prob_O)
  if (use_ipcw) {
    if (is.null(O) || is.null(prob_O)) stop("Both O and prob_O are required for IPCW.")
    stopifnot(length(O) == n, length(prob_O) == n)
    O      <- as.integer(O)
    prob_O <- as.numeric(prob_O)
    if (anyNA(O) || any(!(O %in% 0:1)))      stop("O must be in {0,1} with no NA.")
    if (any(is.na(E) & O == 1L))             stop("E must be observed where O == 1.")
    if (any(!(is.na(E) | E %in% c(0, 1))))   stop("E must be 0/1/NA.")
    if (anyNA(prob_O[O == 1L]))              stop("prob_O has NA where O == 1.")

    denom  <- pmax(prob_O, prob_O_min)
    weight <- ifelse(O == 1L, 1 / denom, 0)                 # w_i = delta_i / max(G, gamma)
    ht_contrib <- ifelse(O == 1L, as.numeric(E) / denom, 0) # E_i * w_i
  } else {
    if (anyNA(E) || any(!(E %in% c(0, 1)))) stop("E must be in {0,1} without IPCW inputs.")
    weight     <- rep(1, n)
    ht_contrib <- as.numeric(E)
  }

  # Event-rate estimator among selected, for a given index set
  event_rate <- function(idx) {
    sel <- idx[S[idx] == 1L]
    if (length(sel) == 0L) return(NA_real_)
    if (ipcw_method == "ht") {
      r <- mean(ht_contrib[sel])
    } else {                      # Hajek
      w <- weight[sel]
      if (sum(w) == 0) return(NA_real_)
      r <- sum(w * ifelse(is.na(E[sel]), 0, E[sel])) / sum(w)
    }
    if (clip_pE) clip01(r) else r
  }

  # ---- Point estimates and vacuous rule (n_min) --------------------------
  n_selected <- sum(S == 1L)
  vacuous    <- n_selected < max(1L, n_min)

  pS_hat <- n_selected / n
  pE_hat <- if (vacuous) NA_real_ else event_rate(seq_len(n))

  # ---- Bootstrap of calibration data + Monte Carlo cohort -----------------
  set.seed(seed)
  idx_mat <- matrix(sample.int(n, n * B, replace = TRUE), nrow = B, ncol = n)
  if (!is.null(m)) {
    U_sel <- matrix(stats::runif(B * m), nrow = B, ncol = m)
    U_evt <- matrix(stats::runif(B * m), nrow = B, ncol = m)
  }

  pS_boot <- numeric(B)
  pE_boot <- numeric(B)
  yield_m <- if (!is.null(m)) numeric(B) else NULL
  rate_m  <- if (!is.null(m)) numeric(B) else NULL

  for (b in seq_len(B)) {
    idx <- idx_mat[b, ]
    pS_boot[b] <- mean(S[idx] == 1L)
    pE_boot[b] <- event_rate(idx)

    if (!is.null(m)) {
      selected   <- U_sel[b, ] <= pS_boot[b]
      A          <- sum(selected)
      yield_m[b] <- A / m
      rate_m[b]  <- if (A > 0L && is.finite(pE_boot[b])) {
        sum(selected & (U_evt[b, ] <= pE_boot[b])) / A
      } else NA_real_
    }
  }

  # ---- Pointwise confidence intervals ------------------------------------
  CI_pS <- percentile_interval(pS_boot)
  CI_pE <- if (vacuous) c(0, 1) else percentile_interval(pE_boot)

  out <- list(
    settings = list(alpha = alpha, B = B, seed = seed, m = m,
                    clip_pE = clip_pE, prob_O_min = prob_O_min,
                    ipcw_used = use_ipcw, ipcw_method = ipcw_method,
                    n_min = n_min, pi_min = pi_min,
                    n_selected = n_selected, vacuous = vacuous,
                    fs_smoothing = FALSE),
    estimate = list(pS = pS_hat, pE = pE_hat),
    confint  = list(pS = CI_pS, pE = CI_pE)
  )

  # ---- Pointwise prediction intervals ------------------------------------
  if (!is.null(m)) {
    PI_yield <- percentile_interval(yield_m)

    nonmissing_frac <- mean(is.finite(rate_m))
    PI_rate <- if (vacuous || nonmissing_frac < pi_min) c(0, 1) else percentile_interval(rate_m)

    # Containment: PI must contain CI
    PI_rate <- c(min(PI_rate[1], CI_pE[1]), max(PI_rate[2], CI_pE[2]))

    out$predint <- list(m = m,
                        prop_selected = PI_yield,
                        prop_event_among_selected = PI_rate)
  }

  class(out) <- "selection_event_bounds"
  out
}
