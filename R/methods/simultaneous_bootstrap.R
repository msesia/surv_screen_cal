# =============================================================================
# Simultaneous confidence and prediction bands for yield and event rate
# among selected patients, over a threshold grid and (optionally) K models.
#
# Implements Algorithms 1-3 of the paper. Notation (paper -> code):
#   n_min  : minimum number of selected calibration patients for the event
#            rate to be estimable at a cell (k, lambda); otherwise the cell is
#            "vacuous": point estimate NA, bands [0,1], excluded from the max.
#   pi_min : minimum fraction of simulated cohorts with >= 1 selected patient
#            for the predictive band of the event rate to be reported.
#   gamma  : prob_O_min (IPCW truncation).
#
# Input `scores` is a numeric vector (K = 1) or a list of K numeric vectors,
# one per model, each of length n. All randomness (multipliers, bootstrap
# indices, cohort uniforms) is shared across models and thresholds, and the
# critical values are quantiles of the maximum over all non-vacuous cells.
#
# Output layout is the same as before for K = 1. For K > 1, the same fields
# are returned as lists indexed by model (see `unpack_model()` at the end).
# =============================================================================


# ---- Wrapper: finite-sample smoothing (Algorithm 3) --------------------------

bootstrap_simultaneous <- function(scores, E,
                                   O = NULL, prob_O = NULL, prob_O_min = 1e-6,
                                   m = NULL,
                                   alpha = 0.05, B = 5000, seed = 1,
                                   clip_pE = TRUE,
                                   fs_smoothing = TRUE,
                                   ipcw_method = c("ht", "hajek"),
                                   lambda_seq = NULL,
                                   n_min = 0, pi_min = 0.9) {

  scores <- as_score_list(scores)
  K <- length(scores)
  n <- length(scores[[1]])
  stopifnot(all(lengths(scores) == n), length(E) == n)
  if (!is.null(O))      stopifnot(length(O) == n)
  if (!is.null(prob_O)) stopifnot(length(prob_O) == n)
  ipcw_method <- match.arg(ipcw_method)

  run <- function(scores, E, O, prob_O, prob_O_min) {
    bootstrap_simultaneous_internal(
      scores = scores, E = E, O = O, prob_O = prob_O, prob_O_min = prob_O_min,
      m = m, alpha = alpha, B = B, seed = seed, clip_pE = clip_pE,
      ipcw_method = ipcw_method, lambda_seq = lambda_seq,
      n_min = n_min, pi_min = pi_min
    )
  }

  base <- run(scores, E, O, prob_O, prob_O_min)
  if (!isTRUE(fs_smoothing)) {
    if (K == 1) base <- flatten_single_model(base)
    return(base)
  }

  # Pseudo-patients have score -Inf / +Inf under every model, so they are
  # never / always selected at every threshold.
  #   run A: one never-selected and one always-selected patient, no event
  #   run B: one always-selected patient with event
  prob_O_min2 <- min(1, prob_O_min)
  augment <- function(s, extra) c(s, extra)
  scores_A <- lapply(scores, augment, extra = c(-Inf, Inf))
  scores_B <- lapply(scores, augment, extra = Inf)

  run_A <- run(scores_A, c(E, 0, 0), c(O, 1, 1), c(prob_O, 1, 1), prob_O_min2)
  run_B <- run(scores_B, c(E, 1),    c(O, 1),    c(prob_O, 1),    prob_O_min2)

  envelope <- function(x, y) cbind(lo = pmin(x[, 1], y[, 1]), hi = pmax(x[, 2], y[, 2]))

  out <- base   # point estimates and settings from the unsmoothed run
  for (k in seq_len(K)) {
    out$confint_simul$pS[[k]] <- envelope(run_A$confint_simul$pS[[k]], run_B$confint_simul$pS[[k]])
    out$confint_simul$pE[[k]] <- envelope(run_A$confint_simul$pE[[k]], run_B$confint_simul$pE[[k]])
    if (!is.null(m)) {
      out$predint_simul$prop_selected[[k]] <-
        envelope(run_A$predint_simul$prop_selected[[k]], run_B$predint_simul$prop_selected[[k]])
      out$predint_simul$prop_event_among_selected[[k]] <-
        envelope(run_A$predint_simul$prop_event_among_selected[[k]],
                 run_B$predint_simul$prop_event_among_selected[[k]])
    }
  }
  out$settings$fs_smoothing <- TRUE

  if (K == 1) out <- flatten_single_model(out)
  class(out) <- "selection_event_bands"
  out
}


# ---- Internal worker: Algorithms 1 and 2, no smoothing -----------------------

bootstrap_simultaneous_internal <- function(scores, E,
                                            O = NULL, prob_O = NULL, prob_O_min = 1e-6,
                                            m = NULL,
                                            alpha = 0.05, B = 5000, seed = 1,
                                            clip_pE = TRUE,
                                            ipcw_method = c("ht", "hajek"),
                                            lambda_seq = NULL,
                                            n_min = 0, pi_min = 0.9) {

  scores <- as_score_list(scores)
  K <- length(scores)
  n <- length(scores[[1]])
  ipcw_method <- match.arg(ipcw_method)
  stopifnot(is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1)
  if (is.null(lambda_seq)) lambda_seq <- seq(0, 1, length.out = 10)
  L <- length(lambda_seq)

  clip01 <- function(x) pmin(1, pmax(0, x))
  q_upper <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    as.numeric(stats::quantile(x, probs = 1 - alpha, names = FALSE, type = 7))
  }

  # ---- IPCW contributions (per patient) -----------------------------------
  use_ipcw <- !is.null(O) || !is.null(prob_O)
  if (use_ipcw) {
    if (is.null(O) || is.null(prob_O)) stop("Both O and prob_O are required for IPCW.")
    stopifnot(length(O) == n, length(prob_O) == n)
    O <- as.integer(O); prob_O <- as.numeric(prob_O)
    if (anyNA(O) || any(!(O %in% 0:1)))    stop("O must be in {0,1} with no NA.")
    if (any(is.na(E) & O == 1L))           stop("E must be observed where O == 1.")
    if (any(!(is.na(E) | E %in% c(0, 1)))) stop("E must be 0/1/NA.")
    if (anyNA(prob_O[O == 1L]))            stop("prob_O has NA where O == 1.")
    w <- ifelse(O == 1L, 1 / pmax(prob_O, prob_O_min), 0)
  } else {
    if (anyNA(E) || any(!(E %in% c(0, 1)))) stop("E must be in {0,1} without IPCW inputs.")
    w <- rep(1, n)
  }
  E0 <- ifelse(is.na(E), 0, as.numeric(E))
  # Numerator / denominator contributions of the event-rate ratio estimator
  z1 <- E0 * w
  z0 <- if (ipcw_method == "ht") rep(1, n) else w

  # ---- Cells: columns of an n x (L*K) selection matrix --------------------
  # cell index c = (k - 1) * L + l  <->  (model k, threshold lambda_l)
  S <- do.call(cbind, lapply(scores, function(s) 1.0 * outer(s, lambda_seq, FUN = ">=")))
  C <- L * K
  cell_k <- rep(seq_len(K), each = L)

  # ---- Point estimates (Algorithm 1, step 1) ------------------------------
  pS_hat <- colMeans(S)
  mu1    <- colMeans(S * z1)
  mu0    <- colMeans(S * z0)

  # ---- Vacuous cells (Algorithm 1, step 2) --------------------------------
  n_selected <- colSums(S)
  vac <- n_selected < max(1L, n_min)
  pE_hat <- rep(NA_real_, C)
  pE_hat[!vac] <- mu1[!vac] / mu0[!vac]
  if (clip_pE) pE_hat <- clip01(pE_hat)

  # ---- Influence contributions and standard errors (step 3) ---------------
  psiS <- sweep(S, 2, pS_hat, "-")                                    # n x C
  psiE <- matrix(0, n, C)
  if (any(!vac)) {
    Z1c <- sweep((S * z1)[, !vac, drop = FALSE], 2, mu1[!vac], "-")
    Z0c <- sweep((S * z0)[, !vac, drop = FALSE], 2, mu0[!vac], "-")
    psiE[, !vac] <- sweep(Z1c, 2, 1 / mu0[!vac], "*") -
                    sweep(Z0c, 2, mu1[!vac] / mu0[!vac]^2, "*")
  }
  sdS <- sqrt(colMeans(psiS^2))
  sdE <- sqrt(colMeans(psiE^2)); sdE[vac] <- NA_real_

  okS <- is.finite(sdS) & sdS > 0
  okE <- !vac & is.finite(sdE) & sdE > 0

  # ---- Multiplier bootstrap, shared across all cells (step 4) -------------
  set.seed(seed)
  G <- matrix(stats::rnorm(B * n), B, n)

  sup_abs <- function(psi, sd, ok) {
    if (!any(ok)) return(rep(NA_real_, B))
    Z <- sweep((G %*% psi[, ok, drop = FALSE]) / sqrt(n), 2, sd[ok], "/")
    apply(abs(Z), 1, max)
  }
  cS <- q_upper(sup_abs(psiS, sdS, okS))
  cE <- q_upper(sup_abs(psiE, sdE, okE))

  # ---- Confidence bands (step 5) ------------------------------------------
  CB_S <- cbind(lo = pS_hat, hi = pS_hat)      # degenerate where sd = 0
  if (is.finite(cS)) {
    CB_S[okS, "lo"] <- pS_hat[okS] - cS * sdS[okS] / sqrt(n)
    CB_S[okS, "hi"] <- pS_hat[okS] + cS * sdS[okS] / sqrt(n)
  }
  CB_E <- cbind(lo = rep(0, C), hi = rep(1, C))  # vacuous by default
  if (is.finite(cE)) {
    CB_E[okE, "lo"] <- pE_hat[okE] - cE * sdE[okE] / sqrt(n)
    CB_E[okE, "hi"] <- pE_hat[okE] + cE * sdE[okE] / sqrt(n)
  }
  CB_S[] <- clip01(CB_S); CB_E[] <- clip01(CB_E)

  out <- list(
    settings = list(alpha = alpha, B = B, seed = seed, n = n, m = m, K = K,
                    prob_O_min = prob_O_min, ipcw_used = use_ipcw,
                    ipcw_method = ipcw_method, lambda_seq = lambda_seq,
                    n_min = n_min, pi_min = pi_min, fs_smoothing = FALSE,
                    critical_value_yield = cS, critical_value_er = cE),
    estimate = list(lambda = lambda_seq,
                    pS = split_by_model(pS_hat, cell_k),
                    pE = split_by_model(pE_hat, cell_k)),
    confint_simul = list(lambda = lambda_seq,
                         pS = split_by_model(CB_S, cell_k),
                         pE = split_by_model(CB_E, cell_k))
  )

  # ---- Prediction bands (Algorithm 2) -------------------------------------
  if (!is.null(m)) {

    set.seed(seed)
    idx_mat <- matrix(sample.int(n, n * B, replace = TRUE), B, n)
    U_sel   <- matrix(stats::runif(B * m), B, m)
    U_evt   <- matrix(stats::runif(B * m), B, m)

    yield_m <- matrix(NA_real_, B, C)   # realized yield in simulated cohort
    rate_m  <- matrix(NA_real_, B, C)   # realized event rate (NA if none selected)

    for (b in seq_len(B)) {
      idx  <- idx_mat[b, ]
      S_b  <- S[idx, , drop = FALSE]
      pS_b <- colMeans(S_b)
      mu0b <- colMeans(S_b * z0[idx])
      mu1b <- colMeans(S_b * z1[idx])
      pE_b <- ifelse(mu0b > 0, mu1b / mu0b, NA_real_)
      if (clip_pE) pE_b <- clip01(pE_b)

      sel <- outer(U_sel[b, ], pS_b, "<=")            # m x C
      A   <- colSums(sel)
      yield_m[b, ] <- A / m
      evt <- outer(U_evt[b, ], pE_b, "<=")
      ok  <- A > 0 & is.finite(pE_b)
      rate_m[b, ok] <- colSums(sel & evt)[ok] / A[ok]
    }

    # Studentisation scales (Algorithm 2, step 3)
    cp_bounds <- t(sapply(n_selected, function(k)
        bound_binomial_exact(k, n, alpha = 0.05, alpha_lower = 0.05, alpha_upper = 0.05)))
    ##sdA <- apply(sqrt(cp_bounds * (1 - cp_bounds) / m), 1, max)
    sdA <- pmax(apply(sqrt(cp_bounds * (1 - cp_bounds) / m), 1, max), 2 / m)
    sdR <- apply(rate_m, 2, stats::sd, na.rm = TRUE)

    # Vacuous cells for the event-rate prediction band (step 4)
    nonmissing_frac <- colMeans(is.finite(rate_m))
    vacR <- vac | !(nonmissing_frac >= pi_min) | !is.finite(sdR) | sdR <= 0

    # Studentised processes, critical values, bands (steps 5-7)
    sup_pred <- function(draws, center, sd, ok) {
      if (!any(ok)) return(rep(NA_real_, B))
      Z <- sweep(sweep(draws[, ok, drop = FALSE], 2, center[ok], "-"), 2, sd[ok], "/")
      apply(abs(Z), 1, max, na.rm = TRUE)
    }
    okA <- is.finite(sdA) & sdA > 0
    qA <- q_upper(sup_pred(yield_m, pS_hat, sdA, okA))
    qR <- q_upper(sup_pred(rate_m,  pE_hat, sdR, !vacR))

    PB_S <- cbind(lo = pS_hat, hi = pS_hat)
    if (is.finite(qA)) {
      PB_S[okA, "lo"] <- pS_hat[okA] - qA * sdA[okA]
      PB_S[okA, "hi"] <- pS_hat[okA] + qA * sdA[okA]
    }
    PB_E <- cbind(lo = rep(0, C), hi = rep(1, C))
    if (is.finite(qR)) {
      PB_E[!vacR, "lo"] <- pE_hat[!vacR] - qR * sdR[!vacR]
      PB_E[!vacR, "hi"] <- pE_hat[!vacR] + qR * sdR[!vacR]
    }

    # Containment (step 8), then clipping
    PB_S[, "lo"] <- pmin(PB_S[, "lo"], CB_S[, "lo"]); PB_S[, "hi"] <- pmax(PB_S[, "hi"], CB_S[, "hi"])
    PB_E[, "lo"] <- pmin(PB_E[, "lo"], CB_E[, "lo"]); PB_E[, "hi"] <- pmax(PB_E[, "hi"], CB_E[, "hi"])
    PB_S[] <- clip01(PB_S); PB_E[] <- clip01(PB_E)

    out$predint_simul <- list(
      lambda = lambda_seq,
      prop_selected             = split_by_model(PB_S, cell_k),
      prop_event_among_selected = split_by_model(PB_E, cell_k),
      critical_value_pred_yield = qA,
      critical_value_pred_er    = qR
    )
  }

  class(out) <- "selection_event_bands"
  out
}


# ---- Helpers ------------------------------------------------------------------

# Accept a numeric vector (K = 1), a list of vectors, or an n x K matrix.
as_score_list <- function(scores) {
  if (is.list(scores)) return(lapply(scores, as.numeric))
  if (is.matrix(scores)) return(lapply(seq_len(ncol(scores)), function(k) scores[, k]))
  list(as.numeric(scores))
}

# Split a length-C vector or C x 2 matrix into a list of K per-model pieces.
split_by_model <- function(x, cell_k) {
  if (is.matrix(x)) lapply(split(seq_len(nrow(x)), cell_k), function(i) x[i, , drop = FALSE])
  else              split(x, cell_k)
}

# For K = 1, unwrap the per-model lists so the output matches the old layout.
flatten_single_model <- function(out) {
  unwrap <- function(x) if (is.list(x) && length(x) == 1 && !is.data.frame(x)) x[[1]] else x
  for (f in c("estimate", "confint_simul", "predint_simul")) {
    if (!is.null(out[[f]])) out[[f]] <- lapply(out[[f]], unwrap)
  }
  out
}

# For K > 1: extract the old single-model layout for model k
# (useful with convert_band_to_long_targets()).
unpack_model <- function(out, k) {
  pick <- function(x) lapply(x, function(v) if (is.list(v) && !is.data.frame(v)) v[[k]] else v)
  res <- list(settings = out$settings,
              estimate = pick(out$estimate),
              confint_simul = pick(out$confint_simul))
  if (!is.null(out$predint_simul)) res$predint_simul <- pick(out$predint_simul)
  class(res) <- "selection_event_bands"
  res
}
