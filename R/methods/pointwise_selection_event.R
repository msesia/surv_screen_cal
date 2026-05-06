print.selection_event_bounds <- function(x, digits = 3, ...) {

  fmt_num <- function(v) {
    if (is.na(v)) return("NA")
    sprintf("%.*f", digits, v)
  }

  fmt_ci <- function(ci) {
    if (all(is.na(ci))) return("[NA, NA]")
    sprintf("[%.*f, %.*f]", digits, ci[1], digits, ci[2])
  }

    ## ---- nominal level text ----
    if (!is.null(x$settings$alpha_lower) &&
        !is.null(x$settings$alpha_upper)) {

        alpha_lower <- x$settings$alpha_lower
        alpha_upper <- x$settings$alpha_upper

        if (is.na(alpha_lower) || is.na(alpha_upper)) {
            level_txt <- "One-sided"
        } else {
            level <- 1 - (alpha_lower + alpha_upper)
            level_txt <- sprintf("%.*f%%", 1, 100 * level)
        }

    } else if (!is.null(x$settings$alpha)) {

        ## finite-sample case
        level_txt <- sprintf("%.*f%%", 1, 100 * (1 - x$settings$alpha))

    } else {
        level_txt <- "Unknown level"
    }


    cat("\nConfidence intervals (nominal level:", level_txt, ")\n")

  if (!is.null(x$estimate)) {
    cat("  Selection proportion (pS):      ",
        fmt_num(x$estimate$pS), " ",
        fmt_ci(x$confint$pS), "\n", sep = "")
    cat("  Event risk among selected (pE): ",
        fmt_num(x$estimate$pE), " ",
        fmt_ci(x$confint$pE), "\n", sep = "")
  }

  if (!is.null(x$predint)) {
    m_txt <- x$settings$m
    cat("\nPredictive intervals (m =", m_txt,
        ", nominal level:", level_txt, ")\n")

    cat("  Proportion selected:          ",
        fmt_num(x$estimate$pS), " ",
        fmt_ci(x$predint$prop_selected), "\n", sep = "")
    cat("  Event risk among selected:    ",
        fmt_num(x$estimate$pE), " ",
        fmt_ci(x$predint$prop_event_among_selected), "\n", sep = "")
  }

  invisible(x)
}

# Horvitz–Thompson–style IPCW estimator
plugin_pE_ipcw_ht <- function(S, E, O, prob_O, prob_O_min = 1e-6) {
  S <- as.integer(S)
  if (sum(S == 1L) == 0L) return(0)
  contrib <- numeric(length(S))
  idx_obs <- which(O == 1L)
  denom <- pmax(prob_O, prob_O_min)
  # E must be observed where O==1 (caller should validate)
  contrib[idx_obs] <- as.numeric(E[idx_obs]) / denom[idx_obs]
  mean(contrib[S == 1L])
}

# Hájek (ratio / normalized-weight) IPCW estimator
plugin_pE_ipcw_hajek <- function(S, E, O, prob_O, prob_O_min = 1e-6) {
  S <- as.integer(S)
  if (sum(S == 1L) == 0L) return(0)
  denom <- pmax(prob_O, prob_O_min)
  # Only observed outcomes contribute
  idx <- which((S == 1L) & (O == 1L))
  if (length(idx) == 0L) return(0)  # or NA_real_ if you prefer "undefined"
  w <- 1 / denom[idx]
  sum(w * as.numeric(E[idx])) / sum(w)
}

# Unified IPCW estimator for pE among selected
# ipcw_method = "hajek" (ratio/normalized) or "ht" (Horvitz–Thompson / unnormalized mean)
plugin_pE_ipcw <- function(S, E, O, prob_O, prob_O_min = 1e-6, ipcw_method = c("hajek", "ht")) {
  method <- match.arg(ipcw_method)
  if (is.null(O) || is.null(prob_O))
    stop("Both O and prob_O must be supplied for IPCW.")
  switch(
    method,
    hajek = plugin_pE_ipcw_hajek(S, E, O, prob_O, prob_O_min),
    ht    = plugin_pE_ipcw_ht(S, E, O, prob_O, prob_O_min)
  )
}

point_estimate_pS_pE <- function(S, E, O = NULL, prob_O = NULL, prob_O_min = 1e-6, clip_pE = TRUE, ipcw_method = c("hajek", "ht")) {
    S <- as.integer(S)
    K <- sum(S == 1L)
    pS_hat <- K / length(S)
    use_ipcw <- !is.null(O) || !is.null(prob_O)
    ## If no selected, pE is undefined
    if (K == 0L) {
        return(list(pS = pS_hat, pE = NA_real_, ipcw_used = use_ipcw))
    }
    if (!use_ipcw) {
        pE_hat <- mean(as.numeric(E)[S == 1L])
    } else {
        if (is.null(O) || is.null(prob_O))
            stop("If using IPCW, both O and prob_O must be supplied.")
        pE_hat <- plugin_pE_ipcw(S, E, O, prob_O, prob_O_min = prob_O_min, ipcw_method = match.arg(ipcw_method)
        )
    }
    if (clip_pE && !is.na(pE_hat))
        pE_hat <- pmin(1, pmax(0, pE_hat))
    list(pS = pS_hat, pE = pE_hat, ipcw_used = use_ipcw)
}


hybrid_pointwise <- function(
  S,
  E,
  O = NULL,
  prob_O = NULL,
  prob_O_min = 1e-6,
  m = NULL,
  alpha = 0.05,
  alpha_lower = alpha / 2,
  alpha_upper = alpha / 2,
  B = 5000,
  seed = 1,
  ratio_if_A0 = NA_real_,
  clip_pE = TRUE,
  fs_smoothing = TRUE,
  ipcw_method = c("hajek", "ht"),
  method_pE_fs = "empirical_bernstein",
  fs_threshold = 10L,
  fs_threshold_n = fs_threshold,
  fs_threshold_sel = fs_threshold
) {
  stopifnot(length(S) == length(E))
  if (!is.null(O))      stopifnot(length(O) == length(S))
  if (!is.null(prob_O)) stopifnot(length(prob_O) == length(S))

  ipcw_method <- match.arg(ipcw_method)

  S_int <- as.integer(S)
  n   <- length(S_int)
  sel <- sum(S_int == 1L)

  use_fs_pS <- (n   < as.integer(fs_threshold_n))
  use_fs_pE <- (sel < as.integer(fs_threshold_sel))

  # ---- If both small → pure finite-sample ----
  if (use_fs_pS && use_fs_pE) {
    out <- finite_sample_pointwise(
      S = S, E = E, O = O, prob_O = prob_O,
      prob_O_min = prob_O_min,
      m = m,
      alpha = alpha,
      ratio_if_A0 = ratio_if_A0,
      clip_pE = clip_pE,
      method_pE = method_pE_fs
    )

    out$settings$hybrid <- TRUE
    out$settings$hybrid_rule <- sprintf(
      "FS used for both: n < %d AND selected < %d",
      as.integer(fs_threshold_n), as.integer(fs_threshold_sel)
    )

    class(out) <- "selection_event_bounds"
    return(out)
  }

  # ---- Otherwise run bootstrap once ----
  boot_out <- bootstrap_bounds_selection_and_event_rate(
    S = S, E = E, O = O, prob_O = prob_O,
    prob_O_min = prob_O_min,
    m = m,
    alpha = alpha,
    alpha_lower = alpha_lower,
    alpha_upper = alpha_upper,
    B = B,
    seed = seed,
    ratio_if_A0 = ratio_if_A0,
    clip_pE = clip_pE,
    fs_smoothing = fs_smoothing,
    ipcw_method = ipcw_method
  )

  out <- boot_out

  # ---- Replace components selectively with FS ----
  if (use_fs_pS || use_fs_pE) {
    fs_out <- finite_sample_pointwise(
      S = S, E = E, O = O, prob_O = prob_O,
      prob_O_min = prob_O_min,
      m = m,
      alpha = alpha,
      ratio_if_A0 = ratio_if_A0,
      clip_pE = clip_pE,
      method_pE = method_pE_fs
    )

    if (use_fs_pS) out$confint$pS <- fs_out$confint$pS
    if (use_fs_pE) out$confint$pE <- fs_out$confint$pE

    if (!is.null(m)) {
      tmp <- predint_wc_yield_and_event_rate(
        m,
        out$confint$pS[1], out$confint$pS[2],
        out$confint$pE[1], out$confint$pE[2],
        alpha = alpha,
        ratio_if_A0 = ratio_if_A0
      )
      out$predint <- tmp$predint
    }
  }

  out$settings$hybrid <- TRUE
  out$settings$fs_threshold_n <- as.integer(fs_threshold_n)
  out$settings$fs_threshold_sel <- as.integer(fs_threshold_sel)
  out$settings$used_fs_pS <- use_fs_pS
  out$settings$used_fs_pE <- use_fs_pE

  class(out) <- "selection_event_bounds"
  out
}
