print.selection_event_bands <- function(x, digits = 3, n_show = 7, ...) {

  fmt_num <- function(v) {
    if (is.na(v)) return("NA")
    sprintf("%.*f", digits, v)
  }

  fmt_ci <- function(lo, hi) {
    if (is.na(lo) && is.na(hi)) return("[NA, NA]")
    sprintf("[%.*f, %.*f]", digits, lo, digits, hi)
  }

  # ---- nominal level text ----
  if (!is.null(x$settings$alpha_lower) && !is.null(x$settings$alpha_upper)) {
    alpha_lower <- x$settings$alpha_lower
    alpha_upper <- x$settings$alpha_upper
    if (is.na(alpha_lower) || is.na(alpha_upper)) {
      level_txt <- "One-sided"
    } else {
      level <- 1 - (alpha_lower + alpha_upper)
      level_txt <- sprintf("%.*f%%", 1, 100 * level)
    }
  } else if (!is.null(x$settings$alpha)) {
    level_txt <- sprintf("%.*f%%", 1, 100 * (1 - x$settings$alpha))
  } else {
    level_txt <- "Unknown level"
  }

  # ---- extract ----
  lam <- x$estimate$lambda
  L <- length(lam)

  pS <- x$estimate$pS
  pE <- x$estimate$pE

  pS_ci <- x$confint_simul$pS  # L x 2 matrix [lo, hi]
  pE_ci <- x$confint_simul$pE

  # ---- pick λ indices to display ----
  n_show <- max(3L, as.integer(n_show))
  if (L <= n_show) {
    idx_show <- seq_len(L)
  } else {
    probs <- seq(0, 1, length.out = n_show)
    idx_show <- unique(pmax(1L, pmin(L, round(probs * (L - 1L) + 1L))))
  }

  # ---- helpers for summaries ----
  band_width <- function(ci_mat) {
    w <- ci_mat[, 2] - ci_mat[, 1]
    w[!is.finite(w)] <- NA_real_
    w
  }

  pS_w <- band_width(pS_ci)
  pE_w <- band_width(pE_ci)

  cat("\nConfidence bands over lambda (nominal level:", level_txt, ")\n")
  cat("  Grid: L =", L,
      "; lambda in [", fmt_num(min(lam)), ", ", fmt_num(max(lam)), "]\n", sep = "")

  # ---- summary over λ ----
  cat("\nSummary over lambda\n")
  cat("  pS:  min/med/max =", fmt_num(min(pS, na.rm = TRUE)), "/",
      fmt_num(stats::median(pS, na.rm = TRUE)), "/",
      fmt_num(max(pS, na.rm = TRUE)), "\n", sep = "")
  cat("  pE:  min/med/max =", fmt_num(min(pE, na.rm = TRUE)), "/",
      fmt_num(stats::median(pE, na.rm = TRUE)), "/",
      fmt_num(max(pE, na.rm = TRUE)), "\n", sep = "")

  if (any(is.finite(pS_w))) {
    cat("  pS band width: min/med/max =",
        fmt_num(min(pS_w, na.rm = TRUE)), "/",
        fmt_num(stats::median(pS_w, na.rm = TRUE)), "/",
        fmt_num(max(pS_w, na.rm = TRUE)), "\n", sep = "")
  }
  if (any(is.finite(pE_w))) {
    cat("  pE band width: min/med/max =",
        fmt_num(min(pE_w, na.rm = TRUE)), "/",
        fmt_num(stats::median(pE_w, na.rm = TRUE)), "/",
        fmt_num(max(pE_w, na.rm = TRUE)), "\n", sep = "")
  }

  # ---- small table ----
  cat("\nSelected lambda points\n")
  tab <- data.frame(
    lambda = lam[idx_show],
    pS_hat = pS[idx_show],
    pS_lo  = pS_ci[idx_show, 1],
    pS_hi  = pS_ci[idx_show, 2],
    pE_hat = pE[idx_show],
    pE_lo  = pE_ci[idx_show, 1],
    pE_hi  = pE_ci[idx_show, 2]
  )

  # format as fixed-width text
  header <- sprintf(
    "%10s  %8s  %17s  %8s  %17s",
    "lambda", "pS", "pS band", "pE", "pE band"
  )
  cat(header, "\n", sep = "")
  for (j in seq_along(idx_show)) {
    i <- idx_show[j]
    cat(sprintf(
      "%10s  %8s  %17s  %8s  %17s\n",
      fmt_num(lam[i]),
      fmt_num(pS[i]),
      fmt_ci(pS_ci[i, 1], pS_ci[i, 2]),
      fmt_num(pE[i]),
      fmt_ci(pE_ci[i, 1], pE_ci[i, 2])
    ))
  }

  # ---- prediction band (if present) ----
  if (!is.null(x$predint_simul) && !is.null(x$predint_simul$prop_selected)) {
    cat("\nPrediction band for yield A/m")
    if (!is.null(x$settings$m)) cat(" (m =", x$settings$m, ")")
    cat(" , nominal level:", level_txt, "\n")

    pS_pred <- x$predint_simul$prop_selected
    pS_pred_w <- band_width(pS_pred)
    if (any(is.finite(pS_pred_w))) {
      cat("  width min/med/max =",
          fmt_num(min(pS_pred_w, na.rm = TRUE)), "/",
          fmt_num(stats::median(pS_pred_w, na.rm = TRUE)), "/",
          fmt_num(max(pS_pred_w, na.rm = TRUE)), "\n", sep = "")
    }

    cat("\n  Yield prediction band at selected lambda points\n")
    header2 <- sprintf("%10s  %17s", "lambda", "A/m band")
    cat(header2, "\n", sep = "")
    for (j in seq_along(idx_show)) {
      i <- idx_show[j]
      cat(sprintf(
        "%10s  %17s\n",
        fmt_num(lam[i]),
        fmt_ci(pS_pred[i, 1], pS_pred[i, 2])
      ))
    }
  }

  invisible(x)
}


as_tibble.selection_event_bands <- function(x, ...) {

  lam  <- x$estimate$lambda
  L    <- length(lam)

  pS_hat <- x$estimate$pS
  pE_hat <- x$estimate$pE

  pS_ci <- x$confint_simul$pS
  pE_ci <- x$confint_simul$pE

  # Optional prediction band
  has_pred <- !is.null(x$predint_simul) &&
              !is.null(x$predint_simul$prop_selected)

  if (has_pred) {
    pS_pred <- x$predint_simul$prop_selected
    pE_pred <- x$predint_simul$prop_event_among_selected
  }

    tibble::tibble(
                lambda = lam,
                num_test = x$settings$m %||% NA_integer_,   # optional if you stored m
                num_cal = x$settings$n %||% NA_integer_,   # optional if you stored n
                num_sel = round(pS_hat * (x$settings$n %||% NA_real_)),                
                pS_hat = pS_hat,
                pE_hat = pE_hat,
                pS_lo  = pS_ci[, 1],
                pS_hi  = pS_ci[, 2],
                pE_lo  = pE_ci[, 1],
                pE_hi  = pE_ci[, 2],
                pS_pred_lo = if (has_pred) pS_pred[, 1] else NA_real_,
                pS_pred_hi = if (has_pred) pS_pred[, 2] else NA_real_,
                pE_pred_lo = if (has_pred) pE_pred[, 1] else NA_real_,
                pE_pred_hi = if (has_pred) pE_pred[, 2] else NA_real_
            ) %>%
        mutate(yield_hat = pS_hat, ppv_hat=1-pE_hat, yield_pop_lo = pS_lo, yield_pop_hi = pS_hi,
               ppv_pop_lo = 1-pE_hi, ppv_pop_hi = 1-pE_lo,
               selected_test_lo = floor(num_test * pS_pred_lo), selected_test_hi = ceiling(num_test * pS_pred_hi),
               ppv_test_lo = 1-pE_pred_hi, ppv_test_hi = 1-pE_pred_lo) %>%
        select(-pS_hat, -pE_hat, -pS_lo, -pS_hi, -pE_lo, -pE_hi, -pS_pred_lo, -pS_pred_hi, -pE_pred_lo, -pE_pred_hi)
      
}
