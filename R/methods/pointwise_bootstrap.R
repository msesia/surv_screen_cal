 #' Bootstrap intervals for selection rate and event risk among selected
#'
#' Wrapper around `bootstrap_pointwise_internal()` that optionally
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
bootstrap_pointwise <- function(
  S,
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
  ipcw_method = c("ht", "hajek")
) {
  stopifnot(length(S) == length(E))
  if (!is.null(O))      stopifnot(length(O) == length(S))
  if (!is.null(prob_O)) stopifnot(length(prob_O) == length(S))

  ## Conservatively merge bands
  envelope <- function(x, y) {c(min(x[1], y[1]), max(x[2], y[2]))}

  ## One place to call the internal function with consistent args.
  run_internal <- function(S, E, O, prob_O, prob_O_min) {
    bootstrap_pointwise_internal(
      S = S, E = E, O = O, prob_O = prob_O, prob_O_min = prob_O_min,
      m = m, alpha = alpha, B = B, seed = seed, clip_pE = clip_pE, ipcw_method=ipcw_method
    )
  }

  ## No smoothing: just delegate.
  if (!isTRUE(fs_smoothing)) {
    return(run_internal(S, E, O, prob_O, prob_O_min))
  }

  ## Keep prob_O_min scalar and <= 1 (matches your current behavior)
  prob_O_min2 <- min(1, prob_O_min)

  ## Two runs differ only in pseudo E*: pessimistic 0 vs optimistic 1
  estim_lower <- run_internal(c(S,0,1), c(E,0,0), c(O,1,1), c(prob_O,1,1), prob_O_min2)
  estim_higher  <- run_internal(c(S,1), c(E, 1), c(O,1), c(prob_O,1), prob_O_min2)
  
  ## ---- Merge into an envelope CI ----
  out <- estim_lower
  out$confint$pS <- envelope(estim_lower$confint$pS, estim_higher$confint$pS)
  out$confint$pE <- envelope(estim_lower$confint$pE, estim_higher$confint$pE)
  
  ### Pointwise estimate without smoothing
  out$estimate <- point_estimate_pS_pE(S, E, O, prob_O, prob_O_min, clip_pE, ipcw_method = ipcw_method)
  
  if (!is.null(m) && !is.null(estim_lower$predint) && !is.null(estim_higher$predint)) {
    out$predint$m <- estim_lower$predint$m
    out$predint$prop_selected <- envelope(estim_lower$predint$prop_selected, estim_higher$predint$prop_selected)
    out$predint$prop_event_among_selected <- envelope(estim_lower$predint$prop_event_among_selected,
                                                      estim_higher$predint$prop_event_among_selected)
  }

  ## Record smoothing metadata in settings.
  out$settings$fs_smoothing <- TRUE
  out$settings$fs_smoothing_detail <- "Envelope CI: lower from pessimistic (E*=0), upper from optimistic (E*=1)"
  out$settings$seed <- seed
  out$settings$B <- B

  class(out) <- "selection_event_bounds"
  
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
bootstrap_pointwise_internal <- function(
                                         S,
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
                                         A_min = 10,
                                         min_nonNA_frac = 0.9
                                         ) {
    stopifnot(length(S) == length(E))
    S <- as.integer(S)
    n <- length(S)
    ipcw_method <- match.arg(ipcw_method)

    if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
        stop("alpha must be a single number in (0,1).")
    }

    ## Symmetric percentile bounds
    q_bounds <- function(x) {
        x <- as.numeric(x)
        x <- x[is.finite(x)]
        if (length(x) == 0L) return(c(NA_real_, NA_real_))
        c(
            as.numeric(stats::quantile(x, probs = alpha/2,     names = FALSE, type = 7)),
            as.numeric(stats::quantile(x, probs = 1 - alpha/2, names = FALSE, type = 7))
        )
    }

    ## ---- construct per-subject contribution for the event risk estimator ----
    use_ipcw <- !is.null(O) || !is.null(prob_O)
    if (use_ipcw) {
        if (is.null(O) || is.null(prob_O))
            stop("If using IPCW, both O and prob_O must be supplied.")
        stopifnot(length(O) == n, length(prob_O) == n)

        O <- as.integer(O)
        prob_O <- as.numeric(prob_O)

        if (any(is.na(O))) stop("O contains NA.")
        if (any(!(O %in% c(0L, 1L)))) stop("O must be in {0,1}.")
        if (any(is.na(E) & O == 1L))
            stop("E has NA where O==1; E must be observed (0/1) when O==1.")
        if (any(!(is.na(E) | E %in% c(0, 1))))
            stop("E must be 0/1/NA.")
        if (any(is.na(prob_O[O == 1L]))) stop("prob_O has NA where O==1.")
        if (any(prob_O[O == 1L] < 0)) stop("prob_O must be >= 0 where O==1.")

        denom <- pmax(prob_O, prob_O_min)

        ## HT contribution vector (only used when ipcw_method=="ht")
        contribution <- numeric(n)
        idx_obs <- which(O == 1L)
        contribution[idx_obs] <- as.numeric(E)[idx_obs] / denom[idx_obs]
    } else {
        if (any(is.na(E)))
            stop("E contains NA but IPCW inputs (O, prob_O) were not provided.")
        if (any(!(E %in% c(0, 1))))
            stop("E must be in {0,1}.")
        contribution <- as.numeric(E)
        denom <- NULL
    }

    ## ---- bootstrap resampling indices ----
    set.seed(seed)
    idx_mat <- matrix(sample.int(n, n * B, replace = TRUE), nrow = B, ncol = n)

    ## Store bootstrap estimator draws (for pointwise CI)
    pS_hat <- rep(NA_real_, B)
    pE_hat <- rep(NA_real_, B)

    ## Predictive draws (uniform-style)
    prop_sel <- if (!is.null(m)) rep(NA_real_, B) else NULL
    ratio    <- if (!is.null(m)) rep(NA_real_, B) else NULL
    A_draw   <- if (!is.null(m)) integer(B) else NULL

    ## Shared uniforms for prediction (uniform-style)
    U_sel_mat <- if (!is.null(m)) matrix(stats::runif(B * m), nrow = B, ncol = m) else NULL
    U_evt_mat <- if (!is.null(m)) matrix(stats::runif(B * m), nrow = B, ncol = m) else NULL

    for (b in 1:B) {
        idx <- idx_mat[b, ]
        Sb  <- S[idx]
        K   <- sum(Sb == 1L)
        pS_b <- K / n

        ## pE_b (bootstrap estimator)
        if (!use_ipcw) {
            Eb <- contribution[idx]
            pE_b <- if (K == 0L) NA_real_ else mean(Eb[Sb == 1L])
        } else {
            Ob <- O[idx]
            denom_b <- denom[idx]

            if (K == 0L) {
                pE_b <- NA_real_
            } else if (ipcw_method == "ht") {
                Eb <- contribution[idx]
                pE_b <- mean(Eb[Sb == 1L])
            } else { # Hajek
                sel_obs <- (Sb == 1L) & (Ob == 1L)
                if (!any(sel_obs)) {
                    pE_b <- NA_real_
                } else {
                    w <- 1 / denom_b[sel_obs]
                    y <- as.numeric(E)[idx][sel_obs]
                    pE_b <- sum(w * y) / sum(w)
                }
            }
        }

        if (clip_pE && !is.na(pE_b)) pE_b <- max(0, min(1, pE_b))

        pS_hat[b] <- pS_b
        pE_hat[b] <- pE_b

        if (!is.null(m)) {
            U_sel <- U_sel_mat[b, ]
            sel  <- (U_sel <= pS_b)
            A    <- sum(sel)

            A_draw[b] <- A
            prop_sel[b] <- A / m

            if (A > 0L && is.finite(pE_b)) {
                U_evt <- U_evt_mat[b, ]
                evt <- (U_evt <= pE_b)
                Bcnt <- sum(sel & evt)
                ratio[b] <- Bcnt / A
            } else {
                ratio[b] <- NA_real_
            }
        }
    }

    ## ---- point estimates ----
    point_est <- point_estimate_pS_pE(S, E, O, prob_O, prob_O_min, clip_pE, ipcw_method = ipcw_method)
    pS_hat0 <- point_est$pS
    pE_hat0 <- point_est$pE

    ## ---- pointwise percentile CIs ----
    CI_pS <- q_bounds(pS_hat)
    CI_pE <- q_bounds(pE_hat)
    CI_pS[1] <- max(0, CI_pS[1]); CI_pS[2] <- min(1, CI_pS[2])
    CI_pE[1] <- max(0, CI_pE[1]); CI_pE[2] <- min(1, CI_pE[2])

    out <- list(
        settings = list(
            alpha = alpha,
            B = B,
            seed = seed,
            m = m,
            clip_pE = clip_pE,
            prob_O_min = prob_O_min,
            ipcw_used = use_ipcw,
            ipcw_method = ipcw_method,
            A_min = A_min,
            min_nonNA_frac = min_nonNA_frac
        ),
        estimate = list(pS = pS_hat0, pE = pE_hat0),
        confint = list(pS = CI_pS, pE = CI_pE)
    )

    ## ---- pointwise predictive percentile intervals ----
    if (!is.null(m)) {
        PI_prop <- q_bounds(prop_sel)
        PI_prop[1] <- max(0, PI_prop[1]); PI_prop[2] <- min(1, PI_prop[2])

        A_expect <- m * pS_hat0
        ok_ratio <- is.finite(ratio)
        nonNA_frac <- mean(ok_ratio)

        if (!is.finite(A_expect) || A_expect < A_min || nonNA_frac < min_nonNA_frac) {
            PI_ratio <- c(0, 1)
        } else {
            PI_ratio <- q_bounds(ratio[ok_ratio])
            if (any(!is.finite(PI_ratio))) PI_ratio <- c(0, 1)
            PI_ratio[1] <- max(0, PI_ratio[1]); PI_ratio[2] <- min(1, PI_ratio[2])
        }

        ## Optional: enforce PI contains CI
        CI_ratio <- out$confint$pE
        if (all(is.finite(CI_ratio)) && all(is.finite(PI_ratio))) {
            PI_ratio[1] <- min(PI_ratio[1], CI_ratio[1])
            PI_ratio[2] <- max(PI_ratio[2], CI_ratio[2])
        }

        PI_ratio[1] <- max(0, PI_ratio[1]); PI_ratio[2] <- min(1, PI_ratio[2])

        out$predint <- list(
            m = m,
            prop_selected = PI_prop,
            prop_event_among_selected = PI_ratio
        )
    }

    out
}
