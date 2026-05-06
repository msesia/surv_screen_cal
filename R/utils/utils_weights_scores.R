## =========================================================
## Fast survival prediction, IPCW weights, and score utilities
##   - Scores:  "survival" or "one_minus_survival"
##   - IPCW:    "et" (event time) or "ft" (fixed time)
##   - 'times' can be scalar or length-n vector
## =========================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

## ---------------------------------------------------------
## Fast survival prediction at arbitrary times (per-row times)
## Expects: model$predict(data, time_grid)$predictions -> n x |time_grid| matrix
## ---------------------------------------------------------
fast_predict_at_times <- function(model, data, target_times, time_grid = NULL) {
    target_times <- as.numeric(target_times)
    n <- nrow(data)
    if (length(target_times) == 1L) target_times <- rep(target_times, n)
    stopifnot(length(target_times) == n)

    time_grid <- if (is.null(time_grid)) sort(unique(target_times)) else sort(unique(as.numeric(time_grid)))
    pred <- model$predict(data, time_grid)
    surv_curves <- pred$predictions  ## n x |time_grid|
    if (all(target_times %in% time_grid)) {
        idx_time <- match(target_times, time_grid)
        surv_curves[cbind(seq_len(n), idx_time)]
    } else {
        vapply(seq_len(n), function(i) {
            stats::approx(time_grid, surv_curves[i, ], xout = target_times[i], rule = 2)$y
        }, numeric(1))
    }
}


## ---------------------------------------------------------
## IPCW weights (model-agnostic)
##   - et: w_i = 1 / Ĝ(Y_i^-|Z_i) for events; non-events -> NA
##   - ft: w_i = 1 / Ĝ(t0_i|Z_i)  for all i  (t0: scalar or length-n)
## ---------------------------------------------------------
build_censor_time_grid <- function(time, status, t0 = NULL) {
    out <- sort(unique(as.numeric(time[status == 1L])))
    if (!is.null(t0)) out <- sort(unique(c(out, as.numeric(t0))))
    out
}

## ---------------------------------------------------------
## Compute survival function of censoring distribution
## for IPCW using winsorization
##   method:
##     - "et" (event-time):  w_i = 1 / Ĝ(Y_i^-|Z_i) for events; non-events -> NA
##     - "ft" (fixed-time):  w_i = 1 / Ĝ(t0_i|Z_i)  for all i
##   trim: fraction of the *upper tail* to cap (default 0.01 → 99th pct cap)
## ---------------------------------------------------------
compute_censoring_survival <- function(
                                       data, cens_model, time, status,
                                       method = c("et","ft"),
                                       t0 = NULL,                      ## scalar or length-n if "ft"
                                       fast = TRUE, time_grid = NULL,
                                       trim = 0.01                     ## reasonable default; set 0 to disable
                                       ) {
    method <- match.arg(method)
    n <- nrow(data)
    time <- as.numeric(time); status <- as.integer(status)
    stopifnot(length(time) == n, length(status) == n)

    if (method == "ft") {
        if (is.null(t0)) stop("t0 must be provided for method='ft'.")
        t0 <- as.numeric(t0); if (length(t0) == 1L) t0 <- rep(t0, n)
        stopifnot(length(t0) == n)
    }

    ## Predict censoring survival Ĝ
    if (method == "et") {
        prob <- rep(NA_real_, n)
        ev_idx <- which(status == 1L)
        if (length(ev_idx)) {
            Ghat <- if (fast) {
                        fast_predict_at_times(
                            cens_model, data[ev_idx, , drop = FALSE], time[ev_idx],
                            time_grid = if (is.null(time_grid)) build_censor_time_grid(time, status) else time_grid
                        )
                    } else {
                        vapply(ev_idx, function(i) {
                            as.numeric(cens_model$predict(data[i, , drop = FALSE], time[i])$predictions)
                        }, numeric(1))
                    }
            prob[ev_idx] <- pmax(Ghat, 1e-12)  ## tiny floor for stability
        }
    } else { ## "ft"
        Ghat <- if (fast) {
                    fast_predict_at_times(
                        cens_model, data, t0,
                        time_grid = if (is.null(time_grid)) build_censor_time_grid(time, status, t0) else time_grid
                    )
                } else {
                    vapply(seq_len(n), function(i) {
                        as.numeric(cens_model$predict(data[i, , drop = FALSE], t0[i])$predictions)
                    }, numeric(1))
                }
        prob <- pmax(Ghat, 1e-12)
    }

    ## Winsorize extreme probabilities (top 'trim' fraction)
    if (is.finite(trim) && trim > 0) {
        cap <- stats::quantile(prob[is.finite(prob)], probs = trim, na.rm = TRUE)
        prob[prob < cap] <- cap
    }
    prob
}

## Convenience wrappers
compute_censoring_survival_et <- function(data, cens_model, time, status,
                                          fast = TRUE, time_grid = NULL, trim = NULL) {
    compute_censoring_survival(data, cens_model, time, status, method = "et",
                               t0 = NULL, fast = fast, time_grid = time_grid, trim = trim)
}
compute_censoring_survival_ft <- function(data, cens_model, time, status, t0,
                                          fast = TRUE, time_grid = NULL, trim = NULL) {
    compute_censoring_survival(data, cens_model, time, status, method = "ft",
                               t0 = t0, fast = fast, time_grid = time_grid, trim = trim)
}

## ---------------------------------------------------------
## Scores (two types only): "survival" or "one_minus_survival"
##   - For low-risk screening:  score_type = "survival"        (select high scores)
##   - For high-risk screening: score_type = "one_minus_survival" (select high scores)
## ---------------------------------------------------------
compute_scores_from_model <- function(
                                      data, surv_model, times,                    ## scalar or length-n
                                      score_type = c("survival","one_minus_survival"),
                                      fast = TRUE, time_grid = NULL
                                      ) {
    score_type <- match.arg(score_type)
    n <- nrow(data)
    times <- as.numeric(times); if (length(times) == 1L) times <- rep(times, n)
    stopifnot(length(times) == n)

    S <- if (fast) {
             fast_predict_at_times(surv_model, data, pmax(0, times), time_grid = time_grid)
         } else {
             vapply(seq_len(n), function(i) {
                 as.numeric(surv_model$predict(data[i, , drop = FALSE], times[i])$predictions)
             }, numeric(1))
         }
    if (score_type == "survival") S else 1 - S
}

compute_scores_from_curves <- function(
                                       surv_curves, time_grid, times,             ## times: scalar or length-n
                                       score_type = c("survival","one_minus_survival")
                                       ) {
    score_type <- match.arg(score_type)
    n <- nrow(surv_curves)
    times <- as.numeric(times); if (length(times) == 1L) times <- rep(times, n)
    stopifnot(length(times) == n)

    S <- vapply(seq_len(n), function(i) {
        stats::approx(time_grid, surv_curves[i, ], xout = pmax(0, times[i]), rule = 2)$y
    }, numeric(1))
    if (score_type == "survival") S else 1 - S
}
