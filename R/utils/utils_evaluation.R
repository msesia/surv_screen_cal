evaluate_yield_and_event_rate <- function(
  data.test,
  selected_idx,
  t0,
  time_col = "time",
  status_col = "status",
  event_time_col = "event_time",
  status_event_value = 1L
) {
  stopifnot(is.data.frame(data.test))
  n <- nrow(data.test)

  if (is.null(selected_idx)) selected_idx <- integer(0)
  selected_idx <- unique(as.integer(selected_idx))
  selected_idx <- selected_idx[selected_idx >= 1L & selected_idx <= n]

  if (!(time_col %in% names(data.test))) stop("Missing column: ", time_col)
  if (!(status_col %in% names(data.test))) stop("Missing column: ", status_col)

  time <- data.test[[time_col]]
  status <- data.test[[status_col]]

  if (!is.numeric(time)) stop(time_col, " must be numeric.")
  if (any(is.na(time))) stop(time_col, " contains NA.")
  if (any(time < 0)) stop(time_col, " must be >= 0.")

  status <- as.integer(status)
  if (any(is.na(status))) stop(status_col, " contains NA.")
  if (any(!(status %in% c(0L, 1L)))) stop(status_col, " must be in {0,1}.")

  if (!is.numeric(t0) || length(t0) != 1L || is.na(t0) || t0 < 0) {
    stop("t0 must be a single nonnegative number.")
  }

  A <- length(selected_idx)
  yield_prop <- if (n == 0L) NA_real_ else A / n

  ## If nobody selected, define event rate as NA (undefined)
  if (A == 0L) {
    return(list(
      n = n,
      num_selected = 0L,
      yield = yield_prop,
      event_rate = NA_real_,
      event_rate_lo = NA_real_,
      event_rate_hi = NA_real_,
      event_rate_type = "undefined (no selected)"
    ))
  }

  ## Helper: among selected, compute optimistic/pessimistic bounds for P(event by t0)
  ## Using (time, status) with right-censoring:
  ##  - Known event by t0: status==1 and time<=t0
  ##  - Known event-free at t0: time>t0 (regardless of status; if time>t0 then not failed before t0)
  ##  - Ambiguous: status==0 and time<=t0 (censored before t0)
  sel_time <- time[selected_idx]
  sel_status <- status[selected_idx]

  n_event_known <- sum(sel_status == status_event_value & sel_time <= t0)
  n_eventfree_known <- sum(sel_time > t0)
  n_ambig <- A - n_event_known - n_eventfree_known

  ## Lower bound: treat ambiguous as no-event by t0
  lo <- n_event_known / A

  ## Upper bound: treat ambiguous as event by t0
  hi <- (n_event_known + n_ambig) / A

  ## If event_time column exists and is usable, compute exact event rate by t0
  exact_available <- (event_time_col %in% names(data.test)) &&
                     all(!is.na(data.test[[event_time_col]][selected_idx]))

  if (exact_available) {
    event_time <- data.test[[event_time_col]]
    if (!is.numeric(event_time)) stop(event_time_col, " must be numeric.")
    if (any(event_time[selected_idx] < 0)) stop(event_time_col, " must be >= 0.")

    ## Exact indicator: event by t0
    ## (Assumes event_time is the actual failure time for everyone, including censored rows.)
    exact <- mean(event_time[selected_idx] <= t0)

    return(list(
      n = n,
      num_selected = A,
      yield = yield_prop,
      event_rate = exact,
      event_rate_lo = exact,
      event_rate_hi = exact,
      event_rate_type = "exact (event_time available)"
    ))
  }

  list(
    n = n,
    num_selected = A,
    yield = yield_prop,
    event_rate = NA_real_,
    event_rate_lo = lo,
    event_rate_hi = hi,
    event_rate_type = "bounds (optimistic/pessimistic censoring)"
  )
}
