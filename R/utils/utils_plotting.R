plot_calibration_bands <- function(band_list,
                                   test_cohort_size = NULL,
                                   plot_title = NULL,
                                   plot_subtitle = NULL) {
  stopifnot(is.list(band_list), length(band_list) > 0)
  if (is.null(names(band_list)) || any(names(band_list) == "")) {
    stop("band_list must be a named list. Names will be used as Method labels.")
  }

  req_cols <- c(
    "lambda",
    "yield_hat", "ppv_hat",
    "yield_pop_lo", "yield_pop_hi",
    "ppv_pop_lo", "ppv_pop_hi"
  )

  # ---- Combine results ----
  band_df <- dplyr::bind_rows(
    lapply(names(band_list), function(nm) {
      df <- band_list[[nm]]
      missing_cols <- setdiff(req_cols, colnames(df))
      if (length(missing_cols) > 0) {
        stop(sprintf("Method '%s' is missing columns: %s",
                     nm, paste(missing_cols, collapse = ", ")))
      }
      dplyr::mutate(df, Method = nm)
    })
  )

  band_df$Method <- factor(band_df$Method, levels = names(band_list))

  # ---- Decide whether prediction bands are available ----
  has_pred_band <- all(c("selected_test_lo", "selected_test_hi", "ppv_test_lo", "ppv_test_hi") %in% colnames(band_df))

  # Determine label for prediction band
  if (is.null(test_cohort_size) && "num_test" %in% colnames(band_df)) {
    nt <- unique(band_df$num_test)
    if (length(nt) == 1 && !is.na(nt)) test_cohort_size <- nt
  }
  pred_lab <- if (has_pred_band) {
    if (!is.null(test_cohort_size)) paste0("Prediction (m=", test_cohort_size, ")") else "Prediction"
  } else {
    NULL
  }

  # ---- Convert to long format ----
  conf_long <- dplyr::bind_rows(
    band_df |>
      dplyr::transmute(
        lambda, Method,
        metric = "Yield",
        band   = "Confidence",
        lo     = yield_pop_lo,
        hi     = yield_pop_hi,
        hat    = yield_hat
      ),
    band_df |>
      dplyr::transmute(
        lambda, Method,
        metric = "PPV",
        band   = "Confidence",
        lo     = ppv_pop_lo,
        hi     = ppv_pop_hi,
        hat    = ppv_hat
      )
  )

  if (has_pred_band) {
    if (!("num_test" %in% colnames(band_df))) {
      stop("Prediction bands require 'num_test' to convert selected_test_lo/hi to yield.")
    }

    pred_long <- dplyr::bind_rows(
      band_df |>
        dplyr::transmute(
          lambda, Method,
          metric = "Yield",
          band   = pred_lab,
          lo     = selected_test_lo / num_test,
          hi     = selected_test_hi / num_test,
          hat    = yield_hat
        ),
      band_df |>
        dplyr::transmute(
          lambda, Method,
          metric = "PPV",
          band   = pred_lab,
          lo     = ppv_test_lo,
          hi     = ppv_test_hi,
          hat    = ppv_hat
        )
    )
    ribbon_long <- dplyr::bind_rows(conf_long, pred_long)
  } else {
    ribbon_long <- conf_long
  }

  line_long <- dplyr::bind_rows(
    band_df |> dplyr::transmute(lambda, Method, metric = "Yield", hat = yield_hat),
    band_df |> dplyr::transmute(lambda, Method, metric = "PPV",   hat = ppv_hat)
  )

  # ---- Plot ----
  band_levels <- c("Confidence", pred_lab)
  band_levels <- band_levels[!is.na(band_levels)]
  ribbon_long$band <- factor(ribbon_long$band, levels = band_levels)

  fill_vals  <- c("Confidence" = "blue")
  alpha_vals <- c("Confidence" = 0.30)
  if (!is.null(pred_lab)) {
    fill_vals[pred_lab]  <- "red"
    alpha_vals[pred_lab] <- 0.18
  }

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = ribbon_long,
      ggplot2::aes(x = lambda, ymin = lo, ymax = hi, fill = band, alpha = band)
    ) +
    ggplot2::geom_line(
      data = line_long,
      ggplot2::aes(x = lambda, y = hat, group = Method),
      linewidth = 1
    ) +
    ggplot2::facet_wrap(
      metric ~ Method,
      scales = "free_x",
      labeller = ggplot2::labeller(
        metric = ggplot2::as_labeller(c(
          Yield = "Yield: P(selected)",
          PPV   = "PPV: P(no event by t | selected)"
        )),
        Method = ggplot2::label_value
      )
    ) +
    ggplot2::scale_fill_manual(name = "Band type", values = fill_vals) +
    ggplot2::scale_alpha_manual(name = "Band type", values = alpha_vals) +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Selection threshold \u03bb",
      y = NULL
    ) +
    ggplot2::theme_minimal()
}
