split_data_n <- function(data, n_train, n_cal, n_test, shuffle = TRUE) {
  n_total <- nrow(data)
  stopifnot(n_train + n_cal + n_test <= n_total)

  if (shuffle) {
    idx <- sample(seq_len(n_total))
  } else {
    idx <- seq_len(n_total)
  }

  idx_train <- idx[1:n_train]
  idx_cal <- idx[(n_train + 1):(n_train + n_cal)]
  idx_test <- idx[(n_train + n_cal + 1):(n_train + n_cal + n_test)]

  return(list(
    data.train = data[idx_train, , drop = FALSE],
    data.cal   = data[idx_cal,   , drop = FALSE],
    data.test  = data[idx_test,  , drop = FALSE]
  ))
}
