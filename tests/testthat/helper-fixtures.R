# Fixtures shared across test files.
#
# They are deterministic on purpose: the snapshot tests in
# `test-selection-snapshot.R` pin fitted quantities, which only means anything
# if the data behind them is fixed. Changing a seed, an `n`, or an `sd` here
# will invalidate every snapshot, and that is a data change rather than a
# regression — accept the new snapshots deliberately.

# Exact (point) distances, half-normal in shape.
exact_data <- function(n = 200, sd = 120, seed = 1) {
  set.seed(seed)
  d <- data.frame(object = seq_len(n), distance = abs(stats::rnorm(n, 0, sd)))
  d[d$distance < 500, , drop = FALSE]
}

# The same distances as intervals, which is the shape a `STRIP` scheme gives.
binned_data <- function(width = 100, ...) {
  d <- exact_data(...)
  d$distbegin <- floor(d$distance / width) * width
  d$distend <- d$distbegin + width
  d$distance <- NULL
  d
}

# One detection covariate, with a real effect on the scale so that the
# covariate model has something to find.
covariate_data <- function(n = 200, seed = 1) {
  set.seed(seed)
  bf <- stats::runif(n, 0, 5)
  d <- data.frame(
    object = seq_len(n),
    bf = bf,
    distance = abs(stats::rnorm(n, 0, 160 - 12 * bf))
  )
  d[d$distance < 500, , drop = FALSE]
}
