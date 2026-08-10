# Regression snapshots over the selection table.
#
# The sweep *selects* a model; these pin the selection, so that an `mrds`
# upgrade cannot quietly change which model wins or by how much. That
# requirement is why this layer is a package rather than a script.
#
# The whole table is pinned rather than the winner alone. A swap between models
# three and four is the same underlying change arriving early, while it is still
# cheap to look at.
#
# When one of these fails the question is not "is the new number reasonable?"
# but "what changed underneath?". Check `news(package = "mrds")` across the
# versions either side before running `testthat::snapshot_accept()`.

skip_if_not_installed("mrds")

# Rounded to a precision an optimiser can be expected to reproduce across
# platforms and BLAS implementations. Absolute `aic` is kept alongside
# `delta_aic`: a likelihood change that moves every model by the same amount
# leaves the ranking untouched and shows up here alone.
snap_table <- function(sw) {
  tab <- selection_table(sw, converged_only = FALSE)
  data.frame(
    model     = tab$model_id,
    converged = tab$converged,
    n_par     = tab$n_par,
    aic       = round(tab$aic, 2),
    delta_aic = round(tab$delta_aic, 2),
    p         = round(tab$p, 4),
    p_cv      = round(tab$p_cv, 3),
    esw       = round(tab$esw, 1),
    cvm_p     = round(tab$cvm_p, 3),
    chisq_p   = round(tab$chisq_p, 3),
    row.names = NULL
  )
}

test_that("the ranking over key functions is stable", {
  sw <- sweep_models(
    exact_data(),
    model_set(key = c("hn", "hr", "gamma")),
    truncation = 400,
    quiet = TRUE
  )
  expect_snapshot(snap_table(sw))
})

test_that("the ranking over adjustment series and orders is stable", {
  # Some of these do not converge, and `mrds` warns that others are not
  # monotonic. Both are pinned here: which candidates in a set fail is part of
  # what an upgrade can change, and the failed rows have to stay on their own
  # rows rather than shifting the models below them up.
  sw <- suppressWarnings(sweep_models(
    exact_data(),
    model_set(key = c("hn", "unif"), adjustment = c("cos", "herm"),
              order = 2:3),
    truncation = 400,
    quiet = TRUE
  ))
  expect_snapshot(snap_table(sw))
})

test_that("the binned likelihood is stable", {
  # Binned and exact fits of the same distances are different likelihoods, so
  # this pins a quantity the exact snapshots cannot reach.
  d <- binned_data()
  sw <- sweep_models(
    d,
    model_set(key = c("hn", "hr")),
    truncation = max(d$distend),
    quiet = TRUE
  )
  expect_snapshot(snap_table(sw))
})

test_that("covariate models are stable, and mcds is the path taken", {
  sw <- sweep_models(
    covariate_data(),
    model_set(key = "hn", formula = list(~1, ~bf)),
    truncation = 400,
    quiet = TRUE
  )
  expect_snapshot(snap_table(sw))
})

test_that("left truncation is stable", {
  # The other route to a blind spot beneath the platform. Its own snapshot
  # because `left` changes both the data and the scaling of g(x).
  sw <- sweep_models(
    exact_data(),
    model_set(key = c("hn", "hr")),
    truncation = 400,
    left = 30,
    quiet = TRUE
  )
  expect_snapshot(snap_table(sw))
})

test_that("what the sweep reports about itself is stable", {
  # The printed object and the progress report are the parts a user reads. They
  # are pinned here so a change to either is a visible diff rather than a
  # surprise at the console.
  d <- rbind(exact_data(), data.frame(object = 9001:9003,
                                      distance = c(450, 800, NA)))
  expect_snapshot(
    sw <- sweep_models(d, model_set(key = c("hn", "hr", "gamma")),
                       truncation = 400)
  )
  expect_snapshot(print(sw))
})
