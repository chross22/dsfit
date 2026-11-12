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
  # Adjustment models are fitted under the monotonicity constraint, which mrds
  # solves with a different optimiser (SLSQP). Its answer differs between
  # platforms in the fourth decimal: on CI, Linux and Windows gave p = 0.3337
  # where macOS gave 0.3338, and hn+cos(2) and hn+herm(4), 0.0006 AIC apart,
  # swapped places. A text snapshot cannot absorb that, so this one test
  # compares with a tolerance instead.
  #
  # Rows are ordered by name, so a near-tie cannot reorder them; delta_aic
  # still pins the ranking. The reference values are from the macOS run of
  # 2026-11-12. As with the snapshots, if this fails, find out what changed
  # in mrds before updating the values.
  sw <- suppressWarnings(sweep_models(
    exact_data(),
    model_set(key = c("hn", "unif"), adjustment = c("cos", "herm"),
              order = 2:4),
    truncation = 400,
    quiet = TRUE
  ))
  tab <- selection_table(sw, converged_only = FALSE)
  # Radix ordering is the C locale on every platform; the default follows the
  # session's locale, which can put "(2)" and "(2,3,4)" in either order.
  tab <- tab[order(tab$model_id, method = "radix"), , drop = FALSE]

  expect_equal(tab$model_id, c(
    "hn+cos(2)", "hn+cos(2,3)", "hn+cos(2,3,4)", "hn+herm(4)",
    "unif+cos(1,2)", "unif+cos(1,2,3)", "unif+cos(1,2,3,4)",
    "unif+herm(2)", "unif+herm(2,4)"
  ))
  expect_true(all(tab$converged))
  expect_equal(tab$n_par, c(2, 3, 4, 2, 2, 3, 4, 1, 2))
  expect_equal(min(tab$aic), 2178.678, tolerance = 1e-6)

  tol <- 1e-3
  expect_equal(tab$delta_aic, c(0.376908, 1.02739, 2.03485, 0.377467, 0,
                                1.81085, 1.67281, 93.6607, 20.7909),
               tolerance = tol)
  expect_equal(tab$p, c(0.348919, 0.324286, 0.340222, 0.349886, 0.370677,
                        0.361466, 0.333757, 0.667925, 0.486222),
               tolerance = tol)
  expect_equal(tab$p_cv, c(0.0838921, 0.108031, 0.127339, 0.0725726,
                           0.0395222, 0.1044, 0.0928966, 0.0817352,
                           0.065755), tolerance = tol)
  expect_equal(tab$esw, c(139.5675, 129.7143, 136.0886, 139.9544, 148.2707,
                          144.5863, 133.5026, 267.1698, 194.4889),
               tolerance = tol)
  expect_equal(tab$cvm_p, c(0.931393, 0.908625, 0.97522, 0.925273, 0.52662,
                            0.743714, 0.95862, 0.00146519, 4.68338e-06),
               tolerance = tol)
  expect_equal(tab$chisq_p, c(0.529661, 0.597656, 0.590906, 0.530986,
                              0.461795, 0.434848, 0.601305, 7.32747e-14,
                              0.00364331), tolerance = tol)
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
  expect_match(paste(deparse(sw$fits[["hn bf"]]$dsmodel), collapse = ""),
               "mcds")
  expect_match(paste(deparse(sw$fits[["hn"]]$dsmodel), collapse = ""), "cds\\(")
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
