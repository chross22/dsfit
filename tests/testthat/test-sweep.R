skip_if_not_installed("mrds")

exact_data <- function(n = 200, sd = 120, seed = 1) {
  set.seed(seed)
  d <- data.frame(object = seq_len(n), distance = abs(stats::rnorm(n, 0, sd)))
  d[d$distance < 500, , drop = FALSE]
}

sweep_fx <- function(...) suppressMessages(sweep_models(...))

test_that("a sweep fits its model set and ranks the result", {
  sw <- sweep_fx(exact_data(), model_set(c("hn", "hr")), truncation = 400)

  expect_s3_class(sw, "dsfit_sweep")
  tab <- selection_table(sw)
  expect_equal(nrow(tab), 2)
  expect_true(all(tab$converged))

  # Sorted by AIC, with the best at zero.
  expect_false(is.unsorted(tab$aic))
  expect_equal(min(tab$delta_aic), 0)
})

test_that("the table carries what abundance actually depends on", {
  tab <- selection_table(sweep_fx(exact_data(), truncation = 400))
  for (nm in c("aic", "delta_aic", "p", "p_se", "p_cv", "esw", "cvm_p")) {
    expect_true(nm %in% names(tab), info = nm)
  }
  # Effective strip half-width is p times the truncation width.
  expect_equal(tab$esw, tab$p * 400)
  expect_equal(tab$p_cv, tab$p_se / tab$p)
})

test_that("truncation is required, because AIC depends on it", {
  expect_error(sweep_models(exact_data()), "truncation` is required")
  expect_error(sweep_models(exact_data(), truncation = -1))
})

test_that("the gamma key fits, and is not available from Distance::ds()", {
  sw <- sweep_fx(exact_data(), model_set(c("hn", "gamma")), truncation = 400)
  tab <- selection_table(sw)
  expect_true("gamma" %in% tab$key)
  expect_true(all(tab$converged))

  # The point of reporting esw: the two keys imply different strip widths even
  # when their AICs are close.
  expect_false(isTRUE(all.equal(tab$esw[1], tab$esw[2])))
})

test_that("gamma together with left truncation warns about double counting", {
  expect_warning(
    suppressMessages(sweep_models(exact_data(), model_set(c("hn", "gamma")),
                                  truncation = 400, left = 20)),
    "twice"
  )
  # Either alone is fine.
  expect_silent(suppressMessages(
    sweep_models(exact_data(), model_set("hn"), truncation = 400, left = 20)
  ))
})

test_that("binned data sweeps through the binned likelihood", {
  d <- exact_data()
  d$distbegin <- floor(d$distance / 100) * 100
  d$distend <- d$distbegin + 100
  d$distance <- NULL

  top <- max(d$distend)
  sw <- sweep_fx(d, model_set("hn"), truncation = top)
  expect_true(sw$settings$binned)
  expect_equal(sw$settings$breaks, seq(0, top, by = 100))
  expect_true(all(selection_table(sw)$converged))

  # Bins that stop short of the truncation are unaccounted effort.
  expect_error(sweep_models(d, model_set("hn"), truncation = top + 100),
               "unaccounted effort")
})

test_that("a model that fails to fit is recorded, not fatal", {
  sw <- suppressWarnings(sweep_fx(
    exact_data(n = 12, seed = 3),
    model_set("hn", adjustment = "cos", order = 8),
    truncation = 400
  ))
  tab <- sw$table
  expect_equal(nrow(tab), 1)
  # Whether it converges is data-dependent; what matters is that the sweep
  # returned a table either way rather than stopping.
  expect_true(is.logical(tab$converged))
})

test_that("converged_only controls what the table shows", {
  sw <- sweep_fx(exact_data(), model_set(c("hn", "hr")), truncation = 400)
  expect_equal(nrow(selection_table(sw, converged_only = FALSE)), nrow(sw$table))
})

test_that("the sweep reports what it did and what it assumes", {
  # Include rows beyond the truncation and one with no distance, so there is
  # something to report as dropped.
  d <- exact_data()
  d <- rbind(d, data.frame(object = 9001:9003, distance = c(450, 800, NA)))
  msg <- tryCatch(
    sweep_models(d, model_set("hn"), truncation = 400),
    message = conditionMessage
  )
  expect_match(msg, "Fitted 1 of 1")
  expect_match(msg, "dropped")
  expect_match(msg, "g\\(0\\) = 1 is assumed")
})

test_that("printing shows the ranking and the caveat", {
  sw <- sweep_fx(exact_data(), model_set(c("hn", "hr")), truncation = 400)
  out <- paste(utils::capture.output(print(sw)), collapse = "\n")
  expect_match(out, "dsfit_sweep")
  expect_match(out, "esw")
  expect_match(out, "g\\(0\\) = 1 is assumed")
})
