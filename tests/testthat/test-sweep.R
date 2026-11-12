skip_if_not_installed("mrds")

# `exact_data()` and `binned_data()` come from helper-fixtures.R.
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
  for (nm in c("aic", "delta_aic", "p", "p_se", "p_cv", "esw", "cvm_p",
               "chisq_p")) {
    expect_true(nm %in% names(tab), info = nm)
  }
  # Effective strip half-width is p times the truncation width.
  expect_equal(tab$esw, tab$p * 400)
  expect_equal(tab$p_cv, tab$p_se / tab$p)

  # Numeric columns come out unnamed, so pulling one out gives a bare vector.
  expect_null(names(tab$aic))
  expect_null(names(tab$esw))
})

test_that("goodness of fit follows the data: CvM exact, chi-square binned", {
  exact <- selection_table(sweep_fx(exact_data(), model_set("hn"),
                                    truncation = 400))
  expect_false(is.na(exact$cvm_p))

  # A binned fit has no empirical distribution for Cramer-von Mises to test,
  # and mrds does not compute one. Without chisq_p the goodness-of-fit column
  # would be silently empty for every binned sweep.
  d <- binned_data()
  binned <- selection_table(sweep_fx(d, model_set("hn"),
                                     truncation = max(d$distend)))
  expect_true(is.na(binned$cvm_p))
  expect_false(is.na(binned$chisq_p))
})

test_that("printing shows whichever goodness-of-fit test applies", {
  d <- binned_data()
  out <- paste(utils::capture.output(print(
    sweep_fx(d, model_set("hn"), truncation = max(d$distend))
  )), collapse = "\n")
  expect_match(out, "chisq_p")
  expect_no_match(out, "CvM_p")
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
  d <- binned_data()

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
  # A covariate that is not in the data fails inside mrds every time, so the
  # failure is certain rather than dependent on the optimiser's start.
  sw <- suppressWarnings(sweep_fx(
    exact_data(), model_set("hn", formula = list(~1, ~absent)),
    truncation = 400
  ))
  tab <- sw$table
  expect_equal(nrow(tab), 2)
  expect_equal(tab$converged[tab$model_id == "hn absent"], FALSE)
  expect_true(tab$converged[tab$model_id == "hn"])
  # The reason is recorded under the model's ID, and is not blank.
  expect_named(sw$settings$failures, "hn absent")
  expect_true(nzchar(sw$settings$failures[["hn absent"]]))
})

test_that("a failed fit does not slide the other models onto wrong rows", {
  # The failures come first, so every model that fits sits after one: a
  # `fits[[i]] <- NULL` shortens the list instead of emptying it, and every
  # later model inherits the previous one's fit. A covariate that is not in
  # the data fails inside mrds every time, which makes the failure certain.
  sw <- suppressWarnings(sweep_fx(
    exact_data(),
    model_set(key = c("hn", "hr"), adjustment = "cos", order = 2:3,
              formula = list(~absent, ~1)),
    truncation = 400
  ))
  tab <- sw$table

  expect_equal(nrow(tab), 8)
  expect_equal(length(sw$fits), 8)
  expect_equal(sum(!tab$converged), 4)
  expect_true(all(grepl("absent", tab$model_id[!tab$converged])))
  expect_true(all(is.na(tab$aic[!tab$converged])))

  # Each row's AIC is the AIC of the fit stored under that row's model_id.
  for (id in tab$model_id[tab$converged]) {
    expect_equal(tab$aic[tab$model_id == id], sw$fits[[id]]$criterion,
                 info = id)
  }
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

test_that("a fit mrds returns unconverged is a failure, not a crash", {
  # mrds warns rather than errors when the optimiser does not converge, and
  # hands back an object whose summary has no average.p. That used to stop
  # the whole sweep inside the selection table. mrds draws random starting
  # values, so whether three points converge depends on the seed; this one
  # does not.
  d <- data.frame(object = 1:3, distance = c(1, 2, 3))
  set.seed(1)
  sw <- suppressWarnings(sweep_fx(d, model_set("hr"), truncation = 400))
  tab <- sw$table

  expect_false(tab$converged)
  expect_true(is.na(tab$aic))
  expect_true(is.na(tab$delta_aic))
  expect_null(sw$fits[["hr"]])
  expect_equal(unname(sw$settings$failures), "did not converge")
})

test_that("one unconverged model leaves the others in the table", {
  d <- data.frame(object = 1:3, distance = c(1, 2, 3))
  set.seed(1)
  sw <- suppressWarnings(sweep_fx(d, model_set(c("hn", "hr")),
                                  truncation = 400))
  tab <- sw$table
  expect_equal(nrow(tab), 2)
  expect_equal(tab$converged[tab$model_id == "hr"], FALSE)
})

test_that("a cumulative adjustment fits the model mrds gives for those orders", {
  d <- exact_data()
  sw <- sweep_fx(d, model_set("hn", adjustment = "cos", order = 3),
                 truncation = 400)
  direct <- mrds::ddf(
    dsmodel = ~cds(key = "hn", formula = ~1, adj.series = "cos",
                   adj.order = c(2, 3)),
    data = d, method = "ds", meta.data = list(width = 400)
  )
  expect_equal(sw$table$aic, direct$criterion)
  expect_equal(sw$table$n_par, 3)
})

test_that("esw is p times the searched width when left truncation is set", {
  # average.p is the mean of g(x) over (left, width), so the strip it scales
  # is width - left. Multiplying by width alone overstates esw, and density
  # with it.
  sw <- sweep_fx(exact_data(), model_set("hn"), truncation = 400, left = 50)
  tab <- selection_table(sw)
  expect_equal(tab$esw, tab$p * (400 - 50))
})

test_that("adjustment models are fitted monotone, as Distance::ds() fits them", {
  # Unconstrained, hr+cos(2) put g(x) at 1.09 here, and its AIC and p differed
  # from ds() with its defaults (3399.27 against 3401.50).
  set.seed(1)
  x <- abs(stats::rnorm(900, 0, 150))
  x <- x[x < 500][1:300]
  sw <- suppressWarnings(sweep_fx(data.frame(distance = x),
                                  model_set("hr", adjustment = "cos", order = 2),
                                  truncation = 400))
  expect_equal(sw$table$aic, 3401.498, tolerance = 1e-6)
  expect_true(sw$fits[[1]]$meta.data$mono.strict)
  est <- effect_estimates_ddf(sw$fits[[1]], n = 501)
  expect_lte(max(est$.estimate), 1 + 1e-8)
})

test_that("keys alone, covariate models and gamma are not constrained", {
  m <- model_set(c("hn", "gamma"), adjustment = "cos", formula = list(~1, ~bf))
  for (i in seq_len(nrow(m))) {
    want <- m$key[i] != "gamma" && m$formula[i] == "~1"
    expect_equal(monotonicity(m[i, ])$mono, want, info = m$model_id[i])
  }
  expect_false(monotonicity(model_set("hn")[1, ])$mono)
})

test_that("a failure reason is the first non-empty line of the error", {
  expect_equal(first_line("\nSome values of object are duplicates\nmore"),
               "Some values of object are duplicates")
  expect_equal(first_line("plain"), "plain")
})

test_that("models must be a model set", {
  expect_error(sweep_fx(exact_data(), data.frame(key = "hn"), truncation = 400),
               "built by `model_set\\(\\)`")
  expect_error(sweep_fx(exact_data(), model_set()[0, ], truncation = 400),
               "no candidates")
})
