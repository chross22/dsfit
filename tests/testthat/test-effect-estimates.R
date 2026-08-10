skip_if_not_installed("mrds")

fit_hn <- function(n = 200, seed = 1, ...) {
  set.seed(seed)
  d <- data.frame(object = seq_len(n), distance = abs(stats::rnorm(n, 0, 120)))
  d <- d[d$distance < 400, , drop = FALSE]
  mrds::ddf(dsmodel = ~cds(key = "hn", formula = ~1), data = d,
            method = "ds", meta.data = list(width = 400, ...))
}

fit_cov <- function(n = 200, seed = 1) {
  set.seed(seed)
  d <- data.frame(object = seq_len(n), distance = abs(stats::rnorm(n, 0, 120)))
  d <- d[d$distance < 400, , drop = FALSE]
  d$bf <- stats::runif(nrow(d), 0, 5)
  mrds::ddf(dsmodel = ~mcds(key = "hn", formula = ~bf), data = d,
            method = "ds", meta.data = list(width = 400))
}

test_that("the frame matches the contract fancyfx documents", {
  est <- effect_estimates_ddf(fit_hn())
  expect_named(est, c(".x", ".estimate", ".lower", ".upper"))
  expect_equal(nrow(est), 100)
  expect_type(attr(est, "quantity"), "character")
  expect_match(attr(est, "quantity"), "Detection probability")
})

test_that("the curve is a detection function", {
  est <- effect_estimates_ddf(fit_hn(), n = 200)

  # Half-normal: certain on the track-line, monotone decreasing, a probability.
  expect_equal(est$.estimate[1], 1)
  expect_false(is.unsorted(rev(est$.estimate)))
  expect_true(all(est$.estimate >= 0 & est$.estimate <= 1))

  # The grid spans the truncation width.
  expect_equal(range(est$.x), c(0, 400))
})

trapz <- function(x, y) {
  sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}

test_that("for an intercept-only model the curve integrates to average.p", {
  # Mean over animals of the integral is the integral of the mean, and with a
  # constant detection probability that is also what mrds reports.
  f <- fit_hn()
  est <- effect_estimates_ddf(f, n = 2000)
  p_bar <- trapz(est$.x, est$.estimate) / diff(range(est$.x))
  expect_equal(p_bar, summary(f)$average.p, tolerance = 1e-4)
})

test_that("for a covariate model it integrates to the arithmetic mean, not average.p", {
  # mrds reports average.p as a Horvitz-Thompson mean, n / sum(1/p), which
  # weights each animal by the inverse of its own detection probability. This
  # curve is the plain mean over the animals as observed. They coincide only
  # when p is constant, and the difference is real rather than numerical.
  f <- fit_cov()
  est <- effect_estimates_ddf(f, n = 2000)
  p_bar <- trapz(est$.x, est$.estimate) / diff(range(est$.x))

  expect_equal(p_bar, mean(f$fitted), tolerance = 1e-4)
  expect_equal(summary(f)$average.p, length(f$fitted) / sum(1 / f$fitted),
               tolerance = 1e-8)
  expect_false(isTRUE(all.equal(p_bar, summary(f)$average.p, tolerance = 1e-5)))
})

test_that("intervals are bounded by the unit interval and widen for ci", {
  f <- fit_hn()
  se <- effect_estimates_ddf(f, interval = "se")
  ci <- effect_estimates_ddf(f, interval = "ci", level = 0.95)

  expect_true(all(se$.lower >= 0) && all(se$.upper <= 1))
  expect_true(all(ci$.lower <= se$.lower + 1e-12))
  expect_true(all(ci$.upper >= se$.upper - 1e-12))

  # A wider level is a wider interval.
  wide <- effect_estimates_ddf(f, interval = "ci", level = 0.99)
  expect_gte(sum(wide$.upper - wide$.lower), sum(ci$.upper - ci$.lower))
})

test_that("covariates are averaged over, and `at` fixes them", {
  f <- fit_cov()
  avg <- effect_estimates_ddf(f, n = 300)
  at3 <- effect_estimates_ddf(f, at = list(bf = 3), n = 300)

  expect_false(isTRUE(all.equal(avg$.estimate, at3$.estimate)))
  expect_equal(at3$.estimate[1], 1)

  # Fixing a covariate at two values gives two different curves.
  at1 <- effect_estimates_ddf(f, at = list(bf = 1), n = 300)
  expect_false(isTRUE(all.equal(at1$.estimate, at3$.estimate)))
})

test_that("`at` is validated against the model's own covariates", {
  f <- fit_cov()
  expect_error(effect_estimates_ddf(f, at = list(nope = 1)), "nope")
  expect_error(effect_estimates_ddf(f, at = list(3)), "named list")
  expect_error(effect_estimates_ddf(f, at = "bf"), "named list")
})

test_that("left truncation moves the start of the grid", {
  f <- fit_hn(left = 40)
  est <- effect_estimates_ddf(f)
  expect_equal(min(est$.x), 40)
  expect_equal(max(est$.x), 400)
})

test_that("only distance is supported, and the message says what is missing", {
  f <- fit_cov()
  expect_error(effect_estimates_ddf(f, var = "bf"), "must be \"distance\"")
  expect_error(effect_estimates_ddf(f, var = "bf"), "not built yet")
})

test_that("the generic's whole argument vocabulary is accepted", {
  f <- fit_hn()
  # fancyfx forwards its own unresolved defaults, so "auto" arrives routinely
  # and a length > 1 vector is normal. match.arg() would reject both.
  expect_s3_class(effect_estimates_ddf(f, interval = "auto"), "data.frame")
  expect_s3_class(effect_estimates_ddf(f, scale = "auto"), "data.frame")
  expect_s3_class(
    effect_estimates_ddf(f, scale = c("auto", "link", "response"),
                         interval = c("auto", "se", "ci")),
    "data.frame"
  )

  # "auto" resolves to the one-standard-error ribbon.
  expect_equal(effect_estimates_ddf(f, interval = "auto"),
               effect_estimates_ddf(f, interval = "se"))

  # A credible interval is not on offer from a frequentist fit.
  expect_error(effect_estimates_ddf(f, interval = "cri"), "no posterior")
  expect_error(effect_estimates_ddf(f, interval = "bootstrap"), "Unknown interval")
})

test_that("a fit from a sweep works without unwrapping", {
  set.seed(1)
  d <- data.frame(object = 1:200, distance = abs(stats::rnorm(200, 0, 120)))
  d <- d[d$distance < 400, ]
  sw <- suppressMessages(sweep_models(d, model_set("hn"), truncation = 400))
  expect_s3_class(effect_estimates_ddf(sw$fits[["hn"]]), "data.frame")
})

test_that("the fancyfx generic dispatches here", {
  skip_if_not_installed("fancyfx")
  est <- fancyfx::effect_estimates(fit_hn(), "distance")
  expect_named(est, c(".x", ".estimate", ".lower", ".upper"))
  expect_match(attr(est, "quantity"), "Detection probability")
})

test_that("plotEffects draws a detection function with a rug", {
  skip_if_not_installed("fancyfx")
  skip_if_not_installed("ggplot2")

  f <- fit_hn()
  p <- fancyfx::plotEffects(f, dat = f$data, var = "distance")
  expect_no_error(ggplot2::ggplot_build(p))
})
