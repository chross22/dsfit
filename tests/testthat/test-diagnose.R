report <- function(...) {
  paste(utils::capture.output(suppressWarnings(diagnose_sweep(...))),
        collapse = "\n")
}

survey <- function(n = 300, blind = 0, sd = 800, ...) {
  set.seed(4)
  d <- abs(stats::rnorm(n, 0, sd)) + blind
  data.frame(object = seq_len(n), distance = d, ...)
}

test_that("it reports rather than fits", {
  out <- report(survey(), model_set(c("hn", "hr")), truncation = 3000)

  expect_match(out, "dsfit sweep diagnosis")
  expect_match(out, "Nothing was fitted")
  expect_match(out, "2 candidates: hn, hr")
})

test_that("the return value is what was reached, for picking up from", {
  res <- suppressWarnings(
    diagnose_sweep(survey(), model_set("hn"), truncation = 3000)
  )
  expect_named(res, c("structure", "prepared", "models"))
  expect_s3_class(res$structure, "dsfit_structure")
  expect_true(is.data.frame(res$models))
})

test_that("truncation is required, as it is for a sweep", {
  expect_error(diagnose_sweep(survey()), "`truncation` is required")
})

test_that("a guard failure stops the run and says which guard", {
  both <- survey(50)
  both$distbegin <- 0
  both$distend <- 100

  out <- report(both, truncation = 400)
  expect_match(out, "FAIL")
  expect_match(out, "different likelihoods")
  expect_match(out, "Stopped")

  # And it hands back what it reached, so the structure is still inspectable.
  res <- suppressWarnings(diagnose_sweep(both, truncation = 400))
  expect_s3_class(res$structure, "dsfit_structure")
  expect_null(res$prepared)
})

test_that("attrition is broken down by reason", {
  d <- survey(200)
  d$distance[1:20] <- NA
  out <- report(d, model_set("hn"), truncation = 500, left = 50)

  expect_match(out, "20 with no distance")
  expect_match(out, "beyond the truncation")
  expect_match(out, "inside `left`")
})

test_that("a truncation that throws away most of the survey is flagged", {
  # Legal, and rarely intended - usually a units mismatch.
  out <- report(survey(300, sd = 2000), model_set("hn"), truncation = 300)
  expect_match(out, "WARN")
  expect_match(out, "% of the rows")
  expect_match(out, "units")
})

test_that("too few detections to fit is a failure, not a warning", {
  out <- report(survey(30), model_set("hn"), truncation = 5000)
  expect_match(out, "FAIL")
  expect_match(out, "shape rather than evidence")

  # And the middle band is a warning.
  mid <- report(survey(50), model_set("hn"), truncation = 5000)
  expect_match(mid, "WARN")
  expect_match(mid, "below the 60-80")
})

test_that("a covariate formula naming an absent column is caught", {
  # This otherwise fails inside mrds, with a message that does not name the
  # column that is missing.
  out <- report(survey(200), model_set("hn", formula = ~beaufort),
                truncation = 3000)
  expect_match(out, "FAIL")
  expect_match(out, "not in `data`: beaufort")
})

test_that("a covariate that cannot inform detection is flagged", {
  constant <- survey(200, beaufort = 2)
  out <- report(constant, model_set("hn", formula = ~beaufort),
                truncation = 3000)
  expect_match(out, "single value")

  gappy <- survey(200)
  gappy$beaufort <- c(rep(NA, 10), rep(1:4, length.out = 190))
  out2 <- report(gappy, model_set("hn", formula = ~beaufort),
                 truncation = 3000)
  expect_match(out2, "missing values")
})

test_that("an abrupt edge is not treated as handled just because gamma is there", {
  # The correction that matters. Gamma models a gradual reduction toward the
  # trackline; a geometric cutoff is a discontinuity, and no key function is
  # discontinuous. Calling it handled was false reassurance in exactly the case
  # the vignette uses to demonstrate every candidate failing its fit test.
  blind <- survey(300, blind = 100)

  with_gamma <- report(blind, model_set(c("hn", "gamma")), truncation = 3000)
  expect_match(with_gamma, "WARN")
  expect_match(with_gamma, "geometric edge")
  expect_match(with_gamma, "including gamma")
  expect_match(with_gamma, "fail its goodness-of-fit test")

  # Without gamma, the same edge, without the clause about gamma.
  without <- report(blind, model_set(c("hn", "hr")), truncation = 3000)
  expect_match(without, "geometric edge")
  expect_no_match(without, "including gamma")

  # Left truncation is the treatment that matches an edge, and removes it.
  with_left <- report(blind, model_set("hn"), truncation = 3000, left = 100)
  expect_match(with_left, "removed by left truncation")
  expect_no_match(with_left, "geometric edge")
})

test_that("onset shape separates an edge from a rise, without fitting", {
  set.seed(3)
  abrupt <- data.frame(distance = 100 + abs(stats::rnorm(1500, 0, 700)))
  expect_equal(onset_shape(abrupt, 100, 2500), "abrupt")

  x <- seq(100, 2500, length.out = 2000)
  gradual <- data.frame(
    distance = sample(x, 1500, replace = TRUE,
                      prob = stats::dnorm(x, mean = 900, sd = 500))
  )
  expect_equal(onset_shape(gradual, 100, 2500), "gradual")

  # Too few detections to judge the shape of anything.
  expect_equal(onset_shape(data.frame(distance = 1:10), 1, 100), "unknown")
})

test_that("a rise into the empty strip is what gamma is for", {
  set.seed(3)
  x <- seq(120, 3000, length.out = 3000)
  rising <- data.frame(
    object = 1:1500,
    distance = sample(x, 1500, replace = TRUE,
                      prob = stats::dnorm(x, mean = 1000, sd = 550))
  )
  out <- report(rising, model_set(c("hn", "gamma")), truncation = 3000)
  expect_match(out, "detections rise into the empty strip")
  expect_no_match(out, "geometric edge")
})

test_that("gamma and left together are still double counting", {
  out <- report(survey(300, blind = 100), model_set(c("hn", "gamma")),
                truncation = 3000, left = 100)
  expect_match(out, "blind spot twice")
})

test_that("absent structure is noted, ambiguous structure is warned about", {
  # A survey that had one observer team is not misconfigured, so saying so on
  # every dataset as a WARN would bury the cases that are. Half-present
  # double-observer structure is the opposite - it gets mistaken for the real
  # thing - so that one warns.
  plain <- report(survey(200), model_set("hn"), truncation = 3000)
  expect_match(plain, "note  estimate perception bias")
  expect_no_match(plain, "WARN  estimate perception bias")

  half <- survey(200)
  half$observer <- rep(1:2, 100)
  ambiguous <- report(half, model_set("hn"), truncation = 3000)
  expect_match(ambiguous, "WARN  estimate perception bias")
})

test_that("universal and merely-absent findings are not warned about", {
  # Availability is unsupported for every table and having no covariates is a
  # fact rather than a fault. Warning about either trains the reader to skim.
  out <- report(survey(200), model_set("hn"), truncation = 3000)
  expect_no_match(out, "estimate availability")
  expect_no_match(out, "WARN  fit covariate models")
})

test_that("an ordinary smallest distance is not mistaken for a blind spot", {
  # Every continuous distance has a minimum, so "nearest > 0" would fire on
  # every survey ever flown. Only an empty strip wide against the truncation
  # counts.
  ordinary <- report(survey(300), model_set(c("hn", "hr")), truncation = 3000)
  expect_no_match(ordinary, "reproduces a cliff")
})

test_that("the assumption no check can make for you is stated", {
  out <- report(survey(200), model_set("hn"), truncation = 3000)
  expect_match(out, "g\\(0\\) = 1")
  expect_match(out, "ranking looks untouched")
})

test_that("a clean survey reports no problems", {
  out <- report(survey(300, size = 2L), model_set(c("hn", "hr")),
                truncation = 3000, left = NULL)
  expect_match(out, "No problems found")
})
