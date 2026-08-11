verdict <- function(x, check) {
  x$table$supported[x$table$check == check]
}
detail <- function(x, check) {
  x$table$detail[x$table$check == check]
}

exact <- function(n = 200, ...) {
  set.seed(1)
  data.frame(object = seq_len(n), distance = abs(stats::rnorm(n, 0, 120)), ...)
}

test_that("exact distances admit Cramer-von Mises, binned do not", {
  e <- detection_structure(exact())
  expect_true(verdict(e, "test fit with Cramer-von Mises"))

  b <- detection_structure(binned_data())
  expect_false(verdict(b, "test fit with Cramer-von Mises"))
  expect_match(detail(b, "test fit with Cramer-von Mises"), "empirical")

  # Chi-square goes the other way: available for both, but weaker for exact
  # because mrds picks the cutpoints.
  expect_true(verdict(b, "test fit with chi-square"))
  expect_match(detail(b, "test fit with chi-square"), "survey's own bins")
  expect_match(detail(e, "test fit with chi-square"), "weaker")
})

test_that("a single-observer table cannot estimate perception, and says why", {
  s <- detection_structure(exact())
  expect_false(verdict(s, "estimate perception bias"))
  expect_match(detail(s, "estimate perception bias"), "no `observer`")
  # The reason matters more than the verdict: this is not fixable with more data.
  expect_match(detail(s, "estimate perception bias"), "at any sample size")
})

test_that("real double-observer structure is recognised", {
  io <- data.frame(
    object = rep(1:50, each = 2),
    distance = rep(abs(stats::rnorm(50, 0, 120)), each = 2),
    observer = rep(1:2, 50),
    detected = rep(1L, 100)
  )
  expect_true(verdict(detection_structure(io), "estimate perception bias"))
})

test_that("partial double-observer structure is neither yes nor no", {
  # `observer` without `detected` is the shape that most invites being treated
  # as double-observer data when it is not.
  half <- exact(50)
  half$observer <- rep(1:2, 25)
  s <- detection_structure(half)
  expect_true(is.na(verdict(s, "estimate perception bias")))
  expect_match(detail(s, "estimate perception bias"), "`detected` is missing")

  # And the reverse.
  other <- exact(50)
  other$detected <- 1L
  expect_true(is.na(verdict(detection_structure(other),
                            "estimate perception bias")))
})

test_that("one observer is not two", {
  one <- exact(50)
  one$observer <- 1L
  one$detected <- 1L
  s <- detection_structure(one)
  expect_false(verdict(s, "estimate perception bias"))
  expect_match(detail(s, "estimate perception bias"), "two independent teams")
})

test_that("observers who never share a sighting are flagged", {
  # Both columns present and two observers, but no object seen twice - so
  # nothing to mark and recapture.
  d <- exact(50)
  d$observer <- rep(1:2, 25)
  d$detected <- 1L
  s <- detection_structure(d)
  expect_true(is.na(verdict(s, "estimate perception bias")))
  expect_match(detail(s, "estimate perception bias"), "no `object` appears twice")
})

test_that("availability is never supported, for any table", {
  # Not a defect in a dataset: it cannot be estimated from distances at all.
  for (d in list(exact(), binned_data(), exact(20))) {
    s <- detection_structure(d)
    expect_false(verdict(s, "estimate availability"))
    expect_match(detail(s, "estimate availability"), "availability\\(\\)")
  }
})

test_that("sample size is reported against the usual guideline", {
  expect_match(detail(detection_structure(exact(200)),
                      "fit a detection function"), "at or above")
  expect_match(detail(detection_structure(exact(50)),
                      "fit a detection function"), "below the 60-80")
  expect_match(detail(detection_structure(exact(20)),
                      "fit a detection function"), "below the 40")
})

test_that("covariate candidates exclude the structural columns", {
  d <- exact(100)
  d$beaufort <- sample(0:4, 100, replace = TRUE)
  d$size <- 1L
  d$Region.Label <- "a"
  s <- detection_structure(d)

  expect_true(verdict(s, "fit covariate models"))
  expect_equal(s$summary$covariates, "beaufort")
  # `size`, `object` and the flatfile columns are structure, not covariates.
  expect_false(any(c("size", "object", "Region.Label") %in% s$summary$covariates))
})

test_that("a constant column is not a covariate candidate", {
  d <- exact(100)
  d$always_three <- 3
  expect_false(verdict(detection_structure(d), "fit covariate models"))
})

test_that("mixed and missing distances are reported, not fitted", {
  both <- exact(50)
  both$distbegin <- 0
  both$distend <- 100
  s <- detection_structure(both)
  expect_true(is.na(verdict(s, "fit a detection function")))
  expect_match(detail(s, "fit a detection function"), "different")

  none <- data.frame(object = 1:10, beaufort = 1:10)
  n <- detection_structure(none)
  expect_false(verdict(n, "fit a detection function"))
  expect_equal(n$summary$type, "none")
})

test_that("group size is reported because abundance needs it", {
  expect_false(verdict(detection_structure(exact()),
                       "carry group size to abundance"))
  with_size <- exact(50)
  with_size$size <- 2L
  expect_true(verdict(detection_structure(with_size),
                      "carry group size to abundance"))
})

test_that("it reports rather than enforces", {
  # prepare_distance_data() errors on data this only describes, which is the
  # division of labour: one is a gate, the other is a briefing.
  both <- exact(50)
  both$distbegin <- 0
  both$distend <- 100

  expect_s3_class(detection_structure(both), "dsfit_structure")
  expect_error(prepare_distance_data(both, truncation = 400))
})

test_that("printing groups by verdict and keeps the caveat", {
  d <- exact(200)
  d$beaufort <- sample(0:4, 200, replace = TRUE)
  out <- paste(utils::capture.output(print(detection_structure(d))),
               collapse = "\n")

  expect_match(out, "dsfit_structure")
  expect_match(out, "can:")
  expect_match(out, "cannot:")
  expect_match(out, "nearest detection")
  expect_match(out, "not a reason to do it")
})

test_that("bad input is refused", {
  expect_error(detection_structure("not a table"), "data frame")
})

test_that("NARWC identifier columns are not offered as detection covariates", {
  skip_if_not_installed("narwcr")

  # Before narwcr's vocabulary was consulted, this reported FILEID, EVENTNO,
  # LATITUDE and LONGITUDE as things to model detection on. A file identifier
  # and an event counter are not survey conditions.
  d <- exact(120)
  d$FILEID <- rep(c("A", "B"), 60)
  d$EVENTNO <- seq_len(120)
  d$LATITUDE <- seq(41, 43, length.out = 120)
  d$LONGITUDE <- seq(-71, -69, length.out = 120)
  d$SPECCODE <- rep(c("RIWH", "HUWH"), 60)
  d$BEAUFORT <- rep(0:4, 24)

  s <- detection_structure(d)
  expect_false(any(c("FILEID", "EVENTNO", "LATITUDE", "LONGITUDE") %in%
                     s$summary$covariates))
  expect_true("BEAUFORT" %in% s$summary$covariates)
})

test_that("a programme's spelling of a condition is still recognised", {
  skip_if_not_installed("narwcr")

  # narwcr's aliases carry SEASTATE, SEA_STATE and BFT onto BEAUFORT. Without
  # them each spelling would be an unrecognised column.
  d <- exact(120)
  d$SEASTATE <- rep(0:4, 24)
  s <- detection_structure(d)

  expect_true("SEASTATE" %in% s$summary$covariates)
  expect_match(detail(s, "fit covariate models"), "survey conditions: SEASTATE")
})

test_that("columns narwcr does not know are still offered", {
  # An unrecognised column may well be a real covariate, so the vocabulary
  # narrows the false positives without narrowing the true ones.
  d <- exact(120)
  d$boatname <- rep(c("x", "y"), 60)
  expect_true("boatname" %in% detection_structure(d)$summary$covariates)
})

test_that("classification works without narwcr installed", {
  # The fallback is the mrds and flatfile names alone, which is right for
  # non-NARWC data and merely less informed for NARWC data.
  vocab <- narwc_vocabulary()
  expect_type(vocab, "list")
  expect_named(vocab, c("structural", "conditions", "aliases"))

  d <- exact(120)
  d$beaufort <- rep(0:4, 24)
  d$size <- 1L
  s <- detection_structure(d)
  expect_true("beaufort" %in% s$summary$covariates)
  expect_false("size" %in% s$summary$covariates)
})
