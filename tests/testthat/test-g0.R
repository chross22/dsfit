avail <- function(value = 0.5, se = 0.04, key = NA_character_) {
  data.frame(key = key, component = "availability", value = value, se = se)
}
percep <- function(value = 0.7, se = 0.06, key = NA_character_) {
  data.frame(key = key, component = "perception", value = value, se = se)
}

test_that("components multiply and their squared CVs add", {
  g <- g0(avail(0.5, 0.04), percep(0.7, 0.06))

  expect_s3_class(g, "dsfit_g0")
  expect_equal(g$table$value, 0.5 * 0.7)

  # CV(g0)^2 = sum of component CV^2, under independence.
  expected_cv <- sqrt((0.04 / 0.5)^2 + (0.06 / 0.7)^2)
  expect_equal(g$table$cv, expected_cv)
  expect_equal(g$table$se, g$table$value * expected_cv)
})

test_that("the least precise component sets the floor", {
  # A very precise availability cannot rescue a vague perception.
  vague <- g0(avail(0.5, 0.001), percep(0.7, 0.21))
  expect_gt(vague$table$cv, 0.29)

  # Tightening the already-precise one barely moves it.
  tighter <- g0(avail(0.5, 0.0001), percep(0.7, 0.21))
  expect_equal(vague$table$cv, tighter$table$cv, tolerance = 0.01)
})

test_that("a correction without a standard error is refused", {
  no_se <- data.frame(key = NA_character_, component = "availability",
                      value = 0.5, se = NA_real_)
  expect_error(suppressWarnings(g0(no_se)), "without one reports a precision")

  # availability() with no interval SEs produces exactly that, so the two
  # rules meet: it will not invent a precision, and this will not accept none.
  a <- availability(surface = 60, dive = 240, window = 24)
  expect_error(suppressWarnings(g0(a)), "No standard error")
})

test_that("an absent component is warned about by name, not assumed silently", {
  expect_warning(g0(avail()), "\"perception\"")
  expect_warning(g0(percep()), "\"availability\"")
  # Both present, nothing to warn about.
  expect_silent(g0(avail(), percep()))
})

test_that("there is no empty default", {
  expect_error(g0(), "no default")
})

test_that("components must be named", {
  unnamed <- data.frame(key = NA_character_, component = NA_character_,
                        value = 0.5, se = 0.04)
  expect_error(suppressWarnings(g0(unnamed)), "needs a name")

  blank <- data.frame(key = NA_character_, component = "  ",
                      value = 0.5, se = 0.04)
  expect_error(suppressWarnings(g0(blank)), "needs a name")
})

test_that("impossible component values are refused", {
  expect_error(suppressWarnings(g0(avail(1.2, 0.04))), "must be in \\(0, 1\\]")
  expect_error(suppressWarnings(g0(avail(0, 0.04))), "must be in \\(0, 1\\]")
  expect_error(suppressWarnings(g0(avail(0.5, -1))), "cannot be negative")
  # Certain detection is a legitimate value, if an unusual one.
  expect_equal(suppressWarnings(g0(avail(1, 0.0)))$table$value, 1)
})

test_that("a table missing columns says which", {
  expect_error(g0(data.frame(component = "availability", value = 0.5)),
               "missing: key, se")
  expect_error(g0("not a table"), "not a data frame")
})

test_that("an unkeyed component broadcasts across the keyed ones", {
  months <- c("Jan", "Feb", "Mar")
  a <- avail(c(0.3, 0.5, 0.8), c(0.03, 0.04, 0.05), key = months)
  g <- g0(a, percep(0.7, 0.06))

  expect_equal(nrow(g$table), 3)
  expect_equal(g$table$key, months)
  expect_equal(g$table$value, c(0.3, 0.5, 0.8) * 0.7)
})

test_that("components keyed on different things are refused, not joined", {
  # Availability by month, perception by year: the honest answer is a value
  # per month-year, and pairing them positionally would be silent nonsense.
  a <- avail(c(0.3, 0.5, 0.8), 0.04, key = c("Jan", "Feb", "Mar"))
  p <- percep(c(0.6, 0.7), 0.06, key = c("2015", "2016"))

  expect_error(g0(a, p), "no correct join")
  expect_error(g0(a, p), "month-year")
})

test_that("half-keyed and duplicated components are refused", {
  half <- avail(c(0.3, 0.5), 0.04, key = c("Jan", NA))
  expect_error(suppressWarnings(g0(half)), "mixes keyed and unkeyed")

  dup <- avail(c(0.3, 0.5), 0.04, key = c("Jan", "Jan"))
  expect_error(suppressWarnings(g0(dup)), "more than one row for the same key")
})

test_that("the components survive in the object, so it can be taken apart", {
  g <- g0(avail(0.5, 0.04), percep(0.7, 0.06))
  expect_equal(nrow(g$components), 2)
  expect_setequal(g$components$component, c("availability", "perception"))
  expect_equal(g$components$value, c(0.5, 0.7))
})

test_that("printing shows what was applied and what was assumed", {
  out <- paste(utils::capture.output(print(g0(avail(), percep()))),
               collapse = "\n")
  expect_match(out, "dsfit_g0")
  expect_match(out, "availability, perception")
  expect_match(out, "abundance step")

  # A missing component is visible in the printout, not just at construction.
  partial <- paste(utils::capture.output(
    print(suppressWarnings(g0(avail())))
  ), collapse = "\n")
  expect_match(partial, "assumed 1:\\s+perception")
})

test_that("availability() feeds g0() directly", {
  a <- availability(
    surface = example_dive_intervals$surface,
    dive = example_dive_intervals$dive,
    window = view_window_aerial(trackline = 51.22, speed = 51.4, slope = 0.03),
    se_surface = example_dive_intervals$se_surface,
    se_dive = example_dive_intervals$se_dive,
    key = as.character(example_dive_intervals$month)
  )
  g <- g0(a, percep(0.68, 0.09))

  expect_equal(nrow(g$table), nrow(example_dive_intervals))
  # Every g(0) is below its availability, because perception can only lower it.
  expect_true(all(g$table$value < a$value))
  # And the combined CV exceeds either component's.
  expect_true(all(g$table$cv > a$se / a$value))
})
