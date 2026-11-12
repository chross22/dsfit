test_that("the default set is the two usual keys", {
  ms <- model_set()
  expect_setequal(ms$key, c("hn", "hr"))
  expect_true(all(is.na(ms$adjustment)))
  expect_equal(nrow(ms), 2)
})

test_that("keys, adjustments and formulas expand into a grid", {
  ms <- model_set(key = c("hn", "hr"), adjustment = "cos", order = 2:3,
                  formula = list(~1, ~beaufort))
  expect_equal(nrow(ms), 2 * 1 * 2 * 2)
  expect_setequal(unique(ms$order), 2:3)
  expect_setequal(unique(ms$formula), c("~1", "~beaufort"))
})

test_that("covariates switch cds for mcds", {
  plain <- model_set(key = "hn")
  cov <- model_set(key = "hn", formula = ~beaufort)

  expect_match(plain$dsmodel, "^~cds\\(")
  expect_match(cov$dsmodel, "^~mcds\\(")
  expect_match(cov$dsmodel, "formula = ~beaufort")
})

test_that("adjustment terms reach the dsmodel string", {
  ms <- model_set(key = "hn", adjustment = "herm", order = 4)
  expect_match(ms$dsmodel, 'adj.series = "herm"')
  expect_match(ms$dsmodel, "adj.order = 4")
  expect_equal(ms$model_id, "hn+herm(4)")
})

test_that("adjustment orders are cumulative, as in Distance::ds()", {
  # order is the highest order; every term of the series below it comes too.
  ms <- model_set(key = "hn", adjustment = "cos", order = 2:4)
  expect_equal(ms$model_id, c("hn+cos(2)", "hn+cos(2,3)", "hn+cos(2,3,4)"))
  expect_match(ms$dsmodel[3], "adj.order = c\\(2, 3, 4\\)")
  expect_equal(ms$order, 2:4)

  # The series depend on the key and the adjustment: Distance's
  # get_adj_orders().
  expect_equal(model_set("unif", adjustment = "cos", order = 3)$model_id,
               "unif+cos(1,2,3)")
  expect_equal(model_set("hr", adjustment = "herm", order = 8)$model_id,
               "hr+herm(4,6,8)")
  expect_equal(model_set("unif", adjustment = "poly", order = 4)$model_id,
               "unif+poly(2,4)")
})

test_that("an order outside its series is dropped, with a reason", {
  # Cosine order 1 on a half-normal, and an odd Hermite order, failed inside
  # mrds before. They are not models it can fit.
  expect_warning(ms <- model_set(c("hn", "unif"), adjustment = "cos", order = 1),
                 "outside their series: hn\\+cos1")
  expect_equal(ms$model_id, "unif+cos(1)")
  expect_warning(ms <- model_set("hn", adjustment = "herm", order = c(4, 5)),
                 "hn\\+herm5")
  expect_equal(ms$model_id, "hn+herm(4)")
  expect_error(model_set("hn", adjustment = "cos", order = 2.5), "whole numbers")
})

test_that("an adjustment with no order defaults rather than erroring", {
  # One term: the first of the series, which depends on the key.
  ms <- model_set(key = c("hn", "unif"), adjustment = c("cos", "herm"))
  expect_setequal(ms$model_id, c("hn+cos(2)", "unif+cos(1)", "hn+herm(4)",
                                 "unif+herm(2)"))
})

test_that("the gamma key is available, unlike in Distance::ds()", {
  ms <- model_set(key = c("hn", "gamma"))
  expect_true("gamma" %in% ms$key)
  expect_match(ms$dsmodel[ms$key == "gamma"], 'key = "gamma"')
})

test_that("a uniform key with no adjustment is dropped, with a reason", {
  # A flat detection function is a strip transect, not a model to fit.
  expect_warning(ms <- model_set(key = c("hn", "unif")), "strip transect")
  expect_setequal(ms$key, "hn")

  # With an adjustment it is a real candidate and is kept.
  expect_silent(keep <- model_set(key = "unif", adjustment = "cos"))
  expect_equal(keep$key, "unif")
})

test_that("a uniform key with covariates is dropped, with a reason", {
  # Covariates act on the scale, which the uniform key lacks. mrds ignores the
  # formula, so the candidate would duplicate unif+cos1 under a covariate label.
  expect_warning(
    ms <- model_set(key = c("hn", "unif"), adjustment = "cos", order = 2,
                    formula = list(~1, ~bf)),
    "does not have"
  )
  expect_setequal(ms$model_id, c("hn+cos(2)", "hn+cos(2) bf", "unif+cos(1,2)"))

  # Without covariates it stays.
  expect_silent(model_set(key = "unif", adjustment = "cos", formula = ~1))
  expect_error(
    suppressWarnings(model_set(key = "unif", adjustment = "cos", formula = ~bf)),
    "No candidate models"
  )
})

test_that("nothing is left is an error, not an empty tibble", {
  expect_error(suppressWarnings(model_set(key = "unif")), "No candidate models")
})

test_that("unknown keys and adjustments are refused by name", {
  expect_error(model_set(key = "hn2"), "hn2")
  expect_error(model_set(key = "hn2"), "\"gamma\"")
  expect_error(model_set(key = "hn", adjustment = "spline"), "spline")
  expect_error(model_set(key = "hn", formula = "beaufort"), "formula")
})

test_that("model ids are unique and readable", {
  ms <- model_set(key = c("hn", "hr"), adjustment = c("cos", "poly"),
                  order = 4, formula = list(~1, ~x))
  expect_equal(anyDuplicated(ms$model_id), 0L)
  expect_true(all(nchar(ms$model_id) > 0))
})

test_that("duplicate candidates are dropped, with a reason", {
  expect_warning(ms <- model_set(key = c("hn", "hn"), formula = list(~1, ~1)),
                 "duplicate candidate\\(s\\): hn")
  expect_equal(ms$model_id, "hn")
  expect_warning(model_set("hn", "cos", order = c(2, 2)), "hn\\+cos\\(2\\)")
})

test_that("two-sided formulas are refused", {
  expect_error(model_set("hn", formula = distance ~ x), "one-sided")
})

test_that("factor keys and adjustments come out as character", {
  ms <- model_set(key = factor(c("hn", "hr")), adjustment = factor("cos"))
  expect_type(ms$key, "character")
  expect_type(ms$adjustment, "character")
})
