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
  expect_equal(ms$model_id, "hn+herm4")
})

test_that("an adjustment with no order defaults rather than erroring", {
  ms <- model_set(key = "hn", adjustment = "cos")
  expect_equal(ms$order, 2L)
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
                  order = 2, formula = list(~1, ~x))
  expect_equal(anyDuplicated(ms$model_id), 0L)
  expect_true(all(nchar(ms$model_id) > 0))
})
