test_that("point distances are kept and counted", {
  d <- data.frame(object = 1:5, distance = c(10, 50, 120, NA, 900))
  p <- prepare_distance_data(d, truncation = 400)

  expect_false(p$binned)
  expect_equal(nrow(p$data), 3)
  # One missing (a sample with no detection) and one beyond truncation.
  expect_equal(p$n_dropped, 2)
})

test_that("left truncation removes the inner distances", {
  d <- data.frame(object = 1:4, distance = c(10, 100, 200, 300))
  p <- prepare_distance_data(d, truncation = 400, left = 150)
  expect_equal(nrow(p$data), 2)
  expect_true(all(p$data$distance >= 150))

  expect_error(prepare_distance_data(d, truncation = 100, left = 200),
               "smaller than")
})

test_that("point and interval distances together are refused", {
  d <- data.frame(object = 1:4, distance = c(10, NA, 30, NA),
                  distbegin = c(NA, 0, NA, 100), distend = c(NA, 50, NA, 200))
  expect_error(prepare_distance_data(d, truncation = 400), "cannot rank both")
  expect_error(prepare_distance_data(d, truncation = 400), "distance-source")
})

test_that("data with no distances at all is refused", {
  expect_error(
    prepare_distance_data(data.frame(object = 1:3), truncation = 400),
    "No distances found"
  )
  expect_error(
    prepare_distance_data(data.frame(object = 1:3, distance = NA), truncation = 400),
    "No distances found"
  )
})

test_that("binned data derives its cutpoints", {
  d <- data.frame(
    object = 1:6,
    distbegin = c(0, 0, 100, 100, 200, 200),
    distend   = c(100, 100, 200, 200, 400, 400)
  )
  p <- prepare_distance_data(d, truncation = 400)
  expect_true(p$binned)
  expect_equal(p$breaks, c(0, 100, 200, 400))
})

test_that("an unbounded top bin is refused rather than fitted", {
  d <- data.frame(object = 1:4, distbegin = c(0, 100, 200, 400),
                  distend = c(100, 200, 400, Inf))
  expect_error(prepare_distance_data(d, truncation = 1e6), "unbounded top bin")
})

test_that("bins that do not tile are refused", {
  # A gap between 100 and 150 means these intervals do not define one set of
  # cutpoints, and deriving breaks anyway would misallocate detections.
  d <- data.frame(object = 1:4, distbegin = c(0, 0, 150, 150),
                  distend = c(100, 100, 300, 300))
  expect_error(prepare_distance_data(d, truncation = 400), "do not tile")

  # Explicit breaks let the caller override, so long as they reach the
  # truncation.
  expect_silent(
    prepare_distance_data(d, truncation = 300, breaks = c(0, 100, 150, 300))
  )
})

test_that("truncating everything away is an error, not an empty fit", {
  d <- data.frame(object = 1:3, distance = c(500, 600, 700))
  expect_error(prepare_distance_data(d, truncation = 100), "No detections left")
})

test_that("an object column is supplied when absent", {
  d <- data.frame(distance = c(10, 20))
  expect_true("object" %in% names(prepare_distance_data(d, truncation = 100)$data))
})

test_that("the tiling error names the interval before the gap", {
  d <- data.frame(object = 1:5, distbegin = c(0, 100, 250, 300, 400),
                  distend = c(100, 200, 300, 400, 500))
  expect_error(prepare_distance_data(d, truncation = 500),
               "interval 2 ends at 200 and the next begins at 250")
})

test_that("bins that start beyond zero are refused without `left`", {
  # mrds integrates from zero, so 0-100 would be effort with no bin. This is
  # what derived breaks give when the innermost bin had no detections.
  d <- data.frame(object = 1:6, distbegin = rep(c(100, 200, 300), 2),
                  distend = rep(c(200, 300, 400), 2))
  expect_error(prepare_distance_data(d, truncation = 400),
               "lowest bin begins at 100")

  # Either remedy the message names works.
  expect_silent(prepare_distance_data(d, truncation = 400, left = 100))
  p <- prepare_distance_data(d, truncation = 400,
                             breaks = c(0, 100, 200, 300, 400))
  expect_equal(p$breaks, c(0, 100, 200, 300, 400))
})

test_that("a `left` that falls inside a bin is refused", {
  # The 0-100 bin is dropped for beginning inside `left`, which would leave
  # 50-100 searched but unbinned.
  d <- data.frame(object = 1:6, distbegin = rep(c(0, 100, 200), 2),
                  distend = rep(c(100, 200, 300), 2))
  expect_error(prepare_distance_data(d, truncation = 300, left = 50),
               "lowest bin begins at 100 but the fit starts at 50")
  expect_silent(prepare_distance_data(d, truncation = 300, left = 100))
})
