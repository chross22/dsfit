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

# A flatfile as `distsamp::detection_data(include_circling = TRUE)` leaves it:
# mostly angle distances, with some from circling.
with_sources <- function(sources) {
  d <- data.frame(object = seq_along(sources),
                  distance = seq(10, 300, length.out = length(sources)))
  d$distance_source <- sources
  d
}

test_that("break-off distances are refused: they are not perpendicular", {
  d <- with_sources(c(rep("angle", 8), "break_off", "break_off"))
  expect_error(prepare_distance_data(d, truncation = 400),
               "2 distance\\(s\\) have `distance_source == \"break_off\"`")
  expect_error(prepare_distance_data(d, truncation = 400), "not perpendicular")
})

test_that("circling distances are admitted with a warning that counts them", {
  d <- with_sources(c(rep("angle", 8), "circling", NA))
  expect_warning(p <- prepare_distance_data(d, truncation = 400),
                 "1 distance\\(s\\) have `distance_source == \"circling\"`")
  expect_equal(nrow(p$data), 10)
})

test_that("a circling row with no distance is dropped, not warned about", {
  d <- with_sources(c(rep("angle", 8), "break_off", "circling"))
  d$distance[9:10] <- NA
  expect_silent(p <- prepare_distance_data(d, truncation = 400))
  expect_equal(p$dropped[["no_distance"]], 2)
})

test_that("ordinary distance sources pass without comment", {
  d <- with_sources(c(rep("angle", 5), rep("exact", 5)))
  expect_silent(prepare_distance_data(d, truncation = 400))
})

binned_with_empty_distance <- function() {
  d <- binned_data()
  d$distance <- NA_real_
  d
}

test_that("a binned flatfile with an all-NA distance column is fitted binned", {
  # The flatfile shape: all three columns present. mrds read the NA column as
  # the distances and had "no data to analyze" for every candidate.
  d <- binned_with_empty_distance()
  p <- prepare_distance_data(d, truncation = max(d$distend))
  expect_true(p$binned)
  expect_false("distance" %in% names(p$data))
})

test_that("empty bin columns on exact data are dropped", {
  d <- exact_data()
  d$distbegin <- NA_real_
  d$distend <- NA_real_
  p <- prepare_distance_data(d, truncation = 400)
  expect_false(any(c("distbegin", "distend") %in% names(p$data)))
})

test_that("distances must be numeric and not negative", {
  d <- data.frame(object = 1:4, distance = c("50", "120", "500", "30"))
  expect_error(prepare_distance_data(d, truncation = 400), "must be numeric")
  d$distance <- factor(d$distance)
  expect_error(prepare_distance_data(d, truncation = 400), "as.character")
  d$distance <- c(-50, 120, 300, 30)
  expect_error(prepare_distance_data(d, truncation = 400), "1 negative")
})

test_that("truncation must be finite", {
  expect_error(prepare_distance_data(exact_data(), truncation = Inf))
})

test_that("duplicated object IDs are refused before mrds sees them", {
  d <- exact_data()
  d$object[1:3] <- 1L
  expect_error(prepare_distance_data(d, truncation = 400), "1 value\\(s\\) repeat")
})

test_that("the no-distances message says what is actually wrong", {
  expect_error(prepare_distance_data(exact_data()[0, ], truncation = 400),
               "no rows")
  expect_error(prepare_distance_data(data.frame(object = 1:3, distance = NA),
                                     truncation = 400),
               "every value is missing")
})

test_that("bins converted between units still tile", {
  d <- data.frame(object = 1:4, distbegin = c(0, .1, .1 * 3, .3),
                  distend = c(.1, .1 * 3, .6, .6))
  p <- prepare_distance_data(d, truncation = .6)
  expect_equal(length(p$breaks), 4)
})

test_that("supplied breaks must be sorted and contain the recorded bins", {
  d <- data.frame(object = 1:6, distbegin = rep(c(0, 100, 200), 2),
                  distend = rep(c(100, 200, 300), 2))
  expect_error(prepare_distance_data(d, truncation = 300,
                                     breaks = c(300, 0, 100, 200)),
               "strictly increasing")
  expect_error(prepare_distance_data(d, truncation = 300,
                                     breaks = c(0, 150, 300)),
               "does not contain the bin edge\\(s\\) 100, 200")
  expect_silent(prepare_distance_data(d, truncation = 300,
                                      breaks = c(0, 100, 200, 300)))
})

test_that("breaks given with exact distances are ignored, with a warning", {
  expect_warning(prepare_distance_data(exact_data(), truncation = 400,
                                       breaks = c(0, 200, 400)), "ignored")
})
