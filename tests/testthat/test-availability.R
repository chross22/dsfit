test_that("the two limits of the formula are what they should be", {
  s <- 60
  d <- 240

  # An instantaneous window is the plain proportion of time at the surface.
  a0 <- availability(surface = s, dive = d, window = 0)
  expect_equal(a0$value, s / (s + d))

  # Watch forever and everything surfaces eventually.
  a_inf <- availability(surface = s, dive = d, window = Inf)
  expect_equal(a_inf$value, 1)

  # In between, monotonically increasing in the window.
  w <- c(0, 10, 30, 60, 120, 600)
  a <- availability(surface = s, dive = d, window = w)$value
  expect_false(is.unsorted(a))
  expect_true(all(a >= s / (s + d) & a <= 1))
})

test_that("availability rises with surface time and falls with dive time", {
  base <- availability(surface = 60, dive = 240, window = 24)$value
  expect_gt(availability(surface = 90, dive = 240, window = 24)$value, base)
  expect_lt(availability(surface = 60, dive = 400, window = 24)$value, base)
})

test_that("the result is the shape a g(0) correction stacks from", {
  a <- availability(surface = c(45, 60, 90), dive = c(300, 240, 150),
                    window = 24, key = c("Feb", "Mar", "Apr"))

  expect_equal(nrow(a), 3)
  expect_equal(names(a), c("key", "component", "value", "se"))
  expect_equal(a$key, c("Feb", "Mar", "Apr"))
  expect_true(all(a$component == "availability"))
  # No standard errors were supplied, so none are invented.
  expect_true(all(is.na(a$se)))
})

test_that("a monthly range like Ganley's comes out of monthly inputs", {
  # The published finding is that availability varies by month between about
  # 0.27 and 0.85. The point of this test is not the numbers but that the
  # function spans that kind of range rather than hovering near one value —
  # a single representative figure is the error this exists to prevent.
  deep <- availability(surface = 30, dive = 600, window = 24)$value
  shallow <- availability(surface = 200, dive = 60, window = 24)$value

  expect_lt(deep, 0.35)
  expect_gt(shallow, 0.75)
})

test_that("standard errors propagate, and are not invented", {
  with_se <- availability(surface = 60, dive = 240, window = 24,
                          se_surface = 8, se_dive = 25)
  expect_false(is.na(with_se$se))
  expect_gt(with_se$se, 0)

  # More uncertain inputs, more uncertain output.
  wider <- availability(surface = 60, dive = 240, window = 24,
                        se_surface = 20, se_dive = 60)
  expect_gt(wider$se, with_se$se)

  # Certain inputs, certain output.
  exact <- availability(surface = 60, dive = 240, window = 24,
                        se_surface = 0, se_dive = 0)
  expect_equal(exact$se, 0)
})

test_that("the delta method agrees with the analytic derivative", {
  # a = (s + d(1 - exp(-w/d))) / (s + d). Differentiating with respect to s:
  #   da/ds = d * exp(-w/d) / (s + d)^2
  s <- 60; d <- 240; w <- 24
  se_s <- 8

  analytic <- (d * exp(-w / d) / (s + d)^2) * se_s
  numeric <- availability(surface = s, dive = d, window = w,
                          se_surface = se_s, se_dive = 0)$se
  expect_equal(numeric, analytic, tolerance = 1e-6)
})

test_that("half a variance is refused rather than reported as all of it", {
  expect_error(
    availability(surface = 60, dive = 240, window = 24, se_surface = 8),
    "both `se_surface` and `se_dive`"
  )
  expect_error(
    availability(surface = 60, dive = 240, window = 24, se_dive = 25),
    "both `se_surface` and `se_dive`"
  )
})

test_that("covariance between the interval means is carried", {
  indep <- availability(surface = 60, dive = 240, window = 24,
                        se_surface = 8, se_dive = 25)$se
  correlated <- availability(surface = 60, dive = 240, window = 24,
                             se_surface = 8, se_dive = 25,
                             cov_surface_dive = 150)$se
  # The two partial derivatives have opposite signs, so a positive covariance
  # reduces the variance of the result rather than inflating it.
  expect_false(isTRUE(all.equal(indep, correlated)))
})

test_that("impossible intervals are refused", {
  expect_error(availability(surface = 0, dive = 240, window = 24), "positive")
  expect_error(availability(surface = 60, dive = -1, window = 24), "positive")
  expect_error(availability(surface = 60, dive = 240, window = -1), "negative")
  expect_error(
    availability(surface = 60, dive = 240, window = 24,
                 se_surface = -1, se_dive = 2),
    "cannot be negative"
  )
})

test_that("a key that does not line up is refused", {
  expect_error(
    availability(surface = c(45, 60), dive = 240, window = 24,
                 key = c("Feb", "Mar", "Apr")),
    "does not line up"
  )
})

test_that("missing inputs give missing availability, not a wrong one", {
  a <- availability(surface = c(60, NA), dive = 240, window = 24)
  expect_false(is.na(a$value[1]))
  expect_true(is.na(a$value[2]))
})

test_that("the view window is a chord, and closes at the edge of view", {
  # On the trackline the chord is the full diameter.
  expect_equal(view_window(radius = 300, speed = 50), 2 * 300 / 50)

  # It narrows away from the trackline.
  w <- view_window(radius = 300, speed = 50, distance = c(0, 150, 290))
  expect_true(all(diff(w) < 0))

  # Beyond the edge of view there is no window at all.
  expect_equal(view_window(radius = 300, speed = 50, distance = 400), 0)
  expect_equal(view_window(radius = 300, speed = 50, distance = 300), 0)

  # Twice the speed, half the time in view.
  expect_equal(view_window(radius = 300, speed = 100),
               view_window(radius = 300, speed = 50) / 2)
})

test_that("view_window feeds availability, and distance lowers it", {
  r <- 300; v <- 50
  on_line <- availability(surface = 60, dive = 240,
                          window = view_window(r, v, 0))$value
  off_line <- availability(surface = 60, dive = 240,
                           window = view_window(r, v, 250))$value
  # A shorter window means less chance a diving animal surfaces in time.
  expect_lt(off_line, on_line)
})

test_that("impossible platform geometry is refused", {
  expect_error(view_window(radius = 0, speed = 50), "positive")
  expect_error(view_window(radius = 300, speed = 0), "point transect")
  expect_error(view_window(radius = 300, speed = 50, distance = -1), "negative")
})

test_that("an aircraft's window grows with distance, unlike a circular one", {
  # The whole point of the aerial geometry: a wedge running forward and aft is
  # wider further out, so an animal off the trackline stays in it longer. The
  # circular geometry does the opposite, and using it for an aircraft gets the
  # sign of the distance effect backwards.
  x <- c(0, 500, 1500, 3000)
  aerial <- view_window_aerial(trackline = 51.22, speed = 51.4, slope = 0.03,
                               distance = x)
  circular <- view_window(radius = 3000, speed = 51.4, distance = x)

  expect_true(all(diff(aerial) > 0))
  expect_true(all(diff(circular) < 0))

  # On the trackline it is the trackline time, by construction.
  expect_equal(aerial[1], 51.22)
})

test_that("the aerial window reproduces Ganley et al.'s measured platform", {
  # Their Skymaster: 51.22 s in view at the effective trackline, rising at
  # about 0.03 s per metre, giving roughly 130-150 s out at 3 km (their Fig. 1).
  w <- view_window_aerial(trackline = 51.22, speed = 51.4, slope = 0.03,
                          distance = 3000)
  expect_gt(w, 120)
  expect_lt(w, 160)
})

test_that("a half-angle and a calibrated slope are two routes to one window", {
  speed <- 51.4
  angle <- 57
  by_angle <- view_window_aerial(trackline = 51.22, speed = speed,
                                 angle = angle, distance = 1000)
  by_slope <- view_window_aerial(trackline = 51.22, speed = speed,
                                 slope = tan(angle * pi / 180) / speed,
                                 distance = 1000)
  expect_equal(by_angle, by_slope)
})

test_that("the aerial window refuses what it cannot mean", {
  expect_error(view_window_aerial(trackline = 51, speed = 51.4),
               "exactly one of")
  expect_error(
    view_window_aerial(trackline = 51, speed = 51.4, angle = 45, slope = 0.03),
    "exactly one of"
  )
  expect_error(view_window_aerial(trackline = 51, speed = 51.4, angle = 90),
               "between 0 and 90")
  # A window that shrinks with distance is the circular case, not this one.
  expect_error(view_window_aerial(trackline = 51, speed = 51.4, slope = -0.01),
               "use `view_window\\(\\)`")
})

test_that("Ganley's percent surface time is NOT the ratio of their means", {
  # Worth a test because it is the obvious wrong assumption, and the one this
  # package made until Table S1 arrived. The reported percentage is a mean of
  # per-follow percentages; the interval columns are means of intervals. A mean
  # of ratios is not a ratio of means, and here they differ by up to 35 points.
  from_intervals <- with(ganley_availability,
                         mean_surface / (mean_surface + mean_dive) * 100)
  reported <- ganley_availability$percent_surface_time

  expect_false(isTRUE(all.equal(from_intervals, reported)))
  expect_gt(max(abs(from_intervals - reported)), 30)

  # January's listed 16% against 8.3% from its intervals, April's 55% against
  # 90% - the gap runs in both directions, so it is not a constant offset.
  expect_lt(from_intervals[1], reported[1])
  expect_gt(from_intervals[4], reported[4])
})

test_that("Ganley's availability is not recoverable at one common window", {
  # Their a(x) is a bootstrap median and varies with distance through the time
  # in view. Solving each month for the window it implies shows they were not
  # evaluated at a shared one - January needs about two minutes, April a few
  # seconds. This is why the values ship as measurements to use rather than
  # outputs to recompute.
  implied <- vapply(seq_len(nrow(ganley_availability)), function(i) {
    s <- ganley_availability$mean_surface[i]
    d <- ganley_availability$mean_dive[i]
    target <- ganley_availability$availability[i]
    f <- function(w) availability(s, d, w)$value - target
    if (f(0) > 0) return(NA_real_)          # already above at an instant
    stats::uniroot(f, c(0, 1e5))$root
  }, numeric(1))

  expect_gt(implied[1], 100)                 # January, about two minutes
  expect_lt(implied[4], 20)                  # April, seconds
  expect_gt(max(implied, na.rm = TRUE) / min(implied, na.rm = TRUE), 5)
})

test_that("the shipped tables are the ones the supplement prints", {
  expect_equal(ganley_availability$percent_surface_time, c(16, 34, 31, 55))
  expect_equal(ganley_availability$mean_dive, c(533, 256, 219, 88))
  expect_equal(ganley_availability$mean_surface, c(48, 226, 67, 801))
  expect_equal(ganley_availability$availability, c(0.27, 0.52, 0.52, 0.91))
  expect_equal(ganley_availability$n_follows, c(7L, 9L, 22L, 48L))

  # The pooled row is an attribute, not a fifth month, so it cannot be summed
  # or plotted alongside them by accident.
  expect_equal(nrow(ganley_availability), 4L)
  expect_equal(attr(ganley_availability, "all_sightings")$n_follows, 86L)

  expect_equal(nrow(ganley_detection), 20L)
  expect_equal(range(ganley_detection$year), c(1998, 2017))
  expect_equal(range(ganley_detection$p), c(0.431, 0.866))
  # 2003 is the year that fails Kolmogorov-Smirnov and has the worst precision.
  worst <- ganley_detection[which.max(ganley_detection$p_se), ]
  expect_equal(worst$year, 2003)
  expect_lt(worst$ks_p, 0.05)
})

test_that("the measured season is threefold, which is the whole argument", {
  a <- ganley_availability$availability
  expect_equal(min(a), 0.27)
  expect_equal(max(a), 0.91)
  expect_gt(max(a) / min(a), 3)
})
