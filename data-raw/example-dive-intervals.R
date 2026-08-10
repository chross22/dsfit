# Generates `example_dive_intervals`.
#
# THESE NUMBERS ARE INVENTED. They are not measurements, and in particular they
# are not Ganley et al. (2019)'s Cape Cod Bay values, which are behind a
# paywall and were never available to this script.
#
# What is real is the *pattern* they were chosen to show. Ganley et al. (2019)
# found right whale availability in Cape Cod Bay varying by month between 0.27
# and 0.85, tracking the depth of the copepod layer: whales feeding deep make
# long dives and are rarely at the surface, whales feeding near the surface are
# available most of the time. These intervals reproduce that shape - dive time
# falling and surface time rising through the season - so that a worked example
# shows why availability cannot be a constant.
#
# To use real numbers, replace this file with a reader for the real table and
# say in the documentation where it came from.

example_dive_intervals <- tibble::tibble(
  month = factor(
    c("Dec", "Jan", "Feb", "Mar", "Apr", "May"),
    levels = c("Dec", "Jan", "Feb", "Mar", "Apr", "May")
  ),
  # Mean surfacing and diving intervals, seconds.
  surface = c(35, 40, 50, 65, 95, 130),
  dive = c(155, 140, 120, 95, 60, 35),
  # Standard errors, loosely scaled to the means so the example can show
  # variance propagation doing something.
  se_surface = c(6, 7, 8, 10, 14, 19),
  se_dive = c(24, 21, 18, 14, 9, 6)
)

usethis::use_data(example_dive_intervals, overwrite = TRUE)
