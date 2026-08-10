# Generates `ganley_availability` and `ganley_detection`.
#
# Transcribed from the supplement to Ganley, L.C., Brault, S. and Mayo, C.A.
# (2019) "What we see is not what there is: estimating North Atlantic right
# whale Eubalaena glacialis local abundance", Endangered Species Research
# 38:101-113, doi:10.3354/esr00938. Open access under CC-BY.
#
# Tables S1 and S2, verbatim. Three things about them are worth recording here
# rather than discovering later.
#
# 1. PERCENT SURFACE TIME IS NOT THE RATIO OF THE TABULATED MEANS. January is
#    listed at 16% surface time with a mean surface interval of 48 s and a mean
#    dive of 533 s; 48/(48+533) is 8.3%, not 16%. The same gap appears in every
#    month, in both directions (April: 801/(801+88) = 90%, listed as 55%). The
#    reported percentage is evidently a mean of per-follow percentages, while
#    the interval columns are means of intervals - a mean of ratios against a
#    ratio of means, which are different quantities. So `percent_surface_time`
#    is NOT E(s)/(E(s)+E(d)) and must not be used as the instantaneous
#    availability.
#
# 2. a(x) IS NOT REPRODUCIBLE FROM THE OTHER COLUMNS AT ANY SINGLE WINDOW. The
#    reported availability is a median over bootstrap replicates, and a(x)
#    depends on distance through the time in view. Working backwards from the
#    tabulated means, January's 0.27 implies a window near 122 s while April's
#    0.91 implies one near 8 s. They are not evaluated at a common window, so
#    availability() will not reproduce them from these columns - which is why
#    the values are shipped as measurements rather than recomputed.
#
# 3. THE VARIANCE COLUMN IS AMBIGUOUS. It is labelled "a(x) variance" and runs
#    0.04-0.09. Read literally, January's SE is sqrt(0.04) = 0.2 against an
#    estimate of 0.27 - a CV of 74%, which would dominate any correction it
#    entered. Read as a standard error instead, the CV is 15%. Neither reading
#    matches the very small error bars in the paper's Fig. 4A. It is shipped
#    under the paper's own label and NOT converted; anyone feeding it to `g0()`
#    has to decide what it is and be able to defend the choice.
#
# Two internal inconsistencies in the source, noted so they are not mistaken
# for transcription errors:
#   - February's sample size is 9 in Table S1 (summing to 86) but 10 in the
#     Fig. 2 caption (summing to 87, which is the figure the Methods gives).
#   - April's availability is 0.91 here and in Section 3.1, but the abstract
#     gives the seasonal range as 0.27-0.85.
# Table S1 is followed here, being the tabulated source.

# Table S1. Mean dive times, mean surface times, and median availability
# estimates a(x), from focal follows in Cape Cod Bay in 2016 and 2017.
ganley_availability <- tibble::tibble(
  month = factor(c("Jan", "Feb", "Mar", "Apr"),
                 levels = c("Jan", "Feb", "Mar", "Apr")),
  percent_surface_time = c(16, 34, 31, 55),
  mean_dive = c(533, 256, 219, 88),
  mean_surface = c(48, 226, 67, 801),
  availability = c(0.27, 0.52, 0.52, 0.91),
  availability_variance = c(0.04, 0.09, 0.06, 0.05),
  n_follows = c(7L, 9L, 22L, 48L),
  hours_followed = c(3.6, 4.0, 16.2, 37.9)
)

# The "All sightings" row of Table S1, kept as an attribute rather than a fifth
# row: it is a pooled summary, not another month, and would otherwise be summed
# or plotted alongside them by accident.
attr(ganley_availability, "all_sightings") <- list(
  percent_surface_time = 44, mean_dive = 175, mean_surface = 492,
  availability = 0.63, availability_variance = 0.07,
  n_follows = 86L, hours_followed = 61.8
)

# Table S2. Goodness-of-fit and detection probabilities for the annual
# detection functions, from the distance sampling analyses.
ganley_detection <- tibble::tibble(
  year = 1998:2017,
  cvm_p = c(0.751, 0.521, 0.531, 0.867, 0.970, 0.103, 0.969, 0.985, 0.431,
            0.866, 0.970, 0.581, 0.668, 0.847, 0.991, 0.980, 0.930, 0.583,
            0.115, 0.341),
  ks_p = c(0.564, 0.426, 0.248, 0.802, 0.965, 0.048, 0.965, 0.995, 0.392,
           0.777, 0.977, 0.507, 0.785, 0.737, 0.982, 0.993, 0.858, 0.815,
           0.065, 0.129),
  p = c(0.663, 0.601, 0.750, 0.640, 0.708, 0.760, 0.655, 0.866, 0.703, 0.712,
        0.595, 0.591, 0.599, 0.604, 0.611, 0.576, 0.507, 0.541, 0.431, 0.518),
  p_se = c(0.156, 0.094, 0.041, 0.070, 0.104, 0.333, 0.142, 0.059, 0.163,
           0.065, 0.050, 0.038, 0.046, 0.043, 0.062, 0.030, 0.035, 0.028,
           0.038, 0.028)
)

stopifnot(
  nrow(ganley_detection) == 20L,
  !anyNA(ganley_detection),
  all(vapply(ganley_detection[-1], function(x) all(x > 0 & x <= 1), logical(1)))
)

usethis::use_data(ganley_availability, overwrite = TRUE)
usethis::use_data(ganley_detection, overwrite = TRUE)
