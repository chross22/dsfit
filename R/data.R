# Documentation for the datasets shipped in data/.
#
#   example_dive_intervals  invented monthly dive intervals, for examples only
#   ganley_availability     real monthly availability, Ganley et al. (2019) Table S1
#   ganley_detection        real annual detection probability, Table S2


#' Example surfacing and diving intervals, by month (invented)
#'
#' Invented mean surfacing and diving intervals, with standard errors, for six
#' months. Use them to try [availability()]. They are not measurements and
#' must not be used to correct anything.
#'
#' @section What the numbers show:
#' The values are invented, and are not Ganley et al. (2019)'s Cape Cod Bay
#' values, which ship as [ganley_availability].
#'
#' They copy a real pattern. Ganley et al. (2019) found right whale
#' availability in Cape Cod Bay varying by month between 0.27 and 0.85 (the
#' range in their abstract), following the depth of the copepods the whales
#' fed on. Deep feeding means long dives and little surface time; shallow
#' feeding means the whales are available most of the time. With a 12 s window
#' (`view_window(radius = 300, speed = 50)`), [availability()] on these
#' intervals ranges from about 0.25 to 0.85. A single constant suited to one
#' month is wrong by a factor of three or more for another.
#'
#' Real dive intervals come from focal follows, tags or drones, and depend on
#' place, season and behaviour. Substitute your own.
#'
#' @format A tibble with 6 rows and 5 columns:
#' \describe{
#'   \item{month}{Month, as an ordered factor from December to May.}
#'   \item{surface}{Mean surfacing interval, seconds.}
#'   \item{dive}{Mean diving interval, seconds.}
#'   \item{se_surface}{Standard error of `surface`, seconds.}
#'   \item{se_dive}{Standard error of `dive`, seconds.}
#' }
#'
#' @references
#' Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what there
#' is: estimating North Atlantic right whale *Eubalaena glacialis* local
#' abundance. *Endangered Species Research* 38:101-113.
#' \doi{10.3354/esr00938} The monthly pattern these invented values copy.
#'
#' @seealso [availability()], [view_window()]
#'
#' @examples
#' example_dive_intervals
#'
#' # Availability for the season, from these intervals and a 12 s window
#' availability(
#'   surface    = example_dive_intervals$surface,
#'   dive       = example_dive_intervals$dive,
#'   window     = view_window(radius = 300, speed = 50),
#'   se_surface = example_dive_intervals$se_surface,
#'   se_dive    = example_dive_intervals$se_dive,
#'   key        = as.character(example_dive_intervals$month)
#' )
"example_dive_intervals"


#' Right whale availability in Cape Cod Bay, by month
#'
#' Table S1 of Ganley et al. (2019): mean dive and surface intervals, percent
#' time at the surface, and median availability, from 86 focal follows of North
#' Atlantic right whales in Cape Cod Bay in 2016 and 2017. These are real
#' measurements, unlike [example_dive_intervals].
#'
#' Availability runs from 0.27 in January to 0.91 in April, as the copepods
#' the whales feed on move up through the water column. This threefold change
#' within one season is measured evidence against a constant `g(0)`.
#'
#' @section Three things not to do with these columns:
#' \describe{
#'   \item{Do not treat `percent_surface_time` as
#'     \eqn{E(s)/(E(s)+E(d))}}{The two do not match. January lists 16%, but
#'     its mean intervals (48 s surface, 533 s dive) give
#'     \eqn{48/(48+533)} = 8.3%. The mismatch occurs every month and in both
#'     directions: April's intervals give 90% against a listed 55%. The
#'     percentage appears to be a mean of per-follow percentages, while the
#'     interval columns are means of intervals (a mean of ratios against a
#'     ratio of means). Using it as instantaneous availability would be wrong
#'     by up to 35 percentage points.}
#'   \item{Do not expect [availability()] to reproduce
#'     `availability`}{The reported value is a median over bootstrap
#'     replicates, and availability depends on the time in view, which varies
#'     with distance. Working back from the tabulated means, January's 0.27
#'     implies a window near 122 s and April's 0.91 a window near 8 s, so the
#'     months were not evaluated at a common window. Use the values as
#'     measurements; do not recompute them.}
#'   \item{Do not pass `availability_variance` to [g0()] without deciding what
#'     it is}{It is labelled a variance and ranges from 0.04 to 0.09. Read as
#'     a variance, January's standard error is \eqn{\sqrt{0.04} = 0.2} against
#'     an estimate of 0.27, a CV of 74%, which would dominate any correction.
#'     Read as a standard error, the CV is 15%. Neither matches the very small
#'     error bars in the paper's Fig. 4A. The column keeps the paper's label
#'     and is not converted.}
#' }
#'
#' @section The survey platform:
#' A Cessna 336/337 Skymaster at 185 km/h and 228 or 304 m altitude, surveying
#' January to May. Time in view at the effective trackline was 51.22 s, rising
#' with perpendicular distance (see [view_window_aerial()]). The effective
#' trackline is at 100 m, not 0 m, because the aircraft's flat windows leave a
#' blind spot beneath it, and the surveys were left-truncated at 100 m. At the
#' fitting end the same blind spot is handled by `sweep_models(left = )` or the
#' gamma key.
#'
#' @section Two inconsistencies in the source:
#' These are in the paper, not transcription errors.
#' \itemize{
#'   \item February's sample size is 9 in Table S1 (summing to the table's
#'     total of 86) but 10 in the Fig. 2 caption (summing to the 87 stated in
#'     the Methods).
#'   \item April's availability is 0.91 in Table S1 and Section 3.1, but the
#'     abstract gives the seasonal range as 0.27-0.85.
#' }
#' Table S1 is followed here, as the tabulated source.
#'
#' @format A tibble with 4 rows and 8 columns:
#' \describe{
#'   \item{month}{Month, as an ordered factor from January to April.}
#'   \item{percent_surface_time}{Percent of time at the surface. Not the ratio
#'     of the two interval columns; see above.}
#'   \item{mean_dive}{Mean diving interval, seconds.}
#'   \item{mean_surface}{Mean surfacing interval, seconds.}
#'   \item{availability}{Median availability, \eqn{a(x)}.}
#'   \item{availability_variance}{As labelled in the source; see above.}
#'   \item{n_follows}{Number of focal follows. January has 7 and April 48:
#'     the month with the lowest availability is also the least measured,
#'     because the weather that keeps the food deep also keeps observers
#'     on the ground.}
#'   \item{hours_followed}{Total follow time, hours.}
#' }
#'
#' The pooled "All sightings" row of Table S1 is stored as the
#' `all_sightings` attribute, not as a fifth row, so it is not summed or
#' plotted with the months by mistake.
#'
#' @source Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not
#'   what there is: estimating North Atlantic right whale *Eubalaena glacialis*
#'   local abundance. *Endangered Species Research* 38:101-113.
#'   \doi{10.3354/esr00938}, Table S1. Open access under CC-BY.
#'
#' @seealso [availability()], [g0()], [view_window_aerial()],
#'   [ganley_detection]
#'
#' @examples
#' ganley_availability
#'
#' # The monthly values as a g(0) component. A standard error must be chosen
#' # first (see the variance column above). This example reads the column as a
#' # variance; decide which reading is right before relying on it.
#' avail <- data.frame(
#'   key = as.character(ganley_availability$month),
#'   component = "availability",
#'   value = ganley_availability$availability,
#'   se = sqrt(ganley_availability$availability_variance)
#' )
#' suppressWarnings(g0(avail))
#'
#' # The pooled row
#' attr(ganley_availability, "all_sightings")$availability
"ganley_availability"


#' Right whale annual detection probability in Cape Cod Bay
#'
#' Table S2 of Ganley et al. (2019): goodness-of-fit p-values and average
#' detection probability for a detection function fitted separately to each
#' year from 1998 to 2017, from aerial line-transect surveys of Cape Cod Bay.
#'
#' @section This is `p`, not perception:
#' `p` is the average detection probability of the fitted detection function,
#' the same quantity [selection_table()] reports as `p` (with standard error
#' `p_se`). It is not a `g(0)` component. Ganley et al. state that perception
#' bias was not addressed in their study, because estimating it needs a
#' second observer team. Correcting for one while believing you have corrected
#' for the other biases abundance by a factor. [g0()] names any standard
#' component that is missing.
#'
#' @section What it is useful for:
#' It is a long, real example that detection probability is not a constant of
#' a survey. Twenty years of the same programme, aircraft and bay give `p`
#' between 0.431 and 0.866, with standard errors from 0.028 to 0.333.
#'
#' The goodness-of-fit columns show that a model can rank well and still fail
#' a fit test. Ganley et al. took p > 0.05 as adequate fit. 2003 fails
#' Kolmogorov-Smirnov (0.048) but passes Cramér-von Mises (0.103), and has the
#' largest standard error on `p`. This is why [selection_table()] reports
#' `cvm_p` alongside AIC.
#'
#' @format A tibble with 20 rows and 5 columns:
#' \describe{
#'   \item{year}{Survey year, 1998 to 2017.}
#'   \item{cvm_p}{Cramér-von Mises p-value for the year's detection function.}
#'   \item{ks_p}{Kolmogorov-Smirnov p-value.}
#'   \item{p}{Average detection probability, \eqn{\hat{P}_a}.}
#'   \item{p_se}{Standard error of `p`.}
#' }
#'
#' @source Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not
#'   what there is: estimating North Atlantic right whale *Eubalaena glacialis*
#'   local abundance. *Endangered Species Research* 38:101-113.
#'   \doi{10.3354/esr00938}, Table S2. Open access under CC-BY.
#'
#' @seealso [selection_table()], [ganley_availability]
#'
#' @examples
#' ganley_detection
#'
#' # Detection probability varies from year to year
#' range(ganley_detection$p)
#'
#' # The year that fails Kolmogorov-Smirnov also has the least precise p
#' ganley_detection[ganley_detection$ks_p < 0.05, ]
"ganley_detection"
