#' Example surfacing and diving intervals, by month
#'
#' A worked example for [availability()]: mean surfacing and diving intervals
#' over a six-month season, with standard errors.
#'
#' @section These numbers are invented:
#' They are not measurements. In particular they are **not** Ganley et al.
#' (2019)'s Cape Cod Bay values, which are not open access. Do not use them to
#' correct anything.
#'
#' What is real is the pattern they were built to show. Ganley et al. (2019)
#' found right whale availability in Cape Cod Bay varying by month between
#' **0.27 and 0.85**, tracking the depth of the copepod layer the whales were
#' feeding on: feeding deep means long dives and little time at the surface,
#' feeding shallow means available most of the time. These intervals reproduce
#' that shape, and running [availability()] over them spans roughly 0.25 to
#' 0.85 — which is the point of shipping them. A single representative
#' availability is a real number for one month and wrong by a factor of two for
#' others.
#'
#' Real dive parameters come from focal follows, tagging, or drone work, and are
#' specific to a place, a season and a behaviour. Substitute your own.
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
#' \doi{10.3354/esr00938} The monthly pattern these values imitate, and the real
#' measurements they are not.
#'
#' @seealso [availability()], [view_window()]
#'
#' @examples
#' example_dive_intervals
#'
#' # The season's availability, from these intervals and a platform's geometry
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
#' surface time, and median availability, from 86 focal follows of North
#' Atlantic right whales in Cape Cod Bay during 2016 and 2017. Unlike
#' [example_dive_intervals], these are real measurements.
#'
#' Availability runs from **0.27 in January to 0.91 in April**, as the copepods
#' the whales feed on move up through the water column. That threefold swing
#' across one season is the case against a constant `g(0)`, in measurements.
#'
#' @section Three things not to do with these columns:
#' \describe{
#'   \item{Do not treat `percent_surface_time` as
#'     \eqn{E(s)/(E(s)+E(d))}}{It is not. January is listed at 16% with a mean
#'     surface interval of 48 s and a mean dive of 533 s, and
#'     \eqn{48/(48+533)} is 8.3%. The gap appears every month and in both
#'     directions — April's intervals give 90% against a listed 55%. The
#'     percentage is evidently a mean of per-follow percentages while the
#'     interval columns are means of intervals: a mean of ratios against a
#'     ratio of means. So it is **not** the instantaneous availability, and
#'     using it as one would be wrong by up to 35 percentage points.}
#'   \item{Do not expect [availability()] to reproduce
#'     `availability`}{The reported figure is a median over bootstrap
#'     replicates, and \eqn{a(x)} varies with distance through the time in
#'     view. Working backwards from the tabulated means, January's 0.27
#'     implies a window near 122 s and April's 0.91 one near 8 s — they are not
#'     evaluated at a common window. The values are measurements to be used,
#'     not outputs to be recomputed.}
#'   \item{Do not feed `availability_variance` to [g0()] without
#'     deciding what it is}{It is labelled a variance and runs 0.04–0.09. Read
#'     literally, January's standard error is \eqn{\sqrt{0.04} = 0.2} against
#'     an estimate of 0.27 — a CV of 74%, which would dominate any correction
#'     it entered. Read as a standard error instead, the CV is 15%. Neither
#'     matches the very small error bars in the paper's Fig. 4A. It is shipped
#'     under the paper's own label and deliberately not converted.}
#' }
#'
#' @section The platform:
#' Cessna 336/337 Skymaster at 185 km/h and 228 or 304 m altitude, surveying
#' January to May. Time in view at the effective trackline was **51.22 s**,
#' rising with perpendicular distance — see [view_window_aerial()]. Their
#' trackline is at 100 m rather than 0 m because the aircraft's flat windows
#' leave a blind spot beneath it, and the surveys were left-truncated there
#' accordingly — the same blind spot `sweep_models(left = )` and the gamma key
#' handle at the fitting end.
#'
#' @section Two inconsistencies in the source:
#' Recorded so they are not mistaken for transcription errors. February's
#' sample size is 9 in Table S1, summing to the 86 its own total row gives, but
#' 10 in the Fig. 2 caption, summing to the 87 the Methods states. And April's
#' availability is 0.91 here and in Section 3.1, while the abstract gives the
#' seasonal range as 0.27–0.85. Table S1 is followed here, being the tabulated
#' source.
#'
#' @format A tibble with 4 rows and 8 columns:
#' \describe{
#'   \item{month}{Month, as an ordered factor from January to April.}
#'   \item{percent_surface_time}{Percent of time at the surface. See above for
#'     why this is not the ratio of the two interval columns.}
#'   \item{mean_dive}{Mean diving interval, seconds.}
#'   \item{mean_surface}{Mean surfacing interval, seconds.}
#'   \item{availability}{Median availability, \eqn{a(x)}.}
#'   \item{availability_variance}{As labelled in the source; see above.}
#'   \item{n_follows}{Focal follows behind each row. January's 7 against
#'     April's 48 is worth noticing — the month with the lowest availability is
#'     also the one measured least well, because the weather that keeps the food
#'     deep also keeps observers on the ground.}
#'   \item{hours_followed}{Total follow time, hours.}
#' }
#'
#' The pooled "All sightings" row of Table S1 is kept as the `all_sightings`
#' attribute rather than a fifth row, so it cannot be summed or plotted
#' alongside the months by accident.
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
#' # The measured seasonal swing, as a g(0) component. A standard error has to
#' # be decided on first - see above on the variance column - so this uses a
#' # deliberately explicit placeholder rather than a silent conversion.
#' avail <- data.frame(
#'   key = as.character(ganley_availability$month),
#'   component = "availability",
#'   value = ganley_availability$availability,
#'   se = sqrt(ganley_availability$availability_variance)
#' )
#' suppressWarnings(g0(avail))
#'
#' attr(ganley_availability, "all_sightings")$availability
"ganley_availability"


#' Right whale annual detection probability in Cape Cod Bay
#'
#' Table S2 of Ganley et al. (2019): the goodness-of-fit and average detection
#' probability of a separately fitted detection function for each year from 1998
#' to 2017, from aerial line-transect surveys of Cape Cod Bay.
#'
#' @section This is `p`, not perception bias:
#' `p` here is the average detection probability of the fitted detection
#' function — the same quantity [selection_table()] reports as `p`, with its
#' standard error as `p_se`. It is **not** a `g(0)` component. Ganley et al.
#' say so directly: perception bias "was not addressed directly in this study",
#' because estimating it needs a second observer team they did not have.
#'
#' The distinction matters because the two are easy to conflate and correcting
#' for one while believing you have corrected for the other is how these
#' estimates go wrong by a factor. See [g0()], which will name any component
#' you have not supplied.
#'
#' @section What it is good for:
#' It is a long, real example of the thing [sweep_models()] is built around:
#' **detection probability is not a constant of the survey**. Twenty years of
#' the same programme, the same aircraft and the same bay give `p` between
#' 0.431 and 0.866 — a twofold range — with standard errors spanning an order
#' of magnitude, from 0.028 to 0.333.
#'
#' The goodness-of-fit columns are worth reading next to it. Ganley et al. took
#' p > 0.05 as adequate fit, and 2003 fails on Kolmogorov-Smirnov (0.048) while
#' passing Cramér-von Mises (0.103); it also carries much the largest standard
#' error on `p`. A model can top a ranking and still fail a fit test, which is
#' why [selection_table()] carries `cvm_p` alongside the AIC.
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
#' # Detection probability is not a constant of a survey programme
#' range(ganley_detection$p)
#'
#' # The year that fails Kolmogorov-Smirnov also has the worst precision on p
#' ganley_detection[ganley_detection$ks_p < 0.05, ]
"ganley_detection"
