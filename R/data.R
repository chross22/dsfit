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
