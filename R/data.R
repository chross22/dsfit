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


#' Right whale percent surface time in Cape Cod Bay, by month
#'
#' Measured percent surface time from 87 focal follows of North Atlantic right
#' whales in Cape Cod Bay, from Ganley et al. (2019). Unlike
#' [example_dive_intervals], these are real measurements.
#'
#' @section What percent surface time is:
#' It is \eqn{E(s) / (E(s) + E(d))} — which is exactly the \eqn{w \to 0} limit
#' of [availability()], the availability of an animal glimpsed instantaneously.
#' Any real window can only raise it, since an animal that was down when the
#' window opened may surface before it closes.
#'
#' So this column is a **floor** on availability, and the seasonal shape of the
#' problem in its rawest form: 16% of the time at the surface in January against
#' 55% in April, as the copepods the whales feed on move up through the water
#' column.
#'
#' @section Why availability itself is not here:
#' Ganley et al. report January's availability exactly (0.27) and give the range
#' in their abstract (0.27–0.85), but the remaining monthly values appear only
#' in their Fig. 4A and in Table S1 of the supplement. Reading bar heights off a
#' figure is not a measurement, so this dataset stops where the text does. See
#' `data-raw/ganley-surface-time.R` for how to add them from the supplement.
#'
#' For the same reason there are no dive times here. The paper gives their
#' ranges across months — mean surface times 0.8 to 13.35 min, mean dive times
#' 1.30 to 8.83 min — but the per-month values are in Table S1.
#'
#' @section The platform:
#' Cessna 336/337 Skymaster at 185 km/h and 228 or 304 m altitude, surveying
#' January to May, 1998–2017. Time in view at the effective trackline was
#' **51.22 s**, rising with perpendicular distance — see [view_window_aerial()].
#' Their trackline is at 100 m rather than 0 m because the aircraft's flat
#' windows leave a blind spot beneath it, and the surveys were left-truncated
#' there accordingly.
#'
#' @format A tibble with 4 rows and 3 columns:
#' \describe{
#'   \item{month}{Month, as an ordered factor from January to April.}
#'   \item{percent_surface_time}{Percent of time at the surface, from focal
#'     follows. Divide by 100 for the instantaneous availability.}
#'   \item{n_follows}{Number of focal follows behind each figure. January's 7
#'     is worth noticing next to April's 48 — the month with the lowest
#'     availability is also the month it is least well measured, because the
#'     weather that keeps whales' food deep also keeps observers on the ground.}
#' }
#'
#' @source Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not
#'   what there is: estimating North Atlantic right whale *Eubalaena glacialis*
#'   local abundance. *Endangered Species Research* 38:101-113.
#'   \doi{10.3354/esr00938}. Section 3.2 and the Fig. 2 caption. Open access
#'   under CC-BY.
#'
#' @seealso [availability()], [view_window_aerial()]
#'
#' @examples
#' ganley_surface_time
#'
#' # Percent surface time is availability at an instantaneous window
#' pst <- ganley_surface_time$percent_surface_time / 100
#'
#' # A real window can only raise it. Ganley et al. measured 51.22 s at the
#' # trackline; dive times are theirs only in range, so this is illustrative.
#' availability(
#'   surface = pst,
#'   dive    = 1 - pst,
#'   window  = 51.22 / (6 * 60),   # window as a fraction of a 6 min cycle
#'   key     = as.character(ganley_surface_time$month)
#' )
"ganley_surface_time"
