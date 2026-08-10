#' Availability: the probability an animal was at the surface to be seen
#'
#' Computes the probability that an animal was available to be detected while it
#' was within view, from its surfacing and diving intervals and the time the
#' survey platform kept it in view. This is the availability component of
#' `g(0)`, and it is **computed from external dive data, not estimated from the
#' survey**.
#'
#' @section Why this is a calculation and not an estimator:
#' An animal submerged while the aircraft passes is missed at every
#' perpendicular distance equally. To first order availability is a pure scale
#' factor on \eqn{g(x)}: it does not dent the near-zero end of the distance
#' distribution, does not change its shape, and leaves no signature for a
#' likelihood to find. That is exactly why a mis-specified `g(0)` shifts every
#' candidate in a [selection_table()] by the same factor and leaves the ranking
#' looking untouched.
#'
#' So there is nothing in the distances to estimate it from, and none of the
#' inputs here come from the survey being corrected. `surface` and `dive` come
#' from focal follows, tagging, or drone observation; `window` comes from the
#' platform's geometry and speed. Both arrive with their own uncertainty, which
#' is the point: it is visible rather than absorbed into a constant.
#'
#' @section The formula:
#' Following Laake and Borchers (2004), for a mean surfacing interval
#' \eqn{E(s)}, a mean diving interval \eqn{E(d)}, and a window \eqn{w} during
#' which the animal is in view:
#'
#' \deqn{a = \frac{E(s) + E(d)\left(1 - e^{-w/E(d)}\right)}{E(s) + E(d)}}
#'
#' The two limits are worth holding onto as a check on any number this returns.
#' As \eqn{w \to 0} the window is instantaneous and \eqn{a} becomes
#' \eqn{E(s)/(E(s) + E(d))}, the plain proportion of time spent at the surface.
#' As \eqn{w \to \infty} the platform watches forever, every animal surfaces
#' eventually, and \eqn{a \to 1}. The middle term is the extra chance that an
#' animal which was down when the window opened comes up before it closes.
#'
#' @section There is no single availability, and that is the finding:
#' Ganley et al. (2019) measured this for right whales in Cape Cod Bay and found
#' availability varying by month between **0.27 and 0.85**, tracking the depth of
#' the copepod layer the whales were feeding on, with perception varying
#' separately by year between 0.43 and 0.87.
#'
#' A plausible-looking single figure is therefore the most dangerous input this
#' package accepts. A value near 0.83 is a real number for some month, and wrong
#' by a factor of two for others — and because it scales every model equally, no
#' amount of model selection reveals it. Compute one per key and pass a vector,
#' rather than picking a representative value.
#'
#' @section How the standard error is obtained:
#' Availability is deterministic given `surface`, `dive` and `window`, so the
#' uncertainty comes from them. Given `se_surface` and `se_dive`, the Jacobian is
#' taken numerically and combined with their covariance matrix — the same
#' delta-method treatment [effect_estimates_ddf()] uses.
#'
#' Where the raw focal follows are in hand, resampling them is better than this:
#' Ganley et al. (2019) and others bootstrap the follows and take the standard
#' deviation of the resulting availabilities, which does not assume the interval
#' means are jointly normal. Do that and pass the result on as `se_surface` and
#' `se_dive`, or bypass this and build the rows directly.
#'
#' Without `se_surface` and `se_dive` the standard error is `NA`, deliberately:
#' this function will not invent a precision for a number whose CV routinely
#' dominates the CV of abundance.
#'
#' @param surface Mean surfacing interval, in seconds. A vector gives one result
#'   per element, which is how a per-month table is built.
#' @param dive Mean diving interval, in seconds. Recycled against `surface`.
#' @param window Time the animal is within view, in seconds. See [view_window()]
#'   to derive it from platform geometry. Recycled against `surface`.
#' @param se_surface,se_dive Standard errors of `surface` and `dive`. Both are
#'   needed for a standard error on the result; either alone is an error, since
#'   propagating one source of variance and silently dropping the other
#'   understates the total.
#' @param cov_surface_dive Covariance between the two interval means. Zero by
#'   default, which is right when they are estimated from separate follows and
#'   optimistic when they are not.
#' @param key Optional labels — months, platforms, years — carried through to the
#'   result so the rows stay attached to what they apply to.
#'
#' @return A tibble with one row per element of `surface`: `key`, `component`
#'   (always `"availability"`), `value`, and `se`. That is the shape a `g(0)`
#'   correction is assembled from, so availability and perception rows stack.
#'
#' @references
#' Laake, J.L. and Borchers, D.L. (2004) Methods for incomplete detection at
#' distance zero. In *Advanced Distance Sampling*, pp. 108-189. Oxford University
#' Press. The formula implemented here.
#'
#' Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what there
#' is: estimating North Atlantic right whale *Eubalaena glacialis* local
#' abundance. *Endangered Species Research* 38:101-113.
#' \doi{10.3354/esr00938} Focal follows and aircraft field of view applied to
#' right whales, and the monthly variation quoted above.
#'
#' Roberts, J.J., Yack, T.M., Fujioka, E., Halpin, P.N., Baumgartner, M.F. and
#' others (2024) North Atlantic right whale density surface model for the US
#' Atlantic evaluated with passive acoustic monitoring. *Marine Ecology Progress
#' Series* 732:167-192. \doi{10.3354/meps14547} Corrects perception and
#' availability per platform, team and conditions across 11 institutions, which
#' is the scale at which these corrections actually vary.
#'
#' @seealso [view_window()] for the window, [sweep_models()] for why none of this
#'   belongs in the detection function fit.
#'
#' @examples
#' # An instantaneous window is the proportion of time at the surface
#' availability(surface = 60, dive = 240, window = 0)
#'
#' # A real window lifts it: some animals that were down come up in time
#' availability(surface = 60, dive = 240, window = 30)
#'
#' # One row per month, which is how it is actually used
#' availability(
#'   surface = c(45, 60, 90),
#'   dive    = c(300, 240, 150),
#'   window  = 24,
#'   key     = c("Feb", "Mar", "Apr")
#' )
#'
#' # With uncertainty on the intervals
#' availability(surface = 60, dive = 240, window = 24,
#'              se_surface = 8, se_dive = 25)
#'
#' @export
availability <- function(surface, dive, window, se_surface = NULL,
                         se_dive = NULL, cov_surface_dive = 0, key = NULL) {
  stopifnot(is.numeric(surface), is.numeric(dive), is.numeric(window))

  n <- max(length(surface), length(dive), length(window))
  if (n == 0L) rlang::abort("`surface`, `dive` and `window` are all empty.")
  surface <- rep_len(surface, n)
  dive <- rep_len(dive, n)
  window <- rep_len(window, n)

  if (any(!is.na(surface) & surface <= 0) || any(!is.na(dive) & dive <= 0)) {
    rlang::abort(paste0(
      "`surface` and `dive` are interval durations in seconds and must be ",
      "positive. A zero surfacing interval is an animal that is never ",
      "available, which is not a correction but a survey that cannot work."
    ))
  }
  if (any(!is.na(window) & window < 0)) {
    rlang::abort("`window` is a duration in seconds and cannot be negative.")
  }

  # One of se_surface/se_dive without the other would propagate half the
  # variance and report it as though it were all of it.
  has_se <- !is.null(se_surface) || !is.null(se_dive)
  if (has_se && (is.null(se_surface) || is.null(se_dive))) {
    rlang::abort(paste0(
      "Give both `se_surface` and `se_dive`, or neither. Propagating one and ",
      "dropping the other reports a standard error smaller than the truth, ",
      "which is worse than reporting none."
    ))
  }
  if (has_se) {
    se_surface <- rep_len(se_surface, n)
    se_dive <- rep_len(se_dive, n)
    cov_surface_dive <- rep_len(cov_surface_dive, n)
    if (any(!is.na(se_surface) & se_surface < 0) ||
        any(!is.na(se_dive) & se_dive < 0)) {
      rlang::abort("Standard errors cannot be negative.")
    }
  }

  if (!is.null(key) && length(key) != n) {
    rlang::abort(paste0(
      "`key` has length ", length(key), " but there are ", n,
      " availabilities. A key that does not line up with its values is worse ",
      "than none."
    ))
  }

  value <- vapply(seq_len(n), function(i) {
    availability_one(c(surface[i], dive[i]), window[i])
  }, numeric(1))

  se <- rep(NA_real_, n)
  if (has_se) {
    se <- vapply(seq_len(n), function(i) {
      availability_se(
        c(surface[i], dive[i]), window[i],
        matrix(c(se_surface[i]^2, cov_surface_dive[i],
                 cov_surface_dive[i], se_dive[i]^2), nrow = 2L)
      )
    }, numeric(1))
  }

  tibble::tibble(
    key = if (is.null(key)) rep(NA_character_, n) else as.character(key),
    component = "availability",
    value = value,
    se = se
  )
}


#' The time an animal stays within view of the platform
#'
#' Turns platform geometry into the `window` [availability()] needs: how long a
#' point on the water remains inside the observers' field of view as the platform
#' passes.
#'
#' @section The geometry:
#' For a circular field of view of radius \eqn{r} and a platform travelling at
#' speed \eqn{v}, a point at perpendicular distance \eqn{x} is crossed along a
#' chord of length \eqn{2\sqrt{r^2 - x^2}}, so
#'
#' \deqn{w(x) = \frac{2\sqrt{r^2 - x^2}}{v}}
#'
#' The window therefore shrinks with perpendicular distance, and is zero beyond
#' the edge of the view — an animal outside it is never available, which this
#' returns as `0` rather than as an error.
#'
#' This is the simplest useful geometry and it will not fit every platform. A
#' forward-looking window, an obscured belly, or a viewing area that is not
#' circular all give a different \eqn{w}, and `availability()` takes `window`
#' directly so that a better one can be substituted. Measuring it, as Ganley et
#' al. (2019) did for the aircraft they flew, beats deriving it.
#'
#' @param radius Radius of the viewing area, in metres.
#' @param speed Platform ground speed, in metres per second.
#' @param distance Perpendicular distance, in metres. `0` (default) is on the
#'   trackline, which gives the widest window.
#'
#' @return A numeric vector of window durations in seconds, `0` where `distance`
#'   exceeds `radius`.
#'
#' @seealso [availability()]
#'
#' @examples
#' # A 300 m viewing radius at 50 m/s, on the trackline
#' view_window(radius = 300, speed = 50)
#'
#' # The window narrows away from the trackline, and closes at the edge
#' view_window(radius = 300, speed = 50, distance = c(0, 150, 290, 400))
#'
#' @export
view_window <- function(radius, speed, distance = 0) {
  stopifnot(is.numeric(radius), is.numeric(speed), is.numeric(distance))
  if (any(!is.na(radius) & radius <= 0)) {
    rlang::abort("`radius` must be positive.")
  }
  if (any(!is.na(speed) & speed <= 0)) {
    rlang::abort(paste0(
      "`speed` must be positive. A stationary platform keeps an animal in ",
      "view indefinitely, which is a point transect rather than a line one."
    ))
  }
  if (any(!is.na(distance) & distance < 0)) {
    rlang::abort("`distance` is a perpendicular distance and cannot be negative.")
  }

  half_chord <- sqrt(pmax(radius^2 - distance^2, 0))
  2 * half_chord / speed
}


# a = (s + d(1 - exp(-w/d))) / (s + d), from Laake and Borchers (2004).
# `par` is c(surface, dive) so that numDeriv can differentiate it as a vector.
availability_one <- function(par, window) {
  s <- par[[1]]
  d <- par[[2]]
  if (is.na(s) || is.na(d) || is.na(window)) return(NA_real_)
  (s + d * (1 - exp(-window / d))) / (s + d)
}

# Delta method: the Jacobian with respect to (surface, dive), combined with
# their covariance matrix. Numerical rather than hand-derived, matching how
# effect_estimates_ddf() handles the same problem.
availability_se <- function(par, window, vcov) {
  if (anyNA(par) || is.na(window) || anyNA(vcov)) return(NA_real_)

  g <- try(
    numDeriv::grad(availability_one, par, window = window),
    silent = TRUE
  )
  if (inherits(g, "try-error") || anyNA(g)) return(NA_real_)

  v <- drop(t(g) %*% vcov %*% g)
  if (is.na(v) || v < 0) return(NA_real_)
  sqrt(v)
}
