# Availability and view windows.
#
# Availability is the probability that an animal is at the surface, and so can
# be seen, while it is within the observers' view. It is one component of g(0),
# the probability of detecting an animal on the trackline. It is computed here
# from dive data collected outside the survey, never estimated from the survey
# distances. g0() in R/g0.R combines it with the other components.
#
# Contents, in file order:
#   availability()        availability from surfacing and diving intervals
#   view_window()         time in view for a circular viewing area
#   view_window_aerial()  time in view for an aircraft side window
#   view_window_vessel()  time in view for a ship's forward search arc
#   check_recycle(), availability_one(), availability_se()   internal helpers


#' Availability of an animal at the surface
#'
#' Computes availability: the probability that an animal is at the surface,
#' and so can be seen, while it is within the observers' view. Use it to build
#' the availability component of a `g(0)` correction, one row per key (month,
#' platform, year), and pass the result to [g0()].
#'
#' @section Why availability is computed, not estimated:
#' An animal that is submerged while the platform passes is missed at every
#' perpendicular distance equally. Availability therefore scales the detection
#' function \eqn{g(x)} without changing its shape, and the distances carry no
#' information about it. The inputs come from outside the survey: `surface`
#' and `dive` from focal follows, tags or drones, and `window` from the
#' platform's geometry and speed.
#'
#' @section How it works:
#' Following Laake et al. (1997), for a mean surfacing interval \eqn{E(s)}, a
#' mean diving interval \eqn{E(d)}, and a window \eqn{w} (the time the animal
#' is in view):
#'
#' \deqn{a = \frac{E(s) + E(d)\left(1 - e^{-w/E(d)}\right)}{E(s) + E(d)}}
#'
#' At \eqn{w = 0}, \eqn{a} is \eqn{E(s)/(E(s) + E(d))}, the proportion of time
#' at the surface. As \eqn{w} grows, \eqn{a} approaches 1. The middle term is
#' the chance that an animal that was submerged when the window opened
#' surfaces before it closes.
#'
#' 1. The inputs are checked and recycled to a common length. Each must have
#'    length 1 or the common length.
#' 2. `surface` and `dive` must be positive and `window` must not be negative.
#' 3. `se_surface` and `se_dive` must be given together or not at all.
#' 4. \eqn{a} is computed for each element.
#' 5. If standard errors were given, the standard error of \eqn{a} is computed
#'    by the delta method: the gradient of \eqn{a} with respect to
#'    \eqn{(E(s), E(d))} is taken numerically and combined with their
#'    covariance matrix, as in [effect_estimates_ddf()].
#'
#' Without `se_surface` and `se_dive` the standard error is `NA`. The function
#' does not invent a precision, because the CV of availability often dominates
#' the CV of abundance. If the raw focal follows are available, bootstrapping
#' them (as Ganley et al. 2019 did) avoids assuming the interval means are
#' jointly normal; the result can be passed as `se_surface` and `se_dive`, or
#' the component rows can be built directly.
#'
#' @section Use one value per key, not one constant:
#' Availability varies. Ganley et al. (2019) found right whale availability in
#' Cape Cod Bay ranging from 0.27 in January to 0.91 in April, following the
#' depth of the copepods the whales fed on (see [ganley_availability]). A
#' single constant suited to one month can be wrong by a factor of three for
#' another. Because it scales every candidate model equally, model selection
#' cannot detect the error. Compute one value per key and pass a vector.
#'
#' This function gives availability only. Perception (whether an animal at the
#' surface and in view is seen) is a separate component and needs a second
#' observer team to estimate. [g0()] warns if it is missing.
#'
#' @param surface Numeric. Mean surfacing interval, in seconds. Must be
#'   positive. A vector gives one result per element.
#' @param dive Numeric. Mean diving interval, in seconds. Must be positive.
#'   Length 1 or the length of the longest input.
#' @param window Numeric. Time the animal is in view, in seconds. Must not be
#'   negative. Derive it with [view_window_aerial()], [view_window_vessel()] or
#'   [view_window()]. Length 1 or the length of the longest input.
#' @param se_surface,se_dive Numeric or `NULL` (default). Standard errors of
#'   `surface` and `dive`, in seconds. Give both or neither: propagating one
#'   and dropping the other would understate the total variance.
#' @param cov_surface_dive Numeric. Covariance between the two interval means,
#'   in seconds squared. Default `0`, which is right when they come from
#'   separate follows and optimistic otherwise. Used only with standard errors.
#' @param key Optional vector of labels (months, platforms, years), one per
#'   result. Must have the same length as the result. Stored as character.
#'
#' @return A tibble with one row per element of the longest input, and columns:
#'   \describe{
#'     \item{key}{The label from `key`, or `NA`.}
#'     \item{component}{Always `"availability"`.}
#'     \item{value}{Availability, between 0 and 1.}
#'     \item{se}{Standard error of `value`, or `NA` without `se_surface` and
#'       `se_dive`.}
#'   }
#'   This is the component shape [g0()] accepts, so availability and perception
#'   rows can be stacked.
#'
#' @references
#' Laake, J.L., Calambokidis, J., Osmek, S.D. and Rugh, D.J. (1997) Probability
#' of detecting harbor porpoise from aerial surveys: estimating g(0). *The
#' Journal of Wildlife Management* 61:63-75. \doi{10.2307/3802415} Source of
#' the formula.
#'
#' Laake, J.L. and Borchers, D.L. (2004) Methods for incomplete detection at
#' distance zero. In *Advanced Distance Sampling*, pp. 108-189. Oxford University
#' Press. Why perception and availability need separate data.
#'
#' Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what there
#' is: estimating North Atlantic right whale *Eubalaena glacialis* local
#' abundance. *Endangered Species Research* 38:101-113.
#' \doi{10.3354/esr00938} Focal follows and monthly availability for right
#' whales.
#'
#' Roberts, J.J., Yack, T.M., Fujioka, E., Halpin, P.N., Baumgartner, M.F. and
#' others (2024) North Atlantic right whale density surface model for the US
#' Atlantic evaluated with passive acoustic monitoring. *Marine Ecology Progress
#' Series* 732:167-192. \doi{10.3354/meps14547} Corrections for perception and
#' availability by platform, team and conditions.
#'
#' @seealso [view_window_aerial()], [view_window_vessel()] and [view_window()]
#'   for the window; [g0()] for combining components; [sweep_models()] for the
#'   detection function fit, which this does not enter.
#'
#' @examples
#' # With a zero window, availability is the proportion of time at the surface
#' availability(surface = 60, dive = 240, window = 0)
#'
#' # A longer window raises it: some submerged animals surface in time
#' availability(surface = 60, dive = 240, window = 30)
#'
#' # One row per month
#' availability(
#'   surface = c(45, 60, 90),
#'   dive    = c(300, 240, 150),
#'   window  = 24,
#'   key     = c("Feb", "Mar", "Apr")
#' )
#'
#' # With standard errors on the intervals
#' availability(surface = 60, dive = 240, window = 24,
#'              se_surface = 8, se_dive = 25)
#'
#' @export
availability <- function(surface, dive, window, se_surface = NULL,
                         se_dive = NULL, cov_surface_dive = 0, key = NULL) {
  stopifnot(is.numeric(surface), is.numeric(dive), is.numeric(window))

  # Recycle the three inputs to a common length n, allowing only length 1 or n.
  n <- max(length(surface), length(dive), length(window))
  if (n == 0L) rlang::abort("`surface`, `dive` and `window` are all empty.")
  check_recycle(list(surface = surface, dive = dive, window = window), n)
  surface <- rep_len(surface, n)
  dive <- rep_len(dive, n)
  window <- rep_len(window, n)

  # Check ranges. NA is allowed and passes through to an NA result.
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

  # Require both standard errors or neither. One alone would propagate half
  # the variance and report it as the whole.
  has_se <- !is.null(se_surface) || !is.null(se_dive)
  if (has_se && (is.null(se_surface) || is.null(se_dive))) {
    rlang::abort(paste0(
      "Give both `se_surface` and `se_dive`, or neither. Propagating one and ",
      "dropping the other reports a standard error smaller than the truth, ",
      "which is worse than reporting none."
    ))
  }
  # Recycle and check the standard errors and covariance.
  if (has_se) {
    check_recycle(list(se_surface = se_surface, se_dive = se_dive,
                       cov_surface_dive = cov_surface_dive), n)
    se_surface <- rep_len(se_surface, n)
    se_dive <- rep_len(se_dive, n)
    cov_surface_dive <- rep_len(cov_surface_dive, n)
    if (any(!is.na(se_surface) & se_surface < 0) ||
        any(!is.na(se_dive) & se_dive < 0)) {
      rlang::abort("Standard errors cannot be negative.")
    }
  }

  # `key` is not recycled: it must label every row.
  if (!is.null(key) && length(key) != n) {
    rlang::abort(paste0(
      "`key` has length ", length(key), " but there are ", n,
      " availabilities. A key that does not line up with its values is worse ",
      "than none."
    ))
  }

  # Availability for each element.
  value <- vapply(seq_len(n), function(i) {
    availability_one(c(surface[i], dive[i]), window[i])
  }, numeric(1))

  # Delta-method standard error for each element, using a 2 x 2 covariance
  # matrix of (surface, dive). NA when no standard errors were given.
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

  # Return in the component shape g0() accepts.
  tibble::tibble(
    key = if (is.null(key)) rep(NA_character_, n) else as.character(key),
    component = "availability",
    value = value,
    se = se
  )
}


#' Time in view for a circular viewing area
#'
#' Computes the `window` for [availability()] when the observers watch a
#' circular area around the platform, fore and aft alike. The window is the
#' time a point on the water spends inside that circle as the platform passes.
#' For an aircraft side window use [view_window_aerial()]; for a ship's
#' forward search use [view_window_vessel()].
#'
#' @section How it works:
#' For a viewing radius \eqn{r} and platform speed \eqn{v}, a point at
#' perpendicular distance \eqn{x} is crossed along a chord of length
#' \eqn{2\sqrt{r^2 - x^2}}, so
#'
#' \deqn{w(x) = \frac{2\sqrt{r^2 - x^2}}{v}}
#'
#' The window is widest on the trackline, shrinks with distance, and is `0`
#' beyond the radius.
#'
#' @section Which platform this fits:
#' Only a platform that watches a circular area. An aircraft observer looking
#' through a side window sees a wedge running forward and aft, so time in view
#' grows with perpendicular distance instead of shrinking; this function gets
#' that direction wrong. Ganley et al. (2019) measured time in view rising
#' with distance for their aircraft, at about 0.03 s per metre. A window
#' measured in the field, as they did by timing a buoy through the field of
#' view, is better than any of these formulas.
#'
#' @param radius Numeric. Radius of the viewing area, in metres. Must be
#'   positive.
#' @param speed Numeric. Platform ground speed, in metres per second. Must be
#'   positive.
#' @param distance Numeric. Perpendicular distance, in metres. Default `0`, the
#'   trackline, which gives the widest window. Must not be negative.
#'
#' @return A numeric vector of windows in seconds, `0` where `distance` exceeds
#'   `radius`.
#'
#' @references
#' Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what there
#' is: estimating North Atlantic right whale *Eubalaena glacialis* local
#' abundance. *Endangered Species Research* 38:101-113.
#' \doi{10.3354/esr00938} Time in view measured for an aircraft, rising with
#' perpendicular distance.
#'
#' @seealso [view_window_aerial()] for an aircraft, [view_window_vessel()] for
#'   a ship, and [availability()] for what the window feeds.
#'
#' @examples
#' # A 300 m viewing radius at 50 m/s, on the trackline
#' view_window(radius = 300, speed = 50)
#'
#' # The window narrows away from the trackline and is 0 beyond the radius
#' view_window(radius = 300, speed = 50, distance = c(0, 150, 290, 400))
#'
#' @export
view_window <- function(radius, speed, distance = 0) {
  # Check inputs.
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

  # Chord length divided by speed. Outside the radius the chord is zero.
  half_chord <- sqrt(pmax(radius^2 - distance^2, 0))
  2 * half_chord / speed
}


#' Time in view for an aircraft
#'
#' Computes the `window` for [availability()] for an observer looking through
#' the side window of an aircraft. The field of view is a wedge running forward
#' and aft of perpendicular, so time in view grows with perpendicular distance.
#' This is the usual geometry for an aerial line-transect survey.
#'
#' @section How it works:
#' Following the linear form fitted by Ganley et al. (2019), for a time in view
#' \eqn{\alpha} at the trackline, a viewing half-angle \eqn{\theta} forward of
#' perpendicular and the same aft, and ground speed \eqn{v}:
#'
#' \deqn{t(x) = \alpha + \frac{2 x \tan\theta}{v}}
#'
#' A point at perpendicular distance \eqn{x} enters the wedge \eqn{x \tan\theta}
#' ahead of abeam and leaves \eqn{x \tan\theta} behind, so it is in view for an
#' extra \eqn{2 x \tan\theta} of track.
#'
#' 1. Exactly one of `angle` and `slope` must be given.
#' 2. If `angle` is given, the slope is computed as
#'    \eqn{2 \tan\theta / v} seconds per metre.
#' 3. The window is `trackline + distance * slope`.
#'
#' @section Calibrating the inputs:
#' Measured values are better than derived ones. Ganley et al. (2019) flew past
#' a navigation buoy and timed it through the field of view. For a Cessna
#' Skymaster at 185 km/h and 228 or 304 m altitude, they found 51.22 s in view
#' at the effective trackline and a slope of about 0.03 s per metre. The two
#' altitudes did not differ detectably.
#'
#' Their effective trackline is at 100 m, not 0 m, because the aircraft's flat
#' windows leave a blind spot beneath it; their data were left-truncated at
#' 100 m for that reason. At the fitting end the same blind spot is handled by
#' `sweep_models(left = )` or the gamma key.
#'
#' @param trackline Numeric. Time in view at the trackline, \eqn{\alpha}, in
#'   seconds. Must not be negative.
#' @param speed Numeric. Ground speed, in metres per second. Must be positive.
#' @param angle Numeric or `NULL` (default). Viewing half-angle \eqn{\theta}, in
#'   degrees from perpendicular: the view runs `angle` forward and `angle` aft.
#'   Must be above 0 and below 90. Give this or `slope`, not both.
#' @param slope Numeric or `NULL` (default). Seconds in view gained per metre
#'   of perpendicular distance, \eqn{2 \tan\theta / v}, if calibrated directly.
#'   Must not be negative. Give this or `angle`, not both.
#' @param distance Numeric. Perpendicular distance, in metres. Default `0`.
#'   Must not be negative.
#'
#' @return A numeric vector of windows in seconds.
#'
#' @references
#' Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what there
#' is: estimating North Atlantic right whale *Eubalaena glacialis* local
#' abundance. *Endangered Species Research* 38:101-113.
#' \doi{10.3354/esr00938}
#'
#' @seealso [availability()], [view_window()], [view_window_vessel()],
#'   [ganley_availability]
#'
#' @examples
#' # Ganley et al.'s Skymaster, calibrated directly: 51.22 s at the trackline,
#' # plus about 0.03 s per metre
#' view_window_aerial(trackline = 51.22, speed = 51.4, slope = 0.03,
#'                    distance = c(0, 1000, 3000))
#'
#' # The same from a half-angle: 2 * tan(37.6 deg) / 51.4 is about 0.03 s/m
#' view_window_aerial(trackline = 51.22, speed = 51.4, angle = 37.6,
#'                    distance = c(0, 1000, 3000))
#'
#' @export
view_window_aerial <- function(trackline, speed, angle = NULL, slope = NULL,
                               distance = 0) {
  # Check inputs: exactly one of angle and slope, and valid ranges.
  stopifnot(is.numeric(trackline), is.numeric(speed), is.numeric(distance))
  if (is.null(angle) == is.null(slope)) {
    rlang::abort(paste0(
      "Give exactly one of `angle` (the viewing half-angle in degrees) or ",
      "`slope` (seconds of view per metre, if it was calibrated directly)."
    ))
  }
  if (any(!is.na(trackline) & trackline < 0)) {
    rlang::abort("`trackline` is a duration in seconds and cannot be negative.")
  }
  if (any(!is.na(speed) & speed <= 0)) {
    rlang::abort("`speed` must be positive.")
  }
  if (any(!is.na(distance) & distance < 0)) {
    rlang::abort("`distance` is a perpendicular distance and cannot be negative.")
  }

  # Get the slope in seconds per metre, from the angle or as given.
  if (is.null(slope)) {
    if (any(!is.na(angle) & (angle <= 0 | angle >= 90))) {
      rlang::abort(paste0(
        "`angle` is a half-angle in degrees from perpendicular, so it must be ",
        "between 0 and 90. At 90 the field of view never closes."
      ))
    }
    # The wedge runs forward and aft, so a point at distance x is in view for
    # 2 * x * tan(angle) of track.
    slope <- 2 * tan(angle * pi / 180) / speed
  } else if (any(!is.na(slope) & slope < 0)) {
    rlang::abort(paste0(
      "`slope` cannot be negative: an aircraft's field of view widens with ",
      "perpendicular distance, so time in view rises with it. For a circular ",
      "viewing area, which narrows, use `view_window()`."
    ))
  }

  # Linear in distance.
  trackline + distance * slope
}


#' Time in view for a vessel
#'
#' Computes the `window` for [availability()] for observers on a ship who
#' search a sector ahead of the bow. The window runs from the moment the animal
#' crosses the search radius to the moment its bearing passes the edge of the
#' search angle (abeam, for a full forward search).
#'
#' @section Which platform this fits:
#' A ship's observers search forward only. An animal behind the beam has
#' already been passed, and counting it would count it twice. This differs
#' from [view_window()], which watches a full circle, and from
#' [view_window_aerial()], whose window grows with perpendicular distance.
#'
#' @section How it works:
#' For a search radius \eqn{r}, a search half-angle \eqn{\theta} measured from
#' the bow, speed \eqn{v} and perpendicular distance \eqn{x}:
#'
#' \deqn{w(x) = \frac{\sqrt{r^2 - x^2} - x/\tan\theta}{v}}
#'
#' 1. The animal crosses the search radius when it is \eqn{\sqrt{r^2 - x^2}}
#'    ahead of the vessel (zero if \eqn{x > r}).
#' 2. Its bearing passes the edge of the search angle when it is
#'    \eqn{x/\tan\theta} ahead. At `angle = 90` this is set to exactly zero
#'    (abeam).
#' 3. The window is the track between the two, divided by speed, and `0` when
#'    the second point is further ahead than the first.
#'
#' At `angle = 90` the window is exactly half of [view_window()]'s: the same
#' chord, with the aft half never searched.
#'
#' @section Calibrating the inputs:
#' `radius` is the effective search radius, the distance at which observers
#' actually detected animals, not the reach of the binoculars. A truncation
#' distance from the fitted detection function is a better starting point than
#' a nominal one. This window does not correct for animals moving toward or
#' away from the ship before detection; that is not an availability problem.
#' See [g0()] for how components combine.
#'
#' @param radius Numeric. Effective search radius ahead of the vessel, in
#'   metres. Must be positive.
#' @param speed Numeric. Vessel speed over ground, in metres per second. Must be
#'   positive.
#' @param angle Numeric. Search half-angle from the bow, in degrees. Default
#'   `90`, a full beam-to-beam forward search. Must be above 0 and at most 90.
#' @param distance Numeric. Perpendicular distance, in metres. Default `0`, the
#'   trackline, which gives the widest window. Must not be negative.
#'
#' @return A numeric vector of windows in seconds, `0` where the animal is never
#'   inside the searched sector. `view_window_ship()` is an alias.
#'
#' @references
#' Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L.
#' and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford University
#' Press. Shipboard search geometry, forward of the beam.
#'
#' Laake, J.L., Calambokidis, J., Osmek, S.D. and Rugh, D.J. (1997) Probability
#' of detecting harbor porpoise from aerial surveys: estimating g(0). *The
#' Journal of Wildlife Management* 61:63-75. \doi{10.2307/3802415} The
#' availability formula this window feeds.
#'
#' @seealso [availability()] for what the window feeds, [view_window_aerial()]
#'   for an aircraft, and [view_window()] for a circular area.
#'
#' @examples
#' # A 6 km effective search radius at 10 knots (5.1 m/s), on the trackline
#' view_window_vessel(radius = 6000, speed = 5.1)
#'
#' # The window shrinks with perpendicular distance and is 0 beyond the radius
#' view_window_vessel(radius = 6000, speed = 5.1,
#'                    distance = c(0, 2000, 5000, 7000))
#'
#' # A narrower search arc shortens the window more at larger distances
#' view_window_vessel(radius = 6000, speed = 5.1, angle = 45,
#'                    distance = c(0, 2000, 5000))
#'
#' # A full forward search gives half the circular window
#' view_window_vessel(radius = 300, speed = 50) == view_window(300, 50) / 2
#'
#' @export
view_window_vessel <- function(radius, speed, angle = 90, distance = 0) {
  # Check inputs.
  stopifnot(is.numeric(radius), is.numeric(speed), is.numeric(distance),
            is.numeric(angle))
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
  if (any(!is.na(angle) & (angle <= 0 | angle > 90))) {
    rlang::abort(paste0(
      "`angle` is a half-angle in degrees from the bow, so it must be above 0 ",
      "and at most 90. Above 90 the search reaches abaft the beam, where an ",
      "animal has already been passed - for a view that keeps watching aft, ",
      "use `view_window()`."
    ))
  }

  # Recycle angle and distance to a common length so the subsetting below
  # pairs them element by element.
  n <- max(length(radius), length(speed), length(angle), length(distance))
  angle <- rep_len(angle, n)
  distance <- rep_len(distance, n)

  # Distance ahead at which the animal crosses the search radius.
  enters_radius <- sqrt(pmax(radius^2 - distance^2, 0))

  # Distance ahead at which its bearing passes the edge of the search angle.
  # tan(90 degrees) is not exactly infinite in floating point, so the
  # full-forward case is set to exactly zero. Subsetting rather than ifelse()
  # keeps each angle paired with its own distance.
  leaves_arc <- distance / tan(angle * pi / 180)
  leaves_arc[angle == 90] <- 0

  # Track between the two points, as time.
  pmax(enters_radius - leaves_arc, 0) / speed
}


#' @rdname view_window_vessel
#' @export
view_window_ship <- view_window_vessel


# Check that each argument has length 1 or `n` before it is recycled.
# Input: a named list of arguments and the target length `n`. Returns nothing;
# aborts naming the offending arguments, reported from the caller's frame.
# rep_len() alone would stretch a length-2 vector over 3 rows and pair values
# with rows they were never meant for.
check_recycle <- function(args, n, call = rlang::caller_env()) {
  bad <- names(args)[!lengths(args) %in% c(1L, n)]
  if (length(bad)) {
    rlang::abort(paste0(
      "Lengths must be 1 or ", n, ": ",
      paste0("`", bad, "` has ", lengths(args)[bad], collapse = ", "), "."
    ), call = call)
  }
}

# Availability for one element: a = (s + d(1 - exp(-w/d))) / (s + d), from
# Laake et al. (1997) as written out by Ganley et al. (2019).
# Input: `par` = c(surface, dive) in seconds, and `window` in seconds.
# Returns a single number, or NA if any input is NA.
# `par` is a vector so that numDeriv::grad() can differentiate it.
availability_one <- function(par, window) {
  s <- par[[1]]
  d <- par[[2]]
  if (is.na(s) || is.na(d) || is.na(window)) return(NA_real_)
  (s + d * (1 - exp(-window / d))) / (s + d)
}

# Delta-method standard error of availability for one element.
# Input: `par` = c(surface, dive), `window`, and `vcov`, the 2 x 2 covariance
# matrix of (surface, dive). Returns a single number, or NA if any input is NA
# or the computation fails.
# The gradient is numerical rather than hand-derived, as in
# effect_estimates_ddf().
availability_se <- function(par, window, vcov) {
  if (anyNA(par) || is.na(window) || anyNA(vcov)) return(NA_real_)

  # Gradient of availability with respect to (surface, dive).
  g <- try(
    numDeriv::grad(availability_one, par, window = window),
    silent = TRUE
  )
  if (inherits(g, "try-error") || anyNA(g)) return(NA_real_)

  # Variance = g' V g; return its square root if it is valid.
  v <- drop(t(g) %*% vcov %*% g)
  if (is.na(v) || v < 0) return(NA_real_)
  sqrt(v)
}
