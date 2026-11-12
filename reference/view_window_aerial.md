# Time in view for an aircraft

Computes the `window` for
[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
for an observer looking through the side window of an aircraft. The
field of view is a wedge running forward and aft of perpendicular, so
time in view grows with perpendicular distance. This is the usual
geometry for an aerial line-transect survey.

## Usage

``` r
view_window_aerial(trackline, speed, angle = NULL, slope = NULL, distance = 0)
```

## Arguments

- trackline:

  Numeric. Time in view at the trackline, \\\alpha\\, in seconds. Must
  not be negative.

- speed:

  Numeric. Ground speed, in metres per second. Must be positive.

- angle:

  Numeric or `NULL` (default). Viewing half-angle \\\theta\\, in degrees
  from perpendicular: the view runs `angle` forward and `angle` aft.
  Must be above 0 and below 90. Give this or `slope`, not both.

- slope:

  Numeric or `NULL` (default). Seconds in view gained per metre of
  perpendicular distance, \\2 \tan\theta / v\\, if calibrated directly.
  Must not be negative. Give this or `angle`, not both.

- distance:

  Numeric. Perpendicular distance, in metres. Default `0`. Must not be
  negative.

## Value

A numeric vector of windows in seconds.

## How it works

Following the linear form fitted by Ganley et al. (2019), for a time in
view \\\alpha\\ at the trackline, a viewing half-angle \\\theta\\
forward of perpendicular and the same aft, and ground speed \\v\\:

\$\$t(x) = \alpha + \frac{2 x \tan\theta}{v}\$\$

A point at perpendicular distance \\x\\ enters the wedge \\x
\tan\theta\\ ahead of abeam and leaves \\x \tan\theta\\ behind, so it is
in view for an extra \\2 x \tan\theta\\ of track.

1.  Exactly one of `angle` and `slope` must be given.

2.  If `angle` is given, the slope is computed as \\2 \tan\theta / v\\
    seconds per metre.

3.  The window is `trackline + distance * slope`.

## Calibrating the inputs

Measured values are better than derived ones. Ganley et al. (2019) flew
past a navigation buoy and timed it through the field of view. For a
Cessna Skymaster at 185 km/h and 228 or 304 m altitude, they found 51.22
s in view at the effective trackline and a slope of about 0.03 s per
metre. The two altitudes did not differ detectably.

Their effective trackline is at 100 m, not 0 m, because the aircraft's
flat windows leave a blind spot beneath it; their data were
left-truncated at 100 m for that reason. At the fitting end the same
blind spot is handled by `sweep_models(left = )` or the gamma key.

## References

Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what
there is: estimating North Atlantic right whale *Eubalaena glacialis*
local abundance. *Endangered Species Research* 38:101-113.
[doi:10.3354/esr00938](https://doi.org/10.3354/esr00938)

## See also

[`availability()`](https://camilleross.org/dsfit/reference/availability.md),
[`view_window()`](https://camilleross.org/dsfit/reference/view_window.md),
[`view_window_vessel()`](https://camilleross.org/dsfit/reference/view_window_vessel.md),
[ganley_availability](https://camilleross.org/dsfit/reference/ganley_availability.md)

## Examples

``` r
# Ganley et al.'s Skymaster, calibrated directly: 51.22 s at the trackline,
# plus about 0.03 s per metre
view_window_aerial(trackline = 51.22, speed = 51.4, slope = 0.03,
                   distance = c(0, 1000, 3000))
#> [1]  51.22  81.22 141.22

# The same from a half-angle: 2 * tan(37.6 deg) / 51.4 is about 0.03 s/m
view_window_aerial(trackline = 51.22, speed = 51.4, angle = 37.6,
                   distance = c(0, 1000, 3000))
#> [1]  51.22000  81.18512 141.11537
```
