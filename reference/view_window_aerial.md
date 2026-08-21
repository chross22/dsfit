# The time an animal stays in view of an aircraft

The window
[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
needs, for an observer looking through a side window of an aircraft: a
field of view running forward and aft, rather than a circular patch.

## Usage

``` r
view_window_aerial(trackline, speed, angle = NULL, slope = NULL, distance = 0)
```

## Arguments

- trackline:

  Time in view at the trackline, in seconds — \\\alpha\\.

- speed:

  Ground speed, in metres per second.

- angle:

  Viewing half-angle forward and aft, in **degrees** from perpendicular.
  Give this or `slope`, not both.

- slope:

  Seconds of view gained per metre of perpendicular distance,
  \\\tan\theta / v\\, if it was calibrated directly. Give this or
  `angle`.

- distance:

  Perpendicular distance, in metres.

## Value

A numeric vector of window durations in seconds.

## Why this grows with distance

The field of view is an angular wedge. An animal further off the
trackline sits inside that wedge for **longer**, because the wedge is
wider out there. This is the opposite of
[`view_window()`](https://camilleross.org/dsfit/reference/view_window.md)'s
circular geometry, and it is the case that applies to a line-transect
aerial survey.

Following Ganley et al. (2019), for a time in view \\\alpha\\ at the
trackline, a viewing half-angle \\\theta\\ forward and aft, and ground
speed \\v\\:

\$\$t(x) = \alpha + \frac{x \tan\theta}{v}\$\$

## Calibrating it, rather than deriving it

Ganley et al. (2019) obtained \\\alpha\\ and \\\theta\\ by flying
transects past a navigation buoy and timing it through the field of
view, which is worth far more than a derivation from window dimensions.
For the Cessna Skymaster over Cape Cod Bay, at 185 km/h and 228 or 304 m
altitude, they found time in view at the effective trackline to be
**51.22 s**, and a slope of about 0.03 s per metre of perpendicular
distance. Altitude made no detectable difference across the two they
flew.

Their trackline is at 100 m rather than 0 m, because the aircraft's flat
windows leave a blind spot directly beneath it — the same blind spot
that `sweep_models(left = )` and the gamma key handle at the fitting
end. Their surveys were left-truncated at 100 m for exactly that reason.

## References

Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what
there is: estimating North Atlantic right whale *Eubalaena glacialis*
local abundance. *Endangered Species Research* 38:101-113.
[doi:10.3354/esr00938](https://doi.org/10.3354/esr00938)

## See also

[`availability()`](https://camilleross.org/dsfit/reference/availability.md),
[`view_window()`](https://camilleross.org/dsfit/reference/view_window.md),
[ganley_availability](https://camilleross.org/dsfit/reference/ganley_availability.md)

## Examples

``` r
# Ganley et al.'s Skymaster, calibrated directly: 51.22 s at the trackline,
# gaining about 0.03 s per metre out
view_window_aerial(trackline = 51.22, speed = 51.4, slope = 0.03,
                   distance = c(0, 1000, 3000))
#> [1]  51.22  81.22 141.22

# The same thing from a viewing half-angle
view_window_aerial(trackline = 51.22, speed = 51.4, angle = 57,
                   distance = c(0, 1000, 3000))
#> [1]  51.22000  81.17846 141.09539
```
