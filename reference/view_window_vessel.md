# Time in view for a vessel

Computes the `window` for
[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
for observers on a ship who search a sector ahead of the bow. The window
runs from the moment the animal crosses the search radius to the moment
its bearing passes the edge of the search angle (abeam, for a full
forward search).

## Usage

``` r
view_window_vessel(radius, speed, angle = 90, distance = 0)

view_window_ship(radius, speed, angle = 90, distance = 0)
```

## Arguments

- radius:

  Numeric. Effective search radius ahead of the vessel, in metres. Must
  be positive.

- speed:

  Numeric. Vessel speed over ground, in metres per second. Must be
  positive.

- angle:

  Numeric. Search half-angle from the bow, in degrees. Default `90`, a
  full beam-to-beam forward search. Must be above 0 and at most 90.

- distance:

  Numeric. Perpendicular distance, in metres. Default `0`, the
  trackline, which gives the widest window. Must not be negative.

## Value

A numeric vector of windows in seconds, `0` where the animal is never
inside the searched sector. `view_window_ship()` is an alias.

## Which platform this fits

A ship's observers search forward only. An animal behind the beam has
already been passed, and counting it would count it twice. This differs
from
[`view_window()`](https://camilleross.org/dsfit/reference/view_window.md),
which watches a full circle, and from
[`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md),
whose window grows with perpendicular distance.

## How it works

For a search radius \\r\\, a search half-angle \\\theta\\ measured from
the bow, speed \\v\\ and perpendicular distance \\x\\:

\$\$w(x) = \frac{\sqrt{r^2 - x^2} - x/\tan\theta}{v}\$\$

1.  The animal crosses the search radius when it is \\\sqrt{r^2 - x^2}\\
    ahead of the vessel (zero if \\x \> r\\).

2.  Its bearing passes the edge of the search angle when it is
    \\x/\tan\theta\\ ahead. At `angle = 90` this is set to exactly zero
    (abeam).

3.  The window is the track between the two, divided by speed, and `0`
    when the second point is further ahead than the first.

At `angle = 90` the window is exactly half of
[`view_window()`](https://camilleross.org/dsfit/reference/view_window.md)'s:
the same chord, with the aft half never searched.

## Calibrating the inputs

`radius` is the effective search radius, the distance at which observers
actually detected animals, not the reach of the binoculars. A truncation
distance from the fitted detection function is a better starting point
than a nominal one. This window does not correct for animals moving
toward or away from the ship before detection; that is not an
availability problem. See
[`g0()`](https://camilleross.org/dsfit/reference/g0.md) for how
components combine.

## References

Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers,
D.L. and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford
University Press. Shipboard search geometry, forward of the beam.

Laake, J.L., Calambokidis, J., Osmek, S.D. and Rugh, D.J. (1997)
Probability of detecting harbor porpoise from aerial surveys: estimating
g(0). *The Journal of Wildlife Management* 61:63-75.
[doi:10.2307/3802415](https://doi.org/10.2307/3802415) The availability
formula this window feeds.

## See also

[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
for what the window feeds,
[`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md)
for an aircraft, and
[`view_window()`](https://camilleross.org/dsfit/reference/view_window.md)
for a circular area.

## Examples

``` r
# A 6 km effective search radius at 10 knots (5.1 m/s), on the trackline
view_window_vessel(radius = 6000, speed = 5.1)
#> [1] 1176.471

# The window shrinks with perpendicular distance and is 0 beyond the radius
view_window_vessel(radius = 6000, speed = 5.1,
                   distance = c(0, 2000, 5000, 7000))
#> [1] 1176.4706 1109.1871  650.3186    0.0000

# A narrower search arc shortens the window more at larger distances
view_window_vessel(radius = 6000, speed = 5.1, angle = 45,
                   distance = c(0, 2000, 5000))
#> [1] 1176.4706  717.0302    0.0000

# A full forward search gives half the circular window
view_window_vessel(radius = 300, speed = 50) == view_window(300, 50) / 2
#> [1] TRUE
```
