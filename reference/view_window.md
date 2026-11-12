# Time in view for a circular viewing area

Computes the `window` for
[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
when the observers watch a circular area around the platform, fore and
aft alike. The window is the time a point on the water spends inside
that circle as the platform passes. For an aircraft side window use
[`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md);
for a ship's forward search use
[`view_window_vessel()`](https://camilleross.org/dsfit/reference/view_window_vessel.md).

## Usage

``` r
view_window(radius, speed, distance = 0)
```

## Arguments

- radius:

  Numeric. Radius of the viewing area, in metres. Must be positive.

- speed:

  Numeric. Platform ground speed, in metres per second. Must be
  positive.

- distance:

  Numeric. Perpendicular distance, in metres. Default `0`, the
  trackline, which gives the widest window. Must not be negative.

## Value

A numeric vector of windows in seconds, `0` where `distance` exceeds
`radius`.

## How it works

For a viewing radius \\r\\ and platform speed \\v\\, a point at
perpendicular distance \\x\\ is crossed along a chord of length
\\2\sqrt{r^2 - x^2}\\, so

\$\$w(x) = \frac{2\sqrt{r^2 - x^2}}{v}\$\$

The window is widest on the trackline, shrinks with distance, and is `0`
beyond the radius.

## Which platform this fits

Only a platform that watches a circular area. An aircraft observer
looking through a side window sees a wedge running forward and aft, so
time in view grows with perpendicular distance instead of shrinking;
this function gets that direction wrong. Ganley et al. (2019) measured
time in view rising with distance for their aircraft, at about 0.03 s
per metre. A window measured in the field, as they did by timing a buoy
through the field of view, is better than any of these formulas.

## References

Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what
there is: estimating North Atlantic right whale *Eubalaena glacialis*
local abundance. *Endangered Species Research* 38:101-113.
[doi:10.3354/esr00938](https://doi.org/10.3354/esr00938) Time in view
measured for an aircraft, rising with perpendicular distance.

## See also

[`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md)
for an aircraft,
[`view_window_vessel()`](https://camilleross.org/dsfit/reference/view_window_vessel.md)
for a ship, and
[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
for what the window feeds.

## Examples

``` r
# A 300 m viewing radius at 50 m/s, on the trackline
view_window(radius = 300, speed = 50)
#> [1] 12

# The window narrows away from the trackline and is 0 beyond the radius
view_window(radius = 300, speed = 50, distance = c(0, 150, 290, 400))
#> [1] 12.000000 10.392305  3.072458  0.000000
```
