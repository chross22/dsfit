# The time an animal stays within view of the platform

Turns platform geometry into the `window`
[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
needs: how long a point on the water remains inside the observers' field
of view as the platform passes.

## Usage

``` r
view_window(radius, speed, distance = 0)
```

## Arguments

- radius:

  Radius of the viewing area, in metres.

- speed:

  Platform ground speed, in metres per second.

- distance:

  Perpendicular distance, in metres. `0` (default) is on the trackline,
  which gives the widest window.

## Value

A numeric vector of window durations in seconds, `0` where `distance`
exceeds `radius`.

## The geometry

For a circular field of view of radius \\r\\ and a platform travelling
at speed \\v\\, a point at perpendicular distance \\x\\ is crossed along
a chord of length \\2\sqrt{r^2 - x^2}\\, so

\$\$w(x) = \frac{2\sqrt{r^2 - x^2}}{v}\$\$

The window therefore shrinks with perpendicular distance, and is zero
beyond the edge of the view — an animal outside it is never available,
which this returns as `0` rather than as an error.

## This is the wrong geometry for most aerial surveys

A circular viewing area is what a platform watching a fixed patch has.
An **aircraft observer looking through a side window does not have
one**: the field of view is an angular wedge running forward and aft, so
an animal further off the trackline sits in that wedge *longer*, not
less. The window **grows** with perpendicular distance, and this
function has the sign of that effect backwards for such a platform.

Ganley et al. (2019) measured it for the Skymaster flown over Cape Cod
Bay and found time in view rising from about 50 s near the trackline to
about 130 s at 3 km, a slope of roughly 0.03 s per metre. Use
[`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md)
there. This function is for a genuinely circular view, and is kept
because that case exists — not because it is the default worth reaching
for.

Measuring the window, as Ganley et al. did by timing a navigation buoy
through the field of view, beats deriving it from either formula.

## References

Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what
there is: estimating North Atlantic right whale *Eubalaena glacialis*
local abundance. *Endangered Species Research* 38:101-113.
[doi:10.3354/esr00938](https://doi.org/10.3354/esr00938) Time in view
measured for an aerial platform, and why it rises rather than falls with
perpendicular distance.

## See also

[`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md)
for an aircraft's forward-and-aft field of view, which is the usual case
here.
[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
for what the window feeds.

## Examples

``` r
# A 300 m viewing radius at 50 m/s, on the trackline
view_window(radius = 300, speed = 50)
#> [1] 12

# The window narrows away from the trackline, and closes at the edge
view_window(radius = 300, speed = 50, distance = c(0, 150, 290, 400))
#> [1] 12.000000 10.392305  3.072458  0.000000
```
