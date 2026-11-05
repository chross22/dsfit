# The time an animal stays in view of a vessel

The window
[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
needs, for observers searching a forward arc from a ship: a sector ahead
of the bow, which a point on the water enters when it comes inside the
search angle and leaves when it draws abeam.

## Usage

``` r
view_window_vessel(radius, speed, angle = 90, distance = 0)

view_window_ship(radius, speed, angle = 90, distance = 0)
```

## Arguments

- radius:

  Effective search radius ahead of the vessel, in metres.

- speed:

  Vessel speed over ground, in metres per second.

- angle:

  Search half-angle from the bow, in **degrees**. `90` (default) is a
  full beam-to-beam forward search; a narrower arc is narrower.

- distance:

  Perpendicular distance, in metres. `0` (default) is on the trackline,
  which gives the widest window.

## Value

A numeric vector of window durations in seconds, `0` where the animal is
never inside the searched sector.

## Why a ship is neither of the other two

[`view_window()`](https://camilleross.org/dsfit/reference/view_window.md)
is a circular patch, watched fore and aft alike, so its window is the
full chord.
[`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md)
is a wedge through a side window, so its window *grows* with
perpendicular distance. A shipboard observer does neither: the search is
**forward only**. Behind the beam there is nothing to see that has not
already been counted, and a survey that counted it would be counting the
same animal twice.

So the window runs from the moment the animal enters the searched sector
to the moment it draws abeam, and no further. For a search radius \\r\\,
a search half-angle \\\theta\\ measured from the bow, a speed \\v\\ and
a perpendicular distance \\x\\:

\$\$w(x) = \frac{\sqrt{r^2 - x^2} - x/\tan\theta}{v}\$\$

The first term is where the animal crosses the radius; the second is
where it crosses the edge of the search angle. Whichever it meets later
is where the window opens, so the window is zero when the second exceeds
the first — an animal too far off the bow to be searched before it
passes.

At `angle = 90` the sector is the full forward half-disc,
\\x/\tan\theta\\ is zero, and the window is exactly **half**
[`view_window()`](https://camilleross.org/dsfit/reference/view_window.md)'s:
same chord, half of it aft of the beam and never searched.

## Calibrating it, rather than deriving it

The same caution applies as for an aircraft. `radius` is an *effective*
search radius — how far the observers were actually detecting, not how
far the big-eyes reach — and it is worth measuring rather than assuming.
A truncation distance from the fitted detection function is a better
starting point than a nominal one.

Nothing here corrects for a ship's own effect on the animals. Vessel
avoidance and attraction both move animals before they are counted, and
neither is an availability problem: availability asks whether an animal
that was there surfaced in time to be seen. See
[`g0()`](https://camilleross.org/dsfit/reference/g0.md) for where the
components combine.

## References

Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers,
D.L. and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford
University Press. Shipboard line transect search geometry, and why the
search is forward of the beam.

Laake, J.L., Calambokidis, J., Osmek, S.D. and Rugh, D.J. (1997)
Probability of detecting harbor porpoise from aerial surveys: estimating
g(0). *The Journal of Wildlife Management* 61:63-75.
[doi:10.2307/3802415](https://doi.org/10.2307/3802415) The availability
formula this window feeds.

## See also

[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
for what the window feeds,
[`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md)
for an aircraft's forward-and-aft wedge, and
[`view_window()`](https://camilleross.org/dsfit/reference/view_window.md)
for a circular patch.

## Examples

``` r
# A 6 km effective search radius at 10 knots (5.1 m/s), on the trackline
view_window_vessel(radius = 6000, speed = 5.1)
#> [1] 1176.471

# The window closes as the animal approaches the edge of the sector
view_window_vessel(radius = 6000, speed = 5.1,
                   distance = c(0, 2000, 5000, 7000))
#> [1] 1176.4706 1109.1871  650.3186    0.0000

# A narrower search arc loses the animals furthest off the bow first
view_window_vessel(radius = 6000, speed = 5.1, angle = 45,
                   distance = c(0, 2000, 5000))
#> [1] 1176.4706 1109.1871  650.3186

# A full forward search is half a circular patch: the aft half is never
# searched
view_window_vessel(radius = 300, speed = 50) == view_window(300, 50) / 2
#> [1] TRUE
```
