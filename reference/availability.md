# Availability of an animal at the surface

Computes availability: the probability that an animal is at the surface,
and so can be seen, while it is within the observers' view. Use it to
build the availability component of a `g(0)` correction, one row per key
(month, platform, year), and pass the result to
[`g0()`](https://camilleross.org/dsfit/reference/g0.md).

## Usage

``` r
availability(
  surface,
  dive,
  window,
  se_surface = NULL,
  se_dive = NULL,
  cov_surface_dive = 0,
  key = NULL
)
```

## Arguments

- surface:

  Numeric. Mean surfacing interval, in seconds. Must be positive. A
  vector gives one result per element.

- dive:

  Numeric. Mean diving interval, in seconds. Must be positive. Length 1
  or the length of the longest input.

- window:

  Numeric. Time the animal is in view, in seconds. Must not be negative.
  Derive it with
  [`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md),
  [`view_window_vessel()`](https://camilleross.org/dsfit/reference/view_window_vessel.md)
  or
  [`view_window()`](https://camilleross.org/dsfit/reference/view_window.md).
  Length 1 or the length of the longest input.

- se_surface, se_dive:

  Numeric or `NULL` (default). Standard errors of `surface` and `dive`,
  in seconds. Give both or neither: propagating one and dropping the
  other would understate the total variance.

- cov_surface_dive:

  Numeric. Covariance between the two interval means, in seconds
  squared. Default `0`, which is right when they come from separate
  follows and optimistic otherwise. Used only with standard errors.

- key:

  Optional vector of labels (months, platforms, years), one per result.
  Must have the same length as the result. Stored as character.

## Value

A tibble with one row per element of the longest input, and columns:

- key:

  The label from `key`, or `NA`.

- component:

  Always `"availability"`.

- value:

  Availability, between 0 and 1.

- se:

  Standard error of `value`, or `NA` without `se_surface` and `se_dive`.

This is the component shape
[`g0()`](https://camilleross.org/dsfit/reference/g0.md) accepts, so
availability and perception rows can be stacked.

## Why availability is computed, not estimated

An animal that is submerged while the platform passes is missed at every
perpendicular distance equally. Availability therefore scales the
detection function \\g(x)\\ without changing its shape, and the
distances carry no information about it. The inputs come from outside
the survey: `surface` and `dive` from focal follows, tags or drones, and
`window` from the platform's geometry and speed.

## How it works

Following Laake et al. (1997), for a mean surfacing interval \\E(s)\\, a
mean diving interval \\E(d)\\, and a window \\w\\ (the time the animal
is in view):

\$\$a = \frac{E(s) + E(d)\left(1 - e^{-w/E(d)}\right)}{E(s) + E(d)}\$\$

At \\w = 0\\, \\a\\ is \\E(s)/(E(s) + E(d))\\, the proportion of time at
the surface. As \\w\\ grows, \\a\\ approaches 1. The middle term is the
chance that an animal that was submerged when the window opened surfaces
before it closes.

1.  The inputs are checked and recycled to a common length. Each must
    have length 1 or the common length.

2.  `surface` and `dive` must be positive and `window` must not be
    negative.

3.  `se_surface` and `se_dive` must be given together or not at all.

4.  \\a\\ is computed for each element.

5.  If standard errors were given, the standard error of \\a\\ is
    computed by the delta method: the gradient of \\a\\ with respect to
    \\(E(s), E(d))\\ is taken numerically and combined with their
    covariance matrix, as in
    [`effect_estimates_ddf()`](https://camilleross.org/dsfit/reference/effect_estimates_ddf.md).

Without `se_surface` and `se_dive` the standard error is `NA`. The
function does not invent a precision, because the CV of availability
often dominates the CV of abundance. If the raw focal follows are
available, bootstrapping them (as Ganley et al. 2019 did) avoids
assuming the interval means are jointly normal; the result can be passed
as `se_surface` and `se_dive`, or the component rows can be built
directly.

## Use one value per key, not one constant

Availability varies. Ganley et al. (2019) found right whale availability
in Cape Cod Bay ranging from 0.27 in January to 0.91 in April, following
the depth of the copepods the whales fed on (see
[ganley_availability](https://camilleross.org/dsfit/reference/ganley_availability.md)).
A single constant suited to one month can be wrong by a factor of three
for another. Because it scales every candidate model equally, model
selection cannot detect the error. Compute one value per key and pass a
vector.

This function gives availability only. Perception (whether an animal at
the surface and in view is seen) is a separate component and needs a
second observer team to estimate.
[`g0()`](https://camilleross.org/dsfit/reference/g0.md) warns if it is
missing.

## References

Laake, J.L., Calambokidis, J., Osmek, S.D. and Rugh, D.J. (1997)
Probability of detecting harbor porpoise from aerial surveys: estimating
g(0). *The Journal of Wildlife Management* 61:63-75.
[doi:10.2307/3802415](https://doi.org/10.2307/3802415) Source of the
formula.

Laake, J.L. and Borchers, D.L. (2004) Methods for incomplete detection
at distance zero. In *Advanced Distance Sampling*, pp. 108-189. Oxford
University Press. Why perception and availability need separate data.

Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what
there is: estimating North Atlantic right whale *Eubalaena glacialis*
local abundance. *Endangered Species Research* 38:101-113.
[doi:10.3354/esr00938](https://doi.org/10.3354/esr00938) Focal follows
and monthly availability for right whales.

Roberts, J.J., Yack, T.M., Fujioka, E., Halpin, P.N., Baumgartner, M.F.
and others (2024) North Atlantic right whale density surface model for
the US Atlantic evaluated with passive acoustic monitoring. *Marine
Ecology Progress Series* 732:167-192.
[doi:10.3354/meps14547](https://doi.org/10.3354/meps14547) Corrections
for perception and availability by platform, team and conditions.

## See also

[`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md),
[`view_window_vessel()`](https://camilleross.org/dsfit/reference/view_window_vessel.md)
and
[`view_window()`](https://camilleross.org/dsfit/reference/view_window.md)
for the window; [`g0()`](https://camilleross.org/dsfit/reference/g0.md)
for combining components;
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
for the detection function fit, which this does not enter.

## Examples

``` r
# With a zero window, availability is the proportion of time at the surface
availability(surface = 60, dive = 240, window = 0)
#> # A tibble: 1 × 4
#>   key   component    value    se
#>   <chr> <chr>        <dbl> <dbl>
#> 1 NA    availability   0.2    NA

# A longer window raises it: some submerged animals surface in time
availability(surface = 60, dive = 240, window = 30)
#> # A tibble: 1 × 4
#>   key   component    value    se
#>   <chr> <chr>        <dbl> <dbl>
#> 1 NA    availability 0.294    NA

# One row per month
availability(
  surface = c(45, 60, 90),
  dive    = c(300, 240, 150),
  window  = 24,
  key     = c("Feb", "Mar", "Apr")
)
#> # A tibble: 3 × 4
#>   key   component    value    se
#>   <chr> <chr>        <dbl> <dbl>
#> 1 Feb   availability 0.197    NA
#> 2 Mar   availability 0.276    NA
#> 3 Apr   availability 0.467    NA

# With standard errors on the intervals
availability(surface = 60, dive = 240, window = 24,
             se_surface = 8, se_dive = 25)
#> # A tibble: 1 × 4
#>   key   component    value     se
#>   <chr> <chr>        <dbl>  <dbl>
#> 1 NA    availability 0.276 0.0297
```
