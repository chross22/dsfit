# Example surfacing and diving intervals, by month

A worked example for
[`availability()`](https://camilleross.org/dsfit/reference/availability.md):
mean surfacing and diving intervals over a six-month season, with
standard errors.

## Usage

``` r
example_dive_intervals
```

## Format

A tibble with 6 rows and 5 columns:

- month:

  Month, as an ordered factor from December to May.

- surface:

  Mean surfacing interval, seconds.

- dive:

  Mean diving interval, seconds.

- se_surface:

  Standard error of `surface`, seconds.

- se_dive:

  Standard error of `dive`, seconds.

## These numbers are invented

They are not measurements. In particular they are **not** Ganley et al.
(2019)'s Cape Cod Bay values, which are not open access. Do not use them
to correct anything.

What is real is the pattern they were built to show. Ganley et al.
(2019) found right whale availability in Cape Cod Bay varying by month
between **0.27 and 0.85**, tracking the depth of the copepod layer the
whales were feeding on: feeding deep means long dives and little time at
the surface, feeding shallow means available most of the time. These
intervals reproduce that shape, and running
[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
over them spans roughly 0.25 to 0.85 — which is the point of shipping
them. A single representative availability is a real number for one
month and wrong by a factor of two for others.

Real dive parameters come from focal follows, tagging, or drone work,
and are specific to a place, a season and a behaviour. Substitute your
own.

## References

Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what
there is: estimating North Atlantic right whale *Eubalaena glacialis*
local abundance. *Endangered Species Research* 38:101-113.
[doi:10.3354/esr00938](https://doi.org/10.3354/esr00938) The monthly
pattern these values imitate, and the real measurements they are not.

## See also

[`availability()`](https://camilleross.org/dsfit/reference/availability.md),
[`view_window()`](https://camilleross.org/dsfit/reference/view_window.md)

## Examples

``` r
example_dive_intervals
#> # A tibble: 6 × 5
#>   month surface  dive se_surface se_dive
#>   <fct>   <dbl> <dbl>      <dbl>   <dbl>
#> 1 Dec        35   155          6      24
#> 2 Jan        40   140          7      21
#> 3 Feb        50   120          8      18
#> 4 Mar        65    95         10      14
#> 5 Apr        95    60         14       9
#> 6 May       130    35         19       6

# The season's availability, from these intervals and a platform's geometry
availability(
  surface    = example_dive_intervals$surface,
  dive       = example_dive_intervals$dive,
  window     = view_window(radius = 300, speed = 50),
  se_surface = example_dive_intervals$se_surface,
  se_dive    = example_dive_intervals$se_dive,
  key        = as.character(example_dive_intervals$month)
)
#> # A tibble: 6 × 4
#>   key   component    value     se
#>   <chr> <chr>        <dbl>  <dbl>
#> 1 Dec   availability 0.245 0.0388
#> 2 Jan   availability 0.286 0.0431
#> 3 Feb   availability 0.361 0.0483
#> 4 Mar   availability 0.477 0.0525
#> 5 Apr   availability 0.683 0.0481
#> 6 May   availability 0.849 0.0339
```
