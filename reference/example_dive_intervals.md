# Example surfacing and diving intervals, by month (invented)

Invented mean surfacing and diving intervals, with standard errors, for
six months. Use them to try
[`availability()`](https://camilleross.org/dsfit/reference/availability.md).
They are not measurements and must not be used to correct anything.

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

## What the numbers show

The values are invented, and are not Ganley et al. (2019)'s Cape Cod Bay
values, which ship as
[ganley_availability](https://camilleross.org/dsfit/reference/ganley_availability.md).

They copy a real pattern. Ganley et al. (2019) found right whale
availability in Cape Cod Bay varying by month between 0.27 and 0.85 (the
range in their abstract), following the depth of the copepods the whales
fed on. Deep feeding means long dives and little surface time; shallow
feeding means the whales are available most of the time. With a 12 s
window (`view_window(radius = 300, speed = 50)`),
[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
on these intervals ranges from about 0.25 to 0.85. A single constant
suited to one month is wrong by a factor of three or more for another.

Real dive intervals come from focal follows, tags or drones, and depend
on place, season and behaviour. Substitute your own.

## References

Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what
there is: estimating North Atlantic right whale *Eubalaena glacialis*
local abundance. *Endangered Species Research* 38:101-113.
[doi:10.3354/esr00938](https://doi.org/10.3354/esr00938) The monthly
pattern these invented values copy.

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

# Availability for the season, from these intervals and a 12 s window
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
