# Right whale availability in Cape Cod Bay, by month

Table S1 of Ganley et al. (2019): mean dive and surface intervals,
percent surface time, and median availability, from 86 focal follows of
North Atlantic right whales in Cape Cod Bay during 2016 and 2017. Unlike
[example_dive_intervals](https://camilleross.org/dsfit/reference/example_dive_intervals.md),
these are real measurements.

## Usage

``` r
ganley_availability
```

## Format

A tibble with 4 rows and 8 columns:

- month:

  Month, as an ordered factor from January to April.

- percent_surface_time:

  Percent of time at the surface. See above for why this is not the
  ratio of the two interval columns.

- mean_dive:

  Mean diving interval, seconds.

- mean_surface:

  Mean surfacing interval, seconds.

- availability:

  Median availability, \\a(x)\\.

- availability_variance:

  As labelled in the source; see above.

- n_follows:

  Focal follows behind each row. January's 7 against April's 48 is worth
  noticing — the month with the lowest availability is also the one
  measured least well, because the weather that keeps the food deep also
  keeps observers on the ground.

- hours_followed:

  Total follow time, hours.

The pooled "All sightings" row of Table S1 is kept as the
`all_sightings` attribute rather than a fifth row, so it cannot be
summed or plotted alongside the months by accident.

## Source

Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what
there is: estimating North Atlantic right whale *Eubalaena glacialis*
local abundance. *Endangered Species Research* 38:101-113.
[doi:10.3354/esr00938](https://doi.org/10.3354/esr00938) , Table S1.
Open access under CC-BY.

## Details

Availability runs from **0.27 in January to 0.91 in April**, as the
copepods the whales feed on move up through the water column. That
threefold swing across one season is the case against a constant `g(0)`,
in measurements.

## Three things not to do with these columns

- Do not treat `percent_surface_time` as \\E(s)/(E(s)+E(d))\\:

  It is not. January is listed at 16% with a mean surface interval of 48
  s and a mean dive of 533 s, and \\48/(48+533)\\ is 8.3%. The gap
  appears every month and in both directions — April's intervals give
  90% against a listed 55%. The percentage is evidently a mean of
  per-follow percentages while the interval columns are means of
  intervals: a mean of ratios against a ratio of means. So it is **not**
  the instantaneous availability, and using it as one would be wrong by
  up to 35 percentage points.

- Do not expect
  [`availability()`](https://camilleross.org/dsfit/reference/availability.md)
  to reproduce `availability`:

  The reported figure is a median over bootstrap replicates, and
  \\a(x)\\ varies with distance through the time in view. Working
  backwards from the tabulated means, January's 0.27 implies a window
  near 122 s and April's 0.91 one near 8 s — they are not evaluated at a
  common window. The values are measurements to be used, not outputs to
  be recomputed.

- Do not feed `availability_variance` to
  [`g0()`](https://camilleross.org/dsfit/reference/g0.md) without
  deciding what it is:

  It is labelled a variance and runs 0.04–0.09. Read literally,
  January's standard error is \\\sqrt{0.04} = 0.2\\ against an estimate
  of 0.27 — a CV of 74%, which would dominate any correction it entered.
  Read as a standard error instead, the CV is 15%. Neither matches the
  very small error bars in the paper's Fig. 4A. It is shipped under the
  paper's own label and deliberately not converted.

## The platform

Cessna 336/337 Skymaster at 185 km/h and 228 or 304 m altitude,
surveying January to May. Time in view at the effective trackline was
**51.22 s**, rising with perpendicular distance — see
[`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md).
Their trackline is at 100 m rather than 0 m because the aircraft's flat
windows leave a blind spot beneath it, and the surveys were
left-truncated there accordingly — the same blind spot
`sweep_models(left = )` and the gamma key handle at the fitting end.

## Two inconsistencies in the source

Recorded so they are not mistaken for transcription errors. February's
sample size is 9 in Table S1, summing to the 86 its own total row gives,
but 10 in the Fig. 2 caption, summing to the 87 the Methods states. And
April's availability is 0.91 here and in Section 3.1, while the abstract
gives the seasonal range as 0.27–0.85. Table S1 is followed here, being
the tabulated source.

## See also

[`availability()`](https://camilleross.org/dsfit/reference/availability.md),
[`g0()`](https://camilleross.org/dsfit/reference/g0.md),
[`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md),
[ganley_detection](https://camilleross.org/dsfit/reference/ganley_detection.md)

## Examples

``` r
ganley_availability
#> # A tibble: 4 × 8
#>   month percent_surface_time mean_dive mean_surface availability
#>   <fct>                <dbl>     <dbl>        <dbl>        <dbl>
#> 1 Jan                     16       533           48         0.27
#> 2 Feb                     34       256          226         0.52
#> 3 Mar                     31       219           67         0.52
#> 4 Apr                     55        88          801         0.91
#> # ℹ 3 more variables: availability_variance <dbl>, n_follows <int>,
#> #   hours_followed <dbl>

# The measured seasonal swing, as a g(0) component. A standard error has to
# be decided on first - see above on the variance column - so this uses a
# deliberately explicit placeholder rather than a silent conversion.
avail <- data.frame(
  key = as.character(ganley_availability$month),
  component = "availability",
  value = ganley_availability$availability,
  se = sqrt(ganley_availability$availability_variance)
)
suppressWarnings(g0(avail))
#> <dsfit_g0>
#>   assumed 1:   perception
#>   components:
#>     availability  source not recorded
#> 
#>  key   g0     se    cv
#>  Jan 0.27 0.2000 0.741
#>  Feb 0.52 0.3000 0.577
#>  Mar 0.52 0.2449 0.471
#>  Apr 0.91 0.2236 0.246
#> 
#>   Divide abundance by g0, and propagate cv. This package does not
#>   apply it: correction happens at the abundance step.

attr(ganley_availability, "all_sightings")$availability
#> [1] 0.63
```
