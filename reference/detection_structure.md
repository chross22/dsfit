# What a set of detections can and cannot support

Reads a table of detections and reports which analyses it admits, and
which it does not, with the reason. Nothing is fitted and nothing is
estimated.

## Usage

``` r
detection_structure(data)
```

## Arguments

- data:

  A data frame with one row per detection, as
  [`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
  takes: `distance` for exact distances, or `distbegin` and `distend`
  for binned ones.

## Value

An object of class `dsfit_structure`: a list with `table` (one row per
question, with `check`, `supported` and `detail`) and `summary` (the
counts and distance type behind it). `supported` is `TRUE`, `FALSE`, or
`NA` where the structure is present but incomplete.

## Why this exists

Most of what decides whether an analysis is possible is structural, and
none of it is announced by the data. Whether distances are exact or
binned decides which goodness-of-fit test can be computed. Whether a
survey ran two independent observer teams decides whether perception
bias is estimable at all — and that is a property of the survey
programme, not of the archive its data ends up in, so a pooled extract
may or may not carry it.

The failure this guards against is not an error but a silence: fitting a
single-observer dataset and reporting the result as though perception
had been handled. That produces a number, and the number is wrong by a
factor. Asking here turns "you had to know that" into something the
package says.

## What it does not tell you

That an analysis is *supported* is a statement about structure, not
about whether it is a good idea. Enough detections to fit a detection
function is not enough detections to fit one well, and a covariate being
present is not a reason to put it in a model.

Availability is reported as unsupported for every table, which is not a
defect in any particular dataset: it cannot be estimated from sighting
distances by construction, because an animal submerged for the whole
pass is missed at every distance equally and leaves no signature. It is
computed from dive data instead — see
[`availability()`](https://camilleross.org/dsfit/reference/availability.md).

## See also

[`prepare_distance_data()`](https://camilleross.org/dsfit/reference/prepare_distance_data.md),
which enforces what this only reports.
[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
and [`g0()`](https://camilleross.org/dsfit/reference/g0.md) for the
components this cannot supply.

## Examples

``` r
set.seed(1)
d <- data.frame(
  object = 1:200,
  distance = abs(rnorm(200, 0, 120)),
  beaufort = sample(0:4, 200, replace = TRUE)
)
detection_structure(d)
#> <dsfit_structure>
#>   200 rows, 200 exact distances
#>   nearest detection at 0.1326 - a blind spot beneath the platform would show here
#> 
#>   can:
#>     fit a detection function        200 exact distances; at or above the 60-80 usually suggested
#>     test fit with Cramer-von Mises  exact distances have an empirical distribution to test
#>     test fit with chi-square        over cutpoints mrds chooses, which makes it the weaker test here
#>     fit covariate models            candidates: beaufort
#> 
#>   cannot:
#>     estimate perception bias       no double-observer structure: no `observer` or `detected` column. Perception needs two independent teams and cannot be recovered from a single-observer survey at any sample size
#>     estimate availability          not estimable from distances by construction: a submerged animal is missed at every distance equally. Compute it from dive data with availability()
#>     carry group size to abundance  no `size` column; abundance would count groups, not individuals
#> 
#>   Structure only. That something is supported is not a reason to do it.

# Binned distances admit a different goodness-of-fit test
b <- data.frame(object = 1:60, distbegin = rep(c(0, 100, 200), 20),
                distend = rep(c(100, 200, 300), 20))
detection_structure(b)
#> <dsfit_structure>
#>   60 rows, 60 binned distances
#> 
#>   can:
#>     fit a detection function  60 binned distances; at or above the 60-80 usually suggested
#>     test fit with chi-square  over the survey's own bins
#> 
#>   cannot:
#>     test fit with Cramer-von Mises  needs exact distances; binned fits have no empirical distribution
#>     fit covariate models            no columns beyond the structural ones vary
#>     estimate perception bias        no double-observer structure: no `observer` or `detected` column. Perception needs two independent teams and cannot be recovered from a single-observer survey at any sample size
#>     estimate availability           not estimable from distances by construction: a submerged animal is missed at every distance equally. Compute it from dive data with availability()
#>     carry group size to abundance   no `size` column; abundance would count groups, not individuals
#> 
#>   Structure only. That something is supported is not a reason to do it.
```
