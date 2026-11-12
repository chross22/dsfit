# Report which analyses a detection table supports

Reads a table of detections and reports which analyses its structure
allows, and the reason for each answer. Nothing is fitted or estimated.
Use it on a new dataset to see what can be done with it before choosing
models.

## Usage

``` r
detection_structure(data)
```

## Arguments

- data:

  Data frame, one row per detection, in the form
  [`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
  takes: `distance` for exact distances, or `distbegin` and `distend`
  for binned ones.

## Value

An object of class `dsfit_structure`, a list with:

- `table`:

  A tibble with one row per check and columns `check` (the question),
  `supported` (`TRUE`, `FALSE`, or `NA` where the structure is ambiguous
  or incomplete) and `detail` (the reason).

- `summary`:

  A list: `n` (rows), `n_distances` (non-missing distances), `type`
  (`"exact"`, `"binned"`, `"both"` or `"none"`), `covariates` (candidate
  covariate names), `duplicates` (named list of column-name pairs naming
  one variable) and `nearest` (smallest distance, or smallest
  `distbegin` for binned data; `NA` if none).

## How it works

1.  Classify the distances as `"exact"` (a non-missing `distance`),
    `"binned"` (a non-missing `distbegin`, with `distend` present),
    `"both"` or `"none"`.

2.  Answer seven questions, one row each (table below).

3.  Find columns that name the same variable twice (see below).

4.  Record the nearest detection distance, where a blind spot beneath
    the platform would show.

The seven checks, in the order they appear in `table`:

|  |  |
|----|----|
| Check | Supported when |
| fit a detection function | Exact or binned distances. `NA` if both are present. The detail compares the count with the 40 workable and 60-80 suggested detections (Buckland et al. 2001). |
| test fit with Cramer-von Mises | Exact distances only. Binned data have no empirical distribution to test. |
| test fit with chi-square | Exact or binned distances. For binned data it uses the survey's bins; for exact data, cutpoints `mrds` chooses. |
| fit covariate models | At least one candidate covariate column (see below). |
| estimate perception bias | A double-observer structure (see below). |
| estimate availability | Never (see below). |
| carry group size to abundance | A `size` column. |

## Why structure matters

Whether an analysis is possible depends on the survey design, and the
data do not announce it. Exact or binned distances decide which
goodness-of-fit test can be computed. Two independent observer teams are
needed to estimate *perception bias* (animals at the surface but missed
by the observers). That depends on the survey programme, so a pooled
extract may or may not carry it.

The risk is silent: a single-observer dataset fits without error, and a
result reported as if perception had been handled is wrong by a factor.
This function states that limit explicitly.

## Candidate covariates

A column is a candidate covariate if all of the following hold:

- It is not a structural column: `object`, `distance`, `distbegin`,
  `distend`, `observer`, `detected`, `size`, `Region.Label`, `Area`,
  `Sample.Label`, `Effort`, `detected_by`, `sample`, `distance_source`,
  or a per-sighting NARWC column (identifiers, positions, dates, raw
  angle and strip fields). Names are matched without regard to case,
  after resolving NARWC aliases.

- It is not a superseded `<NAME>_ORIGINAL` column whose stem is a known
  variable. `narwcr::read_narwc(prefer_source = TRUE)` keeps a replaced
  column under that name for traceability.

- It is numeric, character, factor or logical, and has more than one
  distinct non-missing value.

Candidates that NARWC records once per survey leg (sea state,
visibility, glare, altitude, platform) are listed as survey conditions.
The NARWC vocabulary comes from the `narwcr` package when it is
installed. Without it, only the structural names above are used.

## Duplicate columns

Two columns that name one variable, after case-folding and alias
resolution, are reported in `summary$duplicates`. An example is a
MEMDR-era NARWC file with both `LEGTYPE` and `LEGTYPE_BK`. The pair is
ambiguous: nothing here can tell which one a fit should use.
`narwcr::read_narwc(prefer_source = TRUE)` resolves it.

## Perception bias

Mark-recapture distance sampling (`mrds` methods `"io"` and `"trial"`)
needs one row per observer per object, with an `observer` column and a
`detected` column (the detection indicator). Column values are not
checked. The check gives:

- `FALSE` if neither `observer` nor `detected` is present, or `observer`
  has fewer than two values.

- `NA` if only one of the two columns is present, or if `object` is
  present but no value repeats (no animal recorded by both teams).

- `TRUE` otherwise. Without an `object` column, repeats are not checked.

## What it does not tell you

"Supported" is a statement about structure only. Enough detections to
fit a detection function is not enough to fit one well. A covariate
being present is not a reason to put it in a model.

Availability is unsupported for every table. *Availability* is the
probability that an animal is at the surface and can be seen. An animal
submerged for the whole pass is missed at every distance equally, so
distances carry no information about it. It is computed from dive data
with
[`availability()`](https://camilleross.org/dsfit/reference/availability.md).

## References

Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers,
D.L. and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford
University Press. Source of the 40 and 60-80 detection guidelines.

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
#>     fit covariate models            candidates: beaufort; recorded as survey conditions: beaufort
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
