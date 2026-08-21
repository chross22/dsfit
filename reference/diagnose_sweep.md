# Diagnose common reasons a sweep fails, or succeeds meaninglessly

Runs the input guards, the truncation, and the model-set expansion —
never [`mrds::ddf()`](https://rdrr.io/pkg/mrds/man/ddf.html) itself —
and reports the most common ways a sweep goes wrong before it reaches
the fitting: a truncation that throws away most of the survey, a
covariate formula naming a column that is not there, a model set that
double-counts the blind spot, or too few detections left to fit anything
worth reading.

## Usage

``` r
diagnose_sweep(data, models = NULL, truncation, left = NULL, breaks = NULL)
```

## Arguments

- data:

  A data frame with one row per detection, as
  [`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
  takes.

- models:

  A
  [`model_set()`](https://camilleross.org/dsfit/reference/model_set.md),
  or `NULL` for the default set.

- truncation:

  Right truncation distance. Required, as it is for a sweep.

- left:

  Left truncation distance, or `NULL`.

- breaks:

  Bin cutpoints for interval data, or `NULL` to derive them.

## Value

Invisibly, `list(structure, prepared, models)` — whichever were reached
before a fatal problem stopped the checks, so investigation can pick up
from there.

## Details

Meant to run *before*
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md),
so a misconfiguration is caught in seconds rather than after a sweep
that either errors from inside `mrds` or returns a table that looks fine
and is not.

Every check here reports a problem rather than fixing it. This function
never modifies the data, the model set, or anything else.

## See also

[`detection_structure()`](https://camilleross.org/dsfit/reference/detection_structure.md),
which this uses and which answers the different question of what the
data could support in principle.
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md),
which this is meant to run ahead of.

## Examples

``` r
set.seed(1)
d <- data.frame(
  object = 1:300,
  distance = c(abs(rnorm(295, 0, 800)), rep(NA, 5)),
  beaufort = sample(0:4, 300, replace = TRUE)
)
diagnose_sweep(d, model_set(c("hn", "hr")), truncation = 2000)
#> dsfit sweep diagnosis
#> 
#> == Toolchain ==
#>   ok    mrds 3.0.1 is installed
#> 
#> == Data ==
#>   ok    300 rows, 295 exact distances
#>   note  estimate perception bias: no double-observer structure: no `observer` or `detected` column
#>   note  carry group size to abundance: no `size` column; abundance would count groups, not individuals
#> 
#> == Truncation ==
#>   ok    293 of 300 rows kept at truncation 2000
#>         dropped: 5 with no distance, 2 beyond the truncation
#>   ok    293 detections is at or above the 60-80 usually suggested
#> 
#> == Model set ==
#>   ok    2 candidates: hn, hr
#> 
#> == What is still assumed ==
#>   g(0) = 1, unless a correction is applied at the abundance step.
#>   Nothing here can detect a wrong one: it scales every candidate
#>   equally, so the ranking looks untouched. See ?g0.
#> 
#> No problems found. Nothing was fitted.
```
