# Check and prepare distance data for fitting

Validates the table a sweep is about to be fitted to and puts it in the
shape `mrds` wants. Called by
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md);
exported because the checks are worth running on their own before
committing to a model set.

## Usage

``` r
prepare_distance_data(data, truncation, left = NULL, breaks = NULL)
```

## Arguments

- data:

  A data frame with `distance`, or with `distbegin` and `distend`.

- truncation:

  Right truncation distance.

- left:

  Left truncation distance, or `NULL`.

- breaks:

  Bin cutpoints for interval data. Derived from the data when `NULL`.

## Value

A list with `data` (ready for `mrds`), `binned`, `breaks`, and
`n_dropped`, and `dropped` — the attrition split by reason, so a
truncation that trimmed a tail is distinguishable from one that threw
away half the survey.

## What it refuses

- Both point and interval distances:

  `STRIP`-derived distances are intervals and are fitted binned; angle-
  and position-derived distances are points. The likelihoods differ, so
  one sweep cannot rank both. A table containing both marks a survey-era
  boundary, and should be split on its provenance column and swept
  separately.

- An unbounded top bin:

  `distend` of `Inf` cannot be fitted. The open top bin of every `STRIP`
  scheme has to be dropped or closed before fitting.

- Bins that do not tile:

  Interval data whose `distbegin`/`distend` pairs leave gaps or overlap
  does not define a set of cutpoints, and any `breaks` derived from it
  would silently misallocate detections.

- Bins that stop short of the truncation:

  A binned fit integrates the detection function over the bins. If the
  top bin ends before the truncation width, the strip between them is
  unaccounted effort and the fit describes a narrower survey than the
  one flown.

## What it drops, and reports

Rows with no distance at all — which is how a flatfile records a sample
that produced no detections — and rows beyond `truncation` or inside
`left`. Dropping is counted, never silent.

## See also

[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)

## Examples

``` r
d <- data.frame(object = 1:5, distance = c(10, 50, 120, NA, 900))
prep <- prepare_distance_data(d, truncation = 400)
prep$n_dropped
#> [1] 2
prep$data
#>   object distance
#> 1      1       10
#> 2      2       50
#> 3      3      120
```
