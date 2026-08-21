# Fit a set of detection functions and compare them

Fits every candidate in a
[`model_set()`](https://camilleross.org/dsfit/reference/model_set.md) to
the same data, at the same truncation, and returns a selection table
ranked on more than AIC.

## Usage

``` r
sweep_models(
  data,
  models = NULL,
  truncation,
  left = NULL,
  breaks = NULL,
  quiet = FALSE
)
```

## Arguments

- data:

  A data frame with one row per detection. Needs `distance` for exact
  distances, or `distbegin` and `distend` for binned ones. The flatfile
  from `distsamp::detection_data()` is this shape; rows with no
  detection are dropped here.

- models:

  A
  [`model_set()`](https://camilleross.org/dsfit/reference/model_set.md),
  or `NULL` for the default set.

- truncation:

  Right truncation distance. Required: it defines the data every model
  in the sweep is fitted to.

- left:

  Left truncation distance, for a platform that cannot see beneath
  itself. `NULL` (default) for none.

- breaks:

  Bin cutpoints, for binned data. Taken from `distbegin`/ `distend` when
  not given.

- quiet:

  Suppress the progress and failure report. Default `FALSE`.

## Value

An object of class `dsfit_sweep`: a list with `table` (the selection
table), `fits` (the `ddf` objects, named by `model_id`), and `settings`.

## Why truncation is one argument and not a column

AIC compares models fitted to the same data. Changing the truncation
changes which detections are in the likelihood, so AICs either side of
that are not on the same scale and ranking them together is meaningless.

Truncation is therefore a single value for a whole sweep, by
construction. To compare truncations, run a sweep for each and compare
the *sweeps* — which is a different question, answered by looking at
whether the chosen model and its `p̄` are stable, not by reading one AIC
column.

## Read `p_cv` and `esw`, not just `delta_aic`

What propagates into abundance is the effective strip half-width. Two
models within 2 AIC of one another can imply materially different values
of it, and the difference goes straight into density. The table
therefore carries `p` (average detection probability), its standard
error and CV, and `esw` alongside the AIC, and
[`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md)
sorts on AIC only because something has to be first.

A goodness-of-fit statistic is reported too, whose null is that the
fitted function describes the observed distances. A small p-value is
evidence against the model, and a model can be top of the AIC ranking
and still fail it.

Which test to read follows the data rather than a preference. Exact
distances get Cramér-von Mises (`cvm_p`), which tests their empirical
distribution directly. Binned distances have no such distribution, and
get a chi-square over the survey's own bins (`chisq_p`) — without which
a binned sweep would carry no goodness-of-fit at all. `chisq_p` is
filled in for exact fits too, but over cutpoints `mrds` chooses; a test
whose result moves with an arbitrary binning is the weaker one, and
`cvm_p` is what to read there. Both columns are named for their test, so
a p-value is never read without knowing where it came from, and
[`print()`](https://rdrr.io/r/base/print.html) shows the one that
applies.

A chi-square over few bins runs out of degrees of freedom quickly —
three bins and a two-parameter key leave none — and `mrds` returns `NA`
rather than a p-value when it does. That is a converged model with no
goodness-of-fit test available, not a failed one.

## Binned and exact distances cannot be swept together

`STRIP`-derived distances are intervals and are fitted binned; angle-
and position-derived distances are points. The likelihoods differ, so
their AICs are not comparable, and a table containing both would be a
survey-era boundary rather than a model set. Data carrying both is
refused.

## Left truncation and the gamma key are alternatives

Both describe an aircraft that cannot see beneath itself. Applying
`left` to a model set that contains the gamma key accounts for the blind
spot twice, and warns.

## This does not estimate g(0)

Every fit here conditions on the animal having been available and seen.
For a deep-diving species that assumption biases density low, and no
amount of model selection detects it — a mis-specified `g(0)` shifts
every candidate in the set by the same factor, so the ranking looks
untouched. Correct for it at the abundance step, from external sources,
with its standard error propagated.

## References

Miller, D.L., Rexstad, E., Thomas, L., Marshall, L. and Laake, J.L.
(2019) Distance sampling in R. *Journal of Statistical Software*
89(1):1-28.
[doi:10.18637/jss.v089.i01](https://doi.org/10.18637/jss.v089.i01)

Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers,
D.L. and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford
University Press.

## See also

[`model_set()`](https://camilleross.org/dsfit/reference/model_set.md),
[`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md)

## Examples

``` r
set.seed(1)
d <- data.frame(object = 1:200, distance = abs(rnorm(200, 0, 120)))
d <- d[d$distance < 400, ]

sw <- sweep_models(d, model_set(c("hn", "hr", "gamma")), truncation = 400)
#> Fitted 3 of 3 candidate models to 200 detections, truncation 400.
#>   g(0) = 1 is assumed; see `?sweep_models`.
sw
#> <dsfit_sweep>
#>   detections:  200
#>   truncation:  400
#>   models:      3 of 3 fitted
#> 
#>  model  dAIC      p  p_cv   esw CvM_p
#>     hn  0.00 0.3494 0.052 139.8 0.928
#>  gamma  9.92 0.3010 0.074 120.4 0.498
#>     hr 10.18 0.3781 0.064 151.3 0.786
#> 
#>   Rank on esw and p_cv as well as dAIC: models within 2 AIC can
#>   imply materially different abundance. g(0) = 1 is assumed.

# AIC alone would not have told you this
selection_table(sw)[, c("model_id", "delta_aic", "p", "p_cv", "esw")]
#> # A tibble: 3 × 5
#>   model_id delta_aic     p   p_cv   esw
#>   <chr>        <dbl> <dbl>  <dbl> <dbl>
#> 1 hn            0    0.349 0.0525  140.
#> 2 gamma         9.92 0.301 0.0737  120.
#> 3 hr           10.2  0.378 0.0638  151.
```
