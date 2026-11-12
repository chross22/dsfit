# Fit a set of detection functions and compare them

Fits every candidate in a
[`model_set()`](https://camilleross.org/dsfit/reference/model_set.md) to
the same data at the same truncation, using
[`mrds::ddf()`](https://rdrr.io/pkg/mrds/man/ddf.html), and returns a
selection table. Use it to choose a detection function for a survey.

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

  Data frame with one row per detection. Needs a `distance` column for
  exact distances, or `distbegin` and `distend` for binned ones. The
  flatfile from `distsamp::detection_data()` has this shape. Rows with
  no distance (samples without a detection) are dropped.

- models:

  A
  [`model_set()`](https://camilleross.org/dsfit/reference/model_set.md).
  `NULL` (default) uses
  [`model_set()`](https://camilleross.org/dsfit/reference/model_set.md).

- truncation:

  Number. Right truncation distance, in the units of the distances.
  Required, with no default.

- left:

  Number. Left truncation distance, for a platform that cannot see
  beneath itself. `NULL` (default) for none.

- breaks:

  Numeric vector of bin cutpoints, for binned data. `NULL` (default)
  derives them from `distbegin` and `distend`. Ignored, with a warning,
  for exact distances.

- quiet:

  Logical. `TRUE` suppresses the summary message. Default `FALSE`.

## Value

An object of class `dsfit_sweep`, a list with:

- `table`:

  The selection table, one row per candidate including failures. See
  [`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md)
  for the columns.

- `fits`:

  List of `ddf` objects named by `model_id`. `NULL` for a candidate that
  failed.

- `settings`:

  List with `truncation`, `left`, `binned`, `breaks`, `n` (detections
  fitted), `n_dropped` (rows removed during preparation) and `failures`
  (first line of the error, or `"did not converge"`, named by
  `model_id`).

## How it works

1.  Checks that `mrds` is installed, that `models` is a model set, and
    that `truncation` is a single positive number.

2.  Prepares the data with
    [`prepare_distance_data()`](https://camilleross.org/dsfit/reference/prepare_distance_data.md):
    validates the distances, drops rows with no distance or outside the
    truncation, and decides whether the fit is binned.

3.  Warns if `left` is set and the model set contains the gamma key.

4.  Fits each candidate with
    [`mrds::ddf()`](https://rdrr.io/pkg/mrds/man/ddf.html). A candidate
    that errors or does not converge is recorded as a failure and the
    sweep continues.

5.  Builds the selection table, sorted by AIC, and prints a short report
    unless `quiet = TRUE`.

## Truncation is one value per sweep

AIC only compares models fitted to the same data. A different truncation
keeps a different set of detections, so the AICs are not comparable. To
compare truncations, run one sweep per truncation and check whether the
chosen model and its `p` stay stable.

## What the table reports besides AIC

The effective strip half-width (esw) is the width of a strip in which
detection is certain that would yield the same expected number of
detections. It goes directly into density. Two models within 2 AIC of
each other can have quite different esw, so the table reports `p`
(average detection probability), its standard error and CV, and `esw`
next to AIC.

A goodness-of-fit p-value is also reported. The null hypothesis is that
the fitted function describes the observed distances, so a small p-value
is evidence against the model. The top AIC model can still fail this
test.

## Which goodness-of-fit test to read

Exact distances get a Cramer-von Mises test (`cvm_p`), which tests their
empirical distribution. Binned distances have no empirical distribution,
so they get a chi-square test over the survey's own bins (`chisq_p`).
`chisq_p` is also filled in for exact distances, but over cutpoints
`mrds` chooses; read `cvm_p` there. Printing a sweep shows the test that
applies.

A chi-square test over few bins can have no degrees of freedom left (for
example three bins and a two-parameter key). `mrds` then returns `NA`.
The model converged; it just has no chi-square test.

## Binned and exact distances cannot be mixed

Distances from `STRIP` records are intervals and are fitted binned.
Distances from angles or positions are exact. The two likelihoods
differ, so their AICs are not comparable. Data containing both is
refused; split it and run one sweep for each.

## Left truncation and the gamma key

Both model a platform that cannot see directly beneath itself. Using
both accounts for the blind spot twice, so setting `left` for a model
set that contains the gamma key gives a warning.

## Monotonicity constraint

As in `Distance::ds()`, a candidate with adjustment terms and no
covariates is fitted with a strict monotonicity constraint
(`mono = TRUE, mono.strict = TRUE` in `mrds`). This stops g(x), the
detection probability at distance x, from exceeding 1 or turning upward.
AIC and `p` then match `ds()` with its defaults. Covariate models are
unconstrained, as in `ds()`. The gamma key is also unconstrained,
because it rises to a peak by design.

## g(0) is assumed to be 1

Every fit assumes an animal on the track-line is always detected. A
wrong g(0) scales every candidate by the same factor, so model selection
cannot detect it. Correct for it at the abundance step, from external
sources, and propagate its standard error (see
[`g0()`](https://camilleross.org/dsfit/reference/g0.md)).

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

# Compare esw and p_cv as well as delta AIC
selection_table(sw)[, c("model_id", "delta_aic", "p", "p_cv", "esw")]
#> # A tibble: 3 × 5
#>   model_id delta_aic     p   p_cv   esw
#>   <chr>        <dbl> <dbl>  <dbl> <dbl>
#> 1 hn            0    0.349 0.0525  140.
#> 2 gamma         9.92 0.301 0.0737  120.
#> 3 hr           10.2  0.378 0.0638  151.
```
