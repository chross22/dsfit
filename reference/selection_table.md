# Get the selection table from a sweep

Returns the table of fitted candidates from a
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
result, sorted by AIC (lowest first).

## Usage

``` r
selection_table(x, converged_only = TRUE)
```

## Arguments

- x:

  A `dsfit_sweep` object from
  [`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md).

- converged_only:

  Logical. `TRUE` (default) drops candidates that failed to fit or did
  not converge.

## Value

A tibble with one row per candidate and columns:

- `model_id`, `key`, `adjustment`, `order`, `formula`:

  The candidate, as in
  [`model_set()`](https://camilleross.org/dsfit/reference/model_set.md).

- `converged`:

  `TRUE` if the fit succeeded.

- `n_par`:

  Number of estimated parameters.

- `aic`, `delta_aic`:

  AIC, and its difference from the lowest AIC.

- `p`, `p_se`, `p_cv`:

  Average detection probability over the searched strip, its standard
  error, and its coefficient of variation.

- `esw`:

  Effective strip half-width: `p` times the width of the searched strip
  (`truncation` minus `left`).

- `cvm_p`:

  Cramer-von Mises goodness-of-fit p-value. `NA` for binned distances.

- `chisq_p`:

  Chi-square goodness-of-fit p-value. The test to read for binned
  distances.

Columns describing a fit are `NA` for a candidate that failed.

## See also

[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)

## Examples

``` r
set.seed(1)
d <- data.frame(object = 1:150, distance = abs(rnorm(150, 0, 100)))
d <- d[d$distance < 300, ]
selection_table(sweep_models(d, truncation = 300, quiet = TRUE))
#> # A tibble: 2 × 15
#>   model_id key   adjustment order formula converged n_par   aic delta_aic     p
#>   <chr>    <chr> <chr>      <int> <chr>   <lgl>     <dbl> <dbl>     <dbl> <dbl>
#> 1 hn       hn    NA            NA ~1      TRUE          1 1570.      0    0.378
#> 2 hr       hr    NA            NA ~1      TRUE          2 1574.      4.47 0.388
#> # ℹ 5 more variables: p_se <dbl>, p_cv <dbl>, esw <dbl>, cvm_p <dbl>,
#> #   chisq_p <dbl>
```
