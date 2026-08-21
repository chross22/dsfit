# The ranked selection table from a sweep

The ranked selection table from a sweep

## Usage

``` r
selection_table(x, converged_only = TRUE)
```

## Arguments

- x:

  A `dsfit_sweep` from
  [`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md).

- converged_only:

  Drop models that failed to fit. Default `TRUE`.

## Value

A tibble, sorted by AIC, with `model_id`, `key`, `adjustment`, `order`,
`formula`, `converged`, `n_par`, `aic`, `delta_aic`, `p`, `p_se`,
`p_cv`, `esw`, and a goodness-of-fit p-value in `cvm_p` (exact
distances) or `chisq_p` (binned).

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
