# Build a set of candidate detection functions

Expands key functions, adjustment series and covariate formulas into a
grid of candidate models, one row per candidate. Pass the result to
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
to fit them all.

## Usage

``` r
model_set(key = c("hn", "hr"), adjustment = NULL, order = NULL, formula = ~1)
```

## Arguments

- key:

  Character vector of key functions: any of `"hn"`, `"hr"`, `"unif"`,
  `"gamma"`. Default `c("hn", "hr")`.

- adjustment:

  Character vector of adjustment series: `"cos"` (cosine), `"herm"`
  (Hermite polynomial) or `"poly"` (simple polynomial). One candidate
  per series. `NULL` (default) fits no adjustment terms.

- order:

  Integer vector of highest adjustment orders, one candidate per
  element. Orders are cumulative, as in `Distance::ds()`: a candidate
  carries every term of its series up to `order`. For example,
  `order = 3` on `"hn"` with `"cos"` fits `adj.order = c(2, 3)`. Cosine
  orders run 2, 3, 4, ... (1, 2, 3, ... for `"unif"`). Hermite and
  polynomial orders run 4, 6, 8, ... (2, 4, 6, ... for `"unif"`). `NULL`
  (default) gives the first term of the series only. An order outside
  its series is dropped with a warning. Ignored when `adjustment` is
  `NULL`.

- formula:

  One-sided formula, or a list of them, giving detection covariates.
  `~1` (default) means no covariates and fits with `cds`; any other
  formula fits with `mcds`.

## Value

A tibble with one row per candidate and columns:

- `model_id`:

  Label, e.g. `"hn"`, `"hr+cos(2,3)"`, `"hn wt_beaufort"`.

- `key`:

  Key function.

- `adjustment`:

  Adjustment series, or `NA`.

- `order`:

  Highest adjustment order, or `NA`.

- `formula`:

  Covariate formula as a string.

- `dsmodel`:

  The model formula string passed to
  [`mrds::ddf()`](https://rdrr.io/pkg/mrds/man/ddf.html).

## How it works

1.  Checks that every key, adjustment series, order and formula is
    valid.

2.  Crosses keys, series, orders and formulas into a grid.

3.  Drops, with a warning, candidates that cannot be fitted
    meaningfully: a uniform key with no adjustment, a uniform key with
    covariates, and an adjustment order outside its series.

4.  Builds each candidate's `mrds` model string and its `model_id`
    label.

5.  Drops, with a warning, duplicate candidates.

## Key functions

- `"hn"`:

  Half-normal. A shoulder at zero and a smooth fall-off. The usual
  default.

- `"hr"`:

  Hazard-rate. A broader shoulder and a heavier tail, with a second
  parameter.

- `"unif"`:

  Uniform. Detection is certain out to the truncation distance, which is
  a strip transect. It is kept only with adjustment terms. It takes no
  covariates, because it has no scale parameter for them to act on.

- `"gamma"`:

  Gamma. Unimodal, with its peak away from zero, so detection is lowest
  on the track-line. `mrds` scales it to 1 at its peak, not at zero. Not
  available in `Distance::ds()`.

## Gamma key and left truncation

An aircraft cannot see the water directly beneath it, and a ship cannot
see past its bow wave. For such a platform the gamma key fits the near
distances better than a half-normal, which must peak at zero. Left
truncation models the same blind spot by discarding the near distances.
Use one or the other, not both.
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
takes `left` and warns if it is set for a model set that contains the
gamma key.

## References

Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers,
D.L. and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford
University Press.

Laake, J.L. and Borchers, D.L. (2004) Methods for incomplete detection
at distance zero. In *Advanced Distance Sampling*, Oxford University
Press.

Becker, E.A. and Quang, P.X. (2009) A gamma-shaped detection function
for line-transect surveys with mark-recapture and covariate data.
*Journal of Agricultural, Biological and Environmental Statistics*
14:207-223.
[doi:10.1198/jabes.2009.0013](https://doi.org/10.1198/jabes.2009.0013)
Source of the gamma key, developed for aerial surveys.

## See also

[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)

## Examples

``` r
# Default: half-normal and hazard-rate, no adjustments or covariates
model_set()
#> # A tibble: 2 × 6
#>   model_id key   adjustment order formula dsmodel                           
#>   <chr>    <chr> <chr>      <int> <chr>   <chr>                             
#> 1 hn       hn    NA            NA ~1      "~cds(key = \"hn\", formula = ~1)"
#> 2 hr       hr    NA            NA ~1      "~cds(key = \"hr\", formula = ~1)"

# Cosine adjustments up to orders 2 and 3
model_set(key = c("hn", "hr"), adjustment = "cos", order = 2:3)
#> # A tibble: 4 × 6
#>   model_id    key   adjustment order formula dsmodel                            
#>   <chr>       <chr> <chr>      <int> <chr>   <chr>                              
#> 1 hn+cos(2)   hn    cos            2 ~1      "~cds(key = \"hn\", formula = ~1, …
#> 2 hr+cos(2)   hr    cos            2 ~1      "~cds(key = \"hr\", formula = ~1, …
#> 3 hn+cos(2,3) hn    cos            3 ~1      "~cds(key = \"hn\", formula = ~1, …
#> 4 hr+cos(2,3) hr    cos            3 ~1      "~cds(key = \"hr\", formula = ~1, …

# With and without a detection covariate
model_set(key = "hn", formula = list(~1, ~wt_beaufort))
#> # A tibble: 2 × 6
#>   model_id       key   adjustment order formula      dsmodel                    
#>   <chr>          <chr> <chr>      <int> <chr>        <chr>                      
#> 1 hn             hn    NA            NA ~1           "~cds(key = \"hn\", formul…
#> 2 hn wt_beaufort hn    NA            NA ~wt_beaufort "~mcds(key = \"hn\", formu…

# The gamma key, for a platform that cannot see beneath itself
model_set(key = c("hn", "gamma"))
#> # A tibble: 2 × 6
#>   model_id key   adjustment order formula dsmodel                              
#>   <chr>    <chr> <chr>      <int> <chr>   <chr>                                
#> 1 hn       hn    NA            NA ~1      "~cds(key = \"hn\", formula = ~1)"   
#> 2 gamma    gamma NA            NA ~1      "~cds(key = \"gamma\", formula = ~1)"
```
