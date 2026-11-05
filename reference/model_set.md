# Build a set of candidate detection functions

Expands key functions, adjustment terms, and covariate formulas into a
grid of candidate models, one row each, ready for
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md).

## Usage

``` r
model_set(key = c("hn", "hr"), adjustment = NULL, order = NULL, formula = ~1)
```

## Arguments

- key:

  Character vector of key functions: any of `"hn"`, `"hr"`, `"unif"`,
  `"gamma"`.

- adjustment:

  Adjustment series: `NULL` (none), `"cos"`, `"herm"`, or `"poly"`. A
  vector expands to one candidate per series.

- order:

  Integer vector of highest adjustment orders, one candidate per
  element. Orders are cumulative, as in `Distance::ds()`: each candidate
  carries every term of its series up to `order`, so `order = 3` on `hn`
  with `"cos"` fits `adj.order = c(2, 3)`. The series are cosine 2, 3,
  4, ... (1, 2, 3, ... for `unif`) and Hermite or simple polynomial 4,
  6, 8, ... (2, 4, 6, ... for `unif`). `NULL` (default) gives one term,
  the first of the series. An order outside a key's series is dropped
  with a warning. Ignored when `adjustment` is `NULL`.

- formula:

  One-sided formula, or a list of them, giving detection covariates.
  `~1` (default) is no covariates. Anything else fits through `mcds`
  rather than `cds`.

## Value

A tibble with one row per candidate: `model_id`, `key`, `adjustment`,
`order`, `formula`, and the `dsmodel` string `mrds` is given.

## The key functions

- `"hn"`:

  Half-normal. A shoulder at zero and a smooth fall-off; the usual
  default.

- `"hr"`:

  Hazard-rate. A broader shoulder and a heavier tail, at the cost of a
  second parameter.

- `"unif"`:

  Uniform. Detection certain out to the truncation distance — a strip
  transect. Only interesting with adjustment terms, and takes no
  covariates, having no scale parameter for them to act on. Both cases
  are dropped with a warning.

- `"gamma"`:

  Gamma. **Unimodal, with its peak away from zero**, so detection is
  lowest on the track-line. Not available in `Distance::ds()`, which is
  one reason this package fits through `mrds` directly.

## When gamma is the right shape, and when it is a trap

An aircraft cannot see the water directly beneath it, and neither can a
ship see past its own bow wave. Where that is true the gamma key
describes the survey rather than distorting it, and a half-normal forced
through a peak at zero will misfit the near distances.

But a gamma key and a **left truncation model the same phenomenon**.
Using both accounts for the blind spot twice. Pick one deliberately and
record which —
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
takes `left` for the other route, and will warn if a model set
containing gamma is fitted with one.

Note also that gamma does not assume `g(0) = 1`, so its `p̄` is not
comparable to a half-normal's in the way a shared assumption would make
it. That is a reason to read the whole selection table rather than its
first row.

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
The gamma key function, and the aerial-survey case it was developed for.

## See also

[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)

## Examples

``` r
model_set()
#> # A tibble: 2 × 6
#>   model_id key   adjustment order formula dsmodel                           
#>   <chr>    <chr> <chr>      <int> <chr>   <chr>                             
#> 1 hn       hn    NA            NA ~1      "~cds(key = \"hn\", formula = ~1)"
#> 2 hr       hr    NA            NA ~1      "~cds(key = \"hr\", formula = ~1)"

# With adjustments
model_set(key = c("hn", "hr"), adjustment = "cos", order = 2:3)
#> # A tibble: 4 × 6
#>   model_id    key   adjustment order formula dsmodel                            
#>   <chr>       <chr> <chr>      <int> <chr>   <chr>                              
#> 1 hn+cos(2)   hn    cos            2 ~1      "~cds(key = \"hn\", formula = ~1, …
#> 2 hr+cos(2)   hr    cos            2 ~1      "~cds(key = \"hr\", formula = ~1, …
#> 3 hn+cos(2,3) hn    cos            3 ~1      "~cds(key = \"hn\", formula = ~1, …
#> 4 hr+cos(2,3) hr    cos            3 ~1      "~cds(key = \"hr\", formula = ~1, …

# With a detection covariate
model_set(key = "hn", formula = list(~1, ~wt_beaufort))
#> # A tibble: 2 × 6
#>   model_id       key   adjustment order formula      dsmodel                    
#>   <chr>          <chr> <chr>      <int> <chr>        <chr>                      
#> 1 hn             hn    NA            NA ~1           "~cds(key = \"hn\", formul…
#> 2 hn wt_beaufort hn    NA            NA ~wt_beaufort "~mcds(key = \"hn\", formu…

# The unimodal key, for a platform with a blind spot beneath it
model_set(key = c("hn", "gamma"))
#> # A tibble: 2 × 6
#>   model_id key   adjustment order formula dsmodel                              
#>   <chr>    <chr> <chr>      <int> <chr>   <chr>                                
#> 1 hn       hn    NA            NA ~1      "~cds(key = \"hn\", formula = ~1)"   
#> 2 gamma    gamma NA            NA ~1      "~cds(key = \"gamma\", formula = ~1)"
```
