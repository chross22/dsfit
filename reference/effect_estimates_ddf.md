# Extract a fitted detection function as a plottable data frame

Evaluates a fitted `mrds` detection function \\g(x)\\ on a grid of
distances, with a delta-method interval, and returns it in the format
that
[`fancyfx::plotEffects()`](https://camilleross.org/fancyfx/reference/plotEffects.html)
and
[`fancyfx::plotRugs()`](https://camilleross.org/fancyfx/reference/plotRugs.html)
draw. The detection function \\g(x)\\ is the probability that an animal
at perpendicular distance \\x\\ from the trackline is detected.

## Usage

``` r
# S3 method for class 'ddf'
effect_estimates(
  model,
  var = "distance",
  scale = c("auto", "link", "response"),
  interval = c("auto", "se", "ci"),
  level = 0.95,
  n = 100,
  at = NULL,
  ...
)

effect_estimates_ddf(
  model,
  var = "distance",
  scale = c("auto", "link", "response"),
  interval = c("auto", "se", "ci"),
  level = 0.95,
  n = 100,
  at = NULL,
  ...
)
```

## Arguments

- model:

  A fitted `ddf` object, from
  [`mrds::ddf()`](https://rdrr.io/pkg/mrds/man/ddf.html) or
  [`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md).

- var:

  Must be `"distance"`. Plotting detectability against a covariate is
  not supported.

- scale:

  Accepted for compatibility with the generic and ignored. The result is
  always on the probability scale.

- interval:

  `"se"` for plus and minus one standard error, `"ci"` for a confidence
  interval at `level`, or `"auto"` (the default) for `"se"`. `"cri"` is
  an error, because there is no posterior.

- level:

  Confidence level, used only when `interval = "ci"`.

- n:

  Number of distances at which to evaluate \\g(x)\\.

- at:

  Named list of covariate values to fix, or `NULL` (default) to average
  over the observed covariate values.

- ...:

  Ignored.

## Value

A data frame with one row per grid distance and columns `.x` (distance),
`.estimate` (\\g(x)\\), `.lower` and `.upper` (interval bounds). Its
`"quantity"` attribute is `"Detection probability g(x)"`, which
`fancyfx` uses as the axis label. This is the format
[`fancyfx::effect_estimates()`](https://camilleross.org/fancyfx/reference/effect_estimates.html)
specifies.

## How it works

1.  `interval` and `scale` are checked. `"auto"` becomes `"se"`, and
    `"cri"` is refused.

2.  `var` is checked: only `"distance"` is accepted.

3.  A grid of `n` distances is built from the left truncation distance
    (0 if none) to the right truncation width.

4.  If `at` is given, the named covariates are set to those values.

5.  \\g(x)\\ is evaluated at every grid distance for every observed
    covariate combination, and the curves are averaged.

6.  The standard error at each distance is computed by the delta method
    from the fitted parameters and their inverse Hessian.

7.  The lower and upper bounds are the estimate minus and plus one
    standard error (`"se"`) or a normal quantile times it (`"ci"`),
    clamped to `[0, 1]`.

## What the curve is

The curve is \\g(x)\\ conditional on the animal being available to be
seen. No `g(0)` correction is applied; see
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
and [`g0()`](https://camilleross.org/dsfit/reference/g0.md). Plotting it
over a rug of the observed distances is the standard check of a
detection function's fit.

## Covariate models

A model with detection covariates has one \\g(x)\\ per covariate
combination. This function averages those curves over the covariate
values of the detected animals. It does not evaluate one curve at the
mean covariate values, which is a different and less meaningful
quantity.

`at` fixes named covariates at chosen values. Covariates not named in
`at` are still averaged over their observed values.

## Relation to `esw` and average detection probability

The average detection probability \\p\\ is the mean of \\g(x)\\ over the
truncation interval. The effective strip half-width (`esw`) in
[`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md)
is \\p\\ times the width of that interval.

For an intercept-only model, the mean of the curve returned here equals
the fitted \\p\\. The package tests check this.

For a covariate model the two differ slightly. `mrds` reports
`average.p` as a Horvitz-Thompson mean, \\n / \sum 1/p_i\\, which
weights each animal by the inverse of its own detection probability.
This curve is the plain mean over the detected animals, \\(1/n) \sum g(x
\mid z_i)\\, which is the right curve to compare with a histogram of
their distances. On the package's covariate test fixture (half-normal,
`~bf`, width 400) the plain mean is 0.401 and `average.p` is 0.388, a
difference of 3.5%. The area under this curve therefore does not
reproduce `esw` exactly for a covariate model.

## The interval

\\g(x)\\ is a fixed function of the fitted parameters, so its
uncertainty comes from theirs. The Jacobian of \\g\\ with respect to the
parameters is computed numerically and combined with the inverse Hessian
(the delta method). Bounds are clamped to `[0, 1]`, so a bound at 0 or 1
is a limit, not an estimate. If the Hessian cannot be inverted, a
warning is given and the bounds are `NA`.

`fancyfx` sends models it does not recognise to `marginaleffects`, which
finds no predictors in a `ddf` object. This method is what lets a `ddf`
be plotted.

## References

Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers,
D.L. and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford
University Press.

## See also

[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)

## Examples

``` r
set.seed(1)
d <- data.frame(object = 1:200, distance = abs(rnorm(200, 0, 120)))
d <- d[d$distance < 400, ]
fit <- mrds::ddf(dsmodel = ~cds(key = "hn", formula = ~1), data = d,
                 method = "ds", meta.data = list(width = 400))

est <- effect_estimates_ddf(fit)
head(est)
#>          .x .estimate    .lower    .upper
#> 1  0.000000 1.0000000 1.0000000 1.0000000
#> 2  4.040404 0.9993443 0.9992752 0.9994134
#> 3  8.080808 0.9973797 0.9971037 0.9976556
#> 4 12.121212 0.9941140 0.9934951 0.9947328
#> 5 16.161616 0.9895599 0.9884647 0.9906551
#> 6 20.202020 0.9837353 0.9820342 0.9854364
attr(est, "quantity")
#> [1] "Detection probability g(x)"

# The mean of the curve is the fitted average detection probability
mean(est$.estimate) # approximately summary(fit)$average.p
#> [1] 0.3509216
```
