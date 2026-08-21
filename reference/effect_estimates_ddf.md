# Extract a fitted detection function as a tidy frame

An
[`fancyfx::effect_estimates()`](https://rdrr.io/pkg/fancyfx/man/effect_estimates.html)
method for `mrds` detection functions, so that a fitted `ddf` can be
plotted with
[`fancyfx::plotEffects()`](https://rdrr.io/pkg/fancyfx/man/plotEffects.html)
and
[`fancyfx::plotRugs()`](https://rdrr.io/pkg/fancyfx/man/plotRugs.html)
like any other model.

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

  Must be `"distance"`. A detection function has one predictor; plotting
  detectability against a covariate is a different quantity and is not
  built yet.

- scale:

  Accepted for compatibility with the generic, and ignored. A detection
  function has no separable linear predictor, so every setting gives
  \\g(x)\\ on the probability scale.

- interval:

  `"se"` for a one-standard-error ribbon, `"ci"` for a confidence
  interval at `level`, or `"auto"` (which `fancyfx` forwards by default)
  for `"se"`. `"cri"` is refused: there is no posterior here to
  summarise.

- level:

  Confidence level, used when `interval = "ci"`.

- n:

  Number of distances at which to evaluate the function.

- at:

  Named list fixing covariate values, or `NULL` (default) to average
  over the covariates as observed.

- ...:

  Ignored.

## Value

A data frame with `.x`, `.estimate`, `.lower`, `.upper`, and a
`"quantity"` attribute, as
[`fancyfx::effect_estimates()`](https://rdrr.io/pkg/fancyfx/man/effect_estimates.html)
specifies.

## Why this method has to exist

`fancyfx` sends anything that is not a GAM to `marginaleffects`, which
cannot introspect a `ddf` — it reports no predictors at all. The generic
is designed to be extended from outside the package, and this is that
extension. Nothing downstream needs to change: the transforms, the
ribbon, the rug, and the panel arranging all work off the frame returned
here.

## What it returns, and what it does not

The detection function \\g(x)\\: the probability that an animal at
perpendicular distance \\x\\ was detected, given that it was available
to be. It is **not** the probability that an animal at that distance was
there and seen, and no `g(0)` correction is applied — see
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md).

A rug of the observed distances beneath this curve is the standard
detection-function diagnostic, which is `fancyfx`'s premise applied
almost exactly.

## Covariate models

the curve is averaged, not evaluated at a mean: With detection
covariates there is one \\g(x)\\ per covariate combination. This
averages those curves over the covariate values actually observed,
rather than evaluating the function once at mean covariates. The two are
different, and averaging the curves is the one that has a meaning: it is
the mean detection function of the animals in the sample.

`at` fixes named covariates at chosen values; anything not named is
still averaged over as observed. `at = NULL` averages over everything.

## Why the area under this curve is not `esw`

For an intercept-only model, the mean of this curve over the truncation
width is exactly the fitted average detection probability, and the area
under it is the effective strip half-width in
[`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md).
The tests assert that.

For a **covariate** model the two differ slightly, and it is worth
knowing why rather than discovering it. `mrds` reports `average.p` as a
Horvitz-Thompson mean, \\n / \sum 1/p_i\\, which weights each animal by
the inverse of its own detection probability — correcting for the fact
that a sample of detections over-represents the covariate values that
are easy to detect. The curve here is the plain arithmetic mean over the
animals as observed, \\(1/n) \sum g(x \mid z_i)\\, which is what you
want to look at next to a histogram of those same animals' distances.

The two coincide when detection probability is constant, and differ by a
fraction of a percent otherwise. `esw` in the selection table uses the
Horvitz-Thompson version, so eyeballing the area under this curve will
not reproduce it exactly for a covariate model.

## The interval is a delta-method approximation

\\g(x)\\ is deterministic given the fitted parameters, so the
uncertainty comes from them: the Jacobian of \\g\\ with respect to the
parameter vector is taken numerically and combined with the inverse
Hessian. The result is clamped to `[0, 1]`, since a probability cannot
leave it — which means an interval touching 0 or 1 is a bound, not an
estimate. Where the Hessian cannot be inverted the estimates are
returned with missing bounds rather than no estimates.

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

# The integral of the curve is the fitted average detection probability
mean(est$.estimate) # approximately summary(fit)$average.p
#> [1] 0.3509216
```
