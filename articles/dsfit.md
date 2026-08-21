# Fitting and comparing detection functions

This walks one survey from a table of distances to a `g(0)` correction,
stopping where the package stops. The through-line is that **a detection
function is chosen for its consequences, not its AIC**, and that most of
the ways this goes wrong are silent.

``` r

library(dsfit)
```

## A survey to work with

No detection data ships with the package, so here is a simulated aerial
survey. It has a blind spot beneath the aircraft, which is what makes it
worth using: that single feature drives most of the decisions below.

``` r

set.seed(20260810)

# True perpendicular distances of animals present, out to 3 km
present <- abs(rnorm(4000, 0, 900))
present <- present[present < 3000]

# Sea state, which makes animals harder to see
beaufort <- sample(0:4, length(present), replace = TRUE)

# Detection falls off with distance, falls off further in a rough sea, and is
# zero in the 100 m the aircraft cannot see beneath itself
sigma <- 800 - 60 * beaufort
p <- exp(-present^2 / (2 * sigma^2))
p[present < 100] <- 0

seen <- runif(length(present)) < p

detections <- data.frame(
  object   = seq_len(sum(seen)),
  distance = present[seen],
  beaufort = beaufort[seen],
  size     = sample(1:3, sum(seen), replace = TRUE, prob = c(.6, .3, .1))
)

nrow(detections)
#> [1] 1999
```

## Before fitting: what does this data support?

[`detection_structure()`](https://camilleross.org/dsfit/reference/detection_structure.md)
reads the table and reports which analyses it admits. Nothing is fitted.

``` r

detection_structure(detections)
#> <dsfit_structure>
#>   1999 rows, 1999 exact distances
#>   nearest detection at 100.1 - a blind spot beneath the platform would show here
#> 
#>   can:
#>     fit a detection function        1999 exact distances; at or above the 60-80 usually suggested
#>     test fit with Cramer-von Mises  exact distances have an empirical distribution to test
#>     test fit with chi-square        over cutpoints mrds chooses, which makes it the weaker test here
#>     fit covariate models            candidates: beaufort
#>     carry group size to abundance   `size` present
#> 
#>   cannot:
#>     estimate perception bias  no double-observer structure: no `observer` or `detected` column. Perception needs two independent teams and cannot be recovered from a single-observer survey at any sample size
#>     estimate availability     not estimable from distances by construction: a submerged animal is missed at every distance equally. Compute it from dive data with availability()
#> 
#>   Structure only. That something is supported is not a reason to do it.
```

Two lines there matter more than the rest.

**The nearest detection is not at zero.** That is the blind spot, and it
decides the model set below.

**Perception bias is not estimable.** This survey has one observer team,
so no amount of data recovers it. That is a property of the survey
programme rather than of any archive the data ends up in, and it does
not improve with sample size. The point of asking now is that fitting
anyway produces a number, and the number is wrong by a factor.

## The model set

An aircraft that cannot see beneath itself has two honest treatments,
and they are alternatives rather than a pair:

- a **gamma** key, which is unimodal with its peak away from zero, or
- **left truncation**, which removes the blind spot from the data.

Using both counts the blind spot twice, and
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
warns when a model set containing gamma is fitted with `left`.

``` r

models <- model_set(key = c("hn", "hr", "gamma"))
models[, c("model_id", "key", "formula")]
#> # A tibble: 3 × 3
#>   model_id key   formula
#>   <chr>    <chr> <chr>  
#> 1 hn       hn    ~1     
#> 2 hr       hr    ~1     
#> 3 gamma    gamma ~1
```

## The sweep

Truncation is one value for the whole sweep, by construction. AIC
compares models fitted to the same data, so a table mixing truncations
is not a ranking.

``` r

sw <- sweep_models(detections, models, truncation = 2500)
#> Fitted 3 of 3 candidate models to 1999 detections, truncation 2500.
#>   g(0) = 1 is assumed; see `?sweep_models`.
```

``` r

sw
#> <dsfit_sweep>
#>   detections:  1999
#>   truncation:  2500
#>   models:      3 of 3 fitted
#> 
#>  model   dAIC      p  p_cv   esw CvM_p
#>  gamma   0.00 0.2546 0.016 636.4 0.004
#>     hr 377.36 0.3549 0.015 887.2 0.000
#>     hn 407.33 0.2957 0.018 739.3 0.001
#> 
#>   Rank on esw and p_cv as well as dAIC: models within 2 AIC can
#>   imply materially different abundance. g(0) = 1 is assumed.
```

## Read the table across, not down

The AIC column answers “which model best describes these distances”. It
does not answer “what does this imply about how many animals there
were”, and the two can disagree.

``` r

selection_table(sw)[, c("model_id", "delta_aic", "p", "p_cv", "esw", "cvm_p")]
#> # A tibble: 3 × 6
#>   model_id delta_aic     p   p_cv   esw      cvm_p
#>   <chr>        <dbl> <dbl>  <dbl> <dbl>      <dbl>
#> 1 gamma           0  0.255 0.0164  636. 0.00371   
#> 2 hr            377. 0.355 0.0154  887. 0.00000308
#> 3 hn            407. 0.296 0.0177  739. 0.00132
```

`esw` is the effective strip half-width — the width of the strip that
would have produced this many detections under certain detection. It is
what abundance is divided by, so a spread in `esw` is a spread in
density, directly and proportionally. Across these three it runs from
roughly 640 to 890 m, so the same survey implies densities differing by
about 40% depending on which model you take.

Now read `cvm_p`, whose null is that the fitted function describes the
observed distances. **Every model in this table fails it**, the winner
included. The AIC column ranked them against each other and had no way
to say that none of them fit.

For binned data the column to read is `chisq_p` instead, because
Cramér-von Mises tests an empirical distribution that binned distances
do not have.

### What the failure is telling you

It is not that a fourth key would fit better. It is that a **hard
geometric cutoff is not a detection function**. The gamma key models
detection *falling off* toward the trackline, smoothly; the aircraft
here sees nothing at all below 100 m. No key function reproduces a
cliff, so all of them misfit the near distances and the goodness-of-fit
test says so.

The fix is not a better model. It is to stop asking the model to
describe geometry.

## Left truncation, the other route

Remove the blind spot from the data rather than modelling it. Gamma is
left out of this set deliberately: it and `left` are alternative
treatments of the same phenomenon, and applying both counts the blind
spot twice —
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
warns if you try.

``` r

sw_left <- sweep_models(detections, model_set(c("hn", "hr")),
                        truncation = 2500, left = 100, quiet = TRUE)
selection_table(sw_left)[, c("model_id", "delta_aic", "p", "esw", "cvm_p")]
#> # A tibble: 2 × 5
#>   model_id delta_aic     p   esw    cvm_p
#>   <chr>        <dbl> <dbl> <dbl>    <dbl>
#> 1 hn             0   0.244  609. 0.697   
#> 2 hr            71.8 0.295  737. 0.000352
```

The half-normal now passes comfortably. Same survey, same animals, and
the difference between a model set that fits and one that does not was a
statement about the aircraft rather than about the whales.

Note that `hr` still fails, so the goodness-of-fit column is doing work
here beyond rubber-stamping the winner.

Comparing these AICs with the ones in the previous table would be a
mistake — different data, different likelihood. What *is* comparable is
what the two routes imply for `esw`, and that is a question about the
survey rather than about the models.

## Covariates

A covariate model fits a separate detection function per covariate
value, so `p̄` becomes a mean over the animals seen.

``` r

sw_cov <- sweep_models(
  detections,
  model_set(key = "hn", formula = list(~1, ~beaufort)),
  truncation = 2500, left = 100, quiet = TRUE
)
selection_table(sw_cov)[, c("model_id", "delta_aic", "p", "p_cv", "esw")]
#> # A tibble: 2 × 5
#>   model_id    delta_aic     p   p_cv   esw
#>   <chr>           <dbl> <dbl>  <dbl> <dbl>
#> 1 hn beaufort       0   0.241 0.0179  601.
#> 2 hn               26.2 0.244 0.0177  609.
```

Sea state was built into the simulation, so a model that knows about it
should win here. On real data it often does not, and a covariate that
does not earn its parameter is worth dropping.

## Plotting

[`effect_estimates.ddf()`](https://camilleross.org/dsfit/reference/effect_estimates_ddf.md)
makes a fitted detection function plottable with
[fancyfx](https://github.com/chross22/fancyfx), giving the standard
diagnostic — fitted g(x) with a delta-method ribbon and a rug of the
observed distances.

``` r

fancyfx::plotEffects(sw$fits[["hn"]], dat = detections, var = "distance")
```

For a **covariate** model, the area under that curve is not the `esw` in
the table. `mrds` reports `average.p` as a Horvitz-Thompson mean,
weighting each animal by the inverse of its own detection probability;
the plotted curve is the plain mean over the animals as observed. They
coincide only when detection probability is constant.

## g(0): what the fit assumed

Every fit above conditions on the animal having been available and seen.
A mis-specified `g(0)` scales every candidate in the set by the same
factor, so the ranking looks untouched while the density is wrong —
which is exactly why it cannot be caught by model selection and has to
be handled deliberately.

Three separate things present as `g(0) < 1`:

| component | what it is | where it comes from |
|----|----|----|
| blind spot | the aircraft cannot see beneath itself | handled above, by gamma or `left` |
| availability | the animal was submerged | dive data — external |
| perception | the animal was up and missed | a second observer team |

### Availability is computed, not estimated

An animal submerged for the whole pass is missed at *every* distance
equally, so availability leaves no signature in the distance
distribution. There is nothing to estimate it from. It is computed
instead, from dive intervals and how long the platform keeps a point in
view.

``` r

ganley_availability
#> # A tibble: 4 × 8
#>   month percent_surface_time mean_dive mean_surface availability
#>   <fct>                <dbl>     <dbl>        <dbl>        <dbl>
#> 1 Jan                     16       533           48         0.27
#> 2 Feb                     34       256          226         0.52
#> 3 Mar                     31       219           67         0.52
#> 4 Apr                     55        88          801         0.91
#> # ℹ 3 more variables: availability_variance <dbl>, n_follows <int>,
#> #   hours_followed <dbl>
```

Those are measurements from Cape Cod Bay, and the swing is the point:
**0.27 in January against 0.91 in April**, as the copepods the whales
feed on move up through the water column. A single representative figure
would be a real number for one month and wrong by a factor of three for
another.

[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
computes it from your own platform and dive data. An aircraft’s field of
view is a wedge running forward and aft, so time in view *grows* with
perpendicular distance —
[`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md)
is that geometry, and
[`view_window()`](https://camilleross.org/dsfit/reference/view_window.md)
is the circular case, which narrows.

``` r

avail <- availability(
  surface    = c(48, 226, 67, 801),   # mean surfacing interval, seconds
  dive       = c(533, 256, 219, 88),  # mean diving interval, seconds
  window     = view_window_aerial(trackline = 51.22, speed = 51.4, slope = 0.03),
  se_surface = c(6, 20, 8, 60),
  se_dive    = c(40, 25, 20, 9),
  key        = c("Jan", "Feb", "Mar", "Apr")
)
avail
#> # A tibble: 4 × 4
#>   key   component    value      se
#>   <chr> <chr>        <dbl>   <dbl>
#> 1 Jan   availability 0.167 0.0141 
#> 2 Feb   availability 0.565 0.0337 
#> 3 Mar   availability 0.394 0.0310 
#> 4 Apr   availability 0.945 0.00918
```

### Assembling the correction

[`g0()`](https://camilleross.org/dsfit/reference/g0.md) stacks the
components, multiplies them and propagates the variance. It builds an
object and hands it off; this package does not compute abundance, so it
does not apply it.

``` r

avail$source <- "focal follows, this survey"

perception <- data.frame(
  key = NA_character_, component = "perception", value = 0.68, se = 0.09,
  source = "borrowed from a double-observer programme - not measured here"
)

g0(avail, perception)
#> <dsfit_g0>
#>   components:
#>     availability  focal follows, this survey
#>     perception    borrowed from a double-observer programme - not measured here
#> 
#>  key     g0     se    cv
#>  Jan 0.1133 0.0178 0.157
#>  Feb 0.3843 0.0558 0.145
#>  Mar 0.2679 0.0412 0.154
#>  Apr 0.6424 0.0853 0.133
#> 
#>   Divide abundance by g0, and propagate cv. This package does not
#>   apply it: correction happens at the abundance step.
```

Three things that object refuses to let you do quietly:

**It will not default to 1.** Leave perception out and it says so, in a
warning and again every time it prints:

``` r

g0(avail)
#> Warning: No "perception" component: it is being left at 1. If that is
#> deliberate, say so where this correction is used - correcting for one component
#> while believing you have corrected for another is how these estimates go wrong
#> by a factor rather than a percentage.
#> <dsfit_g0>
#>   assumed 1:   perception
#>   components:
#>     availability  focal follows, this survey
#> 
#>  key     g0     se    cv
#>  Jan 0.1667 0.0141 0.085
#>  Feb 0.5652 0.0337 0.060
#>  Mar 0.3940 0.0310 0.079
#>  Apr 0.9447 0.0092 0.010
#> 
#>   Divide abundance by g0, and propagate cv. This package does not
#>   apply it: correction happens at the abundance step.
```

**It will not accept a component without a standard error.** The CV of
`g(0)` routinely dominates the CV of abundance, so dropping it inverts
which uncertainty matters.
[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
returns `se = NA` when you give it no interval standard errors, and
[`g0()`](https://camilleross.org/dsfit/reference/g0.md) refuses that —
the two rules close a loop.

**It will not join components keyed on different things.** Availability
by month and perception by year is the common case, and it has no
correct join: the answer is a value per month-year. Build that cross
product deliberately.

Under independence the squared CVs add, `CV(g0)² = Σ CV(xᵢ)²`, so the
least precise component sets the floor. A perception estimate with a 30%
CV cannot be rescued by an availability estimate with a 2% one.

## Where this stops

Divide abundance by `g(0)` and carry its CV through — at the abundance
step, in the analysis layer, not here. `dsfit` fits detection functions
and reports what they imply; density surface models and their covariate
stack live in their own package.

And it does not estimate `g(0)`. Correcting for one component while
believing you have corrected for another is how these estimates go wrong
by a factor rather than a percentage, which is why the components are
named separately all the way through.
