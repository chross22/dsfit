# Fitting and comparing detection functions

This guide walks through one basic dsfit run, from a table of detection
distances to a g(0) correction. dsfit fits a set of candidate detection
functions with [`mrds::ddf()`](https://rdrr.io/pkg/mrds/man/ddf.html)
and puts them in one table for comparison.

Some terms used throughout:

- The **detection function** g(x) is the probability that an animal at
  perpendicular distance x from the trackline is detected.
- The **average detection probability** `p` is the mean of g(x) over the
  truncation interval: the fraction of animals in the strip that were
  seen.
- The **effective strip half-width** (`esw`) is `p` times the width of
  the truncation interval. It is the half-width of a strip in which
  every animal would be seen, giving the same expected number of
  detections. Density is estimated by dividing by it.
- **g(0)** is the probability of detecting an animal on the trackline.
  Every fit here assumes g(0) = 1. Section 7 covers how to correct for
  that.

The same steps are available as numbered example scripts, one per step:
`01-load-data.R`, `02-check-data.R`, `03-fit-detection-functions.R`,
`04-inspect-and-plot.R` and `05-g0-correction.R`. They are installed
with the package, in the directory returned by
`system.file("scripts", package = "dsfit")`.

``` r

library(dsfit)
```

## 1. The data and its required columns

dsfit takes a data frame with one row per detection. The columns it uses
are:

| column | required | meaning |
|----|----|----|
| `distance` | yes, for exact distances | perpendicular distance from the trackline |
| `distbegin`, `distend` | yes, for binned distances (instead of `distance`) | the distance bin each detection was recorded in |
| `object` | no | a unique row identifier; added if missing |
| `size` | no | group size, needed to estimate individuals rather than groups |
| covariates | no | any column named in a model formula, e.g. `beaufort` |

A row with a missing distance records a sample with no detections. These
rows are dropped before fitting.

No detection data ships with the package, so this guide uses a simulated
aerial survey. Detection falls off with distance, falls off faster in
rougher seas (higher Beaufort sea state), and is zero within 100 m of
the trackline. That 100 m gap is a blind spot: the observers cannot see
directly beneath the aircraft.

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

## 2. Checking the data

[`detection_structure()`](https://camilleross.org/dsfit/reference/detection_structure.md)
reads the table and reports which analyses it supports. It fits nothing.

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
#>     fit covariate models            candidates: beaufort; recorded as survey conditions: beaufort
#>     carry group size to abundance   `size` present
#> 
#>   cannot:
#>     estimate perception bias  no double-observer structure: no `observer` or `detected` column. Perception needs two independent teams and cannot be recovered from a single-observer survey at any sample size
#>     estimate availability     not estimable from distances by construction: a submerged animal is missed at every distance equally. Compute it from dive data with availability()
#> 
#>   Structure only. That something is supported is not a reason to do it.
```

Two lines in this report affect the analysis:

- **The nearest detection is at 100 m, not near 0.** This is the blind
  spot. It decides the model set in Section 3.
- **Perception bias cannot be estimated.** Perception bias is the chance
  of missing an animal that was at the surface and in view. Estimating
  it needs two independent observer teams, and this survey has one. More
  data does not change this. It is handled in Section 7.

[`diagnose_sweep()`](https://camilleross.org/dsfit/reference/diagnose_sweep.md)
is the second check. It takes the data and a model set (Section 3) and
runs the checks
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
runs before fitting, without fitting anything. It reports problems such
as a formula naming a missing column, a truncation that drops most of
the data, or a model set that handles the blind spot twice. The
availability vignette shows its output
([`vignette("availability", package = "dsfit")`](https://camilleross.org/dsfit/articles/availability.md)).

## 3. Building a model set

[`model_set()`](https://camilleross.org/dsfit/reference/model_set.md)
lists the candidate models: one row per combination of key function,
adjustment terms and covariate formula. The key function is the basic
shape of g(x): half-normal (`"hn"`), hazard-rate (`"hr"`), uniform
(`"unif"`) or gamma (`"gamma"`).

A blind spot can be handled in one of two ways:

- the **gamma** key, whose g(x) rises to a peak away from zero, or
- **left truncation** (`left =` in
  [`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)),
  which removes detections closer than a given distance before fitting.

Use one or the other. Using both accounts for the blind spot twice, and
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
warns if a model set containing gamma is fitted with `left`. The first
set includes gamma and no left truncation:

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

## 4. Running the sweep

[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
fits every model in the set to the same data. The truncation distance is
a single value for the whole sweep, because AIC can only compare models
fitted to the same data.

``` r

sw <- sweep_models(detections, models, truncation = 2500)
#> Fitted 3 of 3 candidate models to 1999 detections, truncation 2500.
#>   g(0) = 1 is assumed; see `?sweep_models`.
```

Printing the result gives a short summary of the fits:

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

## 5. Reading the selection table

[`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md)
returns one row per fitted model, sorted by AIC.

``` r

selection_table(sw)[, c("model_id", "delta_aic", "p", "p_cv", "esw", "cvm_p")]
#> # A tibble: 3 × 6
#>   model_id delta_aic     p   p_cv   esw      cvm_p
#>   <chr>        <dbl> <dbl>  <dbl> <dbl>      <dbl>
#> 1 gamma           0  0.255 0.0164  636. 0.00371   
#> 2 hr            377. 0.355 0.0154  887. 0.00000308
#> 3 hn            407. 0.296 0.0177  739. 0.00132
```

| column | meaning |
|----|----|
| `model_id` | the model’s name from [`model_set()`](https://camilleross.org/dsfit/reference/model_set.md) |
| `key`, `adjustment`, `order`, `formula` | the model specification |
| `converged`, `n_par`, `aic` | whether the fit converged, number of parameters, AIC |
| `delta_aic` | AIC minus the lowest AIC in the table; 0 for the best model |
| `p`, `p_se`, `p_cv` | average detection probability, its standard error, and its coefficient of variation (CV, the standard error divided by the estimate) |
| `esw` | effective strip half-width, in the units of `distance` |
| `cvm_p` | p-value of the Cramér-von Mises goodness-of-fit test, for exact distances |
| `chisq_p` | p-value of the chi-square goodness-of-fit test; use this one for binned distances, where `cvm_p` is `NA` |

Read three things for each model, not only `delta_aic`:

1.  **`delta_aic`** ranks the models against each other. It does not say
    whether any of them fits.
2.  **`cvm_p`** (or `chisq_p` for binned data) tests whether the model
    fits the observed distances. A small value (below 0.05) means it
    does not.
3.  **`esw`** is what abundance is divided by, so a difference in `esw`
    between models is a proportional difference in estimated density.
    Here `esw` runs from about 640 m to about 890 m, so the choice of
    model changes the density estimate by about 40%.

### Goodness of fit and the blind spot

All three models fail the goodness-of-fit test, including the best by
AIC. The blind spot in this survey is a sharp cutoff at 100 m. No key
function, gamma included, has a sharp cutoff, so all of them misfit the
near distances. The fix is to remove the cutoff from the data with left
truncation rather than to try more keys.

### Refitting with left truncation

Gamma is left out of this set because left truncation already handles
the blind spot.

``` r

sw_left <- sweep_models(detections, model_set(c("hn", "hr")),
                        truncation = 2500, left = 100, quiet = TRUE)
selection_table(sw_left)[, c("model_id", "delta_aic", "p", "esw", "cvm_p")]
#> # A tibble: 2 × 5
#>   model_id delta_aic     p   esw    cvm_p
#>   <chr>        <dbl> <dbl> <dbl>    <dbl>
#> 1 hn             0   0.244  584. 0.697   
#> 2 hr            71.8 0.295  707. 0.000352
```

The half-normal now passes the goodness-of-fit test. The hazard-rate
still fails, so the test is distinguishing between models, not only
approving the best one.

The AIC values in this table cannot be compared with those in the
previous table, because the two sweeps were fitted to different data.
The `esw` values can be compared.

### Adding a covariate

A model with a detection covariate, such as `~beaufort`, fits a
different g(x) for each covariate value. Its `p` is then an average over
the detected animals.

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
#> 1 hn beaufort       0   0.241 0.0179  577.
#> 2 hn               26.2 0.244 0.0177  584.
```

Sea state affects detection in the simulation, and the model with
`beaufort` has the lower AIC. On real data, drop a covariate that does
not improve AIC.

## 6. Plotting the chosen detection function

The fitted models are stored in `sw$fits`, named by `model_id`. dsfit
provides an `effect_estimates()` method for these fits so that
[fancyfx](https://github.com/chross22/fancyfx) can plot them. The plot
shows the fitted g(x) with a delta-method standard-error band, over a
rug of the observed distances.

``` r

fancyfx::plotEffects(sw$fits[["hn"]], dat = detections, var = "distance")
```

To get the curve as a data frame instead, use
[`effect_estimates_ddf()`](https://camilleross.org/dsfit/reference/effect_estimates_ddf.md).

For a covariate model, the area under the plotted curve does not exactly
equal the `esw` in the table. The curve is the plain mean of g(x) over
the detected animals. `mrds` computes `p`, and so `esw`, as a
Horvitz-Thompson mean, which weights each animal by the inverse of its
detection probability. The two are equal only when detection probability
is the same for every animal.

## 7. Handing off g(0)

Every fit above assumes g(0) = 1: an animal on the trackline is always
seen. If that is wrong, every model’s `p` is wrong by the same factor.
The AIC ranking does not change, so model selection cannot detect it.
The correction is built separately and applied when abundance is
calculated.

Three separate causes make g(0) less than 1:

| component | what it is | how it is handled |
|----|----|----|
| blind spot | the observers cannot see beneath the aircraft | gamma key or `left`, as in Sections 3 to 5 |
| availability | the animal was underwater as the aircraft passed | computed from dive data with [`availability()`](https://camilleross.org/dsfit/reference/availability.md) |
| perception | the animal was at the surface and missed | estimated from a second observer team, elsewhere |

### Availability

An animal that is underwater for the whole pass is missed at every
distance equally. This leaves no trace in the distances, so availability
cannot be estimated from them. It is computed instead from dive
intervals and the time the aircraft keeps a point in view.

The package includes availability measured for right whales in Cape Cod
Bay by Ganley et al. (2019):

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

Availability was 0.27 in January and 0.91 in April. A single value for
the whole season would be wrong by up to a factor of three in some
months.

[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
computes availability from mean surfacing and diving intervals and a
time-in-view window. The window here is from
[`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md),
which models an aircraft’s field of view; at the trackline it is 51.22
s. The intervals are Ganley et al.’s means. The standard errors are for
illustration.

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

These values differ from Ganley et al.’s published figures. Theirs are
bootstrap medians that were not computed at a common window; the
availability vignette explains this.

### Assembling the correction

[`g0()`](https://camilleross.org/dsfit/reference/g0.md) combines the
components by multiplying them, and propagates their standard errors. It
returns a correction object for use at the abundance step. dsfit does
not compute abundance, so it does not apply the correction. The `source`
column records where each component came from and is printed with the
result.

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

[`g0()`](https://camilleross.org/dsfit/reference/g0.md) enforces three
rules:

1.  **A missing component is named, not set silently to 1.** Leaving out
    perception gives a warning, and the printed object lists it under
    “assumed 1”:

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

2.  **Every component needs a standard error.** The CV of g(0) is often
    the largest part of the CV of abundance.
    [`availability()`](https://camilleross.org/dsfit/reference/availability.md)
    returns `se = NA` when no interval standard errors are given, and
    [`g0()`](https://camilleross.org/dsfit/reference/g0.md) refuses
    that.

3.  **Components must share the same keys.** For example, availability
    by month and perception by year cannot be joined. Build a value per
    month-year for both instead.

If the components are independent, their squared CVs add: CV(g0)² = Σ
CV(xᵢ)². The least precise component therefore sets the minimum CV. A
perception estimate with a 30% CV keeps the CV of g(0) at 30% or more,
however precise availability is.

## Where dsfit stops

dsfit fits detection functions, reports what they imply, and builds the
g(0) correction. Dividing abundance by g(0) and carrying its CV through
is done at the abundance step, outside this package. dsfit does not
estimate g(0) itself. See
[`vignette("availability", package = "dsfit")`](https://camilleross.org/dsfit/articles/availability.md)
for the details of availability and g(0).

## References

Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what
there is: estimating North Atlantic right whale *Eubalaena glacialis*
local abundance. *Endangered Species Research* 38:101-113.
<https://doi.org/10.3354/esr00938>
