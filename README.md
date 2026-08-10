# dsfit

<!-- badges: start -->
[![R-CMD-check](https://github.com/chross22/dsfit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chross22/dsfit/actions/workflows/R-CMD-check.yaml)
[![check-citations](https://github.com/chross22/dsfit/actions/workflows/check-citations.yaml/badge.svg)](https://github.com/chross22/dsfit/actions/workflows/check-citations.yaml)
<!-- badges: end -->

Fit and compare detection functions for line-transect distance sampling — and
compare them on the quantities that decide abundance, not on AIC alone.

Takes the flatfile produced by [`distsamp`](https://github.com/chross22/distsamp)
(or any equivalent table), fits a set of candidate detection functions through
`mrds`, and returns a selection table carrying effective strip half-width and its
CV alongside each AIC.

> **Status: skeleton.** The model set, the sweep, the input guards, and
> snapshot regression tests pinning the selection are built and tested. The
> `g(0)` correction slot and the `dsm` handoff are not yet written. See
> [docs/01-plan.md](docs/01-plan.md).

## Install

```r
# install.packages("remotes")
remotes::install_github("chross22/dsfit")
```

## Use

```r
library(dsfit)

sw <- sweep_models(
  detections,                              # needs `distance`, or distbegin/distend
  model_set(key = c("hn", "hr", "gamma")),
  truncation = 400
)

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
```

Read that table across, not down. `hn` wins on AIC, but `esw` spans 120 to 151
across the three — a 26% spread in the strip width that abundance is divided by.

## Three things it insists on

**Truncation is one value per sweep.** AIC compares models fitted to the same
data, and changing the truncation changes which detections are in the
likelihood. So it is an argument, not a column. To compare truncations, compare
*sweeps*, and look at whether the chosen model and its `p̄` are stable.

**Binned and exact distances cannot be ranked together.** `STRIP`-derived
distances are intervals; angle- and position-derived distances are points. The
likelihoods differ, so a table containing both is a survey-era boundary rather
than a model set. Refused, with the reason.

**Effective strip half-width sits next to the AIC.** It is what propagates into
density. Two models within 2 AIC of one another can imply materially different
abundance, and nothing in an AIC column shows that.

## Goodness of fit, and which test you get

A model can top the AIC ranking and still fail a goodness-of-fit test, so one
sits in the table. Which one follows the data rather than a preference:
Cramér-von Mises (`cvm_p`) for exact distances, and a chi-square over the
survey's own bins (`chisq_p`) for binned ones, which have no empirical
distribution for CvM to test. `print()` shows whichever applies.

Both columns are named for their test, because a p-value read without knowing
what produced it is worse than no p-value. `chisq_p` is filled in for exact fits
too, but over cutpoints `mrds` chooses rather than bins anyone surveyed — read
`cvm_p` there.

## Snapshot tests, and why they are the point

The sweep *selects* a model. The snapshots in
`tests/testthat/test-selection-snapshot.R` pin the selection — the whole table,
not just the winner — so that an `mrds` upgrade cannot quietly change which
model wins or by how much. A swap between models three and four is the same
change arriving early enough to look at.

This is the requirement that made the middle layer a package. Writing them
turned up a live bug in which a candidate that failed to fit shifted every model
below it onto the wrong row, producing a wrong selection table rather than an
error.

## The gamma key

`Distance::ds()` offers `hn`, `hr`, and `unif`. The gamma key — **unimodal, with
its peak away from zero** — is available in `mrds` alone, which is one reason
this package fits through `mrds` directly rather than wrapping `ds()`.

It is the natural shape for an aerial platform that cannot see the water directly
beneath it. But a gamma key and a **left truncation model the same phenomenon**,
so applying both accounts for the blind spot twice. `sweep_models()` warns when
you do.

## Plotting

`effect_estimates.ddf()` makes a fitted detection function plottable with
[fancyfx](https://github.com/chross22/fancyfx):

```r
fancyfx::plotEffects(sw$fits[["hn"]], dat = detections, var = "distance")
```

which gives the fitted g(x) with a delta-method ribbon and a rug of the observed
distances above it — the standard detection-function diagnostic, and fancyfx's
premise applied almost exactly. `fancyfx` is a `Suggests`, so it is optional.

One thing worth knowing: for a **covariate** model the area under that curve is
not the `esw` in the selection table. `mrds` computes `average.p` as a
Horvitz-Thompson mean, weighting each animal by the inverse of its own detection
probability; the plotted curve is the plain mean over the animals as observed.
They coincide only when detection probability is constant.

## Availability

`availability()` computes the availability component of `g(0)` — the probability
an animal was at the surface to be seen while the platform had it in view.

```r
availability(
  surface = example_dive_intervals$surface,   # mean surfacing interval, s
  dive    = example_dive_intervals$dive,      # mean diving interval, s
  window  = view_window(radius = 300, speed = 50),
  key     = as.character(example_dive_intervals$month)
)
#> # A tibble: 6 × 4
#>   key   component    value    se
#>   <chr> <chr>        <dbl> <dbl>
#> 1 Dec   availability 0.245    NA
#> 2 Jan   availability 0.286    NA
#> 3 Feb   availability 0.361    NA
#> 4 Mar   availability 0.477    NA
#> 5 Apr   availability 0.683    NA
#> 6 May   availability 0.849    NA
```

That spread — 0.25 to 0.85 across one season — is the entire argument against a
constant.

### The equation

Following [Laake and Borchers (2004)](https://global.oup.com/), for a mean
surfacing interval `E(s)`, a mean diving interval `E(d)`, and a window `w`
during which the animal is in view:

```
        E(s) + E(d) · (1 − exp(−w / E(d)))
  a =  ────────────────────────────────────
                  E(s) + E(d)
```

Two limits are worth holding onto, and both are pinned as tests:

- **`w → 0`** — an instantaneous glance. The exponential term vanishes and
  `a = E(s) / (E(s) + E(d))`, the plain proportion of time spent at the surface.
- **`w → ∞`** — watch forever, and every animal surfaces eventually, so `a → 1`.

The middle term is the extra chance that an animal which was *down* when the
window opened comes up before it closes. That is why a slower aircraft or a
wider field of view raises availability without the whales changing anything.

The window itself comes from geometry. `view_window()` treats the viewing area
as a circle of radius `r` crossed at speed `v`, so a point at perpendicular
distance `x` is in view along a chord:

```
  w(x) = 2 · sqrt(r² − x²) / v
```

which narrows with distance off the trackline and closes at the edge of view.
It is the simplest useful geometry and will not fit every platform — a
forward-looking window or an obscured belly gives a different `w` — so
`availability()` takes `window` directly and a measured one can be substituted.
Measuring it, as Ganley did for the aircraft they flew, beats deriving it.

### It computes; it does not estimate

An animal submerged while the aircraft passes is missed at **every perpendicular
distance equally**. To first order availability is a pure scale factor on g(x):
it does not dent the near-zero end of the distance distribution, does not change
its shape, and leaves no signature for a likelihood to find. That is the same
reason a mis-specified `g(0)` shifts every candidate in a selection table by the
same factor and leaves the ranking looking untouched.

So there is nothing in the distances to estimate it from, and **none of the
inputs come from the survey being corrected**. The intervals come from focal
follows, tagging or drone work; the window comes from the platform. Both arrive
with their own uncertainty, which is the point — it stays visible instead of
being absorbed into a constant.

### Why it is a function and not a number

[Ganley et al. (2019)](https://doi.org/10.3354/esr00938) measured this for right
whales in Cape Cod Bay with focal follows and the aircraft's field of view:
**availability varied by month between 0.27 and 0.85**, tracking the depth of
the copepod layer the whales were feeding on, while detection probability varied
separately by year between 0.43 and 0.87.
[Roberts et al. (2024)](https://doi.org/10.3354/meps14547) arrive at the same
place from the other direction, correcting perception and availability per
platform, per team and per conditions across 11 institutions and 2.9 million km
of effort.

A single plausible-looking figure is therefore the most dangerous input this
package accepts. Something like 0.83 sits comfortably inside *both* of those
ranges, which is exactly what makes it treacherous: it is one component's value
under one set of conditions. If availability is ~0.5 in a given month and
perception ~0.7, the combined `g(0)` is **~0.35** — using 0.83 there does not
overestimate density by a fifth, it understates the correction by more than a
factor of two. And because it scales every candidate equally, no amount of model
selection reveals it.

### The output shape, and the standard error

The result is a `key` / `component` / `value` / `se` table — the shape a `g(0)`
correction stacks from, so availability rows and perception rows sit together and
it stays visible which have been applied.

Availability is deterministic given its inputs, so its uncertainty comes from
them: pass `se_surface` and `se_dive` and they propagate by delta method.
Passing one without the other is an **error**, not a half-propagated variance —
reporting a standard error smaller than the truth is worse than reporting none.
Passing neither gives `se = NA` rather than an invented precision, because this
is the quantity whose CV routinely dominates the CV of abundance.

Where the raw focal follows are in hand, bootstrapping them beats this, since it
does not assume the interval means are jointly normal. That is what Ganley did.

### The example data is mock

`example_dive_intervals` is **invented**. It is not a measurement, and it is
specifically *not* Ganley's Cape Cod Bay table, which is not open access. What
is real is the pattern it reproduces — dive times falling and surface times
rising through a season, as the food moves up — so that a worked example shows
why availability cannot be a constant. Substitute your own.

## What it will not do

Estimate `g(0)`. Every fit here conditions on the animal having been available
and seen, and a mis-specified `g(0)` shifts every candidate in a set by the same
factor — the ranking looks untouched while the density is wrong. It cannot be
estimated from a standard NARWC extract, which records neither dive data nor the
double-observer structure mark-recapture needs.

`availability()` is not an exception to this: it is a calculation from external
dive data, and it covers one of the three components. Perception needs a
double-observer protocol, and the geometric blind spot is handled by left
truncation or the gamma key.

Correct for it at the abundance step, from external sources, with its standard
error propagated and its components named separately: the geometric blind spot,
availability, and perception are three different things, and correcting for one
while believing you have corrected for another is how these estimates go wrong by
a factor rather than a percentage.

## Where this sits

```
distsamp   NARWC survey data  ->  effort segments + detection distances
   |
dsfit      detection functions: model set, sweep, selection table
   |
analysis   which years, which truncation, which covariates, the report
```

The middle layer is a package because its logic needs tests. The outer layer is
not, because its choices change every run.

## Citing dsfit

```r
citation("dsfit")
```

That returns up to three entries — the package, the `mrds` toolchain every fit
comes off, and the gamma key function if you selected a gamma model:

```
Ross, C. dsfit: Fit and Compare Detection Functions for Line-Transect Survey
Data. R package. https://github.com/chross22/dsfit

Miller, D.L., Rexstad, E., Thomas, L., Marshall, L. and Laake, J.L. (2019)
Distance sampling in R. Journal of Statistical Software 89(1):1-28.
doi:10.18637/jss.v089.i01

Becker, E.A. and Quang, P.X. (2009) A gamma-shaped detection function for
line-transect surveys with mark-recapture and covariate data. Journal of
Agricultural, Biological and Environmental Statistics 14:207-223.
doi:10.1198/jabes.2009.0013
```

The package version comes from `DESCRIPTION` at install time, so it is always
the version you actually have. Record the truncation and the model set as well:
a selection table is not reproducible without them.

## References

What each source is relied on for is in [tools/citations.csv](tools/citations.csv).
They are checked monthly by CI — DOIs still resolve and still describe the paper
cited.

Becker, E.A. and Quang, P.X. (2009) A gamma-shaped detection function for
line-transect surveys with mark-recapture and covariate data. *Journal of
Agricultural, Biological and Environmental Statistics* 14:207–223.
<https://doi.org/10.1198/jabes.2009.0013>
— *the gamma key, and the aerial-survey case it was developed for.*

Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L. and
Thomas, L. (2001) *Introduction to Distance Sampling: Estimating Abundance of
Biological Populations.* Oxford University Press.
— *the standard line-transect reference.*

Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what there is:
estimating North Atlantic right whale *Eubalaena glacialis* local abundance.
*Endangered Species Research* 38:101–113.
<https://doi.org/10.3354/esr00938>
— *availability from focal follows, and the monthly variation that rules out a
constant.*

Laake, J.L. and Borchers, D.L. (2004) Methods for incomplete detection at
distance zero. In *Advanced Distance Sampling*, pp. 108–189. Oxford University
Press.
— *why `g(0)` is not estimable from single-observer data.*

Marques, F.F.C. and Buckland, S.T. (2004) Covariate models for the detection
function. In *Advanced Distance Sampling*, pp. 31–47. Oxford University Press.
— *covariate detection functions, which is what `mcds` fits.*

Miller, D.L., Rexstad, E., Thomas, L., Marshall, L. and Laake, J.L. (2019)
Distance sampling in R. *Journal of Statistical Software* 89(1):1–28.
<https://doi.org/10.18637/jss.v089.i01>
— *the `Distance` and `mrds` toolchain this package fits through.*

Roberts, J.J., Yack, T.M., Fujioka, E., Halpin, P.N., Baumgartner, M.F. and
others (2024) North Atlantic right whale density surface model for the US
Atlantic evaluated with passive acoustic monitoring. *Marine Ecology Progress
Series* 732:167–192. <https://doi.org/10.3354/meps14547>
— *corrections applied per platform, team and conditions, at continental scale.*
