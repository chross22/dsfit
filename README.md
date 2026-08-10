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

## What it will not do

Estimate `g(0)`. Every fit here conditions on the animal having been available
and seen, and a mis-specified `g(0)` shifts every candidate in a set by the same
factor — the ranking looks untouched while the density is wrong. It cannot be
estimated from a standard NARWC extract, which records neither dive data nor the
double-observer structure mark-recapture needs.

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
