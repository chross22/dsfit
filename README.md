# dsfit

Fit and compare detection functions for line-transect distance sampling — and
compare them on the quantities that decide abundance, not on AIC alone.

Takes the flatfile produced by [`distsamp`](https://github.com/chross22/distsamp)
(or any equivalent table), fits a set of candidate detection functions through
`mrds`, and returns a selection table carrying effective strip half-width and its
CV alongside each AIC.

> **Status: skeleton.** The model set, the sweep, and the input guards are built
> and tested. The `g(0)` correction slot, the `dsm` handoff, and the snapshot
> regression tests are not yet written. See [docs/01-plan.md](docs/01-plan.md).

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

## The gamma key

`Distance::ds()` offers `hn`, `hr`, and `unif`. The gamma key — **unimodal, with
its peak away from zero** — is available in `mrds` alone, which is one reason
this package fits through `mrds` directly rather than wrapping `ds()`.

It is the natural shape for an aerial platform that cannot see the water directly
beneath it. But a gamma key and a **left truncation model the same phenomenon**,
so applying both accounts for the blind spot twice. `sweep_models()` warns when
you do.

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
