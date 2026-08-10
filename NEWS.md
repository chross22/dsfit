# dsfit 0.0.0.9000

First skeleton. The detection-function half of the fitting layer described in
`distsamp`'s `docs/07-fitting-architecture.md`.

## Fitting

* `model_set()` — expands key functions, adjustment series and orders, and
  covariate formulas into a grid of candidates. Supports the **gamma** key,
  which `Distance::ds()` does not offer.
* `sweep_models()` — fits every candidate through `mrds::ddf()` at one
  truncation, so every AIC in the table comes off the same likelihood
  machinery.
* `selection_table()` — the ranked result, carrying `p`, `p_se`, `p_cv`, `esw`
  and a Cramér-von Mises p-value next to each AIC.
* `prepare_distance_data()` — the input guards, exported because they are worth
  running before committing to a model set.

## What it refuses, and why

* **Point and interval distances in one sweep.** Binned and exact fits have
  different likelihoods, so their AICs are not comparable.
* **An unbounded top bin.** `distend` of `Inf` cannot be fitted; every `STRIP`
  scheme's top bin is open.
* **Bins that do not tile**, which do not define a set of cutpoints.
* **Bins that stop short of the truncation**, which leave unaccounted effort
  between the last bin and the truncation width.
* **A uniform key with no adjustment terms**, which is a strip transect rather
  than a detection function.

## Warnings rather than refusals

* The gamma key together with left truncation. Both describe a platform that
  cannot see beneath itself, so applying them together counts the blind spot
  twice.

## Plotting, via fancyfx

* `effect_estimates.ddf()` — an [fancyfx](https://github.com/chross22/fancyfx)
  method for `mrds` detection functions, so a fitted `ddf` plots with
  `plotEffects()` and `plotRugs()` like any other model. `fancyfx` sends
  non-GAMs to `marginaleffects`, which cannot introspect a `ddf`; the generic is
  designed to be extended from outside, and this is that extension.

  The result is the standard detection-function diagnostic — fitted g(x) with a
  rug of the observed distances — which is `fancyfx`'s premise almost exactly.

* Intervals are delta-method: the Jacobian of g with respect to the fitted
  parameters, combined with the inverse Hessian, clamped to `[0, 1]`.
* Covariate models average the curves over the covariate values observed rather
  than evaluating g() once at mean covariates. `at` fixes named covariates and
  averages over the rest.
* `effect_estimates_ddf()` is the same function under a plain name, so it can
  be called without `fancyfx` installed.

Two things documented because they would otherwise be discovered the hard way:

* **The area under the curve is not `esw` for a covariate model.** `mrds`
  reports `average.p` as a Horvitz-Thompson mean, `n / sum(1/p)`, which weights
  each animal by the inverse of its own detection probability. The plotted curve
  is the plain arithmetic mean over the animals as observed. They coincide only
  when detection probability is constant.
* `"auto"` is part of `fancyfx`'s argument vocabulary and arrives here
  routinely, since `plotEffects()` forwards its own unresolved defaults.
  `match.arg()` rejects those outright, so the choices are resolved by hand.

## Not yet built

`g(0)` correction, the `dsm` handoff, snapshot regression tests over the
selection table, plotting detectability against a covariate (a different
quantity from g(x), and not the same plot), and a vignette.
