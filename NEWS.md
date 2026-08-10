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
  and a goodness-of-fit p-value next to each AIC.
* `prepare_distance_data()` — the input guards, exported because they are worth
  running before committing to a model set.

## Pinning the selection

* Snapshot regression tests over the whole selection table, which is the
  requirement that made this layer a package. The sweep *selects* a model;
  these pin the selection, so an `mrds` upgrade cannot quietly change which
  model wins or by how much. The whole table is pinned rather than the winner
  alone — a swap between models three and four is the same underlying change
  arriving while it is still cheap to look at. Six of them: key functions,
  adjustment series and orders, the binned likelihood, covariate models, left
  truncation, and what the sweep prints about itself.

  When one fails the question is not whether the new number looks reasonable
  but what changed underneath. The fixtures behind them live in
  `tests/testthat/helper-fixtures.R` and are deliberately fixed; changing a
  seed there invalidates every snapshot, and that is a data change rather than
  a regression.

* Writing them found a real one, which is rather the point. A candidate that
  failed to fit was removed from the list of fits rather than emptied, so every
  model after it in the set inherited the previous one's AIC, `p̄` and `esw` —
  a wrong selection table rather than an error. Fixed, with a test that checks
  each row's AIC against the fit stored under that row's `model_id`.

## Goodness of fit for binned data

* `chisq_p` joins `cvm_p` in the selection table. Cramér-von Mises tests the
  empirical distribution of exact distances, and `mrds` does not compute it for
  a binned fit, which has none — so before this a binned sweep carried no
  goodness-of-fit at all, in a column that looked like it did. Binned fits now
  get a chi-square over the survey's own bins, and `print()` shows whichever
  test applies.

  `chisq_p` is filled in for exact fits too, but over cutpoints `mrds` chooses
  rather than bins the survey defined; a test whose result moves with an
  arbitrary binning is the weaker one, and `cvm_p` is what to read there. Both
  columns are named for their test so a p-value is never separated from what
  produced it.

* A chi-square over few bins runs out of degrees of freedom — three bins and a
  two-parameter key leave none — and reports `NA`. That is a converged model
  with no test available, not a failed one.

* Numeric columns of the selection table are no longer named after the models,
  so pulling one out gives a bare vector.

## The g(0) correction slot

* `g0()` assembles a correction from its named components, multiplies them,
  and propagates their variance. It **builds an object and hands it off**
  rather than applying anything: this package fits detection functions and
  does not compute abundance, so the correction belongs one layer out.

  What it can do from here is make the correction impossible to get wrong
  quietly. The three rules from `docs/01-plan.md` are now enforced rather than
  described:

  * **Never silently 1.** There is no default, and an absent component is
    named in a warning and again in the object's own printout — `assumed 1:
    perception` sits above the table every time it is printed.
  * **Propagate the variance or refuse.** A component without a standard error
    is an error. This meets `availability()`'s refusal to invent one, so the
    two rules close the loop: it will not make up a precision, and this will
    not accept none.
  * **Name the components separately**, so it stays visible which have been
    applied. An unnamed component is an error.

* Components combine multiplicatively and independently — availability comes
  from dive data, perception from a double-observer trial, so they are
  separate studies. Under independence the delta method gives the rule worth
  remembering: **squared CVs add**, `CV(g0)² = Σ CV(xᵢ)²`. The least precise
  component sets the floor, and a perception estimate with a 30% CV cannot be
  rescued by an availability estimate with a 2% one.

* Components keyed on **different things are refused rather than joined**.
  Availability by month and perception by year is the usual way to arrive
  there, and it has no correct join — the answer is a value per month-year.
  Pairing January with 1998 because both come first would be silent nonsense.

## Availability

* `availability()` — the availability component of `g(0)`, computed from mean
  surfacing and diving intervals and the time the platform keeps an animal in
  view, following Laake et al. (1997) equation 5. It returns a `key` /
  `component` / `value` / `se` table, which is the shape a `g(0)` correction
  stacks from.

  It **computes** rather than estimates, and the distinction is the whole
  point. An animal submerged while the aircraft passes is missed at every
  perpendicular distance equally, so to first order availability is a pure
  scale factor on g(x) — no dent in the near-zero end, no change in shape,
  nothing for a likelihood to find. None of the inputs come from the survey
  being corrected.

* Why it is a function rather than a constant: Ganley et al. (2019) measured
  right whale availability in Cape Cod Bay and found it varying by month
  between **0.27 and 0.85**, with perception varying separately by year
  between 0.43 and 0.87. Roberts et al. (2024) correct per platform, team and
  conditions across 11 institutions. A single representative figure is a real
  number for one month and wrong by a factor of two for others — and because
  it scales every candidate equally, model selection never reveals it.

* Standard errors propagate by delta method from `se_surface` and `se_dive`,
  numerically as in `effect_estimates_ddf()`. Supplying one without the other
  is an error rather than a half-propagated variance, and supplying neither
  gives `se = NA` rather than an invented precision.

* `view_window_aerial()`, and a correction. An **aircraft observer's field of
  view is a wedge running forward and aft**, so an animal further off the
  trackline stays in it *longer* — Ganley et al. measured time in view rising
  from about 50 s near the trackline to about 130 s at 3 km. `view_window()`'s
  circular geometry does the opposite, and using it for an aircraft gets the
  sign of the distance effect backwards. Both are kept, because both cases
  exist, but the aerial one is what a line-transect aerial survey wants.

* `ganley_surface_time`: real, cited percent surface time from 87 focal follows
  of right whales in Cape Cod Bay — 16% in January against 55% in April, as the
  copepods move up the water column. Percent surface time is the `w → 0` limit
  of `availability()`, so it is a floor on availability and the seasonal
  problem in its rawest form.

  It carries only what the paper states in text. The monthly availability
  figures are in a figure and a supplementary table, and reading bar heights
  off a chart is not a measurement, so the dataset stops where the text does.

* `example_dive_intervals`, a six-month table of surfacing and diving
  intervals. **The numbers are invented** — not measurements, and not Ganley's
  values. What is real is the pattern: dive times falling and surface times
  rising through the season, which puts `availability()` over them between
  about 0.25 and 0.85. Shipped so the worked example can show the full API,
  standard errors included.

* The test suite now reconciles against a published result. January's 16%
  surface time and reported 0.27 availability, at a measured 51.22 s in view,
  pin the dive time at about 6.1 minutes — inside the 1.30 to 8.83 min range
  the paper gives for monthly mean dive times.

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

`g(0)` correction, the `dsm` handoff, plotting detectability against a
covariate (a different quantity from g(x), and not the same plot), and a
vignette.
