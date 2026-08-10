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

## Not yet built

`g(0)` correction, the `dsm` handoff, snapshot regression tests over the
selection table, and a vignette.
