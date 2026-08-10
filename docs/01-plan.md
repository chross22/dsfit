# dsfit: the plan

Why this package exists, what is built, and what is deliberately left out.

The architecture it implements was settled in `distsamp`'s
[docs/07-fitting-architecture.md](https://github.com/chross22/distsamp/blob/main/docs/07-fitting-architecture.md).
That document is the reasoning; this one is the state.

---

## 1. Three layers

| Layer | Holds | Form |
|---|---|---|
| `distsamp` | NARWC ingest, effort, segmentation, right-angle distances | package |
| **`dsfit`** | the model-selection sweep, gamma support, GOF and `p̄`, `g(0)`, the `dsm` handoff | **package** |
| analysis repo | which years, which truncation, which covariates, the report | `targets` + `renv` |

The middle layer is a package rather than scripts because it is the part with
logic that needs tests. An automated comparison over detection-function
configurations is exactly where a silent error costs a density estimate rather
than throwing, and test infrastructure in a scripts repository does not get run.

The outer layer is not a package because its choices change every run, and
freezing them into one is how a study-area constant outlives its study.

## 2. Why `mrds` and not `Distance::ds()`

`Distance` is the right tool for fitting one model interactively, and this
package does not replace it. Two reasons to go a level down for a sweep:

- **Coverage.** `ds()` accepts `key = c("hn", "hr", "unif")`. The gamma key is
  in `mrds::ddf(key = "gamma")` only, and it is the natural shape for an aerial
  platform with a blind spot beneath it — which is precisely the Skymaster case
  documented in `distsamp`'s `strip_distance()`.
- **Comparability.** Every AIC in a selection table has to come off the same
  likelihood machinery. `ds()` applies its own truncation handling and
  monotonicity constraints; a table mixing `ds()` and `ddf()` fits risks
  presenting those differences as model differences.

## 3. What is built

`model_set()`, `sweep_models()`, `selection_table()`, `prepare_distance_data()`,
and a `print()` method. 137 tests, six of them snapshots over the selection
table.

Left truncation is here too, as `sweep_models(left = )` — the statistically
correct treatment of a blind spot beneath the platform, rather than shifting
bins. `distsamp`'s next-steps list used to carry it; it is fitting-side, so it
belongs to this layer, and that list now points here. The same applies to
truncation and binned fitting, and to the `g(0)` slot below: one list owns each
item.

### The conditions a comparison needs

Each of these is enforced rather than documented, because each one silently
invalidates a ranking rather than failing:

| Condition | Why | Behaviour |
|---|---|---|
| One truncation per sweep | AIC compares models on the same data | `truncation` is an argument, not a column |
| No mixing point and interval distances | Different likelihoods | error |
| No unbounded top bin | Cannot be fitted at all | error |
| Bins must tile | Otherwise they define no cutpoints | error |
| Bins must reach the truncation | Otherwise there is unaccounted effort between the last bin and the width | error |
| Uniform key needs adjustments | A flat detection function is a strip transect | warn and drop |
| Gamma and left truncation are alternatives | Both model the blind spot | warn |

### Ranking on more than AIC

The selection table carries, per model: `aic`, `delta_aic`, `n_par`, `p`,
`p_se`, `p_cv`, `esw`, `cvm_p`, `chisq_p`, `converged`.

`esw` is there because it is what propagates into density, and two models within
2 AIC of one another can imply materially different abundance. On the package's
own example, `hn`, `gamma` and `hr` span 120 to 151 in `esw` — a 26% spread —
while the top two AICs differ by less than the third's advantage suggests.

A goodness-of-fit p-value is reported because a model can top the AIC ranking
and still fail one. Its null is that the fitted function describes the observed
distances, so a small p-value is evidence against the model. Which test that is
follows the data: Cramér-von Mises for exact distances, and a chi-square over
the survey's bins for binned ones, which have no empirical distribution for CvM
to test. Both columns are carried under their own names, because a p-value read
without knowing its test is worse than no p-value.

### The selection is pinned

Snapshot regression tests over the whole table, which is the requirement that
made this a package rather than scripts. The sweep *selects* a model; the
snapshots *pin* the selection, so that an `mrds` or `Distance` upgrade cannot
silently change which model wins. The whole table is pinned rather than the
winner alone — a swap between models 3 and 4 is the same change arriving early.

Six of them, over key functions, adjustment series and orders, the binned
likelihood, covariate models, left truncation, and what the sweep prints. Values
are rounded to a precision an optimiser can be expected to reproduce across
platforms; the fixtures behind them are fixed in
`tests/testthat/helper-fixtures.R`.

Writing them turned up a live bug: a candidate that failed to fit was dropped
from the list of fits rather than emptied, so every model after it inherited the
previous one's AIC and `esw`. A wrong selection table, silently — which is the
class of error this layer exists to catch.

## 4. What is not built

In the order it is worth doing.

1. **The `g(0)` correction slot** — built. `availability()` computes the first
   component and `g0()` assembles them: a **keyed table**, not a value, with
   rows of `component` × `value` × `se`, multiplied together and the variance
   propagated by delta method.

   `g0()` builds an object and hands it off rather than applying anything,
   because this package fits detection functions and does not compute
   abundance — the correction belongs one layer out. What it can do from here
   is make the correction impossible to get wrong quietly: it cannot be built
   without named components, cannot be built without standard errors, and
   names any absent component in a warning and again in its own printout.

   Under independence the delta method gives the rule worth remembering:
   squared CVs add, `CV(g0)² = Σ CV(xᵢ)²`, so the least precise component sets
   the floor. Components keyed on different things — availability by month,
   perception by year — are refused rather than joined, because the honest
   answer there is a value per month-year.

   What is left is **perception**, which needs either a double-observer trial
   (item 2) or values from the literature. `docs/01-plan.md`'s own rules are
   now enforced in code rather than described here.

   Keyed, because a scalar is the error. Ganley et al. (2019) measured right
   whale availability varying by month between 0.27 and 0.85, and perception
   varying by year between 0.43 and 0.87; Roberts et al. (2024) correct per
   platform, team and conditions. A figure like 0.83 sits inside both ranges,
   which is what makes it dangerous — it is one component's value under one set
   of conditions, and using it as `g(0)` understates the correction by more than
   a factor of two when availability is ~0.5 and perception ~0.7.

   The three rules are enforced: no default and never silently 1; propagate the
   variance or refuse the correction, since its CV routinely dominates the CV of
   abundance; and name availability and perception separately, so it is visible
   which have been applied. See section 5 of the architecture document for why
   it cannot be estimated from a NARWC extract.

   Ganley's measurements ship as `ganley_surface_time` — **cited data, not a
   default**. It is Cape Cod Bay 1998–2017, driven by copepod depth in that
   specific bay; as a worked example it shows the shape of the input, and as a
   default it would be the same mistake with better provenance.
2. **An MRDS backend**, conditional on data that actually carries double-observer
   structure — with a guard that errors when it does not, rather than fitting a
   single-observer model and reporting it as though perception bias were
   handled.
3. **The `dsm` handoff**, from `distsamp::segments_as_sf(segs, "midpoints")`.
4. **A vignette**, once there is something end-to-end to walk through.

## 5. Deliberately out of scope

**Estimating `g(0)`.** Not a gap to be filled later — it cannot be estimated
from the data this layer receives. Three distinct things present as `g(0) < 1`:
the geometric blind spot, which is a property of the aircraft and is handled by
left truncation or a gamma key; availability, which needs dive data from tagging
studies; and perception, which needs a double-observer protocol. Correcting for
one while believing you have corrected for another is how these estimates go
wrong by a factor rather than a percentage.

**Wrapping `Distance::ds()`.** It is mature and well documented. A wrapper hides
its options and rots against it.
