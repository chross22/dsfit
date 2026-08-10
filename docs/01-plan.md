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
and a `print()` method. 77 tests.

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
`p_se`, `p_cv`, `esw`, `cvm_p`, `converged`.

`esw` is there because it is what propagates into density, and two models within
2 AIC of one another can imply materially different abundance. On the package's
own example, `hn`, `gamma` and `hr` span 120 to 151 in `esw` — a 26% spread —
while the top two AICs differ by less than the third's advantage suggests.

Cramér-von Mises is reported because a model can top the AIC ranking and still
fail a goodness-of-fit test. Its null is that the fitted function describes the
observed distances, so a small p-value is evidence against the model.

## 4. What is not built

In the order it is worth doing.

1. **Snapshot regression tests over the selection table.** The requirement that
   made this a package. The sweep *selects* a model; regression tests *pin* the
   selection so that an `mrds` or `Distance` upgrade cannot silently change
   which model wins. Snapshot the whole table, not just the winner — a change in
   the ranking of models 3 and 4 is early warning.
2. **The `g(0)` correction slot.** A value *and* its standard error, applied at
   the abundance step, with three rules: no default and never silently 1;
   propagate the variance or refuse the correction, since its CV routinely
   dominates the CV of abundance; and name availability and perception
   separately, so it is visible which have been applied. See section 5 of the
   architecture document for why it cannot be estimated from a NARWC extract.
3. **An MRDS backend**, conditional on data that actually carries double-observer
   structure — with a guard that errors when it does not, rather than fitting a
   single-observer model and reporting it as though perception bias were
   handled.
4. **The `dsm` handoff**, from `distsamp::segments_as_sf(segs, "midpoints")`.
5. **A vignette**, once there is something end-to-end to walk through.

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
