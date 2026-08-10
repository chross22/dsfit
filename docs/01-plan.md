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
| **`dsfit`** | the model-selection sweep, gamma support, GOF and `p̄`, `g(0)` | **package** |
| DSM layer | covariates, density surface models, the `dsm` handoff | separate package |
| analysis repo | which years, which truncation, which covariates, the report | `targets` + `renv` |

The `dsm` handoff was listed here until it was looked at properly. It needs both
`distsamp`'s segment structures and `dsm`'s, and this package depends on
neither; the covariate stack behind a density surface model is heavy (`terra`,
`ncdf4`, `marmap`), which is the same argument `distsamp`'s next-steps already
makes for covariates and DSMs being their own package. Building an adapter here
would drag two dependency trees in for one function. `dsfit` stops at the
detection function and what it implies.

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
and a `print()` method. 274 tests, six of them snapshots over the selection
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

## 4. The `g(0)` correction

### The `g(0)` correction slot — built

`availability()` computes the availability component from dive intervals and
the platform's time in view, and `g0()` assembles the components: a **keyed
table**, not a value, with rows of `component` x `value` x `se`, multiplied
together and the variance propagated by delta method.

Keyed, because a scalar is the error. Ganley et al. (2019) measured right whale
availability in Cape Cod Bay varying by month from 0.27 in January to 0.91 in
April; Roberts et al. (2024) correct per platform, team and conditions. A figure
like 0.83 sits inside that range, which is what makes it dangerous — it is one
component's value under one set of conditions, and using it as `g(0)`
understates the correction by more than a factor of two when availability is
~0.5 and perception ~0.7.

`g0()` builds an object and hands it off rather than applying anything, because
this package fits detection functions and does not compute abundance — the
correction belongs one layer out. What it does from here is make the correction
impossible to get wrong quietly. The three rules are enforced rather than
described: no default and never silently 1, with any absent component named in
a warning and again in the printout; propagate the variance or refuse, since
`g(0)`'s CV routinely dominates the CV of abundance; and name the components
separately, so it stays visible which have been applied. A component may also
record a `source`, so a borrowed value cannot pass for a local one.

Under independence the delta method gives the rule worth remembering: squared
CVs add, `CV(g0)² = Σ CV(xᵢ)²`, so the least precise component sets the floor.
Components keyed on different things — availability by month, perception by
year — are refused rather than joined, because the honest answer there is a
value per month-year.

Ganley's measurements ship as `ganley_availability`, and their annual detection
functions as `ganley_detection` — **cited data, not defaults**. Note that
`ganley_detection`'s 0.43–0.87 figures are `p`, the detection function's own
average probability, which `selection_table()` already reports; they are not
perception, which Ganley et al. did not estimate.

### Perception, and why there is no estimator here

True double-observer data is not the norm in NARWC datasets. Whether perception
is estimable is a property of the contributing survey programme rather than of
the archive — AMAPPS ran two independent teams on both platforms, the Cape Cod
Bay programme did not — and the programmes that did are the exception. An MRDS
backend would be a large piece of machinery, needing maintenance against `mrds`
forever, aimed at a case that mostly does not arrive.

That is not a gap. For nearly every dataset this layer will see, perception has
two honest treatments and `g0()` implements both: left at 1 with the assumption
printed, or borrowed from a programme that could measure it, with `source`
recording that it is on loan. Roberts et al. (2024) took the second route for
ten of their eleven institutions and said so in print. That is the state of the
art, not a shortcut around it.

## 5. What is left

In the order it is worth doing.

1. ~~**A structure detector.**~~ Built, as `detection_structure()`: it reads a
   set of detections and reports what they can and cannot support, with the
   reason. Three verdicts, `can`, `cannot` and `partly` — the last for an
   `observer` column with no `detected` indicator, which looks like
   double-observer data and is not. It reports; `prepare_distance_data()`
   enforces.
2. ~~**A vignette.**~~ Built, as `vignettes/dsfit.Rmd`: one simulated aerial
   survey walked from a table of distances to a `g(0)` correction.

   Its spine is a live demonstration rather than an assertion. The survey has a
   blind spot beneath the aircraft, so the first sweep's AIC winner **fails
   goodness-of-fit, along with every other candidate** — because no key function
   reproduces a hard geometric cutoff. Left truncating fixes it, and the point
   lands on its own: the fix was a statement about the aircraft, not a better
   model, and the AIC column had no way to say so.

The **`dsm` handoff** has moved out of this package — see section 1.

Nothing else is planned. What would come next is use: running this against a
real NARWC extract and finding out which of its guards fire.

## 6. Deliberately out of scope

**Estimating `g(0)`.** Not a gap to be filled later — it cannot be estimated
from the data this layer receives. Three distinct things present as `g(0) < 1`:
the geometric blind spot, which is a property of the aircraft and is handled by
left truncation or a gamma key; availability, which needs dive data from tagging
studies; and perception, which needs a double-observer protocol. Correcting for
one while believing you have corrected for another is how these estimates go
wrong by a factor rather than a percentage.

**Wrapping `Distance::ds()`.** It is mature and well documented. A wrapper hides
its options and rots against it.
