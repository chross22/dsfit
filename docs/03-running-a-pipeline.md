# Running a basic dsfit pipeline

*C. H. Ross. 2026-11-12. Kept in the package repository; not part of the built package.*

## 1. Scope

This chapter describes how a basic analysis is run with dsfit, from a table of
detections to a ranked set of detection functions, a plotted detection
function, and a g(0) correction ready for the abundance step. Each step is a
numbered script in `inst/scripts/`. The scripts run as shipped on simulated
data; each has a Settings block at the top for real data.

Abundance and density are not computed by dsfit. The chosen detection
function and the g(0) correction are handed to the abundance step, which is
outside this package (`Distance::dht()`, `dsm`, or an analysis repository).

Terms used throughout:

| Term | Meaning |
|---|---|
| detection | One animal or group seen, with its perpendicular distance from the trackline. |
| detection function g(x) | The probability that an animal at perpendicular distance x is detected. |
| truncation | The distance beyond which detections are discarded. One value per sweep. |
| `left` | The distance below which detections are discarded, for a platform that cannot see directly beneath itself. |
| candidate | One model to fit: a key function, optional adjustment terms, and an optional covariate formula. |
| sweep | Every candidate fitted to the same data at the same truncation. |
| `p` | Average detection probability within the truncation distance. |
| esw | Effective strip half-width, `p * (truncation - left)`: the width of a strip in which every animal would have given the same number of detections. Abundance is divided by it. |
| g(0) | The probability of detecting an animal on the trackline. Every fit assumes 1. |

## 2. The scripts

| Script | Does | Reads | Writes |
|---|---|---|---|
| `01-load-data.R` | Loads the detections | your data | `01-detections.rds` |
| `02-check-data.R` | Checks the data; chooses truncation and covariates | `01` | `02-settings.rds` |
| `03-fit-detection-functions.R` | Fits every candidate; ranks them | `01`, `02` | `03-sweep.rds`, `03-selection-table.csv` |
| `04-inspect-and-plot.R` | Chooses a model; plots g(x) | `01`, `03` | `04-chosen-model.rds`, `04-detection-function.png` |
| `05-g0-correction.R` | Builds the g(0) correction | dive data | `05-g0.rds` |

All outputs are written to `dsfit-output/` in the working directory. The
scripts are copied out of the installed package and run in order:

```r
scripts <- system.file("scripts", package = "dsfit")
file.copy(list.files(scripts, full.names = TRUE), ".")
for (f in sort(list.files(".", pattern = "^0[1-5].*\\.R$"))) source(f)
```

Each script can also be run alone once the scripts before it have run.

## 3. Step 1: load the detections

The input is one row per detection. It must have either:

- `distance`: exact perpendicular distances, or
- `distbegin` and `distend`: binned distances, one interval per row.

A table with both is refused, because the two are fitted with different
likelihoods and their AICs cannot be compared. Detection covariates (sea
state, glare, group size) are further columns. An `object` column is added if
absent; its values must be unique.

Rows with no distance are allowed. A flatfile uses them for survey samples
with no detection, and they are dropped and counted before fitting.

The flatfile from `distsamp::detection_data()` is in this shape. The
`read.csv()` and `distsamp` lines at the top of `01-load-data.R` replace the
simulated example.

## 4. Step 2: check the data

Two functions are run, and nothing is fitted.

`detection_structure()` reports what the table can support: exact or binned
distances and how many, which columns could be covariates, and whether
perception bias or group size can be handled.

`diagnose_sweep()` runs the same input checks as the sweep, and reports:

- how many rows truncation drops, and why;
- whether there are enough detections (about 60-80 is usually suggested);
- whether the nearest detections start well away from the trackline (a blind
  spot, handled by `left` or the gamma key);
- whether each covariate exists, varies, and has no missing values.

The truncation is set here. By default the script suggests the distance that
keeps the nearest 95% of detections, following Buckland et al. (2001). For
binned data the truncation must equal a bin edge. The truncation, `left` and
the covariate list are saved for step 3.

## 5. Step 3: fit and rank

`model_set()` builds the candidates. The script crosses the key functions
with every combination of the covariates from step 2, including none: three
covariates give 8 formulas, and two keys give 16 candidates. Adjustment terms
can be added; they are fitted without covariates, as `Distance::ds()` fits
them.

`sweep_models()` fits every candidate through `mrds::ddf()` with the same
truncation, so their AICs are comparable. A candidate that fails or does not
converge is recorded with its reason rather than stopping the sweep.

`selection_table()` returns the candidates sorted by AIC. Its columns:

| Column | Meaning |
|---|---|
| `delta_aic` | AIC difference from the best model |
| `p`, `p_cv` | Average detection probability, and its coefficient of variation |
| `esw` | Effective strip half-width |
| `cvm_p` | Cramér-von Mises goodness-of-fit p-value (exact distances) |
| `chisq_p` | Chi-square goodness-of-fit p-value (binned distances) |
| `converged` | `FALSE` for a candidate that failed |

AIC alone is not enough to choose. Two models within 2 AIC can have
different `esw`, and abundance changes by the same proportion. A small
goodness-of-fit p-value is evidence against a model, whatever its AIC.

## 6. Step 4: choose and plot

The script takes the lowest-AIC model unless one is named in its Settings
block. It prints the models within 2 AIC for comparison.

`effect_estimates_ddf()` returns g(x) on a grid of distances with a
confidence band. For a covariate model, the curve is averaged over the
covariate values observed; `at = list(beaufort = 4)` fixes a covariate
instead. The script draws the curve with base graphics and marks the observed
distances along the axis. `plot(fit)` gives the `mrds` plot, and
`fancyfx::plotEffects(fit, dat = detections, var = "distance")` gives the
curve with a histogram rug.

## 7. Step 5: build the g(0) correction

Every fit assumes g(0) = 1. For a diving animal this is false, and abundance
is biased low by the missing factor. The correction is built from separate
components:

| Component | Question | Source |
|---|---|---|
| availability | Was the animal at the surface while in view? | `availability()`, from dive intervals and the time in view |
| perception | Was an available animal seen? | A double-observer study, often from another survey |

The time in view comes from the platform geometry: `view_window_aerial()` for
an aircraft side window, `view_window_vessel()` for a ship's forward search
arc, `view_window()` for a circular view. `g0()` multiplies the components
and combines their uncertainty (squared CVs add). A component left out is
warned about, and each component's source is printed.

The correction is saved, not applied. At the abundance step, abundance is
divided by the g(0) value for the matching key, and the g(0) CV is added in
quadrature to the abundance CV.

## 8. Limitations

- The shipped scripts run on simulated detections and invented dive
  intervals. Their output demonstrates the pipeline and is not a result.
- The default truncation (95% of detections) is a starting point, not a
  decision. It is reviewed against the diagnosis in step 2.
- The default model choice (lowest AIC) is reviewed against goodness of fit
  and `esw` in step 4.
- dsfit fits conventional and covariate detection functions only. Double-
  observer (mark-recapture) models are not fitted.

## 9. Verification performed

| Check | Result |
|---|---|
| Scripts 01-05 run in order from an empty directory against the current package | All seven output files written |
| Step 2 on the simulated data | Truncation 240; 284 of 320 rows kept (20 with no distance, 16 beyond truncation) |
| Step 3 on the simulated data | 16 of 16 candidates fitted; best `hn beaufort` |

## 10. References

Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L. and
Thomas, L. (2001) *Introduction to Distance Sampling: Estimating Abundance of
Biological Populations.* Oxford University Press. Truncation, sample size and
model selection.

Laake, J.L., Calambokidis, J., Osmek, S.D. and Rugh, D.J. (1997) Probability
of detecting harbor porpoise from aerial surveys: estimating g(0). *The
Journal of Wildlife Management* 61:63–75. doi:10.2307/3802415. The
availability formula.

Miller, D.L., Rexstad, E., Thomas, L., Marshall, L. and Laake, J.L. (2019)
Distance sampling in R. *Journal of Statistical Software* 89(1):1–28.
doi:10.18637/jss.v089.i01. The `mrds` fitting interface.
