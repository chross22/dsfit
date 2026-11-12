# dsfit

<!-- badges: start -->
[![R-CMD-check](https://github.com/chross22/dsfit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chross22/dsfit/actions/workflows/R-CMD-check.yaml)
[![check-citations](https://github.com/chross22/dsfit/actions/workflows/check-citations.yaml/badge.svg)](https://github.com/chross22/dsfit/actions/workflows/check-citations.yaml)
<!-- badges: end -->

Fits and compares line-transect detection functions through `mrds`.

Input: one row per detection, with `distance` (exact) or `distbegin` and
`distend` (binned). The flatfile from
[distsamp](https://github.com/chross22/distsamp) has this shape.

## Install

```r
remotes::install_github("chross22/dsfit")
```

Plotting uses [fancyfx](https://github.com/chross22/fancyfx), optional.

## Pipeline

A basic analysis is five numbered scripts, installed with the package. Each
reads the previous step's output from `dsfit-output/`.

| Script | Step |
|---|---|
| `01-load-data.R` | Load the detections |
| `02-check-data.R` | Check the data; choose truncation and covariates |
| `03-fit-detection-functions.R` | Fit every candidate and rank them |
| `04-inspect-and-plot.R` | Choose a model and plot g(x) |
| `05-g0-correction.R` | Build the g(0) correction for the abundance step |

```r
scripts <- system.file("scripts", package = "dsfit")
file.copy(list.files(scripts, full.names = TRUE), ".")
for (f in sort(list.files(".", pattern = "^0[1-5].*\\.R$"))) source(f)
```

They run as shipped on simulated data. Edit the Settings block at the top of
each for real data. Step by step:
[docs/03-running-a-pipeline.md](docs/03-running-a-pipeline.md).

The sections below are the same steps as single function calls.

## Example data

```r
library(dsfit)

set.seed(1)
n <- 300
detections <- data.frame(
  object   = 1:n,
  beaufort = sample(0:4, n, replace = TRUE),
  glare    = sample(c("none", "some"), n, replace = TRUE),
  size     = rpois(n, 1) + 1
)
detections$distance <- abs(rnorm(n, 0, 150 - 15 * detections$beaufort))
```

Synthetic.

## Check

```r
detection_structure(detections)
diagnose_sweep(detections, model_set(c("hn", "hr")), truncation = 400)
```

Nothing is fitted.

## Fit

```r
sw <- sweep_models(detections, model_set(c("hn", "hr", "gamma")), truncation = 400)

sw
selection_table(sw)
sw$fits[["hn"]]     # the mrds fit
```

`truncation` is required. `left = ` sets left truncation.

## Model sets

```r
model_set(key = c("hn", "hr", "gamma"))
model_set(key = "hn", adjustment = "cos", order = 3)   # cos orders 2 and 3
model_set(key = "hn", formula = list(~1, ~beaufort))
```

## All covariate combinations

```r
covs   <- c("beaufort", "glare", "size")
combos <- unlist(lapply(seq_along(covs), function(k) combn(covs, k, simplify = FALSE)),
                 recursive = FALSE)
forms  <- c(list(~1), lapply(combos, reformulate))

sw_cov <- sweep_models(detections, model_set(c("hn", "hr"), formula = forms),
                       truncation = 400)
selection_table(sw_cov)
```

## Plot

```r
fit <- sw$fits[["hn"]]

fancyfx::plotEffects(fit, dat = detections, var = "distance")

est <- effect_estimates_ddf(fit)    # .x, .estimate, .lower, .upper
plot(est$.x, est$.estimate, type = "l", ylim = c(0, 1))

effect_estimates_ddf(sw_cov$fits[["hn beaufort"]], at = list(beaufort = 4))
```

## Binned distances

```r
binned <- transform(detections, distbegin = floor(distance / 100) * 100)
binned <- transform(binned, distend = distbegin + 100, distance = NULL)

sweep_models(binned[binned$distend <= 400, ], model_set("hn"), truncation = 400)
```

## Availability and g(0)

```r
avail <- availability(surface = 60, dive = 240, window = 24,
                      se_surface = 8, se_dive = 25)

perception <- data.frame(key = NA, component = "perception", value = 0.68,
                         se = 0.09, source = "borrowed")
g0(avail, perception)
```

Not applied to the fit. Divide abundance by it.

## Selection table

| Column | Meaning |
|---|---|
| `delta_aic` | AIC difference from the lowest |
| `p`, `p_cv` | Average detection probability and its CV |
| `esw` | Effective strip half-width, `p * (truncation - left)` |
| `cvm_p` / `chisq_p` | Goodness of fit, exact / binned |
| `converged` | `FALSE` for a failed fit |

## Other functions

```r
view_window_aerial(trackline = 51.22, speed = 51.4, slope = 0.03)   # time in view
view_window_vessel(radius = 6000, speed = 5.1)
prepare_distance_data(detections, truncation = 400)                 # the input guards
```

## Documentation

| Topic | File |
|---|---|
| Running a pipeline | [docs/03-running-a-pipeline.md](docs/03-running-a-pipeline.md) |
| Code guide | [docs/02-code-guide.md](docs/02-code-guide.md) ([PDF](docs/02-code-guide.pdf)) |
| Design and guards | [docs/01-plan.md](docs/01-plan.md) |

```r
vignette("dsfit", package = "dsfit")
vignette("availability", package = "dsfit")
```

## Citing dsfit

```r
citation("dsfit")
```

Returns the package, `mrds` (Miller et al. 2019) and the gamma key (Becker and
Quang 2009). Report the truncation and the model set.
