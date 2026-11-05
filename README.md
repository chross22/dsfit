# dsfit

<!-- badges: start -->
[![R-CMD-check](https://github.com/chross22/dsfit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chross22/dsfit/actions/workflows/R-CMD-check.yaml)
[![check-citations](https://github.com/chross22/dsfit/actions/workflows/check-citations.yaml/badge.svg)](https://github.com/chross22/dsfit/actions/workflows/check-citations.yaml)
<!-- badges: end -->

Fits a set of line-transect detection functions through `mrds` and ranks them
in one table: AIC, average detection probability `p`, its CV, effective strip
half-width `esw`, and a goodness-of-fit p-value.

Input is one row per detection: `distance` (exact), or `distbegin` and
`distend` (binned). The flatfile from
[`distsamp`](https://github.com/chross22/distsamp) has this shape.

## Install

```r
# install.packages("remotes")
remotes::install_github("chross22/dsfit")
```

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

## Check the data

```r
detection_structure(detections)

diagnose_sweep(detections, model_set(c("hn", "hr")), truncation = 400)
```

## Fit a sweep

```r
sw <- sweep_models(detections, model_set(c("hn", "hr", "gamma")), truncation = 400)

sw
selection_table(sw)
```

## Build a model set

```r
model_set()                                             # hn, hr
model_set(key = c("hn", "hr", "gamma"))                 # add the gamma key
model_set(key = "hn", adjustment = "cos", order = 2:3)  # adjustment terms
model_set(key = "hn", formula = list(~1, ~beaufort))    # covariates
```

## Sweep every covariate combination

```r
covs <- c("beaufort", "glare", "size")

combos <- unlist(
  lapply(seq_along(covs), function(k) combn(covs, k, simplify = FALSE)),
  recursive = FALSE
)
formulas <- c(list(~1), lapply(combos, reformulate))

ms <- model_set(key = c("hn", "hr"), formula = formulas)   # 2 keys x 8 formulas
sw_cov <- sweep_models(detections, ms, truncation = 400)

selection_table(sw_cov)[, c("model_id", "delta_aic", "p", "p_cv", "esw", "cvm_p")]
```

## Plot a detection function

```r
fit <- sw$fits[["hn"]]

# g(x) with a delta-method band and a histogram rug (needs fancyfx)
fancyfx::plotEffects(fit, dat = detections, var = "distance")

# The same curve, base graphics
est <- effect_estimates_ddf(fit)
plot(est$.x, est$.estimate, type = "l", ylim = c(0, 1),
     xlab = "Distance (m)", ylab = "g(x)")
lines(est$.x, est$.lower, lty = 2)
lines(est$.x, est$.upper, lty = 2)

# A covariate model, with one covariate fixed and the rest averaged
est_bf <- effect_estimates_ddf(sw_cov$fits[["hn beaufort"]], at = list(beaufort = 4))

# mrds's own plot
plot(fit)
```

## Binned distances

```r
binned <- detections
binned$distbegin <- floor(binned$distance / 100) * 100
binned$distend   <- binned$distbegin + 100
binned$distance  <- NULL
binned <- binned[binned$distend <= 400, ]

sweep_models(binned, model_set("hn"), truncation = 400)
```

## Left truncation

```r
sweep_models(detections, model_set(c("hn", "hr")), truncation = 400, left = 20)
```

## Availability and g(0)

```r
avail <- availability(surface = 60, dive = 240, window = 24,
                      se_surface = 8, se_dive = 25)

view_window(radius = 300, speed = 50)                     # circular view
view_window_aerial(trackline = 51.22, speed = 51.4,       # aircraft side window
                   slope = 0.03, distance = c(0, 500, 1000))
view_window_vessel(radius = 6000, speed = 5.1)            # ship, forward arc

perception <- data.frame(key = NA, component = "perception", value = 0.68,
                         se = 0.09, source = "borrowed from another survey")
g0(avail, perception)
```

`g0()` builds the correction. It is applied at the abundance step, not here.

## Rules enforced

| Condition | Behaviour |
|---|---|
| `truncation` missing | error |
| Exact and binned distances in one table | error |
| Open top bin (`distend = Inf`) | error |
| Bins that do not tile | error |
| Bins not spanning `left` (or 0) to `truncation` | error |
| `distance_source == "break_off"` | error |
| `distance_source == "circling"` | warning |
| Gamma key with `left` | warning |
| Uniform key with no adjustment | warning, dropped |
| Uniform key with covariates | warning, dropped |
| Candidate that fails or does not converge | recorded, `converged = FALSE` |

## Outputs

| Column | Meaning |
|---|---|
| `aic`, `delta_aic` | AIC, and difference from the lowest |
| `p`, `p_se`, `p_cv` | Average detection probability, SE, CV |
| `esw` | Effective strip half-width, `p * (truncation - left)` |
| `cvm_p` | Cramér-von Mises p-value (exact distances) |
| `chisq_p` | Chi-square p-value (binned distances) |
| `converged` | Whether the candidate fitted |

All fits assume `g(0) = 1`.

## Documentation

| | |
|---|---|
| [Getting started](vignettes/dsfit.Rmd) | model set, sweep, reading the table |
| [Availability and g(0)](vignettes/availability.Rmd) | availability, view windows, `g0()` |
| [Plan](docs/01-plan.md) | architecture and guards |
| [Code guide](docs/02-code-guide.md) ([PDF](docs/02-code-guide.pdf)) | the code, chunk by chunk |

## Citation

```r
citation("dsfit")
```

## References

Becker, E.A. and Quang, P.X. (2009) *JABES* 14:207–223.
<https://doi.org/10.1198/jabes.2009.0013> — gamma key.

Buckland, S.T. et al. (2001) *Introduction to Distance Sampling.* Oxford
University Press.

Ganley, L.C., Brault, S. and Mayo, C.A. (2019) *Endangered Species Research*
38:101–113. <https://doi.org/10.3354/esr00938> — availability, aerial view
window.

Laake, J.L., Calambokidis, J., Osmek, S.D. and Rugh, D.J. (1997) *Journal of
Wildlife Management* 61:63–75. <https://doi.org/10.2307/3802415> —
availability equation.

Miller, D.L., Rexstad, E., Thomas, L., Marshall, L. and Laake, J.L. (2019)
*Journal of Statistical Software* 89(1):1–28.
<https://doi.org/10.18637/jss.v089.i01> — `mrds`.
