# Availability: the probability an animal was at the surface to be seen

Computes the probability that an animal was available to be detected
while it was within view, from its surfacing and diving intervals and
the time the survey platform kept it in view. This is the availability
component of `g(0)`, and it is **computed from external dive data, not
estimated from the survey**.

## Usage

``` r
availability(
  surface,
  dive,
  window,
  se_surface = NULL,
  se_dive = NULL,
  cov_surface_dive = 0,
  key = NULL
)
```

## Arguments

- surface:

  Mean surfacing interval, in seconds. A vector gives one result per
  element, which is how a per-month table is built.

- dive:

  Mean diving interval, in seconds. Recycled against `surface`.

- window:

  Time the animal is within view, in seconds. See
  [`view_window()`](https://camilleross.org/dsfit/reference/view_window.md)
  to derive it from platform geometry. Recycled against `surface`.

- se_surface, se_dive:

  Standard errors of `surface` and `dive`. Both are needed for a
  standard error on the result; either alone is an error, since
  propagating one source of variance and silently dropping the other
  understates the total.

- cov_surface_dive:

  Covariance between the two interval means. Zero by default, which is
  right when they are estimated from separate follows and optimistic
  when they are not.

- key:

  Optional labels — months, platforms, years — carried through to the
  result so the rows stay attached to what they apply to.

## Value

A tibble with one row per element of `surface`: `key`, `component`
(always `"availability"`), `value`, and `se`. That is the shape a `g(0)`
correction is assembled from, so availability and perception rows stack.

## Why this is a calculation and not an estimator

An animal submerged while the aircraft passes is missed at every
perpendicular distance equally. To first order availability is a pure
scale factor on \\g(x)\\: it does not dent the near-zero end of the
distance distribution, does not change its shape, and leaves no
signature for a likelihood to find. That is exactly why a mis-specified
`g(0)` shifts every candidate in a
[`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md)
by the same factor and leaves the ranking looking untouched.

So there is nothing in the distances to estimate it from, and none of
the inputs here come from the survey being corrected. `surface` and
`dive` come from focal follows, tagging, or drone observation; `window`
comes from the platform's geometry and speed. Both arrive with their own
uncertainty, which is the point: it is visible rather than absorbed into
a constant.

## The formula

Following Laake et al. (1997), for a mean surfacing interval \\E(s)\\, a
mean diving interval \\E(d)\\, and a window \\w\\ during which the
animal is in view:

\$\$a = \frac{E(s) + E(d)\left(1 - e^{-w/E(d)}\right)}{E(s) + E(d)}\$\$

The two limits are worth holding onto as a check on any number this
returns. As \\w \to 0\\ the window is instantaneous and \\a\\ becomes
\\E(s)/(E(s) + E(d))\\, the plain proportion of time spent at the
surface. As \\w \to \infty\\ the platform watches forever, every animal
surfaces eventually, and \\a \to 1\\. The middle term is the extra
chance that an animal which was down when the window opened comes up
before it closes.

## There is no single availability, and that is the finding

Ganley et al. (2019) measured this for right whales in Cape Cod Bay and
found availability varying by month from **0.27 in January to 0.91 in
April**, tracking the depth of the copepod layer the whales were feeding
on. Those measurements ship as
[ganley_availability](https://camilleross.org/dsfit/reference/ganley_availability.md).

A plausible-looking single figure is therefore the most dangerous input
this package accepts. A value near 0.83 is a real number for some month,
and wrong by a factor of three for others — and because it scales every
model equally, no amount of model selection reveals it. Compute one per
key and pass a vector, rather than picking a representative value.

Note that this is availability alone. Perception is a separate component
and a separate measurement, which Ganley et al. did not make —
estimating it needs a second observer team.
[`g0()`](https://camilleross.org/dsfit/reference/g0.md) will say so if
you leave it out.

## How the standard error is obtained

Availability is deterministic given `surface`, `dive` and `window`, so
the uncertainty comes from them. Given `se_surface` and `se_dive`, the
Jacobian is taken numerically and combined with their covariance matrix
— the same delta-method treatment
[`effect_estimates_ddf()`](https://camilleross.org/dsfit/reference/effect_estimates_ddf.md)
uses.

Where the raw focal follows are in hand, resampling them is better than
this: Ganley et al. (2019) and others bootstrap the follows and take the
standard deviation of the resulting availabilities, which does not
assume the interval means are jointly normal. Do that and pass the
result on as `se_surface` and `se_dive`, or bypass this and build the
rows directly.

Without `se_surface` and `se_dive` the standard error is `NA`,
deliberately: this function will not invent a precision for a number
whose CV routinely dominates the CV of abundance.

## References

Laake, J.L., Calambokidis, J., Osmek, S.D. and Rugh, D.J. (1997)
Probability of detecting harbor porpoise from aerial surveys: estimating
g(0). *The Journal of Wildlife Management* 61:63-75.
[doi:10.2307/3802415](https://doi.org/10.2307/3802415) The formula
implemented here.

Laake, J.L. and Borchers, D.L. (2004) Methods for incomplete detection
at distance zero. In *Advanced Distance Sampling*, pp. 108-189. Oxford
University Press. Why perception needs a second observer team, and
availability needs something outside the survey entirely.

Ganley, L.C., Brault, S. and Mayo, C.A. (2019) What we see is not what
there is: estimating North Atlantic right whale *Eubalaena glacialis*
local abundance. *Endangered Species Research* 38:101-113.
[doi:10.3354/esr00938](https://doi.org/10.3354/esr00938) Focal follows
and aircraft field of view applied to right whales, and the monthly
variation quoted above.

Roberts, J.J., Yack, T.M., Fujioka, E., Halpin, P.N., Baumgartner, M.F.
and others (2024) North Atlantic right whale density surface model for
the US Atlantic evaluated with passive acoustic monitoring. *Marine
Ecology Progress Series* 732:167-192.
[doi:10.3354/meps14547](https://doi.org/10.3354/meps14547) Corrects
perception and availability per platform, team and conditions across 11
institutions, which is the scale at which these corrections actually
vary.

## See also

[`view_window()`](https://camilleross.org/dsfit/reference/view_window.md)
for the window,
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
for why none of this belongs in the detection function fit.

## Examples

``` r
# An instantaneous window is the proportion of time at the surface
availability(surface = 60, dive = 240, window = 0)
#> # A tibble: 1 × 4
#>   key   component    value    se
#>   <chr> <chr>        <dbl> <dbl>
#> 1 NA    availability   0.2    NA

# A real window lifts it: some animals that were down come up in time
availability(surface = 60, dive = 240, window = 30)
#> # A tibble: 1 × 4
#>   key   component    value    se
#>   <chr> <chr>        <dbl> <dbl>
#> 1 NA    availability 0.294    NA

# One row per month, which is how it is actually used
availability(
  surface = c(45, 60, 90),
  dive    = c(300, 240, 150),
  window  = 24,
  key     = c("Feb", "Mar", "Apr")
)
#> # A tibble: 3 × 4
#>   key   component    value    se
#>   <chr> <chr>        <dbl> <dbl>
#> 1 Feb   availability 0.197    NA
#> 2 Mar   availability 0.276    NA
#> 3 Apr   availability 0.467    NA

# With uncertainty on the intervals
availability(surface = 60, dive = 240, window = 24,
             se_surface = 8, se_dive = 25)
#> # A tibble: 1 × 4
#>   key   component    value     se
#>   <chr> <chr>        <dbl>  <dbl>
#> 1 NA    availability 0.276 0.0297
```
