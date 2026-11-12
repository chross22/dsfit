# Check a sweep's inputs before fitting

Runs the checks a sweep would run, without calling
[`mrds::ddf()`](https://rdrr.io/pkg/mrds/man/ddf.html), and prints a
report of problems found. Use it before
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
to catch a misconfiguration in seconds, instead of after a sweep that
errors inside `mrds` or returns a table that looks valid but is not.

## Usage

``` r
diagnose_sweep(data, models = NULL, truncation, left = NULL, breaks = NULL)
```

## Arguments

- data:

  Data frame, one row per detection, in the form
  [`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
  takes.

- models:

  A
  [`model_set()`](https://camilleross.org/dsfit/reference/model_set.md),
  or `NULL` (default) for the default set.

- truncation:

  Number. Right truncation distance, in the units of the distances.
  Required, as for a sweep.

- left:

  Number or `NULL` (default). Left truncation distance.

- breaks:

  Numeric vector or `NULL` (default). Bin cutpoints for binned data;
  derived from the bins when `NULL`.

## Value

Invisibly, a list with the results reached before any stop:

- `structure`:

  The
  [`detection_structure()`](https://camilleross.org/dsfit/reference/detection_structure.md)
  result, or `NULL` if the data could not be read.

- `prepared`:

  The
  [`prepare_distance_data()`](https://camilleross.org/dsfit/reference/prepare_distance_data.md)
  result, or `NULL` if it failed. Absent if the report stopped earlier.

- `models`:

  The model set checked. Absent if the report stopped earlier.

The report itself is printed to the console.

## Details

The function only reports. It never changes the data or the model set.

## How it works

Each finding is printed as one line marked `ok`, `note`, `WARN` or
`FAIL`. A `note` is a fact about how the survey was flown (for example,
a single observer team); it does not count as a problem. `WARN` and
`FAIL` lines count as problems. A failure marked "stops" below ends the
report early.

1.  Refuse to run if `truncation` is missing (an error, not a report
    line).

2.  Toolchain: check `mrds` is installed (`WARN` if not).

3.  Data: read the table with
    [`detection_structure()`](https://camilleross.org/dsfit/reference/detection_structure.md).

    - `FAIL` and stop if it cannot be read.

    - `FAIL` if it has no distances; otherwise report the count and
      type.

    - `WARN` for each pair of columns that name the same variable.

    - `WARN` for each analysis whose structure is incomplete
      (`supported` is `NA`), such as a partial double-observer setup.

    - `note` if perception bias cannot be estimated or there is no
      `size` column.

4.  Truncation: run
    [`prepare_distance_data()`](https://camilleross.org/dsfit/reference/prepare_distance_data.md).

    - Its warnings are printed as `WARN` lines. Its errors are printed
      as a `FAIL` line and stop the report.

    - Report how many rows were kept, and how many were dropped for each
      reason.

    - `WARN` if truncation removes more than 50% of the rows that have a
      distance.

5.  Sample size: `FAIL` below 40 detections after truncation, `WARN`
    below 60, `ok` otherwise (Buckland et al. 2001).

6.  Model set: build the default
    [`model_set()`](https://camilleross.org/dsfit/reference/model_set.md)
    or check the one given. `FAIL` and stop if that errors; otherwise
    list the candidates.

7.  Blind spot: see the section below.

8.  Covariates: `FAIL` if a formula names a column not in `data`.
    Otherwise `WARN` for a covariate with only one value, and for one
    with missing values (mrds drops those rows).

9.  Print a reminder that g(0) = 1 is assumed unless corrected at the
    abundance step. g(0) is the detection probability on the trackline.
    An error in it scales every candidate equally, so no check here can
    detect it. See
    [`g0()`](https://camilleross.org/dsfit/reference/g0.md).

10. Print a final line: no problems found, or problems above.

The constant and missing-value covariate checks use the rows that will
be fitted, after truncation, since mrds never sees the rest.

## Blind-spot check

A *blind spot* is a strip next to the trackline where the platform
cannot see, such as beneath an aircraft. It shows as no detections close
to the line.

The check fires when the nearest detection in the input data is more
than 2% of `truncation` from the line. Every dataset has a smallest
distance, so a nonzero nearest distance alone is not evidence of a blind
spot. The 2% threshold is a heuristic.

There are two treatments, and they are alternatives: left truncation
(`left`), or the gamma key, which models detection rising from zero at
the line. The check reports:

|  |  |
|----|----|
| Situation | Result |
| Gamma key in the set and `left` set (blind spot or not) | `WARN`: the blind spot is counted twice. |
| Blind spot, `left` set, no gamma key | `ok` |
| Blind spot, abrupt onset | `WARN`: a hard edge, which `left` removes and no key function fits, gamma included. |
| Blind spot, gradual onset, gamma key in set | `ok` |
| Blind spot, too few detections (under 40) to classify the onset, gamma key in set | `note` |
| Blind spot, neither `left` nor gamma | `WARN` |

*Abrupt onset* means detections start at close to their peak rate just
beyond the empty strip. The range from the nearest detection to
`truncation` is cut into 20 equal bins. The onset is abrupt when the
first bin's count is more than 60% of the largest bin's count. A
geometric edge (a window frame, a hull) gives this shape. A gradual
decline in detectability toward the line instead gives a rising first
few bins, which is what the gamma key fits. The test needs at least 40
detections; with fewer, the onset is not classed as abrupt.

## References

Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers,
D.L. and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford
University Press. Source of the 40 and 60-80 detection guidelines and of
trimming roughly the outer 5% of distances.

## See also

[`detection_structure()`](https://camilleross.org/dsfit/reference/detection_structure.md),
which reports what the data could support in principle.
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md),
which this is meant to run before.

## Examples

``` r
set.seed(1)
d <- data.frame(
  object = 1:300,
  distance = c(abs(rnorm(295, 0, 800)), rep(NA, 5)),
  beaufort = sample(0:4, 300, replace = TRUE)
)
diagnose_sweep(d, model_set(c("hn", "hr")), truncation = 2000)
#> dsfit sweep diagnosis
#> 
#> == Toolchain ==
#>   ok    mrds 3.0.1 is installed
#> 
#> == Data ==
#>   ok    300 rows, 295 exact distances
#>   note  estimate perception bias: no double-observer structure: no `observer` or `detected` column
#>   note  carry group size to abundance: no `size` column; abundance would count groups, not individuals
#> 
#> == Truncation ==
#>   ok    293 of 300 rows kept at truncation 2000
#>         dropped: 5 with no distance, 2 beyond the truncation
#>   ok    293 detections is at or above the 60-80 usually suggested
#> 
#> == Model set ==
#>   ok    2 candidates: hn, hr
#> 
#> == What is still assumed ==
#>   g(0) = 1, unless a correction is applied at the abundance step.
#>   Nothing here can detect a wrong one: it scales every candidate
#>   equally, so the ranking looks untouched. See ?g0.
#> 
#> No problems found. Nothing was fitted.
```
