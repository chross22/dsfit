# Check and prepare distance data for fitting

Checks a table of detections and returns it in the form
[`mrds::ddf()`](https://rdrr.io/pkg/mrds/man/ddf.html) expects.
[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
calls it before fitting. It is exported so the checks can be run on
their own before a sweep.

## Usage

``` r
prepare_distance_data(data, truncation, left = NULL, breaks = NULL)
```

## Arguments

- data:

  Data frame, one row per detection. Needs a numeric `distance` column
  (exact) or numeric `distbegin` and `distend` columns (binned).
  Optional columns: `object` (detection ID) and `distance_source`.

- truncation:

  Number. Right truncation distance, in the units of the distances.
  Required.

- left:

  Number or `NULL` (default). Left truncation distance. `NULL` means no
  left truncation.

- breaks:

  Numeric vector or `NULL` (default). Bin cutpoints for binned data,
  strictly increasing. Derived from the bins when `NULL`. Ignored, with
  a warning, for exact data.

## Value

A list:

- `data`:

  The prepared data frame, ready for
  [`mrds::ddf()`](https://rdrr.io/pkg/mrds/man/ddf.html).

- `binned`:

  `TRUE` for binned data, `FALSE` for exact.

- `breaks`:

  The cutpoints used (binned), or `NULL` (exact).

- `n_dropped`:

  Number of input rows not in `data`.

- `dropped`:

  Named integer vector splitting the drops by reason: `no_distance`,
  `beyond_truncation`, `inside_left`. It shows whether truncation
  trimmed a tail or removed most of the survey.

## Details

Distances are either *exact* (one `distance` per detection) or *binned*
(an interval `distbegin` to `distend` per detection, as recorded by
strip or band surveys). The two are fitted with different likelihoods.

## How it works

The steps run in this order. The first refusal stops the function with
an error.

1.  Check the arguments. `data` must be a data frame. `truncation` must
    be one positive finite number. `left`, if given, must be one
    non-negative finite number smaller than `truncation`.

2.  Check the distance columns. Any of `distance`, `distbegin` and
    `distend` that is present must be numeric (an all-`NA` column is
    allowed) and must not be negative.

3.  Decide the distance type. Refuse a table with no distances, and
    refuse a table with both exact and binned distances.

4.  Check `distance_source`, if present, on the rows with an exact
    `distance`. Refuse `"break_off"` distances. Warn about `"circling"`
    distances. Only rows that will be fitted (within `left` to
    `truncation`) are counted.

5.  Count, by reason, the rows that will be dropped: no distance, beyond
    `truncation`, or inside `left`.

6.  For binned data:

    1.  Drop rows with a missing `distbegin` or `distend`.

    2.  Refuse any `distend` of `Inf`.

    3.  Drop bins that end beyond `truncation` or begin inside `left`.

    4.  Refuse if no rows remain.

    5.  Derive `breaks` from the bins, refusing bins that do not tile;
        or, if `breaks` is supplied, check that it contains every
        recorded bin edge.

    6.  Refuse if the highest break is not `truncation`.

    7.  Refuse if the lowest break is not `left` (or 0 when `left` is
        `NULL`).

    8.  Remove the `distance` column.

7.  For exact data:

    1.  Warn that `breaks` is ignored, if it was supplied.

    2.  Remove the `distbegin` and `distend` columns.

    3.  Drop rows with no distance, beyond `truncation`, or inside
        `left`.

    4.  Refuse if no rows remain.

8.  Add an `object` column numbering the rows if there is none. Refuse
    duplicated `object` values.

Rows with no distance are how a flatfile records a sample with no
detections, so dropping them is expected. Every dropped row is counted
in the return value.

## What it refuses and why

|  |  |
|----|----|
| Condition | Reason |
| Non-numeric distances | Text or factor codes would be compared to `truncation` as text or level numbers, keeping the wrong rows. |
| Negative distances | A perpendicular distance has no side. Use [`abs()`](https://rdrr.io/r/base/MathFun.html) if the sign records the side of the line. |
| Exact and binned distances together | The two likelihoods differ, so their AICs cannot be compared in one sweep. Split the table on its distance-source column and sweep each part. |
| `distance_source == "break_off"` | A break-off distance runs from where the aircraft left the line to the logged position. It is not perpendicular to the line. |
| `distend` of `Inf` | An open top bin cannot be fitted. Drop those rows or close the bin. This applies even if truncation would drop the bin. |
| Bins that do not tile (derived `breaks`) | If bins leave gaps or overlap, they do not define one set of cutpoints. |
| `breaks` missing a recorded bin edge | Goodness of fit is computed over `breaks`, so `breaks` must be the survey's own cutpoints. |
| Highest break not equal to `truncation` | A binned fit integrates over the bins. A gap between the top bin and `truncation` is effort with no bin. |
| Lowest break not equal to `left` (or 0) | The same gap at the near end. It occurs when the innermost bin had no detections, or when `left` falls inside a bin. Supply `breaks` to fix it. |
| No rows left after truncation | There is nothing to fit. |
| Duplicated `object` | `mrds` identifies detections by `object`. |

## What it warns about

- `distance_source == "circling"`: the animal's position while the
  aircraft circled, projected onto the line. It is a perpendicular
  distance, but measured minutes after the detection.
  `distsamp::detection_data()` excludes these rows unless
  `include_circling = TRUE`. They are kept, and the warning gives their
  count.

- `breaks` supplied with exact distances: it is ignored.

## See also

[`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)

## Examples

``` r
d <- data.frame(object = 1:5, distance = c(10, 50, 120, NA, 900))
prep <- prepare_distance_data(d, truncation = 400)
prep$n_dropped
#> [1] 2
prep$dropped
#>       no_distance beyond_truncation       inside_left 
#>                 1                 1                 0 
prep$data
#>   object distance
#> 1      1       10
#> 2      2       50
#> 3      3      120
```
