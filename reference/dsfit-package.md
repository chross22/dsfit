# dsfit: Fit and Compare Detection Functions for Line-Transect Survey Data

Fits a set of candidate detection functions to line-transect distance
data and compares them on the quantities that matter for abundance, not
on AIC alone. Models are fitted through 'mrds' so that every candidate
comes off the same likelihood machinery and the gamma key function is
available alongside half-normal, hazard-rate and uniform. Enforces the
conditions under which a comparison is meaningful: truncation must be
held constant within a model set, and binned and exact distances cannot
be ranked together. Reports effective strip half-width and its
coefficient of variation next to each fit, since two models within two
AIC of one another can imply materially different abundance. Takes the
flatfile produced by 'distsamp' or any equivalent table.

## Main functions

- [`model_set()`](https://camilleross.org/dsfit/reference/model_set.md):

  Builds the grid of candidate detection functions.

- [`diagnose_sweep()`](https://camilleross.org/dsfit/reference/diagnose_sweep.md):

  Checks data and a model set for common problems before anything is
  fitted.

- [`prepare_distance_data()`](https://camilleross.org/dsfit/reference/prepare_distance_data.md):

  Validates distances and applies the truncation. Called by
  [`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md);
  also usable on its own.

- [`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md):

  Fits every candidate with
  [`mrds::ddf()`](https://rdrr.io/pkg/mrds/man/ddf.html) and builds a
  selection table.

- [`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md):

  Returns the ranked table from a sweep.

- [`detection_structure()`](https://camilleross.org/dsfit/reference/detection_structure.md):

  Reports what a table of detections can support, such as exact or
  binned distances and double-observer data.

- [`availability()`](https://camilleross.org/dsfit/reference/availability.md),
  [`g0()`](https://camilleross.org/dsfit/reference/g0.md):

  Build a `g(0)` correction from external data, for use at the abundance
  step.

## Where the package sits

`dsfit` is the middle of three layers. `distsamp` turns NARWC-format
survey data into effort segments and detection distances. `dsfit` fits
and compares detection functions. A separate analysis repository holds
the per-run choices: years, truncation, covariates. Logic that can be
tested lives in the package; project-specific choices do not.

## Why models are fitted through mrds

`Distance::ds()` is a good tool for fitting one model interactively, and
`dsfit` does not replace it. `ds()` offers only the `hn`, `hr` and
`unif` keys. The gamma key, which suits a platform that cannot see
directly beneath itself, is available only in `mrds`. Fitting every
candidate through [`mrds::ddf()`](https://rdrr.io/pkg/mrds/man/ddf.html)
also means every AIC in a selection table comes from the same likelihood
code, with the same handling of truncation and monotonicity.

## g(0) is not estimated

`g(0)` is the probability of detecting an animal on the track-line.
Every fit assumes `g(0) = 1`. A wrong `g(0)` scales every candidate by
the same factor, so the model ranking is unchanged while density is
biased. `g(0)` cannot be estimated from sighting distances alone.
Availability needs dive data, and perception needs a double-observer
design that some NARWC programmes recorded and many did not;
[`detection_structure()`](https://camilleross.org/dsfit/reference/detection_structure.md)
reports which a given extract has. Apply the correction at the abundance
step, from external sources, and propagate its standard error.
[`availability()`](https://camilleross.org/dsfit/reference/availability.md)
and [`g0()`](https://camilleross.org/dsfit/reference/g0.md) build that
correction.

## References

Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers,
D.L. and Thomas, L. (2001) *Introduction to Distance Sampling:
Estimating Abundance of Biological Populations.* Oxford University
Press, New York, NY.

Miller, D.L., Rexstad, E., Thomas, L., Marshall, L. and Laake, J.L.
(2019) Distance sampling in R. *Journal of Statistical Software*
89(1):1-28.
[doi:10.18637/jss.v089.i01](https://doi.org/10.18637/jss.v089.i01)

Laake, J.L. and Borchers, D.L. (2004) Methods for incomplete detection
at distance zero. In S.T. Buckland et al. (eds) *Advanced Distance
Sampling*, pp. 108-189. Oxford University Press.

Marques, F.F.C. and Buckland, S.T. (2004) Covariate models for the
detection function. In S.T. Buckland et al. (eds) *Advanced Distance
Sampling*, pp. 31-47. Oxford University Press.

## See also

Useful links:

- <https://github.com/chross22/dsfit>

- <https://camilleross.org/dsfit/>

- Report bugs at <https://github.com/chross22/dsfit/issues>

## Author

**Maintainer**: Camille Ross <camille.ross@maine.edu>
([ORCID](https://orcid.org/0000-0002-1428-2294))

Authors:

- Camille Ross <camille.ross@maine.edu>
  ([ORCID](https://orcid.org/0000-0002-1428-2294))
