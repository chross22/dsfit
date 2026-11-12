#' @keywords internal
#'
#' @section Main functions:
#' \describe{
#'   \item{[model_set()]}{Builds the grid of candidate detection functions.}
#'   \item{[diagnose_sweep()]}{Checks data and a model set for common problems
#'     before anything is fitted.}
#'   \item{[prepare_distance_data()]}{Validates distances and applies the
#'     truncation. Called by [sweep_models()]; also usable on its own.}
#'   \item{[sweep_models()]}{Fits every candidate with `mrds::ddf()` and builds
#'     a selection table.}
#'   \item{[selection_table()]}{Returns the ranked table from a sweep.}
#'   \item{[detection_structure()]}{Reports what a table of detections can
#'     support, such as exact or binned distances and double-observer data.}
#'   \item{[availability()], [g0()]}{Build a `g(0)` correction from external
#'     data, for use at the abundance step.}
#' }
#'
#' @section Where the package sits:
#' `dsfit` is the middle of three layers. `distsamp` turns NARWC-format survey
#' data into effort segments and detection distances. `dsfit` fits and compares
#' detection functions. A separate analysis repository holds the per-run
#' choices: years, truncation, covariates. Logic that can be tested lives in
#' the package; project-specific choices do not.
#'
#' @section Why models are fitted through mrds:
#' `Distance::ds()` is a good tool for fitting one model interactively, and
#' `dsfit` does not replace it. `ds()` offers only the `hn`, `hr` and `unif`
#' keys. The gamma key, which suits a platform that cannot see directly beneath
#' itself, is available only in `mrds`. Fitting every candidate through
#' `mrds::ddf()` also means every AIC in a selection table comes from the same
#' likelihood code, with the same handling of truncation and monotonicity.
#'
#' @section g(0) is not estimated:
#' `g(0)` is the probability of detecting an animal on the track-line. Every
#' fit assumes `g(0) = 1`. A wrong `g(0)` scales every candidate by the same
#' factor, so the model ranking is unchanged while density is biased.
#' `g(0)` cannot be estimated from sighting distances alone. Availability needs
#' dive data, and perception needs a double-observer design that some NARWC
#' programmes recorded and many did not; [detection_structure()] reports which
#' a given extract has. Apply the correction at the abundance step, from
#' external sources, and propagate its standard error. [availability()] and
#' [g0()] build that correction.
#'
#' @references
#' Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L.
#' and Thomas, L. (2001) *Introduction to Distance Sampling: Estimating
#' Abundance of Biological Populations.* Oxford University Press, New York, NY.
#'
#' Miller, D.L., Rexstad, E., Thomas, L., Marshall, L. and Laake, J.L. (2019)
#' Distance sampling in R. *Journal of Statistical Software* 89(1):1-28.
#' \doi{10.18637/jss.v089.i01}
#'
#' Laake, J.L. and Borchers, D.L. (2004) Methods for incomplete detection at
#' distance zero. In S.T. Buckland et al. (eds) *Advanced Distance Sampling*,
#' pp. 108-189. Oxford University Press.
#'
#' Marques, F.F.C. and Buckland, S.T. (2004) Covariate models for the detection
#' function. In S.T. Buckland et al. (eds) *Advanced Distance Sampling*,
#' pp. 31-47. Oxford University Press.
#'
"_PACKAGE"
