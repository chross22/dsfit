#' @keywords internal
#'
#' @section Where this sits:
#' `dsfit` is the middle of three layers. `distsamp` turns NARWC-format survey
#' data into effort segments and detection distances; `dsfit` fits and compares
#' detection functions; an analysis repository holds the choices that change
#' every run — which years, which truncation, which covariates.
#'
#' The split exists so that the part with testable logic is a package and the
#' part with project-specific choices is not.
#'
#' @section Why not `Distance::ds()`:
#' `Distance` is the right tool for fitting one model interactively, and this
#' package does not replace it. But `ds()` offers `hn`, `hr`, and `unif` only,
#' and the gamma key — the natural model for an aerial platform that cannot see
#' beneath itself — is available in `mrds` alone. Fitting the whole model set
#' through `mrds::ddf()` also keeps every AIC in a selection table coming off
#' the same likelihood machinery, rather than mixing two wrappers' handling of
#' truncation and monotonicity and presenting the difference as a model
#' difference.
#'
#' @section What it will not do:
#' Estimate `g(0)`. Every fit conditions on the animal having been available and
#' seen, and a mis-specified `g(0)` shifts every candidate in a model set by the
#' same factor — so the ranking looks untouched while the density is wrong.
#' It cannot be estimated from a standard NARWC extract, which records neither
#' dive data nor the double-observer structure mark-recapture needs. Correct for
#' it at the abundance step, from external sources, with its standard error
#' propagated.
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
