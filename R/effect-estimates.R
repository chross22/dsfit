#' Extract a fitted detection function as a plottable data frame
#'
#' Evaluates a fitted `mrds` detection function \eqn{g(x)} on a grid of
#' distances, with a delta-method interval, and returns it in the format that
#' [fancyfx::plotEffects()] and [fancyfx::plotRugs()] draw. The detection
#' function \eqn{g(x)} is the probability that an animal at perpendicular
#' distance \eqn{x} from the trackline is detected.
#'
#' @section How it works:
#' 1. `interval` and `scale` are checked. `"auto"` becomes `"se"`, and
#'    `"cri"` is refused.
#' 2. `var` is checked: only `"distance"` is accepted.
#' 3. A grid of `n` distances is built from the left truncation distance (0 if
#'    none) to the right truncation width.
#' 4. If `at` is given, the named covariates are set to those values.
#' 5. \eqn{g(x)} is evaluated at every grid distance for every observed
#'    covariate combination, and the curves are averaged.
#' 6. The standard error at each distance is computed by the delta method
#'    from the fitted parameters and their inverse Hessian.
#' 7. The lower and upper bounds are the estimate minus and plus one standard
#'    error (`"se"`) or a normal quantile times it (`"ci"`), clamped to
#'    `[0, 1]`.
#'
#' @section What the curve is:
#' The curve is \eqn{g(x)} conditional on the animal being available to be
#' seen. No `g(0)` correction is applied; see [sweep_models()] and [g0()].
#' Plotting it over a rug of the observed distances is the standard check of a
#' detection function's fit.
#'
#' @section Covariate models:
#' A model with detection covariates has one \eqn{g(x)} per covariate
#' combination. This function averages those curves over the covariate values
#' of the detected animals. It does not evaluate one curve at the mean
#' covariate values, which is a different and less meaningful quantity.
#'
#' `at` fixes named covariates at chosen values. Covariates not named in `at`
#' are still averaged over their observed values.
#'
#' @section Relation to `esw` and average detection probability:
#' The average detection probability \eqn{p} is the mean of \eqn{g(x)} over the
#' truncation interval. The effective strip half-width (`esw`) in
#' [selection_table()] is \eqn{p} times the width of that interval.
#'
#' For an intercept-only model, the mean of the curve returned here equals the
#' fitted \eqn{p}. The package tests check this.
#'
#' For a covariate model the two differ slightly. `mrds` reports `average.p` as
#' a Horvitz-Thompson mean, \eqn{n / \sum 1/p_i}, which weights each animal by
#' the inverse of its own detection probability. This curve is the plain mean
#' over the detected animals, \eqn{(1/n) \sum g(x \mid z_i)}, which is the right
#' curve to compare with a histogram of their distances. On the package's
#' covariate test fixture (half-normal, `~bf`, width 400) the plain mean is
#' 0.401 and `average.p` is 0.388, a difference of 3.5%. The area under this
#' curve therefore does not reproduce `esw` exactly for a covariate model.
#'
#' @section The interval:
#' \eqn{g(x)} is a fixed function of the fitted parameters, so its uncertainty
#' comes from theirs. The Jacobian of \eqn{g} with respect to the parameters is
#' computed numerically and combined with the inverse Hessian (the delta
#' method). Bounds are clamped to `[0, 1]`, so a bound at 0 or 1 is a limit,
#' not an estimate. If the Hessian cannot be inverted, a warning is given and
#' the bounds are `NA`.
#'
#' `fancyfx` sends models it does not recognise to `marginaleffects`, which
#' finds no predictors in a `ddf` object. This method is what lets a `ddf` be
#' plotted.
#'
#' @param model A fitted `ddf` object, from `mrds::ddf()` or [sweep_models()].
#' @param var Must be `"distance"`. Plotting detectability against a covariate
#'   is not supported.
#' @param scale Accepted for compatibility with the generic and ignored. The
#'   result is always on the probability scale.
#' @param interval `"se"` for plus and minus one standard error, `"ci"` for a
#'   confidence interval at `level`, or `"auto"` (the default) for `"se"`.
#'   `"cri"` is an error, because there is no posterior.
#' @param level Confidence level, used only when `interval = "ci"`.
#' @param n Number of distances at which to evaluate \eqn{g(x)}.
#' @param at Named list of covariate values to fix, or `NULL` (default) to
#'   average over the observed covariate values.
#' @param ... Ignored.
#'
#' @return A data frame with one row per grid distance and columns `.x`
#'   (distance), `.estimate` (\eqn{g(x)}), `.lower` and `.upper` (interval
#'   bounds). Its `"quantity"` attribute is `"Detection probability g(x)"`,
#'   which `fancyfx` uses as the axis label. This is the format
#'   [fancyfx::effect_estimates()] specifies.
#'
#' @references
#' Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L.
#' and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford University
#' Press.
#'
#' @seealso [sweep_models()]
#'
#' @examplesIf requireNamespace("mrds", quietly = TRUE)
#' set.seed(1)
#' d <- data.frame(object = 1:200, distance = abs(rnorm(200, 0, 120)))
#' d <- d[d$distance < 400, ]
#' fit <- mrds::ddf(dsmodel = ~cds(key = "hn", formula = ~1), data = d,
#'                  method = "ds", meta.data = list(width = 400))
#'
#' est <- effect_estimates_ddf(fit)
#' head(est)
#' attr(est, "quantity")
#'
#' # The mean of the curve is the fitted average detection probability
#' mean(est$.estimate) # approximately summary(fit)$average.p
#'
#' @rdname effect_estimates_ddf
#' @exportS3Method fancyfx::effect_estimates
effect_estimates.ddf <- function(model, var = "distance",
                                 scale = c("auto", "link", "response"),
                                 interval = c("auto", "se", "ci"),
                                 level = 0.95,
                                 n = 100,
                                 at = NULL,
                                 ...) {
  # Step 1: check `interval` and `scale`. match.arg() is not used because
  # plotEffects() forwards its own unresolved default vectors, which
  # match.arg() rejects. "cri" is accepted here so it can be refused below
  # with a clear message.
  interval <- resolve_choice(interval, c("auto", "se", "ci", "cri"), "interval")
  resolve_choice(scale, c("auto", "link", "response"), "scale")
  if (identical(interval, "auto")) interval <- "se"
  if (identical(interval, "cri")) {
    rlang::abort(paste0(
      "`interval = \"cri\"` asks for a credible interval, and this is a ",
      "frequentist fit with no posterior to summarise. Use \"se\" or \"ci\"."
    ))
  }
  check_mrds()

  # Step 2: only distance is supported as the predictor.
  if (!identical(as.character(var), "distance")) {
    rlang::abort(paste0(
      "`var` must be \"distance\". A detection function has one predictor, and ",
      "the effect this returns is g(x) against it. Plotting detectability ",
      "against a covariate is a different quantity - average detection ",
      "probability as the covariate varies - and is not built yet."
    ))
  }

  # Step 3: grid of distances from the left truncation (0 if none) to the
  # right truncation width.
  width <- model$meta.data$width
  left <- model$meta.data$left
  if (is.null(left) || !is.finite(left)) left <- 0
  x <- seq(left, width, length.out = n)

  # Step 4: apply `at`. Step 5: build g(x), averaged over covariate values,
  # as a function of the parameter vector.
  ddfobj <- covariate_object(model, at)
  gfun <- averaged_g(ddfobj, x, width)

  # Step 6: estimate at the fitted parameters, and delta-method standard
  # errors.
  est <- gfun(model$par)
  se <- delta_se(gfun, model)

  # Step 7: interval half-width multiplier, then bounds.
  mult <- if (interval == "se") 1 else stats::qnorm(1 - (1 - level) / 2)
  out <- data.frame(
    .x = x,
    .estimate = est,
    # Bounds are clamped to [0, 1] but never moved past the estimate. An
    # adjustment model fitted without a monotonicity constraint can put g(x)
    # above 1, and the band must still contain the curve.
    .lower = pmin(est, pmax(0, est - mult * se)),
    .upper = pmax(est, pmin(1, est + mult * se))
  )
  attr(out, "quantity") <- "Detection probability g(x)"
  out
}

#' @rdname effect_estimates_ddf
#' @export
effect_estimates_ddf <- function(model, var = "distance",
                                 scale = c("auto", "link", "response"),
                                 interval = c("auto", "se", "ci"),
                                 level = 0.95, n = 100, at = NULL, ...) {
  # A plain exported name for the S3 method, so it can be called without
  # fancyfx attached.
  effect_estimates.ddf(model, var = var, scale = scale, interval = interval,
                       level = level, n = n, at = at, ...)
}

# averaged_g(): build g(x), averaged over covariate values, as a function of
# the parameter vector.
#
# Inputs: `ddfobj`, the detection-function object from a fit (possibly with
# covariates changed by covariate_object()); `x`, the grid of distances;
# `width`, the right truncation width.
#
# Returns: a function of the parameter vector `par` that returns g(x) at each
# grid distance, averaged over the detected animals' covariate values.
#
# Why it is built this way:
# - mrds::detfct() pairs each distance with one row of the design matrix,
#   element by element, and recycles the shorter of the two. It does not
#   evaluate every distance for every row. So the design matrix is expanded
#   here to one row per (covariate combination, distance) pair.
# - Hazard-rate and gamma keys also have a shape design matrix, with one row
#   per detection. It is expanded the same way; otherwise it would not match
#   the expanded scale matrix.
# - Only unique covariate combinations are expanded, each weighted by how
#   often it occurs. An intercept-only model then needs one curve, not one per
#   detection.
# - Returning a function of `par` lets delta_se() differentiate it
#   numerically.
averaged_g <- function(ddfobj, x, width) {
  dm <- ddfobj$scale$dm
  shape_dm <- ddfobj$shape$dm
  # The uniform key has no design matrices: one curve applies to every
  # detection.
  if (is.null(dm) && is.null(shape_dm)) {
    return(function(par) {
      as.numeric(mrds::detfct(x, mrds::assign.par(ddfobj, par), width = width))
    })
  }

  # Identify unique covariate combinations (rows of the scale and shape design
  # matrices together), the first row of each, and how many detections share
  # it.
  design <- as.data.frame(dm)
  if (!is.null(shape_dm)) design <- cbind(design, as.data.frame(shape_dm))
  key <- do.call(paste, c(unname(design), sep = "\r"))
  ukey <- unique(key)
  first <- match(ukey, key)
  weight <- tabulate(match(key, ukey), nbins = length(ukey))

  # Expand: each unique combination is repeated once per grid distance.
  n_x <- length(x)
  n_u <- length(ukey)
  rows <- rep(first, each = n_x)

  obj <- ddfobj
  if (!is.null(dm)) obj$scale$dm <- dm[rows, , drop = FALSE]
  if (!is.null(shape_dm)) {
    obj$shape$dm <- shape_dm[rows, , drop = FALSE]
  }
  if (!is.null(ddfobj$xmat)) {
    obj$xmat <- ddfobj$xmat[rows, , drop = FALSE]
  }
  # The distance grid, repeated once per unique combination, to match.
  xs <- rep(x, times = n_u)

  # Evaluate all pairs, reshape to one column per combination, and take the
  # weighted mean across columns.
  function(par) {
    fitted <- mrds::detfct(xs, mrds::assign.par(obj, par), width = width)
    as.numeric(matrix(fitted, nrow = n_x) %*% weight) / sum(weight)
  }
}

# resolve_choice(): pick one value from an argument and check it is allowed.
#
# Inputs: `value`, the argument as passed (possibly a whole default vector);
# `allowed`, the accepted values; `label`, the argument name for the error.
#
# Returns: the single chosen value. Errors if it is not in `allowed`.
#
# Why: match.arg() accepts a vector only when it is identical to the
# function's own default. A default forwarded from another function, such as
# fancyfx::plotEffects(), is not, so the first element is taken here instead.
resolve_choice <- function(value, allowed, label) {
  if (length(value) > 1) value <- value[1]
  if (!value %in% allowed) {
    rlang::abort(paste0(
      "Unknown ", label, ": \"", value, "\". Use one of ",
      paste0("\"", allowed, "\"", collapse = ", "), "."
    ))
  }
  value
}

# covariate_object(): the fit's detection-function object, with covariates
# fixed at the values in `at`.
#
# Inputs: `model`, a fitted ddf; `at`, a named list of covariate values or
# NULL.
#
# Returns: the ddfobj, with its scale design matrix rebuilt from the data
# after the `at` values are substituted. Unchanged when `at` is NULL.
#
# Why: covariates not named in `at` keep their observed values, so they are
# still averaged over by averaged_g() rather than set to an arbitrary value.
covariate_object <- function(model, at) {
  ddfobj <- model$ds$aux$ddfobj
  if (is.null(at)) {
    return(ddfobj)
  }

  # Check `at` is a named list of columns that exist in the model data.
  if (!is.list(at) || is.null(names(at)) || any(names(at) == "")) {
    rlang::abort("`at` must be a named list of covariate values.")
  }

  dat <- model$data
  unknown <- setdiff(names(at), names(dat))
  if (length(unknown)) {
    rlang::abort(paste0(
      "Not a covariate in this model: ",
      paste0("`", unknown, "`", collapse = ", "), "."
    ))
  }

  # Replace each named column with the requested value, in the column's type.
  for (nm in names(at)) dat[[nm]] <- at_value(dat[[nm]], at[[nm]], nm)

  # Rebuild the scale design matrix from the model's scale formula.
  form <- ddfobj$scale$formula
  if (is.character(form)) form <- stats::as.formula(form)
  ddfobj$scale$dm <- stats::model.matrix(form, data = dat)
  ddfobj
}

# at_value(): convert a value from `at` to the type of the column it replaces.
#
# Inputs: `observed`, the column in the model data; `value`, the requested
# value; `name`, the column name for the error message.
#
# Returns: `value` unchanged for a numeric column; for a factor or character
# column, a factor with the observed levels. Errors if `value` is not one of
# those levels.
#
# Why: a level given as a string (e.g. "B") would otherwise become a
# one-level character column, which model.matrix() rejects with an error about
# contrasts. Keeping the observed levels gives the design matrix the same
# columns as the fitted model.
at_value <- function(observed, value, name) {
  if (!is.factor(observed) && !is.character(observed)) return(value)
  lv <- if (is.factor(observed)) levels(observed) else
    sort(unique(observed[!is.na(observed)]))
  bad <- setdiff(as.character(value), lv)
  if (length(bad)) {
    rlang::abort(paste0(
      "`at$", name, "` = ", paste0("\"", bad, "\"", collapse = ", "),
      " is not a level of `", name, "`. Levels: ",
      paste0("\"", lv, "\"", collapse = ", "), "."
    ))
  }
  factor(as.character(value), levels = lv)
}

# delta_se(): delta-method standard error of g(x) at each grid distance.
#
# Inputs: `gfun`, the function of the parameter vector from averaged_g();
# `model`, the fitted ddf (for its parameters and Hessian).
#
# Returns: a numeric vector of standard errors, one per grid distance. All NA,
# with a warning, if the Hessian cannot be inverted.
#
# Why: g(x) depends on the data only through the fitted parameters, so
# Var(g) is approximately J V J', where J is the Jacobian of g with respect to
# the parameters and V is the inverse Hessian. Returning NA rather than an
# error means the curve can still be plotted without a band.
delta_se <- function(gfun, model) {
  n_x <- length(gfun(model$par))

  # Parameter variance-covariance matrix from the inverse Hessian.
  vcov <- try(solve(model$hessian), silent = TRUE)
  if (inherits(vcov, "try-error")) {
    rlang::warn(paste0(
      "The Hessian could not be inverted, so no interval is available for ",
      "this fit. The estimates are returned without one."
    ))
    return(rep(NA_real_, n_x))
  }

  # Numerical Jacobian (one row per distance), then the diagonal of J V J'.
  jac <- numDeriv::jacobian(gfun, model$par)
  variance <- rowSums((jac %*% vcov) * jac)
  # Rounding can give tiny negative variances; floor them at 0.
  sqrt(pmax(0, variance))
}
