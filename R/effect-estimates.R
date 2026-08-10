#' Extract a fitted detection function as a tidy frame
#'
#' An [fancyfx::effect_estimates()] method for `mrds` detection functions, so
#' that a fitted `ddf` can be plotted with [fancyfx::plotEffects()] and
#' [fancyfx::plotRugs()] like any other model.
#'
#' @section Why this method has to exist:
#' `fancyfx` sends anything that is not a GAM to `marginaleffects`, which cannot
#' introspect a `ddf` — it reports no predictors at all. The generic is designed
#' to be extended from outside the package, and this is that extension. Nothing
#' downstream needs to change: the transforms, the ribbon, the rug, and the
#' panel arranging all work off the frame returned here.
#'
#' @section What it returns, and what it does not:
#' The detection function \eqn{g(x)}: the probability that an animal at
#' perpendicular distance \eqn{x} was detected, given that it was available to
#' be. It is **not** the probability that an animal at that distance was there
#' and seen, and no `g(0)` correction is applied — see [sweep_models()].
#'
#' A rug of the observed distances beneath this curve is the standard
#' detection-function diagnostic, which is `fancyfx`'s premise applied almost
#' exactly.
#'
#' @section Covariate models: the curve is averaged, not evaluated at a mean:
#' With detection covariates there is one \eqn{g(x)} per covariate combination.
#' This averages those curves over the covariate values actually observed,
#' rather than evaluating the function once at mean covariates. The two are
#' different, and averaging the curves is the one that has a meaning: it is the
#' mean detection function of the animals in the sample.
#'
#' `at` fixes named covariates at chosen values; anything not named is still
#' averaged over as observed. `at = NULL` averages over everything.
#'
#' @section Why the area under this curve is not `esw`:
#' For an intercept-only model, the mean of this curve over the truncation width
#' is exactly the fitted average detection probability, and the area under it is
#' the effective strip half-width in [selection_table()]. The tests assert that.
#'
#' For a **covariate** model the two differ slightly, and it is worth knowing
#' why rather than discovering it. `mrds` reports `average.p` as a
#' Horvitz-Thompson mean, \eqn{n / \sum 1/p_i}, which weights each animal by the
#' inverse of its own detection probability — correcting for the fact that a
#' sample of detections over-represents the covariate values that are easy to
#' detect. The curve here is the plain arithmetic mean over the animals as
#' observed, \eqn{(1/n) \sum g(x \mid z_i)}, which is what you want to look at
#' next to a histogram of those same animals' distances.
#'
#' The two coincide when detection probability is constant, and differ by a
#' fraction of a percent otherwise. `esw` in the selection table uses the
#' Horvitz-Thompson version, so eyeballing the area under this curve will not
#' reproduce it exactly for a covariate model.
#'
#' @section The interval is a delta-method approximation:
#' \eqn{g(x)} is deterministic given the fitted parameters, so the uncertainty
#' comes from them: the Jacobian of \eqn{g} with respect to the parameter vector
#' is taken numerically and combined with the inverse Hessian. The result is
#' clamped to `[0, 1]`, since a probability cannot leave it — which means an
#' interval touching 0 or 1 is a bound, not an estimate. Where the Hessian
#' cannot be inverted the estimates are returned with missing bounds rather than
#' no estimates.
#'
#' @param model A fitted `ddf` object, from `mrds::ddf()` or [sweep_models()].
#' @param var Must be `"distance"`. A detection function has one predictor;
#'   plotting detectability against a covariate is a different quantity and is
#'   not built yet.
#' @param scale Accepted for compatibility with the generic, and ignored. A
#'   detection function has no separable linear predictor, so every setting
#'   gives \eqn{g(x)} on the probability scale.
#' @param interval `"se"` for a one-standard-error ribbon, `"ci"` for a
#'   confidence interval at `level`, or `"auto"` (which `fancyfx` forwards by
#'   default) for `"se"`. `"cri"` is refused: there is no posterior here to
#'   summarise.
#' @param level Confidence level, used when `interval = "ci"`.
#' @param n Number of distances at which to evaluate the function.
#' @param at Named list fixing covariate values, or `NULL` (default) to average
#'   over the covariates as observed.
#' @param ... Ignored.
#'
#' @return A data frame with `.x`, `.estimate`, `.lower`, `.upper`, and a
#'   `"quantity"` attribute, as [fancyfx::effect_estimates()] specifies.
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
#' # The integral of the curve is the fitted average detection probability
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
  # Resolved by hand rather than with match.arg(), and against the generic's
  # whole vocabulary rather than this method's preferences. `plotEffects()`
  # forwards its own unresolved defaults, so "auto" arrives here routinely and
  # match.arg() rejects the vector outright.
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

  if (!identical(as.character(var), "distance")) {
    rlang::abort(paste0(
      "`var` must be \"distance\". A detection function has one predictor, and ",
      "the effect this returns is g(x) against it. Plotting detectability ",
      "against a covariate is a different quantity - average detection ",
      "probability as the covariate varies - and is not built yet."
    ))
  }

  width <- model$meta.data$width
  left <- model$meta.data$left
  if (is.null(left) || !is.finite(left)) left <- 0
  x <- seq(left, width, length.out = n)

  ddfobj <- covariate_object(model, at)
  gfun <- averaged_g(ddfobj, x, width)

  est <- gfun(model$par)
  se <- delta_se(gfun, model)

  mult <- if (interval == "se") 1 else stats::qnorm(1 - (1 - level) / 2)
  out <- data.frame(
    .x = x,
    .estimate = est,
    .lower = pmax(0, est - mult * se),
    .upper = pmin(1, est + mult * se)
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
  effect_estimates.ddf(model, var = var, scale = scale, interval = interval,
                       level = level, n = n, at = at, ...)
}

# g(x) averaged over the covariate values in a ddfobj, as a function of the
# parameter vector.
#
# `detfct()` pairs its distance argument with the object's scale vector
# *element-wise*, recycling the shorter one - it does not evaluate every
# distance at every covariate row. Handing it a grid and a design matrix of
# different lengths therefore returns whichever recycling happens to produce,
# which is only the intended answer when the scale is constant. The design
# matrix is expanded here so that every pair is explicit.
#
# Rows are expanded once per *unique* covariate combination and weighted by how
# often it occurs, so an intercept-only model costs one column rather than one
# per detection.
averaged_g <- function(ddfobj, x, width) {
  dm <- ddfobj$scale$dm
  key <- do.call(paste, c(as.data.frame(dm), sep = "\r"))
  ukey <- unique(key)
  first <- match(ukey, key)
  weight <- tabulate(match(key, ukey), nbins = length(ukey))

  n_x <- length(x)
  n_u <- length(ukey)
  rows <- rep(first, each = n_x)

  obj <- ddfobj
  obj$scale$dm <- dm[rows, , drop = FALSE]
  if (!is.null(ddfobj$xmat)) {
    obj$xmat <- ddfobj$xmat[rows, , drop = FALSE]
  }
  xs <- rep(x, times = n_u)

  function(par) {
    fitted <- mrds::detfct(xs, mrds::assign.par(obj, par), width = width)
    as.numeric(matrix(fitted, nrow = n_x) %*% weight) / sum(weight)
  }
}

# Take the first element when handed an unresolved default, and check
# membership. match.arg() only tolerates a vector when it is identical to the
# formal default, which is not true of a value forwarded from another function.
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

# The ddfobj with its scale design matrix rebuilt at the requested covariate
# values. Covariates not named in `at` keep the values they were observed at,
# so they are averaged over rather than pinned to something arbitrary.
covariate_object <- function(model, at) {
  ddfobj <- model$ds$aux$ddfobj
  if (is.null(at)) {
    return(ddfobj)
  }
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
  for (nm in names(at)) dat[[nm]] <- at[[nm]]

  form <- ddfobj$scale$formula
  if (is.character(form)) form <- stats::as.formula(form)
  ddfobj$scale$dm <- stats::model.matrix(form, data = dat)
  ddfobj
}

# Delta method on the fitted parameters. Returns zeros with a warning when the
# Hessian will not invert, so a curve is still plotted rather than nothing.
delta_se <- function(gfun, model) {
  n_x <- length(gfun(model$par))
  vcov <- try(solve(model$hessian), silent = TRUE)
  if (inherits(vcov, "try-error")) {
    rlang::warn(paste0(
      "The Hessian could not be inverted, so no interval is available for ",
      "this fit. The estimates are returned without one."
    ))
    return(rep(NA_real_, n_x))
  }

  jac <- numDeriv::jacobian(gfun, model$par)
  variance <- rowSums((jac %*% vcov) * jac)
  sqrt(pmax(0, variance))
}
