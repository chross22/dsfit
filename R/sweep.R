#' Fit a set of detection functions and compare them
#'
#' Fits every candidate in a [model_set()] to the same data, at the same
#' truncation, and returns a selection table ranked on more than AIC.
#'
#' @section Why truncation is one argument and not a column:
#' AIC compares models fitted to the same data. Changing the truncation changes
#' which detections are in the likelihood, so AICs either side of that are not
#' on the same scale and ranking them together is meaningless.
#'
#' Truncation is therefore a single value for a whole sweep, by construction.
#' To compare truncations, run a sweep for each and compare the *sweeps* —
#' which is a different question, answered by looking at whether the chosen
#' model and its `p̄` are stable, not by reading one AIC column.
#'
#' @section Read `p_cv` and `esw`, not just `delta_aic`:
#' What propagates into abundance is the effective strip half-width. Two models
#' within 2 AIC of one another can imply materially different values of it, and
#' the difference goes straight into density. The table therefore carries `p`
#' (average detection probability), its standard error and CV, and `esw`
#' alongside the AIC, and [selection_table()] sorts on AIC only because
#' something has to be first.
#'
#' A goodness-of-fit statistic is reported too: Cramér-von Mises, whose null is
#' that the fitted function describes the observed distances. A small p-value is
#' evidence against the model, and a model can be top of the AIC ranking and
#' still fail it.
#'
#' @section Binned and exact distances cannot be swept together:
#' `STRIP`-derived distances are intervals and are fitted binned; angle- and
#' position-derived distances are points. The likelihoods differ, so their AICs
#' are not comparable, and a table containing both would be a survey-era
#' boundary rather than a model set. Data carrying both is refused.
#'
#' @section Left truncation and the gamma key are alternatives:
#' Both describe an aircraft that cannot see beneath itself. Applying `left`
#' to a model set that contains the gamma key accounts for the blind spot twice,
#' and warns.
#'
#' @section This does not estimate g(0):
#' Every fit here conditions on the animal having been available and seen. For a
#' deep-diving species that assumption biases density low, and no amount of
#' model selection detects it — a mis-specified `g(0)` shifts every candidate in
#' the set by the same factor, so the ranking looks untouched. Correct for it at
#' the abundance step, from external sources, with its standard error
#' propagated.
#'
#' @param data A data frame with one row per detection. Needs `distance` for
#'   exact distances, or `distbegin` and `distend` for binned ones. The flatfile
#'   from `distsamp::detection_data()` is this shape; rows with no detection are
#'   dropped here.
#' @param models A [model_set()], or `NULL` for the default set.
#' @param truncation Right truncation distance. Required: it defines the data
#'   every model in the sweep is fitted to.
#' @param left Left truncation distance, for a platform that cannot see beneath
#'   itself. `NULL` (default) for none.
#' @param breaks Bin cutpoints, for binned data. Taken from `distbegin`/
#'   `distend` when not given.
#' @param quiet Suppress the progress and failure report. Default `FALSE`.
#'
#' @return An object of class `dsfit_sweep`: a list with `table` (the selection
#'   table), `fits` (the `ddf` objects, named by `model_id`), and `settings`.
#'
#' @references
#' Miller, D.L., Rexstad, E., Thomas, L., Marshall, L. and Laake, J.L. (2019)
#' Distance sampling in R. *Journal of Statistical Software* 89(1):1-28.
#' \doi{10.18637/jss.v089.i01}
#'
#' Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L.
#' and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford University
#' Press.
#'
#' @seealso [model_set()], [selection_table()]
#'
#' @examplesIf requireNamespace("mrds", quietly = TRUE)
#' set.seed(1)
#' d <- data.frame(object = 1:200, distance = abs(rnorm(200, 0, 120)))
#' d <- d[d$distance < 400, ]
#'
#' sw <- sweep_models(d, model_set(c("hn", "hr", "gamma")), truncation = 400)
#' sw
#'
#' # AIC alone would not have told you this
#' selection_table(sw)[, c("model_id", "delta_aic", "p", "p_cv", "esw")]
#'
#' @export
sweep_models <- function(data, models = NULL, truncation, left = NULL,
                         breaks = NULL, quiet = FALSE) {
  check_mrds()
  if (is.null(models)) models <- model_set()
  if (missing(truncation)) {
    rlang::abort(paste0(
      "`truncation` is required. It defines the data every model in the sweep ",
      "is fitted to, and AIC only compares models fitted to the same data."
    ))
  }
  stopifnot(is.numeric(truncation), length(truncation) == 1L, truncation > 0)

  prepared <- prepare_distance_data(data, truncation = truncation, left = left,
                                    breaks = breaks)

  if (any(models$key == "gamma") && !is.null(left)) {
    rlang::warn(paste0(
      "The model set contains the gamma key and `left` is set. Both describe a ",
      "platform that cannot see beneath itself, so applying them together ",
      "accounts for the blind spot twice. Choose one deliberately."
    ))
  }

  meta <- list(width = truncation)
  if (!is.null(left)) meta$left <- left
  if (prepared$binned) {
    meta$binned <- TRUE
    meta$breaks <- prepared$breaks
  }

  fits <- vector("list", nrow(models))
  names(fits) <- models$model_id
  failures <- character(0)

  for (i in seq_len(nrow(models))) {
    fit <- try(
      mrds::ddf(
        dsmodel = stats::as.formula(models$dsmodel[i]),
        data = prepared$data,
        method = "ds",
        meta.data = meta
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      failures <- c(failures, stats::setNames(
        sub("\n.*", "", conditionMessage(attr(fit, "condition"))),
        models$model_id[i]
      ))
      fits[[i]] <- NULL
    } else {
      fits[[i]] <- fit
    }
  }

  out <- list(
    table = build_selection_table(models, fits),
    fits = fits,
    settings = list(
      truncation = truncation, left = left, binned = prepared$binned,
      breaks = prepared$breaks, n = nrow(prepared$data),
      n_dropped = prepared$n_dropped, failures = failures
    )
  )
  class(out) <- "dsfit_sweep"

  if (!quiet) report_sweep(out)
  out
}


#' The ranked selection table from a sweep
#'
#' @param x A `dsfit_sweep` from [sweep_models()].
#' @param converged_only Drop models that failed to fit. Default `TRUE`.
#'
#' @return A tibble, sorted by AIC, with `model_id`, `key`, `adjustment`,
#'   `order`, `formula`, `converged`, `n_par`, `aic`, `delta_aic`, `p`, `p_se`,
#'   `p_cv`, `esw`, and `cvm_p`.
#'
#' @seealso [sweep_models()]
#'
#' @examplesIf requireNamespace("mrds", quietly = TRUE)
#' set.seed(1)
#' d <- data.frame(object = 1:150, distance = abs(rnorm(150, 0, 100)))
#' d <- d[d$distance < 300, ]
#' selection_table(sweep_models(d, truncation = 300, quiet = TRUE))
#'
#' @export
selection_table <- function(x, converged_only = TRUE) {
  stopifnot(inherits(x, "dsfit_sweep"))
  out <- x$table
  if (converged_only) out <- out[out$converged, , drop = FALSE]
  out
}

# One row per candidate, whether or not it fitted.
build_selection_table <- function(models, fits) {
  n <- nrow(models)
  grab <- function(f, what) {
    if (is.null(f)) return(NA_real_)
    what(f)
  }

  aic <- vapply(fits, function(f) if (is.null(f)) NA_real_ else f$criterion,
                numeric(1))
  stats_list <- lapply(fits, fit_summary)

  out <- models[, c("model_id", "key", "adjustment", "order", "formula")]
  out$converged <- !vapply(fits, is.null, logical(1))
  out$n_par <- vapply(stats_list, function(s) s$n_par, numeric(1))
  out$aic <- aic
  out$delta_aic <- aic - min(aic, na.rm = TRUE)
  out$p <- vapply(stats_list, function(s) s$p, numeric(1))
  out$p_se <- vapply(stats_list, function(s) s$p_se, numeric(1))
  out$p_cv <- out$p_se / out$p
  out$esw <- vapply(stats_list, function(s) s$esw, numeric(1))
  out$cvm_p <- vapply(stats_list, function(s) s$cvm_p, numeric(1))

  out[order(out$aic, na.last = TRUE), , drop = FALSE]
}

# What a fitted ddf gives back, in one shape. NA all through when it failed.
fit_summary <- function(f) {
  blank <- list(n_par = NA_real_, p = NA_real_, p_se = NA_real_,
                esw = NA_real_, cvm_p = NA_real_)
  if (is.null(f)) return(blank)

  s <- try(summary(f), silent = TRUE)
  if (inherits(s, "try-error")) return(blank)

  # Cramer-von Mises: the null is that the fitted function describes the
  # observed distances, so a small p-value is evidence against the model.
  cvm <- try(mrds::ddf.gof(f, qq = FALSE)$dsgof$CvM$p, silent = TRUE)
  if (inherits(cvm, "try-error") || is.null(cvm)) cvm <- NA_real_

  list(
    n_par = length(f$par),
    p = s$average.p,
    p_se = s$average.p.se,
    # Effective strip half-width: the width of the strip that would have given
    # the same number of detections under certain detection. This is what
    # propagates into density, which is why it sits next to the AIC.
    esw = s$average.p * s$width,
    cvm_p = as.numeric(cvm)
  )
}

report_sweep <- function(x) {
  tab <- x$table
  ok <- sum(tab$converged)
  lines <- paste0(
    "Fitted ", ok, " of ", nrow(tab), " candidate models to ",
    x$settings$n, " detections",
    if (x$settings$binned) " (binned)" else "",
    ", truncation ", x$settings$truncation,
    if (!is.null(x$settings$left)) paste0(", left ", x$settings$left) else "",
    "."
  )
  if (x$settings$n_dropped > 0) {
    lines <- c(lines, paste0(
      "  ", x$settings$n_dropped,
      " row(s) dropped: no distance, or beyond the truncation."
    ))
  }
  if (length(x$settings$failures)) {
    lines <- c(lines, paste0(
      "  did not fit: ",
      paste(names(x$settings$failures), collapse = ", ")
    ))
  }
  lines <- c(lines, "  g(0) = 1 is assumed; see `?sweep_models`.")
  rlang::inform(paste(lines, collapse = "\n"))
}

check_mrds <- function(call = rlang::caller_env()) {
  if (!requireNamespace("mrds", quietly = TRUE)) {
    rlang::abort(
      paste0(
        "The `mrds` package is required to fit detection functions. Install it ",
        "with `install.packages(\"mrds\")`."
      ),
      call = call
    )
  }
}
