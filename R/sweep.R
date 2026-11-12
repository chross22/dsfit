#' Fit a set of detection functions and compare them
#'
#' Fits every candidate in a [model_set()] to the same data at the same
#' truncation, using `mrds::ddf()`, and returns a selection table. Use it to
#' choose a detection function for a survey.
#'
#' @section How it works:
#' 1. Checks that `mrds` is installed, that `models` is a model set, and that
#'    `truncation` is a single positive number.
#' 2. Prepares the data with [prepare_distance_data()]: validates the
#'    distances, drops rows with no distance or outside the truncation, and
#'    decides whether the fit is binned.
#' 3. Warns if `left` is set and the model set contains the gamma key.
#' 4. Fits each candidate with `mrds::ddf()`. A candidate that errors or does
#'    not converge is recorded as a failure and the sweep continues.
#' 5. Builds the selection table, sorted by AIC, and prints a short report
#'    unless `quiet = TRUE`.
#'
#' @section Truncation is one value per sweep:
#' AIC only compares models fitted to the same data. A different truncation
#' keeps a different set of detections, so the AICs are not comparable. To
#' compare truncations, run one sweep per truncation and check whether the
#' chosen model and its `p` stay stable.
#'
#' @section What the table reports besides AIC:
#' The effective strip half-width (esw) is the width of a strip in which
#' detection is certain that would yield the same expected number of
#' detections. It goes directly into density. Two models within 2 AIC of each
#' other can have quite different esw, so the table reports `p` (average
#' detection probability), its standard error and CV, and `esw` next to AIC.
#'
#' A goodness-of-fit p-value is also reported. The null hypothesis is that the
#' fitted function describes the observed distances, so a small p-value is
#' evidence against the model. The top AIC model can still fail this test.
#'
#' @section Which goodness-of-fit test to read:
#' Exact distances get a Cramer-von Mises test (`cvm_p`), which tests their
#' empirical distribution. Binned distances have no empirical distribution, so
#' they get a chi-square test over the survey's own bins (`chisq_p`).
#' `chisq_p` is also filled in for exact distances, but over cutpoints `mrds`
#' chooses; read `cvm_p` there. Printing a sweep shows the test that applies.
#'
#' A chi-square test over few bins can have no degrees of freedom left (for
#' example three bins and a two-parameter key). `mrds` then returns `NA`. The
#' model converged; it just has no chi-square test.
#'
#' @section Binned and exact distances cannot be mixed:
#' Distances from `STRIP` records are intervals and are fitted binned.
#' Distances from angles or positions are exact. The two likelihoods differ, so
#' their AICs are not comparable. Data containing both is refused; split it and
#' run one sweep for each.
#'
#' @section Left truncation and the gamma key:
#' Both model a platform that cannot see directly beneath itself. Using both
#' accounts for the blind spot twice, so setting `left` for a model set that
#' contains the gamma key gives a warning.
#'
#' @section Monotonicity constraint:
#' As in `Distance::ds()`, a candidate with adjustment terms and no covariates
#' is fitted with a strict monotonicity constraint (`mono = TRUE,
#' mono.strict = TRUE` in `mrds`). This stops g(x), the detection probability
#' at distance x, from exceeding 1 or turning upward. AIC and `p` then match
#' `ds()` with its defaults. Covariate models are unconstrained, as in `ds()`.
#' The gamma key is also unconstrained, because it rises to a peak by design.
#'
#' @section g(0) is assumed to be 1:
#' Every fit assumes an animal on the track-line is always detected. A wrong
#' g(0) scales every candidate by the same factor, so model selection cannot
#' detect it. Correct for it at the abundance step, from external sources, and
#' propagate its standard error (see [g0()]).
#'
#' @param data Data frame with one row per detection. Needs a `distance`
#'   column for exact distances, or `distbegin` and `distend` for binned ones.
#'   The flatfile from `distsamp::detection_data()` has this shape. Rows with
#'   no distance (samples without a detection) are dropped.
#' @param models A [model_set()]. `NULL` (default) uses `model_set()`.
#' @param truncation Number. Right truncation distance, in the units of the
#'   distances. Required, with no default.
#' @param left Number. Left truncation distance, for a platform that cannot
#'   see beneath itself. `NULL` (default) for none.
#' @param breaks Numeric vector of bin cutpoints, for binned data. `NULL`
#'   (default) derives them from `distbegin` and `distend`. Ignored, with a
#'   warning, for exact distances.
#' @param quiet Logical. `TRUE` suppresses the summary message. Default
#'   `FALSE`.
#'
#' @return An object of class `dsfit_sweep`, a list with:
#' \describe{
#'   \item{`table`}{The selection table, one row per candidate including
#'     failures. See [selection_table()] for the columns.}
#'   \item{`fits`}{List of `ddf` objects named by `model_id`. `NULL` for a
#'     candidate that failed.}
#'   \item{`settings`}{List with `truncation`, `left`, `binned`, `breaks`, `n`
#'     (detections fitted), `n_dropped` (rows removed during preparation) and
#'     `failures` (first line of the error, or `"did not converge"`, named
#'     by `model_id`).}
#' }
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
#' # Compare esw and p_cv as well as delta AIC
#' selection_table(sw)[, c("model_id", "delta_aic", "p", "p_cv", "esw")]
#'
#' @export
sweep_models <- function(data, models = NULL, truncation, left = NULL,
                         breaks = NULL, quiet = FALSE) {
  # 1. Check inputs: mrds is installed, the model set is valid, and there is
  #    one positive, finite truncation distance.
  check_mrds()
  if (is.null(models)) models <- model_set()
  check_model_set(models)
  if (missing(truncation)) {
    rlang::abort(paste0(
      "`truncation` is required. It defines the data every model in the sweep ",
      "is fitted to, and AIC only compares models fitted to the same data."
    ))
  }
  stopifnot(is.numeric(truncation), length(truncation) == 1L,
            is.finite(truncation), truncation > 0)

  # 2. Validate and truncate the distances, and find out whether they are
  #    binned. Every candidate is fitted to this one data set.
  prepared <- prepare_distance_data(data, truncation = truncation, left = left,
                                    breaks = breaks)

  # 3. Warn when the blind spot under the platform is modelled twice.
  if (any(models$key == "gamma") && !is.null(left)) {
    rlang::warn(paste0(
      "The model set contains the gamma key and `left` is set. Both describe a ",
      "platform that cannot see beneath itself, so applying them together ",
      "accounts for the blind spot twice. Choose one deliberately."
    ))
  }

  # 4. Build the mrds meta.data shared by every candidate: truncation, left
  #    truncation, and the bins for a binned fit.
  meta <- list(width = truncation)
  if (!is.null(left)) meta$left <- left
  if (prepared$binned) {
    meta$binned <- TRUE
    meta$breaks <- prepared$breaks
  }

  # 5. Fit each candidate. `fits` keeps one slot per candidate, in model-set
  #    order, so the selection table can match fits to rows by position.
  fits <- vector("list", nrow(models))
  names(fits) <- models$model_id
  failures <- character(0)

  for (i in seq_len(nrow(models))) {
    # The monotonicity setting differs by candidate, so it is added here.
    fit <- try(
      mrds::ddf(
        dsmodel = stats::as.formula(models$dsmodel[i]),
        data = prepared$data,
        method = "ds",
        meta.data = c(meta, monotonicity(models[i, , drop = FALSE]))
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      # An error: record its first line and leave the slot empty.
      failures <- c(failures, stats::setNames(
        first_line(conditionMessage(attr(fit, "condition"))),
        models$model_id[i]
      ))
      # `fits[[i]] <- NULL` would delete the element and shift every later
      # fit onto the wrong row of the table. `fits[i] <- list(NULL)` empties
      # the slot instead.
      fits[i] <- list(NULL)
    } else if (!identical(as.numeric(fit$ds$converge), 0)) {
      # No convergence: mrds warns instead of erroring, and returns a fit
      # whose summary has no average.p. Treat it as a failure so it is not
      # reported as converged.
      failures <- c(failures, stats::setNames("did not converge",
                                              models$model_id[i]))
      fits[i] <- list(NULL)
    } else {
      fits[[i]] <- fit
    }
  }

  # 6. Assemble the result: table, fits, and the settings the sweep used.
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

  # 7. Report what was fitted, unless asked to be quiet.
  if (!quiet) report_sweep(out)
  out
}


#' Get the selection table from a sweep
#'
#' Returns the table of fitted candidates from a [sweep_models()] result,
#' sorted by AIC (lowest first).
#'
#' @param x A `dsfit_sweep` object from [sweep_models()].
#' @param converged_only Logical. `TRUE` (default) drops candidates that
#'   failed to fit or did not converge.
#'
#' @return A tibble with one row per candidate and columns:
#' \describe{
#'   \item{`model_id`, `key`, `adjustment`, `order`, `formula`}{The candidate,
#'     as in [model_set()].}
#'   \item{`converged`}{`TRUE` if the fit succeeded.}
#'   \item{`n_par`}{Number of estimated parameters.}
#'   \item{`aic`, `delta_aic`}{AIC, and its difference from the lowest AIC.}
#'   \item{`p`, `p_se`, `p_cv`}{Average detection probability over the
#'     searched strip, its standard error, and its coefficient of variation.}
#'   \item{`esw`}{Effective strip half-width: `p` times the width of the
#'     searched strip (`truncation` minus `left`).}
#'   \item{`cvm_p`}{Cramer-von Mises goodness-of-fit p-value. `NA` for binned
#'     distances.}
#'   \item{`chisq_p`}{Chi-square goodness-of-fit p-value. The test to read for
#'     binned distances.}
#' }
#' Columns describing a fit are `NA` for a candidate that failed.
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

# Build the selection table from a model set and its fits.
# Inputs: the model set, and the list of ddf fits (NULL for failures) in the
# same order. Returns a tibble with one row per candidate, fitted or not,
# sorted by AIC with failures last. Columns are listed in ?selection_table.
build_selection_table <- function(models, fits) {
  # 1. Pull AIC and the summary statistics out of each fit. `fits` is named
  #    by model_id and vapply() would copy those names onto every column, so
  #    they are removed here once.
  aic <- unname(vapply(fits, function(f) if (is.null(f)) NA_real_
                       else f$criterion, numeric(1)))
  stats_list <- lapply(fits, fit_summary)
  stat <- function(what) unname(vapply(stats_list, function(s) s[[what]],
                                       numeric(1)))

  # 2. Start from the candidate description and add one column per statistic.
  out <- models[, c("model_id", "key", "adjustment", "order", "formula")]
  out$converged <- unname(!vapply(fits, is.null, logical(1)))
  out$n_par <- stat("n_par")
  out$aic <- aic
  # If every candidate failed there is no minimum; min() would warn and
  # return Inf.
  out$delta_aic <- if (all(is.na(aic))) aic else aic - min(aic, na.rm = TRUE)
  out$p <- stat("p")
  out$p_se <- stat("p_se")
  out$p_cv <- out$p_se / out$p
  out$esw <- stat("esw")
  out$cvm_p <- stat("cvm_p")
  out$chisq_p <- stat("chisq_p")

  # 3. Sort by AIC, failed candidates last.
  out[order(out$aic, na.last = TRUE), , drop = FALSE]
}

# Summarise one fitted ddf object.
# Input: a ddf fit, or NULL for a failed candidate. Returns a list with
# n_par, p, p_se, esw, cvm_p and chisq_p. Every element is NA when the fit is
# NULL or its summary cannot be computed, so the table always gets one row
# of the same shape.
fit_summary <- function(f) {
  # 1. Return all-NA for a failed fit or a summary that errors.
  blank <- list(n_par = NA_real_, p = NA_real_, p_se = NA_real_,
                esw = NA_real_, cvm_p = NA_real_, chisq_p = NA_real_)
  if (is.null(f)) return(blank)

  s <- try(summary(f), silent = TRUE)
  if (inherits(s, "try-error")) return(blank)

  # 2. Goodness of fit.
  gof <- fit_gof(f)

  # 3. Average detection probability and esw. average.p is the mean of g(x)
  #    over the strip searched, which is (left, width) under left truncation,
  #    so esw multiplies it by width minus left.
  left <- if (is.null(s$left)) 0 else s$left
  p <- scalar(s$average.p)

  list(
    n_par = length(f$par),
    p = p,
    p_se = scalar(s$average.p.se),
    esw = p * (scalar(s$width) - left),
    cvm_p = gof$cvm_p,
    chisq_p = gof$chisq_p
  )
}

# Coerce an mrds summary field to a single number.
# Input: any value. Returns it as numeric if it has length 1, otherwise NA
# (mrds leaves some fields empty or absent).
scalar <- function(x) {
  if (length(x) != 1L) return(NA_real_)
  as.numeric(x)
}

# Compute goodness-of-fit p-values for one fit.
# Input: a ddf fit. Returns a list with cvm_p (Cramer-von Mises) and chisq_p
# (chi-square), each NA when mrds does not provide it or ddf.gof() errors.
# mrds computes Cramer-von Mises only for exact distances. Chi-square is
# always computed, over the survey's bins for binned data and over cutpoints
# mrds chooses for exact data. Both are kept under their own names so the
# reader knows which test a p-value came from.
fit_gof <- function(f) {
  out <- list(cvm_p = NA_real_, chisq_p = NA_real_)
  g <- try(mrds::ddf.gof(f, qq = FALSE), silent = TRUE)
  if (inherits(g, "try-error")) return(out)

  cvm <- g$dsgof$CvM$p
  if (!is.null(cvm)) out$cvm_p <- as.numeric(cvm)

  chisq <- g$chisquare$chi1$p
  if (!is.null(chisq)) out$chisq_p <- as.numeric(chisq)

  out
}

# Print the summary message after a sweep.
# Input: a dsfit_sweep. Emits a message (via rlang::inform) giving how many
# candidates fitted, the number of detections and truncation, rows dropped,
# which candidates failed, and a reminder that g(0) = 1 is assumed. Returns
# nothing useful.
report_sweep <- function(x) {
  # 1. Headline: fits, detections, truncation.
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
  # 2. Rows dropped during preparation, if any.
  if (x$settings$n_dropped > 0) {
    lines <- c(lines, paste0(
      "  ", x$settings$n_dropped,
      " row(s) dropped: no distance, or beyond the truncation",
      if (!is.null(x$settings$left)) ", or inside `left`" else "", "."
    ))
  }
  # 3. Candidates that failed.
  if (length(x$settings$failures)) {
    lines <- c(lines, paste0(
      "  did not fit: ",
      paste(names(x$settings$failures), collapse = ", ")
    ))
  }
  # 4. The g(0) reminder, then emit.
  lines <- c(lines, "  g(0) = 1 is assumed; see `?sweep_models`.")
  rlang::inform(paste(lines, collapse = "\n"))
}

# Choose the mrds monotonicity setting for one candidate, matching the
# Distance::ds() defaults.
# Input: one row of a model set. Returns list(mono, mono.strict), both TRUE
# for a candidate with adjustment terms, no covariates and a key other than
# gamma; both FALSE otherwise.
# Adjustment terms can bend g(x) above 1 or back upward, and without the
# constraint AIC and p differ from the same model fitted in Distance. A key
# alone is already monotone. The gamma key rises to a peak by design.
monotonicity <- function(model) {
  constrain <- !is.na(model$adjustment) && model$key != "gamma" &&
    identical(trimws(model$formula), "~1")
  list(mono = constrain, mono.strict = constrain)
}

# Get the first non-empty line of an error message.
# Input: a message string. Returns its first non-blank line, or the whole
# message if every line is blank. Some mrds errors start with a newline, so
# taking the text before the first newline would give "".
first_line <- function(msg) {
  lines <- trimws(strsplit(msg, "\n")[[1]])
  lines <- lines[nzchar(lines)]
  if (length(lines)) lines[1] else msg
}

# Check that `models` is a model set from model_set().
# Inputs: the object, and the calling environment for the error. Returns
# `models` invisibly, or errors if it is not a data frame with the model-set
# columns, or has no rows. Without this check a wrong object fails later with
# an error that does not name the argument.
check_model_set <- function(models, call = rlang::caller_env()) {
  need <- c("model_id", "key", "adjustment", "order", "formula", "dsmodel")
  if (!is.data.frame(models) || !all(need %in% names(models))) {
    rlang::abort("`models` must be a model set built by `model_set()`.",
                 call = call)
  }
  if (!nrow(models)) {
    rlang::abort("`models` has no candidates.", call = call)
  }
  invisible(models)
}

# Check that the mrds package is installed.
# Input: the calling environment for the error. Errors with an install hint
# if mrds is missing; returns nothing otherwise. mrds is in Suggests, not
# Imports, so it may be absent.
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
