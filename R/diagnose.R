#' Check a sweep's inputs before fitting
#'
#' Runs the checks a sweep would run, without calling `mrds::ddf()`, and prints
#' a report of problems found. Use it before [sweep_models()] to catch a
#' misconfiguration in seconds, instead of after a sweep that errors inside
#' `mrds` or returns a table that looks valid but is not.
#'
#' The function only reports. It never changes the data or the model set.
#'
#' @section How it works:
#' Each finding is printed as one line marked `ok`, `note`, `WARN` or `FAIL`.
#' A `note` is a fact about how the survey was flown (for example, a single
#' observer team); it does not count as a problem. `WARN` and `FAIL` lines
#' count as problems. A failure marked "stops" below ends the report early.
#'
#' 1. Refuse to run if `truncation` is missing (an error, not a report line).
#' 2. Toolchain: check `mrds` is installed (`WARN` if not).
#' 3. Data: read the table with [detection_structure()].
#'    * `FAIL` and stop if it cannot be read.
#'    * `FAIL` if it has no distances; otherwise report the count and type.
#'    * `WARN` for each pair of columns that name the same variable.
#'    * `WARN` for each analysis whose structure is incomplete (`supported`
#'      is `NA`), such as a partial double-observer setup.
#'    * `note` if perception bias cannot be estimated or there is no `size`
#'      column.
#' 4. Truncation: run [prepare_distance_data()].
#'    * Its warnings are printed as `WARN` lines. Its errors are printed as a
#'      `FAIL` line and stop the report.
#'    * Report how many rows were kept, and how many were dropped for each
#'      reason.
#'    * `WARN` if truncation removes more than 50% of the rows that have a
#'      distance.
#' 5. Sample size: `FAIL` below 40 detections after truncation, `WARN`
#'    below 60, `ok` otherwise (Buckland et al. 2001).
#' 6. Model set: build the default [model_set()] or check the one given.
#'    `FAIL` and stop if that errors; otherwise list the candidates.
#' 7. Blind spot: see the section below.
#' 8. Covariates: `FAIL` if a formula names a column not in `data`.
#'    Otherwise `WARN` for a covariate with only one value, and for one with
#'    missing values (mrds drops those rows).
#' 9. Print a reminder that g(0) = 1 is assumed unless corrected at the
#'    abundance step. g(0) is the detection probability on the trackline. An
#'    error in it scales every candidate equally, so no check here can detect
#'    it. See [g0()].
#' 10. Print a final line: no problems found, or problems above.
#'
#' The constant and missing-value covariate checks use the rows that will be
#' fitted, after truncation, since mrds never sees the rest.
#'
#' @section Blind-spot check:
#' A *blind spot* is a strip next to the trackline where the platform cannot
#' see, such as beneath an aircraft. It shows as no detections close to the
#' line.
#'
#' The check fires when the nearest detection in the input data is more than
#' 2% of `truncation` from the line. Every dataset has a smallest distance, so
#' a nonzero nearest distance alone is not evidence of a blind spot. The 2%
#' threshold is a heuristic.
#'
#' There are two treatments, and they are alternatives: left truncation
#' (`left`), or the gamma key, which models detection rising from zero at the
#' line. The check reports:
#'
#' | Situation | Result |
#' |---|---|
#' | Gamma key in the set and `left` set (blind spot or not) | `WARN`: the blind spot is counted twice. |
#' | Blind spot, `left` set, no gamma key | `ok` |
#' | Blind spot, abrupt onset | `WARN`: a hard edge, which `left` removes and no key function fits, gamma included. |
#' | Blind spot, gradual onset, gamma key in set | `ok` |
#' | Blind spot, too few detections (under 40) to classify the onset, gamma key in set | `note` |
#' | Blind spot, neither `left` nor gamma | `WARN` |
#'
#' *Abrupt onset* means detections start at close to their peak rate just
#' beyond the empty strip. The range from the nearest detection to
#' `truncation` is cut into 20 equal bins. The onset is abrupt when the first
#' bin's count is more than 60% of the largest bin's count. A geometric edge
#' (a window frame, a hull) gives this shape. A gradual decline in
#' detectability toward the line instead gives a rising first few bins, which
#' is what the gamma key fits. The test needs at least 40 detections; with
#' fewer, the onset is not classed as abrupt.
#'
#' @param data Data frame, one row per detection, in the form
#'   [sweep_models()] takes.
#' @param models A [model_set()], or `NULL` (default) for the default set.
#' @param truncation Number. Right truncation distance, in the units of the
#'   distances. Required, as for a sweep.
#' @param left Number or `NULL` (default). Left truncation distance.
#' @param breaks Numeric vector or `NULL` (default). Bin cutpoints for binned
#'   data; derived from the bins when `NULL`.
#'
#' @return Invisibly, a list with the results reached before any stop:
#'   \describe{
#'     \item{`structure`}{The [detection_structure()] result, or `NULL` if the
#'       data could not be read.}
#'     \item{`prepared`}{The [prepare_distance_data()] result, or `NULL` if it
#'       failed. Absent if the report stopped earlier.}
#'     \item{`models`}{The model set checked. Absent if the report stopped
#'       earlier.}
#'   }
#'   The report itself is printed to the console.
#'
#' @references
#' Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L.
#' and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford University
#' Press. Source of the 40 and 60-80 detection guidelines and of trimming
#' roughly the outer 5% of distances.
#'
#' @seealso [detection_structure()], which reports what the data could
#'   support in principle. [sweep_models()], which this is meant to run
#'   before.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   object = 1:300,
#'   distance = c(abs(rnorm(295, 0, 800)), rep(NA, 5)),
#'   beaufort = sample(0:4, 300, replace = TRUE)
#' )
#' diagnose_sweep(d, model_set(c("hn", "hr")), truncation = 2000)
#'
#' @export
diagnose_sweep <- function(data, models = NULL, truncation, left = NULL,
                           breaks = NULL) {
  # Report-line printers. `ok` becomes FALSE on any WARN or FAIL and decides
  # the final line.
  ok <- TRUE
  header <- function(x) cat("\n== ", x, " ==\n", sep = "")
  pass <- function(...) cat("  ok    ", ..., "\n", sep = "")
  # A note is a fact about how the survey was flown, which reconfiguring
  # cannot change, so it does not count as a problem. Otherwise every
  # single-observer survey would end in warnings and readers would learn to
  # skim them.
  note <- function(...) cat("  note  ", ..., "\n", sep = "")
  warn <- function(...) {
    cat("  WARN  ", ..., "\n", sep = "")
    ok <<- FALSE
  }
  fail <- function(...) {
    cat("  FAIL  ", ..., "\n", sep = "")
    ok <<- FALSE
  }

  cat("dsfit sweep diagnosis\n")

  # 1. Require `truncation`: it defines the data the models would be fitted to.
  if (missing(truncation)) {
    rlang::abort(paste0(
      "`truncation` is required here for the same reason it is required by ",
      "`sweep_models()`: it defines the data the models would be fitted to."
    ))
  }

  # 2. Toolchain ------------------------------------------------------------
  header("Toolchain")
  if (requireNamespace("mrds", quietly = TRUE)) {
    pass("mrds ", as.character(utils::packageVersion("mrds")), " is installed")
  } else {
    warn("mrds is not installed, so fitting will fail even if everything ",
         "below passes. install.packages(\"mrds\")")
  }

  # 3. Data -----------------------------------------------------------------
  header("Data")
  # 3a. Read the structure; stop if the table cannot be read at all.
  struct <- tryCatch(detection_structure(data), error = function(e) {
    fail("could not read the detections: ", conditionMessage(e))
    NULL
  })
  if (is.null(struct)) {
    cat("\nStopped: fix the error above before continuing.\n")
    return(invisible(list(structure = NULL)))
  }
  # 3b. Report the number and type of distances.
  if (identical(struct$summary$type, "none")) {
    fail(struct$summary$n, " rows, and no distances in them")
  } else {
    pass(struct$summary$n, " rows, ", struct$summary$n_distances, " ",
         struct$summary$type, " distances")
  }

  # 3c. Two columns naming one variable. This is ambiguous, not just
  # redundant: nothing here can tell which one a fit should use.
  for (k in names(struct$summary$duplicates)) {
    warn("two columns name ", k, ": ",
         paste(struct$summary$duplicates[[k]], collapse = " and "),
         ". Read the extract with narwcr::read_narwc(prefer_source = TRUE), ",
         "which decides between them and keeps the loser as ", k, "_ORIGINAL")
  }

  # 3d. Unsupported analyses. Incomplete structure (supported = NA) is a WARN,
  # because a half-present double-observer setup is easily mistaken for a
  # complete one. Absent structure is a note, and only for the two checks
  # below: a survey with one observer team is not misconfigured. Other
  # unsupported checks (e.g. availability) are not printed.
  reportable <- c("estimate perception bias", "carry group size to abundance")
  for (i in seq_len(nrow(struct$table))) {
    row <- struct$table[i, ]
    if (isTRUE(row$supported)) next
    # Print only the first sentence of the detail.
    short <- sub("[.].*$", "", row$detail)
    if (is.na(row$supported)) {
      warn(row$check, ": ", short)
    } else if (row$check %in% reportable) {
      note(row$check, ": ", short)
    }
  }

  # 4. Truncation -----------------------------------------------------------
  header("Truncation")
  # 4a. Run the same preparation as a sweep. Its warnings become WARN lines
  # here, instead of R warnings printed after the report. An error stops the
  # report.
  prepared <- tryCatch(
    withCallingHandlers(
      prepare_distance_data(data, truncation = truncation, left = left,
                            breaks = breaks),
      warning = function(w) {
        warn(conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      fail(conditionMessage(e))
      NULL
    }
  )
  if (is.null(prepared)) {
    cat("\nStopped: the data would not reach `mrds`. Fix the above.\n")
    return(invisible(list(structure = struct, prepared = NULL)))
  }

  # 4b. Report rows kept and rows dropped by reason.
  n_kept <- nrow(prepared$data)
  n_in <- nrow(data)
  pass(n_kept, " of ", n_in, " rows kept at truncation ", truncation,
       if (!is.null(left)) paste0(", left ", left) else "")

  d <- prepared$dropped
  if (sum(d) > 0) {
    parts <- c(
      if (d[["no_distance"]] > 0) paste0(d[["no_distance"]], " with no distance"),
      if (d[["beyond_truncation"]] > 0)
        paste0(d[["beyond_truncation"]], " beyond the truncation"),
      if (d[["inside_left"]] > 0) paste0(d[["inside_left"]], " inside `left`")
    )
    cat("        dropped: ", paste(parts, collapse = ", "), "\n", sep = "")
  }

  # 4c. Warn if truncation removes most of the survey. That is allowed but
  # rarely intended; a common cause is `truncation` in the wrong units. The
  # share is taken over rows that have a distance, because rows without one
  # are a flatfile's empty samples and truncation did not remove them.
  with_distance <- n_in - d[["no_distance"]]
  lost <- (d[["beyond_truncation"]] + d[["inside_left"]]) /
    max(with_distance, 1L)
  if (lost > 0.5) {
    warn("truncation drops ", round(100 * lost), "% of the detections. Buckland et ",
         "al. suggest trimming roughly the outer 5%; check `truncation` is in ",
         "the units the distances are in")
  }

  # 5. Sample size: 40 workable, 60-80 suggested (Buckland et al. 2001).
  # Printed under the Truncation header.
  if (n_kept < 40) {
    fail(n_kept, " detections is below the 40 usually workable; a fitted ",
         "detection function here would be shape rather than evidence")
  } else if (n_kept < 60) {
    warn(n_kept, " detections is workable but below the 60-80 usually ",
         "suggested")
  } else {
    pass(n_kept, " detections is at or above the 60-80 usually suggested")
  }

  # 6. Model set ------------------------------------------------------------
  header("Model set")
  # Build the default set or validate the given one; stop on error.
  models <- tryCatch(
    if (is.null(models)) model_set() else check_model_set(models),
    error = function(e) {
      fail("could not build the model set: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(models)) {
    cat("\nStopped: fix the model set above.\n")
    return(invisible(list(structure = struct, prepared = prepared)))
  }
  pass(nrow(models), " candidate", if (nrow(models) != 1) "s" else "", ": ",
       paste(models$model_id, collapse = ", "))

  # 7. Blind spot. Left truncation and the gamma key are alternative
  # treatments, so this check needs the model set and runs after it.
  has_gamma <- any(models$key == "gamma")
  nearest <- struct$summary$nearest

  # 7a. Decide whether there is a blind spot. Every dataset has a smallest
  # distance, so "nearest > 0" alone would fire on every survey. A blind spot
  # is an empty strip that is wide against the truncation; 2% is a heuristic.
  # `nearest` comes from the input data, before truncation.
  blind <- !is.na(nearest) && is.finite(nearest) &&
    nearest > 0.02 * truncation

  # 7b. Report, by treatment in use. Gamma with `left` is warned about even
  # without a detected blind spot.
  if (has_gamma && !is.null(left)) {
    warn("the set contains the gamma key and `left` is set. Both describe a ",
         "platform that cannot see beneath itself, so together they count the ",
         "blind spot twice")
  } else if (blind && !is.null(left)) {
    pass("empty strip to ", signif(nearest, 4), " removed by left truncation at ",
         left)
  } else if (blind) {
    # Classify the onset on the truncated data: "abrupt", "gradual", or
    # "unknown" when there are too few detections to tell.
    onset <- onset_shape(prepared$data, nearest, truncation)
    abrupt <- identical(onset, "abrupt")
    where <- paste0("no detections inside ", signif(nearest, 4), ", which is ",
                    round(100 * nearest / truncation),
                    "% of the truncation width")

    if (abrupt) {
      # A geometric cutoff is a discontinuity, and every key function is
      # continuous, so gamma cannot fit it either.
      warn(where, ", and detections begin at close to their peak rate rather ",
           "than rising into it. That is the shape of a geometric edge, which ",
           "`left` removes and no key function reproduces",
           if (has_gamma) paste0(
             " - including gamma, which models a gradual reduction toward the ",
             "trackline rather than an edge. Expect it to fail its ",
             "goodness-of-fit test too") else "")
    } else if (has_gamma && identical(onset, "unknown")) {
      # Too few detections to say whether gamma's shape fits the edge.
      note(where, ". The gamma key is in the set, but there are too few ",
           "detections to tell whether they rise into the strip (gamma's ",
           "shape) or start abruptly (an edge, which `left` handles)")
    } else if (has_gamma) {
      pass("detections rise into the empty strip to ", signif(nearest, 4),
           " rather than starting abruptly, which is the shape the gamma key ",
           "is for")
    } else {
      warn(where, ", and neither `left` nor the gamma key is in use")
    }
  }

  # 8. Covariates. A formula naming a missing column fails inside mrds with a
  # message that does not name the column, so check names here.
  wanted <- unique(unlist(lapply(models$formula, formula_vars)))
  wanted <- setdiff(wanted, "")
  absent <- setdiff(wanted, names(data))
  if (length(absent)) {
    fail("covariate", if (length(absent) > 1) "s" else "", " not in `data`: ",
         paste(absent, collapse = ", "))
  } else if (length(wanted)) {
    pass("covariates present: ", paste(wanted, collapse = ", "))
    # A covariate with one distinct value cannot explain variation in
    # detection.
    # Checked on the rows that will be fitted: rows with no distance, or
    # beyond the truncation, are dropped before mrds sees their covariates.
    kept <- prepared$data
    constant <- wanted[vapply(kept[wanted], function(x)
      length(unique(x[!is.na(x)])) < 2L, logical(1))]
    if (length(constant)) {
      warn("covariate", if (length(constant) > 1) "s" else "", " with a single ",
           "value, which cannot inform detection: ",
           paste(constant, collapse = ", "))
    }
    # mrds drops rows with a missing covariate, shrinking the sample silently.
    missing_cov <- wanted[vapply(kept[wanted], anyNA, logical(1))]
    if (length(missing_cov)) {
      warn("covariate", if (length(missing_cov) > 1) "s" else "",
           " with missing values, which mrds drops rows for: ",
           paste(missing_cov, collapse = ", "))
    }
  }

  # 9. The g(0) assumption, which no check on the data can test.
  header("What is still assumed")
  cat("  g(0) = 1, unless a correction is applied at the abundance step.\n")
  cat("  Nothing here can detect a wrong one: it scales every candidate\n")
  cat("  equally, so the ranking looks untouched. See ?g0.\n")

  # 10. Summary line and return.
  cat("\n", if (ok) "No problems found. Nothing was fitted." else
      "Problems above. Nothing was fitted.", "\n", sep = "")

  invisible(list(structure = struct, prepared = prepared, models = models))
}


# Classify the near edge of the distance distribution as "abrupt", "gradual"
# or "unknown".
# Inputs: prepared data (`distance`, or `distbegin` for binned data), the
# nearest detection distance `nearest`, the truncation, and the number of bins.
# Returns one of the three strings.
#
# Why: the two shapes need different treatments, and the difference is
# visible without fitting. A geometric cutoff (an aircraft's window frame, a
# hull blocking the view) is a discontinuity: detection is zero, then
# immediately near its usual rate, so the first bin past the edge is already
# near the modal count. A gradual decline in detectability toward the
# trackline ramps up instead, which is the shape the gamma key fits. No key
# function is discontinuous, so none fits a hard edge.
onset_shape <- function(data, nearest, truncation, bins = 20L) {
  # 1. Take the distances; for binned data, the bin starts.
  x <- if ("distance" %in% names(data)) data$distance else data$distbegin
  x <- x[!is.na(x)]
  # 2. Too few detections, or no usable edge: cannot classify.
  if (length(x) < 40L || !is.finite(nearest) || nearest >= truncation) {
    return("unknown")
  }

  # 3. Count detections in equal-width bins from the edge to the truncation.
  edges <- seq(nearest, truncation, length.out = bins + 1L)
  counts <- table(cut(x, breaks = edges, include.lowest = TRUE))
  if (!length(counts) || max(counts) == 0) return("unknown")

  # 4. Abrupt if the first bin holds more than 60% of the largest bin's count.
  if (counts[[1]] / max(counts) > 0.6) "abrupt" else "gradual"
}


# List the variable names in a one-sided formula string such as "~bf + size".
# Input: a formula string. Returns a character vector, empty if the string
# does not parse.
formula_vars <- function(x) {
  f <- try(stats::as.formula(x), silent = TRUE)
  if (inherits(f, "try-error")) return(character(0))
  all.vars(f)
}
