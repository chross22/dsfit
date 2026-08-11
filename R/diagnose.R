#' Diagnose common reasons a sweep fails, or succeeds meaninglessly
#'
#' Runs the input guards, the truncation, and the model-set expansion — never
#' `mrds::ddf()` itself — and reports the most common ways a sweep goes wrong
#' before it reaches the fitting: a truncation that throws away most of the
#' survey, a covariate formula naming a column that is not there, a model set
#' that double-counts the blind spot, or too few detections left to fit
#' anything worth reading.
#'
#' Meant to run *before* [sweep_models()], so a misconfiguration is caught in
#' seconds rather than after a sweep that either errors from inside `mrds` or
#' returns a table that looks fine and is not.
#'
#' Every check here reports a problem rather than fixing it. This function never
#' modifies the data, the model set, or anything else.
#'
#' @param data A data frame with one row per detection, as [sweep_models()]
#'   takes.
#' @param models A [model_set()], or `NULL` for the default set.
#' @param truncation Right truncation distance. Required, as it is for a sweep.
#' @param left Left truncation distance, or `NULL`.
#' @param breaks Bin cutpoints for interval data, or `NULL` to derive them.
#'
#' @return Invisibly, `list(structure, prepared, models)` — whichever were
#'   reached before a fatal problem stopped the checks, so investigation can
#'   pick up from there.
#'
#' @seealso [detection_structure()], which this uses and which answers the
#'   different question of what the data could support in principle.
#'   [sweep_models()], which this is meant to run ahead of.
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
  ok <- TRUE
  header <- function(x) cat("\n== ", x, " ==\n", sep = "")
  pass <- function(...) cat("  ok    ", ..., "\n", sep = "")
  # A fact about how the survey was flown, which no amount of reconfiguring
  # changes. Reported, but it does not make the run "a problem" - otherwise
  # every single-observer survey ends in warnings and the reader learns to skim.
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

  if (missing(truncation)) {
    rlang::abort(paste0(
      "`truncation` is required here for the same reason it is required by ",
      "`sweep_models()`: it defines the data the models would be fitted to."
    ))
  }

  # --- Toolchain ----------------------------------------------------------
  header("Toolchain")
  if (requireNamespace("mrds", quietly = TRUE)) {
    pass("mrds ", as.character(utils::packageVersion("mrds")), " is installed")
  } else {
    warn("mrds is not installed, so fitting will fail even if everything ",
         "below passes. install.packages(\"mrds\")")
  }

  # --- What the data is ---------------------------------------------------
  header("Data")
  struct <- tryCatch(detection_structure(data), error = function(e) {
    fail("could not read the detections: ", conditionMessage(e))
    NULL
  })
  if (is.null(struct)) {
    cat("\nStopped: fix the error above before continuing.\n")
    return(invisible(list(structure = NULL)))
  }
  pass(struct$summary$n, " rows, ", struct$summary$n_distances, " ",
       struct$summary$type, " distances")

  # Absent structure is noted; ambiguous structure is warned about. A survey
  # that simply had one observer team is not misconfigured, and saying so as a
  # WARN on every dataset would bury the cases that are. A half-present
  # double-observer structure is the opposite: it is the shape that gets
  # mistaken for the real thing.
  reportable <- c("estimate perception bias", "carry group size to abundance")
  for (i in seq_len(nrow(struct$table))) {
    row <- struct$table[i, ]
    if (isTRUE(row$supported)) next
    short <- sub("[.].*$", "", row$detail)
    if (is.na(row$supported)) {
      warn(row$check, ": ", short)
    } else if (row$check %in% reportable) {
      note(row$check, ": ", short)
    }
  }

  # --- The guards, and what truncation costs -------------------------------
  header("Truncation")
  prepared <- tryCatch(
    prepare_distance_data(data, truncation = truncation, left = left,
                          breaks = breaks),
    error = function(e) {
      fail(conditionMessage(e))
      NULL
    }
  )
  if (is.null(prepared)) {
    cat("\nStopped: the data would not reach `mrds`. Fix the above.\n")
    return(invisible(list(structure = struct, prepared = NULL)))
  }

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

  # Losing most of the survey to truncation is legal, and rarely intended.
  lost <- 1 - n_kept / max(n_in, 1L)
  if (lost > 0.5) {
    warn("truncation drops ", round(100 * lost), "% of the rows. Buckland et ",
         "al. suggest trimming roughly the outer 5%; check `truncation` is in ",
         "the units the distances are in")
  }

  # --- Enough to fit ------------------------------------------------------
  if (n_kept < 40) {
    fail(n_kept, " detections is below the 40 usually workable; a fitted ",
         "detection function here would be shape rather than evidence")
  } else if (n_kept < 60) {
    warn(n_kept, " detections is workable but below the 60-80 usually ",
         "suggested")
  } else {
    pass(n_kept, " detections is at or above the 60-80 usually suggested")
  }

  # --- The model set ------------------------------------------------------
  header("Model set")
  models <- tryCatch(
    if (is.null(models)) model_set() else models,
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

  # The blind spot has exactly two honest treatments, and they are
  # alternatives. Checked here rather than earlier because whether the model
  # set already handles it decides whether there is anything to say.
  has_gamma <- any(models$key == "gamma")
  near <- struct$summary$nearest

  # Any continuous distance has a smallest value, so "nearest > 0" is not a
  # blind spot - it would fire on every survey ever flown. What marks one is an
  # empty strip wide enough to matter against the truncation. Two percent is a
  # heuristic and is stated as one below.
  blind <- !is.na(near) && is.finite(near) && near > 0.02 * truncation

  if (has_gamma && !is.null(left)) {
    warn("the set contains the gamma key and `left` is set. Both describe a ",
         "platform that cannot see beneath itself, so together they count the ",
         "blind spot twice")
  } else if (blind && !is.null(left)) {
    pass("empty strip to ", signif(near, 4), " removed by left truncation at ",
         left)
  } else if (blind) {
    abrupt <- identical(onset_shape(prepared$data, near, truncation), "abrupt")
    where <- paste0("no detections inside ", signif(near, 4), ", which is ",
                    round(100 * near / truncation),
                    "% of the truncation width")

    if (abrupt) {
      # The distinction gamma cannot cross. A geometric cutoff is a
      # discontinuity, and every key function is continuous.
      warn(where, ", and detections begin at close to their peak rate rather ",
           "than rising into it. That is the shape of a geometric edge, which ",
           "`left` removes and no key function reproduces",
           if (has_gamma) paste0(
             " - including gamma, which models a gradual reduction toward the ",
             "trackline rather than an edge. Expect it to fail its ",
             "goodness-of-fit test too") else "")
    } else if (has_gamma) {
      pass("detections rise into the empty strip to ", signif(near, 4),
           " rather than starting abruptly, which is the shape the gamma key ",
           "is for")
    } else {
      warn(where, ", and neither `left` nor the gamma key is in use")
    }
  }

  # A formula naming a column that is not there fails inside mrds, with a
  # message that does not name the column.
  wanted <- unique(unlist(lapply(models$formula, formula_vars)))
  wanted <- setdiff(wanted, "")
  absent <- setdiff(wanted, names(data))
  if (length(absent)) {
    fail("covariate", if (length(absent) > 1) "s" else "", " not in `data`: ",
         paste(absent, collapse = ", "))
  } else if (length(wanted)) {
    pass("covariates present: ", paste(wanted, collapse = ", "))
    constant <- wanted[vapply(data[wanted], function(x)
      length(unique(x[!is.na(x)])) < 2L, logical(1))]
    if (length(constant)) {
      warn("covariate", if (length(constant) > 1) "s" else "", " with a single ",
           "value, which cannot inform detection: ",
           paste(constant, collapse = ", "))
    }
    missing_cov <- wanted[vapply(data[wanted], anyNA, logical(1))]
    if (length(missing_cov)) {
      warn("covariate", if (length(missing_cov) > 1) "s" else "",
           " with missing values, which mrds drops rows for: ",
           paste(missing_cov, collapse = ", "))
    }
  }

  # --- The assumption no check can make for you ----------------------------
  header("What is still assumed")
  cat("  g(0) = 1, unless a correction is applied at the abundance step.\n")
  cat("  Nothing here can detect a wrong one: it scales every candidate\n")
  cat("  equally, so the ranking looks untouched. See ?g0.\n")

  cat("\n", if (ok) "No problems found. Nothing was fitted." else
      "Problems above. Nothing was fitted.", "\n", sep = "")

  invisible(list(structure = struct, prepared = prepared, models = models))
}


# Does the near edge of the distance distribution begin abruptly, or rise into
# itself? The two have different treatments and the difference is visible
# without fitting anything.
#
# A geometric cutoff - an aircraft's flat windows, a hull blocking the view -
# is a discontinuity: detection is zero, then immediately whatever it would
# have been. So the first bin above the edge already sits near the modal rate.
# A genuine decline in detectability toward the trackline instead ramps up, and
# that is the unimodal shape the gamma key exists to fit.
#
# It matters because no key function is discontinuous. Gamma models a gradual
# reduction, not an edge, so against a hard cutoff it misfits like the others.
onset_shape <- function(data, near, truncation, bins = 20L) {
  x <- if ("distance" %in% names(data)) data$distance else data$distbegin
  x <- x[!is.na(x)]
  if (length(x) < 40L || !is.finite(near) || near >= truncation) {
    return("unknown")
  }

  edges <- seq(near, truncation, length.out = bins + 1L)
  counts <- table(cut(x, breaks = edges, include.lowest = TRUE))
  if (!length(counts) || max(counts) == 0) return("unknown")

  if (counts[[1]] / max(counts) > 0.6) "abrupt" else "gradual"
}


# Variable names on the right of a one-sided formula string like "~bf + size".
formula_vars <- function(x) {
  f <- try(stats::as.formula(x), silent = TRUE)
  if (inherits(f, "try-error")) return(character(0))
  all.vars(f)
}
