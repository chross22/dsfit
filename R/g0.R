#' Assemble a g(0) correction from its named components
#'
#' Stacks the components of `g(0)` — availability, perception, anything else
#' measured separately — multiplies them, and propagates their variance. The
#' result is an object to hand to the abundance step, not a number applied here.
#'
#' @section Why this is an object and not an argument:
#' `dsfit` fits detection functions; it does not compute abundance. The
#' correction is applied one layer out, where density is calculated, and that
#' layer is a `targets` repository rather than a package.
#'
#' What this can do from here is make the correction impossible to get wrong
#' quietly. A `dsfit_g0` cannot be built without naming its components, cannot
#' be built without standard errors, and cannot be built at all by accident —
#' so whatever consumes it downstream is handed a documented object rather than
#' a bare number of unknown provenance.
#'
#' @section The three rules:
#' From `docs/01-plan.md`, and each one is enforced rather than documented:
#'
#' \describe{
#'   \item{Never silently 1}{There is no default. A component that is absent is
#'     warned about by name, because assuming perception is 1 while believing
#'     you have corrected for it is how these estimates go wrong by a factor.}
#'   \item{Propagate the variance or refuse}{A component without a standard
#'     error is an error, not a point estimate. The CV of `g(0)` routinely
#'     dominates the CV of abundance, so dropping it inverts which uncertainty
#'     matters.}
#'   \item{Name the components separately}{Availability, perception and the
#'     geometric blind spot are three different things measured three different
#'     ways. They are kept as rows, not collapsed, so it stays visible which
#'     have been applied.}
#' }
#'
#' @section How the components combine:
#' Multiplicatively, and independently — availability comes from dive data and
#' perception from a double-observer trial, so they are estimated from separate
#' studies and the covariance between them is not merely unknown but usually
#' undefined.
#'
#' Under independence the delta method gives a result worth remembering: the
#' squared CVs add.
#'
#' \deqn{CV(g_0)^2 = \sum_i CV(x_i)^2}
#'
#' So the least precise component sets the floor. A perception estimate with a
#' 30% CV cannot be rescued by an availability estimate with a 2% one.
#'
#' @section Keys, and the trap in them:
#' Components are matched on `key`. A component with a single `NA` key applies
#' to every key — a perception estimate that does not vary by month, say,
#' combined with availability that does.
#'
#' Components keyed on **different things** are refused rather than joined.
#' Availability by month and perception by year is the common case, and it has
#' no correct join: the answer is a value per month-year, which means building
#' the cross product yourself and keying on it. Silently pairing January with
#' 1998 because both are first in their vectors is exactly the failure this
#' refuses to commit.
#'
#' @section Recording where a component came from:
#' Components may carry an optional `source` column, and it is worth using.
#'
#' Perception is the reason. Estimating it needs two independent observer teams,
#' and whether a survey programme ran them is a property of that programme
#' rather than of the archive its data ends up in. A NARWC extract may or may
#' not carry the structure, so for many datasets the only available perception
#' estimate is one **borrowed from a different programme**. Roberts et al.
#' (2024) did exactly this — perception was estimable only from NOAA's AMAPPS
#' surveys, and they applied those corrections to the other ten institutions'
#' data, cautioning in print that it could have biased their density estimates.
#'
#' A borrowed correction that looks local is the failure this guards against. So
#' `source` travels with the component into the object and is printed every
#' time, and a component with none is shown as `source not recorded` rather than
#' silently blank.
#'
#' @param ... One or more component tables, each with `key`, `component`,
#'   `value` and `se` columns, and optionally `source`. [availability()] returns
#'   this shape.
#'
#' @return An object of class `dsfit_g0`: a list with `table` (one row per key,
#'   carrying `value`, `se` and `cv`) and `components` (the rows it was built
#'   from, `source` included, kept so the correction can always be taken apart
#'   again).
#'
#' @references
#' Laake, J.L. and Borchers, D.L. (2004) Methods for incomplete detection at
#' distance zero. In *Advanced Distance Sampling*, pp. 108-189. Oxford
#' University Press. Why the components are separate quantities.
#'
#' @seealso [availability()] for the availability component. [sweep_models()]
#'   for why none of this enters the detection function fit.
#'
#' @examples
#' avail <- availability(surface = 60, dive = 240, window = 24,
#'                       se_surface = 8, se_dive = 25)
#' avail$source <- "focal follows, this survey"
#'
#' # A perception estimate that is not this survey's, said so
#' perception <- data.frame(
#'   key = NA_character_, component = "perception", value = 0.68, se = 0.09,
#'   source = "AMAPPS double-observer, borrowed - not measured on this survey"
#' )
#'
#' g0(avail, perception)
#'
#' # Without perception at all, which is said too
#' suppressWarnings(g0(avail))
#'
#' @export
g0 <- function(...) {
  parts <- list(...)
  if (!length(parts)) {
    rlang::abort(paste0(
      "A g(0) correction needs at least one component. There is no default: ",
      "a silent 1 is the assumption this object exists to prevent."
    ))
  }

  cols <- c("key", "component", "value", "se")
  parts <- lapply(seq_along(parts), function(i) {
    p <- parts[[i]]
    if (!is.data.frame(p)) {
      rlang::abort(paste0(
        "Component ", i, " is not a data frame. Each component needs `key`, ",
        "`component`, `value` and `se` columns - the shape `availability()` ",
        "returns."
      ))
    }
    missing <- setdiff(cols, names(p))
    if (length(missing)) {
      rlang::abort(paste0(
        "Component ", i, " is missing: ", paste(missing, collapse = ", "), "."
      ))
    }
    # `source` is optional, and absent means unrecorded rather than local.
    p$source <- if ("source" %in% names(p)) as.character(p$source) else
      NA_character_
    p <- p[, c(cols, "source"), drop = FALSE]
    p$key <- if (all(is.na(p$key))) NA_character_ else as.character(p$key)
    p$component <- as.character(p$component)
    p
  })

  comp <- do.call(rbind, lapply(parts, as.data.frame))
  comp <- tibble::as_tibble(comp)

  if (anyNA(comp$component) || any(!nzchar(trimws(comp$component)))) {
    rlang::abort(paste0(
      "Every component needs a name. Availability, perception and the ",
      "geometric blind spot are measured three different ways, and a ",
      "correction that does not say which it contains cannot be checked."
    ))
  }

  # Rule two: propagate the variance or refuse the correction.
  bad_se <- is.na(comp$se)
  if (any(bad_se)) {
    rlang::abort(paste0(
      "No standard error for: ",
      paste(unique(comp$component[bad_se]), collapse = ", "),
      ". The CV of g(0) routinely dominates the CV of abundance, so a ",
      "correction without one reports a precision it does not have. Supply ",
      "`se`, or do not apply the correction."
    ))
  }
  if (anyNA(comp$value)) {
    rlang::abort("Every component needs a `value`; some are missing.")
  }
  if (any(comp$value <= 0 | comp$value > 1)) {
    rlang::abort(paste0(
      "Component values are probabilities and must be in (0, 1]. A value ",
      "above 1 would be a correction that finds more animals than were ",
      "available to be seen."
    ))
  }
  if (any(comp$se < 0)) {
    rlang::abort("Standard errors cannot be negative.")
  }

  keys <- unique(stats::na.omit(comp$key))
  if (!length(keys)) keys <- NA_character_

  # A component keyed on one thing and another keyed on something else have no
  # correct join, so this refuses rather than guessing.
  for (nm in unique(comp$component)) {
    rows <- comp[comp$component == nm, , drop = FALSE]
    if (nrow(rows) == 1L && is.na(rows$key[1])) next
    if (anyNA(rows$key)) {
      rlang::abort(paste0(
        "Component \"", nm, "\" mixes keyed and unkeyed rows. Give it one row ",
        "with no key, to apply everywhere, or one row per key."
      ))
    }
    if (anyDuplicated(rows$key)) {
      rlang::abort(paste0(
        "Component \"", nm, "\" has more than one row for the same key."
      ))
    }
    if (!setequal(rows$key, keys)) {
      rlang::abort(paste0(
        "Component \"", nm, "\" is keyed on something different from the ",
        "others: it has ", nrow(rows), " key(s) against ", length(keys),
        ". Availability by month and perception by year is the usual way to ",
        "arrive here, and it has no correct join - the answer is a value per ",
        "month-year. Build that cross product and key both components on it."
      ))
    }
  }

  named <- unique(comp$component)
  for (expected in c("availability", "perception")) {
    if (!expected %in% named) {
      rlang::warn(paste0(
        "No \"", expected, "\" component: it is being left at 1. If that is ",
        "deliberate, say so where this correction is used - correcting for ",
        "one component while believing you have corrected for another is how ",
        "these estimates go wrong by a factor rather than a percentage."
      ))
    }
  }

  rows <- lapply(keys, function(k) {
    part <- comp[is.na(comp$key) | comp$key == k, , drop = FALSE]
    value <- prod(part$value)
    # Independent components, so the squared CVs add.
    cv <- sqrt(sum((part$se / part$value)^2))
    tibble::tibble(key = k, value = value, se = value * cv, cv = cv)
  })

  out <- list(
    table = do.call(rbind, rows),
    components = comp
  )
  class(out) <- "dsfit_g0"
  out
}


#' @export
print.dsfit_g0 <- function(x, ...) {
  cat("<dsfit_g0>\n")

  missing <- setdiff(c("availability", "perception"),
                     unique(x$components$component))
  if (length(missing)) {
    cat("  assumed 1:   ", paste(missing, collapse = ", "), "\n", sep = "")
  }

  # Where each component came from. Printed every time, with anything
  # unrecorded said so, because a borrowed correction that looks local is the
  # mistake this is here to prevent.
  named <- unique(x$components$component)
  width <- max(nchar(named))
  cat("  components:\n")
  for (nm in named) {
    src <- unique(x$components$source[x$components$component == nm])
    src <- src[!is.na(src)]
    label <- if (!length(src)) "source not recorded" else
      paste(src, collapse = "; ")
    cat("    ", formatC(nm, width = -width), "  ", label, "\n", sep = "")
  }

  tab <- x$table
  cat("\n")
  show <- data.frame(
    key = if (all(is.na(tab$key))) "-" else tab$key,
    g0 = round(tab$value, 4),
    se = round(tab$se, 4),
    cv = round(tab$cv, 3)
  )
  print(utils::head(show, 12), row.names = FALSE)
  if (nrow(show) > 12) cat("  ... and ", nrow(show) - 12, " more\n", sep = "")

  cat("\n  Divide abundance by g0, and propagate cv. This package does not\n")
  cat("  apply it: correction happens at the abundance step.\n")
  invisible(x)
}
