# Combining g(0) components.
#
# g(0) is the probability of detecting an animal on the trackline. It is the
# product of separately measured components: availability (whether the animal
# is at the surface while in view; see R/availability.R), perception (whether
# an available animal is seen; needs a second observer team), and any others.
# g0() checks and combines them into a `dsfit_g0` object for the abundance
# step. dsfit does not apply the correction itself.
#
# Contents, in file order:
#   g0()             build a dsfit_g0 object from component tables
#   print.dsfit_g0() print method


#' Combine g(0) components into one correction
#'
#' Combines the separately measured components of `g(0)` (availability,
#' perception, and any others), multiplies them per key, and propagates their
#' standard errors. Use it to hand a documented correction to the abundance
#' step; `dsfit` does not apply the correction itself.
#'
#' @section Why the result is an object:
#' `dsfit` fits detection functions and does not compute abundance. The
#' correction is applied where density is calculated, outside this package.
#' A `dsfit_g0` object records which components went in, where each came
#' from, and its standard error, so that what reaches the abundance step is
#' not a bare number of unknown origin.
#'
#' @section Rules enforced:
#' These follow `docs/01-plan.md`.
#' \describe{
#'   \item{No silent default of 1}{At least one component is required. A
#'     missing `"availability"` or `"perception"` component is named in a
#'     warning and in the printout. Assuming a component is 1 while believing
#'     it was corrected for biases abundance by a factor.}
#'   \item{A standard error is required}{A component with an `NA` standard
#'     error is an error. The CV of `g(0)` often dominates the CV of
#'     abundance, so leaving it out understates the uncertainty.}
#'   \item{Components are kept separate}{Availability, perception and the
#'     geometric blind spot are measured in different ways. Each is kept as
#'     its own rows, so it stays visible which have been applied.}
#' }
#'
#' @section How it works:
#' Components are assumed independent, because they come from separate
#' studies (dive data for availability, double-observer trials for
#' perception). Under independence the delta method gives the coefficient of
#' variation (CV) of the product as
#'
#' \deqn{CV(g_0)^2 = \sum_i CV(x_i)^2}
#'
#' so the least precise component sets the floor: a perception CV of 30%
#' cannot be offset by an availability CV of 2%.
#'
#' 1. Each component table is checked for the columns `key`, `component`,
#'    `value` and `se`, and at least one row. A missing `source` column is
#'    filled with `NA`.
#' 2. The tables are stacked. Every row must have a component name, a
#'    standard error, and a value in (0, 1].
#' 3. Keys are checked (see below).
#' 4. A warning is given for each of `"availability"` and `"perception"`
#'    that is absent.
#' 5. For each key, the matching rows are multiplied to give `value`, and
#'    \eqn{CV = \sqrt{\sum (se/value)^2}} gives `cv`; `se` is
#'    `value * cv`.
#'
#' @section Keys:
#' Components are matched on `key` (for example a month). A component given as
#' a single row with `key = NA` applies to every key, such as one perception
#' estimate combined with monthly availability. Otherwise each component must
#' have exactly one row for each key in the full set of keys.
#'
#' Components keyed on different things are refused, not joined. Availability
#' by month and perception by year have no correct join: the answer needs one
#' value per month-year, so build that cross product and key both components
#' on it. Also refused: several rows with no keys, a mix of keyed and unkeyed
#' rows, and duplicate keys within a component.
#'
#' @section Recording where a component came from:
#' A component may carry a `source` column, which is kept in the object and
#' printed. A component without one is printed as `source not recorded`.
#'
#' This matters most for perception. Estimating it needs two independent
#' observer teams, and whether a programme ran them varies by programme, not
#' by the archive (such as NARWC) that holds its data. Perception is often
#' borrowed from another programme. Roberts et al. (2024) estimated it only
#' from NOAA's AMAPPS surveys, applied it to data from ten other
#' institutions, and noted this could bias their density estimates. Recording
#' the source keeps a borrowed correction from looking local.
#'
#' @param ... One or more data frames of components, each with columns `key`
#'   (character, or `NA` to apply to every key), `component` (name, such as
#'   `"availability"`), `value` (probability in (0, 1]) and `se` (standard
#'   error, finite and not negative), and optionally `source` (character).
#'   [availability()] returns this shape.
#'
#' @return An object of class `dsfit_g0`: a list with
#'   \describe{
#'     \item{table}{A tibble with one row per key: `key`, `value` (the
#'       combined `g(0)`), `se` and `cv`. If no component is keyed there is one
#'       row with `key = NA`.}
#'     \item{components}{A tibble of the input rows (`key`, `component`,
#'       `value`, `se`, `source`), kept so the correction can be taken apart
#'       again.}
#'   }
#'
#' @references
#' Laake, J.L. and Borchers, D.L. (2004) Methods for incomplete detection at
#' distance zero. In *Advanced Distance Sampling*, pp. 108-189. Oxford
#' University Press. Why the components are separate quantities.
#'
#' Roberts, J.J., Yack, T.M., Fujioka, E., Halpin, P.N., Baumgartner, M.F. and
#' others (2024) North Atlantic right whale density surface model for the US
#' Atlantic evaluated with passive acoustic monitoring. *Marine Ecology Progress
#' Series* 732:167-192. \doi{10.3354/meps14547} Perception borrowed across
#' programmes.
#'
#' @seealso [availability()] for the availability component, [sweep_models()]
#'   for the detection function fit, which this does not enter.
#'
#' @examples
#' avail <- availability(surface = 60, dive = 240, window = 24,
#'                       se_surface = 8, se_dive = 25)
#' avail$source <- "focal follows, this survey"
#'
#' # A perception estimate borrowed from another programme, recorded as such
#' perception <- data.frame(
#'   key = NA_character_, component = "perception", value = 0.68, se = 0.09,
#'   source = "AMAPPS double-observer, borrowed - not measured on this survey"
#' )
#'
#' g0(avail, perception)
#'
#' # Without perception: a warning, and "assumed 1" in the printout
#' suppressWarnings(g0(avail))
#'
#' @export
g0 <- function(...) {
  # Rule one: at least one component, with no default.
  parts <- list(...)
  if (!length(parts)) {
    rlang::abort(paste0(
      "A g(0) correction needs at least one component. There is no default: ",
      "a silent 1 is the assumption this object exists to prevent."
    ))
  }

  # Check each component table and reduce it to the standard columns.
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
    if (!nrow(p)) {
      rlang::abort(paste0("Component ", i, " has no rows."))
    }
    # `source` is optional. Absent means unrecorded, not local.
    p$source <- if ("source" %in% names(p)) as.character(p$source) else
      NA_character_
    if (!is.numeric(p$value) || !is.numeric(p$se)) {
      rlang::abort(paste0(
        "Component ", i, ": `value` and `se` must be numeric."
      ))
    }
    # Standardise columns and types. An all-NA key column becomes character
    # NA; an empty-string key is refused because NA already means "every key".
    p <- p[, c(cols, "source"), drop = FALSE]
    p$key <- if (all(is.na(p$key))) NA_character_ else as.character(p$key)
    if (any(!is.na(p$key) & !nzchar(trimws(p$key)))) {
      rlang::abort(paste0(
        "Component ", i, " has an empty `key`. Use `NA` for a row that ",
        "applies to every key."
      ))
    }
    p$component <- as.character(p$component)
    p
  })

  # Stack all components into one table.
  comp <- do.call(rbind, lapply(parts, as.data.frame))
  comp <- tibble::as_tibble(comp)

  # Rule three: every row names its component.
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
  # Values are probabilities; standard errors are finite and not negative.
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
  if (any(comp$se < 0) || any(!is.finite(comp$se))) {
    rlang::abort("Standard errors must be finite and not negative.")
  }

  # The full set of keys across all components. If none are keyed, there is
  # a single unkeyed result.
  keys <- unique(stats::na.omit(comp$key))
  if (!length(keys)) keys <- NA_character_

  # Check each component's keys. A component is either one unkeyed row that
  # applies everywhere, or exactly one row per key in the full set. Components
  # keyed on different things have no correct join, so they are refused.
  for (nm in unique(comp$component)) {
    rows <- comp[comp$component == nm, , drop = FALSE]
    if (nrow(rows) == 1L && is.na(rows$key[1])) next
    # Several rows and no keys, e.g. monthly availability built without
    # `key`: nothing says which row is which.
    if (all(is.na(rows$key))) {
      rlang::abort(paste0(
        "Component \"", nm, "\" has ", nrow(rows), " rows and no `key`. ",
        "Give each row a key (a month, a platform), or give one row with no ",
        "key to apply everywhere."
      ))
    }
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

  # Warn, by name, for each standard component that is absent.
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

  # Combine per key: unkeyed rows apply to every key. Values multiply.
  rows <- lapply(keys, function(k) {
    part <- comp[is.na(comp$key) | comp$key == k, , drop = FALSE]
    value <- prod(part$value)
    # Independent components, so the squared CVs add.
    cv <- sqrt(sum((part$se / part$value)^2))
    tibble::tibble(key = k, value = value, se = value * cv, cv = cv)
  })

  # Keep the input rows alongside the result so it can be taken apart.
  out <- list(
    table = do.call(rbind, rows),
    components = comp
  )
  class(out) <- "dsfit_g0"
  out
}


# Print method for dsfit_g0. Shows, in order: components assumed to be 1,
# the source of each component, the per-key table (first 12 rows), and a
# reminder that the correction is applied at the abundance step.
# Returns `x` invisibly.
#' @export
print.dsfit_g0 <- function(x, ...) {
  cat("<dsfit_g0>\n")

  # Standard components that are absent, and so effectively 1.

  missing <- setdiff(c("availability", "perception"),
                     unique(x$components$component))
  if (length(missing)) {
    cat("  assumed 1:   ", paste(missing, collapse = ", "), "\n", sep = "")
  }

  # Where each component came from, printed every time so a borrowed
  # correction does not look local. Unrecorded sources are said so.
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

  # The combined correction per key, rounded for display.
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

  # How to use it.
  cat("\n  Divide abundance by g0, and propagate cv. This package does not\n")
  cat("  apply it: correction happens at the abundance step.\n")
  invisible(x)
}
