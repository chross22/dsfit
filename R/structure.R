#' Report which analyses a detection table supports
#'
#' Reads a table of detections and reports which analyses its structure
#' allows, and the reason for each answer. Nothing is fitted or estimated. Use
#' it on a new dataset to see what can be done with it before choosing models.
#'
#' @section How it works:
#' 1. Classify the distances as `"exact"` (a non-missing `distance`),
#'    `"binned"` (a non-missing `distbegin`, with `distend` present), `"both"`
#'    or `"none"`.
#' 2. Answer seven questions, one row each (table below).
#' 3. Find columns that name the same variable twice (see below).
#' 4. Record the nearest detection distance, where a blind spot beneath the
#'    platform would show.
#'
#' The seven checks, in the order they appear in `table`:
#'
#' | Check | Supported when |
#' |---|---|
#' | fit a detection function | Exact or binned distances. `NA` if both are present. The detail compares the count with the 40 workable and 60-80 suggested detections (Buckland et al. 2001). |
#' | test fit with Cramer-von Mises | Exact distances only. Binned data have no empirical distribution to test. |
#' | test fit with chi-square | Exact or binned distances. For binned data it uses the survey's bins; for exact data, cutpoints `mrds` chooses. |
#' | fit covariate models | At least one candidate covariate column (see below). |
#' | estimate perception bias | A double-observer structure (see below). |
#' | estimate availability | Never (see below). |
#' | carry group size to abundance | A `size` column. |
#'
#' @section Why structure matters:
#' Whether an analysis is possible depends on the survey design, and the data
#' do not announce it. Exact or binned distances decide which goodness-of-fit
#' test can be computed. Two independent observer teams are needed to
#' estimate *perception bias* (animals at the surface but missed by the
#' observers). That depends on the survey programme, so a pooled extract may
#' or may not carry it.
#'
#' The risk is silent: a single-observer dataset fits without error, and a
#' result reported as if perception had been handled is wrong by a factor.
#' This function states that limit explicitly.
#'
#' @section Candidate covariates:
#' A column is a candidate covariate if all of the following hold:
#' * It is not a structural column: `object`, `distance`, `distbegin`,
#'   `distend`, `observer`, `detected`, `size`, `Region.Label`, `Area`,
#'   `Sample.Label`, `Effort`, `detected_by`, `sample`, `distance_source`, or
#'   a per-sighting NARWC column (identifiers, positions, dates, raw angle and
#'   strip fields). Names are matched without regard to case, after
#'   resolving NARWC aliases.
#' * It is not a superseded `<NAME>_ORIGINAL` column whose stem is a known
#'   variable. `narwcr::read_narwc(prefer_source = TRUE)` keeps a replaced
#'   column under that name for traceability.
#' * It is numeric, character, factor or logical, and has more than one
#'   distinct non-missing value.
#'
#' Candidates that NARWC records once per survey leg (sea state, visibility,
#' glare, altitude, platform) are listed as survey conditions. The NARWC
#' vocabulary comes from the `narwcr` package when it is installed. Without
#' it, only the structural names above are used.
#'
#' @section Duplicate columns:
#' Two columns that name one variable, after case-folding and alias
#' resolution, are reported in `summary$duplicates`. An example is a
#' MEMDR-era NARWC file with both `LEGTYPE` and `LEGTYPE_BK`. The pair is
#' ambiguous: nothing here can tell which one a fit should use.
#' `narwcr::read_narwc(prefer_source = TRUE)` resolves it.
#'
#' @section Perception bias:
#' Mark-recapture distance sampling (`mrds` methods `"io"` and `"trial"`)
#' needs one row per observer per object, with an `observer` column and a
#' `detected` column (the detection indicator). Column values are not
#' checked. The check gives:
#' * `FALSE` if neither `observer` nor `detected` is present, or `observer`
#'   has fewer than two values.
#' * `NA` if only one of the two columns is present, or if `object` is present
#'   but no value repeats (no animal recorded by both teams).
#' * `TRUE` otherwise. Without an `object` column, repeats are not checked.
#'
#' @section What it does not tell you:
#' "Supported" is a statement about structure only. Enough detections to fit a
#' detection function is not enough to fit one well. A covariate being present
#' is not a reason to put it in a model.
#'
#' Availability is unsupported for every table. *Availability* is the
#' probability that an animal is at the surface and can be seen. An animal
#' submerged for the whole pass is missed at every distance equally, so
#' distances carry no information about it. It is computed from dive data
#' with [availability()].
#'
#' @param data Data frame, one row per detection, in the form
#'   [sweep_models()] takes: `distance` for exact distances, or `distbegin`
#'   and `distend` for binned ones.
#'
#' @return An object of class `dsfit_structure`, a list with:
#'   \describe{
#'     \item{`table`}{A tibble with one row per check and columns `check`
#'       (the question), `supported` (`TRUE`, `FALSE`, or `NA` where the
#'       structure is ambiguous or incomplete) and `detail` (the reason).}
#'     \item{`summary`}{A list: `n` (rows), `n_distances` (non-missing
#'       distances), `type` (`"exact"`, `"binned"`, `"both"` or `"none"`),
#'       `covariates` (candidate covariate names), `duplicates` (named list of
#'       column-name pairs naming one variable) and `nearest` (smallest
#'       distance, or smallest `distbegin` for binned data; `NA` if none).}
#'   }
#'
#' @references
#' Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L.
#' and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford University
#' Press. Source of the 40 and 60-80 detection guidelines.
#'
#' @seealso [prepare_distance_data()], which enforces what this only reports.
#'   [availability()] and [g0()] for the components this cannot supply.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(
#'   object = 1:200,
#'   distance = abs(rnorm(200, 0, 120)),
#'   beaufort = sample(0:4, 200, replace = TRUE)
#' )
#' detection_structure(d)
#'
#' # Binned distances admit a different goodness-of-fit test
#' b <- data.frame(object = 1:60, distbegin = rep(c(0, 100, 200), 20),
#'                 distend = rep(c(100, 200, 300), 20))
#' detection_structure(b)
#'
#' @export
detection_structure <- function(data) {
  # 1. Validate the input.
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame with one row per detection.")
  }
  nm <- names(data)
  n <- nrow(data)

  # 2. Classify the distances. A column counts only if it has a non-missing
  # value.
  has_exact <- "distance" %in% nm && any(!is.na(data$distance))
  has_binned <- all(c("distbegin", "distend") %in% nm) &&
    any(!is.na(data$distbegin))

  type <- if (has_exact && has_binned) "both" else if (has_exact) "exact" else
    if (has_binned) "binned" else "none"

  # Rows of the result table are collected here by add().
  rows <- list()
  add <- function(check, supported, detail) {
    rows[[length(rows) + 1L]] <<- tibble::tibble(
      check = check, supported = supported, detail = detail
    )
  }

  # 3. Can a detection function be fitted at all? -------------------------
  # Count the distances; for "both", the exact ones are counted.
  n_dist <- if (has_exact) sum(!is.na(data$distance)) else
    if (has_binned) sum(!is.na(data$distbegin)) else 0L

  if (type == "none") {
    add("fit a detection function", FALSE,
        if (!n) "no rows" else
          if (any(c("distance", "distbegin") %in% nm))
            "the distance columns are present but every value is missing"
          else "no `distance`, and no `distbegin`/`distend` pair")
  } else if (type == "both") {
    add("fit a detection function", NA,
        paste0("both exact and binned distances present; they have different ",
               "likelihoods and cannot be swept together"))
  } else {
    # Compare with Buckland et al. (2001): 60-80 detections suggested, 40
    # workable.
    detail <- paste0(n_dist, " ", type, " distances")
    detail <- paste0(detail, if (n_dist < 40) {
      "; below the 40 that is usually workable"
    } else if (n_dist < 60) {
      "; workable, but below the 60-80 usually suggested"
    } else {
      "; at or above the 60-80 usually suggested"
    })
    add("fit a detection function", TRUE, detail)
  }

  # 4. Goodness-of-fit tests available ---------------------------------------
  add("test fit with Cramer-von Mises", isTRUE(type == "exact"),
      if (identical(type, "exact")) "exact distances have an empirical distribution to test"
      else "needs exact distances; binned fits have no empirical distribution")

  add("test fit with chi-square", type %in% c("exact", "binned"),
      if (identical(type, "binned")) "over the survey's own bins"
      else if (identical(type, "exact"))
        "over cutpoints mrds chooses, which makes it the weaker test here"
      else "needs distances")

  # 5. Candidate covariates --------------------------------------------------
  # 5a. Structural names are never covariates: mrds and flatfile columns,
  # plus NARWC's per-sighting columns when narwcr is installed.
  vocab <- narwc_vocabulary()
  structural <- c("object", "distance", "distbegin", "distend", "observer",
                  "detected", "size", "Region.Label", "Area", "Sample.Label",
                  "Effort", "detected_by", "sample", "distance_source",
                  vocab$structural)

  # 5b. Map each column to its canonical NARWC name, if it is an alias.
  canonical <- vapply(nm, function(x) {
    hit <- vocab$aliases[[toupper(x)]]
    if (is.null(hit)) x else hit
  }, character(1))
  # Compare names in upper case. NARWC names are upper case and a flatfile's
  # often are not: `latitude` is LATITUDE, and BEAUFORT beside beaufort is one
  # variable named twice.
  canon_up <- toupper(canonical)

  # 5c. Exclude superseded `<NAME>_ORIGINAL` columns. read_narwc(prefer_source
  # = TRUE) keeps a replaced column under that name for traceability (e.g.
  # LEGTYPE_ORIGINAL when a MEMDR-era file had both LEGTYPE and LEGTYPE_BK).
  # Only names whose stem is a known variable are excluded, so an unrelated
  # column ending in _ORIGINAL is kept.
  known <- unique(c(structural, vocab$conditions, unlist(vocab$aliases)))
  superseded <- grepl("_ORIGINAL$", toupper(nm)) &
    sub("_ORIGINAL$", "", toupper(nm)) %in% toupper(known)

  # 5d. Keep non-structural, non-superseded columns of a usable type that
  # vary (more than one distinct non-missing value).
  keep <- !(canon_up %in% toupper(structural)) & !superseded
  candidates <- nm[keep]
  candidates <- candidates[vapply(data[candidates], function(x) {
    (is.numeric(x) || is.character(x) || is.factor(x) || is.logical(x)) &&
      length(unique(x[!is.na(x)])) > 1L
  }, logical(1))]

  # 5e. Flag the candidates NARWC records once per leg (survey conditions).
  conditions <- candidates[canon_up[match(candidates, nm)] %in%
                             toupper(vocab$conditions)]

  # 5f. Describe the candidates (first eight named) and add the row.
  detail <- if (!length(candidates)) {
    "no columns beyond the structural ones vary"
  } else {
    d <- paste0("candidates: ", paste(utils::head(candidates, 8), collapse = ", "),
                if (length(candidates) > 8) ", ..." else "")
    if (length(conditions)) {
      d <- paste0(d, "; recorded as survey conditions: ",
                  paste(conditions, collapse = ", "))
    }
    d
  }
  add("fit covariate models", length(candidates) > 0L, detail)

  # 6. Duplicate columns: two or more columns with the same canonical name,
  # across all columns (structural ones included). Named by the canonical
  # name. Unresolved, the pair is ambiguous: nothing here can tell which one
  # a fit should use.
  dup_up <- unique(canon_up[duplicated(canon_up)])
  dups <- lapply(dup_up, function(k) nm[canon_up == k])
  names(dups) <- canonical[match(dup_up, canon_up)]

  # 7. Perception bias: needs a double-observer structure -------------------
  obs <- perception_structure(data)
  add("estimate perception bias", obs$supported, obs$detail)

  # 8. Availability: never estimable from distances -------------------------
  add("estimate availability", FALSE,
      paste0("not estimable from distances by construction: a submerged animal ",
             "is missed at every distance equally. Compute it from dive data ",
             "with availability()"))

  # 9. Group size, needed to convert group abundance to individuals ---------
  add("carry group size to abundance", "size" %in% nm,
      if ("size" %in% nm) "`size` present"
      else "no `size` column; abundance would count groups, not individuals")

  # 10. Assemble the result.
  tab <- do.call(rbind, rows)

  out <- list(
    table = tab,
    summary = list(n = n, n_distances = n_dist, type = type,
                   covariates = candidates, duplicates = dups,
                   nearest = nearest_distance(data, has_exact, has_binned))
  )
  class(out) <- "dsfit_structure"
  out
}


# Get NARWC's column vocabulary from narwcr, if it is installed.
# No inputs. Returns a list: `structural` (per-sighting columns), `conditions`
# (per-leg columns) and `aliases` (alternative name -> canonical name). All
# empty when narwcr is missing or its functions fail.
#
# The split is narwcr's. `narwc_never_fill()` lists columns recorded per
# sighting (identifiers, positions, dates, raw angle and strip fields); none
# is a detection covariate. `narwc_fill_columns()` lists columns recorded once
# per leg, which are survey conditions (sea state, visibility, glare,
# altitude, platform). Without narwcr, classification uses the mrds and
# flatfile names alone. That is correct for non-NARWC data and less informed
# for NARWC data.
narwc_vocabulary <- function() {
  # 1. Without narwcr, return the empty vocabulary.
  empty <- list(structural = character(0), conditions = character(0),
                aliases = list())
  if (!requireNamespace("narwcr", quietly = TRUE)) return(empty)

  # 2. Read the two column lists, treating any error as an empty list.
  grab <- function(f) {
    out <- try(f(), silent = TRUE)
    if (inherits(out, "try-error")) character(0) else as.character(out)
  }
  structural <- grab(narwcr::narwc_never_fill)
  conditions <- grab(narwcr::narwc_fill_columns)
  if (!length(structural) && !length(conditions)) return(empty)

  # 3. Read the alias map, if the schema provides one.
  aliases <- try(narwcr::narwc_schema()$aliases, silent = TRUE)
  aliases <- if (inherits(aliases, "try-error") || is.null(aliases)) list() else
    as.list(aliases)

  list(structural = structural, conditions = conditions, aliases = aliases)
}


# Check for the double-observer structure mrds needs for method = "io" or
# "trial": one row per observer per object, with a detection indicator.
# Input: the detection table. Returns list(supported = TRUE/FALSE/NA,
# detail = reason string).
perception_structure <- function(data) {
  nm <- names(data)
  has_observer <- "observer" %in% nm
  has_detected <- "detected" %in% nm
  has_object <- "object" %in% nm

  # 1. Neither column: a single-observer survey. No sample size fixes this.
  if (!has_observer && !has_detected) {
    return(list(supported = FALSE, detail = paste0(
      "no double-observer structure: no `observer` or `detected` column. ",
      "Perception needs two independent teams and cannot be recovered from a ",
      "single-observer survey at any sample size"
    )))
  }
  # 2. Only one of the two columns: incomplete structure.
  if (!has_observer || !has_detected) {
    missing <- if (has_observer) "`detected`" else "`observer`"
    return(list(supported = NA, detail = paste0(
      "partial double-observer structure: ", missing, " is missing. ",
      "Mark-recapture distance sampling needs one row per observer per object, ",
      "with a detection indicator"
    )))
  }

  # 3. Fewer than two observer values: only one team.
  n_obs <- length(unique(data$observer[!is.na(data$observer)]))
  if (n_obs < 2L) {
    return(list(supported = FALSE, detail = paste0(
      "`observer` takes only ", n_obs, " value; two independent teams are needed"
    )))
  }
  # 4. No repeated object: no animal is recorded by both teams. Skipped when
  # there is no `object` column.
  if (has_object && !anyDuplicated(data$object)) {
    return(list(supported = NA, detail = paste0(
      "`observer` and `detected` present, but no `object` appears twice, so ",
      "no animal was seen by both teams as recorded"
    )))
  }

  # 5. Otherwise supported.
  list(supported = TRUE, detail = paste0(
    n_obs, " observers with a detection indicator; mark-recapture distance ",
    "sampling is possible"
  ))
}


# Find the nearest detection, where a blind spot beneath the platform shows.
# Inputs: the table and the two distance-type flags from
# detection_structure(). Returns the smallest `distance` (exact, or both),
# else the smallest `distbegin` (binned), else NA.
nearest_distance <- function(data, has_exact, has_binned) {
  if (has_exact) return(suppressWarnings(min(data$distance, na.rm = TRUE)))
  if (has_binned) return(suppressWarnings(min(data$distbegin, na.rm = TRUE)))
  NA_real_
}


# Print method for detection_structure() results. Prints the counts, the
# nearest detection (if above zero), any duplicate columns, and the checks
# grouped as "can", "partly" and "cannot". Returns `x` invisibly.
#' @export
print.dsfit_structure <- function(x, ...) {
  # 1. Header: rows, distance count and type.
  s <- x$summary
  cat("<dsfit_structure>\n")
  cat("  ", s$n, " rows, ", s$n_distances, " ", s$type, " distances\n",
      sep = "")
  # 2. Nearest detection, where a blind spot would show.
  if (!is.na(s$nearest) && is.finite(s$nearest) && s$nearest > 0) {
    cat("  nearest detection at ", signif(s$nearest, 4),
        " - a blind spot beneath the platform would show here\n", sep = "")
  }
  # 3. Duplicate columns.
  for (k in names(s$duplicates)) {
    cat("  two columns name ", k, ": ",
        paste(s$duplicates[[k]], collapse = " and "),
        " - narwcr::read_narwc(prefer_source = TRUE) resolves this\n", sep = "")
  }

  # 4. The checks, grouped by answer, with check names padded to align.
  tab <- x$table
  show_block <- function(label, rows) {
    if (!nrow(rows)) return(invisible(NULL))
    cat("\n  ", label, ":\n", sep = "")
    width <- max(nchar(rows$check))
    for (i in seq_len(nrow(rows))) {
      cat("    ", formatC(rows$check[i], width = -width), "  ",
          rows$detail[i], "\n", sep = "")
    }
  }

  show_block("can", tab[!is.na(tab$supported) & tab$supported, , drop = FALSE])
  show_block("partly", tab[is.na(tab$supported), , drop = FALSE])
  show_block("cannot", tab[!is.na(tab$supported) & !tab$supported, , drop = FALSE])

  cat("\n  Structure only. That something is supported is not a reason to do it.\n")
  invisible(x)
}
