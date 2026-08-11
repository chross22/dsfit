#' What a set of detections can and cannot support
#'
#' Reads a table of detections and reports which analyses it admits, and which
#' it does not, with the reason. Nothing is fitted and nothing is estimated.
#'
#' @section Why this exists:
#' Most of what decides whether an analysis is possible is structural, and none
#' of it is announced by the data. Whether distances are exact or binned decides
#' which goodness-of-fit test can be computed. Whether a survey ran two
#' independent observer teams decides whether perception bias is estimable at
#' all — and that is a property of the survey programme, not of the archive its
#' data ends up in, so a pooled extract may or may not carry it.
#'
#' The failure this guards against is not an error but a silence: fitting a
#' single-observer dataset and reporting the result as though perception had
#' been handled. That produces a number, and the number is wrong by a factor.
#' Asking here turns "you had to know that" into something the package says.
#'
#' @section What it does not tell you:
#' That an analysis is *supported* is a statement about structure, not about
#' whether it is a good idea. Enough detections to fit a detection function is
#' not enough detections to fit one well, and a covariate being present is not a
#' reason to put it in a model.
#'
#' Availability is reported as unsupported for every table, which is not a
#' defect in any particular dataset: it cannot be estimated from sighting
#' distances by construction, because an animal submerged for the whole pass is
#' missed at every distance equally and leaves no signature. It is computed from
#' dive data instead — see [availability()].
#'
#' @param data A data frame with one row per detection, as [sweep_models()]
#'   takes: `distance` for exact distances, or `distbegin` and `distend` for
#'   binned ones.
#'
#' @return An object of class `dsfit_structure`: a list with `table` (one row
#'   per question, with `check`, `supported` and `detail`) and `summary` (the
#'   counts and distance type behind it). `supported` is `TRUE`, `FALSE`, or
#'   `NA` where the structure is present but incomplete.
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
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame with one row per detection.")
  }
  nm <- names(data)
  n <- nrow(data)

  has_exact <- "distance" %in% nm && any(!is.na(data$distance))
  has_binned <- all(c("distbegin", "distend") %in% nm) &&
    any(!is.na(data$distbegin))

  type <- if (has_exact && has_binned) "both" else if (has_exact) "exact" else
    if (has_binned) "binned" else "none"

  rows <- list()
  add <- function(check, supported, detail) {
    rows[[length(rows) + 1L]] <<- tibble::tibble(
      check = check, supported = supported, detail = detail
    )
  }

  # --- Fitting a detection function at all -------------------------------
  n_dist <- if (has_exact) sum(!is.na(data$distance)) else
    if (has_binned) sum(!is.na(data$distbegin)) else 0L

  if (type == "none") {
    add("fit a detection function", FALSE,
        "no `distance`, and no `distbegin`/`distend` pair")
  } else if (type == "both") {
    add("fit a detection function", NA,
        paste0("both exact and binned distances present; they have different ",
               "likelihoods and cannot be swept together"))
  } else {
    # Buckland et al. (2001) suggest 60-80 detections, with 40 workable.
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

  # --- Goodness of fit ----------------------------------------------------
  add("test fit with Cramer-von Mises", isTRUE(type == "exact"),
      if (identical(type, "exact")) "exact distances have an empirical distribution to test"
      else "needs exact distances; binned fits have no empirical distribution")

  add("test fit with chi-square", type %in% c("exact", "binned"),
      if (identical(type, "binned")) "over the survey's own bins"
      else if (identical(type, "exact"))
        "over cutpoints mrds chooses, which makes it the weaker test here"
      else "needs distances")

  # --- Covariates ---------------------------------------------------------
  vocab <- narwc_vocabulary()
  structural <- c("object", "distance", "distbegin", "distend", "observer",
                  "detected", "size", "Region.Label", "Area", "Sample.Label",
                  "Effort", "detected_by", "sample", vocab$structural)

  canonical <- vapply(nm, function(x) {
    hit <- vocab$aliases[[toupper(x)]]
    if (is.null(hit)) x else hit
  }, character(1))

  keep <- !(canonical %in% structural)
  candidates <- nm[keep]
  candidates <- candidates[vapply(data[candidates], function(x) {
    (is.numeric(x) || is.character(x) || is.factor(x) || is.logical(x)) &&
      length(unique(x[!is.na(x)])) > 1L
  }, logical(1))]

  conditions <- candidates[canonical[match(candidates, nm)] %in% vocab$conditions]

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

  # --- Perception: the one that is usually the answer ---------------------
  obs <- perception_structure(data)
  add("estimate perception bias", obs$supported, obs$detail)

  # --- Availability: never, and not a defect ------------------------------
  add("estimate availability", FALSE,
      paste0("not estimable from distances by construction: a submerged animal ",
             "is missed at every distance equally. Compute it from dive data ",
             "with availability()"))

  # --- Group size, which abundance needs downstream ------------------------
  add("carry group size to abundance", "size" %in% nm,
      if ("size" %in% nm) "`size` present"
      else "no `size` column; abundance would count groups, not individuals")

  tab <- do.call(rbind, rows)

  out <- list(
    table = tab,
    summary = list(n = n, n_distances = n_dist, type = type,
                   covariates = candidates,
                   nearest = nearest_distance(data, has_exact, has_binned))
  )
  class(out) <- "dsfit_structure"
  out
}


# NARWC's own column vocabulary, when narwcr is installed to supply it.
#
# The distinction that matters here is narwcr's, not ours. `narwc_never_fill()`
# is the columns recorded per sighting - identifiers, positions, dates, the raw
# angle and strip fields - none of which is a detection covariate.
# `narwc_fill_columns()` is the columns recorded once per leg, which is what a
# survey condition is: sea state, visibility, glare, altitude, platform.
#
# Without narwcr this returns empty and the classification falls back to the
# `mrds` and flatfile names alone, which is right for non-NARWC data and merely
# less informed for NARWC data.
narwc_vocabulary <- function() {
  empty <- list(structural = character(0), conditions = character(0),
                aliases = list())
  if (!requireNamespace("narwcr", quietly = TRUE)) return(empty)

  grab <- function(f) {
    out <- try(f(), silent = TRUE)
    if (inherits(out, "try-error")) character(0) else as.character(out)
  }
  structural <- grab(narwcr::narwc_never_fill)
  conditions <- grab(narwcr::narwc_fill_columns)
  if (!length(structural) && !length(conditions)) return(empty)

  aliases <- try(narwcr::narwc_schema()$aliases, silent = TRUE)
  aliases <- if (inherits(aliases, "try-error") || is.null(aliases)) list() else
    as.list(aliases)

  list(structural = structural, conditions = conditions, aliases = aliases)
}


# Double-observer structure, as mrds would need it for method = "io"/"trial":
# one row per observer per object, with a detection indicator.
perception_structure <- function(data) {
  nm <- names(data)
  has_observer <- "observer" %in% nm
  has_detected <- "detected" %in% nm
  has_object <- "object" %in% nm

  if (!has_observer && !has_detected) {
    return(list(supported = FALSE, detail = paste0(
      "no double-observer structure: no `observer` or `detected` column. ",
      "Perception needs two independent teams and cannot be recovered from a ",
      "single-observer survey at any sample size"
    )))
  }
  if (!has_observer || !has_detected) {
    missing <- if (has_observer) "`detected`" else "`observer`"
    return(list(supported = NA, detail = paste0(
      "partial double-observer structure: ", missing, " is missing. ",
      "Mark-recapture distance sampling needs one row per observer per object, ",
      "with a detection indicator"
    )))
  }

  n_obs <- length(unique(data$observer[!is.na(data$observer)]))
  if (n_obs < 2L) {
    return(list(supported = FALSE, detail = paste0(
      "`observer` takes only ", n_obs, " value; two independent teams are needed"
    )))
  }
  if (has_object && !anyDuplicated(data$object)) {
    return(list(supported = NA, detail = paste0(
      "`observer` and `detected` present, but no `object` appears twice, so ",
      "no animal was seen by both teams as recorded"
    )))
  }

  list(supported = TRUE, detail = paste0(
    n_obs, " observers with a detection indicator; mark-recapture distance ",
    "sampling is possible"
  ))
}


# The closest detection, which is where a blind spot would show.
nearest_distance <- function(data, has_exact, has_binned) {
  if (has_exact) return(suppressWarnings(min(data$distance, na.rm = TRUE)))
  if (has_binned) return(suppressWarnings(min(data$distbegin, na.rm = TRUE)))
  NA_real_
}


#' @export
print.dsfit_structure <- function(x, ...) {
  s <- x$summary
  cat("<dsfit_structure>\n")
  cat("  ", s$n, " rows, ", s$n_distances, " ", s$type, " distances\n",
      sep = "")
  if (!is.na(s$nearest) && is.finite(s$nearest) && s$nearest > 0) {
    cat("  nearest detection at ", signif(s$nearest, 4),
        " - a blind spot beneath the platform would show here\n", sep = "")
  }

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
