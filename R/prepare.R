#' Check and prepare distance data for fitting
#'
#' Checks a table of detections and returns it in the form `mrds::ddf()`
#' expects. [sweep_models()] calls it before fitting. It is exported so the
#' checks can be run on their own before a sweep.
#'
#' Distances are either *exact* (one `distance` per detection) or *binned*
#' (an interval `distbegin` to `distend` per detection, as recorded by strip
#' or band surveys). The two are fitted with different likelihoods.
#'
#' @section How it works:
#' The steps run in this order. The first refusal stops the function with an
#' error.
#'
#' 1. Check the arguments. `data` must be a data frame. `truncation` must be
#'    one positive finite number. `left`, if given, must be one non-negative
#'    finite number smaller than `truncation`.
#' 2. Check the distance columns. Any of `distance`, `distbegin` and `distend`
#'    that is present must be numeric (an all-`NA` column is allowed) and must
#'    not be negative.
#' 3. Decide the distance type. Refuse a table with no distances, and refuse a
#'    table with both exact and binned distances.
#' 4. Check `distance_source`, if present, on the rows with an exact
#'    `distance`. Refuse `"break_off"` distances. Warn about `"circling"`
#'    distances. Only rows that will be fitted (within `left` to `truncation`)
#'    are counted.
#' 5. Count, by reason, the rows that will be dropped: no distance, beyond
#'    `truncation`, or inside `left`.
#' 6. For binned data:
#'    1. Drop rows with a missing `distbegin` or `distend`.
#'    2. Refuse any `distend` of `Inf`.
#'    3. Drop bins that end beyond `truncation` or begin inside `left`.
#'    4. Refuse if no rows remain.
#'    5. Derive `breaks` from the bins, refusing bins that do not tile; or,
#'       if `breaks` is supplied, check that it contains every recorded bin
#'       edge.
#'    6. Refuse if the highest break is not `truncation`.
#'    7. Refuse if the lowest break is not `left` (or 0 when `left` is
#'       `NULL`).
#'    8. Remove the `distance` column.
#' 7. For exact data:
#'    1. Warn that `breaks` is ignored, if it was supplied.
#'    2. Remove the `distbegin` and `distend` columns.
#'    3. Drop rows with no distance, beyond `truncation`, or inside `left`.
#'    4. Refuse if no rows remain.
#' 8. Add an `object` column numbering the rows if there is none. Refuse
#'    duplicated `object` values.
#'
#' Rows with no distance are how a flatfile records a sample with no
#' detections, so dropping them is expected. Every dropped row is counted in
#' the return value.
#'
#' @section What it refuses and why:
#' | Condition | Reason |
#' |---|---|
#' | Non-numeric distances | Text or factor codes would be compared to `truncation` as text or level numbers, keeping the wrong rows. |
#' | Negative distances | A perpendicular distance has no side. Use `abs()` if the sign records the side of the line. |
#' | Exact and binned distances together | The two likelihoods differ, so their AICs cannot be compared in one sweep. Split the table on its distance-source column and sweep each part. |
#' | `distance_source == "break_off"` | A break-off distance runs from where the aircraft left the line to the logged position. It is not perpendicular to the line. |
#' | `distend` of `Inf` | An open top bin cannot be fitted. Drop those rows or close the bin. This applies even if truncation would drop the bin. |
#' | Bins that do not tile (derived `breaks`) | If bins leave gaps or overlap, they do not define one set of cutpoints. |
#' | `breaks` missing a recorded bin edge | Goodness of fit is computed over `breaks`, so `breaks` must be the survey's own cutpoints. |
#' | Highest break not equal to `truncation` | A binned fit integrates over the bins. A gap between the top bin and `truncation` is effort with no bin. |
#' | Lowest break not equal to `left` (or 0) | The same gap at the near end. It occurs when the innermost bin had no detections, or when `left` falls inside a bin. Supply `breaks` to fix it. |
#' | No rows left after truncation | There is nothing to fit. |
#' | Duplicated `object` | `mrds` identifies detections by `object`. |
#'
#' @section What it warns about:
#' * `distance_source == "circling"`: the animal's position while the aircraft
#'   circled, projected onto the line. It is a perpendicular distance, but
#'   measured minutes after the detection. `distsamp::detection_data()`
#'   excludes these rows unless `include_circling = TRUE`. They are kept, and
#'   the warning gives their count.
#' * `breaks` supplied with exact distances: it is ignored.
#'
#' @param data Data frame, one row per detection. Needs a numeric `distance`
#'   column (exact) or numeric `distbegin` and `distend` columns (binned).
#'   Optional columns: `object` (detection ID) and `distance_source`.
#' @param truncation Number. Right truncation distance, in the units of the
#'   distances. Required.
#' @param left Number or `NULL` (default). Left truncation distance. `NULL`
#'   means no left truncation.
#' @param breaks Numeric vector or `NULL` (default). Bin cutpoints for binned
#'   data, strictly increasing. Derived from the bins when `NULL`. Ignored,
#'   with a warning, for exact data.
#'
#' @return A list:
#'   \describe{
#'     \item{`data`}{The prepared data frame, ready for `mrds::ddf()`.}
#'     \item{`binned`}{`TRUE` for binned data, `FALSE` for exact.}
#'     \item{`breaks`}{The cutpoints used (binned), or `NULL` (exact).}
#'     \item{`n_dropped`}{Number of input rows not in `data`.}
#'     \item{`dropped`}{Named integer vector splitting the drops by reason:
#'       `no_distance`, `beyond_truncation`, `inside_left`. It shows whether
#'       truncation trimmed a tail or removed most of the survey.}
#'   }
#'
#' @seealso [sweep_models()]
#'
#' @examples
#' d <- data.frame(object = 1:5, distance = c(10, 50, 120, NA, 900))
#' prep <- prepare_distance_data(d, truncation = 400)
#' prep$n_dropped
#' prep$dropped
#' prep$data
#'
#' @export
prepare_distance_data <- function(data, truncation, left = NULL,
                                  breaks = NULL) {
  # 1. Validate the arguments.
  stopifnot(is.data.frame(data))
  stopifnot(is.numeric(truncation), length(truncation) == 1L,
            is.finite(truncation), truncation > 0)
  if (!is.null(left)) {
    stopifnot(is.numeric(left), length(left) == 1L, is.finite(left), left >= 0)
    if (left >= truncation) {
      rlang::abort("`left` must be smaller than `truncation`.")
    }
  }

  # 2. Check the distance columns are numeric and non-negative. A character or
  # factor distance would be compared to the truncation as text or as level
  # codes, keeping and dropping the wrong rows.
  for (col in intersect(c("distance", "distbegin", "distend"), names(data))) {
    x <- data[[col]]
    if (!is.numeric(x) && !all(is.na(x))) {
      rlang::abort(paste0(
        "`", col, "` must be numeric, not ", class(x)[1], ". Convert it with ",
        "`as.numeric()` - for a factor, `as.numeric(as.character(x))`."
      ))
    }
    if (is.numeric(x) && any(x < 0, na.rm = TRUE)) {
      rlang::abort(paste0(
        "`", col, "` has ", sum(x < 0, na.rm = TRUE), " negative value(s). ",
        "A perpendicular distance has no side; take `abs()` if a negative ",
        "sign records which side of the line the animal was on."
      ))
    }
  }

  # 3. Decide whether the distances are exact or binned. A column counts only
  # if it has at least one non-missing value.
  has_point <- "distance" %in% names(data) && any(!is.na(data$distance))
  has_bin <- all(c("distbegin", "distend") %in% names(data)) &&
    any(!is.na(data$distbegin))

  # Neither: say whether the table is empty, the columns are all missing, or
  # the columns are absent.
  if (!has_point && !has_bin) {
    rlang::abort(if (!nrow(data)) "`data` has no rows." else
      if (any(c("distance", "distbegin") %in% names(data))) paste0(
        "No distances found: the distance columns are present but every ",
        "value is missing."
      ) else paste0(
        "No distances found. `data` needs a `distance` column, or `distbegin` ",
        "and `distend` for binned data."
      ))
  }
  # Both: the likelihoods differ, so the AICs cannot be compared.
  if (has_point && has_bin) {
    rlang::abort(paste0(
      "Both point and interval distances are present, and one sweep cannot ",
      "rank both: binned and exact fits have different likelihoods, so their ",
      "AICs are not comparable. Split the table on its distance-source column ",
      "and sweep each separately."
    ))
  }

  # 4. Refuse break-off distances and warn about circling ones.
  check_circling(data, truncation, left)

  n0 <- nrow(data)

  # 5. Count the rows to be dropped, by reason. The three reasons are
  # exclusive, so each dropped row is counted once. A total alone would hide
  # whether truncation trimmed a tail or removed half the survey.
  if (has_bin) {
    gone_missing <- is.na(data$distbegin) | is.na(data$distend)
    gone_beyond <- !gone_missing & data$distend > truncation
    gone_inside <- if (is.null(left)) rep(FALSE, nrow(data)) else
      !gone_missing & !gone_beyond & data$distbegin < left
  } else {
    gone_missing <- is.na(data$distance)
    gone_beyond <- !gone_missing & data$distance > truncation
    gone_inside <- if (is.null(left)) rep(FALSE, nrow(data)) else
      !gone_missing & !gone_beyond & data$distance < left
  }
  dropped <- c(
    no_distance = sum(gone_missing),
    beyond_truncation = sum(gone_beyond),
    inside_left = sum(gone_inside)
  )

  # 6. Binned data.
  if (has_bin) {
    # 6.1 Drop rows with no bin.
    keep <- !is.na(data$distbegin) & !is.na(data$distend)
    data <- data[keep, , drop = FALSE]

    # 6.2 Refuse an open top bin. Checked before truncation, so an open bin is
    # refused even when truncation would have dropped it.
    if (any(is.infinite(data$distend))) {
      rlang::abort(paste0(
        "An unbounded top bin (`distend` of `Inf`) cannot be fitted. Drop ",
        "those detections or close the bin before sweeping."
      ))
    }

    # 6.3 Truncate whole bins: keep a bin only if it lies inside
    # [left, truncation].
    data <- data[data$distend <= truncation, , drop = FALSE]
    if (!is.null(left)) {
      data <- data[data$distbegin >= left, , drop = FALSE]
    }
    # 6.4 Refuse if nothing is left.
    if (!nrow(data)) {
      rlang::abort("No detections left after truncation.")
    }

    # 6.5 Derive the cutpoints from the bins, or check the supplied ones.
    if (is.null(breaks)) {
      breaks <- derive_breaks(data)
    } else {
      check_breaks(breaks, data)
    }

    # 6.6 The top break must equal the truncation. A binned fit integrates the
    # detection function over the bins, so a gap between the top bin and the
    # truncation is effort with no bin, and the fit describes a narrower
    # survey than the one flown.
    top <- max(breaks)
    if (!isTRUE(all.equal(top, truncation))) {
      rlang::abort(paste0(
        "The top bin ends at ", top, " but `truncation` is ", truncation,
        ". A binned fit integrates over the bins, so a gap between the last ",
        "one and the truncation width is unaccounted effort. Either set ",
        "`truncation = ", top, "`, or supply `breaks` reaching ", truncation, "."
      ))
    }

    # 6.7 The bottom break must equal `left`, or 0 without it, because mrds
    # integrates from there. A gap arises from derived breaks when the
    # innermost bin had no detections, or from a `left` inside a bin (whose
    # detections were dropped in 6.3).
    bottom <- min(breaks)
    start <- if (is.null(left)) 0 else left
    if (!isTRUE(all.equal(bottom, start))) {
      rlang::abort(paste0(
        "The lowest bin begins at ", bottom, " but the fit starts at ", start,
        if (is.null(left)) " (no `left` is set)" else " (`left`)",
        ". A binned fit integrates over the bins, so the strip between them is ",
        "unaccounted effort. Either set `left = ", bottom, "`, or supply ",
        "`breaks` starting at ", start, " - an inner bin with no detections ",
        "is still a bin, and derived breaks cannot know it was there."
      ))
    }

    # 6.8 Remove `distance`. A flatfile carries an all-NA `distance` on binned
    # rows, and mrds would read it as the distances and report "no data to
    # analyze".
    data$distance <- NULL
    return(list(data = ensure_object(data), binned = TRUE, breaks = breaks,
                n_dropped = n0 - nrow(data), dropped = dropped))
  }

  # 7. Exact data.
  # 7.1 Cutpoints are only used by a binned fit.
  if (!is.null(breaks)) {
    rlang::warn(paste0(
      "`breaks` is ignored: the distances are exact, and only a binned fit ",
      "uses cutpoints."
    ))
  }
  # 7.2 Remove empty bin columns. Left in, they make mrds warn once per
  # candidate model that the data has bins but is not fitted as binned.
  data$distbegin <- NULL
  data$distend <- NULL

  # 7.3 Drop rows with no distance or outside [left, truncation].
  keep <- !is.na(data$distance) & data$distance <= truncation
  if (!is.null(left)) keep <- keep & data$distance >= left
  data <- data[keep, , drop = FALSE]
  # 7.4 Refuse if nothing is left.
  if (!nrow(data)) {
    rlang::abort("No detections left after truncation.")
  }

  # 8. Number or check `object` (inside ensure_object()) and return.
  list(data = ensure_object(data), binned = FALSE, breaks = NULL,
       n_dropped = n0 - nrow(data), dropped = dropped)
}

# Check the `distance_source` column for distances measured while circling.
# Input: the detection table. Returns NULL invisibly; errors or warns as a
# side effect.
#
# `distsamp` labels two kinds of distance taken after the aircraft broke off
# the line to circle a sighting. `distsamp::detection_data()` excludes both
# unless `include_circling = TRUE`.
# - "break_off": the great-circle distance from where the aircraft left the
#   line to where the record was logged. It is not a perpendicular distance,
#   so a detection function is not defined on it. Refused.
# - "circling": an exact position projected onto the line. It is a
#   perpendicular distance, but measured minutes after the detection, after
#   the animal and aircraft have moved. Kept, with a warning giving the count.
#
# Only rows that would be fitted are counted: a distance within `left` to
# `truncation`. Tables without a `distance` column (binned only) are not
# checked.
check_circling <- function(data, truncation, left = NULL) {
  # 1. Nothing to check without both columns.
  if (!"distance_source" %in% names(data) || !"distance" %in% names(data)) {
    return(invisible(NULL))
  }
  src <- as.character(data$distance_source)
  d <- data$distance
  lower <- if (is.null(left)) 0 else left
  fitted <- !is.na(d) & d >= lower & d <= truncation

  # 2. Refuse break-off distances.
  n_break <- sum(fitted & src %in% "break_off")
  if (n_break > 0) {
    rlang::abort(paste0(
      n_break, " distance(s) have `distance_source == \"break_off\"`: slant ",
      "distances from where the aircraft broke off to circle, not ",
      "perpendicular distances, so not what a detection function is defined ",
      "on. Drop them before fitting (they still count towards abundance in ",
      "`distsamp`'s segments), e.g. ",
      "`data[!data$distance_source %in% \"break_off\", ]`."
    ))
  }

  # 3. Warn about circling distances, so keeping them is a stated choice.
  n_circ <- sum(fitted & src %in% "circling")
  if (n_circ > 0) {
    rlang::warn(paste0(
      n_circ, " distance(s) have `distance_source == \"circling\"`: ",
      "perpendicular, but from a position fixed minutes after the detection. ",
      "`distsamp::detection_data()` excludes them by default. They are the ",
      "close or conspicuous groups worth breaking off for, so keeping or ",
      "dropping them moves the near end of the distance distribution - state ",
      "which in the methods."
    ))
  }
  invisible(NULL)
}

# Test whether two numbers are equal up to floating-point error.
# Inputs: numeric vectors `a` and `b`. Returns a logical vector.
# Needed because bin edges converted between units differ in the last digits
# (0.1 * 3 against 0.3). The tolerance is relative, with a floor of 1e-8.
near <- function(a, b) abs(a - b) <= 1e-8 * pmax(1, abs(a), abs(b))

# Sort and deduplicate a numeric vector, treating near-equal values as equal.
# Input: numeric vector. Returns the sorted unique values, keeping the first
# (smallest) of each near-equal run.
unique_near <- function(x) {
  x <- sort(unique(x))
  if (length(x) < 2L) return(x)
  x[c(TRUE, !near(x[-1], x[-length(x)]))]
}

# Derive bin cutpoints from the recorded bins.
# Input: binned data with `distbegin` and `distend`. Returns the sorted unique
# bin edges. Errors if the bins do not tile, i.e. if any bin does not begin
# where the previous one ends. A gap or overlap means the bins do not define
# one set of cutpoints, and using one anyway would put detections in the
# wrong bin.
derive_breaks <- function(data) {
  # 1. Collect every edge, merging near-equal values.
  edges <- unique_near(c(data$distbegin, data$distend))

  # 2. List the distinct bins. Each edge is first replaced by its
  # representative in `edges`, so one bin written two ways (0.1 * 3 and 0.3)
  # counts once.
  snap <- function(v) edges[vapply(v, function(e) which(near(e, edges))[1],
                                   integer(1))]
  bins <- unique(data.frame(b = snap(data$distbegin), e = snap(data$distend)))
  bins <- bins[order(bins$b), , drop = FALSE]

  # 3. Check each bin begins where the previous one ends.
  if (nrow(bins) > 1) {
    gaps <- !near(bins$b[-1], bins$e[-nrow(bins)])
    if (any(gaps)) {
      first_gap <- which(gaps)[1]
      rlang::abort(paste0(
        "The bins do not tile: interval ", first_gap,
        " ends at ", bins$e[first_gap], " and the next begins at ",
        bins$b[first_gap + 1], ". Supply `breaks` explicitly."
      ))
    }
  }
  edges
}

# Check user-supplied cutpoints against the recorded bins.
# Inputs: `breaks` and binned data. Returns nothing useful; errors on failure.
# `breaks` must be strictly increasing and contain every recorded bin edge.
# Extra breaks are allowed. Why: the likelihood uses each row's own
# distbegin/distend, but goodness of fit uses `breaks`. Cutpoints that cross
# the recorded bins leave AIC unchanged and compute the chi-square over bins
# that were never surveyed.
check_breaks <- function(breaks, data) {
  # 1. Shape: numeric, no NA, at least two values, strictly increasing.
  if (!is.numeric(breaks) || anyNA(breaks) || length(breaks) < 2L ||
      is.unsorted(breaks, strictly = TRUE)) {
    rlang::abort("`breaks` must be at least two strictly increasing numbers.")
  }
  # 2. Every recorded edge must be (nearly) one of the breaks. Up to five
  # missing edges are named in the error.
  edges <- unique(c(data$distbegin, data$distend))
  off <- edges[!vapply(edges, function(e) any(near(e, breaks)), logical(1))]
  if (length(off)) {
    rlang::abort(paste0(
      "`breaks` does not contain the bin edge(s) ",
      paste(utils::head(sort(off), 5), collapse = ", "),
      if (length(off) > 5) ", ..." else "", " that the distances were ",
      "recorded in. Goodness of fit is computed over `breaks`, so they must ",
      "be the survey's own cutpoints."
    ))
  }
}

# Make sure every row has a unique `object` ID.
# Input: the prepared data. Returns it with `object` added (1, 2, ...) if it
# was absent. Errors on duplicated IDs. mrds identifies detections by
# `object` and also refuses duplicates, but with a less helpful message.
ensure_object <- function(data) {
  if (!"object" %in% names(data)) {
    data$object <- seq_len(nrow(data))
  }
  dup <- unique(data$object[duplicated(data$object)])
  if (length(dup)) {
    rlang::abort(paste0(
      "`object` must identify one detection per row; ", length(dup),
      " value(s) repeat, e.g. ", paste(utils::head(dup, 3), collapse = ", "),
      ". Drop `object` to have it numbered, or deduplicate the rows."
    ))
  }
  data
}
