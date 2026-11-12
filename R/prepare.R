#' Check and prepare distance data for fitting
#'
#' Validates the table a sweep is about to be fitted to and puts it in the shape
#' `mrds` wants. Called by [sweep_models()]; exported because the checks are
#' worth running on their own before committing to a model set.
#'
#' @section What it refuses:
#' \describe{
#'   \item{Both point and interval distances}{`STRIP`-derived distances are
#'     intervals and are fitted binned; angle- and position-derived distances
#'     are points. The likelihoods differ, so one sweep cannot rank both. A
#'     table containing both marks a survey-era boundary, and should be split on
#'     its provenance column and swept separately.}
#'   \item{An unbounded top bin}{`distend` of `Inf` cannot be fitted. The open
#'     top bin of every `STRIP` scheme has to be dropped or closed before
#'     fitting.}
#'   \item{Bins that do not tile}{Interval data whose `distbegin`/`distend`
#'     pairs leave gaps or overlap does not define a set of cutpoints, and any
#'     `breaks` derived from it would silently misallocate detections.}
#'   \item{Bins that stop short of the truncation}{A binned fit integrates the
#'     detection function over the bins. If the top bin ends before the
#'     truncation width, the strip between them is unaccounted effort and the
#'     fit describes a narrower survey than the one flown.}
#'   \item{Bins that start beyond `left`}{The same gap at the near end: the
#'     lowest break must equal `left`, or zero without it. Derived breaks fall
#'     short when the innermost bin had no detections, and a `left` inside a
#'     bin drops that bin's detections; supply `breaks` in either case.}
#'   \item{Break-off distances}{`distance_source == "break_off"`, which
#'     `distsamp` gives a circling sighting with no exact position: a slant
#'     distance from where the aircraft left the line, not a perpendicular
#'     one.}
#' }
#'
#' @section What it warns about:
#' `distance_source == "circling"`: a circling sighting's exact position
#' projected onto the census line. The right quantity, but fixed minutes after
#' the detection, and excluded by `distsamp::detection_data()` unless
#' `include_circling = TRUE`. Admitted, with the count, so the choice is made
#' knowingly.
#'
#' @section What it drops, and reports:
#' Rows with no distance at all — which is how a flatfile records a sample that
#' produced no detections — and rows beyond `truncation` or inside `left`.
#' Dropping is counted, never silent.
#'
#' @param data A data frame with `distance`, or with `distbegin` and `distend`.
#' @param truncation Right truncation distance.
#' @param left Left truncation distance, or `NULL`.
#' @param breaks Bin cutpoints for interval data. Derived from the data when
#'   `NULL`.
#'
#' @return A list with `data` (ready for `mrds`), `binned`, `breaks`, and
#'   `n_dropped`, and `dropped` — the attrition split by reason, so a
#'   truncation that trimmed a tail is distinguishable from one that threw away
#'   half the survey.
#'
#' @seealso [sweep_models()]
#'
#' @examples
#' d <- data.frame(object = 1:5, distance = c(10, 50, 120, NA, 900))
#' prep <- prepare_distance_data(d, truncation = 400)
#' prep$n_dropped
#' prep$data
#'
#' @export
prepare_distance_data <- function(data, truncation, left = NULL,
                                  breaks = NULL) {
  stopifnot(is.data.frame(data))
  stopifnot(is.numeric(truncation), length(truncation) == 1L,
            is.finite(truncation), truncation > 0)
  if (!is.null(left)) {
    stopifnot(is.numeric(left), length(left) == 1L, is.finite(left), left >= 0)
    if (left >= truncation) {
      rlang::abort("`left` must be smaller than `truncation`.")
    }
  }

  # A character or factor distance would be compared as text or as level
  # codes against the truncation, keeping and dropping the wrong rows.
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

  has_point <- "distance" %in% names(data) && any(!is.na(data$distance))
  has_bin <- all(c("distbegin", "distend") %in% names(data)) &&
    any(!is.na(data$distbegin))

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
  if (has_point && has_bin) {
    rlang::abort(paste0(
      "Both point and interval distances are present, and one sweep cannot ",
      "rank both: binned and exact fits have different likelihoods, so their ",
      "AICs are not comparable. Split the table on its distance-source column ",
      "and sweep each separately."
    ))
  }

  check_circling(data)

  n0 <- nrow(data)

  # Attrition by reason, partitioned so every dropped row is counted once.
  # A total on its own hides the difference between a truncation that trimmed
  # a tail and one that threw away half the survey.
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

  if (has_bin) {
    keep <- !is.na(data$distbegin) & !is.na(data$distend)
    data <- data[keep, , drop = FALSE]

    if (any(is.infinite(data$distend))) {
      rlang::abort(paste0(
        "An unbounded top bin (`distend` of `Inf`) cannot be fitted. Drop ",
        "those detections or close the bin before sweeping."
      ))
    }

    data <- data[data$distend <= truncation, , drop = FALSE]
    if (!is.null(left)) {
      data <- data[data$distbegin >= left, , drop = FALSE]
    }
    if (!nrow(data)) {
      rlang::abort("No detections left after truncation.")
    }

    if (is.null(breaks)) {
      breaks <- derive_breaks(data)
    } else {
      check_breaks(breaks, data)
    }

    # A binned fit integrates the detection function over the bins. If the top
    # bin stops short of the truncation width there is a strip of unaccounted
    # effort between them, and the fit silently describes a narrower survey
    # than the one flown.
    top <- max(breaks)
    if (!isTRUE(all.equal(top, truncation))) {
      rlang::abort(paste0(
        "The top bin ends at ", top, " but `truncation` is ", truncation,
        ". A binned fit integrates over the bins, so a gap between the last ",
        "one and the truncation width is unaccounted effort. Either set ",
        "`truncation = ", top, "`, or supply `breaks` reaching ", truncation, "."
      ))
    }

    # The same holds at the near end. mrds integrates from `left` (zero without
    # it), so a lowest bin that begins further out leaves a strip of effort
    # with no bin to put detections in. It arises from derived breaks when the
    # innermost bin had no detections, and from a `left` that falls inside a
    # bin, whose detections are then dropped.
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

    # An all-NA `distance` column, which a flatfile carries for binned rows,
    # is read by mrds as the distances and leaves it "no data to analyze".
    data$distance <- NULL
    return(list(data = ensure_object(data), binned = TRUE, breaks = breaks,
                n_dropped = n0 - nrow(data), dropped = dropped))
  }

  if (!is.null(breaks)) {
    rlang::warn(paste0(
      "`breaks` is ignored: the distances are exact, and only a binned fit ",
      "uses cutpoints."
    ))
  }
  # Likewise empty bin columns on exact data, which make mrds warn once per
  # candidate that the data has bins but is not being fitted binned.
  data$distbegin <- NULL
  data$distend <- NULL

  keep <- !is.na(data$distance) & data$distance <= truncation
  if (!is.null(left)) keep <- keep & data$distance >= left
  data <- data[keep, , drop = FALSE]
  if (!nrow(data)) {
    rlang::abort("No detections left after truncation.")
  }

  list(data = ensure_object(data), binned = FALSE, breaks = NULL,
       n_dropped = n0 - nrow(data), dropped = dropped)
}

# Distances `distsamp` derives from a break-off to circle, which
# `detection_data()` excludes unless `include_circling = TRUE`. They differ in
# kind, and are treated accordingly:
#
# - "break_off" is the great-circle distance from where the aircraft left the
#   line to where the record was logged. It is not perpendicular to anything,
#   so it is not the quantity a detection function is defined on. Refused.
# - "circling" is an exact position projected onto the census line: the right
#   quantity, but fixed minutes after the detection, after the animal and the
#   aircraft have both moved. Admitted, with a warning naming the count.
#
# Only distances that would be fitted are counted; a circling row with no
# distance is dropped with the rest of the missing ones.
check_circling <- function(data) {
  if (!"distance_source" %in% names(data) || !"distance" %in% names(data)) {
    return(invisible(NULL))
  }
  src <- as.character(data$distance_source)
  fitted <- !is.na(data$distance)

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

# Two edges that differ only by floating-point error, as bins converted
# between units do (0.1 * 3 against 0.3).
near <- function(a, b) abs(a - b) <= 1e-8 * pmax(1, abs(a), abs(b))

# Sorted unique values, with near-equal ones collapsed onto the first.
unique_near <- function(x) {
  x <- sort(unique(x))
  if (length(x) < 2L) return(x)
  x[c(TRUE, !near(x[-1], x[-length(x)]))]
}

# Cutpoints implied by the bins present. They have to tile: a gap or an overlap
# means the intervals do not define a single set of breaks, and deriving one
# anyway would put detections in the wrong bin.
derive_breaks <- function(data) {
  edges <- unique_near(c(data$distbegin, data$distend))
  # Each recorded edge is replaced by its representative in `edges`, so the
  # same bin written two ways (0.1 * 3 and 0.3) is one bin.
  snap <- function(v) edges[vapply(v, function(e) which(near(e, edges))[1],
                                   integer(1))]
  bins <- unique(data.frame(b = snap(data$distbegin), e = snap(data$distend)))
  bins <- bins[order(bins$b), , drop = FALSE]

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

# Supplied cutpoints must be the bins the distances were recorded in. The
# likelihood uses each row's own distbegin/distend, but goodness of fit uses
# `breaks`, so cutpoints that cut across the recorded bins leave AIC unchanged
# and the chi-square computed over bins nobody surveyed.
check_breaks <- function(breaks, data) {
  if (!is.numeric(breaks) || anyNA(breaks) || length(breaks) < 2L ||
      is.unsorted(breaks, strictly = TRUE)) {
    rlang::abort("`breaks` must be at least two strictly increasing numbers.")
  }
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

# mrds keys detections on `object`, and refuses duplicates with a message that
# begins with a newline.
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
