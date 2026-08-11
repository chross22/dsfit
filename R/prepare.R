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
#' }
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
  stopifnot(is.numeric(truncation), length(truncation) == 1L, truncation > 0)
  if (!is.null(left)) {
    stopifnot(is.numeric(left), length(left) == 1L, left >= 0)
    if (left >= truncation) {
      rlang::abort("`left` must be smaller than `truncation`.")
    }
  }

  has_point <- "distance" %in% names(data) && any(!is.na(data$distance))
  has_bin <- all(c("distbegin", "distend") %in% names(data)) &&
    any(!is.na(data$distbegin))

  if (!has_point && !has_bin) {
    rlang::abort(paste0(
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

    if (is.null(breaks)) breaks <- derive_breaks(data)

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

    return(list(data = ensure_object(data), binned = TRUE, breaks = breaks,
                n_dropped = n0 - nrow(data), dropped = dropped))
  }

  keep <- !is.na(data$distance) & data$distance <= truncation
  if (!is.null(left)) keep <- keep & data$distance >= left
  data <- data[keep, , drop = FALSE]
  if (!nrow(data)) {
    rlang::abort("No detections left after truncation.")
  }

  list(data = ensure_object(data), binned = FALSE, breaks = NULL,
       n_dropped = n0 - nrow(data), dropped = dropped)
}

# Cutpoints implied by the bins present. They have to tile: a gap or an overlap
# means the intervals do not define a single set of breaks, and deriving one
# anyway would put detections in the wrong bin.
derive_breaks <- function(data) {
  edges <- sort(unique(c(data$distbegin, data$distend)))
  bins <- unique(data.frame(b = data$distbegin, e = data$distend))
  bins <- bins[order(bins$b), , drop = FALSE]

  if (nrow(bins) > 1) {
    gaps <- bins$b[-1] != bins$e[-nrow(bins)]
    if (any(gaps)) {
      rlang::abort(paste0(
        "The bins do not tile: interval ", nrow(bins) - sum(gaps),
        " ends at ", bins$e[which(gaps)[1]], " and the next begins at ",
        bins$b[which(gaps)[1] + 1], ". Supply `breaks` explicitly."
      ))
    }
  }
  edges
}

# mrds keys detections on `object`.
ensure_object <- function(data) {
  if (!"object" %in% names(data)) {
    data$object <- seq_len(nrow(data))
  }
  data
}
