#' Build a set of candidate detection functions
#'
#' Expands key functions, adjustment terms, and covariate formulas into a grid
#' of candidate models, one row each, ready for [sweep_models()].
#'
#' @section The key functions:
#' \describe{
#'   \item{`"hn"`}{Half-normal. A shoulder at zero and a smooth fall-off; the
#'     usual default.}
#'   \item{`"hr"`}{Hazard-rate. A broader shoulder and a heavier tail, at the
#'     cost of a second parameter.}
#'   \item{`"unif"`}{Uniform. Detection certain out to the truncation distance —
#'     a strip transect. Only interesting with adjustment terms.}
#'   \item{`"gamma"`}{Gamma. **Unimodal, with its peak away from zero**, so
#'     detection is lowest on the track-line. Not available in
#'     `Distance::ds()`, which is one reason this package fits through `mrds`
#'     directly.}
#' }
#'
#' @section When gamma is the right shape, and when it is a trap:
#' An aircraft cannot see the water directly beneath it. Where that is true the
#' gamma key describes the survey rather than distorting it, and a half-normal
#' forced through a peak at zero will misfit the near distances.
#'
#' But a gamma key and a **left truncation model the same phenomenon**. Using
#' both accounts for the blind spot twice. Pick one deliberately and record
#' which — [sweep_models()] takes `left` for the other route, and will warn if a
#' model set containing gamma is fitted with one.
#'
#' Note also that gamma does not assume `g(0) = 1`, so its `p̄` is not
#' comparable to a half-normal's in the way a shared assumption would make it.
#' That is a reason to read the whole selection table rather than its first row.
#'
#' @param key Character vector of key functions: any of `"hn"`, `"hr"`,
#'   `"unif"`, `"gamma"`.
#' @param adjustment Adjustment series: `NULL` (none), `"cos"`, `"herm"`, or
#'   `"poly"`. A vector expands to one candidate per series.
#' @param order Integer vector of adjustment orders. Ignored when `adjustment`
#'   is `NULL`.
#' @param formula One-sided formula, or a list of them, giving detection
#'   covariates. `~1` (default) is no covariates. Anything else fits through
#'   `mcds` rather than `cds`.
#'
#' @return A tibble with one row per candidate: `model_id`, `key`,
#'   `adjustment`, `order`, `formula`, and the `dsmodel` string `mrds` is given.
#'
#' @references
#' Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L.
#' and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford University
#' Press.
#'
#' Laake, J.L. and Borchers, D.L. (2004) Methods for incomplete detection at
#' distance zero. In *Advanced Distance Sampling*, Oxford University Press.
#'
#' @seealso [sweep_models()]
#'
#' @examples
#' model_set()
#'
#' # With adjustments
#' model_set(key = c("hn", "hr"), adjustment = "cos", order = 2:3)
#'
#' # With a detection covariate
#' model_set(key = "hn", formula = list(~1, ~wt_beaufort))
#'
#' # The unimodal key, for a platform with a blind spot beneath it
#' model_set(key = c("hn", "gamma"))
#'
#' @export
model_set <- function(key = c("hn", "hr"), adjustment = NULL, order = NULL,
                      formula = ~1) {
  known_key <- c("hn", "hr", "unif", "gamma")
  bad <- setdiff(key, known_key)
  if (length(bad)) {
    rlang::abort(paste0(
      "Unknown key function", if (length(bad) > 1) "s" else "", ": ",
      paste0("\"", bad, "\"", collapse = ", "), ". Use one of ",
      paste0("\"", known_key, "\"", collapse = ", "), "."
    ))
  }

  known_adj <- c("cos", "herm", "poly")
  if (!is.null(adjustment)) {
    bad <- setdiff(adjustment, known_adj)
    if (length(bad)) {
      rlang::abort(paste0(
        "Unknown adjustment series: ",
        paste0("\"", bad, "\"", collapse = ", "), ". Use one of ",
        paste0("\"", known_adj, "\"", collapse = ", "), "."
      ))
    }
    if (is.null(order)) order <- 2L
  }

  formulas <- if (inherits(formula, "formula")) list(formula) else formula
  if (!all(vapply(formulas, inherits, logical(1), "formula"))) {
    rlang::abort("`formula` must be a one-sided formula, or a list of them.")
  }

  grid <- expand.grid(
    key = key,
    adjustment = if (is.null(adjustment)) NA_character_ else adjustment,
    order = if (is.null(adjustment)) NA_integer_ else order,
    formula_i = seq_along(formulas),
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS = FALSE
  )

  # A uniform key with no adjustment is a flat line: detection certain out to
  # the truncation distance. That is a strip transect, not a detection
  # function, and fitting it as one is almost never what was meant.
  drop <- grid$key == "unif" & is.na(grid$adjustment)
  if (any(drop)) {
    rlang::warn(paste0(
      "Dropping the uniform key with no adjustment terms: it is a flat ",
      "detection function, which is a strip transect rather than a model to ",
      "fit. Give it an `adjustment` to keep it."
    ))
    grid <- grid[!drop, , drop = FALSE]
  }
  if (!nrow(grid)) {
    rlang::abort("No candidate models left after validation.")
  }

  form_chr <- vapply(formulas, deparse_formula, character(1))
  out <- tibble::tibble(
    key = grid$key,
    adjustment = grid$adjustment,
    order = as.integer(grid$order),
    formula = form_chr[grid$formula_i]
  )
  out$dsmodel <- vapply(seq_len(nrow(out)), function(i) {
    build_dsmodel(out$key[i], out$adjustment[i], out$order[i], out$formula[i])
  }, character(1))
  out$model_id <- vapply(seq_len(nrow(out)), function(i) {
    label_model(out$key[i], out$adjustment[i], out$order[i], out$formula[i])
  }, character(1))

  out[, c("model_id", "key", "adjustment", "order", "formula", "dsmodel")]
}

deparse_formula <- function(f) paste(deparse(f), collapse = " ")

# The dsmodel string mrds is handed. `cds` when there are no covariates,
# `mcds` when there are - they are the same family, but mcds carries a scale
# formula and cds does not.
build_dsmodel <- function(key, adjustment, order, formula) {
  covariate <- !identical(trimws(formula), "~1")
  fn <- if (covariate) "mcds" else "cds"

  parts <- c(
    paste0("key = \"", key, "\""),
    paste0("formula = ", formula)
  )
  if (!is.na(adjustment)) {
    parts <- c(parts,
               paste0("adj.series = \"", adjustment, "\""),
               paste0("adj.order = ", order))
  }
  paste0("~", fn, "(", paste(parts, collapse = ", "), ")")
}

label_model <- function(key, adjustment, order, formula) {
  out <- key
  if (!is.na(adjustment)) out <- paste0(out, "+", adjustment, order)
  if (!identical(trimws(formula), "~1")) {
    out <- paste0(out, " ", sub("^~", "", trimws(formula)))
  }
  out
}
