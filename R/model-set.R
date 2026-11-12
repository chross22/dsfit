#' Build a set of candidate detection functions
#'
#' Expands key functions, adjustment series and covariate formulas into a grid
#' of candidate models, one row per candidate. Pass the result to
#' [sweep_models()] to fit them all.
#'
#' @section How it works:
#' 1. Checks that every key, adjustment series, order and formula is valid.
#' 2. Crosses keys, series, orders and formulas into a grid.
#' 3. Drops, with a warning, candidates that cannot be fitted meaningfully:
#'    a uniform key with no adjustment, a uniform key with covariates, and an
#'    adjustment order outside its series.
#' 4. Builds each candidate's `mrds` model string and its `model_id` label.
#' 5. Drops, with a warning, duplicate candidates.
#'
#' @section Key functions:
#' \describe{
#'   \item{`"hn"`}{Half-normal. A shoulder at zero and a smooth fall-off. The
#'     usual default.}
#'   \item{`"hr"`}{Hazard-rate. A broader shoulder and a heavier tail, with a
#'     second parameter.}
#'   \item{`"unif"`}{Uniform. Detection is certain out to the truncation
#'     distance, which is a strip transect. It is kept only with adjustment
#'     terms. It takes no covariates, because it has no scale parameter for
#'     them to act on.}
#'   \item{`"gamma"`}{Gamma. Unimodal, with its peak away from zero, so
#'     detection is lowest on the track-line. `mrds` scales it to 1 at its
#'     peak, not at zero. Not available in `Distance::ds()`.}
#' }
#'
#' @section Gamma key and left truncation:
#' An aircraft cannot see the water directly beneath it, and a ship cannot see
#' past its bow wave. For such a platform the gamma key fits the near distances
#' better than a half-normal, which must peak at zero. Left truncation models
#' the same blind spot by discarding the near distances. Use one or the other,
#' not both. [sweep_models()] takes `left` and warns if it is set for a model
#' set that contains the gamma key.
#'
#' @param key Character vector of key functions: any of `"hn"`, `"hr"`,
#'   `"unif"`, `"gamma"`. Default `c("hn", "hr")`.
#' @param adjustment Character vector of adjustment series: `"cos"` (cosine),
#'   `"herm"` (Hermite polynomial) or `"poly"` (simple polynomial). One
#'   candidate per series. `NULL` (default) fits no adjustment terms.
#' @param order Integer vector of highest adjustment orders, one candidate per
#'   element. Orders are cumulative, as in `Distance::ds()`: a candidate carries
#'   every term of its series up to `order`. For example, `order = 3` on `"hn"`
#'   with `"cos"` fits `adj.order = c(2, 3)`. Cosine orders run 2, 3, 4, ...
#'   (1, 2, 3, ... for `"unif"`). Hermite and polynomial orders run 4, 6, 8, ...
#'   (2, 4, 6, ... for `"unif"`). `NULL` (default) gives the first term of the
#'   series only. An order outside its series is dropped with a warning.
#'   Ignored when `adjustment` is `NULL`.
#' @param formula One-sided formula, or a list of them, giving detection
#'   covariates. `~1` (default) means no covariates and fits with `cds`; any
#'   other formula fits with `mcds`.
#'
#' @return A tibble with one row per candidate and columns:
#' \describe{
#'   \item{`model_id`}{Label, e.g. `"hn"`, `"hr+cos(2,3)"`, `"hn wt_beaufort"`.}
#'   \item{`key`}{Key function.}
#'   \item{`adjustment`}{Adjustment series, or `NA`.}
#'   \item{`order`}{Highest adjustment order, or `NA`.}
#'   \item{`formula`}{Covariate formula as a string.}
#'   \item{`dsmodel`}{The model formula string passed to `mrds::ddf()`.}
#' }
#'
#' @references
#' Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L.
#' and Thomas, L. (2001) *Introduction to Distance Sampling.* Oxford University
#' Press.
#'
#' Laake, J.L. and Borchers, D.L. (2004) Methods for incomplete detection at
#' distance zero. In *Advanced Distance Sampling*, Oxford University Press.
#'
#' Becker, E.A. and Quang, P.X. (2009) A gamma-shaped detection function for
#' line-transect surveys with mark-recapture and covariate data. *Journal of
#' Agricultural, Biological and Environmental Statistics* 14:207-223.
#' \doi{10.1198/jabes.2009.0013} Source of the gamma key, developed for aerial
#' surveys.
#'
#' @seealso [sweep_models()]
#'
#' @examples
#' # Default: half-normal and hazard-rate, no adjustments or covariates
#' model_set()
#'
#' # Cosine adjustments up to orders 2 and 3
#' model_set(key = c("hn", "hr"), adjustment = "cos", order = 2:3)
#'
#' # With and without a detection covariate
#' model_set(key = "hn", formula = list(~1, ~wt_beaufort))
#'
#' # The gamma key, for a platform that cannot see beneath itself
#' model_set(key = c("hn", "gamma"))
#'
#' @export
model_set <- function(key = c("hn", "hr"), adjustment = NULL, order = NULL,
                      formula = ~1) {
  # 1. Validate the key functions
  key <- as.character(key)
  if (!is.null(adjustment)) adjustment <- as.character(adjustment)
  known_key <- c("hn", "hr", "unif", "gamma")
  bad <- setdiff(key, known_key)
  if (length(bad)) {
    rlang::abort(paste0(
      "Unknown key function", if (length(bad) > 1) "s" else "", ": ",
      paste0("\"", bad, "\"", collapse = ", "), ". Use one of ",
      paste0("\"", known_key, "\"", collapse = ", "), "."
    ))
  }

  # 2. Validate the adjustment series and orders. Orders matter only when
  #    there is a series to apply them to.
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
    if (!is.null(order) &&
        (!is.numeric(order) || anyNA(order) || any(order != round(order)))) {
      rlang::abort("`order` must be whole numbers.")
    }
  }

  # 3. Validate the covariate formulas. A single formula is wrapped in a list
  #    so one and several are handled alike.
  formulas <- if (inherits(formula, "formula")) list(formula) else formula
  if (!all(vapply(formulas, inherits, logical(1), "formula")) ||
      any(lengths(formulas) != 2L)) {
    rlang::abort(paste0(
      "`formula` must be a one-sided formula, or a list of them: `~beaufort`, ",
      "not `distance ~ beaufort`."
    ))
  }

  # 4. Cross every key, series, order and formula. Formulas enter by
  #    position (`formula_i`) so the grid holds only plain vectors.
  grid <- expand.grid(
    key = key,
    adjustment = if (is.null(adjustment)) NA_character_ else adjustment,
    order = if (is.null(adjustment) || is.null(order)) NA_integer_ else order,
    formula_i = seq_along(formulas),
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS = FALSE
  )

  # 5. Drop the uniform key without adjustments. It is a flat line (certain
  #    detection out to the truncation distance), i.e. a strip transect, not
  #    a detection function.
  drop <- grid$key == "unif" & is.na(grid$adjustment)
  if (any(drop)) {
    rlang::warn(paste0(
      "Dropping the uniform key with no adjustment terms: it is a flat ",
      "detection function, which is a strip transect rather than a model to ",
      "fit. Give it an `adjustment` to keep it."
    ))
    grid <- grid[!drop, , drop = FALSE]
  }

  # 6. Drop the uniform key with covariates. Covariates act on the scale
  #    parameter, which the uniform key lacks. mrds ignores the formula, so
  #    `unif+cos(1) bf` would fit as `unif+cos(1)` and appear twice in the
  #    table under different labels.
  form_chr <- vapply(formulas, deparse_formula, character(1))
  covariate <- !vapply(trimws(form_chr), identical, logical(1), "~1")
  drop <- grid$key == "unif" & covariate[grid$formula_i]
  if (any(drop)) {
    rlang::warn(paste0(
      "Dropping the uniform key with a covariate formula: covariates act on ",
      "the scale parameter, which the uniform key does not have, so `mrds` ",
      "ignores the formula and the fit duplicates the uniform model without ",
      "it. Use a key with a scale (\"hn\", \"hr\", \"gamma\") for covariates."
    ))
    grid <- grid[!drop, , drop = FALSE]
  }

  # 7. Expand each highest order into the full set of orders it carries
  #    (cumulative, as in Distance::ds()). Drop orders outside the series,
  #    such as cosine order 1 on a non-uniform key or an odd Hermite order:
  #    mrds cannot fit them and would fail inside the sweep.
  adj <- !is.na(grid$adjustment)
  grid$orders <- vector("list", nrow(grid))
  grid$orders[adj] <- Map(adjustment_orders, grid$key[adj],
                          grid$adjustment[adj], grid$order[adj])
  drop <- adj & vapply(grid$orders, is.null, logical(1))
  if (any(drop)) {
    bad <- unique(paste0(grid$key[drop], "+", grid$adjustment[drop],
                         grid$order[drop]))
    rlang::warn(paste0(
      "Dropping adjustment order(s) outside their series: ",
      paste(bad, collapse = ", "), ". Cosine orders run 2, 3, ... (1, 2, ... ",
      "for unif); Hermite and polynomial orders run 4, 6, ... (2, 4, ... for ",
      "unif)."
    ))
    grid <- grid[!drop, , drop = FALSE]
  }
  if (!nrow(grid)) {
    rlang::abort("No candidate models left after validation.")
  }

  # 8. Build the output: one row per candidate, with its mrds model string
  #    and its label.
  out <- tibble::tibble(
    key = grid$key,
    adjustment = grid$adjustment,
    order = vapply(grid$orders, function(o) {
      if (is.null(o)) NA_integer_ else as.integer(max(o))
    }, integer(1)),
    formula = form_chr[grid$formula_i]
  )
  out$dsmodel <- vapply(seq_len(nrow(out)), function(i) {
    build_dsmodel(out$key[i], out$adjustment[i], grid$orders[[i]],
                  out$formula[i])
  }, character(1))
  out$model_id <- vapply(seq_len(nrow(out)), function(i) {
    label_model(out$key[i], out$adjustment[i], grid$orders[[i]],
                out$formula[i])
  }, character(1))

  # 9. Drop duplicates. Repeated keys, orders or formulas give identical
  #    candidates with the same model_id; the fits list could not tell them
  #    apart, and each would tie at delta AIC 0 with its copy.
  dup <- duplicated(out$model_id)
  if (any(dup)) {
    rlang::warn(paste0(
      "Dropping duplicate candidate(s): ",
      paste(unique(out$model_id[dup]), collapse = ", "), "."
    ))
    out <- out[!dup, , drop = FALSE]
  }

  out[, c("model_id", "key", "adjustment", "order", "formula", "dsmodel")]
}

# Turn a formula into one line of text, e.g. `~wt_beaufort`.
# Input: a formula. Returns: a length-1 character string. deparse() can split
# a long formula over several strings, so they are pasted back together.
deparse_formula <- function(f) paste(deparse(f), collapse = " ")

# Build the `dsmodel` formula string passed to mrds::ddf().
# Inputs: key, adjustment series (or NA), the vector of adjustment orders, and
# the covariate formula as text. Returns a string such as
# `~cds(key = "hn", formula = ~1, adj.series = "cos", adj.order = c(2, 3))`.
# `cds` is used without covariates and `mcds` with them; mcds carries a scale
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
               paste0("adj.order = ", format_orders(order)))
  }
  paste0("~", fn, "(", paste(parts, collapse = ", "), ")")
}

# Build the human-readable model_id for one candidate.
# Inputs: as for build_dsmodel(). Returns a string such as `"hr+cos(2,3)"` or
# `"hn wt_beaufort"`: key, then adjustment orders, then covariates (no `~`).
label_model <- function(key, adjustment, order, formula) {
  out <- key
  if (!is.na(adjustment)) {
    out <- paste0(out, "+", adjustment, "(", paste(order, collapse = ","), ")")
  }
  if (!identical(trimws(formula), "~1")) {
    out <- paste0(out, " ", sub("^~", "", trimws(formula)))
  }
  out
}

# List the adjustment orders a candidate carries, following Distance's
# get_adj_orders().
# Inputs: key, adjustment series, and the highest order `max` (NA for none
# given). Returns an integer vector from the first order of the series up to
# `max`; the first order alone when `max` is NA; NULL when `max` is not a
# member of the series.
# Series: cosine starts at 2 (1 for unif) in steps of 1; Hermite and
# polynomial start at 4 (2 for unif) in steps of 2.
adjustment_orders <- function(key, adjustment, max) {
  if (adjustment == "cos") {
    start <- if (key == "unif") 1L else 2L
    step <- 1L
  } else {
    start <- if (key == "unif") 2L else 4L
    step <- 2L
  }
  if (is.na(max)) return(start)
  if (max < start || (max - start) %% step != 0) return(NULL)
  seq.int(start, max, by = step)
}

# Write adjustment orders as R code for the dsmodel string.
# Input: integer vector. Returns `"2"` for one order, `"c(2, 3)"` for several.
format_orders <- function(order) {
  if (length(order) == 1L) as.character(order) else
    paste0("c(", paste(order, collapse = ", "), ")")
}
