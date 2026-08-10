#' @export
print.dsfit_sweep <- function(x, ...) {
  tab <- x$table
  ok <- tab[tab$converged, , drop = FALSE]

  cat("<dsfit_sweep>\n")
  cat("  detections:  ", x$settings$n,
      if (x$settings$binned) " (binned)" else "", "\n", sep = "")
  cat("  truncation:  ", x$settings$truncation,
      if (!is.null(x$settings$left)) paste0("   left: ", x$settings$left) else "",
      "\n", sep = "")
  cat("  models:      ", nrow(ok), " of ", nrow(tab), " fitted\n", sep = "")

  if (nrow(ok)) {
    show <- utils::head(ok, 5)
    cat("\n")
    out <- data.frame(
      model = show$model_id,
      dAIC = round(show$delta_aic, 2),
      p = round(show$p, 4),
      p_cv = round(show$p_cv, 3),
      esw = round(show$esw, 1),
      CvM_p = round(show$cvm_p, 3)
    )
    print(out, row.names = FALSE)
    if (nrow(ok) > 5) cat("  ... and ", nrow(ok) - 5, " more\n", sep = "")
    cat("\n  Rank on esw and p_cv as well as dAIC: models within 2 AIC can\n")
    cat("  imply materially different abundance. g(0) = 1 is assumed.\n")
  }
  invisible(x)
}
