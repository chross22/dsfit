# Print method for the result of sweep_models().
# Input: a dsfit_sweep. Prints the number of detections, the truncation, how
# many candidates fitted, and the top five converged candidates by AIC with
# delta AIC, p, p_cv, esw and one goodness-of-fit p-value. Returns `x`
# invisibly. The full table is available from selection_table().
#' @export
print.dsfit_sweep <- function(x, ...) {
  tab <- x$table
  ok <- tab[tab$converged, , drop = FALSE]

  # 1. Header: detections, truncation, number of models fitted.
  cat("<dsfit_sweep>\n")
  cat("  detections:  ", x$settings$n,
      if (x$settings$binned) " (binned)" else "", "\n", sep = "")
  cat("  truncation:  ", x$settings$truncation,
      if (!is.null(x$settings$left)) paste0("   left: ", x$settings$left) else "",
      "\n", sep = "")
  cat("  models:      ", nrow(ok), " of ", nrow(tab), " fitted\n", sep = "")

  if (nrow(ok)) {
    # 2. Table of the top five converged candidates, rounded for display.
    show <- utils::head(ok, 5)
    cat("\n")
    out <- data.frame(
      model = show$model_id,
      dAIC = round(show$delta_aic, 2),
      p = round(show$p, 4),
      p_cv = round(show$p_cv, 3),
      esw = round(show$esw, 1)
    )
    # 3. Add one goodness-of-fit column. Binned distances have no
    #    Cramer-von Mises value, so show chi-square when CvM is missing and
    #    chi-square is available; otherwise show CvM.
    if (all(is.na(show$cvm_p)) && any(!is.na(show$chisq_p))) {
      out$chisq_p <- round(show$chisq_p, 3)
    } else {
      out$CvM_p <- round(show$cvm_p, 3)
    }
    print(out, row.names = FALSE)

    # 4. Count of candidates not shown, then a reminder of what to read.
    if (nrow(ok) > 5) cat("  ... and ", nrow(ok) - 5, " more\n", sep = "")
    cat("\n  Rank on esw and p_cv as well as dAIC: models within 2 AIC can\n")
    cat("  imply materially different abundance. g(0) = 1 is assumed.\n")
  }
  invisible(x)
}
