# 02-check-data.R
#
# Step 2 of the basic dsfit pipeline: check the data and choose the truncation.
# Nothing is fitted here.
#
# Input:  dsfit-output/01-detections.rds
# Output: dsfit-output/02-settings.rds  (truncation, left, covariates)

library(dsfit)

# --- Settings ------------------------------------------------------------------

out_dir <- "dsfit-output"

# Right truncation. NULL uses the suggestion below: the distance that keeps the
# nearest 95% of detections (Buckland et al. 2001 suggest trimming about 5%).
truncation <- NULL

# Left truncation, for a platform that cannot see directly beneath itself.
# NULL for none. Use this or the gamma key in step 3, not both.
left <- NULL

# Detection covariates to consider in step 3. Every combination is fitted.
covariates <- c("beaufort", "glare", "size")

# --- 1. Load -------------------------------------------------------------------

detections <- readRDS(file.path(out_dir, "01-detections.rds"))

# --- 2. What the data can support ----------------------------------------------

# Exact or binned distances, how many, which columns could be covariates, and
# whether perception bias or group size can be handled.
print(detection_structure(detections))

# --- 3. Choose the truncation ----------------------------------------------------

if (is.null(truncation)) {
  d <- if ("distance" %in% names(detections) &&
           any(!is.na(detections$distance))) {
    detections$distance
  } else {
    detections$distend
  }
  truncation <- signif(stats::quantile(d, 0.95, na.rm = TRUE), 2)
  truncation <- unname(truncation)
  message("Suggested truncation: ", truncation)
}

# For binned data the truncation must equal a bin edge; set it explicitly above.

# --- 4. Check the planned sweep ---------------------------------------------------

# Runs the same input checks as the sweep, reports what truncation drops,
# looks for a blind spot near the line, and checks the covariates exist.
formulas <- c(list(~1), lapply(covariates, reformulate))
diagnose_sweep(detections, model_set(c("hn", "hr"), formula = formulas),
               truncation = truncation, left = left)

# --- 5. Save the settings for step 3 ----------------------------------------------

saveRDS(list(truncation = truncation, left = left, covariates = covariates),
        file.path(out_dir, "02-settings.rds"))
