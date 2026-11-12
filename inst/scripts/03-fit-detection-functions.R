# 03-fit-detection-functions.R
#
# Step 3 of the basic dsfit pipeline: fit every candidate detection function
# and rank them.
#
# Input:  dsfit-output/01-detections.rds, dsfit-output/02-settings.rds
# Output: dsfit-output/03-sweep.rds            the full sweep (fits included)
#         dsfit-output/03-selection-table.csv  the ranked table

library(dsfit)

# --- Settings ------------------------------------------------------------------

out_dir <- "dsfit-output"

# Key functions: "hn" half-normal, "hr" hazard-rate, "gamma" (peak away from
# the line, for a platform with a blind spot), "unif" (needs adjustments).
keys <- c("hn", "hr")

# Adjustment terms, fitted only without covariates. NULL for none, or e.g.
# "cos" with orders up to 3: adjustment <- "cos"; order <- 2:3
adjustment <- NULL
order <- NULL

# --- 1. Load -------------------------------------------------------------------

detections <- readRDS(file.path(out_dir, "01-detections.rds"))
settings   <- readRDS(file.path(out_dir, "02-settings.rds"))

# --- 2. Build the model set --------------------------------------------------------

# Every combination of the covariates from step 2, plus no covariates: for
# three covariates that is 2^3 = 8 formulas.
covs   <- settings$covariates
combos <- unlist(
  lapply(seq_along(covs), function(k) combn(covs, k, simplify = FALSE)),
  recursive = FALSE
)
formulas <- c(list(~1), lapply(combos, reformulate))

models <- model_set(key = keys, formula = formulas)
if (!is.null(adjustment)) {
  models <- rbind(models, model_set(key = keys, adjustment = adjustment,
                                    order = order))
}
print(models[, c("model_id", "dsmodel")], n = Inf)

# --- 3. Fit --------------------------------------------------------------------

# One truncation for every candidate, so their AICs are comparable. mrds draws
# random starting values, so a seed makes a marginal fit repeat.
set.seed(1)
sw <- sweep_models(detections, models,
                   truncation = settings$truncation, left = settings$left)

# --- 4. Rank -------------------------------------------------------------------

# Sorted by AIC. Read p_cv and esw as well: two models within 2 AIC can imply
# different abundance.
tab <- selection_table(sw)
print(tab[, c("model_id", "delta_aic", "p", "p_cv", "esw", "cvm_p", "chisq_p")],
      n = Inf)

# --- 5. Save -------------------------------------------------------------------

saveRDS(sw, file.path(out_dir, "03-sweep.rds"))
utils::write.csv(selection_table(sw, converged_only = FALSE),
                 file.path(out_dir, "03-selection-table.csv"), row.names = FALSE)
