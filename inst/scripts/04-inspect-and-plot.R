# 04-inspect-and-plot.R
#
# Step 4 of the basic dsfit pipeline: choose a model and plot its detection
# function against the observed distances.
#
# Input:  dsfit-output/03-sweep.rds, dsfit-output/01-detections.rds
# Output: dsfit-output/04-detection-function.png
#         dsfit-output/04-chosen-model.rds  (the chosen mrds fit)

library(dsfit)

# --- Settings ------------------------------------------------------------------

out_dir <- "dsfit-output"

# The model to carry forward. NULL takes the lowest AIC. Check its goodness of
# fit (cvm_p, or chisq_p for binned data) and its esw against the runners-up
# before accepting it.
chosen <- NULL

# --- 1. Load -------------------------------------------------------------------

sw         <- readRDS(file.path(out_dir, "03-sweep.rds"))
detections <- readRDS(file.path(out_dir, "01-detections.rds"))
tab        <- selection_table(sw)

# --- 2. Choose -----------------------------------------------------------------

if (is.null(chosen)) chosen <- tab$model_id[1]
fit <- sw$fits[[chosen]]
print(tab[tab$model_id == chosen, ])

# Models within 2 AIC of the best, for comparison.
print(tab[tab$delta_aic < 2, c("model_id", "delta_aic", "p", "esw", "cvm_p")])

# --- 3. Plot -------------------------------------------------------------------

png(file.path(out_dir, "04-detection-function.png"), width = 1600,
    height = 1200, res = 200)

# g(x) averaged over the observed covariates, with a 95% band. For a covariate
# model, fix a covariate with e.g. at = list(beaufort = 4).
est <- effect_estimates_ddf(fit, interval = "ci")
plot(est$.x, est$.estimate, type = "l", lwd = 2, ylim = c(0, 1),
     xlab = "Perpendicular distance", ylab = "Detection probability g(x)",
     main = chosen)
polygon(c(est$.x, rev(est$.x)), c(est$.lower, rev(est$.upper)),
        col = grDevices::adjustcolor("grey50", 0.3), border = NA)
# Observed distances along the axis (exact distances only; binned fits have
# intervals instead).
if (!is.null(fit$data$distance)) graphics::rug(fit$data$distance)

dev.off()

# Alternatives:
#   plot(fit)                                                   # mrds's own
#   fancyfx::plotEffects(fit, dat = detections, var = "distance")  # with rug

# --- 4. Save -------------------------------------------------------------------

saveRDS(fit, file.path(out_dir, "04-chosen-model.rds"))
