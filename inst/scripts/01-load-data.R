# 01-load-data.R
#
# Step 1 of the basic dsfit pipeline: load the detections.
#
# Input:  one row per detection, with either
#           - `distance`              exact perpendicular distances, or
#           - `distbegin`, `distend`  binned distances (one interval per row),
#         plus any detection covariates (e.g. `beaufort`).
#         The flatfile from `distsamp::detection_data()` has this shape.
# Output: dsfit-output/01-detections.rds
#
# Run the scripts in order, 01 to 05, from the same working directory. Each one
# reads what the one before it saved to `dsfit-output/`.

library(dsfit)

# --- Settings ------------------------------------------------------------------

out_dir <- "dsfit-output"

# --- 1. Read the detections ----------------------------------------------------

# Replace this block with your own data. Two common sources:
#
#   detections <- read.csv("detections.csv")
#
#   segs       <- distsamp::segment_survey(dat, seg_length = 5, seed = 1)
#   detections <- distsamp::detection_data(segs, area = 5811)
#
# The example below simulates a survey so the pipeline runs as shipped: 300
# detections from a half-normal whose scale falls with sea state, plus 20 rows
# for samples with no detection (as a flatfile records them).

set.seed(1)
n <- 300
detections <- data.frame(
  object   = seq_len(n),
  beaufort = sample(0:4, n, replace = TRUE),
  glare    = sample(c("none", "some"), n, replace = TRUE),
  size     = rpois(n, 1) + 1
)
detections$distance <- abs(rnorm(n, mean = 0, sd = 150 - 15 * detections$beaufort))
detections <- rbind(
  detections,
  data.frame(object = n + 1:20, beaufort = NA, glare = NA, size = NA,
             distance = NA)
)

# --- 2. Save -------------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE)
saveRDS(detections, file.path(out_dir, "01-detections.rds"))
message("Saved ", nrow(detections), " rows to ",
        file.path(out_dir, "01-detections.rds"))
