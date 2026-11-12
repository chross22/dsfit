# 05-g0-correction.R
#
# Step 5 of the basic dsfit pipeline: build the g(0) correction to hand to the
# abundance step. dsfit does not apply it.
#
# Every fit in step 3 assumes g(0) = 1: an animal on the trackline is always
# seen. For whales that is false. The correction has separate components:
#   - availability: was the animal at the surface while in view?
#     Computed here from dive data with availability().
#   - perception: was an available animal seen? Needs a double-observer
#     study; usually taken from another source.
#
# Input:  dive intervals and platform geometry (settings below)
# Output: dsfit-output/05-g0.rds

library(dsfit)

# --- Settings ------------------------------------------------------------------

out_dir <- "dsfit-output"

# Dive data: mean surfacing and diving intervals (s) with standard errors, one
# row per key (month here). The shipped example is invented; use your own.
dives <- example_dive_intervals

# Time an animal stays in view (s). Aircraft side window: time on the
# trackline plus seconds gained per metre out, here at 500 m.
window <- view_window_aerial(trackline = 51.22, speed = 51.4, slope = 0.03,
                             distance = 500)

# Perception: value and SE from a double-observer study, and where it came
# from. NA key applies it to every month.
perception <- data.frame(
  key = NA, component = "perception", value = 0.68, se = 0.09,
  source = "replace with your source"
)

# --- 1. Availability -------------------------------------------------------------

avail <- availability(
  surface    = dives$surface,
  dive       = dives$dive,
  window     = window,
  se_surface = dives$se_surface,
  se_dive    = dives$se_dive,
  key        = as.character(dives$month)
)
avail$source <- "example_dive_intervals (invented)"
print(avail)

# --- 2. Combine ----------------------------------------------------------------

# Values multiply; squared CVs add.
correction <- g0(avail, perception)
print(correction)

# --- 3. Save -------------------------------------------------------------------

# At the abundance step: divide abundance by correction$table$value for the
# matching key, and add correction$table$cv in quadrature to abundance's CV.
saveRDS(correction, file.path(out_dir, "05-g0.rds"))
