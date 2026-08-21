# Package index

## Fitting and comparing

Fit a set of candidate keys at once and compare them on effective strip
half-width and its coefficient of variation, rather than on AIC alone.

- [`diagnose_sweep()`](https://camilleross.org/dsfit/reference/diagnose_sweep.md)
  : Diagnose common reasons a sweep fails, or succeeds meaninglessly
- [`model_set()`](https://camilleross.org/dsfit/reference/model_set.md)
  : Build a set of candidate detection functions
- [`selection_table()`](https://camilleross.org/dsfit/reference/selection_table.md)
  : The ranked selection table from a sweep
- [`sweep_models()`](https://camilleross.org/dsfit/reference/sweep_models.md)
  : Fit a set of detection functions and compare them

## Preparing survey data

Getting distance data into the shape a detection function wants, and
describing the window an observer could actually see through.

- [`detection_structure()`](https://camilleross.org/dsfit/reference/detection_structure.md)
  : What a set of detections can and cannot support
- [`prepare_distance_data()`](https://camilleross.org/dsfit/reference/prepare_distance_data.md)
  : Check and prepare distance data for fitting
- [`view_window()`](https://camilleross.org/dsfit/reference/view_window.md)
  : The time an animal stays within view of the platform
- [`view_window_aerial()`](https://camilleross.org/dsfit/reference/view_window_aerial.md)
  : The time an animal stays in view of an aircraft

## Detection on the line

Availability and perception bias: what was there, versus what could have
been seen at all.

- [`availability()`](https://camilleross.org/dsfit/reference/availability.md)
  : Availability: the probability an animal was at the surface to be
  seen
- [`effect_estimates(`*`<ddf>`*`)`](https://camilleross.org/dsfit/reference/effect_estimates_ddf.md)
  [`effect_estimates_ddf()`](https://camilleross.org/dsfit/reference/effect_estimates_ddf.md)
  : Extract a fitted detection function as a tidy frame
- [`g0()`](https://camilleross.org/dsfit/reference/g0.md) : Assemble a
  g(0) correction from its named components

## Bundled data

Dive-interval and detection data shipped with the package, so the
examples and the availability correction run without needing an extract.

- [`example_dive_intervals`](https://camilleross.org/dsfit/reference/example_dive_intervals.md)
  : Example surfacing and diving intervals, by month
- [`ganley_availability`](https://camilleross.org/dsfit/reference/ganley_availability.md)
  : Right whale availability in Cape Cod Bay, by month
- [`ganley_detection`](https://camilleross.org/dsfit/reference/ganley_detection.md)
  : Right whale annual detection probability in Cape Cod Bay
