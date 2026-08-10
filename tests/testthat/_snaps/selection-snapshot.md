# the ranking over key functions is stable

    Code
      snap_table(sw)
    Output
        model converged n_par     aic delta_aic      p  p_cv   esw cvm_p chisq_p
      1    hn      TRUE     1 2177.06      0.00 0.3494 0.052 139.8 0.928   0.613
      2 gamma      TRUE     2 2186.97      9.92 0.3010 0.074 120.4 0.498   0.339
      3    hr      TRUE     2 2187.24     10.18 0.3781 0.064 151.3 0.786   0.126

# the ranking over adjustment series and orders is stable

    Code
      snap_table(sw)
    Output
             model converged n_par     aic delta_aic      p  p_cv   esw cvm_p chisq_p
      1    hn+cos3      TRUE     2 2178.51      0.00 0.3278 0.099 131.1 0.916   0.634
      2    hn+cos2      TRUE     2 2179.05      0.54 0.3488 0.084 139.5 0.932   0.530
      3 unif+herm2      TRUE     1 2272.10     93.59 0.6674 0.082 267.0 0.001   0.000
      4  unif+cos2      TRUE     1 2379.98    201.46 0.7043 0.064 281.7 0.012   0.000
      5  unif+cos3      TRUE     1 2396.71    218.20 0.8776 0.089 351.0 0.020   0.000
      6   hn+herm2     FALSE    NA      NA        NA     NA    NA    NA    NA      NA
      7   hn+herm3     FALSE    NA      NA        NA     NA    NA    NA    NA      NA
      8 unif+herm3     FALSE    NA      NA        NA     NA    NA    NA    NA      NA

# the binned likelihood is stable

    Code
      snap_table(sw)
    Output
        model converged n_par    aic delta_aic      p  p_cv   esw cvm_p chisq_p
      1    hn      TRUE     1 329.85      0.00 0.4565 0.057 137.0    NA   0.471
      2    hr      TRUE     2 331.33      1.48 0.4940 0.077 148.2    NA      NA

# covariate models are stable, and mcds is the path taken

    Code
      snap_table(sw)
    Output
        model converged n_par     aic delta_aic      p  p_cv   esw cvm_p chisq_p
      1 hn bf      TRUE     2 2210.39      0.00 0.3879 0.055 155.1 0.724   0.490
      2    hn      TRUE     1 2220.34      9.95 0.4035 0.047 161.4 0.497   0.386

# left truncation is stable

    Code
      snap_table(sw)
    Output
        model converged n_par     aic delta_aic      p  p_cv   esw cvm_p chisq_p
      1    hn      TRUE     1 1695.39      0.00 0.2959 0.071 118.4 0.809   0.683
      2    hr      TRUE     2 1705.36      9.97 0.3330 0.088 133.2 0.709   0.179

# what the sweep reports about itself is stable

    Code
      sw <- sweep_models(d, model_set(key = c("hn", "hr", "gamma")), truncation = 400)
    Message
      Fitted 3 of 3 candidate models to 200 detections, truncation 400.
        3 row(s) dropped: no distance, or beyond the truncation.
        g(0) = 1 is assumed; see `?sweep_models`.

---

    Code
      print(sw)
    Output
      <dsfit_sweep>
        detections:  200
        truncation:  400
        models:      3 of 3 fitted
      
       model  dAIC      p  p_cv   esw CvM_p
          hn  0.00 0.3494 0.052 139.8 0.928
       gamma  9.92 0.3010 0.074 120.4 0.498
          hr 10.18 0.3781 0.064 151.3 0.786
      
        Rank on esw and p_cv as well as dAIC: models within 2 AIC can
        imply materially different abundance. g(0) = 1 is assumed.

