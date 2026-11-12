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
                    model converged n_par     aic delta_aic      p  p_cv   esw cvm_p
      1     unif+cos(1,2)      TRUE     2 2178.68      0.00 0.3707 0.040 148.3 0.527
      2         hn+cos(2)      TRUE     2 2179.05      0.38 0.3489 0.084 139.6 0.931
      3        hn+herm(4)      TRUE     2 2179.06      0.38 0.3499 0.073 140.0 0.925
      4       hn+cos(2,3)      TRUE     3 2179.71      1.03 0.3243 0.108 129.7 0.909
      5 unif+cos(1,2,3,4)      TRUE     4 2180.35      1.67 0.3338 0.093 133.5 0.959
      6   unif+cos(1,2,3)      TRUE     3 2180.49      1.81 0.3615 0.104 144.6 0.744
      7     hn+cos(2,3,4)      TRUE     4 2180.71      2.03 0.3402 0.127 136.1 0.975
      8    unif+herm(2,4)      TRUE     2 2199.47     20.79 0.4862 0.066 194.5 0.000
      9      unif+herm(2)      TRUE     1 2272.34     93.66 0.6679 0.082 267.2 0.001
        chisq_p
      1   0.462
      2   0.530
      3   0.531
      4   0.598
      5   0.601
      6   0.435
      7   0.591
      8   0.004
      9   0.000

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
      1    hn      TRUE     1 1695.39      0.00 0.2959 0.071 109.5 0.809   0.683
      2    hr      TRUE     2 1705.36      9.97 0.3330 0.088 123.2 0.709   0.179

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

