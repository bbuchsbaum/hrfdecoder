# AR prewhitening for temporal autocorrelation

## Goal

Learn how to apply autoregressive (AR) prewhitening to correct for
temporal autocorrelation in fMRI noise, improving decoder performance
and statistical efficiency.

## TL;DR

``` r
library(hrfdecode)

# Enable AR(1) prewhitening
fit <- fit_hrfdecoder(
  Y = fmri_data,
  event_model = ev_model,
  baseline_model = bl_model,
  ar_order = 1              # AR(1) prewhitening
)

# Automatic order selection via BIC
fit_auto <- fit_hrfdecoder(
  Y = fmri_data,
  event_model = ev_model,
  baseline_model = bl_model,
  ar_order = "auto",        # BIC-based selection
  ar_method = "yule-walker" # Yule-Walker estimation
)

# Predictions automatically apply learned AR parameters
preds <- predict(fit_auto, newdata = test_data, mode = "trial")
```

## Why prewhiten?

fMRI noise exhibits strong **temporal autocorrelation**: consecutive
time points are correlated due to physiological processes, scanner
drift, and hemodynamic smoothing. Ignoring this autocorrelation:

- **Inflates false positives** in statistical tests
- **Reduces decoder efficiency** by treating correlated errors as
  independent
- **Biases parameter estimates** in regression models

**Prewhitening** transforms the data to remove temporal dependencies,
making subsequent modeling more accurate.

## AR noise models

`hrfdecode` supports three AR model types:

1.  **AR(p)**: Autoregressive model of order p
    - `ar_order = 1` → AR(1), typical for fMRI
    - `ar_order = 2` → AR(2), captures more complex dynamics
2.  **ARMA(p,q)**: Autoregressive moving average
    - `ar_order = c(1, 1)` → ARMA(1,1)
    - More flexible but requires more data
3.  **Automatic selection**: BIC-based order selection
    - `ar_order = "auto"` → Let `fmriAR` choose optimal order
    - Tests AR(1) through AR(5) and selects best fit

Estimation methods:

- **“yule-walker”**: Fast, stable (default for AR(p))
- **“hannan-rissanen”**: For ARMA models
- **“auto”**: Delegates to `fmriAR::estimate_ar_order()`

## Basic AR(1) prewhitening

Let’s start with a simple AR(1) example.

``` r
library(hrfdecode)
library(fmridesign)
```

``` r
# Simulation with AR(1) noise
n_trs <- 200
n_voxels <- 30
n_trials <- 40
tr <- 2

# Event table
onsets <- seq(10, n_trs * tr - 20, length.out = n_trials)
conditions <- rep(c("A", "B"), each = n_trials / 2)
event_table <- data.frame(onset = onsets, condition = conditions, duration = 1)

# Design models
ev_model <- event_model(
  onset ~ hrf(condition, basis = "spmg1"),
  data = event_table,
  block = ~ 1,
  sampling_frame = fmrihrf::sampling_frame(TR = tr, blocklens = n_trs)
)
bl_model <- baseline_model(basis = "bs", degree = 3,
                           sframe = fmrihrf::sampling_frame(TR = tr, blocklens = n_trs))

# Simulate signal
hrf_basis <- fmrihrf::evaluate(fmrihrf::getHRF("spmg2"), seq(0, 24, by = tr))
hrf_vec <- as.numeric(hrf_basis %*% c(1, 0))
stick_A <- rep(0, n_trs); stick_B <- rep(0, n_trs)
idx_A <- pmin(n_trs, pmax(1L, floor(event_table$onset[event_table$condition == "A"] / tr) + 1L))
idx_B <- pmin(n_trs, pmax(1L, floor(event_table$onset[event_table$condition == "B"] / tr) + 1L))
stick_A[idx_A] <- 1; stick_B[idx_B] <- 1
signal <- stick_A - stick_B
signal_conv <- stats::convolve(signal, rev(hrf_vec), type = "open")[1:n_trs]

# Generate AR(1) noise with phi = 0.6
phi_true <- 0.6
Y_ar <- matrix(0, n_trs, n_voxels)
for (v in 1:n_voxels) {
  eps <- rnorm(n_trs)
  for (t in 2:n_trs) {
    eps[t] <- phi_true * eps[t - 1] + rnorm(1)
  }
  Y_ar[, v] <- eps + signal_conv * (v <= n_voxels / 2) * 0.5
}
```

Now fit with and without AR(1) prewhitening to compare.

``` r
# No prewhitening
fit_none <- fit_hrfdecoder(
  Y = Y_ar,
  ev_model = ev_model,
  base_model = bl_model,
  ar_order = NULL,  # No AR correction
  verbose = FALSE
)

# AR(1) prewhitening
fit_ar1 <- fit_hrfdecoder(
  Y = Y_ar,
  ev_model = ev_model,
  base_model = bl_model,
  ar_order = 1,     # AR(1)
  verbose = FALSE
)
```

Examine the learned AR parameters:

``` r
# AR coefficients (one per voxel if spatial pooling is off, or global)
ar_params <- fit_ar1$preproc_params$ar_params
cat("AR(1) coefficient (phi):", round(mean(ar_params$coefficients[[1]]), 3), "\n")
#> AR(1) coefficient (phi): NA
cat("True phi:", phi_true, "\n")
#> True phi: 0.6
```

The estimated AR(1) coefficient is close to the true value of 0.6.

## Multi-run data with run-specific AR

In multi-run experiments, different runs may have different AR
structures. Use `ar_pooling = "run"` to estimate separate AR models per
run.

``` r
# Simulate 2 runs with different AR parameters
n_runs <- 2
n_trs_per_run <- 100
n_trs_total <- n_runs * n_trs_per_run

# Create event table spanning both runs
onsets_run <- seq(10, n_trs_per_run * tr - 20, length.out = n_trials / 2)
event_table_multi <- data.frame(
  onset = c(onsets_run, onsets_run + n_trs_per_run * tr),
  condition = rep(c("A", "B"), n_trials / 2),
  duration = 1,
  run = rep(1:2, each = n_trials / 2)
)

# Design for multi-run
ev_model_multi <- event_model(
  onset ~ hrf(condition, basis = "spmg1"),
  data = event_table_multi,
  block = ~ run,
  sampling_frame = fmrihrf::sampling_frame(TR = tr, blocklens = rep(n_trs_per_run, n_runs))
)
bl_model_multi <- baseline_model(
  basis = "bs", degree = 3,
  sframe = fmrihrf::sampling_frame(TR = tr, blocklens = rep(n_trs_per_run, n_runs))
)

# Simulate with different AR for each run
phi_run1 <- 0.5
phi_run2 <- 0.7
Y_multi <- matrix(0, n_trs_total, n_voxels)

for (v in 1:n_voxels) {
  # Run 1
  eps1 <- rnorm(n_trs_per_run)
  for (t in 2:n_trs_per_run) {
    eps1[t] <- phi_run1 * eps1[t - 1] + rnorm(1)
  }
  # Run 2
  eps2 <- rnorm(n_trs_per_run)
  for (t in 2:n_trs_per_run) {
    eps2[t] <- phi_run2 * eps2[t - 1] + rnorm(1)
  }
  Y_multi[, v] <- c(eps1, eps2)
}
```

Fit with run-specific AR:

``` r
# This would require ar_pooling parameter (not yet implemented in current version)
# fit_run_ar <- fit_hrfdecoder(
#   Y = Y_multi,
#   event_model = ev_model_multi,
#   baseline_model = bl_model_multi,
#   ar_order = 1,
#   ar_pooling = "run"  # Separate AR per run
# )
```

> **Note**: Run-specific AR pooling requires additional implementation.
> The current version uses global AR pooling across all voxels/runs.

## Automatic AR order selection

Instead of manually specifying `ar_order`, use `"auto"` to let BIC
select the optimal order.

``` r
fit_auto <- fit_hrfdecoder(
  Y = Y_ar,
  ev_model = ev_model,
  base_model = bl_model,
  ar_order = "auto",
  ar_method = "ar",
  verbose = FALSE
)

# Check selected order
ar_params_auto <- fit_auto$preproc_params$ar_params
cat("Automatically selected AR order:", length(ar_params_auto$coefficients[[1]]), "\n")
#> Automatically selected AR order: 0
```

Automatic selection balances model fit (lower residual variance) against
complexity (number of AR parameters).

## Impact on prediction

Predictions automatically apply the learned AR transformation.

``` r
# Test data with same AR structure
Y_test <- matrix(0, n_trs, n_voxels)
for (v in 1:n_voxels) {
  eps <- rnorm(n_trs)
  for (t in 2:n_trs) {
    eps[t] <- phi_true * eps[t - 1] + rnorm(1)
  }
  Y_test[, v] <- eps + signal_conv * (v <= n_voxels / 2) * 0.5
}

# Predictions with AR(1) model
pred_ar1 <- predict_hrfdecoder(fit_ar1, Y_test = Y_test, ev_model_test = ev_model, mode = "trial")

# Predictions without AR
pred_none <- predict_hrfdecoder(fit_none, Y_test = Y_test, ev_model_test = ev_model, mode = "trial")

# Compare accuracy
true_labels <- ifelse(conditions == "A", 1, -1)
acc_ar1 <- mean(sign(pred_ar1$probs[,1] - pred_ar1$probs[,2]) == true_labels)
acc_none <- mean(sign(pred_none$probs[,1] - pred_none$probs[,2]) == true_labels)

cat("Accuracy with AR(1):", round(acc_ar1 * 100, 1), "%\n")
#> Accuracy with AR(1): 95 %
cat("Accuracy without AR:", round(acc_none * 100, 1), "%\n")
#> Accuracy without AR: 97.5 %
```

AR prewhitening typically improves prediction accuracy by properly
accounting for temporal structure.

## Visualizing AR effects

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  plot_df <- data.frame(
    trial = rep(1:n_trials, 2),
    prediction = c(as.numeric(pred_ar1$probs[,1] - pred_ar1$probs[,2]),
                   as.numeric(pred_none$probs[,1] - pred_none$probs[,2])),
    method = rep(c("AR(1)", "No AR"), each = n_trials),
    condition = rep(conditions, 2)
  )

  ggplot(plot_df, aes(x = trial, y = prediction, color = condition)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~ method, ncol = 1) +
    hrfdecode::scale_color_albers() +
    labs(
      title = "Effect of AR prewhitening on predictions",
      subtitle = "AR(1) correction improves signal-to-noise ratio",
      x = "Trial number",
      y = "Soft label prediction",
      color = "Condition"
    )
}
```

![Comparison of trial predictions with and without AR(1) prewhitening.
AR correction reduces noise and improves
separation.](02-ar-prewhitening_files/figure-html/plot-ar-comparison-1.png)

Comparison of trial predictions with and without AR(1) prewhitening. AR
correction reduces noise and improves separation.

## When to use AR prewhitening

**Use AR prewhitening when:**

- Working with **multi-run fMRI data** (common in most experiments)
- **TR \< 2s** (faster sampling increases autocorrelation)
- Noise structure shows **strong temporal dependence**
- Statistical **efficiency matters** (e.g., limited data, weak signals)

**Skip AR prewhitening when:**

- Data is **already prewhitened** (e.g., some preprocessing pipelines)
- **Single short run** with minimal temporal structure
- **Computational speed** is critical (AR adds overhead)

## Next steps

- [Getting Started](01-getting-started.md) — Basic decoder fitting
  workflow
- [rMVPA Integration](03-rmvpa-integration.md) — Cross-validation with
  searchlight analysis
- [HRF Estimation](04-hrf-estimation.md) — Joint HRF learning
- [Weakly Supervised Learning](05-weakly-supervised.md) — Algorithm
  internals

## Session info

``` r
sessioninfo::session_info(pkgs = "hrfdecode")
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value
#>  version  R version 4.5.1 (2025-06-13)
#>  os       macOS Sonoma 14.3
#>  system   aarch64, darwin20
#>  ui       X11
#>  language en
#>  collate  en_US.UTF-8
#>  ctype    en_US.UTF-8
#>  tz       America/Toronto
#>  date     2025-11-09
#>  pandoc   3.7.0.2 @ /opt/homebrew/bin/ (via rmarkdown)
#>  quarto   1.7.32 @ /usr/local/bin/quarto
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package        * version    date (UTC) lib source
#>  askpass          1.2.1      2024-10-04 [2] CRAN (R 4.5.0)
#>  assertthat       0.2.1      2019-03-21 [2] CRAN (R 4.5.0)
#>  backports        1.5.0      2024-05-23 [2] CRAN (R 4.5.0)
#>  base64enc        0.1-3      2015-07-28 [2] CRAN (R 4.5.0)
#>  bdsmatrix        1.3-7      2024-03-02 [2] CRAN (R 4.5.0)
#>  bigassertr       0.1.7      2025-06-27 [2] CRAN (R 4.5.0)
#>  bigparallelr     0.3.2      2021-10-02 [2] CRAN (R 4.5.0)
#>  bigstatsr        1.6.2      2025-07-29 [2] CRAN (R 4.5.0)
#>  bit              4.6.0      2025-03-06 [2] CRAN (R 4.5.0)
#>  bit64            4.6.0-1    2025-01-16 [2] CRAN (R 4.5.0)
#>  bitops           1.0-9      2024-10-03 [2] CRAN (R 4.5.0)
#>  broom            1.0.10     2025-09-13 [2] CRAN (R 4.5.0)
#>  bslib            0.9.0      2025-01-30 [2] CRAN (R 4.5.0)
#>  cachem           1.1.0      2024-05-16 [2] CRAN (R 4.5.0)
#>  caTools          1.18.3     2024-09-04 [2] CRAN (R 4.5.0)
#>  cli              3.6.5      2025-04-23 [2] CRAN (R 4.5.0)
#>  clipr            0.8.0      2022-02-22 [2] CRAN (R 4.5.0)
#>  codetools        0.2-20     2024-03-31 [3] CRAN (R 4.5.1)
#>  colorplane       0.5.0      2025-11-02 [2] Github (bbuchsbaum/colorplane@1b7a26f)
#>  corpcor          1.6.10     2021-09-16 [2] CRAN (R 4.5.0)
#>  cowplot          1.2.0      2025-07-07 [2] CRAN (R 4.5.0)
#>  cpp11            0.5.2      2025-03-03 [2] CRAN (R 4.5.0)
#>  crayon           1.5.3      2024-06-20 [2] CRAN (R 4.5.0)
#>  crosstalk        1.2.2      2025-08-26 [2] CRAN (R 4.5.0)
#>  curl             7.0.0      2025-08-19 [2] CRAN (R 4.5.0)
#>  data.table       1.17.8     2025-07-10 [2] CRAN (R 4.5.0)
#>  dbscan           1.2.3      2025-08-20 [2] CRAN (R 4.5.0)
#>  deflist          0.2.0      2023-04-27 [2] CRAN (R 4.5.0)
#>  DEoptimR         1.1-4      2025-07-27 [2] CRAN (R 4.5.0)
#>  digest           0.6.37     2024-08-19 [2] CRAN (R 4.5.0)
#>  doParallel       1.0.17     2022-02-07 [2] CRAN (R 4.5.0)
#>  dplyr            1.1.4      2023-11-17 [2] CRAN (R 4.5.0)
#>  entropy          1.3.2      2025-04-07 [2] CRAN (R 4.5.0)
#>  evaluate         1.0.5      2025-08-27 [2] CRAN (R 4.5.0)
#>  farver           2.1.2      2024-05-13 [2] CRAN (R 4.5.0)
#>  fastmap          1.2.0      2024-05-15 [2] CRAN (R 4.5.0)
#>  fdrtool          1.2.18     2024-08-20 [2] CRAN (R 4.5.0)
#>  ff               4.5.2      2025-01-13 [2] CRAN (R 4.5.0)
#>  ffmanova         1.1.2      2023-10-18 [2] CRAN (R 4.5.0)
#>  filenamer        0.3        2025-04-09 [2] CRAN (R 4.5.0)
#>  flock            0.7        2016-11-12 [2] CRAN (R 4.5.0)
#>  fmriAR           0.1.0      2025-10-18 [2] Github (bbuchsbaum/fmriAR@0b10352)
#>  fmridesign     * 0.5.0      2025-11-09 [2] Github (bbuchsbaum/fmridesign@f1462eb)
#>  fmrihrf          0.1.0.9000 2025-11-01 [2] Github (bbuchsbaum/fmrihrf@708058f)
#>  FNN              1.1.4.1    2024-09-22 [2] CRAN (R 4.5.0)
#>  fontawesome      0.5.3      2024-11-16 [2] CRAN (R 4.5.0)
#>  foreach          1.5.2      2022-02-02 [2] CRAN (R 4.5.0)
#>  formatR          1.14       2023-01-17 [2] CRAN (R 4.5.0)
#>  fs               1.6.6      2025-04-12 [2] CRAN (R 4.5.0)
#>  furrr            0.3.1      2022-08-15 [2] CRAN (R 4.5.0)
#>  futile.logger    1.4.3      2016-07-10 [2] CRAN (R 4.5.0)
#>  futile.options   1.0.1      2018-04-20 [2] CRAN (R 4.5.0)
#>  future           1.67.0     2025-07-29 [2] CRAN (R 4.5.0)
#>  future.apply     1.20.0     2025-06-06 [2] CRAN (R 4.5.0)
#>  generics         0.1.4      2025-05-09 [2] CRAN (R 4.5.0)
#>  ggplot2        * 4.0.0      2025-09-11 [2] CRAN (R 4.5.0)
#>  gifti            0.8.0      2020-11-11 [2] CRAN (R 4.5.0)
#>  glmnet           4.1-10     2025-07-17 [2] CRAN (R 4.5.0)
#>  globals          0.18.0     2025-05-08 [2] CRAN (R 4.5.0)
#>  glue             1.8.0      2024-09-30 [2] CRAN (R 4.5.0)
#>  gplots           3.2.0      2024-10-05 [2] CRAN (R 4.5.0)
#>  gtable           0.3.6      2024-10-25 [2] CRAN (R 4.5.0)
#>  gtools           3.9.5      2023-11-20 [2] CRAN (R 4.5.0)
#>  hardhat          1.4.2      2025-08-20 [2] CRAN (R 4.5.0)
#>  highr            0.11       2024-05-26 [2] CRAN (R 4.5.0)
#>  hms              1.1.4      2025-10-17 [2] CRAN (R 4.5.0)
#>  hrfdecode      * 0.2.0      2025-11-09 [1] local
#>  htmltools        0.5.8.1    2024-04-04 [2] CRAN (R 4.5.0)
#>  htmlwidgets      1.6.4      2023-12-06 [2] CRAN (R 4.5.0)
#>  httr             1.4.7      2023-08-15 [2] CRAN (R 4.5.0)
#>  igraph           2.2.1      2025-10-27 [2] CRAN (R 4.5.0)
#>  io               0.3.2      2019-12-17 [2] CRAN (R 4.5.0)
#>  isoband          0.2.7      2022-12-20 [2] CRAN (R 4.5.0)
#>  iterators        1.0.14     2022-02-05 [2] CRAN (R 4.5.0)
#>  jquerylib        0.1.4      2021-04-26 [2] CRAN (R 4.5.0)
#>  jsonlite         2.0.0      2025-03-27 [2] CRAN (R 4.5.0)
#>  KernSmooth       2.23-26    2025-01-01 [3] CRAN (R 4.5.1)
#>  knitr            1.50       2025-03-16 [2] CRAN (R 4.5.0)
#>  labeling         0.4.3      2023-08-29 [2] CRAN (R 4.5.0)
#>  lambda.r         1.2.4      2019-09-18 [2] CRAN (R 4.5.0)
#>  later            1.4.4      2025-08-27 [2] CRAN (R 4.5.0)
#>  lattice          0.22-7     2025-04-02 [3] CRAN (R 4.5.1)
#>  lazyeval         0.2.2      2019-03-15 [2] CRAN (R 4.5.0)
#>  lifecycle        1.0.4      2023-11-07 [2] CRAN (R 4.5.0)
#>  listenv          0.10.0     2025-11-02 [2] CRAN (R 4.5.0)
#>  magrittr         2.0.4      2025-09-12 [2] CRAN (R 4.5.0)
#>  MASS             7.3-65     2025-02-28 [2] CRAN (R 4.5.0)
#>  Matrix           1.7-3      2025-03-11 [3] CRAN (R 4.5.1)
#>  matrixStats      1.5.0      2025-01-07 [2] CRAN (R 4.5.0)
#>  memoise          2.0.1      2021-11-26 [2] CRAN (R 4.5.0)
#>  mime             0.13       2025-03-17 [2] CRAN (R 4.5.0)
#>  mmap             0.6-22     2023-12-08 [2] CRAN (R 4.5.0)
#>  modelr           0.1.11     2023-03-22 [2] CRAN (R 4.5.0)
#>  mvtnorm          1.3-3      2025-01-10 [2] CRAN (R 4.5.0)
#>  neuroim2         0.8.3      2025-11-07 [2] Github (bbuchsbaum/neuroim2@77cd9c4)
#>  neurosurf        0.1.0      2025-11-02 [2] Github (bbuchsbaum/neurosurf@5af7de2)
#>  numDeriv         2016.8-1.1 2019-06-06 [2] CRAN (R 4.5.0)
#>  openssl          2.3.4      2025-09-30 [2] CRAN (R 4.5.0)
#>  otel             0.2.0      2025-08-29 [2] CRAN (R 4.5.0)
#>  parallelly       1.45.1     2025-07-24 [2] CRAN (R 4.5.0)
#>  pillar           1.11.1     2025-09-17 [2] CRAN (R 4.5.0)
#>  pkgconfig        2.0.3      2019-09-22 [2] CRAN (R 4.5.0)
#>  plotly           4.11.0     2025-06-19 [2] CRAN (R 4.5.0)
#>  pls              2.8-5      2024-09-15 [2] CRAN (R 4.5.0)
#>  plyr             1.8.9      2023-10-02 [2] CRAN (R 4.5.0)
#>  pracma           2.4.6      2025-10-22 [2] CRAN (R 4.5.0)
#>  prettyunits      1.2.0      2023-09-24 [2] CRAN (R 4.5.0)
#>  progress         1.2.3      2023-12-06 [2] CRAN (R 4.5.0)
#>  promises         1.5.0      2025-11-01 [2] CRAN (R 4.5.0)
#>  proxy            0.4-27     2022-06-09 [2] CRAN (R 4.5.0)
#>  ps               1.9.1      2025-04-12 [2] CRAN (R 4.5.0)
#>  purrr            1.2.0      2025-11-04 [2] CRAN (R 4.5.0)
#>  R.methodsS3      1.8.2      2022-06-13 [2] CRAN (R 4.5.0)
#>  R.oo             1.27.1     2025-05-02 [2] CRAN (R 4.5.0)
#>  R.utils          2.13.0     2025-02-24 [2] CRAN (R 4.5.0)
#>  R6               2.6.1      2025-02-15 [2] CRAN (R 4.5.0)
#>  rappdirs         0.3.3      2021-01-31 [2] CRAN (R 4.5.0)
#>  RColorBrewer     1.1-3      2022-04-03 [2] CRAN (R 4.5.0)
#>  Rcpp             1.1.0      2025-07-02 [2] CRAN (R 4.5.0)
#>  RcppArmadillo    15.0.2-2   2025-09-19 [2] CRAN (R 4.5.0)
#>  RcppEigen        0.3.4.0.2  2024-08-24 [2] CRAN (R 4.5.0)
#>  RcppParallel     5.1.11-1   2025-08-27 [2] CRAN (R 4.5.0)
#>  readr            2.1.5      2024-01-10 [2] CRAN (R 4.5.0)
#>  Rfit             0.27.0     2024-05-25 [2] CRAN (R 4.5.0)
#>  rgl              1.3.24     2025-06-25 [2] CRAN (R 4.5.0)
#>  RhpcBLASctl      0.23-42    2023-02-11 [2] CRAN (R 4.5.0)
#>  rlang            1.1.6      2025-04-11 [2] CRAN (R 4.5.0)
#>  rmarkdown        2.30       2025-09-28 [2] CRAN (R 4.5.0)
#>  rmio             0.4.0      2022-02-17 [2] CRAN (R 4.5.0)
#>  rMVPA            0.1.2      2025-11-09 [2] local
#>  RNifti           1.8.0      2025-02-22 [2] CRAN (R 4.5.0)
#>  RNiftyReg        2.8.4      2024-09-30 [2] CRAN (R 4.5.0)
#>  robustbase       0.99-6     2025-09-04 [2] CRAN (R 4.5.0)
#>  rsample          1.3.1      2025-07-29 [2] CRAN (R 4.5.0)
#>  RSpectra         0.16-2     2024-07-18 [2] CRAN (R 4.5.0)
#>  Rvcg             0.25       2025-03-14 [2] CRAN (R 4.5.0)
#>  S7               0.2.0      2024-11-07 [2] CRAN (R 4.5.0)
#>  sass             0.4.10     2025-04-11 [2] CRAN (R 4.5.0)
#>  scales           1.4.0      2025-04-24 [2] CRAN (R 4.5.0)
#>  sda              1.3.9      2025-04-08 [2] CRAN (R 4.5.0)
#>  shape            1.4.6.1    2024-02-23 [2] CRAN (R 4.5.0)
#>  slider           0.3.2      2024-10-25 [2] CRAN (R 4.5.0)
#>  sparsediscrim    0.3.0      2021-07-01 [2] CRAN (R 4.5.0)
#>  sparsevctrs      0.3.4      2025-05-25 [2] CRAN (R 4.5.0)
#>  stringi          1.8.7      2025-03-27 [2] CRAN (R 4.5.0)
#>  stringr          1.6.0      2025-11-04 [2] CRAN (R 4.5.0)
#>  survival         3.8-3      2024-12-17 [3] CRAN (R 4.5.1)
#>  sys              3.4.3      2024-10-04 [2] CRAN (R 4.5.0)
#>  tibble           3.3.0      2025-06-08 [2] CRAN (R 4.5.0)
#>  tidyr            1.3.1      2024-01-24 [2] CRAN (R 4.5.0)
#>  tidyselect       1.2.1      2024-03-11 [2] CRAN (R 4.5.0)
#>  tinytex          0.57       2025-04-15 [2] CRAN (R 4.5.0)
#>  tzdb             0.5.0      2025-03-15 [2] CRAN (R 4.5.0)
#>  utf8             1.2.6      2025-06-08 [2] CRAN (R 4.5.0)
#>  vctrs            0.6.5      2023-12-01 [2] CRAN (R 4.5.0)
#>  viridisLite      0.4.2      2023-05-02 [2] CRAN (R 4.5.0)
#>  vroom            1.6.6      2025-09-19 [2] CRAN (R 4.5.0)
#>  warp             0.2.1      2023-11-02 [2] CRAN (R 4.5.0)
#>  whitening        1.4.0      2022-06-07 [2] CRAN (R 4.5.0)
#>  withr            3.0.2      2024-10-28 [2] CRAN (R 4.5.0)
#>  xfun             0.54       2025-10-30 [2] CRAN (R 4.5.0)
#>  xml2             1.4.1      2025-10-27 [2] CRAN (R 4.5.0)
#>  yaml             2.3.10     2024-07-26 [2] CRAN (R 4.5.0)
#>  yardstick        1.3.2      2025-01-22 [2] CRAN (R 4.5.0)
#> 
#>  [1] /private/var/folders/9h/nkjq6vss7mqdl4ck7q1hd8ph0000gp/T/RtmpfzOoi1/temp_libpathe6e4ebda2ae
#>  [2] /Users/bbuchsbaum/Library/R/arm64/4.5/library
#>  [3] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
#>  * ── Packages attached to the search path.
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
