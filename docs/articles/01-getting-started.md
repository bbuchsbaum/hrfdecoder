# Getting started with hrfdecoder

## Goal

Learn how to fit an HRF-aware weakly supervised decoder to fMRI time
series data and make predictions on new data.

## TL;DR

``` r
library(hrfdecode)
library(fmridesign)

# Fit decoder from fMRI timeseries and event onsets
fit <- fit_hrfdecoder(
  Y = fmri_data,           # T × V matrix (time × voxels)
  event_model = ev_model,  # event_model from fmridesign
  baseline_model = bl_model
)

# Predict on new data
preds <- predict(fit, newdata = test_data, mode = "trial")
```

## Introduction

Traditional MVPA requires trial-averaged data or explicit trial labels.
The `hrfdecode` package takes a different approach: it learns directly
from continuous fMRI time series and event timing information. This
“weakly supervised” approach jointly estimates:

1.  **Soft labels** — continuous trial-level predictions
2.  **HRF parameters** — subject-specific hemodynamic response
3.  **Decoder weights** — multivariate classification function

All three components are optimized together using alternating least
squares.

## Setup

``` r
library(hrfdecode)
library(fmridesign)
```

For this tutorial, we’ll create synthetic data that mimics a simple
two-class decoding experiment.

``` r
# Simulation parameters
n_trs <- 200        # Number of time points
n_voxels <- 50      # Number of voxels
n_trials <- 40      # Number of trials
tr <- 2             # TR in seconds

# Create event table (onset times and conditions)
onsets <- seq(10, n_trs * tr - 20, length.out = n_trials)
conditions <- rep(c("A", "B"), each = n_trials / 2)

event_table <- data.frame(
  onset = onsets,
  condition = conditions,
  duration = 1
)

head(event_table)
#>      onset condition duration
#> 1 10.00000         A        1
#> 2 19.48718         A        1
#> 3 28.97436         A        1
#> 4 38.46154         A        1
#> 5 47.94872         A        1
#> 6 57.43590         A        1
```

Next, we define the experimental design using the `fmridesign` package.

``` r
# Create event model from event table
  ev_model <- event_model(
  onset ~ hrf(condition, basis = "spmg1"),
  data = event_table,
  block = ~ 1,
  sampling_frame = fmrihrf::sampling_frame(TR = tr, blocklens = n_trs)
)

# Create baseline model (intercept + linear drift)
  bl_model <- baseline_model(basis = "bs", degree = 3, sframe = fmrihrf::sampling_frame(TR = tr, blocklens = n_trs))
```

Now simulate fMRI data with signal in a subset of voxels.

``` r
# Generate synthetic BOLD data
# First half of voxels discriminate condition A vs B
true_pattern <- c(rep(1, n_voxels / 2), rep(0, n_voxels / 2))

# Build TR-grid stick functions per condition
stick_A <- rep(0, n_trs)
stick_B <- rep(0, n_trs)
idx_A <- pmin(n_trs, pmax(1L, floor(event_table$onset[event_table$condition == "A"] / tr) + 1L))
idx_B <- pmin(n_trs, pmax(1L, floor(event_table$onset[event_table$condition == "B"] / tr) + 1L))
stick_A[idx_A] <- 1
stick_B[idx_B] <- 1
signal <- stick_A - stick_B

# Convolve with HRF (SPMG2 basis) sampled at TR grid
hrf_obj <- fmrihrf::getHRF("spmg2")
span <- attr(hrf_obj, "span") %||% 24
K <- max(1L, ceiling(span / tr))
tgrid <- seq(0, (K - 1L) * tr, by = tr)
hrf_basis <- fmrihrf::evaluate(hrf_obj, tgrid)
hrf_vec <- as.numeric(hrf_basis %*% c(1, 0))
signal_conv <- stats::convolve(signal, rev(hrf_vec), type = "open")[1:n_trs]

# Add to random noise
Y_train <- matrix(rnorm(n_trs * n_voxels, sd = 1), n_trs, n_voxels)
for (v in 1:n_voxels) {
  Y_train[, v] <- Y_train[, v] + signal_conv * true_pattern[v] * 0.5
}

dim(Y_train)
#> [1] 200  50
```

## Fitting the decoder

Now we fit the decoder using
[`fit_hrfdecoder()`](https://bbuchsbaum.github.io/hrfdecoder/reference/fit_hrfdecoder.md).
The key arguments are:

- `Y`: The fMRI data matrix (time × voxels)
- `event_model`: Event design from `fmridesign`
- `baseline_model`: Baseline/nuisance model
- `lambda_W`: Ridge penalty on decoder weights (default: 0.1)
- `max_iter`: Maximum ALS iterations (default: 10)

``` r
fit <- fit_hrfdecoder(
  Y = Y_train,
  ev_model = ev_model,
  base_model = bl_model,
  lambda_W = 0.1,
  max_iter = 10,
  verbose = FALSE
)

# Inspect fit object
class(fit)
#> [1] "hrfdecoder_fit"
names(fit)
#>  [1] "W"           "P"           "b"           "theta"       "hrf"        
#>  [6] "conditions"  "background"  "converged"   "iterations"  "settings"   
#> [11] "train"       "preproc"     "diagnostics"
```

The fitted object contains:

- `W`: Decoder weight vector (V × 1)
- `theta`: HRF basis coefficients
- `y_soft`: Soft labels (continuous trial predictions)
- `preproc_params`: Preprocessing metadata (centering, scaling, AR
  parameters)
- `convergence`: Convergence diagnostics

## Inspecting decoder weights

Let’s examine which voxels contribute most to the decoder.

``` r
# Get top 10 voxels by absolute weight
top_voxels <- order(abs(fit$W), decreasing = TRUE)[1:10]
top_weights <- fit$W[top_voxels]

data.frame(
  voxel = top_voxels,
  weight = round(top_weights, 3),
  true_pattern = true_pattern[top_voxels]
)
#>    voxel weight true_pattern
#> 1     22  0.028            1
#> 2      6  0.026            1
#> 3     15  0.025            1
#> 4     52 -0.023           NA
#> 5     65 -0.023           NA
#> 6     51 -0.022           NA
#> 7      1  0.022            1
#> 8     69 -0.022           NA
#> 9     62 -0.021           NA
#> 10    81 -0.021           NA
```

The decoder correctly identifies signal-bearing voxels (those with
`true_pattern = 1`).

## Making predictions

The [`predict()`](https://rdrr.io/r/stats/predict.html) method supports
two modes:

1.  **“tr”** — TR-level predictions (one per time point)
2.  **“trial”** — Trial-level predictions (aggregated across event
    duration)

``` r
# Create test data (new noise, same signal)
Y_test <- matrix(rnorm(n_trs * n_voxels, sd = 1), n_trs, n_voxels)
for (v in 1:n_voxels) {
  Y_test[, v] <- Y_test[, v] + signal_conv * true_pattern[v] * 0.5
}

# TR-level predictions
pred_tr <- predict(fit, newdata = Y_test, mode = "tr")
nrow(pred_tr)  # One prediction per TR
#> [1] 200
```

``` r
# Trial-level predictions (aggregated within event windows)
pred_trial <- predict(fit, newdata = Y_test, ev_model_test = ev_model, mode = "trial")
nrow(pred_trial$probs)  # One prediction per trial
#> [1] 40
```

Trial-level predictions aggregate the TR-level signal using HRF-weighted
averaging within each event window.

## Evaluating performance

For a two-class problem, we can compute classification accuracy.

``` r
# Convert event-level probabilities to class predictions
true_labels <- ifelse(conditions == "A", 1, -1)
pred_classes <- ifelse(pred_trial$probs[, 1] >= pred_trial$probs[, 2], 1, -1)

# Accuracy
accuracy <- mean(pred_classes == true_labels)
cat("Classification accuracy:", round(accuracy * 100, 1), "%\n")
#> Classification accuracy: 100 %
```

## Visualizing predictions

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  plot_df <- data.frame(
    trial = 1:n_trials,
    prediction = as.numeric(pred_trial$probs[, 1] - pred_trial$probs[, 2]),
    true_label = true_labels,
    condition = conditions
  )

  scl <- if (requireNamespace("albersdown", quietly = TRUE)) {
    albersdown::scale_color_albers(params$family)
  } else {
    ggplot2::scale_color_discrete()
  }

  ggplot(plot_df, aes(x = trial, y = prediction, color = condition)) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scl +
    labs(
      title = "Trial-level decoder predictions",
      subtitle = "Weakly supervised decoder correctly separates conditions",
      x = "Trial number",
      y = "Soft label prediction",
      color = "Condition"
    )
}
```

![Trial-level predictions vs. true condition labels. Positive values
indicate condition A, negative values condition
B.](01-getting-started_files/figure-html/plot-predictions-1.png)

Trial-level predictions vs. true condition labels. Positive values
indicate condition A, negative values condition B.

## Next steps

This tutorial covered the basics of fitting and predicting with
`hrfdecode`. For more advanced topics:

- [AR
  Prewhitening](https://bbuchsbaum.github.io/hrfdecoder/articles/02-ar-prewhitening.md)
  — Handle temporal autocorrelation in multi-run data
- [rMVPA
  Integration](https://bbuchsbaum.github.io/hrfdecoder/articles/03-rmvpa-integration.md)
  — Run searchlight analysis with cross-validation
- [HRF
  Estimation](https://bbuchsbaum.github.io/hrfdecoder/articles/04-hrf-estimation.md)
  — Understand joint HRF learning and event aggregation
- [Weakly Supervised
  Learning](https://bbuchsbaum.github.io/hrfdecoder/articles/05-weakly-supervised.md)
  — Deep dive into the ALS algorithm

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
#>  date     2025-11-10
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
#>  fmridesign     * 0.5.0      2025-11-10 [2] Github (bbuchsbaum/fmridesign@9af622f)
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
#>  hrfdecode      * 0.2.0      2025-11-10 [1] local
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
#>  rMVPA            0.1.2      2025-11-10 [2] Github (bbuchsbaum/rMVPA@8858610)
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
#>  [1] /private/var/folders/9h/nkjq6vss7mqdl4ck7q1hd8ph0000gp/T/RtmpsnT3u0/temp_libpathad1faeea7f2
#>  [2] /Users/bbuchsbaum/Library/R/arm64/4.5/library
#>  [3] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
#>  * ── Packages attached to the search path.
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
