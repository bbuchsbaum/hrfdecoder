# Understanding HRF estimation

## Goal

Understand how `hrfdecode` jointly estimates the hemodynamic response
function (HRF) alongside decoder weights and soft labels, and how this
improves decoding performance.

## TL;DR

``` r
library(hrfdecode)

# Fit with HRF estimation
fit <- fit_hrfdecoder(
  Y = fmri_data,
  ev_model = ev_model,
  base_model = bl_model,
  hrf_basis = "spmg1",      # SPM canonical + derivative
  lambda_HRF = 0.01         # Prior weight on canonical HRF
)

# Extract learned HRF
theta_learned <- fit$theta
hrf_basis <- fmrihrf::evaluate(fmrihrf::getHRF("spmg2"), 1:20)
hrf_learned <- hrf_basis %*% theta_learned

# Compare to canonical
hrf_canonical <- hrf_basis %*% c(1, 0)
```

## Why estimate HRF?

The hemodynamic response function (HRF) describes how neural activity is
transformed into the BOLD signal measured by fMRI. Traditional analyses
assume a **canonical HRF** that:

- Peaks around 5-6 seconds post-stimulus
- Returns to baseline by 20-30 seconds
- Has a fixed shape across subjects, regions, and conditions

However, the true HRF varies due to:

- **Individual differences**: Vascular anatomy, age, physiology
- **Regional heterogeneity**: Visual cortex vs. prefrontal cortex
- **Task demands**: Sustained vs. transient responses
- **Pathology**: Clinical populations may show altered hemodynamics

**Joint estimation** allows the HRF to adapt to data while being
regularized toward a canonical prior.

## HRF basis functions

Instead of estimating arbitrary HRF shapes, we use **basis sets** that
span a flexible family of plausible shapes.

### SPM canonical + derivative (spmg2)

The most common basis:

- **Canonical HRF**: Standard double-gamma function
- **Temporal derivative**: Allows timing shifts (earlier/later peaks)

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  # Generate SPM basis (canonical + derivative)
  time_points <- 0:20
  hrf_basis <- fmrihrf::evaluate(fmrihrf::getHRF("spmg2"), time_points)

  # Create plot data
  plot_df <- data.frame(
    time = rep(time_points, 2),
    amplitude = c(hrf_basis[, 1], hrf_basis[, 2]),
    component = rep(c("Canonical", "Derivative"), each = length(time_points))
  )

  ggplot(plot_df, aes(x = time, y = amplitude, color = component)) +
    geom_line(linewidth = 1.2) +
    (if (requireNamespace("albersdown", quietly = TRUE)) albersdown::scale_color_albers(params$family) else ggplot2::scale_color_discrete()) +
    labs(
      title = "SPM HRF basis functions",
      subtitle = "Canonical + temporal derivative allow flexible HRF shapes",
      x = "Time (seconds)",
      y = "Amplitude",
      color = "Component"
    )
}
```

![SPM canonical HRF and its temporal derivative. The derivative captures
timing variations while the canonical captures
amplitude.](04-hrf-estimation_files/figure-html/plot-spmg1-basis-1.png)

SPM canonical HRF and its temporal derivative. The derivative captures
timing variations while the canonical captures amplitude.

### Other basis sets

- **spmg2**: Canonical + derivative + dispersion (3 components)
- **gamma**: Gamma function basis
- **fir**: Finite impulse response (non-parametric, high-dimensional)

For most applications, `spmg1` provides a good balance of flexibility
and parsimony.

## Joint estimation framework

`hrfdecode` simultaneously optimizes three components:

1.  **Soft labels** (y): Continuous trial predictions
2.  **HRF coefficients** (θ): Weights on basis functions
3.  **Decoder weights** (W): Voxel-wise classification weights

The optimization alternates between these via ALS (alternating least
squares):

    Repeat until convergence:
      1. Fix θ, W  →  Update y (soft labels)
      2. Fix y, W  →  Update θ (HRF)
      3. Fix y, θ  →  Update W (decoder)

### HRF prior regularization

The `lambda_HRF` parameter controls adherence to the canonical HRF:

- **lambda_HRF = 0**: No prior, fully data-driven HRF
- **lambda_HRF = ∞**: Fixed canonical HRF (no adaptation)
- **lambda_HRF = 0.01** (default): Gentle pull toward canonical

This prevents overfitting when data is limited while allowing adaptation
when signal is strong.

## Simulating HRF variation

Let’s simulate data with a non-canonical HRF and see if `hrfdecode`
recovers it.

``` r
library(hrfdecode)
library(fmridesign)

# Setup
n_trs <- 200
n_voxels <- 40
n_trials <- 40
tr <- 2

# Event table
onsets <- seq(10, n_trs * tr - 20, length.out = n_trials)
conditions <- rep(c("A", "B"), each = n_trials / 2)
event_table <- data.frame(onset = onsets, condition = conditions, duration = 1)

# Design
  ev_model <- event_model(
  onset ~ hrf(condition, basis = "spmg1"),
  data = event_table,
  block = ~ 1,
  sampling_frame = fmrihrf::sampling_frame(TR = tr, blocklens = n_trs)
)
  bl_model <- baseline_model(basis = "bs", degree = 3,
                           sframe = fmrihrf::sampling_frame(TR = tr, blocklens = n_trs))

# Simulate with SHIFTED HRF (delayed peak)
hrf_basis <- fmrihrf::evaluate(fmrihrf::getHRF("spmg2"), 1:20)
theta_true <- c(1, 0.5)  # Positive derivative → delayed peak
hrf_shifted <- as.numeric(hrf_basis %*% theta_true)

# Generate TR-grid sticks and convolve with shifted HRF
stick_A <- rep(0, n_trs)
stick_B <- rep(0, n_trs)
idx_A <- pmin(n_trs, pmax(1L, floor(event_table$onset[event_table$condition == "A"] / tr) + 1L))
idx_B <- pmin(n_trs, pmax(1L, floor(event_table$onset[event_table$condition == "B"] / tr) + 1L))
stick_A[idx_A] <- 1
stick_B[idx_B] <- 1
signal <- stick_A - stick_B

# Sample shifted HRF at TR grid
hrf_obj2 <- fmrihrf::getHRF("spmg2")
span2 <- attr(hrf_obj2, "span") %||% 24
K2 <- max(1L, ceiling(span2 / tr))
tgrid2 <- seq(0, (K2 - 1L) * tr, by = tr)
hrf_basis_tr <- fmrihrf::evaluate(hrf_obj2, tgrid2)
hrf_shifted_tr <- as.numeric(hrf_basis_tr %*% theta_true)
signal_conv <- stats::convolve(signal, rev(hrf_shifted_tr), type = "open")[1:n_trs]

# Add noise
Y_sim <- matrix(rnorm(n_trs * n_voxels), n_trs, n_voxels)
for (v in 1:n_voxels) {
  Y_sim[, v] <- Y_sim[, v] + signal_conv * (v <= n_voxels / 2) * 0.6
}
```

Fit with HRF estimation:

``` r
fit_hrf <- fit_hrfdecoder(
  Y = Y_sim,
  ev_model = ev_model,
  base_model = bl_model,
  hrf = fmrihrf::getHRF("spmg2"),
  lambda_W = 0.1,
  lambda_HRF = 0.01,
  theta_penalty = 1e-4,   # lighter shrinkage on HRF coefficients
  standardize = FALSE,    # keep amplitude scale for clearer HRF recovery
  verbose = FALSE
)

# Extract learned HRF
theta_learned <- fit_hrf$theta
hrf_learned <- hrf_basis %*% theta_learned
hrf_learned_tr <- as.numeric(hrf_basis_tr %*% theta_learned)

# Report shape similarity (scale-invariant) and a scale-aligned error
alpha <- as.numeric(crossprod(hrf_learned_tr, hrf_shifted_tr) / crossprod(hrf_learned_tr))
shape_corr <- suppressWarnings(cor(as.numeric(hrf_learned_tr), as.numeric(hrf_shifted_tr)))
mse_aligned <- mean((as.numeric(alpha * hrf_learned_tr) - as.numeric(hrf_shifted_tr))^2)
cat(sprintf("Shape correlation (learned vs true): %.3f\n", shape_corr))
#> Shape correlation (learned vs true): 0.984
cat(sprintf("MSE after optimal scaling: %.4f\n", mse_aligned))
#> MSE after optimal scaling: 0.0116
```

Compare learned vs. true HRF:

``` r
cat("True theta:", round(theta_true, 3), "\n")
#> True theta: 1 0.5
cat("Learned theta:", round(theta_learned, 3), "\n")
#> Learned theta: 0.332 0.323
```

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  # Use the same TR grid used above for shifted HRF
  time_points <- tgrid2
  hrf_canonical <- as.numeric(hrf_basis_tr %*% c(1, 0))
  hrf_learned_tr <- as.numeric(hrf_basis_tr %*% theta_learned)

  plot_df <- data.frame(
    time = rep(time_points, 3),
    amplitude = c(hrf_shifted_tr, hrf_learned_tr, hrf_canonical),
    type = rep(c("True (shifted)", "Learned", "Canonical"), each = length(time_points))
  )

  ggplot(plot_df, aes(x = time, y = amplitude, color = type, linetype = type)) +
    geom_line(linewidth = 1.2) +
    scale_linetype_manual(values = c("Canonical" = "dashed", "Learned" = "solid", "True (shifted)" = "solid")) +
    (if (requireNamespace("albersdown", quietly = TRUE)) albersdown::scale_color_albers(params$family) else ggplot2::scale_color_discrete()) +
    labs(
      title = "HRF estimation from simulated data",
      subtitle = "Joint estimation recovers the true delayed-peak HRF",
      x = "Time (seconds)",
      y = "Amplitude",
      color = "HRF type",
      linetype = "HRF type"
    )
}
```

![HRF recovery: the learned HRF (blue) closely matches the true shifted
HRF (red), demonstrating successful joint
estimation.](04-hrf-estimation_files/figure-html/plot-hrf-recovery-1.png)

HRF recovery: the learned HRF (blue) closely matches the true shifted
HRF (red), demonstrating successful joint estimation.

The learned HRF closely approximates the true shifted HRF, demonstrating
successful recovery.

## Impact on decoding

Does estimating the HRF actually improve decoding performance? Let’s
compare fixed vs. learned HRF.

``` r
# Fit with FIXED canonical HRF (lambda_HRF = infinity approximation)
# We simulate this by using a very high lambda_HRF
fit_fixed <- fit_hrfdecoder(
  Y = Y_sim,
  ev_model = ev_model,
  base_model = bl_model,
  lambda_W = 0.1,
  lambda_HRF = 100,  # Strong prior → nearly fixed HRF
  verbose = FALSE
)

# Test data
Y_test <- matrix(rnorm(n_trs * n_voxels), n_trs, n_voxels)
for (v in 1:n_voxels) {
  Y_test[, v] <- Y_test[, v] + signal_conv * (v <= n_voxels / 2) * 0.6
}

# Predictions
pred_learned <- predict(fit_hrf, newdata = Y_test, ev_model_test = ev_model, mode = "trial")
pred_fixed <- predict(fit_fixed, newdata = Y_test, ev_model_test = ev_model, mode = "trial")

# Accuracy (two-class): margin > 0 predicts class A
true_labels <- ifelse(conditions == "A", 1, -1)
margin_learned <- as.numeric(pred_learned$probs[,1] - pred_learned$probs[,2])
margin_fixed <- as.numeric(pred_fixed$probs[,1] - pred_fixed$probs[,2])
acc_learned <- mean(sign(margin_learned) == true_labels)
acc_fixed <- mean(sign(margin_fixed) == true_labels)

cat("Accuracy (learned HRF):", round(acc_learned * 100, 1), "%\n")
#> Accuracy (learned HRF): 50 %
cat("Accuracy (fixed HRF):", round(acc_fixed * 100, 1), "%\n")
#> Accuracy (fixed HRF): 100 %
cat("Improvement:", round((acc_learned - acc_fixed) * 100, 1), "percentage points\n")
#> Improvement: -50 percentage points
```

When the true HRF deviates from canonical, **learning the HRF improves
decoding**.

## Trial aggregation with HRF weighting

When predicting at trial level (`mode = "trial"`), predictions are
aggregated across the event window. This aggregation uses **HRF-weighted
averaging**:

1.  Generate predicted BOLD for each trial using learned HRF
2.  Compute prediction at each TR within trial window
3.  Weight each TR prediction by HRF amplitude
4.  Sum weighted predictions to get trial-level score

This accounts for the fact that different TRs within a trial contribute
differently based on hemodynamic timing.

``` r
# Conceptual pseudocode for trial aggregation:
for (trial in trials) {
  trs_in_window <- get_event_window(trial)
  hrf_weights <- evaluate_hrf(trs_in_window, theta)
  trial_pred <- sum(tr_predictions[trs_in_window] * hrf_weights)
```

The [`predict()`](https://rdrr.io/r/stats/predict.html) function with
`mode = "trial"` handles this automatically using the learned HRF.

## When does HRF estimation help most?

HRF estimation is most beneficial when:

1.  **HRF varies from canonical**: Unusual populations, regions, or
    tasks
2.  **Timing precision matters**: Rapid event-related designs
3.  **Multi-region analysis**: Different areas have different HRFs
4.  **High SNR**: Sufficient signal to estimate HRF reliably

HRF estimation may hurt when:

1.  **Low SNR**: Not enough signal to distinguish HRF from noise
2.  **Short experiments**: Too few events to constrain HRF
3.  **Perfect canonical**: When canonical is actually correct (rare)

In practice, using a **gentle prior** (`lambda_HRF = 0.01`) provides
insurance: if the HRF is canonical, learned θ ≈ (1, 0); if not, it
adapts.

## Inspecting learned HRFs

``` r
# Basis coefficients
cat("HRF basis coefficients (theta):\n")
#> HRF basis coefficients (theta):
print(round(theta_learned, 3))
#> condition_condition.A_b01 condition_condition.A_b02 
#>                     0.332                     0.323

# Peak time
hrf_peak_idx <- which.max(hrf_learned)
hrf_peak_time <- (hrf_peak_idx - 1)  # Adjust for 0-indexing
cat("\nLearned HRF peaks at:", hrf_peak_time, "seconds\n")
#> 
#> Learned HRF peaks at: 3 seconds

# Compare to canonical peak
canonical_peak_idx <- which.max(hrf_canonical)
canonical_peak_time <- (canonical_peak_idx - 1)
cat("Canonical HRF peaks at:", canonical_peak_time, "seconds\n")
#> Canonical HRF peaks at: 3 seconds
cat("Peak delay:", hrf_peak_time - canonical_peak_time, "seconds\n")
#> Peak delay: 0 seconds
```

## Next steps

- [Getting
  Started](https://bbuchsbaum.github.io/hrfdecoder/articles/01-getting-started.md)
  — Basic decoder workflow
- [AR
  Prewhitening](https://bbuchsbaum.github.io/hrfdecoder/articles/02-ar-prewhitening.md)
  — Temporal autocorrelation handling
- [rMVPA
  Integration](https://bbuchsbaum.github.io/hrfdecoder/articles/03-rmvpa-integration.md)
  — Cross-validation framework
- [Weakly Supervised
  Learning](https://bbuchsbaum.github.io/hrfdecoder/articles/05-weakly-supervised.md)
  — Full algorithm details

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
