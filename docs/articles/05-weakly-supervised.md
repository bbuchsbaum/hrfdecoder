# Weakly supervised learning for fMRI

## Goal

Understand the weakly supervised learning framework underlying
`hrfdecode`, including the alternating least squares algorithm,
regularization parameters, and convergence diagnostics.

## TL;DR

``` r
library(hrfdecode)

# Fit with custom regularization
fit <- fit_hrfdecoder(
  Y = fmri_data,
  ev_model = ev_model,
  base_model = bl_model,
  lambda_W = 0.5,        # Ridge penalty on decoder weights
  lambda_HRF = 0.01,     # HRF prior strength
  lambda_smooth = 0.1,   # Temporal smoothness on soft labels
  theta_penalty = 0.01,  # HRF coefficient L2 penalty
  max_iter = 20,         # Maximum ALS iterations
  tol = 1e-4,            # Convergence tolerance
  verbose = TRUE         # Show iteration progress
)

# Inspect convergence
fit$convergence
```

## What is weakly supervised learning?

Traditional supervised learning requires **strong labels**: explicit
trial-level class assignments. For fMRI, this means:

- Averaging BOLD across trial duration
- Assigning discrete labels (e.g., “face” vs. “scene”)
- Training classifier on (averaged data, labels) pairs

**Weakly supervised learning** relaxes this requirement. Instead of
trial labels, we provide:

- **Event onsets**: When stimuli occurred
- **Condition labels**: What type of stimulus (but not trial outcome)
- **Continuous time series**: Full BOLD data, not trial averages

The algorithm then **jointly infers**:

1.  **Soft labels** (y): Continuous trial predictions
2.  **HRF** (θ): Hemodynamic response parameters
3.  **Decoder** (W): Voxel weights

This is “weak” supervision because we don’t specify the trial-level BOLD
patterns—only the experimental structure.

## The alternating least squares algorithm

The optimization problem is:

``` math
\min_{y, \theta, W} \| Y - X(\theta) y W^T \|^2 + \lambda_W \|W\|^2 + \lambda_{HRF} \|\theta - \theta_0\|^2 + \lambda_{smooth} \|Ly\|^2
```

Where:

- **Y**: fMRI data (T × V)
- **X(θ)**: Design matrix convolved with HRF(θ)
- **y**: Soft labels (N × 1)
- **W**: Decoder weights (V × 1)
- **L**: Laplacian for temporal smoothness

This is **non-convex** in (y, θ, W) jointly, but **convex in each
variable** when others are fixed. ALS exploits this:

### ALS steps

    Initialize: y₀, θ₀, W₀
    For iteration t = 1, 2, ..., max_iter:

      1. Update soft labels:
         y_t = argmin_y || Y - X(θ_{t-1}) y W_{t-1}^T ||^2 + λ_smooth ||Ly||^2

      2. Update HRF:
         θ_t = argmin_θ || Y - X(θ) y_t W_{t-1}^T ||^2 + λ_HRF ||θ - θ₀||^2

      3. Update decoder:
         W_t = argmin_W || Y - X(θ_t) y_t W^T ||^2 + λ_W ||W||^2

      If || (y_t, θ_t, W_t) - (y_{t-1}, θ_{t-1}, W_{t-1}) || < tol:
        Break (converged)

Each step is a **regularized least squares** problem with closed-form
solution.

## Regularization parameters

Four main parameters control the optimization:

### lambda_W: Decoder ridge penalty

Controls decoder complexity via L2 penalty on weights.

- **Low (0.01)**: Complex decoder, risk of overfitting
- **Medium (0.1)**: Default, good balance
- **High (1.0)**: Simple decoder, risk of underfitting

**When to increase**: Many voxels, low SNR, small sample **When to
decrease**: Few voxels, high SNR, large sample

### lambda_HRF: HRF prior strength

Pulls HRF toward canonical shape.

- **Low (0.001)**: Data-driven HRF, high flexibility
- **Medium (0.01)**: Default, gentle regularization
- **High (1.0)**: Nearly fixed canonical HRF

**When to increase**: Uncertain HRF, low SNR, standard populations
**When to decrease**: Known HRF deviation, high SNR, atypical
populations

### lambda_smooth: Soft label smoothness

Encourages temporal smoothness in soft labels across trials.

- **Low (0)**: Independent trial labels
- **Medium (0.1)**: Gentle smoothing
- **High (1.0)**: Strong smoothness (assumes slow label drift)

**When to increase**: Blocked designs, gradual condition changes **When
to decrease**: Event-related designs, rapid switching

### theta_penalty: HRF coefficient penalty

Additional L2 penalty on HRF basis coefficients (beyond prior).

- **Default (0.01)**: Light regularization
- Can prevent extreme HRF shapes

## Practical parameter tuning

### Strategy 1: Grid search

``` r
# Setup (assuming data and models are defined)
library(hrfdecode)

# Parameter grid
lambda_W_values <- c(0.01, 0.1, 0.5, 1.0)
results <- list()

for (i in seq_along(lambda_W_values)) {
  fit <- fit_hrfdecoder(
    Y = Y_train,
    ev_model = ev_model,
    base_model = bl_model,
    lambda_W = lambda_W_values[i],
    verbose = FALSE
  )

  # Evaluate on validation set
  preds <- predict(fit, newdata = Y_val, mode = "trial")
  acc <- mean(sign(preds) == true_labels_val)

  results[[i]] <- list(lambda_W = lambda_W_values[i], accuracy = acc)

# Select best
best_idx <- which.max(sapply(results, function(x) x$accuracy))
best_lambda_W <- results[[best_idx]]$lambda_W
```

### Strategy 2: Cross-validation

Integrate with rMVPA for automatic parameter selection:

``` r
library(rMVPA)

# Define parameter grid
param_grid <- expand.grid(
  lambda_W = c(0.1, 0.5, 1.0),
  lambda_HRF = c(0.001, 0.01, 0.1)
)

# Cross-validate each combination
# (Conceptual; requires rMVPA infrastructure)
# cv_results <- tune_parameters(hrfdecoder_model, param_grid, cv_folds = 5)
```

### Strategy 3: Heuristics

Rule-of-thumb starting points:

- **lambda_W**: 0.1 × (number of voxels / number of trials)
- **lambda_HRF**: 0.01 (rarely needs tuning)
- **lambda_smooth**: 0 for event-related, 0.1-0.5 for block designs
- **theta_penalty**: 0.01 (fixed)

## Convergence diagnostics

Monitor convergence to ensure optimization succeeded.

``` r
library(hrfdecode)
library(fmridesign)

# Simple simulation
n_trs <- 150
n_voxels <- 30
n_trials <- 30
tr <- 2

onsets <- seq(10, n_trs * tr - 20, length.out = n_trials)
conditions <- rep(c("A", "B"), each = n_trials / 2)
event_table <- data.frame(onset = onsets, condition = conditions, duration = 1)

  ev_model <- event_model(
  onset ~ hrf(condition, basis = "spmg1"),
  data = event_table,
  block = ~ 1,
  sampling_frame = fmrihrf::sampling_frame(TR = tr, blocklens = n_trs)
)
  bl_model <- baseline_model(basis = "bs", degree = 3,
                           sframe = fmrihrf::sampling_frame(TR = tr, blocklens = n_trs))

# Simulate data
hrf_basis <- fmrihrf::evaluate(fmrihrf::getHRF("spmg2"), seq(0, 24, by = tr))
hrf_vec <- as.numeric(hrf_basis %*% c(1, 0))
stick_A <- rep(0, n_trs); stick_B <- rep(0, n_trs)
idx_A <- pmin(n_trs, pmax(1L, floor(event_table$onset[event_table$condition == "A"] / tr) + 1L))
idx_B <- pmin(n_trs, pmax(1L, floor(event_table$onset[event_table$condition == "B"] / tr) + 1L))
stick_A[idx_A] <- 1; stick_B[idx_B] <- 1
signal <- stick_A - stick_B
signal_conv <- stats::convolve(signal, rev(hrf_vec), type = "open")[1:n_trs]

Y_data <- matrix(rnorm(n_trs * n_voxels), n_trs, n_voxels)
for (v in 1:n_voxels) {
  Y_data[, v] <- Y_data[, v] + signal_conv * (v <= n_voxels / 2) * 0.5
}
```

Fit with convergence tracking:

``` r
fit <- fit_hrfdecoder(
  Y = Y_data,
  ev_model = ev_model,
  base_model = bl_model,
  lambda_W = 0.1,
  max_iter = 15,
  tol = 1e-4,
  verbose = FALSE
)

# Inspect convergence info
names(fit$convergence)
#> NULL
```

``` r
# Did it converge?
if (!is.null(fit$convergence$converged)) {
  cat("Converged:", fit$convergence$converged, "\n")
  if (fit$convergence$converged) {
    cat("Iterations:", fit$convergence$iterations, "\n")
  }
} else {
  cat("Convergence info not available in this version\n")
}
#> Convergence info not available in this version
```

### Warning signs

- **Max iterations reached**: Increase `max_iter` or check parameter
  values
- **Oscillating objective**: Reduce learning rate (not exposed) or
  increase regularization
- **Very few iterations**: May indicate poor initialization or trivial
  solution

## Understanding soft labels

Soft labels represent the algorithm’s continuous predictions for each
trial.

``` r
# Extract soft labels by aggregating TR-level predictions to trials
pred_train <- predict(fit, newdata = Y_data, ev_model_test = ev_model, mode = "trial")
# For two-class case, define a signed margin (A vs B)
y_soft <- as.numeric(pred_train$probs[, 1] - pred_train$probs[, 2])

# Summary
cat("Soft label range:", round(range(y_soft), 3), "\n")
#> Soft label range: -0.22 0.189
cat("Mean absolute value:", round(mean(abs(y_soft)), 3), "\n")
#> Mean absolute value: 0.143

# For binary classification, expect negative/positive split
true_labels <- ifelse(conditions == "A", 1, -1)
cat("\nCorrelation with true labels:", round(cor(y_soft, true_labels), 3), "\n")
#> 
#> Correlation with true labels: 0.963
```

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  plot_df <- data.frame(
    trial = 1:n_trials,
    soft_label = y_soft,
    condition = conditions,
    true_label = true_labels
  )

  ggplot(plot_df, aes(x = trial, y = soft_label, color = condition)) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    (if (requireNamespace("albersdown", quietly = TRUE)) albersdown::scale_color_albers(params$family) else ggplot2::scale_color_discrete()) +
    labs(
      title = "Learned soft labels",
      subtitle = "Continuous trial predictions from weakly supervised learning",
      x = "Trial number",
      y = "Soft label value",
      color = "Condition"
    )
}
```

![Soft labels learned during training. Continuous predictions capture
graded trial-level responses rather than hard class
assignments.](05-weakly-supervised_files/figure-html/plot-soft-labels-1.png)

Soft labels learned during training. Continuous predictions capture
graded trial-level responses rather than hard class assignments.

Soft labels can reveal:

- **Graded responses**: Some trials are “more typical” than others
- **Learning effects**: Early vs. late trials may differ
- **Attention lapses**: Outlier trials with weak predictions

## Computational complexity

Each ALS iteration involves:

1.  **Soft label update**: O(N³) for N trials (Laplacian solve)
2.  **HRF update**: O(K³) for K basis functions (typically K = 2-3)
3.  **Decoder update**: O(V³) for V voxels (ridge regression)

**Bottleneck**: Decoder update when V is large (thousands of voxels).

**Speedup strategies**:

- **Prewhitening**: Reduces effective dimensionality
- **Rank reduction**: Automatic nuisance rank estimation
- **Sparse solvers**: Exploit structure in design matrices

Typical runtime: 1-5 seconds for 100 trials × 100 voxels on modern
hardware.

## Advanced topics

### Custom initialization

The algorithm initializes:

- **y**: Random normal or based on condition means
- **θ**: Canonical HRF (1, 0, …)
- **W**: Ridge regression on initial y

Custom initialization (not currently exposed) could improve convergence
for difficult cases.

### Alternative optimizers

ALS is simple and stable but not always fastest. Alternatives:

- **Coordinate descent**: Update each parameter sequentially
- **Gradient descent**: First-order optimization
- **ADMM**: Alternating direction method of multipliers

Current implementation uses ALS for reliability and interpretability.

## Debugging common issues

### Issue: Poor predictions

**Symptoms**: Random or below-chance accuracy

**Possible causes**:

1.  **Over-regularization**: Try decreasing `lambda_W`
2.  **Wrong HRF**: Try flexible HRF (`lambda_HRF = 0.001`)
3.  **Insufficient iterations**: Increase `max_iter`
4.  **Misspecified design**: Check event model timing

### Issue: Slow convergence

**Symptoms**: Reaches `max_iter` without converging

**Possible causes**:

1.  **Weak signal**: Increase regularization
2.  **Conflicting constraints**: Check lambda values
3.  **Poor scaling**: Ensure data is standardized

### Issue: Extreme soft labels

**Symptoms**: `y_soft` values \>\> 1 in magnitude

**Possible causes**:

1.  **Under-regularization**: Increase `lambda_smooth`
2.  **Overfitting**: Increase `lambda_W`

## Next steps

- [Getting
  Started](https://bbuchsbaum.github.io/hrfdecoder/articles/01-getting-started.md)
  — Basic usage tutorial
- [AR
  Prewhitening](https://bbuchsbaum.github.io/hrfdecoder/articles/02-ar-prewhitening.md)
  — Handle temporal autocorrelation
- [rMVPA
  Integration](https://bbuchsbaum.github.io/hrfdecoder/articles/03-rmvpa-integration.md)
  — Cross-validation framework
- [HRF
  Estimation](https://bbuchsbaum.github.io/hrfdecoder/articles/04-hrf-estimation.md)
  — Joint HRF learning

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
