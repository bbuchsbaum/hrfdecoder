# Integrating with rMVPA workflows

## Goal

Learn how to integrate `hrfdecode` with the `rMVPA` package for
whole-brain searchlight analysis and cross-validation.

## TL;DR

``` r
library(hrfdecode)
library(rMVPA)
library(fmridesign)
library(fmrihrf)

# Build event_model and design using helper
sframe <- fmrihrf::sampling_frame(TR = 2, blocklens = c(200, 200))
ev_model <- fmridesign::event_model(
  onset ~ fmridesign::hrf(condition, basis = "spmg1"),
  data = events,
  block = ~ run,
  sampling_frame = sframe
)
block_ids <- rep(1:2, each = 200)
design <- hrfdecode::continuous_mvpa_design(
  event_model = ev_model,
  block_var = block_ids,
  design_df_events = events
)

# Create hrfdecoder model specification
model <- hrfdecoder_model(
  design = design,
  lambda_W = 0.1,
  ar_order = 1
)

# Run searchlight with cross-validation
result <- run_searchlight(
  model_spec = model,
  radius = 8,
  method = "randomized"
)
```

## Background

The `rMVPA` package provides a framework for:

- **Searchlight analysis** — sliding-window decoding across the brain
- **ROI-based analysis** — regional classification
- **Cross-validation** — train/test splits with multiple folds

Traditional rMVPA workflows require **trial-averaged data** and
**explicit labels**. The `hrfdecode` integration enables **continuous
MVPA**: decoding directly from time series using weakly supervised
learning.

### Key differences

| Aspect | Traditional rMVPA | hrfdecoder + rMVPA |
|----|----|----|
| **Input data** | Trial-averaged beta maps | Continuous fMRI time series |
| **Labels** | Explicit trial labels | Event onsets + conditions |
| **Temporal structure** | Ignored (trials independent) | Modeled via HRF + soft labels |
| **Cross-validation** | Across trials | Across runs |

## Setup

``` r
library(hrfdecode)
library(fmridesign)
# library(rMVPA)  # Load if available
```

## Creating a continuous design

The
[`continuous_mvpa_design()`](https://bbuchsbaum.github.io/hrfdecoder/reference/continuous_mvpa_design.md)
function wraps
[`fmridesign::event_model()`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html)
to create an rMVPA-compatible design object.

``` r
# Example event table (multi-run experiment)
n_trials_per_run <- 20
tr <- 2
n_trs_per_run <- 200

event_table <- data.frame(
  onset = c(
    seq(10, (n_trs_per_run - 10) * tr, length.out = n_trials_per_run),
    seq(10, (n_trs_per_run - 10) * tr, length.out = n_trials_per_run)
  ),
  condition = rep(c("face", "scene"), n_trials_per_run),
  duration = 1,
  run = rep(1:2, each = n_trials_per_run)
)

head(event_table, 3)
#>      onset condition duration run
#> 1 10.00000      face        1   1
#> 2 29.47368     scene        1   1
#> 3 48.94737      face        1   1
```

Create the continuous design:

``` r
# Build sampling frame and event_model for two runs
sframe <- fmrihrf::sampling_frame(TR = tr, blocklens = c(n_trs_per_run, n_trs_per_run))
ev_model <- fmridesign::event_model(
  onset ~ fmridesign::hrf(condition, basis = "spmg1"),
  data = event_table,
  block = ~ run,
  sampling_frame = sframe
)

# Block/run ids per TR (length equals total TRs across runs)
block_ids <- rep(1:2, each = n_trs_per_run)

# Create hrfdecoder design via helper
design <- hrfdecode::continuous_mvpa_design(
  event_model = ev_model,
  block_var = block_ids,
  design_df_events = event_table
)

# Inspect design object
class(design)
#> [1] "mvpa_design" "list"
```

The design contains:

- `event_model`: `fmridesign` event model with HRF basis
- `baseline_model`: Nuisance/baseline regressor model
- `block_var`: Blocking structure for cross-validation (e.g., runs)

## Creating an hrfdecoder model

Use
[`hrfdecoder_model()`](https://bbuchsbaum.github.io/hrfdecoder/reference/hrfdecoder_model.md)
to create an rMVPA-compatible model specification.

``` r
model <- hrfdecoder_model(
  design = design,
  lambda_W = 0.1,      # Ridge penalty on decoder weights
  lambda_HRF = 0.01,   # Prior adherence to HRF basis
  ar_order = 1,        # AR(1) prewhitening
  max_iter = 10        # ALS iterations
)

class(model)
#> [1] "hrfdecoder_model"
```

This model specification implements the rMVPA interface:

- [`train_model()`](http://bbuchsbaum.github.io/rMVPA/reference/train_model.md):
  Fit decoder on training fold
- [`predict()`](https://rdrr.io/r/stats/predict.html): Apply decoder to
  test fold (now using
  [`format_result()`](http://bbuchsbaum.github.io/rMVPA/reference/format_result.md))
- [`merge_results()`](http://bbuchsbaum.github.io/rMVPA/reference/merge_results.md):
  Combine cross-validation results

## Cross-validation structure

With continuous data, **cross-validation happens across runs**, not
trials. Each fold:

1.  **Training**: All runs except one
2.  **Testing**: Held-out run

This respects temporal dependencies within runs while providing
independent test sets.

``` r
# Conceptual fold structure for 3 runs:
# Fold 1: Train on runs 2,3 | Test on run 1
# Fold 2: Train on runs 1,3 | Test on run 2
# Fold 3: Train on runs 1,2 | Test on run 3
```

## Running searchlight analysis

> **Note**: The following examples are conceptual. Full searchlight
> integration requires an `rMVPA` installation and appropriate data
> structures.

``` r
library(rMVPA)

# Prepare 4D fMRI dataset (X × Y × Z × T)
# Assume fmri_4d is a neuroim2::NeuroVec object
# fmri_4d <- neuroim2::read_vec("func_data.nii.gz")

# Run searchlight
result <- run_searchlight(
  model_spec = model,
  dataset = fmri_4d,
  radius = 8,          # 8mm searchlight sphere
  method = "randomized",
  niter = 4            # Iterations for randomized searchlight
)

# Extract performance map
perf_map <- result$performance

# Threshold at accuracy > 60%
sig_map <- perf_map > 0.60
```

### Searchlight workflow

For each searchlight sphere:

1.  **Extract voxel time series** within radius
2.  **Fit hrfdecoder** using
    [`train_model()`](http://bbuchsbaum.github.io/rMVPA/reference/train_model.md)
    on training folds
3.  **Predict** on test fold using learned decoder, HRF, AR
4.  **Aggregate predictions** to trial level
5.  **Compute accuracy** (or other metric)
6.  **Assign performance** to center voxel

The result is a whole-brain performance map showing where decoding
succeeds.

## ROI-based analysis

For hypothesis-driven analyses, use predefined ROIs instead of
searchlight.

``` r
# Define ROI mask (e.g., fusiform face area)
# roi_mask <- neuroim2::read_vol("ffa_mask.nii.gz")

# Extract ROI time series
# Y_roi <- extract_roi_timeseries(fmri_4d, roi_mask)

# Fit decoder on ROI
fit_roi <- fit_hrfdecoder(
  Y = Y_roi,
  event_model = design$event_model,
  baseline_model = design$baseline_model,
  lambda_W = 0.1,
  ar_order = 1
)

# Cross-validate manually
# (or use rMVPA::crossval_* functions with hrfdecoder_model)
```

## Merging fold results

After cross-validation, results from each fold need to be combined.

``` r
# Conceptual: After running cross-validation
# fold_results <- list(fold1_pred, fold2_pred, fold3_pred)

# Merge using hrfdecoder's merge method
# merged <- merge_results(model, fold_results)

# Compute overall accuracy
# overall_acc <- mean(merged$predictions == merged$true_labels)
```

The
[`merge_results()`](http://bbuchsbaum.github.io/rMVPA/reference/merge_results.md)
method concatenates predictions across folds, respecting the temporal
structure of test sets.

## Handling multi-run data

Multi-run experiments require careful handling:

1.  **Event table**: Include `run` column
2.  **Blocklens**: Specify TR counts per run
3.  **Baseline model**: Separate baseline per run
4.  **AR parameters**: Optionally pool per run (if implemented)

``` r
# Already created above with run column
str(event_table)
#> 'data.frame':    40 obs. of  4 variables:
#>  $ onset    : num  10 29.5 48.9 68.4 87.9 ...
#>  $ condition: chr  "face" "scene" "face" "scene" ...
#>  $ duration : num  1 1 1 1 1 1 1 1 1 1 ...
#>  $ run      : int  1 1 1 1 1 1 1 1 1 1 ...

# Blocklens match number of runs
blocklens <- c(n_trs_per_run, n_trs_per_run)

# Design respects run structure
# Design respects run structure (via rMVPA helper)
sframe_multi <- fmrihrf::sampling_frame(TR = tr, blocklens = blocklens)
ev_model_multi <- fmridesign::event_model(
  onset ~ fmridesign::hrf(condition, basis = "spmg1"),
  data = event_table,
  block = ~ run,
  sampling_frame = sframe_multi
)
block_ids_multi <- rep(1:2, each = n_trs_per_run)
design_multi <- hrfdecode::continuous_mvpa_design(
  event_model = ev_model_multi,
  block_var = block_ids_multi,
  design_df_events = event_table
)
```

## Interpreting results

Searchlight and ROI results provide different insights:

### Searchlight maps

- **High accuracy regions**: Information-bearing areas
- **Cluster extent**: Spatial distribution of coding
- **Peak locations**: Anatomical specificity

### ROI performance

- **Above-chance accuracy**: Evidence for information coding
- **Comparison across ROIs**: Selectivity and specificity
- **Temporal dynamics**: Via TR-level predictions

## Practical tips

1.  **Start with ROIs**: Test decoder on known regions before
    whole-brain
2.  **Check convergence**: Monitor ALS iterations (`verbose = TRUE`)
3.  **Tune regularization**: Cross-validate `lambda_W` and `lambda_HRF`
4.  **Verify predictions**: Inspect `y_soft` and `theta` for sanity
5.  **Use parallel processing**: rMVPA supports parallel searchlight

## Next steps

- [Getting
  Started](https://bbuchsbaum.github.io/hrfdecoder/articles/01-getting-started.md)
  — Basic decoder fitting
- [AR
  Prewhitening](https://bbuchsbaum.github.io/hrfdecoder/articles/02-ar-prewhitening.md)
  — Temporal autocorrelation correction
- [HRF
  Estimation](https://bbuchsbaum.github.io/hrfdecoder/articles/04-hrf-estimation.md)
  — Understanding joint HRF learning
- [Weakly Supervised
  Learning](https://bbuchsbaum.github.io/hrfdecoder/articles/05-weakly-supervised.md)
  — Algorithm details

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
