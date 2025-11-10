# hrfdecoder Integration with rMVPA: Comprehensive Analysis

## Overview

This directory contains detailed architectural analysis of the rMVPA
package with focus on how the **hrfdecoder** continuous-time decoder
model integrates into the framework.

## Documents in This Analysis

1.  **`rMVPA_ARCHITECTURE_REPORT.md`** (main document)
    - Complete 9-section report covering all aspects of rMVPA
      architecture
    - Detailed data flow diagrams
    - Code signatures and patterns
    - Integration recommendations
    - **~3,000 lines of detailed documentation**
2.  **`RMVPA_KEY_FINDINGS.md`** (quick reference)
    - Quick file reference with line numbers
    - Critical integration points
    - Data structure flows
    - Implementation checklist
    - **Best for quick lookups and navigation**

## Key Architectural Insights

### 1. S3 Plugin System

rMVPA extends models via S3 method dispatch. Core operations:

``` r
train_model(model_spec, data, labels, ...)          # S3 generic
format_result(model_spec, fit, error, context, ...) # S3 generic
merge_results(model_spec, fold_results, ...)        # S3 generic
compute_performance(model_spec, classification_result) # S3 generic
```

**For hrfdecoder**: Implement 5 methods to plug into pipeline: -
`train_model.hrfdecoder_model()` - trains solver, returns fit object -
`format_result.hrfdecoder_model()` - predicts + aggregates to events -
`merge_results.hrfdecoder_model()` - combines fold results -
`compute_performance.hrfdecoder_model()` - extracts metrics -
`y_train.hrfdecoder_model()` - returns dummy sequence (no real labels)

### 2. Data Flow: Strict Temporal Ordering

    ROI Extracted (TR x features)
        ↓
    [PREPROCESSING OPPORTUNITY] ← AR whitening goes here
        ↓
    CV Folds Created (based on block_var/runs)
        ↓
    For Each Fold:
      ├─ train_model() [train data only]
      ├─ format_result() [test data, using train params]
      ├─ Result: tibble(class, probs, y_true, test_ind)
        ↓
    merge_results() [combine folds, compute metrics]
        ↓
    Final: classification_result with event-level predictions

**Critical**: Data extraction (step 1) occurs BEFORE fold creation. This
determines where preprocessing can happen.

### 3. Preprocessing Timing (Two Options)

**Option A: Per-ROI, Before Folds (PREFERRED)**

    Extract ROI → AR estimation → Prewhiten ROI → Create folds → Train/test

- Pros: Stable AR estimation, whole-ROI context
- Cons: Must apply params per-fold in train_model

**Option B: Per-Fold, Inside train_model**

    Create folds → For each fold: estimate AR, prewhiten, fit → apply params to test

- Pros: Strict train/test separation
- Cons: Slower (AR estimated per fold), harder parameter tuning

### 4. hrfdecoder-Specific Design

**Special Extension**: `hrfdecoder_design` subclasses `mvpa_design` -
Stores event metadata: `event_model`, `events` - Dummy `y_train` (TR
indices 1:T) for fold construction - Actual targets come from event
table, NOT y_train - Fold assignment via `block_var` (runs), NOT y
values - Ensures temporal structure preserved in cross-validation

**Example Usage**:

``` r
design <- hrfdecoder_design(
  event_model = evm,        # from fmridesign
  events = events_df,       # with onset, condition
  block_var = run_ids       # 1:3, determines folds
)

mspec <- hrfdecoder_model(
  dataset = dset,
  design = design,
  lambda_W = 10,
  window = c(4, 8),
  # ... other parameters
)

results <- run_searchlight(mspec, radius=8)
```

### 5. Where AR Prewhitening Fits

**In `train_model.hrfdecoder_model()`**: 1. Receive training data matrix
(ROI: obs x voxels) 2. Estimate AR model per voxel 3. Apply whitening to
training data 4. Fit decoder on whitened training data 5. Store AR
parameters in fit object

**In `format_result.hrfdecoder_model()`**: 1. Receive test data matrix
2. Apply SAME AR transformation (using training params) 3. Predict on
whitened test data 4. Aggregate TR-level predictions to event level 5.
Return classification_result

**Key**: AR parameters estimated from training data only, applied
identically to test data.

## Integration Checklist

For implementing AR prewhitening in hrfdecoder:

**Add ar_order parameter** to
[`hrfdecoder_model()`](https://bbuchsbaum.github.io/hrfdecoder/reference/hrfdecoder_model.md)
signature

- Default: `ar_order = 0` (no whitening)
- Callable: `hrfdecoder_model(..., ar_order = 2, ...)`

**Implement AR model fitting**

- Per-voxel AR(p) estimation on training data
- Store coefficients in fit object

**Implement whitening application**

- Forward: residuals from AR model (training data)
- Inverse: apply same transformation to test data

**Modify `train_model.hrfdecoder_model()`**

- Lines 175-210 in `/R/hrfdecoder_model.R`
- Add AR estimation before `hrfdecoder_fit()`
- Store `ar_params` in returned fit object

**Modify `format_result.hrfdecoder_model()`**

- Lines 216-297 in `/R/hrfdecoder_model.R`
- Apply AR transformation to test data before prediction
- Use AR params from training fold

**Validation tests**

- Verify train/test independence (no data leakage)
- Cross-validation with/without AR
- Parameter sensitivity analysis

## Key Code Locations

| Task | File | Lines | Function |
|----|----|----|----|
| Model spec creation | `/R/hrfdecoder_model.R` | 97-146 | [`hrfdecoder_model()`](https://bbuchsbaum.github.io/hrfdecoder/reference/hrfdecoder_model.md) |
| Training | `/R/hrfdecoder_model.R` | 175-210 | `train_model.hrfdecoder_model()` |
| Prediction | `/R/hrfdecoder_model.R` | 216-297 | `format_result.hrfdecoder_model()` |
| Fold merging | `/R/hrfdecoder_model.R` | 303-328 | `merge_results.hrfdecoder_model()` |
| CV iteration | `/R/mvpa_iterate.R` | 224-296 | `internal_crossval()` |
| Entry point | `/R/allgeneric.R` | 206-222 | `process_roi.default()` |
| S3 generics | `/R/allgeneric.R` | 351-708 | `train_model`, `format_result`, etc. |
| Model spec factory | `/R/mvpa_model.R` | 196-213 | `create_model_spec()` |
| CV specification | `/R/crossval.R` | 402-546 | [`blocked_cross_validation()`](http://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md) |

## Method Signatures Reference

### train_model

``` r
train_model.hrfdecoder_model <- function(
  obj,          # hrfdecoder_model spec
  train_dat,    # tibble/matrix: obs x features (already split by fold)
  y,            # vector: dummy sequence (ignored by hrfdecoder)
  sl_info,      # list: center_local_id, center_global_id (searchlight info)
  cv_spec,      # cross_validation spec
  indices,      # integer: global voxel indices
  ...           # additional arguments
)
# Returns: object of class "hrfdecoder_fit_wrap" with fit, ar_params, indices
```

### format_result

``` r
format_result.hrfdecoder_model <- function(
  obj,              # hrfdecoder_model spec
  result,           # fit object from train_model
  error_message,    # NULL if no error
  context,          # list with test, ytest, roi, ytrain, train, .id
  ...               # additional arguments
)
# Returns: tibble with class, probs, y_true, test_ind, fit, error, error_message
```

### merge_results

``` r
merge_results.hrfdecoder_model <- function(
  obj,          # hrfdecoder_model spec
  result_set,   # tibble: results from all folds (from format_result)
  indices,      # integer: global voxel indices for ROI
  id,           # identifier for this ROI
  ...           # additional arguments
)
# Returns: tibble with result, indices, performance, id, error, error_message
```

## File Organization

    rMVPA/R/
    ├── allgeneric.R          [944 lines] Generic functions & S3 dispatch
    ├── mvpa_model.R          [387 lines] Standard MVPA model specs
    ├── hrfdecoder_model.R    [336 lines] hrfdecoder adapter ← MODIFY THIS
    ├── hrfdecoder_design.R   [174 lines] hrfdecoder design extension
    ├── mvpa_iterate.R        [585 lines] CV iteration & fold generation
    ├── crossval.R            [903 lines] CV specifications
    ├── common.R              [624 lines] Normalization utilities
    ├── dataset.R             [479 lines] Data structures
    ├── design.R              [344 lines] Design specifications
    ├── classifiers.R         [739 lines] Model registry
    ├── searchlight.R        [1130 lines] Searchlight execution
    └── regional.R            [549 lines] Regional analysis

## Performance Considerations

### Current Bottlenecks

1.  **AR Estimation per Fold**: If AR estimated separately for each
    fold, will be slower
    - Mitigation: Estimate once per ROI, apply per-fold
2.  **AR Parameter Storage**: Need efficient serialization in fit object
    - Current: Simple list assignment (should be fine)
3.  **Data Copying**: Prewhitening creates new data copies
    - Mitigation: Use sparse matrices where possible

### Optimization Opportunities

- Cache AR estimates if they’re stable across folds
- Vectorize per-voxel AR estimation (use RcppArmadillo)
- Parallelize fold-level prewhitening

## Testing Strategy

``` r
# 1. Unit tests for AR functions
test_that("AR whitening preserves shape", { ... })
test_that("AR parameters recoverable", { ... })

# 2. Integration tests for train_model
test_that("train_model returns ar_params", { ... })
test_that("format_result applies ar_params correctly", { ... })

# 3. Cross-validation tests
test_that("train/test data independence", { ... })
test_that("AR estimation only uses training data", { ... })

# 4. End-to-end tests
test_that("searchlight with AR produces valid results", { ... })
test_that("performance metrics computed correctly", { ... })
```

## References

- **rMVPA Main Report**: See `rMVPA_ARCHITECTURE_REPORT.md` for
  comprehensive documentation
- **Quick Reference**: See `RMVPA_KEY_FINDINGS.md` for file locations
  and navigation
- **hrfdecoder Package**: Continuous-time decoding solver with ALS
  optimization
- **fmridesign Package**: Event model construction and HRF convolution
- **neuroim2 Package**: Neuroimaging data structures (NeuroVec,
  NeuroVol)

## Questions & Answers

**Q: Where does preprocessing happen in the pipeline?** A: Two options
documented in the main report (Section 5.1). Either per-ROI before fold
creation (preferred) or per-fold inside
[`train_model()`](http://bbuchsbaum.github.io/rMVPA/reference/train_model.md)
(stricter but slower).

**Q: How are train/test splits handled?** A: By `block_var` (run IDs) in
[`blocked_cross_validation()`](http://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md).
Each run is held out completely, no temporal mixing. AR parameters
estimated from training runs only, applied to held-out runs.

**Q: Can preprocessing be model-specific?** A: Yes! S3 method dispatch
allows `preprocess_roi.hrfdecoder_model()` to exist alongside
`preprocess_roi.mvpa_model()`, etc. Each model can implement its own
preprocessing.

**Q: How to ensure prewhitening happens before train_model?** A: Modify
`process_roi.default()` in `/R/allgeneric.R` to call
`roi <- preprocess_roi(mod_spec, roi)` before `internal_crossval()`.

------------------------------------------------------------------------

**Report Generated**: November 9, 2024 **rMVPA Version Analyzed**:
Current (dev) **Scope**: Complete architectural analysis (7,000+ lines
of code reviewed)
