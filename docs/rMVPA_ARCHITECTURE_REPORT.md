# rMVPA Package Architecture Report: Integration Points for hrfdecoder

## Executive Summary

This report documents the rMVPA plugin architecture with specific focus
on how the **hrfdecoder** decoder model integrates. The analysis covers
(1) the S3 model plugin system, (2) data flow through
searchlight/cross-validation pipelines, (3) existing preprocessing
capabilities, and (4) architectural patterns for extending the
framework.

------------------------------------------------------------------------

## 1. Model Plugin Architecture

### 1.1 S3 Method Dispatch System

rMVPA uses S3 classes for extensibility. Core model operations dispatch
through:

**Generic Functions** (defined in `/R/allgeneric.R`): -
`train_model(obj, ...)` → trains model on ROI/fold data -
`predict_model(object, fit, newdata, ...)` → generates predictions -
`format_result(obj, result, error_message, context, ...)` → formats
fold-level outputs - `merge_results(obj, result_set, indices, id, ...)`
→ aggregates fold results to ROI-level performance -
`compute_performance(obj, result)` → extracts metrics from result

### 1.2 Model Spec Creation: `create_model_spec()` and `mvpa_model()`

**Function**: `/R/mvpa_model.R:196-213` and `/R/mvpa_model.R:260-337`

``` r
create_model_spec(name, dataset, design, return_predictions=FALSE, 
                  compute_performance=FALSE, tune_reps=FALSE, ...)
# Returns: S3 object with class c(name, "model_spec", "list")

mvpa_model(model, dataset, design, model_type=c("classification", "regression"), 
           crossval=NULL, feature_selector=NULL, tune_grid=NULL, ...)
# Returns: S3 object with class c("mvpa_model", "model_spec", "list")
```

**Key Fields in Model Spec**:

``` r
model_spec$dataset          # mvpa_dataset with train_data, test_data, mask
model_spec$design           # mvpa_design with y_train, y_test, block_var
model_spec$crossval         # cross_validation spec (blocked_cross_validation, etc.)
model_spec$compute_performance  # logical: compute metrics after fold merging
model_spec$return_predictions   # logical: return per-fold predictions
model_spec$return_fits      # logical: keep fitted models
model_spec$performance      # function: (result) -> named_metrics
model_spec$has_test_set     # logical: external test set present
```

### 1.3 hrfdecoder Plugin: Specialized Model Spec

**File**: `/R/hrfdecoder_model.R:97-146`

``` r
hrfdecoder_model <- function(
  dataset,           # mvpa_dataset with TR x V data
  design,            # hrfdecoder_design (custom mvpa_design subclass)
  lambda_W = 10,     # ridge penalty for decoder weights
  lambda_HRF = 1,    # HRF conformity penalty
  lambda_smooth = 5, # temporal smoothness penalty
  basis = NULL,      # HRF basis (fmrihrf::spmg1, etc.)
  window = c(4, 8),  # event aggregation window (seconds)
  nonneg = TRUE,     # project soft labels to [0,1]
  max_iter = 10,     # ALS iterations
  tol = 1e-4,        # convergence tolerance
  performance = NULL,
  crossval = NULL,
  return_predictions = TRUE,
  return_fits = FALSE
)
# Returns: S3 object with class c("hrfdecoder_model", "model_spec", "list")
```

**Key Differences from Standard MVPA Models**: - Operates on TR-level
(continuous) data, not trial-level betas -
[`y_train()`](http://bbuchsbaum.github.io/rMVPA/reference/y_train-methods.md)
returns dummy TR sequence (1:T) for fold construction - Actual targets
come from `design$event_model` and `design$events` (ignored by fold
machinery) - Fold assignment determined by `block_var` (run IDs), not
`y_train` values

### 1.4 S3 Method Implementation Pattern

Example from hrfdecoder:

``` r
# STEP 1: y_train extraction (for fold construction)
y_train.hrfdecoder_model <- function(obj) {
  seq_len(nobs(obj$dataset))  # Dummy sequence 1:T
}

# STEP 2: Model training (ignored y parameter, uses design metadata)
train_model.hrfdecoder_model <- function(obj, train_dat, y, sl_info, cv_spec, indices, ...) {
  X <- as.matrix(train_dat)
  evm <- obj$design$event_model
  events <- obj$design$events
  fit <- hrfdecoder::hrfdecoder_fit(X, event_model=evm, basis=obj$basis, ...)
  structure(list(fit=fit, sl_info=sl_info, indices=indices), 
            class="hrfdecoder_fit_wrap")
}

# STEP 3: Format fold result (predict + aggregate to events)
format_result.hrfdecoder_model <- function(obj, result, error_message=NULL, context, ...) {
  Xtest <- as.matrix(context$test)
  Ptest <- predict(result$fit, Xtest)  # TR-level soft labels
  agg <- hrfdecoder::aggregate_events(Ptest, obj$design$events, window=obj$window)
  # Return tibble with event-level probs and predictions
}

# STEP 4: Merge fold results
merge_results.hrfdecoder_model <- function(obj, result_set, indices, id, ...) {
  # Concatenate fold results, build classification_result, compute performance
}

# STEP 5: Compute performance
compute_performance.hrfdecoder_model <- function(obj, result) {
  obj$performance(result)  # Delegates to perf_fun
}
```

------------------------------------------------------------------------

## 2. Data Flow: Searchlight and Cross-Validation

### 2.1 ROI/Searchlight Extraction → Cross-Validation Folds → Training/Testing

**Main Pipeline** (from `/R/mvpa_iterate.R` and `/R/allgeneric.R`):

    process_roi(mod_spec, roi, rnum, center_global_id)
      ├─→ Check for external test set (has_test_set)
      │    └─→ external_crossval(mspec, roi, id, center_global_id)
      ├─→ Check for internal cross-validation (has_crossval)
      │    └─→ internal_crossval(mspec, roi, id, center_global_id)
      └─→ Fallback: no CV
           └─→ process_roi_default(mspec, roi, rnum, center_global_id)

### 2.2 Internal Cross-Validation (Most Common)

**Function**: `/R/mvpa_iterate.R:224-296`

``` r
internal_crossval <- function(mspec, roi, id, center_global_id = NA) {
  # 1. EXTRACT DATA: Get train_roi from searchlight/region
  xtrain <- neuroim2::values(roi$train_roi)  # ROI matrix: obs x features
  ind <- neuroim2::indices(roi$train_roi)     # Global voxel indices
  
  # 2. GENERATE FOLDS: Create train/test splits based on crossval spec
  samples <- crossval_samples(
    mspec$crossval,
    tibble::as_tibble(xtrain),
    y_train(mspec)  # For hrfdecoder: dummy TR sequence
  )
  # Returns tibble with columns: ytrain, ytest, train, test, .id
  
  # 3. ITERATE FOLDS: For each fold
  ret <- samples %>% pmap(function(ytrain, ytest, train, test, .id) {
    # 3a. TRAIN
    result <- train_model(
      mspec, 
      tibble::as_tibble(train),
      ytrain,           # For hrfdecoder: dummy sequence
      sl_info = list(center_local_id, center_global_id),
      cv_spec = mspec$crossval,
      indices = ind
    )
    
    # 3b. FORMAT (predict + aggregate)
    format_result(
      mspec, result, 
      error_message = NULL,
      context = list(
        test = test,      # Test ROI data
        ytest = ytest,    # Test labels/dummy
        ...
      )
    )
  }) %>% bind_rows()
  
  # 4. MERGE FOLDS: Aggregate across folds
  merge_results(mspec, ret, indices=ind, id=id)
  # Returns: tibble with result, indices, performance, id, error
}
```

**Key Insight**: Data extraction (step 1) happens **before** fold
creation. This means: - **Preprocessing MUST happen before fold
iteration** to avoid data leakage - Preprocessing can be done per-ROI,
before step 2 - No per-fold preprocessing without custom handling

### 2.3 Cross-Validation Specification

**File**: `/R/crossval.R:402-546`

``` r
blocked_cross_validation(block_var)
# Standard: leave one block (run) out at a time
# Returns: object with class "blocked_cross_validation"
# Used by: crossval_samples() → creates folds

sequential_blocked_cross_validation(block_var, nfolds=2, nreps=4)
# Advanced: subdivides each block further

custom_cross_validation(sample_set)
# User provides explicit train/test index pairs
```

Implementation:

``` r
crossval_samples.blocked_cross_validation <- function(obj, data, y, ...) {
  crossv_block(data, y, obj$block_var)
  # Returns: tibble(ytrain, ytest, train, test, .id) for each fold
}
```

### 2.4 Data Sample Structure

The extracted ROI data is wrapped in **modelr::resample** objects
(lightweight index-based views):

``` r
# Inside crossval_samples:
list(
  ytrain = y[train_indices],
  ytest = y[test_indices],
  train = modelr::resample(data, train_indices),  # Lazy view
  test = modelr::resample(data, test_indices)      # Lazy view
)
```

Unpacking for model training:

``` r
xtrain <- tibble::as_tibble(train, .name_repair=.name_repair)
# Converts to matrix/tibble of actual values
```

------------------------------------------------------------------------

## 3. Existing Preprocessing Capabilities

### 3.1 Dataset-Level Normalization

**File**: `/R/common.R:52-88`

``` r
# Per-volume z-scoring (before fold creation)
normalize_image_samples <- function(bvec, mask) {
  # Scale each volume independently
  vlist <- vols(bvec) %>% future_map(function(v) {
    scale(v[mask>0])[,1]
  })
  SparseNeuroVec(do.call(cbind, vlist), space(bvec), mask=mask)
}

# Per-block standardization (within-run z-scoring)
standardize_vars <- function(bvec, mask, blockvar) {
  # Split by block, z-score within each block
  vlist <- vectors(bvec, subset=which(mask>0)) %>% future_map(function(v) {
    unlist(map(split(v, blockvar), scale))
  })
  SparseNeuroVec(do.call(cbind, vlist), space(bvec), mask=mask)
}
```

**When Applied**: During dataset initialization in the configuration
pipeline (not fold-aware).

### 3.2 Feature Selection

**File**: `/R/feature_selection.R` and `/R/allgeneric.R:36-64`

``` r
select_features <- function(obj, X, Y, ...)
# Generic called within internal_crossval before training
# Returns: logical vector of selected features

# Implemented in mvpa_model:
select_features.mvpa_model <- function(obj, X, Y,...) {
  if (!is.null(obj$feature_selector)) {
    select_features(obj$feature_selector, X, Y,...)
  } else {
    rep(TRUE, length(X))  # No selection
  }
}
```

**Timing**: Called per-fold **within** `internal_crossval` **after**
fold split. Avoids feature selection bias.

### 3.3 Nuisance Regression & AR Handling

**Current Status**: **NOT IMPLEMENTED** in rMVPA core.

- No built-in AR whitening, ARMA models, or nuisance regression
- Normalization (`normalize_samples` flag in config) is limited to
  z-scoring
- **This is a key integration point for hrfdecoder’s prewhitening**

------------------------------------------------------------------------

## 4. Configuration and Hyperparameters

### 4.1 Model Specification Flow

    mvpa_model()                          # Specify model, dataset, design
        ↓
    create_model_spec()                   # Build S3 object
        ↓
    run_searchlight(model_spec, ...)      # Execute searchlight
        ↓
    process_roi(model_spec, roi, ...)     # Per-ROI entry point
        ↓
    internal_crossval(model_spec, ...)    # CV iteration
        ↓
    train_model(model_spec, ...)          # Delegate to S3 method

### 4.2 Hyperparameter Specification

**For standard models** (e.g., `sda_notune`): - Passed via `model`
object in registry (MVPAModels) -
`fit(x, y, wts, param, lev, last, weights, classProbs, ...)` signature

**For hrfdecoder**: - Passed as explicit function arguments: `lambda_W`,
`lambda_HRF`, `lambda_smooth`, etc. - Stored directly in model spec -
Accessed in `train_model.hrfdecoder_model()` via `obj$lambda_W`, etc.

### 4.3 Cross-Validation Interaction with Preprocessing

**Key Issue**: CV folds created **after** ROI extraction but
preprocessing timing matters:

    ROI Extracted → [PREPROCESSING WINDOW] → Folds Created → Train/Test per Fold

**For AR Prewhitening**: 1. Must happen **before** fold creation to
avoid cross-fold contamination 2. Must use **training data only** per
fold (inside `train_model`) 3. Decoder weights fit on prewhitened
training data 4. Test predictions made on independently prewhitened test
data

------------------------------------------------------------------------

## 5. Integration Points for hrfdecoder

### 5.1 Where Prewhitening Fits

**Option A: Per-ROI, Before Fold Creation** (Preferred for AR)

``` r
process_roi_with_prewhitening <- function(mspec, roi, rnum, ...) {
  # Extract ROI data
  roi_data <- neuroim2::values(roi$train_roi)
  
  # PREWHITEN per ROI (not per fold)
  roi_data_pw <- prewhiten_ar(roi_data, mspec$ar_order)
  
  # Update roi with prewhitened data
  roi$train_roi <- update_roi_data(roi$train_roi, roi_data_pw)
  
  # Proceed with normal flow
  internal_crossval(mspec, roi, rnum, ...)
}
```

**Option B: Per-Fold, Inside `train_model`** (Fold-specific AR params)

``` r
train_model.hrfdecoder_model <- function(obj, train_dat, y, ...) {
  X <- as.matrix(train_dat)
  
  # Estimate AR model and prewhiten TRAINING data only
  X_pw <- prewhiten_ar_estimate(X, obj$ar_order)
  
  # Fit decoder on prewhitened data
  fit <- hrfdecoder::hrfdecoder_fit(X_pw, event_model=obj$design$event_model, ...)
  
  # Store AR parameters for test data
  structure(list(fit=fit, ar_params=attr(X_pw, "ar_params")), 
            class="hrfdecoder_fit_wrap")
}

format_result.hrfdecoder_model <- function(obj, result, error_message=NULL, context, ...) {
  Xtest <- as.matrix(context$test)
  
  # Apply SAME AR whitening to test data
  Xtest_pw <- apply_ar_whitening(Xtest, result$ar_params)
  
  Ptest <- predict(result$fit, Xtest_pw)
  # ... rest of aggregation
}
```

### 5.2 Design Extension: `hrfdecoder_design`

**File**: `/R/hrfdecoder_design.R:60-130`

``` r
hrfdecoder_design <- function(event_model, events, block_var, split_by=NULL) {
  # Creates mvpa_design subclass with additional fields
  
  mvdes <- mvpa_design(
    train_design = data.frame(y = seq_len(Tlen), block = block_var),
    y_train = ~ y,           # Dummy: TR indices
    block_var = ~ block,     # Run IDs - DETERMINES folds
    split_by = split_by
  )
  
  # Attach metadata
  mvdes$event_model <- event_model   # fmridesign object
  mvdes$events <- events              # event-level labels
  
  class(mvdes) <- c("hrfdecoder_design", class(mvdes))
  mvdes
}
```

**Key**: Block-level CV ensures entire runs are held out (respects
temporal structure).

### 5.3 Model-Specific Preprocessing Hook

**New Generic Function** (could be added):

``` r
preprocess_roi <- function(obj, roi, ...) {
  UseMethod("preprocess_roi")
}

preprocess_roi.default <- function(obj, roi, ...) {
  roi  # No preprocessing
}

preprocess_roi.hrfdecoder_model <- function(obj, roi, ...) {
  # Apply AR prewhitening
  roi_data <- neuroim2::values(roi$train_roi)
  roi_data_pw <- apply_prewhitening(roi_data, obj)
  
  # Update roi
  roi$train_roi <- replace_roi_data(roi$train_roi, roi_data_pw)
  roi
}
```

Called from `process_roi.default`:

``` r
process_roi.default <- function(mod_spec, roi, rnum, center_global_id = NA, ...) {
  roi <- preprocess_roi(mod_spec, roi)  # NEW HOOK
  
  xtrain <- neuroim2::values(roi$train_roi)
  # ... rest of pipeline
}
```

------------------------------------------------------------------------

## 6. Key Architectural Patterns

### 6.1 S3 Method Resolution Order

For hrfdecoder model operations:

    train_model(hrfdecoder_model, ...)
      → train_model.hrfdecoder_model()

    format_result(hrfdecoder_model, ...)
      → format_result.hrfdecoder_model()

    merge_results(hrfdecoder_model, ...)
      → merge_results.hrfdecoder_model()

Each method receives the full model_spec object, allowing access to all
parameters.

### 6.2 Preprocessing Timing Patterns

    PATTERN 1: Dataset-level (immediate, no folds)
      mvpa_dataset() → normalize_image_samples() → data ready

    PATTERN 2: Per-ROI (before fold iteration)
      process_roi() → [preprocess_roi() NEW] → internal_crossval()

    PATTERN 3: Per-Fold (inside train_model)
      internal_crossval() → train_model() → [local preprocessing] → fit()

    PATTERN 4: Test Data (must use training params)
      format_result() → [apply training preprocessing] → predict()

### 6.3 Result Aggregation Hierarchy

    Fold-Level Results (format_result):
      class: tibble with columns [class, probs, y_true, test_ind, fit, error, error_message]

    ROI-Level Results (merge_results):
      class: tibble with columns [result, indices, performance, id, error, error_message]
      - Combines all folds into single classification_result
      - Computes metrics via compute_performance()

    Searchlight Results (wrap_out, run_searchlight):
      class: searchlight_result (list of NeuroVol/NeuroSurface maps)

------------------------------------------------------------------------

## 7. File Organization Reference

| File | Purpose | Key Content |
|----|----|----|
| `allgeneric.R` | Generic function definitions | `train_model`, `predict_model`, `format_result`, `merge_results`, `process_roi` |
| `mvpa_model.R` | Standard MVPA model spec | [`mvpa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_model.md), `create_model_spec()`, S3 methods |
| `hrfdecoder_model.R` | hrfdecoder adapter | [`hrfdecoder_model()`](https://bbuchsbaum.github.io/hrfdecoder/reference/hrfdecoder_model.md), `train_model.hrfdecoder_model()`, `format_result.hrfdecoder_model()` |
| `hrfdecoder_design.R` | hrfdecoder design extension | `hrfdecoder_design()`, validation, metadata attachment |
| `mvpa_iterate.R` | CV iteration logic | `internal_crossval()`, `external_crossval()`, fold generation |
| `crossval.R` | CV specification classes | [`blocked_cross_validation()`](http://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md), [`crossval_samples()`](http://bbuchsbaum.github.io/rMVPA/reference/crossval_samples.md) |
| `model_fit.R` | Model tuning utilities | `tune_model()`, caret integration |
| `common.R` | Normalization utilities | `normalize_image_samples()`, `standardize_vars()` |
| `searchlight.R` | Searchlight iteration | [`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md), output wrapping |
| `regional.R` | Regional analysis | [`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md), region mask handling |
| `dataset.R` | Dataset structures | [`mvpa_dataset()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_dataset.md), [`gen_sample_dataset()`](http://bbuchsbaum.github.io/rMVPA/reference/gen_sample_dataset.md) |
| `design.R` | Design specifications | [`mvpa_design()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_design.md), formula parsing |
| `classifiers.R` | Model registry | `MVPAModels` environment, [`load_model()`](http://bbuchsbaum.github.io/rMVPA/reference/load_model.md) |

------------------------------------------------------------------------

## 8. Integration Checklist for AR Prewhitening

Implement `preprocess_roi()` generic function

Add `preprocess_roi.hrfdecoder_model()` method

Store AR parameters in fit object for test data application

Ensure `format_result.hrfdecoder_model()` applies same AR transform to
test data

Validate that prewhitening happens before fold creation (no data
leakage)

Document AR parameter control in
[`hrfdecoder_model()`](https://bbuchsbaum.github.io/hrfdecoder/reference/hrfdecoder_model.md)
signature

Test with cross-validation to ensure train/test independence

Add `ar_order` parameter to
[`hrfdecoder_model()`](https://bbuchsbaum.github.io/hrfdecoder/reference/hrfdecoder_model.md)
spec

Consider model-agnostic interface for other preprocessing (nuisance
regression, etc.)

------------------------------------------------------------------------

## 9. Recommended Preprocessing Architecture for hrfdecoder

``` r
#' Apply AR Prewhitening to ROI Data
#'
#' @param obj hrfdecoder_model spec
#' @param roi ROI data (list with train_roi)
#' @param ... Additional args
#' @return roi with prewhitened train_roi
preprocess_roi.hrfdecoder_model <- function(obj, roi, ...) {
  if (is.null(obj$ar_order) || obj$ar_order == 0) {
    return(roi)  # No whitening
  }
  
  roi_data <- neuroim2::values(roi$train_roi)
  
  # Fit AR model per voxel
  ar_fit <- fit_ar_model(roi_data, obj$ar_order)
  
  # Prewhiten data
  roi_data_pw <- apply_ar_whitening(roi_data, ar_fit)
  
  # Update roi, preserving indices and metadata
  roi$train_roi@data <- roi_data_pw
  attr(roi$train_roi, "ar_fit") <- ar_fit
  
  roi
}

# Inside train_model:
train_model.hrfdecoder_model <- function(obj, train_dat, y, sl_info, cv_spec, indices, ...) {
  X <- as.matrix(train_dat)
  
  # Check if ar_fit is pre-computed (from preprocess_roi)
  ar_fit <- attr(train_dat, "ar_fit")
  if (!is.null(ar_fit)) {
    # Already prewhitened
    X_pw <- X
  } else if (!is.null(obj$ar_order) && obj$ar_order > 0) {
    # Estimate per-fold
    ar_fit <- fit_ar_model(X, obj$ar_order)
    X_pw <- apply_ar_whitening(X, ar_fit)
  } else {
    X_pw <- X
    ar_fit <- NULL
  }
  
  fit <- hrfdecoder::hrfdecoder_fit(X_pw, event_model=obj$design$event_model, ...)
  
  structure(
    list(fit=fit, ar_fit=ar_fit, sl_info=sl_info, indices=indices),
    class="hrfdecoder_fit_wrap"
  )
}

# Inside format_result:
format_result.hrfdecoder_model <- function(obj, result, error_message=NULL, context, ...) {
  Xtest <- as.matrix(context$test)
  
  # Apply same AR whitening to test data
  if (!is.null(result$ar_fit)) {
    Xtest_pw <- apply_ar_whitening(Xtest, result$ar_fit)
  } else {
    Xtest_pw <- Xtest
  }
  
  Ptest <- predict(result$fit, Xtest_pw)
  # ... aggregation and result formatting
}
```

This architecture ensures: 1. **Train/test separation**: Each fold’s AR
model estimated from training data only 2. **Consistency**: Test data
transformed using training fold’s AR parameters 3. **Model
transparency**: AR fit stored in result object for inspection 4.
**Flexibility**: Works with or without prewhitening 5.
**Extensibility**: Same pattern works for other preprocessing (nuisance
regression, etc.)
