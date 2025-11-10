# rMVPA Architecture: Key Findings for hrfdecoder Integration

## Quick File Reference

### Core Plugin Architecture Files

1. **`/R/allgeneric.R` (944 lines)**
   - Generic functions: `train_model()`, `predict_model()`, `format_result()`, `merge_results()`, `compute_performance()`
   - Key class resolution dispatcher
   - `process_roi.default()` entry point (line 206-222) - controls flow to CV methods

2. **`/R/mvpa_model.R` (387 lines)**
   - `mvpa_model()` (line 260-337) - creates standard model spec
   - `create_model_spec()` (line 196-213) - internal spec factory
   - S3 methods for mvpa_model class
   - Shows how to attach metadata to model specs

3. **`/R/hrfdecoder_model.R` (336 lines)** - Already exists!
   - `hrfdecoder_model()` (line 97-146) - creates hrfdecoder spec
   - `train_model.hrfdecoder_model()` (line 175-210) - trains decoder per fold
   - `format_result.hrfdecoder_model()` (line 216-297) - predicts + aggregates to events
   - `merge_results.hrfdecoder_model()` (line 303-328) - folds → classification_result
   - `compute_performance.hrfdecoder_model()` (line 334-336) - delegates to perf_fun

4. **`/R/hrfdecoder_design.R` (174 lines)** - Already exists!
   - `hrfdecoder_design()` (line 60-130) - extends mvpa_design with event metadata
   - Validates events, sampling_frame, TR alignment
   - Stores `event_model` and `events` for training

### Data Flow Architecture

5. **`/R/mvpa_iterate.R` (585 lines)**
   - `internal_crossval()` (line 224-296) - MAIN CV LOOP
     - Line 228-229: Generates folds via `crossval_samples()`
     - Line 263-291: Per-fold training via `train_model()` + `format_result()`
     - Line 295: Merges via `merge_results()`
   - `external_crossval()` (line 89-175) - for external test sets
   - Key insight: Data extraction (line 228) happens BEFORE fold creation
     - **This is where prewhitening must hook in (or inside `train_model`)**

6. **`/R/crossval.R` (903 lines)**
   - `blocked_cross_validation()` (line 402-406)
   - `crossval_samples()` generics (line 867)
   - `crossval_samples.blocked_cross_validation()` (line 572-574)
   - Implementation: `crossv_block()` (line 170-192)

### Preprocessing Architecture

7. **`/R/common.R` (624 lines)**
   - `normalize_image_samples()` (line 52-60) - per-volume z-scoring
   - `standardize_vars()` (line 66-77) - per-block z-scoring
   - `normalize_surface_samples()` (line 82-89)
   - **IMPORTANT**: These are dataset-level (not fold-aware)
   - **No AR/nuisance regression currently**

### Dataset & Design Structures

8. **`/R/dataset.R` (479 lines)**
   - `mvpa_dataset()` (line 193+) - creates dataset spec
   - `gen_sample_dataset()` (line 65-156) - test data generation
   - Structure: `dataset$train_data`, `dataset$test_data`, `dataset$mask`

9. **`/R/design.R` (344 lines)**
   - `mvpa_design()` - creates design spec
   - `y_train()`, `y_test()` generics (line 366-383)

### Model Registry

10. **`/R/classifiers.R` (739 lines)**
    - `MVPAModels` environment - registry of available models
    - Example: `MVPAModels$corclass` (line 134-158)
    - Each has `fit()`, `predict()`, `prob()` methods

---

## Critical Integration Points for AR Prewhitening

### Where Preprocessing Can Hook In

```
OPTION A: Per-ROI, Before Fold Creation (PREFERRED)
  Location: process_roi() → [preprocess_roi() HOOK] → internal_crossval()
  
  Pros:
  ✓ Whole-ROI AR estimation (more stable)
  ✓ Easier to manage parameters
  ✗ Must handle fold-level AR application in train_model
  
OPTION B: Per-Fold, Inside train_model()
  Location: train_model.hrfdecoder_model() → [estimate AR] → train/test
  
  Pros:
  ✓ Fold-specific AR parameters (stricter separation)
  ✓ Simpler logic
  ✗ Repeated AR estimation per fold (slower, harder tuning)
```

### Key Code Locations for Integration

**Entry point modification** (would add hook):
```
/R/allgeneric.R:206-222 (process_roi.default)
  Add: roi <- preprocess_roi(mod_spec, roi)
```

**AR whitening storage**:
```
/R/hrfdecoder_model.R:209 (train_model returns fit object)
  Store: AR params in fit object for test data application
```

**Test data application**:
```
/R/hrfdecoder_model.R:237-243 (format_result.hrfdecoder_model)
  Apply: Same AR params to test data before prediction
```

---

## Data Structure Flow for hrfdecoder

```
hrfdecoder_design()
  ├─ train_design = data.frame(y=1:T, block=block_var)
  ├─ y_train = ~ y  (dummy, just for length)
  ├─ block_var = ~ block  (DETERMINES folds)
  ├─ event_model  (from fmridesign)
  └─ events  (data.frame with onset, condition)

mvpa_dataset()
  ├─ train_data = NeuroVec(T x V)  [TR-level]
  ├─ test_data = NeuroVec(T_test x V)
  └─ mask = NeuroVol

hrfdecoder_model()
  ├─ lambda_W, lambda_HRF, lambda_smooth
  ├─ basis = HRF (fmrihrf)
  ├─ window = event aggregation
  └─ crossval = blocked_cross_validation(block_var)

run_searchlight(spec) → process_roi() → internal_crossval()
  ├─ [EXTRACT ROI: TR x voxels]
  ├─ [PREPROCESS: AR whitening???]
  ├─ [CREATE FOLDS by block_var]
  ├─ For each fold:
  │  ├─ train_model() → hrfdecoder_fit(X_train, event_model)
  │  └─ format_result() → predict(X_test) → aggregate_events()
  └─ merge_results() → classification_result with event-level probs
```

---

## Key Method Signatures

### For implementing train_model with preprocessing:

```r
train_model.hrfdecoder_model <- function(obj, train_dat, y, sl_info, cv_spec, indices, ...) {
  X <- as.matrix(train_dat)
  
  # PREWHITEN training data
  X_pw <- apply_ar(X, obj$ar_order)  # Returns: (X_pw, ar_params)
  ar_params <- attr(X_pw, "ar_params")
  
  # FIT on prewhitened data
  fit <- hrfdecoder::hrfdecoder_fit(X_pw, event_model=obj$design$event_model, ...)
  
  # RETURN with AR params for test data
  structure(
    list(fit=fit, ar_params=ar_params, sl_info=sl_info, indices=indices),
    class="hrfdecoder_fit_wrap"
  )
}

format_result.hrfdecoder_model <- function(obj, result, error_message=NULL, context, ...) {
  Xtest <- as.matrix(context$test)
  
  # APPLY SAME AR to test data using training fold's params
  Xtest_pw <- apply_ar_inverse(Xtest, result$ar_params)
  
  # PREDICT on prewhitened test data
  Ptest <- predict(result$fit, Xtest_pw)
  
  # AGGREGATE to event level
  agg <- hrfdecoder::aggregate_events(Ptest, obj$design$events, window=obj$window)
  
  # RETURN results tibble
  tibble::tibble(
    class = list(pred_class),
    probs = list(probs),
    y_true = list(observed),
    test_ind = list(as.integer(context$test)),
    fit = list(if (obj$return_fits) result$fit else NULL),
    error = FALSE,
    error_message = "~"
  )
}
```

---

## Summary: What's Already in Place vs. What's Needed

### Already Implemented ✓
- [x] hrfdecoder_model() spec creation
- [x] train_model.hrfdecoder_model() solver integration
- [x] format_result.hrfdecoder_model() event aggregation
- [x] merge_results.hrfdecoder_model() fold merging
- [x] hrfdecoder_design() for metadata
- [x] Blocked cross-validation (respects TR structure)
- [x] S3 method dispatch system

### Needs Implementation for AR Prewhitening ✗
- [ ] AR model fitting (per-voxel or global)
- [ ] AR whitening application (forward)
- [ ] AR inverse application (for test data)
- [ ] Storage in fit object (ar_params)
- [ ] Integration point in process_roi() or train_model()
- [ ] Parameter control (ar_order) in hrfdecoder_model()
- [ ] Validation tests (train/test separation)

---

## Files Used in This Analysis

Absolute paths in rMVPA package:

- `/Users/bbuchsbaum/code/rMVPA/R/allgeneric.R` - 944 lines
- `/Users/bbuchsbaum/code/rMVPA/R/mvpa_model.R` - 387 lines
- `/Users/bbuchsbaum/code/rMVPA/R/hrfdecoder_model.R` - 336 lines
- `/Users/bbuchsbaum/code/rMVPA/R/hrfdecoder_design.R` - 174 lines
- `/Users/bbuchsbaum/code/rMVPA/R/mvpa_iterate.R` - 585 lines
- `/Users/bbuchsbaum/code/rMVPA/R/crossval.R` - 903 lines
- `/Users/bbuchsbaum/code/rMVPA/R/common.R` - 624 lines
- `/Users/bbuchsbaum/code/rMVPA/R/dataset.R` - 479 lines
- `/Users/bbuchsbaum/code/rMVPA/R/design.R` - 344 lines
- `/Users/bbuchsbaum/code/rMVPA/R/classifiers.R` - 739 lines
- `/Users/bbuchsbaum/code/rMVPA/R/searchlight.R` - 1130 lines
- `/Users/bbuchsbaum/code/rMVPA/R/regional.R` - 549 lines

**Total analyzed**: ~7,000 lines of rMVPA code focused on architecture
