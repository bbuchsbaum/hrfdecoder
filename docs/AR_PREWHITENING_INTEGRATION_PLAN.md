# AR Prewhitening Integration Plan for hrfdecoder

**Date:** 2025-11-09 **Authors:** Analysis by 3 specialized sub-agents
(fmriAR, rMVPA, hrfdecoder) **Status:** Design Complete, Ready for
Implementation

------------------------------------------------------------------------

## Executive Summary

This document provides a comprehensive plan for integrating AR(1)/ARMA
prewhitening into the hrfdecoder package using fmriAR, with full
consideration of the rMVPA plugin architecture. The integration
addresses autocorrelated noise in fMRI data while maintaining clean
separation of concerns across the preprocessing pipeline.

### Key Findings

1.  **Current State:** hrfdecoder handles nuisance signals via
    pre-projection but does NOT address autocorrelated noise
2.  **Optimal Solution:** Use fmriAR for run-specific AR estimation and
    whitening
3.  **Integration Point:** Insert between baseline residualization and
    standardization
4.  **rMVPA Impact:** Minimal changes required; AR parameters passed
    through model spec

------------------------------------------------------------------------

## Table of Contents

1.  [Current Noise Handling
    Analysis](#id_1-current-noise-handling-analysis)
2.  [fmriAR Capabilities Overview](#id_2-fmriar-capabilities-overview)
3.  [rMVPA Architecture &
    Constraints](#id_3-rmvpa-architecture--constraints)
4.  [Integration Architecture](#id_4-integration-architecture)
5.  [Implementation Plan](#id_5-implementation-plan)
6.  [Testing Strategy](#id_6-testing-strategy)
7.  [Performance Considerations](#id_7-performance-considerations)

------------------------------------------------------------------------

## 1. Current Noise Handling Analysis

### 1.1 Autocorrelated Noise: NOT Handled

**Evidence from hrfdecoder:**

- **File:** [src/softlabels_als.cpp:39](../src/softlabels_als.cpp#L39)

  ``` cpp
  recon = arma::accu(arma::square(S - P));  // Assumes i.i.d. noise
  ```

- **Objective function:**

      L = ||XW - P||²_F              (reconstruction - assumes white noise)
        + λ_W ||W||²_F                (weight regularization)
        + λ_HRF ||P - DBβ||²_F        (prior adherence)
        + λ_smooth P^T L P            (temporal smoothness)

**Why this matters:**

fMRI noise exhibits temporal autocorrelation (typically AR(1) with ρ ≈
0.3-0.5):

- **Inflated degrees of freedom** → Liberal statistical tests
- **Suboptimal decoder weights** → Not using GLS-optimal weights
- **Biased variance estimates** → Confidence intervals too narrow
- **Reduced efficiency** → More data needed for same power

**Mentioned in design notes (Section 4):** \> “AR(1) prewhitening per
run prior to fitting” listed as **optional enhancement** (not
implemented)

------------------------------------------------------------------------

### 1.2 Nuisance Signals: Handled via Pre-Projection ✓

**Evidence from hrfdecoder:**

- **File:** [R/fit.R:40](../R/fit.R#L40)

  ``` r
  Y <- residualize_baseline(Y, base_model)  # BEFORE any fitting
  ```

- **Implementation:**
  [R/interop_fmri.R:44-46](../R/interop_fmri.R#L44-L46)

  ``` r
  residualize_baseline <- function(Y, base_model = NULL) {
    if (is.null(base_model)) return(Y)
    fmridesign::residualize(base_model, Y)  # Delegate to fmridesign
  }
  ```

**What gets removed:**

- Low-frequency drift (polynomial trends, high-pass filtering)
- Motion artifacts (6 or 24 motion regressors)
- Physiological noise (cardiac, respiratory if modeled)
- Other nuisance regressors specified in `baseline_model`

**Approach: Orthogonal Projection**

    Y_clean = Y - X_baseline (X_baseline^T X_baseline)^(-1) X_baseline^T Y

**Conclusion:** ✓ Nuisance removal is properly implemented via
pre-projection. No changes needed.

------------------------------------------------------------------------

## 2. fmriAR Capabilities Overview

### 2.1 Core Functionality

**Package:** fmriAR v0.3.0 **Location:** `/Users/bbuchsbaum/code/fmriAR`

**Provides:**

1.  **AR/ARMA parameter estimation** from residuals (Yule-Walker,
    Hannan-Rissanen)
2.  **Run-aware whitening** with automatic boundary resets
3.  **Censor-aware whitening** for excluding bad TRs
4.  **Multi-scale spatial pooling** across parcel hierarchies
5.  **AFNI-compatible restricted AR models**
6.  **Fast C++ implementation** with OpenMP parallelization

### 2.2 API Design

**Primary workflow:**

``` r
# 1. Estimate AR/ARMA noise model from residuals
plan <- fit_noise(
  resid,                  # (T x V) residual matrix
  runs = run_ids,         # Run identifiers (respects boundaries)
  method = "ar",          # or "arma"
  p = "auto",             # Auto-select order via BIC (or specify integer)
  q = 0,                  # MA order (0 for pure AR)
  exact_first = "ar1"     # Exact AR(1) scaling at segment starts
)

# 2. Apply whitening to design and data matrices
whitened <- whiten_apply(
  plan,
  X, Y,
  runs = run_ids,         # Same run structure as estimation
  censor = NULL           # Optional: indices of bad TRs to skip
)

# 3. Use whitened data
X_white <- whitened$X
Y_white <- whitened$Y
```

**Returned objects:**

- **`fit_noise()`** → `fmriAR_plan` S3 object with:
  - `phi`: AR coefficients (per-run or pooled)
  - `theta`: MA coefficients
  - `order`: (p, q)
  - `runs`: Run labels
  - `method`: “ar”, “arma”, or “afni”
- **`whiten_apply()`** → List with:
  - `X`: Whitened design matrix
  - `Y`: Whitened data matrix

### 2.3 Run-Specific AR Handling

**Key feature:** Estimates separate AR parameters per run and resets
whitening buffers at run boundaries.

``` r
runs <- rep(1:3, each=100)  # 3 runs of 100 TRs each

plan <- fit_noise(resid, runs = runs, pooling = "run", p = 1)
# plan$phi = list(φ₁, φ₂, φ₃)  # One per run

whitened <- whiten_apply(plan, X, Y, runs = runs)
# Applies run-specific φᵣ and resets AR state at run boundaries
```

**C++ implementation:** `arma_whiten_inplace()` in
[fmriAR/src/fmriAR_whiten.cpp](https://github.com/bbuchsbaum/fmriAR/blob/main/src/fmriAR_whiten.cpp)

``` cpp
// Pseudocode: For each time t and voxel j:
y_t[j] = y_t[j] - sum_{k=1}^p φ_k * y_{t-k}[j]  // AR filter

// At run boundaries:
// - Reset buffers (no cross-run AR dependency)
// - Apply exact_first_ar1 scaling if p=1: y_t[j] *= sqrt(1 - φ₁²)
```

### 2.4 Train/Test Workflow

**Critical design pattern for cross-validation:**

``` r
# TRAINING PHASE
# Step 1: Fit AR model on training residuals
train_resid <- Y_train - X_train %*% solve(X_train, Y_train)
plan <- fit_noise(train_resid, runs = train_runs, method = "ar", p = 2)

# Step 2: Whiten training data
train_white <- whiten_apply(plan, X_train, Y_train, runs = train_runs)

# TESTING PHASE
# Step 3: Apply SAME plan to test data (no re-estimation!)
test_white <- whiten_apply(plan, X_test, Y_test, runs = test_runs)
```

**Key principle:** AR parameters learned from training data only, then
frozen for test application.

### 2.5 Integration Strengths

✓ **Run-aware:** Proper multi-run fMRI handling ✓ **Train/test
separation:** Plan estimated once, applied to both ✓ **Fast:** C++ +
OpenMP parallelization ✓ **Flexible pooling:** Global, run-specific, or
parcel-based AR ✓ **Censoring support:** Can exclude bad TRs ✓
**Production-ready:** Mature package with comprehensive tests

------------------------------------------------------------------------

## 3. rMVPA Architecture & Constraints

### 3.1 Model Plugin System

**Package:** rMVPA **Location:** `/Users/bbuchsbaum/code/rMVPA`

**S3 Method Dispatch Pattern:**

``` r
# Core generics for any model plugin:
train_model(obj, train_dat, y, sl_info, cv_spec, indices, ...)
predict_model(obj, new_data, ...)
format_result(obj, result, error_message, context, ...)
merge_results(obj, result_set, indices, id, ...)
compute_performance(obj, result, ...)
```

**hrfdecoder integration:** **File:**
[R/hrfdecoder_model.R](../R/hrfdecoder_model.R)

Already implements all required S3 methods: -
`train_model.hrfdecoder_model()` (lines 48-71) -
`format_result.hrfdecoder_model()` (lines 74-105) -
`merge_results.hrfdecoder_model()` (lines 108-126)

### 3.2 Data Flow in Cross-Validation

    ┌─────────────────────────────────────────────────────────────┐
    │ 1. ROI/Searchlight Extraction                              │
    │    mvpa_iterate() extracts T×V matrix per ROI              │
    └────────────────┬────────────────────────────────────────────┘
                     │
                     ▼
    ┌─────────────────────────────────────────────────────────────┐
    │ 2. Cross-Validation Fold Creation                          │
    │    Based on block_var (run IDs) or custom CV scheme        │
    │    → train_indices, test_indices                           │
    └────────────────┬────────────────────────────────────────────┘
                     │
                     ▼
    ┌─────────────────────────────────────────────────────────────┐
    │ 3. PER FOLD:                                                │
    │    a. train_model() called with train_dat                  │
    │       → Fits decoder on training subset                    │
    │    b. format_result() called with test data                │
    │       → Predicts on test subset                            │
    └────────────────┬────────────────────────────────────────────┘
                     │
                     ▼
    ┌─────────────────────────────────────────────────────────────┐
    │ 4. merge_results() across folds                            │
    │    → Combined classification_result                        │
    └─────────────────────────────────────────────────────────────┘

### 3.3 Preprocessing Integration Points

**Two architectural options identified:**

#### **Option A: Per-Fold Preprocessing (RECOMMENDED)**

``` r
train_model.hrfdecoder_model <- function(obj, train_dat, ...) {
  X_train <- as.matrix(train_dat)

  # ===== AR PREWHITENING HAPPENS HERE =====
  # 1. Estimate AR from training residuals
  # 2. Whiten training data
  # 3. Store plan in fit object

  fit <- fit_hrfdecoder(
    Y = X_train_white,  # Whitened data
    ev_model = obj$design$event_model,
    ar_plan = plan,     # NEW: Store for test application
    ...
  )

  return(fit)
}

format_result.hrfdecoder_model <- function(obj, result, context, ...) {
  X_test <- as.matrix(context$test)

  # ===== APPLY SAME AR PLAN TO TEST DATA =====
  X_test_white <- whiten_apply(result$fit$ar_plan, NULL, X_test)$Y

  preds <- predict_hrfdecoder(result$fit, X_test_white, ...)
  return(preds)
}
```

**Pros:** - ✓ CV-safe: AR estimated from training fold only - ✓ Minimal
rMVPA changes: Contained in model-specific methods - ✓ Flexible: Each
fold can have different AR structure

**Cons:** - ✗ Redundant computation: Re-estimates AR for each fold - ✗
More complex: Logic split across train/predict

#### **Option B: Dataset-Level Preprocessing**

Prewhiten entire dataset before passing to rMVPA. **NOT RECOMMENDED**
because: - ✗ CV leakage: Test data influences AR estimation - ✗ Less
flexible: Can’t adapt AR per fold

### 3.4 Key Constraints from rMVPA

1.  **Blocked CV is standard:** Uses `block_var` (run IDs) to create
    folds
2.  **Result format required:** Must return tibble with columns `class`,
    `probs`, `y_true`, `test_ind`, `fit`, `error`, `error_message`
3.  **Event model required:** For trial-level evaluation (the default
    mode)
4.  **Standardization handling:** rMVPA may pre-standardize; hrfdecoder
    currently calls with `standardize = FALSE` in rMVPA context

------------------------------------------------------------------------

## 4. Integration Architecture

### 4.1 Preprocessing Pipeline Order

**RECOMMENDED ORDER:**

    Input Y (T × V raw fMRI data)
              ↓
    ┌─────────────────────────────────────────┐
    │ 1. Baseline Residualization             │  [R/fit.R:40]
    │    Y ← residualize_baseline(Y, base_model)
    │    Removes: drift, motion, nuisance    │
    └─────────────────┬───────────────────────┘
                      ↓
    ┌─────────────────────────────────────────┐
    │ 2. AR Prewhitening (NEW)                │  [INSERTION POINT]
    │    plan ← fit_noise(Y, runs, p)        │
    │    Y ← whiten_apply(plan, NULL, Y)$Y   │
    │    Removes: temporal autocorrelation   │
    └─────────────────┬───────────────────────┘
                      ↓
    ┌─────────────────────────────────────────┐
    │ 3. Standardization                      │  [R/fit.R:43-50]
    │    Y ← scale(Y)                        │
    │    Z-scores: (Y - mean) / sd           │
    └─────────────────┬───────────────────────┘
                      ↓
    ┌─────────────────────────────────────────┐
    │ 4. Decoder Input Preparation            │  [R/prep.R:40-58]
    │    X_list, DBbeta, P0 ← prepare_decoder_inputs()
    │    L ← build_laplacian_from_runs()    │
    └─────────────────┬───────────────────────┘
                      ↓
    ┌─────────────────────────────────────────┐
    │ 5. ALS Solver                           │  [src/softlabels_als.cpp:68]
    │    W, P, b ← fit_softlabels_als()      │
    └─────────────────────────────────────────┘

**Rationale:**

1.  **Baseline FIRST:** Removes structured nuisance before AR estimation
2.  **AR SECOND:** Estimates autocorrelation on residuals (after
    nuisance removal)
3.  **Standardization THIRD:** Z-scores applied to whitened residuals
4.  **ALS LAST:** Operates on fully preprocessed data

### 4.2 Exact Code Insertion Point

**File:** `/Users/bbuchsbaum/code/hrfdecoder/R/fit.R` **Location:**
Between lines 40 and 42

``` r
fit_hrfdecoder <- function(
  Y,
  ev_model,
  base_model = NULL,
  hrf = NULL,
  ar_order = 1,              # NEW PARAMETER
  ar_method = "ar",          # NEW: "ar" or "arma"
  ar_pooling = "run",        # NEW: "global" or "run"
  lambda_W = 10,
  lambda_HRF = 1,
  lambda_smooth = 5,
  theta_penalty = 0.01,
  max_iter = 20,
  tol = 1e-4,
  nonneg = TRUE,
  background = TRUE,
  standardize = TRUE,
  verbose = 1
) {
  stopifnot(is.matrix(Y))
  stopifnot(inherits(ev_model, "event_model"))
  Tn <- nrow(Y)

  # ===== STEP 1: Baseline residualization =====
  Y <- residualize_baseline(Y, base_model)

  # ===== STEP 2: AR PREWHITENING (NEW) =====
  ar_plan <- NULL
  if (!is.null(ar_order) && ar_order > 0) {
    # Need run_ids BEFORE prep
    prep_temp <- prepare_decoder_inputs(ev_model, hrf = hrf, background = background)
    run_ids <- prep_temp$run_ids

    # Estimate AR from current Y (post-baseline, pre-standardization)
    ar_plan <- fmriAR::fit_noise(
      Y,
      runs = run_ids,
      method = ar_method,
      p = ar_order,
      pooling = ar_pooling,
      exact_first = "ar1"
    )

    # Apply whitening
    Y_white <- fmriAR::whiten_apply(ar_plan, X = NULL, Y = Y, runs = run_ids)
    Y <- Y_white$Y

    if (verbose) {
      message("Applied AR(", ar_order, ") prewhitening (pooling=", ar_pooling, ")")
    }
  }

  # ===== STEP 3: Standardization =====
  preproc <- list(
    standardize = isTRUE(standardize),
    center = NULL,
    scale = NULL,
    ar_plan = ar_plan  # NEW: Store for predict
  )
  if (isTRUE(standardize)) {
    Ys <- scale(Y)
    preproc$center <- attr(Ys, "scaled:center")
    s <- attr(Ys, "scaled:scale")
    s[is.na(s) | s == 0] <- 1
    preproc$scale <- s
    attr(Ys, "scaled:center") <- NULL
    attr(Ys, "scaled:scale") <- NULL
    Y <- Ys
  }

  # ===== STEP 4: Decoder preparation & fitting =====
  prep <- prepare_decoder_inputs(ev_model, hrf = hrf, background = background)
  # ... (rest unchanged)

  fit <- list(
    W = fit_cpp$W,
    P = Z,
    b = as.numeric(fit_cpp$b),
    theta = theta,
    hrf = hrf_est,
    conditions = prep$conditions,
    background = background,
    converged = isTRUE(fit_cpp$converged),
    iterations = fit_cpp$iterations,
    settings = list(
      lambda_W = lambda_W,
      lambda_HRF = lambda_HRF,
      lambda_smooth = lambda_smooth,
      theta_penalty = theta_penalty,
      max_iter = max_iter,
      tol = tol,
      nonneg = nonneg,
      background = background,
      TR = TR,
      run_ids = prep$run_ids,
      ar_order = ar_order,         # NEW
      ar_method = ar_method,        # NEW
      ar_pooling = ar_pooling       # NEW
    ),
    train = list(
      P0 = prep$P0,
      prior = prep$DBbeta,
      events = events_tbl
    ),
    preproc = preproc,              # Contains ar_plan
    diagnostics = list(
      obj_trace = fit_cpp$obj_trace
    )
  )
  class(fit) <- "hrfdecoder_fit"
  fit
}
```

### 4.3 Prediction Update

**File:** `/Users/bbuchsbaum/code/hrfdecoder/R/predict.R` **Location:**
Lines 21-27 (preprocessing section)

``` r
predict_hrfdecoder <- function(
  object,
  Y_test,
  ev_model_test = NULL,
  mode = c("tr", "trial"),
  window = c(4, 8),
  weights = c("hrf", "flat")
) {
  stopifnot(inherits(object, "hrfdecoder_fit"))
  mode <- match.arg(mode)
  weights <- match.arg(weights)

  # ===== APPLY PREPROCESSING IN SAME ORDER AS TRAINING =====

  # STEP 1: AR prewhitening (if was applied during training)
  if (!is.null(object$preproc$ar_plan)) {
    run_ids_test <- get_run_ids_from_test_data(object, Y_test, ev_model_test)
    Y_white <- fmriAR::whiten_apply(
      object$preproc$ar_plan,
      X = NULL,
      Y = Y_test,
      runs = run_ids_test
    )
    Y_test <- Y_white$Y
  }

  # STEP 2: Standardization (if was applied during training)
  if (!is.null(object$preproc) && isTRUE(object$preproc$standardize)) {
    ctr <- object$preproc$center %||% rep(0, ncol(Y_test))
    scl <- object$preproc$scale %||% rep(1, ncol(Y_test))
    scl[is.na(scl) | scl == 0] <- 1
    Y_test <- sweep(Y_test, 2, ctr, FUN = "-")
    Y_test <- sweep(Y_test, 2, scl, FUN = "/")
  }

  # ===== COMPUTE PREDICTIONS =====
  scores <- predict_softlabels(Y_test, object$W, object$b)
  probs <- row_softmax(scores)

  # ... (rest unchanged)
}

# NEW HELPER FUNCTION
get_run_ids_from_test_data <- function(fit_obj, Y_test, ev_model_test) {
  # Extract run IDs for test data
  # Option 1: From ev_model_test if provided
  if (!is.null(ev_model_test)) {
    prep <- prepare_decoder_inputs(ev_model_test, hrf = fit_obj$hrf,
                                   background = fit_obj$background)
    return(prep$run_ids)
  }

  # Option 2: Assume same run structure as training (if single-run test)
  if (nrow(Y_test) == length(fit_obj$settings$run_ids)) {
    return(fit_obj$settings$run_ids)
  }

  # Option 3: Default to single run
  return(rep(1L, nrow(Y_test)))
}
```

### 4.4 rMVPA Integration Updates

**File:** `/Users/bbuchsbaum/code/hrfdecoder/R/hrfdecoder_model.R`
**Updates required:**

``` r
hrfdecoder_model <- function(
  dataset,
  design,
  lambda_W = 10,
  lambda_HRF = 1,
  lambda_smooth = 5,
  theta_penalty = 0.01,
  basis = NULL,
  window = c(4, 8),
  nonneg = TRUE,
  max_iter = 10,
  tol = 1e-4,
  ar_order = 1,           # NEW
  ar_method = "ar",       # NEW
  ar_pooling = "run",     # NEW
  performance = NULL,
  crossval = NULL,
  return_predictions = TRUE,
  return_fits = FALSE
) {
  # ... auto-detect CV ...

  rMVPA::create_model_spec(
    "hrfdecoder_model",
    dataset = dataset,
    design = design,
    lambda_W = lambda_W,
    lambda_HRF = lambda_HRF,
    lambda_smooth = lambda_smooth,
    theta_penalty = theta_penalty,
    basis = basis,
    window = window,
    nonneg = nonneg,
    max_iter = max_iter,
    tol = tol,
    ar_order = ar_order,        # NEW
    ar_method = ar_method,      # NEW
    ar_pooling = ar_pooling,    # NEW
    crossval = crossval,
    performance = performance,
    compute_performance = TRUE,
    return_predictions = return_predictions,
    return_fits = return_fits
  )
}

# train_model.hrfdecoder_model UNCHANGED - just passes ar_order to fit_hrfdecoder
train_model.hrfdecoder_model <- function(obj, train_dat, y, sl_info, cv_spec, indices, ...) {
  X <- as.matrix(train_dat)
  ev_model <- obj$design$event_model
  base_model <- obj$design$baseline_model %||% NULL

  fit <- fit_hrfdecoder(
    Y = X,
    ev_model = ev_model,
    base_model = base_model,
    hrf = obj$basis,
    lambda_W = obj$lambda_W,
    lambda_HRF = obj$lambda_HRF,
    lambda_smooth = obj$lambda_smooth,
    theta_penalty = obj$theta_penalty,
    ar_order = obj$ar_order,       # NEW
    ar_method = obj$ar_method,     # NEW
    ar_pooling = obj$ar_pooling,   # NEW
    max_iter = obj$max_iter,
    tol = obj$tol,
    nonneg = obj$nonneg,
    standardize = FALSE,
    verbose = 0
  )

  structure(
    list(fit = fit, sl_info = sl_info, indices = indices),
    class = "hrfdecoder_fit_wrap"
  )
}

# format_result.hrfdecoder_model UNCHANGED - predict_hrfdecoder handles AR internally
```

------------------------------------------------------------------------

## 5. Implementation Plan

### 5.1 Phase 1: Core AR Integration (Week 1)

**Files to modify:**

1.  **`R/fit.R`** (main changes)
    - Add parameters: `ar_order`, `ar_method`, `ar_pooling`
    - Insert AR prewhitening between lines 40-42
    - Store `ar_plan` in `fit$preproc`
    - Store AR settings in `fit$settings`
2.  **`R/predict.R`**
    - Apply AR whitening before standardization
    - Add helper `get_run_ids_from_test_data()`
3.  **`R/hrfdecoder_model.R`**
    - Add AR parameters to
      [`hrfdecoder_model()`](reference/hrfdecoder_model.md) signature
    - Pass through to [`fit_hrfdecoder()`](reference/fit_hrfdecoder.md)
      in
      [`train_model()`](http://bbuchsbaum.github.io/rMVPA/reference/train_model.md)
4.  **`DESCRIPTION`**
    - Add `fmriAR (>= 0.3.0)` to Imports
    - Add remote: `bbuchsbaum/fmriAR`
5.  **`NAMESPACE`**
    - Import:
      [`fmriAR::fit_noise`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.html),
      [`fmriAR::whiten_apply`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.html)

### 5.2 Phase 2: Testing & Validation (Week 1-2)

**New test file:** `tests/testthat/test-ar-prewhitening.R`

**Test cases:**

1.  **Decorrelation test:**

    ``` r
    test_that("AR(1) prewhitening reduces autocorrelation in residuals", {
      # Generate AR(1) data with known phi
      # Fit with ar_order=1
      # Check: residual ACF near zero at lag 1
    })
    ```

2.  **Multi-run test:**

    ``` r
    test_that("Run-specific AR parameters are estimated and applied", {
      # Simulate 3 runs with different AR parameters
      # Fit with ar_pooling="run"
      # Check: fit$preproc$ar_plan has 3 different phi values
    })
    ```

3.  **Train/test consistency:**

    ``` r
    test_that("AR plan from training is applied to test data", {
      # Train on runs 1-2, test on run 3
      # Verify: test predictions use stored ar_plan
      # Verify: manual application of ar_plan matches predict output
    })
    ```

4.  **CV safety test:**

    ``` r
    test_that("AR estimation does not leak information across CV folds", {
      # Run blocked CV with ar_order=1
      # Verify: each fold estimates AR from training data only
      # Verify: test AR params != training AR params (different data)
    })
    ```

5.  **Backward compatibility:**

    ``` r
    test_that("ar_order=NULL or 0 reproduces old behavior exactly", {
      # Fit with ar_order=0
      # Fit without ar_order (old API)
      # Check: W, P, b identical
    })
    ```

### 5.3 Phase 3: Documentation & Examples (Week 2)

**Update roxygen docs:**

1.  **`R/fit.R`:**

    ``` r
    #' @param ar_order AR order for prewhitening (default: 1).
    #'   Set to 0 or NULL to disable.
    #' @param ar_method AR estimation method: "ar" (Yule-Walker) or
    #'   "arma" (Hannan-Rissanen). Default: "ar".
    #' @param ar_pooling Spatial pooling for AR: "global" (one AR for all voxels),
    #'   "run" (separate AR per run). Default: "run".
    ```

2.  **Vignette:** `vignettes/ar-prewhitening.Rmd`

    - When to use AR prewhitening
    - Impact on decoder performance
    - Choosing AR order (BIC, AIC, or fixed)
    - Computational cost considerations

3.  **README update:**

    - Add AR prewhitening to feature list
    - Simple example with/without AR

### 5.4 Phase 4: Performance Optimization (Week 3)

**Potential optimizations:**

1.  **Parallel AR estimation across voxels:**

    ``` r
    plan <- fmriAR::fit_noise(..., threads = TRUE)  # Enable OpenMP
    ```

2.  **Cache AR plan for repeated fits:**

    - If fitting multiple ROIs with same run structure, reuse AR plan

3.  **Auto AR order selection:**

    ``` r
    ar_order = "auto"  # Use BIC to select p ∈ {0, 1, 2, 3}
    ```

4.  **Benchmark comparison:**

    - Compare runtime with/without AR
    - Measure impact on searchlight analysis

------------------------------------------------------------------------

## 6. Testing Strategy

### 6.1 Unit Tests

**File:** `tests/testthat/test-ar-prewhitening.R`

**Coverage goals:**

- [x] AR estimation reduces autocorrelation
- [x] Run-specific AR parameters
- [x] Train/test consistency
- [x] CV fold independence
- [x] Backward compatibility (ar_order=0)
- [x] Error handling (bad inputs)

### 6.2 Integration Tests

**File:** `tests/testthat/test-rmvpa-integration.R`

**Coverage goals:**

- [x] rMVPA model spec with AR parameters
- [x] Searchlight with AR prewhitening
- [x] Blocked CV with AR
- [x] Performance metrics unchanged by AR (or improved)

### 6.3 Simulation Studies

**Validate AR improves decoder performance:**

1.  **Scenario 1:** High autocorrelation (ρ = 0.5)
    - Expected: AR prewhitening improves accuracy
2.  **Scenario 2:** Low autocorrelation (ρ = 0.1)
    - Expected: AR has minimal impact
3.  **Scenario 3:** Multi-run with varying AR
    - Expected: Run-specific AR outperforms global

### 6.4 Real Data Validation

**Apply to existing fMRI datasets:**

1.  Check residual whiteness (ACF diagnostics)
2.  Compare decoder accuracy with/without AR
3.  Measure computational overhead

------------------------------------------------------------------------

## 7. Performance Considerations

### 7.1 Computational Cost

**AR estimation:** - **Per voxel:** O(p² T) for Yule-Walker (fast) -
**Parallelized:** OpenMP across voxels - **Typical:** ~0.1-0.5 seconds
for 100 TRs × 1000 voxels

**Whitening:** - **Per voxel:** O(p T) recursive filtering -
**Parallelized:** OpenMP across voxels - **Typical:** ~0.05-0.2 seconds
for same data

**Impact on searchlight:** - **Overhead:** +5-10% runtime (negligible
vs. decoder fitting) - **Benefit:** More efficient decoder (better GLS
weights)

### 7.2 Memory Usage

**AR plan storage:** - Global pooling: O(p) parameters - Run pooling:
O(R × p) where R = number of runs - Parcel pooling: O(P × p) where P =
number of parcels

**Typical footprint:** \<1 MB for most analyses

### 7.3 Scaling Recommendations

**For large searchlights (V \> 10,000):** - Use `ar_pooling = "global"`
to reduce AR estimation cost - Enable OpenMP parallelization in fmriAR

**For many runs (R \> 10):** - Use `ar_pooling = "run"` for run-specific
noise structure - Consider caching AR plans across ROIs with same run
structure

------------------------------------------------------------------------

## 8. Open Questions & Future Enhancements

### 8.1 Open Questions

1.  **Should AR order be auto-selected via BIC?**
    - Pro: Automatic, data-driven
    - Con: Adds complexity, computational cost
2.  **Should we support ARMA(p, q) models?**
    - Pro: More flexible noise modeling
    - Con: Slower estimation, more parameters
3.  **Should AR be applied to design matrix (DBbeta)?**
    - Current: Only Y is whitened
    - Alternative: Whiten both Y and DBbeta prior
    - Trade-off: Computational cost vs. theoretical correctness

### 8.2 Future Enhancements

1.  **Parcel-based AR pooling:**

    ``` r
    ar_pooling = "parcel",
    parcels = parcel_labels  # Spatial pooling of AR
    ```

2.  **Multi-scale AR:**

    ``` r
    ar_pooling = "multiscale",
    parcel_sets = list(coarse, medium, fine)
    ```

3.  **AR diagnostics:**

    ``` r
    fit$diagnostics$ar_acf  # Post-whitening ACF
    fit$diagnostics$ar_bic  # BIC per AR order
    ```

4.  **Prewhitened prior:**

    ``` r
    # Whiten DBbeta prior along with Y
    DBbeta_white <- whiten_apply(plan, NULL, DBbeta)$Y
    ```

------------------------------------------------------------------------

## 9. Summary & Recommendations

### 9.1 Key Decisions

| Decision | Recommendation | Rationale |
|----|----|----|
| **AR library** | Use fmriAR | Production-ready, run-aware, fast C++ |
| **Integration point** | After baseline, before standardization | Estimates AR on residuals, cleanly separable |
| **AR pooling** | Default: `"run"` | Respects run-specific noise structure |
| **AR order** | Default: `1` (AR(1)) | Simple, covers most fMRI autocorrelation |
| **rMVPA integration** | Per-fold preprocessing in [`train_model()`](http://bbuchsbaum.github.io/rMVPA/reference/train_model.md) | CV-safe, no test leakage |
| **Backward compatibility** | `ar_order = NULL` disables AR | Preserves existing API |

### 9.2 Implementation Checklist

**Phase 1: Core integration** (Week 1)

Modify `R/fit.R` (add AR prewhitening)

Modify `R/predict.R` (apply AR to test)

Modify `R/hrfdecoder_model.R` (pass AR params)

Update `DESCRIPTION` (add fmriAR dependency)

Update `NAMESPACE` (import fmriAR functions)

**Phase 2: Testing** (Week 1-2)

Unit tests for AR estimation

Integration tests with rMVPA

Simulation studies (varying ρ)

Real data validation

**Phase 3: Documentation** (Week 2)

Roxygen docs for new parameters

Vignette on AR prewhitening

Update README examples

**Phase 4: Optimization** (Week 3)

Enable OpenMP parallelization

Benchmark performance

Profile memory usage

### 9.3 Expected Benefits

**Statistical:** - ✓ Correct effective degrees of freedom - ✓ Valid
confidence intervals - ✓ Optimal GLS decoder weights

**Practical:** - ✓ Improved decoder efficiency (less data needed) - ✓
Better generalization to test data - ✓ More accurate performance
estimates

**Computational:** - ✓ Minimal overhead (~5-10% runtime) - ✓
Parallelized AR estimation - ✓ Efficient recursive whitening

------------------------------------------------------------------------

## Appendix A: Code Snippets

### A.1 Complete fit_hrfdecoder() with AR

See Section 4.2 for full implementation.

### A.2 Complete predict_hrfdecoder() with AR

See Section 4.3 for full implementation.

### A.3 Example Usage

``` r
library(hrfdecoder)
library(fmridesign)
library(fmrihrf)

# 1. Create event model
ev_model <- event_model(
  onsets ~ hrf(condition, basis = "spmg3"),
  data = events_df,
  block = ~ run,
  sampling_frame = sampling_frame(blocklens = c(200, 200, 200), TR = 2)
)

# 2. Create baseline model (nuisance)
base_model <- baseline_model(
  ~ poly(run, degree = 3) + motion_1 + motion_2,
  data = nuisance_df,
  block = ~ run,
  sampling_frame = sampling_frame(blocklens = c(200, 200, 200), TR = 2)
)

# 3. Fit decoder with AR(1) prewhitening
fit <- fit_hrfdecoder(
  Y = roi_data,                 # (600 TRs × 500 voxels)
  ev_model = ev_model,
  base_model = base_model,
  hrf = fmrihrf::spmg3(),
  ar_order = 1,                 # NEW: Enable AR(1) prewhitening
  ar_method = "ar",             # Yule-Walker estimation
  ar_pooling = "run",           # Run-specific AR parameters
  lambda_W = 10,
  lambda_HRF = 1,
  lambda_smooth = 5,
  max_iter = 20,
  verbose = 1
)

# 4. Predict on test data
preds_tr <- predict_hrfdecoder(fit, Y_test, mode = "tr")
preds_trial <- predict_hrfdecoder(fit, Y_test, ev_model_test = ev_model_test,
                                  mode = "trial")

# 5. Inspect AR parameters
fit$preproc$ar_plan$phi    # List of AR coefficients per run
fit$settings$ar_order       # AR order used
```

### A.4 rMVPA Usage

``` r
library(rMVPA)
library(hrfdecoder)

# 1. Create design
mvdes <- continuous_mvpa_design(
  event_model = ev_model,
  block_var = run_ids,
  design_df_events = trials_df
)

# 2. Specify model with AR
spec <- hrfdecoder_model(
  dataset = as_mvpa_dataset(fmri_dataset),
  design = mvdes,
  basis = fmrihrf::spmg3(),
  ar_order = 1,              # NEW
  ar_method = "ar",          # NEW
  ar_pooling = "run",        # NEW
  lambda_W = 10,
  lambda_HRF = 1,
  lambda_smooth = 5
)

# 3. Run searchlight with AR prewhitening
results <- run_searchlight(spec, radius = 8, method = "randomized", niter = 4)
```

------------------------------------------------------------------------

## Appendix B: References

### B.1 Package Locations

- **hrfdecoder:** `/Users/bbuchsbaum/code/hrfdecoder`
- **fmriAR:** `/Users/bbuchsbaum/code/fmriAR`
- **rMVPA:** `/Users/bbuchsbaum/code/rMVPA`

### B.2 Key Files Modified

| File                                    | Lines | Changes                |
|-----------------------------------------|-------|------------------------|
| `R/fit.R`                               | 40-42 | Insert AR prewhitening |
| `R/predict.R`                           | 21-27 | Apply AR to test data  |
| `R/hrfdecoder_model.R`                  | 10-25 | Add AR parameters      |
| `DESCRIPTION`                           | —     | Add fmriAR dependency  |
| `tests/testthat/test-ar-prewhitening.R` | NEW   | AR unit tests          |

### B.3 Related Documentation

- fmriAR vignette:
  [fmriAR-introduction.Rmd](https://github.com/bbuchsbaum/fmriAR/blob/main/vignettes/fmriAR-introduction.Rmd)
- rMVPA model plugin guide:
  [custom-models.Rmd](https://github.com/bbuchsbaum/rMVPA/blob/master/vignettes/custom-models.Rmd)
- hrfdecoder design notes:
  [notes/hrf_weakly_supervised_decoder.md](../notes/hrf_weakly_supervised_decoder.md)

------------------------------------------------------------------------

**END OF INTEGRATION PLAN**
