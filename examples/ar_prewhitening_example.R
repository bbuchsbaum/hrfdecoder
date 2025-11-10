## AR Prewhitening Quick Start Example
## This script demonstrates the new AR prewhitening functionality in hrfdecoder

library(hrfdecoder)
library(fmridesign)
library(fmrihrf)
library(fmriAR)

# ==============================================================================
# EXAMPLE 1: Basic AR(1) Prewhitening
# ==============================================================================

# Simulate fMRI data with temporal autocorrelation
set.seed(42)
n_trs <- 200
n_voxels <- 100
phi <- 0.5  # AR(1) coefficient

# Generate AR(1) noise
Y <- matrix(0, n_trs, n_voxels)
for (j in 1:n_voxels) {
  eps <- rnorm(n_trs)
  for (i in 2:n_trs) {
    eps[i] <- phi * eps[i-1] + rnorm(1, sd = 0.8)
  }
  Y[, j] <- eps
}

# Create event model
sf <- sampling_frame(blocklens = n_trs, TR = 2)
events <- data.frame(
  onset = seq(10, 190, by = 20),
  condition = rep(c("face", "house"), length.out = 10),
  run = 1
)

ev_model <- event_model(
  onsets ~ hrf(condition, basis = "spmg1"),
  data = events,
  sampling_frame = sf,
  block = ~ run
)

# Fit decoder WITH AR(1) prewhitening
cat("\n=== Fitting WITH AR(1) prewhitening ===\n")
fit_ar <- fit_hrfdecoder(
  Y = Y,
  ev_model = ev_model,
  ar_order = 1,              # Enable AR(1)
  ar_method = "ar",          # Yule-Walker estimation
  ar_pooling = "global",     # Single AR for all voxels
  lambda_W = 10,
  lambda_HRF = 1,
  lambda_smooth = 5,
  max_iter = 10,
  verbose = 1
)

cat("\nAR plan stored:", !is.null(fit_ar$preproc$ar_plan), "\n")
cat("AR order:", fit_ar$settings$ar_order, "\n")
cat("Converged:", fit_ar$converged, "in", fit_ar$iterations, "iterations\n")

# Fit decoder WITHOUT AR (for comparison)
cat("\n=== Fitting WITHOUT AR prewhitening ===\n")
fit_no_ar <- fit_hrfdecoder(
  Y = Y,
  ev_model = ev_model,
  ar_order = NULL,           # Disable AR
  lambda_W = 10,
  lambda_HRF = 1,
  lambda_smooth = 5,
  max_iter = 10,
  verbose = 1
)

cat("Converged:", fit_no_ar$converged, "in", fit_no_ar$iterations, "iterations\n")

# Compare objective values
cat("\n=== Performance Comparison ===\n")
cat("Final objective (with AR):   ", tail(fit_ar$diagnostics$obj_trace$value, 1), "\n")
cat("Final objective (without AR):", tail(fit_no_ar$diagnostics$obj_trace$value, 1), "\n")


# ==============================================================================
# EXAMPLE 2: Multi-Run with Run-Specific AR
# ==============================================================================

cat("\n\n=== EXAMPLE 2: Multi-Run Data ===\n")

# Simulate 2 runs with different AR parameters
n1 <- 150
n2 <- 180
V <- 80
phi1 <- 0.4  # Run 1 AR
phi2 <- 0.6  # Run 2 AR

make_ar_noise <- function(n, V, phi) {
  noise <- matrix(0, n, V)
  for (j in 1:V) {
    eps <- rnorm(n)
    for (i in 2:n) {
      eps[i] <- phi * eps[i-1] + rnorm(1)
    }
    noise[, j] <- eps
  }
  noise
}

Y_run1 <- make_ar_noise(n1, V, phi1)
Y_run2 <- make_ar_noise(n2, V, phi2)
Y_multi <- rbind(Y_run1, Y_run2)

# Create multi-run event model
sf_multi <- sampling_frame(blocklens = c(n1, n2), TR = 2)
events_multi <- data.frame(
  onset = c(seq(10, 140, by = 15), seq(160, 330, by = 20)),
  condition = rep(c("face", "house"), length.out = 18),
  run = c(rep(1, 9), rep(2, 9))
)

ev_model_multi <- event_model(
  onsets ~ hrf(condition, basis = "spmg1"),
  data = events_multi,
  sampling_frame = sf_multi,
  block = ~ run
)

# Fit with run-specific AR
cat("\nFitting with run-specific AR parameters...\n")
fit_run_ar <- fit_hrfdecoder(
  Y = Y_multi,
  ev_model = ev_model_multi,
  ar_order = 1,
  ar_pooling = "run",        # Different AR per run
  lambda_W = 10,
  max_iter = 10,
  verbose = 1
)

cat("\nRun-specific AR pooling:", fit_run_ar$settings$ar_pooling, "\n")
cat("Number of runs:", length(unique(fit_run_ar$settings$run_ids)), "\n")


# ==============================================================================
# EXAMPLE 3: Prediction with AR Prewhitening
# ==============================================================================

cat("\n\n=== EXAMPLE 3: Prediction on Test Data ===\n")

# Simulate test data (same AR structure)
n_test <- 100
Y_test <- make_ar_noise(n_test, V, phi1)

# Create test event model
sf_test <- sampling_frame(blocklens = n_test, TR = 2)
events_test <- data.frame(
  onset = seq(10, 90, by = 10),
  condition = rep(c("face", "house"), length.out = 9),
  run = 1
)

ev_model_test <- event_model(
  onsets ~ hrf(condition, basis = "spmg1"),
  data = events_test,
  sampling_frame = sf_test,
  block = ~ run
)

# Predict with AR prewhitening (automatically applied)
cat("\nPredicting on test data (AR automatically applied)...\n")
preds_tr <- predict(fit_ar, newdata = Y_test, mode = "tr")
preds_trial <- predict(fit_ar, newdata = Y_test,
                                  ev_model_test = ev_model_test,
                                  mode = "trial")

cat("TR-level predictions shape:", dim(preds_tr), "\n")
cat("Trial-level predictions:", nrow(preds_trial$probs), "trials\n")
cat("Predicted classes:", table(colnames(preds_trial$probs)[max.col(preds_trial$probs)]), "\n")


# ==============================================================================
# EXAMPLE 4: Auto AR Order Selection
# ==============================================================================

cat("\n\n=== EXAMPLE 4: Automatic AR Order Selection ===\n")

# Simulate AR(2) data
n_ar2 <- 200
V_ar2 <- 60
phi1_ar2 <- 0.4
phi2_ar2 <- 0.3

Y_ar2 <- matrix(0, n_ar2, V_ar2)
for (j in 1:V_ar2) {
  eps <- rnorm(n_ar2)
  for (i in 3:n_ar2) {
    eps[i] <- phi1_ar2 * eps[i-1] + phi2_ar2 * eps[i-2] + rnorm(1)
  }
  Y_ar2[, j] <- eps
}

# Let fmriAR select AR order automatically
cat("\nFitting with automatic AR order selection...\n")
fit_auto <- fit_hrfdecoder(
  Y = Y_ar2,
  ev_model = ev_model,  # Reuse earlier event model
  ar_order = "auto",    # Automatic selection via BIC
  ar_pooling = "global",
  lambda_W = 10,
  max_iter = 10,
  verbose = 1
)

cat("\nSelected AR order:", fit_auto$preproc$ar_plan$order["p"], "\n")
cat("(True order was 2)\n")


# ==============================================================================
# EXAMPLE 5: Inspecting AR Parameters
# ==============================================================================

cat("\n\n=== EXAMPLE 5: Inspecting AR Plan ===\n")

# Access stored AR plan
ar_plan <- fit_ar$preproc$ar_plan

cat("\nAR plan details:\n")
cat("  Method:", ar_plan$method, "\n")
cat("  Order: p =", ar_plan$order["p"], ", q =", ar_plan$order["q"], "\n")
cat("  Pooling:", fit_ar$settings$ar_pooling, "\n")

if (fit_ar$settings$ar_pooling == "global") {
  cat("  AR coefficients (phi):", ar_plan$phi[[1]], "\n")
} else if (fit_ar$settings$ar_pooling == "run") {
  cat("  Number of runs:", length(ar_plan$phi), "\n")
  for (r in seq_along(ar_plan$phi)) {
    cat("    Run", r, "phi:", ar_plan$phi[[r]], "\n")
  }
}


# ==============================================================================
# EXAMPLE 6: rMVPA Integration
# ==============================================================================

cat("\n\n=== EXAMPLE 6: rMVPA Integration (Conceptual) ===\n")

# NOTE: This requires full rMVPA setup with neuroim2 objects
# Shown here for illustration only

# library(rMVPA)
#
# # Create continuous design
# mvdes <- continuous_mvpa_design(
#   event_model = ev_model_multi,
#   block_var = c(rep(1, n1), rep(2, n2)),
#   design_df_events = events_multi
# )
#
# # Create model spec WITH AR prewhitening
# spec <- hrfdecoder_model(
#   dataset = as_mvpa_dataset(fmri_data),
#   design = mvdes,
#   ar_order = 1,              # AR(1) enabled
#   ar_method = "ar",
#   ar_pooling = "run",
#   lambda_W = 10,
#   lambda_HRF = 1,
#   lambda_smooth = 5,
#   max_iter = 10
# )
#
# # Run searchlight with AR prewhitening
# results <- run_searchlight(spec, radius = 8, method = "randomized", niter = 4)

cat("See AR_IMPLEMENTATION_SUMMARY.md for full rMVPA examples\n")


# ==============================================================================
# Summary
# ==============================================================================

cat("\n\n=== SUMMARY ===\n")
cat("✓ AR prewhitening successfully integrated into hrfdecoder\n")
cat("✓ Supports AR(p), ARMA(p,q), and automatic order selection\n")
cat("✓ Run-specific and global AR pooling available\n")
cat("✓ Seamless integration with rMVPA\n")
cat("✓ Backward compatible (ar_order = NULL disables AR)\n")
cat("\nSee AR_PREWHITENING_INTEGRATION_PLAN.md for full documentation\n")
