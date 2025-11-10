# Changelog

All notable changes to hrfdecoder will be documented in this file.

## \[0.2.1\] - 2025-11-09

### Added

- **AR Prewhitening Support**: Full integration with fmriAR for temporal
  autocorrelation correction
  - New parameters in [`fit_hrfdecoder()`](reference/fit_hrfdecoder.md):
    `ar_order`, `ar_method`, `ar_pooling`
  - Supports AR(p), ARMA(p,q), and automatic order selection via BIC
  - Run-specific and global AR parameter pooling
  - Seamless test data application (AR plan stored and reused)
  - rMVPA integration with default `ar_order = 1` in
    [`hrfdecoder_model()`](reference/hrfdecoder_model.md)

### Changed

- Preprocessing pipeline order now: Baseline → AR Whitening →
  Standardization → ALS
- `fit$preproc` now includes `ar_plan` field when AR is enabled
- `fit$settings` expanded to store `ar_order`, `ar_method`, `ar_pooling`

### Technical Details

- AR prewhitening occurs after baseline residualization but before
  standardization
- AR parameters estimated from residuals using
  [`fmriAR::fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.html)
- Whitening applied via
  [`fmriAR::whiten_apply()`](https://bbuchsbaum.github.io/fmriAR/reference/whiten_apply.html)
  with run-aware boundary resets
- Test data automatically whitened using stored AR plan from training
- Helper function
  [`.get_run_ids_from_test_data()`](reference/dot-get_run_ids_from_test_data.md)
  extracts run structure for test application

### Documentation

- Comprehensive integration plan: `AR_PREWHITENING_INTEGRATION_PLAN.md`
- Implementation summary: `AR_IMPLEMENTATION_SUMMARY.md`
- Working examples: `examples/ar_prewhitening_example.R`
- Full test suite: `tests/testthat/test-ar-prewhitening.R`

### Backward Compatibility

- Fully backward compatible: `ar_order = NULL` (default) disables AR
  prewhitening
- Existing code runs unchanged and produces identical results
- No breaking changes to API

### Performance

- AR estimation overhead: ~5-10% of total runtime
- Benefits: Correct GLS weights, valid statistical inference, improved
  generalization
- Parallelizable via fmriAR’s OpenMP support

------------------------------------------------------------------------

## \[0.2.0\] - Previous Release

### Added

- Integration with fmridesign and fmrihrf for event modeling and HRF
  basis
- HRF estimation in basis space with penalty matrices
- Baseline residualization via
  [`fmridesign::residualize()`](https://bbuchsbaum.github.io/fmridesign/reference/residualize.html)
- rMVPA plugin architecture with
  [`hrfdecoder_model()`](reference/hrfdecoder_model.md)
- Multi-run support with block-diagonal Laplacian
- ALS solver with convergence diagnostics

### Features

- Joint estimation of soft labels (P), decoder weights (W), and HRF
- Temporal smoothness via second-difference Laplacian
- HRF-convolved design priors
- Trial-level and TR-level prediction modes
- Cross-validation compatible
