# Fit HRF-aware weakly supervised decoder

Fit HRF-aware weakly supervised decoder

## Usage

``` r
fit_hrfdecoder(
  Y,
  ev_model,
  base_model = NULL,
  hrf = NULL,
  lambda_W = 10,
  lambda_HRF = 1,
  lambda_smooth = 5,
  theta_penalty = 0.01,
  max_iter = 20,
  tol = 1e-04,
  nonneg = TRUE,
  background = TRUE,
  standardize = TRUE,
  ar_order = NULL,
  ar_method = c("ar", "arma"),
  ar_pooling = c("run", "global"),
  verbose = 1
)
```

## Arguments

- Y:

  numeric matrix (T x V) of fMRI data (time by voxel)

- ev_model:

  fmridesign::event_model describing events

- base_model:

  optional fmridesign::baseline_model for nuisance removal

- hrf:

  optional fmrihrf basis object

- lambda_W:

  ridge penalty on decoder weights

- lambda_HRF:

  adherence weight to HRF prior

- lambda_smooth:

  temporal smoothness weight

- theta_penalty:

  ridge penalty on HRF basis coefficients

- max_iter:

  maximum ALS iterations

- tol:

  convergence tolerance on P updates

- nonneg:

  enforce non-negative soft labels

- background:

  include a background column in the soft labels

- standardize:

  z-score Y columns before fitting

- ar_order:

  AR order for prewhitening (default: NULL for no AR prewhitening). Set
  to 1 for AR(1), 2 for AR(2), etc. Use "auto" for automatic BIC-based
  selection.

- ar_method:

  AR estimation method: "ar" (Yule-Walker) or "arma" (Hannan-Rissanen).
  Default: "ar".

- ar_pooling:

  Spatial pooling for AR parameters: "global" (one AR model for all
  voxels) or "run" (separate AR model per run). Default: "run".

- verbose:

  integer verbosity level

## Value

object of class \`hrfdecoder_fit\`
