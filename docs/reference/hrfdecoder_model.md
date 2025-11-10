# Create an hrfdecoder model spec for rMVPA

Create an hrfdecoder model spec for rMVPA

## Usage

``` r
hrfdecoder_model(
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
  tol = 1e-04,
  ar_order = 1,
  ar_method = c("ar", "arma"),
  ar_pooling = c("run", "global"),
  performance = NULL,
  crossval = NULL,
  return_predictions = TRUE,
  return_fits = FALSE
)
```

## Arguments

- ar_order:

  AR order for prewhitening (default: 1 for AR(1)). Set to NULL or 0 to
  disable.

- ar_method:

  AR estimation method: "ar" or "arma". Default: "ar".

- ar_pooling:

  Spatial pooling for AR: "global" or "run". Default: "run".
