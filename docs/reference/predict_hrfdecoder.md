# Predict with an hrfdecoder fit

Predict with an hrfdecoder fit

## Usage

``` r
predict_hrfdecoder(
  object,
  Y_test,
  ev_model_test = NULL,
  mode = c("tr", "trial"),
  window = c(4, 8),
  weights = c("hrf", "flat")
)
```

## Arguments

- object:

  hrfdecoder_fit

- Y_test:

  numeric matrix (T x V)

- ev_model_test:

  optional fmridesign::event_model for trial-level outputs

- mode:

  "tr" or "trial"

- window:

  time window (seconds) relative to onset for aggregation

- weights:

  weighting scheme ("hrf" uses fitted HRF; "flat" uniform)
