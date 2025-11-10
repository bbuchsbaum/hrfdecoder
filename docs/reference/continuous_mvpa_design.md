# Continuous-time MVPA design wrapper

Continuous-time MVPA design wrapper

## Usage

``` r
continuous_mvpa_design(
  event_model,
  block_var,
  design_df_events = NULL,
  split_by = NULL
)
```

## Arguments

- event_model:

  fmridesign::event_model

- block_var:

  run/block ids per TR

- design_df_events:

  optional trial table (defaults to events(event_model))

- split_by:

  optional formula passed to mvpa_design
