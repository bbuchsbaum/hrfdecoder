# Aggregate TR-level soft labels to events

Aggregate TR-level soft labels to events

## Usage

``` r
aggregate_events(
  P,
  events,
  TR,
  conditions,
  window = c(4, 8),
  hrf = NULL,
  normalize = FALSE
)
```

## Arguments

- P:

  matrix (T x K_event)

- events:

  event data.frame (needs columns onset, condition)

- TR:

  TR duration (seconds)

- conditions:

  ordered condition labels

- window:

  time window (s) after onset

- hrf:

  optional fmrihrf HRF object for weighting

## Value

list with probs matrix and y_true factor
