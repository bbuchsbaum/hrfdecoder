#' Convert to rMVPA dataset
#'
#' Generic function to convert various fMRI data objects to rMVPA's mvpa_dataset format
#'
#' @param x object to convert (e.g., fmri_dataset)
#' @param ... additional arguments passed to methods
#' @return an mvpa_dataset object
#' @export
as_mvpa_dataset <- function(x, ...) {
  UseMethod("as_mvpa_dataset")
}

#' @describeIn as_mvpa_dataset Convert fmri_dataset to mvpa_dataset
#' @param mask optional mask to use (defaults to x$mask)
#' @export
as_mvpa_dataset.fmri_dataset <- function(x, mask = NULL, ...) {
  mask <- mask %||% x$mask
  rMVPA::mvpa_dataset(train_data = x$neurovec, mask = mask, ...)
}

#' Continuous-time MVPA design wrapper
#' @param event_model fmridesign::event_model
#' @param block_var run/block ids per TR
#' @param design_df_events optional trial table (defaults to events(event_model))
#' @param split_by optional formula passed to mvpa_design
#' @export
continuous_mvpa_design <- function(event_model,
                                   block_var,
                                   design_df_events = NULL,
                                   split_by = NULL) {
  stopifnot(inherits(event_model, "event_model"))
  block_var <- block_var %||% rep(1, nrow(fmridesign::design_matrix(event_model)))
  y_dummy <- seq_along(block_var)
  design_df <- data.frame(y = y_dummy, block = block_var)
  mvdes <- rMVPA::mvpa_design(
    train_design = design_df,
    y_train = ~ y,
    block_var = ~ block,
    split_by = split_by
  )
  mvdes$event_model <- event_model
  mvdes$events <- design_df_events %||% fmridesign::events(event_model)
  mvdes$block_var <- block_var
  mvdes
}
