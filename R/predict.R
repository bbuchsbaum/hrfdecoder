#' Predict with an hrfdecoder fit
#'
#' @param object hrfdecoder_fit
#' @param Y_test numeric matrix (T x V)
#' @param ev_model_test optional fmridesign::event_model for trial-level outputs
#' @param mode "tr" or "trial"
#' @param window time window (seconds) relative to onset for aggregation
#' @param weights weighting scheme ("hrf" uses fitted HRF; "flat" uniform)
#'
#' @export
predict_hrfdecoder <- function(object, Y_test,
                               ev_model_test = NULL,
                               mode = c("tr", "trial"),
                               window = c(4, 8),
                               weights = c("hrf", "flat")) {
  # Gentle deprecation notice; can be silenced via option
  if (isTRUE(getOption("hrfdecode.deprecate.predict_hrfdecoder", TRUE))) {
    .Deprecated("predict", package = "hrfdecode",
                msg = "predict_hrfdecoder() is soft-deprecated; use predict(object, newdata = ..., ...) instead. Set options(hrfdecode.deprecate.predict_hrfdecoder = FALSE) to silence this message.")
  }
  stopifnot(inherits(object, "hrfdecoder_fit"))
  stopifnot(is.matrix(Y_test))
  mode <- match.arg(mode)
  weights <- match.arg(weights)

  # Apply preprocessing in same order as training

  # STEP 1: AR prewhitening (if was applied during training)
  if (!is.null(object$preproc$ar_plan)) {
    run_ids_test <- .get_run_ids_from_test_data(object, Y_test, ev_model_test)
    Y_white <- fmriAR::whiten_apply(
      object$preproc$ar_plan,
      X = matrix(0, nrow(Y_test), 1L),
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

  scores <- predict_softlabels(Y_test, object$W, object$b)
  probs <- row_softmax(scores)
  if (mode == "tr") {
    colnames(probs) <- c(object$conditions, if (isTRUE(object$background)) "background")
    return(probs)
  }
  stopifnot(inherits(ev_model_test, "event_model"))
  events_tbl <- get_event_table(ev_model_test)
  P_event <- probs[, seq_along(object$conditions), drop = FALSE]
  hrf_weights <- if (weights == "hrf") object$hrf else NULL
  agg <- aggregate_events(
    P = P_event,
    events = events_tbl,
    TR = object$settings$TR,
    conditions = object$conditions,
    window = window,
    hrf = hrf_weights,
    normalize = TRUE
  )
  agg
}

#' Predict for hrfdecoder_fit objects
#'
#' S3 wrapper around `predict_hrfdecoder()` to integrate with base
#' `predict()` workflows. This preserves return types and arguments used
#' in `predict_hrfdecoder()` while exposing the conventional
#' `predict(object, newdata, ...)` interface.
#'
#' @param object hrfdecoder_fit
#' @param newdata numeric matrix (T x V)
#' @param ev_model_test optional fmridesign::event_model for trial-level outputs
#' @param mode "tr" or "trial"
#' @param window time window (seconds) relative to onset for aggregation
#' @param weights weighting scheme ("hrf" uses fitted HRF; "flat" uniform)
#' @param ... unused
#' @return Same as `predict_hrfdecoder()`
#' @export
predict.hrfdecoder_fit <- function(object, newdata,
                                   ev_model_test = NULL,
                                   mode = c("tr", "trial"),
                                   window = c(4, 8),
                                   weights = c("hrf", "flat"),
                                   ...) {
  if (missing(newdata)) stop("newdata must be supplied for predict().")
  Y_test <- as.matrix(newdata)
  predict_hrfdecoder(
    object = object,
    Y_test = Y_test,
    ev_model_test = ev_model_test,
    mode = mode,
    window = window,
    weights = weights
  )
}

#' Aggregate TR-level soft labels to events
#' @param P matrix (T x K_event)
#' @param events event data.frame (needs columns onset, condition)
#' @param TR TR duration (seconds)
#' @param conditions ordered condition labels
#' @param window time window (s) after onset
#' @param hrf optional fmrihrf HRF object for weighting
#' @param normalize logical; if TRUE, normalize per-event probabilities to sum to 1 across conditions (default: FALSE)
#' @return list with probs matrix and y_true factor
#' @export
aggregate_events <- function(P, events, TR, conditions,
                             window = c(4, 8), hrf = NULL, normalize = FALSE) {
  stopifnot(is.matrix(P))
  stopifnot(all(c("onset", "condition") %in% names(events)))
  K <- length(conditions)
  if (ncol(P) != K) {
    stop("Number of columns in P must match conditions.")
  }
  window <- sort(window)
  tgrid <- seq(window[1], window[2], by = TR)
  if (length(tgrid) == 0L) tgrid <- window[2]
  weights <- if (is.null(hrf)) rep(1, length(tgrid)) else as.numeric(fmrihrf::evaluate(hrf, tgrid))
  weights <- weights / sum(weights)
  NE <- nrow(events)
  probs <- matrix(0, nrow = NE, ncol = K)
  colnames(probs) <- conditions

  for (i in seq_len(NE)) {
    onset <- events$onset[i]
    start_idx <- floor((onset + window[1]) / TR) + 1L
    idx <- start_idx + seq_along(weights) - 1L
    idx <- idx[idx >= 1L & idx <= nrow(P)]
    if (length(idx) == 0L) next
    w <- weights[seq_along(idx)]
    slice <- P[idx, , drop = FALSE]
    probs[i, ] <- colSums(slice * w)
  }
  if (isTRUE(normalize)) {
    # Normalize per-event to sum to 1 across conditions
    rs <- rowSums(probs)
    nz <- rs > 0
    if (any(nz)) probs[nz, ] <- probs[nz, , drop = FALSE] / rs[nz]
  }
  y_true <- factor(events$condition, levels = conditions)
  list(probs = probs, y_true = y_true)
}

#' Get run IDs for test data (internal helper)
#' @keywords internal
.get_run_ids_from_test_data <- function(fit_obj, Y_test, ev_model_test) {
  # Option 1: Extract from ev_model_test if provided
  if (!is.null(ev_model_test) && inherits(ev_model_test, "event_model")) {
    prep <- prepare_decoder_inputs(ev_model_test, hrf = fit_obj$hrf,
                                   background = fit_obj$background)
    return(prep$run_ids)
  }

  # Option 2: If test data has same length as training, use training run_ids
  if (nrow(Y_test) == length(fit_obj$settings$run_ids)) {
    return(fit_obj$settings$run_ids)
  }

  # Option 3: Default to single run
  return(rep(1L, nrow(Y_test)))
}
