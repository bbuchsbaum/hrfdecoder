#' Create an hrfdecoder model spec for rMVPA
#' @param dataset optional mvpa_dataset object (can be NULL if using rMVPA's run_searchlight)
#' @param design mvpa_design object containing event_model and baseline_model
#' @param lambda_W ridge penalty on decoder weights (default: 10)
#' @param lambda_HRF prior strength pulling HRF toward canonical shape (default: 1)
#' @param lambda_smooth temporal smoothness penalty on soft labels (default: 5)
#' @param theta_penalty L2 penalty on HRF basis coefficients (default: 0.01)
#' @param basis optional fmrihrf HRF basis specification (default: NULL uses "spmg1")
#' @param window time window in seconds for trial aggregation, e.g., c(4, 8) (default: c(4, 8))
#' @param nonneg logical; if TRUE, enforce non-negative decoder weights (default: TRUE)
#' @param max_iter maximum number of alternating least squares iterations (default: 10)
#' @param tol convergence tolerance for ALS algorithm (default: 1e-4)
#' @param ar_order AR order for prewhitening (default: 1 for AR(1)). Set to NULL or 0 to disable.
#' @param ar_method AR estimation method: "ar" or "arma". Default: "ar".
#' @param ar_pooling Spatial pooling for AR: "global" or "run". Default: "run".
#' @param performance optional custom performance function (default: NULL uses built-in)
#' @param metrics character vector of performance metrics to compute, e.g., c("accuracy", "auc") (default: c("accuracy", "auc"))
#' @param primary_metric character; which metric to use as primary for reporting (default: "auc")
#' @param crossval optional cross-validation specification (default: NULL uses blocked_cross_validation if block_var exists)
#' @param return_predictions logical; if TRUE, return predictions for each fold (default: TRUE)
#' @param return_fits logical; if TRUE, return fitted model objects for each fold (default: FALSE)
#' @export
hrfdecoder_model <- function(dataset = NULL, design,
                             lambda_W = 10,
                             lambda_HRF = 1,
                             lambda_smooth = 5,
                             theta_penalty = 0.01,
                             basis = NULL,
                             window = c(4, 8),
                             nonneg = TRUE,
                             max_iter = 10,
                             tol = 1e-4,
                             ar_order = 1,
                             ar_method = c("ar", "arma"),
                             ar_pooling = c("run", "global"),
                             performance = NULL,
                             metrics = c("accuracy", "auc"),
                             primary_metric = "auc",
                             crossval = NULL,
                             return_predictions = TRUE,
                             return_fits = FALSE) {
  if (is.null(crossval) && !is.null(design$block_var)) {
    crossval <- rMVPA::blocked_cross_validation(design$block_var)
  }
  ar_method <- match.arg(ar_method)
  ar_pooling <- match.arg(ar_pooling)
  obj <- list(
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
    ar_order = ar_order,
    ar_method = ar_method,
    ar_pooling = ar_pooling,
    crossval = crossval,
    performance = performance,
    compute_performance = TRUE,
    return_predictions = return_predictions,
    metrics = metrics,
    primary_metric = primary_metric,
    return_fits = return_fits
  )
  class(obj) <- "hrfdecoder_model"
  obj
}

#' @export
y_train.hrfdecoder_model <- function(obj) {
  n <- nrow(obj$design$train_design)
  seq_len(n)
}

#' @export
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
    max_iter = obj$max_iter,
    tol = obj$tol,
    nonneg = obj$nonneg,
    ar_order = obj$ar_order,
    ar_method = obj$ar_method,
    ar_pooling = obj$ar_pooling,
    standardize = FALSE,
    verbose = 0
  )
  structure(
    list(fit = fit, sl_info = sl_info, indices = indices),
    class = "hrfdecoder_fit_wrap"
  )
}

#' @export
format_result.hrfdecoder_model <- function(obj, result, error_message = NULL, context, ...) {
  if (!is.null(error_message)) {
    return(tibble::tibble(
      class = list(NULL),
      probs = list(NULL),
      y_true = list(NULL),
      test_ind = list(context$test_ind %||% NA_integer_),
      fit = list(NULL),
      error = TRUE,
      error_message = error_message
    ))
  }
  Xtest <- as.matrix(context$test)
  preds <- predict(
    object = result$fit,
    newdata = Xtest,
    ev_model_test = obj$design$event_model,
    mode = "trial",
    window = obj$window
  )
  probs <- preds$probs
  classes <- factor(colnames(probs)[max.col(probs)], levels = colnames(probs))
  tibble::tibble(
    class = list(classes),
    probs = list(probs),
    y_true = list(preds$y_true),
    test_ind = list(context$test_ind %||% seq_len(nrow(Xtest))),
    fit = list(if (isTRUE(obj$return_fits)) result$fit else NULL),
    error = FALSE,
    error_message = NA_character_
  )
}

#' @export
merge_results.hrfdecoder_model <- function(obj, result_set, indices, id, ...) {
  if (any(result_set$error)) {
    emsg <- result_set$error_message[which(result_set$error)[1]]
    return(tibble::tibble(
      result = list(NULL),
      indices = list(indices),
      performance = list(NULL),
      id = id,
      error = TRUE,
      error_message = emsg
    ))
  }
  combined <- wrap_classification_result_from_folds(result_set)
  perf <- compute_acc_perf(combined, metrics = obj$metrics %||% c("accuracy"))
  # Flag primary metric for convenience
  if (!is.null(obj$primary_metric)) {
    perf$primary <- perf$metric == obj$primary_metric
  }
  # Also surface primary metric/value at the top-level for easy tabulation
  pm <- obj$primary_metric %||% (perf$metric[1] %||% NA_character_)
  pv <- tryCatch({ perf$value[match(pm, perf$metric)] }, error = function(e) NA_real_)
  tibble::tibble(
    result = list(combined),
    indices = list(indices),
    performance = list(perf),
    primary_metric = pm,
    primary_value = pv,
    id = id,
    error = FALSE,
    error_message = NA_character_
  )
}

#' @keywords internal
wrap_classification_result_from_folds <- function(result_set) {
  probs <- do.call(rbind, result_set$probs)
  y_true <- do.call(c, lapply(result_set$y_true, as.character))
  y_true <- factor(y_true, levels = colnames(probs))
  classes <- factor(colnames(probs)[max.col(probs)], levels = colnames(probs))
  list(class = classes, probs = probs, y_true = y_true)
}

#' @keywords internal
compute_acc_perf <- function(result, metrics = c("accuracy")) {
  probs <- result$probs
  y_true <- result$y_true
  compute_metrics(probs, y_true, metrics = metrics)
}
