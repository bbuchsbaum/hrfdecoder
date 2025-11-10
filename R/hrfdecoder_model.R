#' Create an hrfdecoder model spec for rMVPA
#' @param ar_order AR order for prewhitening (default: 1 for AR(1)). Set to NULL or 0 to disable.
#' @param ar_method AR estimation method: "ar" or "arma". Default: "ar".
#' @param ar_pooling Spatial pooling for AR: "global" or "run". Default: "run".
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
  preds <- predict_hrfdecoder(
    object = result$fit,
    Y_test = Xtest,
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
  perf <- compute_acc_perf(combined)
  tibble::tibble(
    result = list(combined),
    indices = list(indices),
    performance = list(perf),
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
compute_acc_perf <- function(result) {
  acc <- mean(result$class == result$y_true)
  tibble::tibble(metric = "accuracy", value = acc)
}
