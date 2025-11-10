#' Internal metrics utilities for classification
#' @keywords internal
one_hot <- function(y, levels) {
  y <- factor(y, levels = levels)
  Y <- matrix(0, nrow = length(y), ncol = length(levels))
  colnames(Y) <- levels
  idx <- as.integer(y)
  ok <- !is.na(idx)
  Y[cbind(which(ok), idx[ok])] <- 1
  Y
}

#' @keywords internal
metric_accuracy <- function(probs, y_true) {
  pred <- colnames(probs)[max.col(probs)]
  mean(pred == as.character(y_true))
}

#' Fast binary AUC (Mann–Whitney)
#' @keywords internal
metric_auc_binary <- function(p_pos, y01) {
  # y01 must be 0/1 vector; drop NAs conservatively
  ok <- !is.na(y01) & !is.na(p_pos)
  p_pos <- p_pos[ok]
  y01 <- y01[ok]
  n1 <- sum(y01 == 1L)
  n0 <- sum(y01 == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(p_pos, ties.method = "average")
  sum_r1 <- sum(r[y01 == 1L])
  (sum_r1 - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' @keywords internal
metric_auc <- function(probs, y_true) {
  lev <- colnames(probs)
  K <- ncol(probs)
  if (K == 2L) {
    y <- factor(y_true, levels = lev)
    y01 <- as.integer(y) - 1L
    return(metric_auc_binary(probs[, lev[1L]], y01))
  }
  # Macro one-vs-rest AUC for multi-class
  vals <- vapply(seq_len(K), function(i) {
    y_bin <- as.integer(factor(y_true, levels = lev))
    y01 <- as.integer(y_bin == i)
    metric_auc_binary(probs[, i], y01)
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}

#' Macro one-vs-one AUC (Hand & Till, 2001)
#' For each pair (i, j), compute AUC for class i vs j using
#' P_i restricted to samples of classes i or j, then average 0.5*(A_ij + A_ji)
#' across all i < j.
#' @keywords internal
metric_auc_ovo <- function(probs, y_true) {
  lev <- colnames(probs)
  K <- ncol(probs)
  if (K <= 2L) {
    return(metric_auc(probs, y_true))
  }
  y <- factor(y_true, levels = lev)
  pairs <- combn(seq_len(K), 2)
  vals <- apply(pairs, 2, function(idx) {
    i <- idx[1]; j <- idx[2]
    sel <- (y %in% lev[c(i, j)])
    y_sel <- y[sel]
    p_i <- probs[sel, i]
    # A_ij: class i positive
    y01_ij <- as.integer(y_sel == lev[i])
    a_ij <- metric_auc_binary(p_i, y01_ij)
    # A_ji: class j positive (use P_j)
    p_j <- probs[sel, j]
    y01_ji <- as.integer(y_sel == lev[j])
    a_ji <- metric_auc_binary(p_j, y01_ji)
    0.5 * (a_ij + a_ji)
  })
  mean(vals, na.rm = TRUE)
}

#' @keywords internal
metric_logloss <- function(probs, y_true, eps = 1e-12) {
  lev <- colnames(probs)
  y <- factor(y_true, levels = lev)
  idx <- as.integer(y)
  p_true <- probs[cbind(seq_len(nrow(probs)), pmax(idx, 1L))]
  p_true[is.na(p_true)] <- eps
  -mean(log(pmin(pmax(p_true, eps), 1 - eps)))
}

#' @keywords internal
metric_brier <- function(probs, y_true) {
  lev <- colnames(probs)
  Y <- one_hot(y_true, lev)
  mean(rowSums((probs - Y)^2))
}

#' @keywords internal
metric_balanced_accuracy <- function(probs, y_true) {
  lev <- colnames(probs)
  pred <- colnames(probs)[max.col(probs)]
  y <- factor(y_true, levels = lev)
  p <- factor(pred, levels = lev)
  recalls <- vapply(lev, function(l) {
    idx <- (y == l)
    if (!any(idx)) return(NA_real_)
    mean(p[idx] == l)
  }, numeric(1))
  mean(recalls, na.rm = TRUE)
}

#' Compute a set of metrics from probs/y_true
#' @keywords internal
compute_metrics <- function(probs, y_true, metrics = c("accuracy", "auc")) {
  # Ensure column names and aligned factor levels
  if (is.null(colnames(probs))) {
    lev <- levels(as.factor(y_true))
    if (!is.null(lev) && length(lev) == ncol(probs)) {
      colnames(probs) <- lev
    } else {
      colnames(probs) <- paste0("C", seq_len(ncol(probs)))
    }
  }
  y_true <- factor(y_true, levels = colnames(probs))
  vals <- lapply(metrics, function(m) {
    switch(m,
      accuracy = metric_accuracy(probs, y_true),
      auc = metric_auc(probs, y_true),
      auc_ovr = metric_auc(probs, y_true),
      auc_ovo = metric_auc_ovo(probs, y_true),
      logloss = metric_logloss(probs, y_true),
      brier = metric_brier(probs, y_true),
      bal_acc = metric_balanced_accuracy(probs, y_true),
      NA_real_
    )
  })
  tibble::tibble(metric = metrics, value = unlist(vals))
}

#' Compute metrics from prediction probabilities and true labels
#' @param probs matrix of probabilities (N x K)
#' @param y_true factor of true labels with same levels as colnames(probs)
#' @param metrics character vector of metrics: "accuracy", "auc", "logloss", "brier", "bal_acc"
#' @export
hrf_metrics <- function(probs, y_true, metrics = c("accuracy", "auc")) {
  compute_metrics(probs, y_true, metrics)
}
