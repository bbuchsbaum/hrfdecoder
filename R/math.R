#' @keywords internal
row_softmax <- function(X) {
  X <- as.matrix(X)
  if (nrow(X) == 0L) return(X)
  max_per_row <- apply(X, 1, max)
  exps <- exp(X - max_per_row)
  denom <- rowSums(exps)
  exps / denom
}
