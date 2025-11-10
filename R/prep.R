#' @keywords internal
build_laplacian_from_runs <- function(run_ids) {
  stopifnot(length(run_ids) >= 1L)
  blocks <- split(seq_along(run_ids), run_ids)
  mats <- lapply(blocks, function(idx) {
    n <- length(idx)
    if (n <= 2L) {
      return(Matrix::Diagonal(n))
    }
    rseq <- seq_len(n - 2L)
    rows <- rep(rseq, each = 3L)
    cols <- as.vector(rbind(rseq, rseq + 1L, rseq + 2L))
    vals <- rep(c(1, -2, 1), times = length(rseq))
    D <- Matrix::sparseMatrix(i = rows, j = cols, x = vals, dims = c(n - 2L, n))
    Matrix::crossprod(D)
  })
  Matrix::bdiag(mats)
}

#' @keywords internal
build_softlabel_prior <- function(X_list, theta, background = TRUE) {
  stopifnot(length(X_list) >= 1L)
  Tn <- nrow(X_list[[1L]])
  G <- vapply(
    X_list,
    function(Xc) {
      Xmat <- as.matrix(Xc)
      as.numeric(Xmat %*% theta)
    },
    numeric(Tn)
  )
  if (background) {
    G <- cbind(G, rep(0, Tn))
  }
  list(DBbeta = G, P0 = row_softmax(G))
}

#' Prepare decoder design inputs
#' @keywords internal
prepare_decoder_inputs <- function(ev_model, hrf = NULL, background = TRUE) {
  interop <- build_condition_basis(ev_model, hrf)
  Pdim <- ncol(interop$X_list[[1L]])
  theta0 <- numeric(Pdim)
  theta0[1L] <- 1
  priors <- build_softlabel_prior(interop$X_list, theta0, background = background)
  blocklens <- interop$sframe$blocklens %||% nrow(interop$X_list[[1L]])
  run_ids <- rep(seq_along(blocklens), blocklens)
  list(
    X_list = interop$X_list,
    hrf = interop$hrf,
    sframe = interop$sframe,
    conditions = interop$conditions,
    DBbeta = priors$DBbeta,
    P0 = priors$P0,
    theta = theta0,
    run_ids = run_ids
  )
}

#' Estimate HRF basis coefficients from soft labels
#' @keywords internal
estimate_hrf_theta <- function(X_list, Z_event, hrf_obj, penalty = 0.01) {
  P <- ncol(X_list[[1L]])
  XtX <- matrix(0, P, P)
  Xtz <- numeric(P)
  for (j in seq_along(X_list)) {
    Xc <- as.matrix(X_list[[j]])
    XtX <- XtX + crossprod(Xc)
    Xtz <- Xtz + as.numeric(crossprod(Xc, Z_event[, j]))
  }
  Rmat <- fmrihrf::penalty_matrix(hrf_obj)
  theta <- solve(XtX + penalty * Rmat, Xtz)
  theta
}

#' @keywords internal
sampling_frame_tr <- function(sframe) {
  tr <- sframe$TR %||% attr(sframe, "TR", exact = TRUE)
  if (!is.null(tr)) return(tr)
  sam <- fmrihrf::samples(sframe)
  if (length(sam) > 1L) stats::median(diff(sam)) else 2
}
