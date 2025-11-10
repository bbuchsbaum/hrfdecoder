#' Estimate nuisance rank from whitened residuals
#'
#' Uses a Marčenko–Pastur bulk-edge threshold together with
#' a parallel-analysis null (voxel-wise circular shifts) to select the
#' low-rank nuisance dimension after whitening and optional baseline removal.
#'
#' @param Y Numeric matrix (T x V) of voxel time series.
#' @param X_list List of design blocks (each T x K). If no design exists,
#'   pass a single all-zero column to keep the code path identical to the fitter.
#' @param X_base Optional baseline design (T x B) subtracted after whitening.
#' @param runs Optional integer vector of length T with run identifiers.
#' @param control Optional list with entries `lambda_A_init_ridge` and
#'   `lambda_base` for the initial ridge fits.
#' @param alpha Family-wise level for the parallel-analysis quantile.
#' @param B Number of circular-shift null draws (set to 0 to disable).
#' @param r_max Maximum number of singular values to probe.
#' @param block_shifts Logical; if TRUE, circular shifts are done voxel-wise.
#'
#' @return List with fields `r`, `r_mp`, `r_pa`, `sigma2`, `lambda_plus`, `ev`.
#' @export
estimate_rank_auto <- function(Y,
                               X_list,
                               X_base = NULL,
                               runs = NULL,
                               control = list(),
                               alpha = 0.05,
                               B = 50L,
                               r_max = 50L,
                               block_shifts = TRUE) {
  stopifnot(is.matrix(Y))
  stopifnot(is.list(X_list), length(X_list) >= 1L)
  Tlen <- nrow(Y)
  V <- ncol(Y)
  if (is.null(runs)) runs <- rep(1L, Tlen)
  if (length(runs) != Tlen) stop("runs must have length equal to nrow(Y).")

  # Construct the joint design matrix (task blocks + optional baseline)
  block_lengths <- vapply(X_list, nrow, integer(1))
  if (any(block_lengths != Tlen)) stop("All X_list blocks must match Y rows.")

  design_blocks <- lapply(X_list, function(Xb) {
    if (!is.matrix(Xb)) stop("Each element of X_list must be a matrix.")
    Xb
  })
  X_full <- do.call(cbind, c(design_blocks, if (!is.null(X_base)) list(X_base) else list()))
  if (is.null(X_full) || ncol(X_full) == 0L) {
    X_full <- matrix(0, Tlen, 1L)
    design_blocks <- list(X_full)
  }

  lamA0 <- control$lambda_A_init_ridge
  if (is.null(lamA0)) lamA0 <- 1e-6
  XtX <- crossprod(X_full)
  beta0 <- solve(XtX + lamA0 * diag(ncol(X_full)), crossprod(X_full, Y))
  resid0 <- Y - X_full %*% beta0

  noise_plan <- fmriAR::fit_noise(resid = resid0, runs = runs, method = "ar", p = "auto")
  whitened <- fmriAR::whiten_apply(noise_plan, X_full, Y, runs = runs)
  Yw <- whitened$Y
  Xw_all <- whitened$X

  block_sizes <- vapply(design_blocks, ncol, integer(1))
  base_cols <- if (!is.null(X_base)) ncol(X_base) else 0L
  split_idx <- cumsum(c(0L, block_sizes, base_cols))

  Xw_base <- NULL
  if (base_cols > 0L) {
    rng <- (split_idx[length(block_sizes) + 1L] + 1L):split_idx[length(block_sizes) + 2L]
    Xw_base <- Xw_all[, rng, drop = FALSE]
  }

  if (!is.null(Xw_base)) {
    lamB <- control$lambda_base
    if (is.null(lamB)) lamB <- 1e-3
    Bdim <- ncol(Xw_base)
    Bcoef <- solve(crossprod(Xw_base) + lamB * diag(Bdim), crossprod(Xw_base, Yw))
    E <- Yw - Xw_base %*% Bcoef
  } else {
    E <- Yw
  }

  S <- tcrossprod(E) / V
  ev <- sort(Re(eigen(S, symmetric = TRUE, only.values = TRUE)$values), decreasing = TRUE)
  q <- min(Tlen, V) / max(Tlen, V)
  tail_n <- max(1L, floor(0.3 * length(ev)))
  sigma2_hat <- mean(tail(ev, tail_n))
  lambda_plus <- sigma2_hat * (1 + sqrt(q))^2
  r_mp <- sum(ev > lambda_plus)

  r_probe <- min(r_max, length(ev))
  sv_obs <- sqrt(ev[seq_len(r_probe)] * V)
  thr <- rep(Inf, r_probe)
  if (B > 0L) {
    sv_null <- matrix(NA_real_, r_probe, B)
    for (b in seq_len(B)) {
      Ep <- E
      if (block_shifts) {
        shifts <- sample.int(Tlen, V, replace = TRUE) - 1L
        for (i in seq_len(V)) {
          s <- shifts[i]
          if (s > 0L) {
            Ep[, i] <- c(E[(s + 1L):Tlen, i], E[seq_len(s), i])
          }
        }
      } else {
        perm <- sample.int(Tlen)
        Ep <- Ep[perm, , drop = FALSE]
      }
      Sp <- tcrossprod(Ep) / ncol(Ep)
      evp <- sort(Re(eigen(Sp, symmetric = TRUE, only.values = TRUE)$values), decreasing = TRUE)
      sv_null[, b] <- sqrt(evp[seq_len(r_probe)] * ncol(Ep))
    }
    thr <- apply(sv_null, 1L, stats::quantile, probs = 1 - alpha)
  }
  r_pa <- sum(sv_obs > thr)

  list(
    r = min(r_mp, r_pa),
    r_mp = r_mp,
    r_pa = r_pa,
    sigma2 = sigma2_hat,
    lambda_plus = lambda_plus,
    ev = ev
  )
}
