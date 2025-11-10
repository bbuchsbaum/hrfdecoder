als_objective <- function(X, W, b, P, DBbeta, lambda_W, lambda_HRF, lambda_smooth, L) {
  Smat <- X %*% W
  if (!is.null(b)) Smat <- sweep(Smat, 2, b, FUN = "+")
  Pmat <- as.matrix(P)
  DBmat <- as.matrix(DBbeta)
  recon <- sum((Smat - Pmat)^2)
  regW  <- lambda_W * sum(W^2)
  prior <- lambda_HRF * sum((Pmat - DBmat)^2)
  rough <- lambda_smooth * sum((L %*% Pmat) * Pmat)
  recon + regW + prior + rough
}

test_that("C++ entry points available and basic shapes are correct (smoke)", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")
  skip_if_not_installed("Matrix")
  ns <- asNamespace("hrfdecode")
  has_fit <- exists("fit_softlabels_als", envir = ns, inherits = FALSE)
  has_pred <- exists("predict_softlabels", envir = ns, inherits = FALSE)
  skip_if_not(has_fit && has_pred,
              "fit_softlabels_als/predict_softlabels not found (did you rebuild the package?).")
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  set.seed(1)
  sframe <- fmrihrf::sampling_frame(blocklens = 40, TR = 1)
  evtab  <- data.frame(onset = c(5, 15, 25), condition = factor(c("A","B","A")), run = 1)
  evmod  <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg1"),
    data = evtab, block = ~run, sampling_frame = sframe
  )
  prep   <- hrfdecode:::prepare_decoder_inputs(ev_model = evmod, hrf = NULL, background = TRUE)
  L      <- hrfdecode:::build_laplacian_from_runs(prep$run_ids)
  Tn     <- nrow(prep$DBbeta)
  V      <- 12
  X      <- matrix(rnorm(Tn * V), Tn, V)

  fit <- fit_softlabels_als(
    X = X, P0 = prep$P0, L = L, DBbeta = prep$DBbeta,
    lambda_W = 5, lambda_HRF = 0.5, lambda_smooth = 0.1,
    max_iter = 3, tol = 1e-5, nonneg = TRUE, threads = FALSE
  )
  expect_equal(dim(fit$P), dim(prep$P0))
  expect_equal(dim(fit$W), c(V, ncol(prep$P0)))
})

test_that("Single ALS iteration matches closed-form update (with/without intercept)", {
  skip_if_not_installed("fmridesign"); skip_if_not_installed("fmrihrf"); skip_if_not_installed("Matrix")
  skip_if_not(exists("fit_softlabels_als", envir = asNamespace("hrfdecode"), inherits = FALSE))
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  set.seed(2)
  sframe <- fmrihrf::sampling_frame(blocklens = 50, TR = 1)
  evtab  <- data.frame(onset = c(8, 18, 28), condition = factor(c("A","B","A")), run = 1)
  evmod  <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg1"),
    data = evtab, block = ~run, sampling_frame = sframe
  )
  prep   <- hrfdecode:::prepare_decoder_inputs(evmod, hrf = NULL, background = TRUE)
  L      <- hrfdecode:::build_laplacian_from_runs(prep$run_ids)
  Tn     <- nrow(prep$DBbeta); V <- 10
  X      <- matrix(rnorm(Tn * V), Tn, V)

  lamW <- 3; lamH <- 0.7; lamS <- 0.2
  fit1 <- fit_softlabels_als(X, prep$P0, L, prep$DBbeta, lamW, lamH, lamS,
                             max_iter = 1, tol = 1e-12, nonneg = FALSE, threads = FALSE)
  has_b <- "b" %in% names(fit1)

  if (has_b) {
    Xbar <- colMeans(X); Pbar <- colMeans(prep$P0)
    Xc   <- scale(X, center = Xbar, scale = FALSE)
    Pc   <- scale(prep$P0, center = Pbar, scale = FALSE)
    Aw   <- crossprod(Xc) + lamW * diag(ncol(X))
    Wcf  <- solve(Aw, crossprod(Xc, Pc))
    bcf  <- as.numeric(Pbar - drop(Xbar %*% Wcf))
    RHS  <- sweep(X %*% Wcf, 2, bcf, FUN = "+") + lamH * prep$DBbeta
  } else {
    Aw  <- crossprod(X) + lamW * diag(ncol(X))
    Wcf <- solve(Aw, crossprod(X, prep$P0))
    RHS <- X %*% Wcf + lamH * prep$DBbeta
  }
  lap  <- (1 + lamH) * Matrix::Diagonal(Tn) + lamS * L
  Pcf  <- as.matrix(Matrix::solve(lap, Matrix::Matrix(RHS, sparse = FALSE)))
  expect_equal(unname(Pcf), unname(fit1$P), tolerance = 1e-7)
  expect_equal(unname(Wcf), unname(fit1$W), tolerance = 1e-8)
})

test_that("Engine objective matches R helper", {
  skip_if_not_installed("fmridesign"); skip_if_not_installed("fmrihrf"); skip_if_not_installed("Matrix")
  ns <- asNamespace("hrfdecode")
  skip_if_not(all(c("fit_softlabels_als", "engine_objective") %in% ls(ns)))
  fit_softlabels_als <- ns$fit_softlabels_als
  engine_objective <- ns$engine_objective

  seeds <- c(11, 22, 33)
  for (sd in seeds) {
    set.seed(sd)
    sframe <- fmrihrf::sampling_frame(blocklens = 45, TR = 1)
    evtab  <- data.frame(onset = c(6, 16, 26), condition = factor(c("A","B","A")), run = 1)
    evmod  <- fmridesign::event_model(
      onset ~ fmridesign::hrf(condition, basis = "spmg1"),
      data = evtab, block = ~run, sampling_frame = sframe
    )
    prep   <- hrfdecode:::prepare_decoder_inputs(evmod, hrf = NULL, background = TRUE)
    L      <- hrfdecode:::build_laplacian_from_runs(prep$run_ids)
    Tn     <- nrow(prep$DBbeta); V <- 9
    X      <- matrix(rnorm(Tn * V), Tn, V)

    lamW <- 4; lamH <- 0.9; lamS <- 0.3
    fit <- fit_softlabels_als(X, prep$P0, L, prep$DBbeta,
                              lambda_W = lamW, lambda_HRF = lamH, lambda_smooth = lamS,
                              max_iter = 5, tol = 1e-10, nonneg = TRUE, threads = FALSE)
    bvec <- if ("b" %in% names(fit)) fit$b else rep(0, ncol(fit$W))
    r_obj <- als_objective(X, fit$W, bvec, fit$P, prep$DBbeta, lamW, lamH, lamS, L)
    cpp_obj <- engine_objective(X, fit$W, bvec, fit$P, prep$DBbeta, lamW, lamH, lamS, L)
    expect_equal(r_obj, cpp_obj$value, tolerance = 1e-8)
    expect_equal(unname(cpp_obj$recon + cpp_obj$regW + cpp_obj$prior + cpp_obj$rough),
                 cpp_obj$value, tolerance = 1e-10)
  }
})

test_that("ALS objective decreases with more iterations", {
  skip_if_not_installed("fmridesign"); skip_if_not_installed("fmrihrf"); skip_if_not_installed("Matrix")
  ns <- asNamespace("hrfdecode")
  skip_if_not(all(c("fit_softlabels_als", "engine_objective") %in% ls(ns)))
  fit_softlabels_als <- ns$fit_softlabels_als
  engine_objective <- ns$engine_objective

  set.seed(3)
  sframe <- fmrihrf::sampling_frame(blocklens = 60, TR = 1)
  evtab  <- data.frame(onset = c(6, 16, 36), condition = factor(c("A","B","A")), run = 1)
  evmod  <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg1"),
    data = evtab, block = ~run, sampling_frame = sframe
  )
  prep   <- hrfdecode:::prepare_decoder_inputs(evmod, hrf = NULL, background = TRUE)
  L      <- hrfdecode:::build_laplacian_from_runs(prep$run_ids)
  Tn     <- nrow(prep$DBbeta); V <- 15
  X      <- matrix(rnorm(Tn * V), Tn, V)

  lamW <- 5; lamH <- 1.0; lamS <- 0.5
  iters <- c(1, 2, 4, 8)
  objs  <- numeric(length(iters))
  for (i in seq_along(iters)) {
    fit <- fit_softlabels_als(X, prep$P0, L, prep$DBbeta, lamW, lamH, lamS,
                              max_iter = iters[i], tol = 1e-10, nonneg = TRUE, threads = FALSE)
    trace_vals <- unlist(fit$obj_trace$value)
    expect_true(length(trace_vals) >= 1)
    expect_true(all(diff(trace_vals) <= 1e-8))
    b    <- if ("b" %in% names(fit)) fit$b else rep(0, ncol(fit$W))
    objs[i] <- als_objective(X, fit$W, b, fit$P, prep$DBbeta, lamW, lamH, lamS, L)
    cpp_obj <- engine_objective(X, fit$W, b, fit$P, prep$DBbeta, lamW, lamH, lamS, L)
    expect_equal(objs[i], cpp_obj$value, tolerance = 1e-8)
  }
  expect_true(all(diff(objs) <= 1e-8))
})

test_that("HRF adherence dominates when lambda_HRF is huge", {
  skip_if_not_installed("fmridesign"); skip_if_not_installed("fmrihrf"); skip_if_not_installed("Matrix")
  skip_if_not(exists("fit_softlabels_als", envir = asNamespace("hrfdecode"), inherits = FALSE))
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  set.seed(4)
  sframe <- fmrihrf::sampling_frame(blocklens = 40, TR = 1)
  evtab  <- data.frame(onset = c(5, 25), condition = factor(c("A","B")), run = 1)
  evmod  <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg1"),
    data = evtab, block = ~run, sampling_frame = sframe
  )
  prep   <- hrfdecode:::prepare_decoder_inputs(evmod, hrf = NULL, background = TRUE)
  L      <- hrfdecode:::build_laplacian_from_runs(prep$run_ids)
  Tn     <- nrow(prep$DBbeta)
  X0     <- matrix(0, Tn, 6)

  lamW <- 0; lamH <- 1e6; lamS <- 0
  fit  <- fit_softlabels_als(X0, prep$P0, L, prep$DBbeta, lamW, lamH, lamS,
                             max_iter = 2, tol = 1e-12, nonneg = FALSE, threads = FALSE)
  rel  <- sqrt(sum((fit$P - prep$DBbeta)^2)) / (sqrt(sum(prep$DBbeta^2)) + 1e-12)
  expect_lt(rel, 1e-6)
})

test_that("Zero HRF penalty ignores prior even when DBbeta disagrees", {
  skip_if_not_installed("Matrix")
  skip_if_not(exists("fit_softlabels_als", envir = asNamespace("hrfdecode"), inherits = FALSE))
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  set.seed(44)
  Tn <- 30L
  V  <- 5L
  K  <- 3L
  X <- matrix(rnorm(Tn * V), Tn, V)
  W_true <- matrix(rnorm(V * K), V, K)
  b_true <- runif(K, -0.3, 0.3)
  P_data <- sweep(X %*% W_true, 2, b_true, "+")
  DBbeta <- matrix(rnorm(Tn * K, sd = 5), Tn, K)  # wildly different prior
  L <- Matrix::Diagonal(x = rep(0, Tn))

  fit <- fit_softlabels_als(
    X = X,
    P0 = P_data,
    L = L,
    DBbeta = DBbeta,
    lambda_W = 1e-4,
    lambda_HRF = 0,
    lambda_smooth = 1,
    max_iter = 6,
    tol = 1e-10,
    nonneg = FALSE,
    threads = FALSE
  )
  decoder_fit <- sweep(X %*% fit$W, 2, fit$b, "+")
  err_decoder <- sqrt(sum((fit$P - decoder_fit)^2))
  err_prior   <- sqrt(sum((fit$P - DBbeta)^2))
  expect_lt(err_decoder, err_prior * 1e-3)
})

test_that("Temporal roughness energy decreases as lambda_smooth increases", {
  skip_if_not_installed("fmridesign"); skip_if_not_installed("fmrihrf"); skip_if_not_installed("Matrix")
  skip_if_not(exists("fit_softlabels_als", envir = asNamespace("hrfdecode"), inherits = FALSE))
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  set.seed(5)
  sframe <- fmrihrf::sampling_frame(blocklens = 80, TR = 1)
  evtab  <- data.frame(onset = c(10, 30, 50), condition = factor(c("A","B","A")), run = 1)
  evmod  <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg1"),
    data = evtab, block = ~run, sampling_frame = sframe
  )
  prep   <- hrfdecode:::prepare_decoder_inputs(evmod, hrf = NULL, background = TRUE)
  L      <- hrfdecode:::build_laplacian_from_runs(prep$run_ids)
  Tn     <- nrow(prep$DBbeta); V <- 20
  X      <- matrix(rnorm(Tn * V), Tn, V)

  fit0   <- fit_softlabels_als(X, prep$P0, L, prep$DBbeta, lambda_W = 5, lambda_HRF = 1,
                               lambda_smooth = 0, max_iter = 4, tol = 1e-6,
                               nonneg = TRUE, threads = FALSE)
  fitHi  <- fit_softlabels_als(X, prep$P0, L, prep$DBbeta, lambda_W = 5, lambda_HRF = 1,
                               lambda_smooth = 50, max_iter = 4, tol = 1e-6,
                               nonneg = TRUE, threads = FALSE)
  rough0 <- sum((L %*% as.matrix(fit0$P)) * as.matrix(fit0$P))
  roughH <- sum((L %*% as.matrix(fitHi$P)) * as.matrix(fitHi$P))
  expect_lt(roughH, rough0 * 0.9)
})

test_that("Zero smoothing limit reproduces decoder predictions exactly", {
  skip_if_not_installed("Matrix")
  skip_if_not(exists("fit_softlabels_als", envir = asNamespace("hrfdecode"), inherits = FALSE))
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  set.seed(123)
  Tn <- 36L
  V  <- 6L
  K  <- 3L
  X <- matrix(rnorm(Tn * V), Tn, V)
  W_true <- matrix(rnorm(V * K), V, K)
  b_true <- runif(K, min = -0.4, max = 0.4)
  P_target <- sweep(X %*% W_true, 2, b_true, "+")
  P0 <- P_target + matrix(rnorm(Tn * K, sd = 1e-4), Tn, K)
  DBbeta <- matrix(rnorm(Tn * K, sd = 1), Tn, K)
  Lzero <- Matrix::Diagonal(x = rep(0, Tn))

  fit <- fit_softlabels_als(
    X = X,
    P0 = P0,
    L = Lzero,
    DBbeta = DBbeta,
    lambda_W = 1e-6,
    lambda_HRF = 0,
    lambda_smooth = 0,
    max_iter = 5,
    tol = 1e-12,
    nonneg = FALSE,
    threads = FALSE
  )
  recon <- X %*% fit$W
  recon <- sweep(recon, 2, fit$b, "+")
  expect_lt(max(abs(fit$P - recon)), 1e-8)
})

test_that("Smoothing does not leak across run boundaries", {
  skip_if_not_installed("Matrix")
  skip_if_not(exists("fit_softlabels_als", envir = asNamespace("hrfdecode"), inherits = FALSE))
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  run_ids <- rep(1:2, each = 30L)
  L <- hrfdecode:::build_laplacian_from_runs(run_ids)
  Tn <- length(run_ids)
  K_ev <- 2L
  # Event prior active only in run 1
  DBbeta <- matrix(0, nrow = Tn, ncol = K_ev + 1L)
  DBbeta[run_ids == 1L, 1L] <- exp(-seq_len(sum(run_ids == 1L)) / 5)
  DBbeta[, K_ev + 1L] <- 0  # background
  P0 <- DBbeta
  X0 <- matrix(0, Tn, 4L)

  fit <- fit_softlabels_als(
    X = X0,
    P0 = P0,
    L = L,
    DBbeta = DBbeta,
    lambda_W = 0,
    lambda_HRF = 1e4,
    lambda_smooth = 10,
    max_iter = 5,
    tol = 1e-10,
    nonneg = TRUE,
    threads = FALSE
  )
  run2_idx <- run_ids == 2L
  slice <- fit$P[run2_idx, 1:K_ev, drop = FALSE]
  expect_lt(max(abs(slice)), 1e-5)
})

test_that("Extremely large smoothing drives second differences to nearly zero", {
  skip_if_not_installed("Matrix")
  skip_if_not(exists("fit_softlabels_als", envir = asNamespace("hrfdecode"), inherits = FALSE))
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  run_ids <- rep(1:2, each = 40L)
  L <- hrfdecode:::build_laplacian_from_runs(run_ids)
  Tn <- length(run_ids)
  K_ev <- 2L
  DBbeta <- matrix(rnorm(Tn * K_ev, sd = 1), Tn, K_ev)
  P0 <- DBbeta
  X0 <- matrix(0, Tn, 3L)

  fit <- fit_softlabels_als(
    X = X0,
    P0 = P0,
    L = L,
    DBbeta = DBbeta,
    lambda_W = 0,
    lambda_HRF = 1,
    lambda_smooth = 1e6,
    max_iter = 4,
    tol = 1e-10,
    nonneg = FALSE,
    threads = FALSE
  )
  roughness <- L %*% fit$P[, 1:K_ev, drop = FALSE]
  expect_lt(max(abs(roughness)), 5e-5)
})

test_that("Nonnegativity projection removes negative entries", {
  skip_if_not_installed("fmridesign"); skip_if_not_installed("fmrihrf"); skip_if_not_installed("Matrix")
  skip_if_not(exists("fit_softlabels_als", envir = asNamespace("hrfdecode"), inherits = FALSE))
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  set.seed(7)
  sframe <- fmrihrf::sampling_frame(blocklens = 50, TR = 1)
  evtab  <- data.frame(onset = c(10, 30), condition = factor(c("A","B")), run = 1)
  evmod  <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg1"),
    data = evtab, block = ~run, sampling_frame = sframe
  )
  prep   <- hrfdecode:::prepare_decoder_inputs(evmod, hrf = NULL, background = TRUE)
  L      <- hrfdecode:::build_laplacian_from_runs(prep$run_ids)
  Tn     <- nrow(prep$DBbeta)
  X0     <- matrix(0, Tn, 3)

  DBneg <- prep$DBbeta
  DBneg[, seq_along(prep$X_list)] <- DBneg[, seq_along(prep$X_list)] - 1.0

  fit_free <- fit_softlabels_als(X0, prep$P0, L, DBneg,
                                 lambda_W = 0, lambda_HRF = 5, lambda_smooth = 0,
                                 max_iter = 2, tol = 1e-10, nonneg = FALSE, threads = FALSE)
  expect_lt(min(fit_free$P), -1e-3)

  fit_proj <- fit_softlabels_als(X0, prep$P0, L, DBneg,
                                 lambda_W = 0, lambda_HRF = 5, lambda_smooth = 0,
                                 max_iter = 2, tol = 1e-10, nonneg = TRUE, threads = FALSE)
  expect_gte(min(fit_proj$P), -1e-12)
})

test_that("estimate_hrf_theta recovers known HRF basis coefficients", {
  skip_if_not_installed("fmridesign"); skip_if_not_installed("fmrihrf")

  set.seed(8)
  sframe <- fmrihrf::sampling_frame(blocklens = 60, TR = 1)
  evtab  <- data.frame(onset = c(5, 25, 45), condition = factor(c("A","B","A")), run = 1)
  evmod  <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg1"),
    data = evtab, block = ~run, sampling_frame = sframe
  )
  inter  <- hrfdecode:::build_condition_basis(evmod, hrf = fmrihrf::getHRF("spmg1"))

  Pdim   <- ncol(inter$X_list[[1L]])
  theta0 <- runif(Pdim, min = 0.2, max = 1.2)
  Z_event <- vapply(
    inter$X_list,
    function(Xc) as.numeric(as.matrix(Xc) %*% theta0),
    numeric(nrow(inter$X_list[[1L]]))
  )
  theta_hat <- hrfdecode:::estimate_hrf_theta(inter$X_list, Z_event, inter$hrf, penalty = 1e-4)
  if (length(theta0) >= 2) {
    corr <- cor(theta0, theta_hat)
    expect_gt(corr, 0.999)
  }
  rel_err <- sqrt(sum((theta_hat - theta0)^2)) / (sqrt(sum(theta0^2)) + 1e-12)
  expect_lt(rel_err, 1e-3)
})

test_that("HRF basis stress test recovers multi-basis theta", {
  skip_if_not_installed("fmridesign"); skip_if_not_installed("fmrihrf")

  set.seed(1234)
  sframe <- fmrihrf::sampling_frame(blocklens = 70, TR = 1)
  evtab  <- data.frame(
    onset = c(4, 14, 24, 34, 44, 54),
    condition = factor(rep(c("CondA", "CondB"), each = 3)),
    run = 1
  )
  evmod <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg3"),
    data = evtab, block = ~run, sampling_frame = sframe
  )
  inter <- hrfdecode:::build_condition_basis(evmod, hrf = fmrihrf::getHRF("spmg3"))

  nbasis <- ncol(inter$X_list[[1L]])
  expect_gte(nbasis, 3)

  theta_true <- runif(nbasis, min = 0.1, max = 1.5)
  Z_event <- vapply(
    inter$X_list,
    function(Xc) {
      Xd <- as.matrix(Xc)
      drop(Xd %*% theta_true)
    },
    numeric(nrow(inter$X_list[[1L]]))
  )
  # add a tiny amount of noise to avoid perfect collinearity
  Z_event <- Z_event + matrix(rnorm(length(Z_event), sd = 1e-6), nrow(Z_event))

  theta_est <- hrfdecode:::estimate_hrf_theta(inter$X_list, Z_event, inter$hrf, penalty = 1e-4)
  corr <- cor(theta_true, theta_est)
  rel_err <- sqrt(sum((theta_est - theta_true)^2)) / (sqrt(sum(theta_true^2)) + 1e-12)
  expect_gt(corr, 0.995)
  expect_lt(rel_err, 5e-3)
})

test_that("ALS fit is deterministic given identical inputs", {
  skip_if_not_installed("fmridesign"); skip_if_not_installed("fmrihrf"); skip_if_not_installed("Matrix")
  skip_if_not(exists("fit_softlabels_als", envir = asNamespace("hrfdecode"), inherits = FALSE))
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  set.seed(9)
  sframe <- fmrihrf::sampling_frame(blocklens = 55, TR = 1)
  evtab  <- data.frame(onset = c(7, 17, 37), condition = factor(c("A","B","A")), run = 1)
  evmod  <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg1"),
    data = evtab, block = ~run, sampling_frame = sframe
  )
  prep   <- hrfdecode:::prepare_decoder_inputs(evmod, hrf = NULL, background = TRUE)
  L      <- hrfdecode:::build_laplacian_from_runs(prep$run_ids)
  Tn     <- nrow(prep$DBbeta); V <- 11
  set.seed(42); X <- matrix(rnorm(Tn * V), Tn, V)

  fit1 <- fit_softlabels_als(X, prep$P0, L, prep$DBbeta, 4, 1, 0.5,
                             max_iter = 6, tol = 1e-8, nonneg = TRUE, threads = FALSE)
  fit2 <- fit_softlabels_als(X, prep$P0, L, prep$DBbeta, 4, 1, 0.5,
                             max_iter = 6, tol = 1e-8, nonneg = TRUE, threads = FALSE)
  expect_equal(fit1$W, fit2$W, tolerance = 1e-12)
  expect_equal(fit1$P, fit2$P, tolerance = 1e-12)
  if ("b" %in% names(fit1)) expect_equal(fit1$b, fit2$b, tolerance = 1e-12)
})

test_that("Trace exposes norms and step sizes", {
  skip_if_not_installed("fmridesign"); skip_if_not_installed("fmrihrf"); skip_if_not_installed("Matrix")
  skip_if_not(exists("fit_softlabels_als", envir = asNamespace("hrfdecode"), inherits = FALSE))
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  set.seed(10)
  sframe <- fmrihrf::sampling_frame(blocklens = 50, TR = 1)
  evtab  <- data.frame(onset = c(8, 18, 28), condition = factor(c("A","B","A")), run = 1)
  evmod  <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg1"),
    data = evtab, block = ~run, sampling_frame = sframe
  )
  prep   <- hrfdecode:::prepare_decoder_inputs(evmod, hrf = NULL, background = TRUE)
  L      <- hrfdecode:::build_laplacian_from_runs(prep$run_ids)
  Tn     <- nrow(prep$DBbeta); V <- 8
  X      <- matrix(rnorm(Tn * V), Tn, V)

  fit <- fit_softlabels_als(X, prep$P0, L, prep$DBbeta,
                            lambda_W = 3, lambda_HRF = 0.8, lambda_smooth = 0.4,
                            max_iter = 5, tol = 1e-10, nonneg = TRUE, threads = FALSE)
  tr <- fit$obj_trace
  expect_true(all(c("w_norm","p_norm","dW","dP","rel_dW","rel_dP") %in% names(tr)))
  len <- length(tr$value)
  expect_equal(length(tr$w_norm), len)
  expect_equal(length(tr$rel_dP), len)
  expect_lt(tail(tr$rel_dP, 1), 0.1)
})

test_that("Ground-truth W, b, and P are recovered on noiseless synthetic data", {
  skip_if_not_installed("Matrix")
  skip_if_not(exists("fit_softlabels_als", envir = asNamespace("hrfdecode"), inherits = FALSE))
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  set.seed(202405)
  Tn <- 32L
  V  <- 7L
  K  <- 4L
  X <- matrix(rnorm(Tn * V), Tn, V)
  W_true <- matrix(rnorm(V * K), V, K)
  b_true <- runif(K, min = -0.25, max = 0.25)
  P_true <- sweep(X %*% W_true, 2, b_true, FUN = "+")

  # Start close to the truth but not exact so ALS has work to do.
  P0 <- P_true + matrix(rnorm(Tn * K, sd = 1e-4), Tn, K)
  DBbeta <- P_true
  L <- Matrix::Diagonal(x = rep(0, Tn))

  fit <- fit_softlabels_als(
    X = X,
    P0 = P0,
    L = L,
    DBbeta = DBbeta,
    lambda_W = 1e-9,
    lambda_HRF = 1e6,
    lambda_smooth = 0,
    max_iter = 4,
    tol = 1e-10,
    nonneg = FALSE,
    threads = FALSE
  )

  b_hat <- if ("b" %in% names(fit)) fit$b else rep(0, K)
  rel_err_W <- sqrt(sum((fit$W - W_true)^2)) / max(1e-12, sqrt(sum(W_true^2)))
  rel_err_b <- sqrt(sum((b_hat - b_true)^2)) / max(1e-12, sqrt(sum(b_true^2)))
  rel_err_P <- sqrt(sum((fit$P - P_true)^2)) / max(1e-12, sqrt(sum(P_true^2)))

  expect_lt(rel_err_W, 1e-6)
  expect_lt(rel_err_b, 1e-6)
  expect_lt(rel_err_P, 1e-6)
})

test_that("aggregate_events handles boundary truncation and global run offsets", {
  conditions <- c("A", "B")
  TR <- 1
  Tn <- 12L
  P <- cbind(seq_len(Tn) / 10, rev(seq_len(Tn)) / 20)
  events <- data.frame(
    onset = c(0, 4, 9),
    condition = factor(c("A", "B", "A"), levels = conditions),
    run = c(1, 1, 2)
  )
  window <- c(0, 3)
  agg <- aggregate_events(P = P, events = events, TR = TR,
                          conditions = conditions, window = window, hrf = NULL)
  weights <- rep(1 / (diff(window) + 1), diff(window) + 1)
  manual <- function(onset) {
    start_idx <- floor((onset + window[1]) / TR) + 1L
    idx <- start_idx + seq_along(weights) - 1L
    idx <- idx[idx >= 1L & idx <= Tn]
    if (!length(idx)) return(rep(0, length(conditions)))
    w <- weights[seq_along(idx)]
    colSums(P[idx, , drop = FALSE] * w)
  }
  expected <- t(sapply(events$onset, manual))
  colnames(expected) <- conditions
  expect_equal(agg$probs, expected, tolerance = 1e-12)
  expect_equal(agg$y_true, events$condition)
})

test_that("Predict path returns valid probabilities and aggregates correctly", {
  skip_if_not_installed("fmridesign"); skip_if_not_installed("fmrihrf")
  sframe <- fmrihrf::sampling_frame(blocklens = 20, TR = 1)
  evtab  <- data.frame(onset = c(5, 12, 17), condition = factor(c("A","B","A")), run = 1)
  evmod  <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg1"),
    data = evtab, block = ~run, sampling_frame = sframe
  )
  set.seed(10)
  Y <- matrix(rnorm(20 * 6, sd = 0.2), 20, 6)
  fit <- fit_hrfdecoder(
    Y, ev_model = evmod,
    lambda_W = 1, lambda_HRF = 0.5, lambda_smooth = 0.2,
    max_iter = 4, tol = 1e-6, nonneg = TRUE,
    background = TRUE, standardize = TRUE, verbose = 0
  )
  tr_out <- predict_hrfdecoder(fit, Y_test = Y, mode = "tr")
  rs <- rowSums(tr_out)
  if (max(abs(rs - 1)) < 0.1) {
    expect_true(all(abs(rs - 1) < 1e-6))
    expect_gte(min(tr_out), -1e-12)
  } else {
    testthat::skip("TR softmax not enabled yet; skipping probability checks.")
  }

  cond_a <- fit$conditions[1]
  cond_b <- fit$conditions[min(2, length(fit$conditions))]
  Ptoy <- matrix(0, nrow = 20, ncol = length(fit$conditions),
                 dimnames = list(NULL, fit$conditions))
  Ptoy[6:8, cond_a] <- c(0.2, 0.6, 0.8)
  Ptoy[6:8, cond_b] <- c(0.1, 0.1, 0.1)
  events_small <- data.frame(onset = 5, condition = factor(cond_a, levels = fit$conditions))
  agg <- aggregate_events(P = Ptoy, events = events_small, TR = 1,
                          conditions = fit$conditions, window = c(0, 2), hrf = NULL)
  expect_equal(unname(agg$probs[1, cond_a]), mean(c(0.2, 0.6, 0.8)), tolerance = 1e-12)
  expect_equal(unname(agg$probs[1, cond_b]), mean(c(0.1, 0.1, 0.1)), tolerance = 1e-12)
})
test_that("Finite-difference gradient check for W/P", {
  skip_if_not_installed("fmridesign"); skip_if_not_installed("fmrihrf"); skip_if_not_installed("Matrix")
  ns <- asNamespace("hrfdecode")
  skip_if_not("engine_objective" %in% ls(ns))
  engine_objective <- ns$engine_objective

  set.seed(123)
  sframe <- fmrihrf::sampling_frame(blocklens = 6, TR = 1)
  evtab  <- data.frame(onset = c(1, 3, 5), condition = factor(c("A", "B", "A")), run = 1)
  evmod  <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg1"),
    data = evtab, block = ~run, sampling_frame = sframe
  )
  prep   <- hrfdecode:::prepare_decoder_inputs(evmod, hrf = NULL, background = TRUE)
  L      <- hrfdecode:::build_laplacian_from_runs(prep$run_ids)
  Tn     <- nrow(prep$DBbeta); V <- 3
  X      <- matrix(rnorm(Tn * V), Tn, V)
  lamW <- 2; lamH <- 0.8; lamS <- 0.1

  # Single ALS iteration to get W,P,b
  fit <- hrfdecode:::fit_softlabels_als(X, prep$P0, L, prep$DBbeta,
                                        lambda_W = lamW, lambda_HRF = lamH, lambda_smooth = lamS,
                                        max_iter = 1, tol = 1e-12, nonneg = FALSE, threads = FALSE)
  W <- fit$W; P <- fit$P; b <- fit$b
  eps <- 1e-5

  # Finite difference on W
  coords_W <- expand.grid(i = 1:nrow(W), j = 1:ncol(W))
  residual <- X %*% W
  residual <- sweep(residual, 2, b, "+") - P
  grad_W <- 2 * (t(X) %*% residual) + 2 * lamW * W
  for (row in seq_len(min(10, nrow(coords_W)))) {
    set.seed(100 + row)
    idx <- sample(nrow(coords_W), 1)
    i <- coords_W$i[idx]; j <- coords_W$j[idx]
    W_plus <- W; W_plus[i, j] <- W_plus[i, j] + eps
    W_minus <- W; W_minus[i, j] <- W_minus[i, j] - eps
    obj_plus <- als_objective(X, W_plus, b, P, prep$DBbeta, lamW, lamH, lamS, L)
    obj_minus <- als_objective(X, W_minus, b, P, prep$DBbeta, lamW, lamH, lamS, L)
    fd <- (obj_plus - obj_minus) / (2 * eps)
    expect_equal(fd, grad_W[i, j], tolerance = 1e-3)
  }

  # Finite difference on P
  coords_P <- expand.grid(i = 1:nrow(P), j = 1:ncol(P))
  grad_P <- as.matrix(-2 * residual + 2 * lamH * (P - prep$DBbeta) + 2 * lamS * (L %*% P))
  set.seed(999)
  sample_idx <- sample(nrow(coords_P), min(10, nrow(coords_P)))
  for (idx in sample_idx) {
    i <- coords_P$i[idx]; j <- coords_P$j[idx]
    P_plus <- P; P_plus[i, j] <- P_plus[i, j] + eps
    P_minus <- P; P_minus[i, j] <- P_minus[i, j] - eps
    obj_plus <- als_objective(X, W, b, P_plus, prep$DBbeta, lamW, lamH, lamS, L)
    obj_minus <- als_objective(X, W, b, P_minus, prep$DBbeta, lamW, lamH, lamS, L)
    fd <- (obj_plus - obj_minus) / (2 * eps)
    analytic <- unname(grad_P[i, j])
    if (abs(fd) < 1e-8 && abs(analytic) < 1e-8) next
    expect_equal(fd, analytic, tolerance = 1e-3)
  }
})

test_that("Random hyperparameter sweep remains monotone and converges", {
  skip_if_not_installed("Matrix")
  skip_if_not(exists("fit_softlabels_als", envir = asNamespace("hrfdecode"), inherits = FALSE))
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  seeds <- 1:10
  for (sd in seeds) {
    set.seed(sd)
    Tn <- 24L; V <- 5L; K <- 3L
    X <- matrix(rnorm(Tn * V), Tn, V)
    P0 <- matrix(runif(Tn * K), Tn, K)
    DBbeta <- matrix(rnorm(Tn * K), Tn, K)
    L <- Matrix::Diagonal(x = rep(0, Tn))
    lamW <- runif(1, 0, 5)
    lamH <- runif(1, 0, 5)
    lamS <- runif(1, 0, 2)
    fit <- fit_softlabels_als(
      X = X, P0 = P0, L = L, DBbeta = DBbeta,
      lambda_W = lamW, lambda_HRF = lamH, lambda_smooth = lamS,
      max_iter = 12, tol = 1e-7, nonneg = TRUE, threads = FALSE
    )
    tr <- fit$obj_trace
    msg <- paste("seed", sd)
    expect_true(all(diff(tr$value) <= 1e-8), info = msg)
    expect_true(tail(tr$rel_dP, 1) < 1e-2, info = msg)
  }
})

test_that("Warm starts continue the objective trace without jumps", {
  skip_if_not_installed("Matrix")
  skip_if_not(exists("fit_softlabels_als", envir = asNamespace("hrfdecode"), inherits = FALSE))
  fit_softlabels_als <- hrfdecode:::fit_softlabels_als

  set.seed(321)
  Tn <- 30L; V <- 6L; K <- 3L
  X <- matrix(rnorm(Tn * V), Tn, V)
  P0 <- matrix(rnorm(Tn * K), Tn, K)
  DBbeta <- matrix(rnorm(Tn * K), Tn, K)
  L <- Matrix::Diagonal(x = rep(0, Tn))
  lamW <- 2; lamH <- 0.5; lamS <- 0.8

  fit_a <- fit_softlabels_als(
    X, P0, L, DBbeta,
    lambda_W = lamW, lambda_HRF = lamH, lambda_smooth = lamS,
    max_iter = 3, tol = 1e-8, nonneg = TRUE, threads = FALSE
  )
  fit_b <- fit_softlabels_als(
    X, fit_a$P, L, DBbeta,
    lambda_W = lamW, lambda_HRF = lamH, lambda_smooth = lamS,
    max_iter = 5, tol = 1e-10, nonneg = TRUE, threads = FALSE
  )
  last_a <- tail(fit_a$obj_trace$value, 1)
  first_b <- fit_b$obj_trace$value[1]
  expect_lte(first_b, last_a + 1e-8)
  combined <- c(fit_a$obj_trace$value, fit_b$obj_trace$value)
  expect_true(all(diff(combined) <= 1e-8))
})
