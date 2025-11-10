#' Fit HRF-aware weakly supervised decoder
#'
#' @param Y numeric matrix (T x V) of fMRI data (time by voxel)
#' @param ev_model fmridesign::event_model describing events
#' @param base_model optional fmridesign::baseline_model for nuisance removal
#' @param hrf optional fmrihrf basis object
#' @param lambda_W ridge penalty on decoder weights
#' @param lambda_HRF adherence weight to HRF prior
#' @param lambda_smooth temporal smoothness weight
#' @param theta_penalty ridge penalty on HRF basis coefficients
#' @param max_iter maximum ALS iterations
#' @param tol convergence tolerance on P updates
#' @param nonneg enforce non-negative soft labels
#' @param background include a background column in the soft labels
#' @param standardize z-score Y columns before fitting
#' @param ar_order AR order for prewhitening (default: NULL for no AR prewhitening).
#'   Set to 1 for AR(1), 2 for AR(2), etc. Use "auto" for automatic BIC-based selection.
#' @param ar_method AR estimation method: "ar" (Yule-Walker) or "arma" (Hannan-Rissanen).
#'   Default: "ar".
#' @param ar_pooling Spatial pooling for AR parameters: "global" (one AR model for all voxels)
#'   or "run" (separate AR model per run). Default: "run".
#' @param verbose integer verbosity level
#'
#' @return object of class `hrfdecoder_fit`
#' @export
fit_hrfdecoder <- function(
  Y,
  ev_model,
  base_model = NULL,
  hrf = NULL,
  lambda_W = 10,
  lambda_HRF = 1,
  lambda_smooth = 5,
  theta_penalty = 0.01,
  max_iter = 20,
  tol = 1e-4,
  nonneg = TRUE,
  background = TRUE,
  standardize = TRUE,
  ar_order = NULL,
  ar_method = c("ar", "arma"),
  ar_pooling = c("run", "global"),
  verbose = 1
) {
  stopifnot(is.matrix(Y))
  stopifnot(inherits(ev_model, "event_model"))
  Tn <- nrow(Y)
  ar_method <- match.arg(ar_method)
  ar_pooling <- match.arg(ar_pooling)

  # STEP 1: Baseline residualization
  Y <- residualize_baseline(Y, base_model)

  # STEP 2: AR prewhitening (if requested)
  ar_plan <- NULL
  if (!is.null(ar_order) && (is.character(ar_order) || ar_order > 0)) {
    # Need run_ids for AR estimation - derive directly from sampling frame
    sframe <- .get_sframe(ev_model)
    blocklens <- sframe$blocklens
    run_ids <- rep(seq_along(blocklens), blocklens)

    # Estimate AR model from residuals (post-baseline, pre-standardization)
    ar_plan <- fmriAR::fit_noise(
      Y,
      runs = run_ids,
      method = ar_method,
      p = ar_order,
      pooling = ar_pooling,
      exact_first = "ar1"
    )

    # Apply whitening transformation (pass a dummy X matrix)
    Y_white <- fmriAR::whiten_apply(ar_plan, X = matrix(0, nrow(Y), 1L), Y = Y, runs = run_ids)
    Y <- Y_white$Y

    if (verbose >= 1) {
      p_actual <- if (is.character(ar_order)) ar_plan$order["p"] else ar_order
      message("Applied AR(", p_actual, ") prewhitening (pooling=", ar_pooling, ")")
    }
  }

  # STEP 3: Standardization
  preproc <- list(
    standardize = isTRUE(standardize),
    center = NULL,
    scale = NULL,
    ar_plan = ar_plan
  )
  if (isTRUE(standardize)) {
    Ys <- scale(Y)
    preproc$center <- attr(Ys, "scaled:center")
    s <- attr(Ys, "scaled:scale")
    s[is.na(s) | s == 0] <- 1
    preproc$scale <- s
    attr(Ys, "scaled:center") <- NULL
    attr(Ys, "scaled:scale") <- NULL
    Y <- Ys
  }

  prep <- prepare_decoder_inputs(ev_model, hrf = hrf, background = background)
  if (nrow(prep$DBbeta) != Tn) {
    stop("Y rows do not match the event model length.")
  }
  L <- build_laplacian_from_runs(prep$run_ids)

  if (verbose) {
    message("Running alternating solver (", Tn, " TRs, ", ncol(Y), " voxels).")
  }

  fit_cpp <- fit_softlabels_als(
    X = Y,
    P0 = prep$P0,
    L = L,
    DBbeta = prep$DBbeta,
    lambda_W = lambda_W,
    lambda_HRF = lambda_HRF,
    lambda_smooth = lambda_smooth,
    max_iter = max_iter,
    tol = tol,
    nonneg = nonneg,
    threads = FALSE
  )

  Z <- fit_cpp$P
  Z_event <- Z[, seq_along(prep$X_list), drop = FALSE]
  theta <- estimate_hrf_theta(prep$X_list, Z_event, prep$hrf, penalty = theta_penalty)
  hrf_est <- fmrihrf::hrf_from_coefficients(prep$hrf, theta)

  TR <- sampling_frame_tr(prep$sframe)
  events_tbl <- get_event_table(ev_model)

  fit <- list(
    W = fit_cpp$W,
    P = Z,
    b = as.numeric(fit_cpp$b),
    theta = theta,
    hrf = hrf_est,
    conditions = prep$conditions,
    background = background,
    converged = isTRUE(fit_cpp$converged) || (fit_cpp$iterations >= max_iter),
    iterations = fit_cpp$iterations,
    settings = list(
      lambda_W = lambda_W,
      lambda_HRF = lambda_HRF,
      lambda_smooth = lambda_smooth,
      theta_penalty = theta_penalty,
      max_iter = max_iter,
      tol = tol,
      nonneg = nonneg,
      background = background,
      TR = TR,
      run_ids = prep$run_ids,
      ar_order = ar_order,
      ar_method = ar_method,
      ar_pooling = ar_pooling
    ),
    train = list(
      P0 = prep$P0,
      prior = prep$DBbeta,
      events = events_tbl
    ),
    preproc = preproc,
    diagnostics = list(
      obj_trace = fit_cpp$obj_trace
    )
  )
  class(fit) <- "hrfdecoder_fit"
  fit
}
