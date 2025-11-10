test_that("AUC increases with stronger signal in simulation", {
  set.seed(123)
  n_trs <- 120; n_vox <- 20; n_trials <- 24; tr <- 2
  onsets <- seq(10, n_trs * tr - 20, length.out = n_trials)
  cond <- rep(c("A","B"), each = n_trials/2)
  events <- data.frame(onset = onsets, condition = cond, duration = 1)
  ev_model <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg1"),
    data = events, block = ~ 1,
    sampling_frame = fmrihrf::sampling_frame(TR = tr, blocklens = n_trs)
  )
  hrf_basis <- fmrihrf::evaluate(fmrihrf::getHRF("spmg2"), seq(0, 24, by = tr))
  h <- as.numeric(hrf_basis %*% c(1,0))
  stick_A <- rep(0, n_trs); stick_B <- rep(0, n_trs)
  idx_A <- pmin(n_trs, pmax(1L, floor(onsets[cond=="A"]/tr) + 1L))
  idx_B <- pmin(n_trs, pmax(1L, floor(onsets[cond=="B"]/tr) + 1L))
  stick_A[idx_A] <- 1; stick_B[idx_B] <- 1
  signal <- stick_A - stick_B
  sconv <- stats::convolve(signal, rev(h), type = "open")[1:n_trs]

  sim_dat <- function(scale_sig) {
    Y <- matrix(rnorm(n_trs * n_vox), n_trs, n_vox)
    for (v in 1:n_vox) Y[,v] <- Y[,v] + sconv * (v <= n_vox/2) * scale_sig
    fit <- fit_hrfdecoder(Y, ev_model, lambda_W = 0.1, max_iter = 8, verbose = 0)
    pred <- predict(fit, newdata = Y, ev_model_test = ev_model, mode = "trial")
    hrfdecode:::compute_metrics(pred$probs, pred$y_true, metrics = c("auc"))$value
  }
  auc_low <- sim_dat(0.3)
  auc_high <- sim_dat(0.8)
  expect_gt(auc_high, auc_low + 0.1)
})

