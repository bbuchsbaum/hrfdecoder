test_that("AR(1) prewhitening improves AUC with AR(1) noise", {
  set.seed(42)
  n_trs <- 160; n_vox <- 24; n_trials <- 32; tr <- 2
  onsets <- seq(10, n_trs * tr - 20, length.out = n_trials)
  cond <- rep(c("A","B"), each = n_trials/2)
  events <- data.frame(onset=onsets, condition=cond, duration=1)
  ev <- fmridesign::event_model(onset ~ fmridesign::hrf(condition, basis='spmg1'),
                                data=events, block=~1,
                                sampling_frame=fmrihrf::sampling_frame(TR=tr, blocklens=n_trs))
  hrfb <- fmrihrf::evaluate(fmrihrf::getHRF('spmg2'), seq(0,24,by=tr))
  h <- as.numeric(hrfb %*% c(1,0))
  stick_A <- rep(0, n_trs); stick_B <- rep(0, n_trs)
  idx_A <- pmin(n_trs, pmax(1L, floor(onsets[cond=='A']/tr) + 1L))
  idx_B <- pmin(n_trs, pmax(1L, floor(onsets[cond=='B']/tr) + 1L))
  stick_A[idx_A] <- 1; stick_B[idx_B] <- 1
  sconv <- stats::convolve(stick_A - stick_B, rev(h), type='open')[1:n_trs]

  # AR(1) noise
  phi <- 0.6
  Y <- matrix(0, n_trs, n_vox)
  for (v in 1:n_vox) {
    eps <- rnorm(n_trs)
    for (t in 2:n_trs) eps[t] <- phi * eps[t-1] + rnorm(1)
    Y[,v] <- eps + sconv * (v <= n_vox/2) * 0.6
  }

  # No AR vs AR(1)
  fit_none <- fit_hrfdecoder(Y, ev, ar_order = NULL, lambda_W=0.1, max_iter=8, verbose=0)
  fit_ar1  <- fit_hrfdecoder(Y, ev, ar_order = 1,    lambda_W=0.1, max_iter=8, verbose=0)
  pred_none <- predict(fit_none, newdata=Y, ev_model_test=ev, mode='trial')
  pred_ar1  <- predict(fit_ar1,  newdata=Y, ev_model_test=ev, mode='trial')
  auc_none <- hrfdecode:::metric_auc(pred_none$probs, pred_none$y_true)
  auc_ar1  <- hrfdecode:::metric_auc(pred_ar1$probs,  pred_ar1$y_true)
  expect_gt(auc_ar1, auc_none + 0.05)
})

