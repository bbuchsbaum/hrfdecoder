test_that("fit_hrfdecoder runs on synthetic data", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")
  sframe <- fmrihrf::sampling_frame(blocklens = 60, TR = 1)
  evtab <- data.frame(
    onset = c(5, 15, 35),
    condition = factor(c("A", "B", "A")),
    run = c(1, 1, 1)
  )
  evmod <- fmridesign::event_model(
    onset ~ fmridesign::hrf(condition, basis = "spmg1"),
    data = evtab,
    block = ~run,
    sampling_frame = sframe
  )
  Y <- matrix(rnorm(60 * 20), 60, 20)
  fit <- fit_hrfdecoder(
    Y = Y,
    ev_model = evmod,
    lambda_W = 1,
    lambda_HRF = 0.5,
    lambda_smooth = 0.1,
    max_iter = 5,
    verbose = 0
  )
  expect_equal(nrow(fit$P), 60)
  preds <- predict_hrfdecoder(fit, Y, ev_model_test = evmod, mode = "trial")
  expect_equal(ncol(preds$probs), length(fit$conditions))
})
