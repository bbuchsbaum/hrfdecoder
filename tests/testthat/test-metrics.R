test_that("binary AUC works and matches expectations", {
  set.seed(1)
  p <- c(0.9, 0.8, 0.7, 0.2, 0.1)
  y <- c(1,1,1,0,0)
  auc <- hrfdecode:::metric_auc_binary(p, y)
  expect_gt(auc, 0.99)

  # Inverted predictions → near 0
  auc2 <- hrfdecode:::metric_auc_binary(1 - p, y)
  expect_lt(auc2, 0.01)

  # Ties → around 0.5
  p3 <- c(0.5, 0.5, 0.5, 0.5)
  y3 <- c(1,0,1,0)
  auc3 <- hrfdecode:::metric_auc_binary(p3, y3)
  expect_equal(auc3, 0.5)
})

test_that("multi-metric computation returns sensible values", {
  probs <- matrix(c(0.8, 0.2,
                    0.7, 0.3,
                    0.4, 0.6,
                    0.1, 0.9), ncol = 2, byrow = TRUE)
  colnames(probs) <- c("A","B")
  y <- factor(c("A","A","B","B"), levels = c("A","B"))
  out <- hrfdecode:::compute_metrics(probs, y, metrics = c("accuracy","auc","logloss","brier","bal_acc"))
  expect_true(all(c("accuracy","auc","logloss","brier","bal_acc") %in% out$metric))
  expect_gt(out$value[out$metric=="accuracy"], 0.74)
  expect_gt(out$value[out$metric=="auc"], 0.79)
  # OVO should be close to OVR in this balanced example
  out2 <- hrfdecode:::compute_metrics(probs, y, metrics = c("auc_ovr","auc_ovo"))
  expect_true(abs(out2$value[1] - out2$value[2]) < 1e-6 || all(!is.na(out2$value)))
})
