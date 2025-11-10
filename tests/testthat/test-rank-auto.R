skip_rank_tests <- function() {
  skip_if_not_installed("fmriAR")
  ns <- asNamespace("hrfdecode")
  skip_if_not("estimate_rank_auto" %in% ls(ns))
  ns$estimate_rank_auto
}

test_that("estimate_rank_auto recovers planted rank within ±1", {
  estimate_rank_auto <- skip_rank_tests()
  set.seed(2601)
  Tlen <- 220L
  V <- 260L
  r_true <- 5L
  U_true <- qr.Q(qr(matrix(rnorm(Tlen * r_true), Tlen, r_true)))
  R_true <- matrix(rnorm(V * r_true), V, r_true)
  Y <- U_true %*% t(R_true) + matrix(rnorm(Tlen * V, sd = 0.6), Tlen, V)
  X_list <- list(matrix(0, Tlen, 1))
  est <- estimate_rank_auto(Y, X_list, runs = rep(1L, Tlen),
                            alpha = 0.05, B = 20L, r_max = 15L)
  expect_gte(est$r, r_true - 1)
  expect_lte(est$r, r_true + 1)
})

test_that("estimate_rank_auto stays near zero under null residuals", {
  estimate_rank_auto <- skip_rank_tests()
  set.seed(2701)
  Tlen <- 220L
  V <- 280L
  Y <- matrix(rnorm(Tlen * V), Tlen, V)
  X_list <- list(matrix(0, Tlen, 1))
  est <- estimate_rank_auto(Y, X_list, runs = rep(1L, Tlen),
                            alpha = 0.05, B = 20L, r_max = 15L)
  expect_lte(est$r, 1L)
})

test_that("AR(1) noise still yields correct rank after whitening", {
  estimate_rank_auto <- skip_rank_tests()
  set.seed(2801)
  Tlen <- 240L
  V <- 220L
  r_true <- 4L
  phi <- 0.6
  U_true <- qr.Q(qr(matrix(rnorm(Tlen * r_true), Tlen, r_true)))
  R_true <- matrix(rnorm(V * r_true), V, r_true)
  signal <- U_true %*% t(R_true)
  noise <- matrix(0, Tlen, V)
  innov <- matrix(rnorm(Tlen * V, sd = 1), Tlen, V)
  for (t in 2:Tlen) {
    noise[t, ] <- phi * noise[t - 1L, ] + innov[t, ]
  }
  Y <- signal + 0.4 * noise
  X_list <- list(matrix(0, Tlen, 1))
  est <- estimate_rank_auto(Y, X_list, runs = rep(1L, Tlen),
                            alpha = 0.05, B = 20L, r_max = 15L)
  expect_gte(est$r, r_true - 1)
  expect_lte(est$r, r_true + 1)
})

test_that("Baseline removal does not distort rank selection", {
  estimate_rank_auto <- skip_rank_tests()
  set.seed(2901)
  Tlen <- 200L
  V <- 240L
  r_true <- 3L
  U_true <- qr.Q(qr(matrix(rnorm(Tlen * r_true), Tlen, r_true)))
  R_true <- matrix(rnorm(V * r_true), V, r_true)
  Y <- U_true %*% t(R_true) + matrix(rnorm(Tlen * V, sd = 0.5), Tlen, V)
  Kb <- 6L
  t <- seq_len(Tlen)
  X_base <- scale(sapply(seq_len(Kb), function(k) cos(2 * pi * k * t / Tlen)),
                  center = TRUE, scale = FALSE)
  X_list <- list(matrix(0, Tlen, 1))
  est <- estimate_rank_auto(
    Y, X_list, X_base = X_base, runs = rep(1L, Tlen),
    control = list(lambda_base = 1e-3),
    alpha = 0.05, B = 20L, r_max = 12L
  )
  expect_gte(est$r, r_true - 1)
  expect_lte(est$r, r_true + 1)
})

test_that("Selected rank increases with SNR", {
  estimate_rank_auto <- skip_rank_tests()
  set.seed(3001)
  Tlen <- 200L
  V <- 230L
  r_true <- 6L
  U_true <- qr.Q(qr(matrix(rnorm(Tlen * r_true), Tlen, r_true)))
  R_true <- matrix(rnorm(V * r_true), V, r_true)
  base_signal <- U_true %*% t(R_true)
  X_list <- list(matrix(0, Tlen, 1))
  Y_low <- 0.6 * base_signal + matrix(rnorm(Tlen * V, sd = 1.2), Tlen, V)
  Y_high <- 1.2 * base_signal + matrix(rnorm(Tlen * V, sd = 0.6), Tlen, V)
  est_low <- estimate_rank_auto(Y_low, X_list, runs = rep(1L, Tlen),
                                alpha = 0.05, B = 15L, r_max = 15L)
  est_high <- estimate_rank_auto(Y_high, X_list, runs = rep(1L, Tlen),
                                 alpha = 0.05, B = 15L, r_max = 15L)
  expect_lte(est_low$r, est_high$r)
  expect_gte(est_high$r, r_true - 1L)
})
