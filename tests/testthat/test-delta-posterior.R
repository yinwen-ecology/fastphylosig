test_that("Delta and its MCSE summarize the nonlinear ratio trace", {
  a <- cbind(seq(1, 5, length.out = 100), seq(2, 7, length.out = 100))
  b <- a^2 + .1
  ratio <- b / a
  batch <- floor(sqrt(nrow(a))); nb <- nrow(a) %/% batch
  variance <- sum(vapply(1:2, function(j) {
    means <- colMeans(matrix(ratio[seq_len(nb * batch), j], batch, nb))
    var(means) / nb
  }, numeric(1))) / 4
  cpp <- delta_diagnostics_cpp(a, b)
  r <- .delta_diagnostics_fallback(a, b)
  expect_equal(cpp$MCSE_Delta, sqrt(variance), tolerance = 1e-12)
  expect_equal(r$MCSE_Delta, sqrt(variance), tolerance = 1e-12)
  expect_equal(cpp$ESS_Delta, r$ESS_Delta, tolerance = 1e-10)
  expect_equal(cpp$split_Rhat_Delta, r$split_Rhat_Delta, tolerance = 1e-12)
  expect_gt(abs(mean(ratio) - mean(b) / mean(a)), .1)
  constant <- delta_diagnostics_cpp(a, 2 * a)
  expect_identical(constant$MCSE_Delta, 0)
})

test_that("Delta rejects invalid prior and proposal controls", {
  p <- matrix(c(.2, .8, .4, .6), 2, 2, byrow = TRUE)
  for (bad in c(0, -1, Inf, NA_real_)) {
    expect_error(delta_mcmc_cpp(p, bad, .5, 10L, 1L, 2L, 1L), "finite and positive")
    expect_error(delta_mcmc_cpp(p, .1, bad, 10L, 1L, 2L, 1L), "finite and positive")
  }
})
