test_that("K no-test scalar stays numeric and derives canonical estimate", {
  x <- structure(0.42, class = "phylosig", method = "K", test = FALSE)
  fit <- fastphylosig:::.decorate_fastphylosig_result(
    x, method = "K", vector_input = TRUE
  )

  expect_true(is.numeric(fit))
  expect_equal(
    class(fit), c("fastphylosig_signal", "fastphylosig_result", "phylosig")
  )
  expect_equal(as.data.frame(fit)$estimate, 0.42)
  expect_equal(summary(fit)$estimate, 0.42)
})

test_that("signal list aliases and classes preserve legacy fields", {
  fit <- fastphylosig:::.decorate_fastphylosig_result(
    structure(list(K = 0.3, P = 0.2, MCSE_P = 0.01, nsim = 20,
                   nsim_successful = 18, sim.K = 1:4),
              class = "phylosig", method = "K"),
    method = "K", vector_input = TRUE
  )

  expect_s3_class(fit, "fastphylosig_signal")
  expect_s3_class(fit, "fastphylosig_result")
  expect_s3_class(fit, "phylosig")
  expect_equal(fit$estimate, fit$K)
  expect_equal(fit$p_value, fit$P)
  expect_equal(fit$p_mcse, fit$MCSE_P)
  expect_equal(fit$n_sim_requested, fit$nsim)
  expect_equal(fit$n_sim_successful, fit$nsim_successful)
  expect_true("sim.K" %in% names(fit))
  expect_false(grepl("sim.K", paste(capture.output(print(fit)), collapse = "\n")))
})

test_that("lambda and D aliases remain statistic-specific", {
  lambda <- fastphylosig:::.decorate_fastphylosig_result(
    structure(list(lambda = 0.7, logL = -4, LR = 2, P = 0.1,
                   legacy = "kept"), class = "phylosig", method = "lambda"),
    method = "lambda", vector_input = TRUE
  )
  expect_equal(lambda$estimate, lambda$lambda)
  expect_equal(lambda$logLik, lambda$logL)
  expect_equal(lambda$p_value, lambda$P)
  expect_equal(lambda$legacy, "kept")

  d <- fastphylosig:::.decorate_fastphylosig_result(
    structure(list(DEstimate = 0.5, Pval1 = 0.2, Pval0 = 0.8,
                   MCSE_P_random = 0.1, MCSE_P_Brownian = 0.15),
              class = "phylo.d"),
    method = "D", vector_input = TRUE
  )
  expect_equal(d$estimate, d$DEstimate)
  expect_equal(d$P_random, d$Pval1)
  expect_equal(d$P_Brownian, d$Pval0)
  expect_equal(d$p_mcse_random, d$MCSE_P_random)
  expect_equal(d$p_mcse_brownian, d$MCSE_P_Brownian)
  expect_false("p_value" %in% names(d))
})

test_that("Delta aliases include uncertainty and diagnostics counts", {
  fit <- fastphylosig:::.decorate_fastphylosig_result(
    structure(list(delta = 0.6, P = 0.25, P_MCSE = 0.05,
                   MCSE_Delta = 0.02, ESS_alpha = 30, ESS_beta = 25,
                   Rhat_alpha = 1.01, Rhat_beta = 1.02,
                   requested_simulations = 10,
                   successful_simulations = 9),
              class = "phylo_delta"),
    method = "Delta", vector_input = TRUE
  )
  expect_equal(fit$estimate, fit$delta)
  expect_equal(fit$p_value, fit$P)
  expect_equal(fit$p_mcse, fit$P_MCSE)
  expect_equal(fit$estimate_mcse, fit$MCSE_Delta)
  expect_equal(fit$ess_alpha, fit$ESS_alpha)
  expect_equal(fit$rhat_beta, fit$Rhat_beta)
  expect_equal(fit$ESS, 25)
  expect_equal(fit$Rhat, 1.02)
  expect_equal(fit$n_sim_requested, 10)
  expect_equal(fit$n_sim_successful, 9)
})

test_that("batch decoration adds canonical columns and compact coercion", {
  tab <- data.frame(
    trait = "x", K_fast = 0.4, P_fast = 0.3, MCSE_P = 0.1,
    nsim_requested = 10L, nsim_successful = 9L,
    sim.K_fast = I(list(1:10)), stringsAsFactors = FALSE
  )
  fit <- fastphylosig:::.decorate_fastphylosig_result(
    tab, method = "K", vector_input = FALSE
  )
  expect_equal(
    class(fit), c("fastphylosig_table", "fastphylosig_result", "data.frame")
  )
  expect_equal(fit$estimate, fit$K_fast)
  expect_equal(fit$p_value, fit$P_fast)
  expect_equal(fit$p_mcse, fit$MCSE_P)
  expect_equal(fit$n_sim_requested, fit$nsim_requested)
  expect_equal(fit$n_sim_successful, fit$nsim_successful)
  compact <- as.data.frame(fit)
  expect_false("sim.K_fast" %in% names(compact))
  expect_true(all(c("estimate", "p_value", "p_mcse") %in% names(compact)))
})
