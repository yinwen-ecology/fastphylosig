test_that("lowercase D probabilities select their matching null distributions", {
  for (table in c(FALSE, TRUE)) {
    if (table) {
      fit <- data.frame(trait = "x", D_fast = 0.5, observed = 2,
                        mean_random = 3, mean_brownian = 1,
                        p_random = 0.4, p_brownian = 0.2)
      fit$random_fast <- I(list(c(0, 1, 2, 2, 3)))
      fit$brownian_fast <- I(list(c(1, 2, 2, 3, 4)))
    } else {
      fit <- structure(list(DEstimate = 0.5, binvar = "x",
        Parameters = list(Observed = 2, MeanRandom = 3, MeanBrownian = 1),
        Permutations = list(random = c(0, 1, 2, 2, 3),
                            brownian = c(1, 2, 2, 3, 4)),
        p_random = 0.4, p_brownian = 0.2), class = "phylo.d")
    }
    for (null in c("random", "brownian")) {
      alias <- paste0("p_", null)
      expected_p <- if (null == "random") 0.4 else 0.2
      expected_sim <- (if (null == "random") c(0, 1, 2, 2, 3) else
        c(1, 2, 2, 3, 4)) / 2 - 0.5
      dat <- fastphylosig:::.signal_distribution_data(fit, p_col = alias)
      expect_equal(dat$p_value, expected_p)
      expect_identical(dat$null_model, paste(null, "D"))
      expect_equal(dat$sim[[1]], expected_sim)
      expect_equal(fastphylosig:::.signal_plot_data(fit, p_col = alias)$p_value,
                   expected_p)
      # Explicit null selection wins over an opposite lowercase alias.
      opposite <- paste0("p_", if (null == "random") "brownian" else "random")
      dat <- fastphylosig:::.signal_distribution_data(fit, null = null,
                                                    p_col = opposite)
      expect_equal(dat$p_value, expected_p)
      expect_equal(dat$sim[[1]], expected_sim)
    }
    overlay <- fastphylosig:::.signal_d_overlay_data(fit)
    expect_equal(overlay$Pval1, 0.4)
    expect_equal(overlay$Pval0, 0.2)
  }
})

test_that("one retained iteration preserves separate Delta chains", {
  alpha <- matrix(c(1, 3), nrow = 1)
  beta <- matrix(c(2, 12), nrow = 1)
  for (x in list(alpha, as.data.frame(alpha), list(1, 3))) {
    expect_equal(fastphylosig:::.delta_chain_matrix(x), alpha,
                 ignore_attr = TRUE)
    expect_identical(dim(fastphylosig:::.delta_chain_matrix(x)), c(1L, 2L))
  }
  fallback <- fastphylosig:::.delta_diagnostics_fallback(alpha, beta)
  native <- fastphylosig:::delta_diagnostics_cpp(alpha, beta)
  dispatched <- fastphylosig:::.delta_mcmc_diagnostics(
    list(alpha_chain = alpha, beta_chain = beta))
  expect_equal(fallback$n_iter, 1)
  expect_equal(fallback$n_chains, 2)
  expect_equal(fallback$chain_alpha_mean, c(1, 3))
  expect_equal(fallback$chain_beta_mean, c(2, 12))
  for (nm in c("n_iter", "n_chains", "n_saved", "chain_alpha_mean",
               "chain_beta_mean", "MCSE_Delta", "ESS_alpha", "ESS_beta")) {
    expect_equal(as.numeric(fallback[[nm]]), as.numeric(native[[nm]]), info = nm)
    expect_equal(as.numeric(dispatched[[nm]]), as.numeric(native[[nm]]), info = nm)
  }
  expect_true(is.na(fallback$split_Rhat_Delta))
  expect_equal(fastphylosig:::.delta_chain_matrix(matrix(1:20, 2, 10)),
               t(matrix(as.numeric(1:20), 2, 10)))
})
