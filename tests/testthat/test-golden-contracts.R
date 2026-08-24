# Golden contracts for the public fastphylosig API --------------------------
#
# These fixtures are deliberately small and deterministic.  They are not a
# replacement for the broader reference-package tests in test-fast-signal.R;
# instead they provide a compact release gate for values and object contracts
# that should not move during a performance refactor.

.golden_tree <- function() {
  ape::read.tree(
    text = "((a:0.5,b:0.7):0.8,(c:0.4,d:0.6):0.9);"
  )
}

.golden_continuous <- function() {
  out <- cbind(
    trait_a = c(a = 1.2, b = -0.7, c = 2.3, d = 0.4),
    trait_b = c(a = -2, b = 0.3, c = 1.1, d = 3.5)
  )
  out
}

.golden_k_permutations <- function() {
  rbind(
    c(1L, 2L, 3L, 4L),
    c(2L, 1L, 4L, 3L),
    c(3L, 4L, 1L, 2L),
    c(4L, 3L, 2L, 1L)
  )
}

.golden_d_states <- function() {
  list(
    random = cbind(
      c(0, 1, 0, 1), c(1, 0, 1, 0), c(0, 0, 1, 1),
      c(1, 1, 0, 0), c(0, 1, 1, 0)
    ),
    brownian = cbind(
      c(0, 0, 1, 0), c(1, 1, 0, 1), c(0, 1, 0, 0),
      c(1, 0, 1, 1), c(0, 0, 1, 1)
    )
  )
}

.expect_delta_mcmc_contract <- function(fit) {
  expect_s3_class(fit, "phylo_delta")
  expect_true(is.finite(fit$delta))
  expect_true(is.finite(fit$alpha_mean))
  expect_true(is.finite(fit$beta_mean))
  expect_true(is.finite(fit$n_saved))
  expect_gte(fit$n_saved, 1)
  expect_true(is.list(fit$parameters))

  # Newer builds may expose Monte Carlo uncertainty and chain diagnostics.
  # Keep this contract forward-compatible without requiring fields that are
  # not part of the current public API.
  uncertainty <- intersect(
    c("mc_se", "MCSE", "mcse", "delta_mc_se", "MonteCarloSE",
      "MCSE_Delta"),
    names(fit)
  )
  if (length(uncertainty)) {
    expect_true(all(vapply(fit[uncertainty], is.numeric, logical(1))))
    expect_true(all(vapply(fit[uncertainty], function(z) {
      length(z) == 1L && is.finite(z)
    }, logical(1))))
  }
  chain_fields <- intersect(
    c("alpha_sd", "beta_sd", "ESS_alpha", "ESS_beta",
      "split_Rhat_alpha", "split_Rhat_beta", "alpha_beta_cov"),
    names(fit)
  )
  if (length(chain_fields)) {
    expect_true(all(vapply(fit[chain_fields], is.numeric, logical(1))))
    expect_true(all(vapply(fit[chain_fields], function(z) {
      length(z) == 1L && (is.na(z) || is.finite(z))
    }, logical(1))))
  }
  schema_fields <- intersect(
    c("n_saved_requested", "n_saved_successful", "requested_simulations",
      "successful_simulations", "diagnostics_available"), names(fit)
  )
  if (length(schema_fields)) {
    expect_true(all(vapply(fit[schema_fields], function(z) {
      length(z) == 1L
    }, logical(1))))
  }
  diagnostics <- intersect(
    c("diagnostics", "mcmc_diagnostics", "chain_diagnostics",
      "diagnostics_note", "note"),
    names(fit)
  )
  if (length(diagnostics)) {
    expect_true(all(vapply(fit[diagnostics], function(z) {
      is.list(z) || is.character(z) || is.logical(z)
    }, logical(1))))
  }
  invisible(fit)
}

test_that("golden K batch values and schema are stable", {
  tree <- .golden_tree()
  X <- .golden_continuous()
  fit <- suppressMessages(
    fast_signal(tree, X, method = "K", test = FALSE, verbose = FALSE)
  )

  expect_true(all(c("trait", "K_fast", "n_species", "n_removed_na",
                    "matched_species", "removed_tree_tips",
                    "removed_data_rows", "note") %in% names(fit)))
  expect_equal(fit$trait, c("trait_a", "trait_b"))
  expect_equal(fit$n_species, c(4L, 4L))
  expect_equal(fit$K_fast, c(0.609682598152122, 0.965219909161572),
               tolerance = 1e-8)
  expect_equal(attr(fit, "match_report")$matched_species, 4L)
})

test_that("golden K algebra freezes GLS and quadratic components", {
  tree <- .golden_tree()
  x <- .golden_continuous()[, "trait_a"]
  C <- ape::vcv.phylo(tree)
  invC <- solve(C)
  gls_mean <- sum(invC %*% x) / sum(invC)
  centered <- x - gls_mean
  numerator <- sum(centered^2)
  denominator <- drop(t(centered) %*% invC %*% centered)
  normalization <- (sum(diag(C)) - length(x) / sum(invC)) /
    (length(x) - 1)

  expect_equal(gls_mean, 0.961911874533234, tolerance = 1e-12)
  expect_equal(numerator, 4.92486182045946, tolerance = 1e-12)
  expect_equal(denominator, 7.19219566840926, tolerance = 1e-12)
  expect_equal(sum(invC), 1.79322351680729, tolerance = 1e-12)
  expect_equal(normalization, 1.12312671147623, tolerance = 1e-12)
  expect_equal((numerator / denominator) / normalization,
               0.609682598152122, tolerance = 1e-12)
})

test_that("golden K permutation values are stable and thread invariant", {
  tree <- .golden_tree()
  x <- .golden_continuous()[, "trait_a"]
  perms <- .golden_k_permutations()

  one <- suppressMessages(fast_signal(
    tree, x, method = "K", test = TRUE, nsim = nrow(perms),
    permutations = perms, return_sim = TRUE, verbose = FALSE, ncores = 1
  ))
  two <- suppressMessages(fast_signal(
    tree, x, method = "K", test = TRUE, nsim = nrow(perms),
    permutations = perms, return_sim = TRUE, verbose = FALSE, ncores = 2
  ))

  expect_equal(one$K, 0.609682598152122, tolerance = 1e-8)
  # phytools uses the inclusive upper tail, so permutations tied with the
  # observed K count as exceedances. The K kernel treats only machine-level
  # cross-platform rounding as a tie; this keeps the contract stable without
  # changing materially different permutation values.
  expect_equal(one$P, 1, tolerance = 1e-12)
  expect_equal(one$exceedance_count, 4)
  expect_equal(one$sim.K,
               c(0.609682598152122, 0.619158721940097,
                 0.619158721940097, 0.609682598152122),
               tolerance = 1e-8)
  expect_equal(one$K, two$K, tolerance = 1e-12)
  expect_equal(one$P, two$P, tolerance = 1e-12)
  expect_equal(one$exceedance_count, two$exceedance_count)
  expect_equal(one$sim.K, two$sim.K, tolerance = 1e-12)
})

test_that("golden lambda likelihood and profile contract are stable", {
  fit <- suppressMessages(fast_signal(
    .golden_tree(), .golden_continuous()[, "trait_a"], method = "lambda",
    test = TRUE, lambda_profile_points = 11L, verbose = FALSE
  ))

  expect_s3_class(fit, "phylosig")
  expect_equal(attr(fit, "method"), "lambda")
  expect_equal(fit$lambda, 6.61069613518961e-05, tolerance = 1e-5)
  expect_equal(fit$logL, -6.04621136219611, tolerance = 1e-5)
  expect_equal(fit$logL0, -6.04617125212792, tolerance = 1e-5)
  expect_equal(fit$P, 1, tolerance = 1e-12)
  C <- ape::vcv.phylo(.golden_tree())
  y <- .golden_continuous()[, "trait_a"]
  Cl <- fit$lambda * (C - diag(diag(C))) + diag(diag(C))
  invCl <- solve(Cl)
  gls_mean <- sum(invCl %*% y) / sum(invCl)
  centered <- y - gls_mean
  sigma2 <- drop(t(centered) %*% invCl %*% centered) / length(y)
  expect_equal(gls_mean, 0.867858568028865, tolerance = 1e-10)
  expect_equal(sigma2, 0.861839629072312, tolerance = 1e-10)
  expect_equal(fastphylosig:::.max_lambda(.golden_tree()), 1,
               tolerance = 1e-12)
  expect_equal(max(0, 2 * (fit$logL - fit$logL0)), 0,
               tolerance = 1e-12)
  expect_true(is.data.frame(fit$lambda_profile))
  expect_true(all(c("lambda", "logL") %in% names(fit$lambda_profile)))
  expect_true(any(abs(fit$lambda_profile$lambda) < 1e-12))
  expect_true(any(abs(fit$lambda_profile$lambda - 1) < 1e-12))
  expect_equal(fit$lambda_profile$logL[
    which.min(abs(fit$lambda_profile$lambda - fit$lambda))
  ], fit$logL, tolerance = 1e-10)
})

test_that("prepared-tree K path is a numerical golden equivalent", {
  tree <- .golden_tree()
  X <- .golden_continuous()
  ctx <- prepare_tree(tree, max_cached_subsets = 3L)
  expect_s3_class(ctx, "fastphylosig_tree")

  raw <- suppressMessages(
    fast_signal(tree, X, method = "K", test = FALSE, verbose = FALSE)
  )
  cached <- suppressMessages(
    fast_signal(ctx, X, method = "K", test = FALSE, verbose = FALSE)
  )
  expect_equal(cached$K_fast, raw$K_fast, tolerance = 1e-12)
  expect_equal(cached$n_species, raw$n_species)
  expect_equal(cached$matched_species, raw$matched_species)
  expect_true(length(ls(ctx$cache, all.names = TRUE)) >= 1L)
})

test_that("golden D calibration uses explicit null matrices", {
  states <- .golden_d_states()
  x <- c(a = 0, b = 0, c = 1, d = 1)
  fit <- suppressMessages(fast_d(
    .golden_tree(), x, test = TRUE, nsim = ncol(states$random),
    random_states = states$random, brownian_states = states$brownian,
    return_sim = TRUE, verbose = FALSE
  ))

  expect_s3_class(fit, "phylo.d")
  expect_equal(fit$DEstimate, -1.71830985915493, tolerance = 1e-10)
  expect_equal(fit$Pval1, 0, tolerance = 1e-12)
  expect_equal(fit$Pval0, 0.8, tolerance = 1e-12)
  expect_equal(fit$Parameters$Observed, 1, tolerance = 1e-12)
  expect_equal(fit$Parameters$MeanRandom, 1.64333333333333,
               tolerance = 1e-10)
  expect_equal(fit$Parameters$MeanBrownian, 1.40666666666667,
               tolerance = 1e-10)
  expect_equal(fit$Permutations$random,
               c(2.01666666666667, 2.01666666666667, 1, 1,
                 2.18333333333333), tolerance = 1e-10)
  expect_equal(fit$Permutations$brownian,
               c(1.6, 1.6, 1.41666666666667, 1.41666666666667, 1),
               tolerance = 1e-10)

  labels <- c(a = 10, b = 10, c = 20, d = 20)
  relabelled <- suppressMessages(fast_d(
    .golden_tree(), labels, test = TRUE, nsim = ncol(states$random),
    random_states = states$random, brownian_states = states$brownian,
    return_sim = FALSE, verbose = FALSE
  ))
  expect_equal(relabelled$DEstimate, fit$DEstimate, tolerance = 1e-12)
  expect_equal(relabelled$Pval1, fit$Pval1, tolerance = 1e-12)
  expect_equal(relabelled$Pval0, fit$Pval0, tolerance = 1e-12)
})

test_that("golden ER ancestral likelihoods are stable", {
  x <- factor(c(a = "A", b = "A", c = "B", d = "C"))
  fit <- suppressWarnings(fast_ace(x, .golden_tree(), model = "ER"))

  expect_s3_class(fit, "ace")
  expect_equal(fit$loglik, -3.285635060172, tolerance = 1e-8)
  expect_equal(fit$rates, 0.668752665860438, tolerance = 1e-8)
  expect_equal(fit$lik.anc, matrix(
    c(0.379570009827257, 0.322736327202790, 0.297693662969953,
      0.724615560559838, 0.138848938803113, 0.136535500637049,
      0.156291549659079, 0.506675342183192, 0.337033108157729),
    nrow = 3L, byrow = TRUE,
    dimnames = list(as.character(5:7), c("A", "B", "C"))
  ), tolerance = 1e-7)
})

test_that("golden entropy values and Delta chain-length contract are stable", {
  prob <- matrix(
    c(0.8, 0.1, 0.1,
      0.34, 0.33, 0.33,
      0.2, 0.5, 0.3),
    ncol = 3L, byrow = TRUE
  )
  expect_equal(as.numeric(fastphylosig:::delta_entropy_cpp(prob, 1L)),
               c(0.3, 0.99, 0.75), tolerance = 1e-12)
  expect_equal(as.numeric(fastphylosig:::delta_entropy_cpp(prob, 2L)),
               c(0.581671865717887, 0.999909274984070,
                 0.937230563216130), tolerance = 1e-12)
  expect_equal(as.numeric(fastphylosig:::delta_entropy_cpp(prob, 3L)),
               c(0.51, 0.9999, 0.93), tolerance = 1e-12)

  x <- factor(c(a = "A", b = "A", c = "B", d = "C"))
  set.seed(401)
  short <- suppressWarnings(fast_delta(
    .golden_tree(), x, test = FALSE, mcmc_sim = 400L, thin = 10L,
    burn = 100L, model = "ER", verbose = FALSE
  ))
  set.seed(402)
  long <- suppressWarnings(fast_delta(
    .golden_tree(), x, test = FALSE, mcmc_sim = 800L, thin = 10L,
    burn = 100L, model = "ER", verbose = FALSE
  ))
  .expect_delta_mcmc_contract(short)
  .expect_delta_mcmc_contract(long)
  # This is deliberately a broad Monte Carlo guard, not an equality claim.
  expect_gt(long$n_saved, short$n_saved)
  expect_lt(abs(short$delta - long$delta), 1)
  expect_equal(long$parameters$mcmc_sim, 800L)
  expect_equal(long$parameters$thin, 10L)
  expect_equal(long$parameters$burn, 100L)
})
