test_that("D reports separate null P values and exact Bernoulli MCSE", {
  set.seed(7301)
  tree <- ape::reorder.phylo(ape::rtree(10), "pruningwise")
  x <- setNames(rep(0:1, length.out = 10L), tree$tip.label)
  nsim <- 6L

  # Controlled states force the compatibility path and make the null draws
  # inspectable.  The tail definitions are deliberately different: random
  # uses < observed and Brownian uses > observed.
  random_states <- matrix(rep(x, nsim), nrow = 10L)
  random_states[, 1L] <- rev(x)
  brownian_states <- matrix(rep(x, nsim), nrow = 10L)
  brownian_states[, 1L] <- rev(x)
  names(x) <- tree$tip.label

  fit <- fastphylosig::fast_d(
    tree, x, nsim = nsim, test = TRUE, return_sim = TRUE,
    random_states = random_states, brownian_states = brownian_states,
    verbose = FALSE
  )
  observed <- fit$Parameters$Observed
  random <- fit$Permutations$random
  brownian <- fit$Permutations$brownian
  p_random <- sum(random < observed) / nsim
  p_brownian <- sum(brownian > observed) / nsim

  expect_equal(fit$Pval1, p_random, tolerance = 0)
  expect_equal(fit$Pval0, p_brownian, tolerance = 0)
  expect_equal(fit$P_random, p_random, tolerance = 0)
  expect_equal(fit$P_Brownian, p_brownian, tolerance = 0)
  expect_equal(fit$nsim_requested, nsim)
  expect_equal(fit$nsim_successful_random, nsim)
  expect_equal(fit$nsim_successful_brownian, nsim)
  expect_equal(fit$nsim_failed_random, 0L)
  expect_equal(fit$nsim_failed_brownian, 0L)
  expect_equal(
    fit$MCSE_P_random,
    sqrt(p_random * (1 - p_random) / nsim),
    tolerance = 0
  )
  expect_equal(
    fit$MCSE_P_Brownian,
    sqrt(p_brownian * (1 - p_brownian) / nsim),
    tolerance = 0
  )
})

test_that("D matrix output carries null-specific diagnostics", {
  set.seed(7302)
  tree <- ape::rtree(9)
  x <- setNames(rep(0:1, length.out = 9L), tree$tip.label)
  X <- cbind(trait_a = x, trait_b = 1 - x)
  nsim <- 5L
  random_states <- matrix(rep(x, nsim), nrow = 9L)
  brownian_states <- matrix(rep(1 - x, nsim), nrow = 9L)

  fit <- fastphylosig::fast_d(
    tree, X, nsim = nsim, test = TRUE, return_sim = TRUE,
    random_states = random_states, brownian_states = brownian_states,
    verbose = FALSE
  )
  expected <- c(
    "P_random", "P_Brownian", "MCSE_P_random", "MCSE_P_Brownian",
    "nsim_requested", "nsim_successful_random",
    "nsim_successful_brownian", "nsim_failed_random",
    "nsim_failed_brownian"
  )
  expect_true(all(expected %in% names(fit)))
  expect_equal(fit$nsim_requested, rep(nsim, 2L))
  expect_equal(fit$nsim_successful_random, rep(nsim, 2L))
  expect_equal(fit$nsim_successful_brownian, rep(nsim, 2L))
  expect_equal(fit$nsim_failed_random, rep(0L, 2L))
  expect_equal(fit$nsim_failed_brownian, rep(0L, 2L))
  expect_equal(fit$MCSE_P_random,
               sqrt(fit$P_random * (1 - fit$P_random) / nsim),
               tolerance = 0)
  expect_equal(fit$MCSE_P_Brownian,
               sqrt(fit$P_Brownian * (1 - fit$P_Brownian) / nsim),
               tolerance = 0)
})

test_that("boundary P values have zero MCSE and nonfinite draws are explicit", {
  zero <- fastphylosig:::.d_null_summary(
    c(1, 2, 3, 4), observed = 0, direction = "less"
  )
  one <- fastphylosig:::.d_null_summary(
    c(1, 2, 3, 4), observed = 10, direction = "less"
  )
  failed <- fastphylosig:::.d_null_summary(
    c(1, NA_real_, Inf, 2), observed = 1.5, direction = "less"
  )

  expect_equal(zero$p, 0)
  expect_equal(zero$mcse, 0)
  expect_equal(one$p, 1)
  expect_equal(one$mcse, 0)
  expect_equal(failed$n_successful, 2L)
  expect_equal(failed$n_failed, 2L)
  expect_equal(failed$p, 0.5)
  expect_equal(failed$mcse, sqrt(0.5 * 0.5 / 2), tolerance = 0)
  expect_match(failed$note, "excluded 2 non-finite")
})

test_that("test = FALSE hides P and MCSE but retains simulation counts", {
  set.seed(7303)
  tree <- ape::rtree(8)
  x <- setNames(rep(0:1, length.out = 8L), tree$tip.label)
  fit <- fastphylosig::fast_d(
    tree, x, nsim = 8L, test = FALSE, return_sim = FALSE,
    verbose = FALSE
  )
  expect_true(is.na(fit$P_random))
  expect_true(is.na(fit$P_Brownian))
  expect_true(is.na(fit$MCSE_P_random))
  expect_true(is.na(fit$MCSE_P_Brownian))
  expect_equal(fit$nsim_requested, 8L)
  expect_equal(fit$nsim_successful_random, 8L)
  expect_equal(fit$nsim_successful_brownian, 8L)
})
