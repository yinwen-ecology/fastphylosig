test_that("corrected Delta MCSE has honest finite-sampling boundaries", {
  for (B in c(1, 2, 5, 1000)) {
    expect_true(is.na(.delta_p_mcse(1 / (B + 1), B)))
    expect_true(is.na(.delta_p_mcse(1, B)))
    if (B > 1) for (b in unique(c(1, floor(B / 2), B - 1))) {
      expect_equal(.delta_p_mcse((b + 1) / (B + 1), B),
                   sqrt(B * (b / B) * (1 - b / B)) / (B + 1))
    }
  }
  for (p in c(NA_real_, Inf, -1, 2)) expect_true(is.na(.delta_p_mcse(p, 10)))
  for (B in c(NA_real_, Inf, 0, -1)) expect_true(is.na(.delta_p_mcse(.5, B)))
})

.delta_inference_fixture <- function() {
  tree <- ape::stree(8, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  x <- setNames(rep(c("A", "B"), 4), tree$tip.label)
  list(tree = tree, x = x)
}

.delta_inference_call <- function(f, x = f$x, ...) suppressWarnings(fast_delta(
  f$tree, x, test = TRUE, nsim = 4, mcmc_sim = 40, thin = 5, burn = 10,
  model = "ER", return_sim = TRUE, progress = FALSE, ...
))

test_that("Delta never tests a successful subset of failed permutations", {
  f <- .delta_inference_fixture()
  for (failed in list(2L, c(1L, 3L), 1:4)) local({
    bad <- failed
    testthat::local_mocked_bindings(.delta_permutation_worker = function(i, ...) {
      if (i %in% bad) NA_real_ else 0
    }, .package = "fastphylosig")
    set.seed(991)
    fit <- .delta_inference_call(f)
    expect_true(is.finite(fit$delta))
    expect_true(is.na(fit$P))
    expect_true(is.na(fit$MCSE_P))
    expect_identical(fit$status, "partial")
    expect_equal(fit$n_failed_sim, length(bad))
    expect_equal(fit$successful_simulations, 4 - length(bad))
    expect_equal(length(fit$sim.delta), 4)
    expect_match(fit$message, "P and MCSE_P are unavailable")
    X <- cbind(a = f$x, b = f$x)
    fit <- .delta_inference_call(f, X)
    expect_true(all(is.na(fit$P_fast)))
    expect_true(all(is.na(fit$P_MCSE)))
    expect_true(all(fit$status == "partial"))
    expect_equal(fit$n_failed_sim, rep(length(bad), 2))
  })
})

test_that("Delta rank P includes the observation at both exceedance boundaries", {
  f <- .delta_inference_fixture()
  for (value in c(0, Inf)) local({
    # Use a finite value above every possible observed value, not Inf failure.
    v <- if (is.infinite(value)) .Machine$double.xmax else value
    testthat::local_mocked_bindings(.delta_permutation_worker = function(...) v, .package = "fastphylosig")
    fit <- .delta_inference_call(f)
    expect_equal(fit$P, if (v == 0) 1 / 5 else 1)
    expect_true(is.na(fit$P_MCSE))
  })
})

test_that("random and supplied permutations use the documented rank formula", {
  f <- .delta_inference_fixture()
  set.seed(992)
  perms <- t(replicate(4, sample.int(8)))
  for (cores in c(1L, 2L)) for (custom in c(FALSE, TRUE)) {
    set.seed(993)
    fit <- .delta_inference_call(f, ncores = cores,
      permutations = if (custom) perms else NULL)
    expect_equal(fit$successful_simulations, 4)
    b <- sum(fit$sim.delta >= fit$delta)
    expect_equal(fit$P, (b + 1) / 5)
    expect_equal(fit$P_MCSE, if (b %in% c(0, 4)) NA_real_ else
      sqrt(4 * (b / 4) * (1 - b / 4)) / 5)
  }
})

test_that("retaining Delta chains does not change statistic or RNG", {
  f <- .delta_inference_fixture()
  f$x <- factor(f$x)
  # Match the public Delta path, which does not request ACE rate standard errors.
  workspace <- .fast_ace_workspace(.fast_ace_tree_workspace(f$tree),
                                   model = "ER", nl = nlevels(factor(f$x)))
  args <- list(tree = f$tree, y = f$x, lambda0 = .1, proposal_sd = .5,
    mcmc_sim = 40L, thin = 5L, burn = 10L, entropy_code = 1L,
    model = "ER", ace_engine = "fast", ace_workspace = workspace)
  set.seed(994)
  a <- do.call(.delta_one, c(args, list(keep_chains = TRUE)))
  rng <- .Random.seed
  set.seed(994)
  b <- do.call(.delta_one, c(args, list(keep_chains = FALSE)))
  expect_identical(a$delta, b$delta)
  expect_identical(.Random.seed, rng)
  expect_false(a$calculation_diagnostics$other_warning)
  expect_false(b$calculation_diagnostics$other_warning)
})
test_that("Delta failed-permutation warning is emitted once per batch", {
  f <- .delta_inference_fixture()
  testthat::local_mocked_bindings(.delta_permutation_worker = function(i, ...) {
    if (i == 2L) NA_real_ else 0
  }, .package = "fastphylosig")
  messages <- character()
  fit <- withCallingHandlers(fast_delta(f$tree, cbind(a=f$x,b=f$x),
    test=TRUE, nsim=4, mcmc_sim=40, thin=5, burn=10, model="ER",
    progress=FALSE), warning=function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  expect_equal(sum(grepl("Delta permutation failures for", messages)), 1L)
})
test_that("Delta exceedance counting includes ties and mixes both tails", {
  f <- .delta_inference_fixture()
  set.seed(995)
  observed <- suppressWarnings(fast_delta(f$tree, f$x, mcmc_sim=40,
    thin=5, burn=10, model="ER", progress=FALSE))$delta
  testthat::local_mocked_bindings(.delta_permutation_worker = function(i, ...) {
    c(0, observed, observed, .Machine$double.xmax)[i]
  }, .package="fastphylosig")
  set.seed(995)
  fit <- .delta_inference_call(f)
  expect_identical(fit$delta, observed)
  expect_equal(fit$P, 4 / 5)
  expect_equal(fit$P_MCSE, sqrt(4 * .75 * .25) / 5)
})
test_that("parallel Delta failures also suppress selected-subset inference", {
  f <- .delta_inference_fixture()
  testthat::local_mocked_bindings(.delta_permutation_worker = function(i, ...) {
    if (i == 2L) NA_real_ else 0
  }, .package="fastphylosig")
  fit <- .delta_inference_call(f, ncores=2L)
  expect_equal(fit$n_failed_sim, 1L)
  expect_equal(fit$successful_simulations, 3L)
  expect_identical(fit$status, "partial")
  expect_true(is.na(fit$P))
  expect_true(is.na(fit$P_MCSE))
})