.ace_scale_fixture <- function() {
  set.seed(101)
  tree <- ape::rtree(128L)
  x <- ape::rTraitDisc(tree, model = "ER", k = 2L, rate = .35,
                      states = c("A", "B"))
  list(tree = tree, x = x)
}

test_that("ACE rate and ancestral estimates survive branch-unit changes", {
  z <- .ace_scale_fixture()
  for (model in c("ER", "ARD")) {
    baseline <- suppressWarnings(fast_ace(z$x, z$tree, model = model,
      progress = FALSE))
    for (scale in c(1e-3, 1e3, 1e4)) {
      tree <- z$tree
      tree$edge.length <- tree$edge.length * scale
      set.seed(8361)
      seed <- .Random.seed
      actual <- suppressWarnings(fast_ace(z$x, tree, model = model,
        progress = FALSE))
      expect_identical(.Random.seed, seed)
      expect_identical(actual$convergence, 0L)
      expect_equal(actual$loglik, baseline$loglik, tolerance = 1e-7)
      expect_equal(actual$rates * scale, baseline$rates, tolerance = 1e-5)
      expect_equal(actual$lik.anc, baseline$lik.anc, tolerance = 1e-5)
      expect_true(all(is.finite(actual$se)))
      expect_equal(actual$se * scale, baseline$se, tolerance = .003)
      if (!is.null(actual$diagnostics$rate_scale_recovery)) {
        expect_gt(actual$diagnostics$deviance_improvement, 0)
        expect_true(actual$diagnostics$optimizer_retry$recovered)
        expect_gte(nrow(actual$diagnostics$attempts), 2L)
        expect_equal(actual$diagnostics$attempts$objective[actual$diagnostics$selected_attempt], -2 * actual$loglik)
      }
      # Explicit ip is still a rate in the supplied tree's physical units.
      explicit <- suppressWarnings(fast_ace(z$x, tree, model = model,
        ip = .1 / scale, progress = FALSE))
      expect_null(explicit$diagnostics$rate_scale_recovery)
      expect_equal(explicit$rates * scale, baseline$rates, tolerance = 1e-5)
      expect_equal(explicit$loglik, baseline$loglik, tolerance = 1e-7)
      expect_equal(explicit$se * scale, baseline$se, tolerance = .003)
    }
  }
})

test_that("rate-scale recovery retains the baseline when optimized candidates fail or do not improve", {
  fit <- list(par = .2, objective = 12, convergence = 0L)
  set.seed(8362)
  seed <- .Random.seed
  expect_identical(.fast_ace_recheck_rate_scale(
    function(p) stop("unexpected objective evaluation"), fit, .1,
    list(edge_length = c(.2, .8))), fit)
  expect_identical(.fast_ace_recheck_rate_scale(
    function(p) 14, fit, .1, list(edge_length = c(200, 800))), fit)
  expect_identical(.fast_ace_recheck_rate_scale(
    function(p) stop("failed probe"), fit, .1,
    list(edge_length = c(200, 800))), fit)
  expect_identical(.Random.seed, seed)
})

test_that("Delta uses the recovered ACE fit and exposes optimizer diagnostics", {
  z <- .ace_scale_fixture()
  large <- z$tree
  large$edge.length <- large$edge.length * 1000
  run <- function(tree) {
    set.seed(123)
    result <- suppressWarnings(fast_delta(tree, z$x, model = "ER",
      mcmc_sim = 1000L, thin = 10L, burn = 100L, progress = FALSE))
    list(result = result, seed = .Random.seed)
  }
  original <- run(z$tree)
  actual <- run(large)
  expect_equal(actual$result$delta, original$result$delta, tolerance = 1e-7)
  expect_identical(actual$seed, original$seed)
  expect_identical(actual$result$status, "ok")
  if (!is.null(actual$result$optimizer_retry)) {
    expect_true(actual$result$optimizer_retry)
    expect_true(actual$result$optimizer_recovered)
    expect_false(isTRUE(actual$result$optimizer_warning))
  }
})

test_that("Delta crown and overall fits share ACE branch-scale recovery", {
  z <- .ace_scale_fixture()
  large <- z$tree
  large$edge.length <- large$edge.length * 1000
  run <- function(tree) {
    set.seed(123)
    suppressWarnings(fast_hetero(tree, z$x, method = "Delta",
      nodes = length(tree$tip.label) + 1L, min_tips = 20L, overall = TRUE,
      delta_control = list(model = "ER", mcmc_sim = 1000L,
                           thin = 10L, burn = 100L), progress = FALSE))
  }
  original <- run(z$tree)
  actual <- run(large)
  expect_equal(actual$clades$estimate, original$clades$estimate, tolerance = 1e-7)
  expect_equal(actual$overall$estimate, original$overall$estimate, tolerance = 1e-7)
  expect_true(all(actual$clades$status == "ok"))
})
