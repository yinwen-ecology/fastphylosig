test_that("ARD fixed-rate fits do not depend on state names", {
  original_nlminb <- stats::nlminb
  # Resolve ape's lazy optimizer import before temporarily mocking stats.
  # Otherwise a later ace() call can retain the mock after stats is restored.
  original_ape_nlminb <- get("nlminb", asNamespace("ape"), inherits = TRUE)
  original_recheck <- .fast_ace_recheck_rate_scale
  run_fixture <- function() {
    testthat::local_mocked_bindings(nlminb = function(start, objective, ...) {
      par <- seq_along(start) * .01
      list(par = par, objective = objective(par), convergence = 0L,
        message = "fixed-rate mapping test", iterations = 0L,
        evaluations = c("function" = 1L))
    }, .package = "stats")
    testthat::local_mocked_bindings(.fast_ace_recheck_rate_scale =
      function(dev, first, start, workspace) first)
  for (k in 2:5) {
    set.seed(808)
    tree <- ape::rtree(48)
    # Relabelling needs every state present; stochastic simulation can omit one.
    x <- stats::setNames(factor(rep(letters[seq_len(k)], length.out = 48),
      levels = letters[seq_len(k)]), tree$tip.label)
    expect_identical(sort(unique(as.integer(x))), seq_len(k))
    labels <- levels(x)
    baseline <- suppressWarnings(fast_ace(x, tree, model = "ARD", progress = FALSE))
    for (permutation in list(rev(seq_len(k)), c(seq.int(2L, k), 1L))) {
      y <- stats::setNames(labels[permutation[as.integer(x)]], names(x))
      set.seed(622)
      before <- .Random.seed
      actual <- suppressWarnings(fast_ace(y, tree, model = "ARD", progress = FALSE))
      expect_identical(.Random.seed, before)
      expect_identical(actual$loglik, baseline$loglik)
      expect_identical(unname(actual$lik.anc[, permutation]), unname(baseline$lik.anc))
      to_matrix <- function(z, field) {
        out <- matrix(0, k, k)
        ind <- z$index.matrix
        out[!is.na(ind)] <- z[[field]][ind[!is.na(ind)]]
        out
      }
      expect_identical(to_matrix(actual, "rates")[permutation, permutation],
                       to_matrix(baseline, "rates"))
      expect_identical(to_matrix(actual, "se")[permutation, permutation],
                       to_matrix(baseline, "se"))
      expect_identical(actual$index.matrix, baseline$index.matrix)
    }
  }
  }
  run_fixture()
  expect_identical(stats::nlminb, original_nlminb)
  expect_identical(get("nlminb", asNamespace("ape"), inherits = TRUE),
                   original_ape_nlminb)
  expect_identical(.fast_ace_recheck_rate_scale, original_recheck)
})

test_that("ARD external initial parameters and diagnostic orders are preserved", {
  set.seed(808)
  tree <- ape::rtree(48)
  x <- stats::setNames(factor(rep(letters[1:4], length.out = 48),
    levels = letters[1:4]), tree$tip.label)
  w <- .fast_ace_workspace(.fast_ace_tree_workspace(tree), "ARD", 4L)
  states <- as.integer(x)
  labels <- levels(x)
  # Force nonidentity internal state ordering without changing memberships.
  order0 <- .fast_ace_state_order(states, w)
  permutation <- rev(order0)
  states <- match(states, permutation)
  labels <- labels[permutation]
  internal_order <- .fast_ace_state_order(states, w)
  expect_false(identical(internal_order, 1:4))
  ri <- w$rate$rate_index
  map <- ri[internal_order, internal_order][ri > 0L]
  inverse <- order(map)
  ip <- seq(.01, .12, .01)
  original_impl <- .fast_ace_fit_workspace_impl
  reference <- suppressWarnings(original_impl(
    match(states, internal_order), labels[internal_order], w, ip[map],
    estimate_se = FALSE, fit_strategy = "multistart"))
  run_fixture <- function() {
    seen <- 0L
    testthat::local_mocked_bindings(.fast_ace_fit_workspace_impl =
      function(tip_state, lvls, workspace, ip, ...) {
        seen <<- seen + 1L
        expect_identical(tip_state, match(states, internal_order))
        expect_identical(lvls, labels[internal_order])
        expect_identical(workspace$edge, w$edge)
        expect_identical(workspace$edge_length, w$edge_length)
        expect_identical(workspace$rate$rate_index, w$rate$rate_index)
        expect_identical(ip, ip_external[map])
        reference
      })
    ip_external <- ip
    result <- .fast_ace_fit_workspace(states, labels, w, ip,
      estimate_se = FALSE, fit_strategy = "multistart")
    expect_identical(seen, 1L)
    result
  }
  actual <- run_fixture()
  expect_identical(.fast_ace_fit_workspace_impl, original_impl)
  expect_identical(actual$loglik, reference$loglik)
  expect_identical(actual$rates, reference$rates[inverse])
  expect_identical(actual$lik.anc, reference$lik.anc[, order(internal_order)])
  expect_identical(actual$diagnostics$starts[[1L]], ip)
  expect_identical(actual$diagnostics, .fast_ace_remap_parameters(reference$diagnostics, inverse))
  # Test nested optimizer diagnostic remapping, including warning conditions.
  sample <- list(optimizer_retry = list(original = list(par = 1:12)),
                 starts = list(1:12), attempts = data.frame(start = "baseline"))
  mapped <- .fast_ace_remap_parameters(sample, inverse)
  expect_identical(mapped$optimizer_retry$original$par, (1:12)[inverse])
  expect_identical(mapped$starts[[1L]], (1:12)[inverse])
  expect_identical(mapped$attempts, sample$attempts)
})

test_that("Delta fixed-rate ARD preserves sampling under relabelling", {
  original_nlminb <- stats::nlminb
  # Resolve ape's lazy optimizer import before temporarily mocking stats.
  # Otherwise a later ace() call can retain the mock after stats is restored.
  original_ape_nlminb <- get("nlminb", asNamespace("ape"), inherits = TRUE)
  original_recheck <- .fast_ace_recheck_rate_scale
  run_fixture <- function() {
    testthat::local_mocked_bindings(nlminb = function(start, objective, ...) {
      par <- seq_along(start) * .01
      list(par = par, objective = objective(par), convergence = 0L,
        message = "fixed-rate mapping test", iterations = 0L,
        evaluations = c("function" = 1L))
    }, .package = "stats")
    testthat::local_mocked_bindings(.fast_ace_recheck_rate_scale =
      function(dev, first, start, workspace) first)

  set.seed(808)
  tree <- ape::rtree(48)
  x <- stats::setNames(factor(rep(letters[1:4], length.out = 48),
    levels = letters[1:4]), tree$tip.label)
  y <- stats::setNames(rev(levels(x))[as.integer(x)], names(x))
  fit <- function(z) {
    set.seed(123)
    suppressWarnings(fast_delta(tree, z, model = "ARD", mcmc_sim = 1000,
      thin = 10, burn = 100, progress = FALSE))
  }
  a <- fit(x)
  seed_a <- .Random.seed
  b <- fit(y)
  expect_identical(.Random.seed, seed_a)
  expect_identical(a$delta, b$delta)
  expect_identical(a$status, b$status)
  expect_true(is.finite(a$delta))
  }
  run_fixture()
  expect_identical(stats::nlminb, original_nlminb)
  expect_identical(get("nlminb", asNamespace("ape"), inherits = TRUE),
                   original_ape_nlminb)
  expect_identical(.fast_ace_recheck_rate_scale, original_recheck)
})

test_that("ARD optimizer warning diagnostics use external parameter order", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  workspace <- .fast_ace_workspace(.fast_ace_tree_workspace(tree), "ARD", 4L)
  state <- 4:1
  order <- .fast_ace_state_order(state, workspace)
  ri <- workspace$rate$rate_index
  map <- ri[order, order][ri > 0L]
  inverse <- order(map)
  fake <- function(...) {
    diagnostics <- list(optimizer_retry = list(original = list(par = 1:12)))
    warning(structure(list(message = "test recovery", call = NULL,
      diagnostics = diagnostics), class = c("fastphylosig_optimizer_retry_warning",
                                            "warning", "condition")))
    list(rates = 1:12, se = 1:12, diagnostics = diagnostics)
  }
  testthat::local_mocked_bindings(.fast_ace_fit_workspace_impl = fake)
  caught <- NULL
  result <- withCallingHandlers(.fast_ace_fit_workspace(state, letters[1:4],
    workspace, .1, CI = FALSE, estimate_se = FALSE), warning = function(w) {
      caught <<- w
      invokeRestart("muffleWarning")
    })
  expect_s3_class(caught, "fastphylosig_optimizer_retry_warning")
  expect_identical(caught$diagnostics$optimizer_retry$original$par, (1:12)[inverse])
  expect_identical(caught$diagnostics, result$diagnostics)
})
