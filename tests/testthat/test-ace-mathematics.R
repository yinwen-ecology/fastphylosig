test_that("binary ARD probabilities agree with closed-form latent-state sums", {
  trees <- c("((a:0.2,b:0.7):0.5,(c:0.3,d:0.9):0.4);",
             "(((a:0.2,b:0.7):0.5,c:0.3):0.4,d:0.9);")
  index <- matrix(c(0L, 1L, 2L, 0L), 2L)
  for (text in trees) for (missing_tip in c(FALSE, TRUE)) {
    tree <- ape::read.tree(text = text)
    states <- c(1L, 2L, 1L, 1L)
    if (missing_tip) states[4L] <- NA_integer_
    Q <- matrix(c(-.15, .15, 1.4, -1.4), 2L, byrow = TRUE)
    truth <- .ace_math_enumerate(tree, states, Q)
    w <- .fast_ace_tree_workspace(tree)
    set.seed(8271); seed <- .Random.seed
    for (reuse in c(FALSE, TRUE)) {
      got <- fast_ace_discrete_liks_cpp(w$edge, w$edge_length, states,
        index, c(1.4, .15), FALSE, reuse)
      expect_equal(unname(got$lik.anc), unname(truth$probability), tolerance = 1e-12)
      expect_equal(got$loglik - log(2), truth$loglik, tolerance = 1e-12)
    }
    ptr <- fast_ace_fit_prepare_cpp(w$edge, w$edge_length, states, index, 1L)
    got <- fast_ace_fit_liks_cpp(ptr, c(1.4, .15), FALSE)
    fast_ace_fit_release_cpp(ptr)
    expect_equal(unname(got$lik.anc), unname(truth$probability), tolerance = 1e-12)
    expect_identical(.Random.seed, seed)
  }
})

test_that("three-state ARD and reducible boundaries agree with exhaustive sums", {
  skip_if_not_installed("expm")
  index <- matrix(0L, 3L, 3L)
  index[row(index) != col(index)] <- 1:6
  generators <- list(
    # Includes a complex conjugate eigensystem; do not use ape::matexpo as truth.
    matrix(c(0, .7, .4, .2, 0, 1.7, 1.1, .1, 0), 3L, byrow = TRUE),
    # Reversible with a nonuniform stationary distribution, hence asymmetric P.
    matrix(c(0, .3, .28, .5, 0, .16, .7, .24, 0), 3L, byrow = TRUE),
    # Two closed absorbing classes: impossible parent-state messages are 0/0.
    matrix(c(0, 0, 0, 0, 0, 0, .4, .6, 0), 3L, byrow = TRUE))
  for (Q in generators) {
    diag(Q) <- -rowSums(Q)
    rates <- Q[row(Q) != col(Q)]
    for (text in c("((a:0.2,b:0.7):0.5,(c:0.3,d:0.9):0.4);",
                   "(((a:0.2,b:0.7):0.5,c:0.3):0.4,d:0.9);")) {
      tree <- ape::read.tree(text = text); states <- c(1L, 2L, 3L, 1L)
      truth <- .ace_math_enumerate(tree, states, Q)
      w <- .fast_ace_tree_workspace(tree)
      for (reuse in c(FALSE, TRUE)) {
        got <- fast_ace_discrete_liks_cpp(w$edge, w$edge_length, states,
          index, rates, FALSE, reuse)
        expect_equal(unname(got$lik.anc), unname(truth$probability), tolerance = 1e-11)
        expect_equal(got$loglik - log(3), truth$loglik, tolerance = 1e-11)
      }
      ptr <- fast_ace_fit_prepare_cpp(w$edge, w$edge_length, states, index, 1L)
      got <- fast_ace_fit_liks_cpp(ptr, rates, FALSE)
      fast_ace_fit_release_cpp(ptr)
      expect_equal(unname(got$lik.anc), unname(truth$probability), tolerance = 1e-11)
    }
  }
})

test_that("public ARD probabilities satisfy the full-tree probability model", {
  skip_if_not_installed("expm")
  tree <- ape::read.tree(text = "((a:0.2,b:0.7):0.5,(c:0.3,d:0.9):0.4);")
  x <- c(a = "a", b = "b", c = "c", d = "a")
  set.seed(8272); seed <- .Random.seed
  # This small tree can give boundary rates and undefined numerical rate SEs.
  # Here only the full-tree probabilities are the target; interior SE curvature
  # is independently checked below on larger, identifiable ER fits.
  fit <- suppressWarnings(fast_ace(x, tree, model = "ARD", progress = FALSE))
  expect_identical(fit$convergence, 0L)
  Q <- .ace_math_q(fit$index.matrix, fit$rates)
  truth <- .ace_math_enumerate(tree, as.integer(factor(x)), Q)
  expect_equal(unname(fit$lik.anc), unname(truth$probability[rownames(fit$lik.anc), ]),
               tolerance = 1e-10)
  expect_equal(fit$loglik - log(3), truth$loglik, tolerance = 1e-10)
  expect_identical(.Random.seed, seed)
})

test_that("Delta receives mathematically correct ARD probabilities", {
  skip_if_not_installed("expm")
  tree <- ape::read.tree(text = "((a:0.2,b:0.7):0.5,(c:0.3,d:0.9):0.4);")
  y <- factor(c(a = "a", b = "b", c = "c", d = "a"))
  w <- .fast_ace_workspace(.fast_ace_tree_workspace(tree), "ARD", 3L)
  fit <- .fast_ace_fit_workspace(as.integer(y), levels(y), w, .1,
    estimate_se = FALSE, reuse_buffers = TRUE, recover_ard = TRUE)
  truth <- .ace_math_enumerate(tree, as.integer(y), .ace_math_q(fit$index.matrix, fit$rates))
  sampler <- delta_mcmc_cpp
  seen <- NULL
  local_mocked_bindings(delta_mcmc_cpp = function(probabilities, ...) {
    seen <<- probabilities
    sampler(probabilities, ...)
  }, .package = "fastphylosig")
  .delta_one(tree, y, .1, .5, 1000L, 10L, 100L, 1L, "ARD", "fast",
             ace_workspace = w)
  expect_equal(unname(seen), unname(truth$probability), tolerance = 1e-10)
})

test_that("ACE standard errors use likelihood information, not deviance curvature", {
  for (s in c(101L, 102L, 103L)) {
    set.seed(s); tree <- ape::rtree(128L)
    x <- ape::rTraitDisc(tree, model = "ER", k = 2L, rate = .35,
                        states = c("A", "B"))
    fit <- fast_ace(x, tree, model = "ER", progress = FALSE)
    expect_identical(fit$convergence, 0L)
    expect_gt(fit$rates, 0)
    states <- as.integer(factor(x[tree$tip.label]))
    nll <- function(q) -.ace_math_loglik(tree, states,
      matrix(c(-q, q, q, -q), 2L, byrow = TRUE))
    expect_equal(fit$loglik - log(2), -nll(fit$rates), tolerance = 1e-11)
    for (relative_step in c(1e-3, 3e-4)) {
      h <- fit$rates * relative_step
      information <- (nll(fit$rates+h) - 2*nll(fit$rates) + nll(fit$rates-h)) / h^2
      expect_equal(as.numeric(fit$se), sqrt(1/information), tolerance = .003)
    }
  }
})
