test_that("dimensionless ACE fits recover the analytic four-tip optimum", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- c(a = "A", b = "A", c = "B", d = "B")
  for (scale in c(1e-40, 1e-20, 1e-9, 1, 1e9, 1e20, 1e40)) {
    scaled <- tree
    scaled$edge.length <- tree$edge.length * scale
    fit <- suppressWarnings(fast_ace(x, scaled, progress = FALSE))
    expect_equal(fit$loglik, -log(6), tolerance = 1e-7)
    expect_equal(fit$rates * scale, log(3) / 4, tolerance = 1e-5)
  }
})

test_that("asymmetric fits are optimized before scale candidates are rejected", {
  for (seed in c(8, 11, 101)) {
    set.seed(seed)
    tree <- ape::rtree(48)
    x <- ape::rTraitDisc(tree, model = "ER", k = 2, rate = .2,
                        states = LETTERS[1:2])
    for (model in c("ER", "ARD")) {
      base <- suppressWarnings(fast_ace(x, tree, model = model, progress = FALSE))
      for (scale in c(1e-9, 1e3, 1e9)) {
        scaled <- tree
        scaled$edge.length <- tree$edge.length * scale
        fit <- suppressWarnings(fast_ace(x, scaled, model = model, progress = FALSE))
        expect_equal(fit$loglik, base$loglik, tolerance = 1e-7)
        expect_equal(fit$rates * scale, base$rates, tolerance = 1e-5)
        expect_equal(fit$lik.anc, base$lik.anc, tolerance = 1e-5)
      }
    }
  }
})

test_that("optimized scale candidates recover an invalid baseline", {
  first <- list(par = .1, objective = 1e50, convergence = 0L)
  fit <- .fast_ace_recheck_rate_scale(function(p) (p - .3)^2,
    first, .1, list(edge_length = c(1, 1)))
  expect_equal(fit$par, .3, tolerance = 1e-6)
  expect_gte(nrow(fit$diagnostics$attempts), 4L)
  expect_equal(fit$diagnostics$attempts$objective[
    fit$diagnostics$selected_attempt], fit$objective)
})
test_that("multistart uses dimensionless optimization and physical rate units", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- factor(c("A", "A", "B", "B"))
  run <- function(scale, model) {
    scaled <- tree
    scaled$edge.length <- tree$edge.length * scale
    workspace <- .fast_ace_workspace(.fast_ace_tree_workspace(scaled), model, 2L)
    suppressWarnings(.fast_ace_fit_workspace(as.integer(x), levels(x),
      workspace, ip = .1 / scale, CI = TRUE, estimate_se = TRUE,
      fit_strategy = "multistart"))
  }
  for (model in c("ER", "ARD")) {
    baseline <- run(1, model)
    for (scale in c(1e-9, 1, 1e9)) {
      fit <- run(scale, model)
      expect_equal(fit$loglik, baseline$loglik, tolerance = 1e-7)
      expect_equal(fit$rates * scale, baseline$rates, tolerance = 1e-5)
      expect_true(all(is.finite(fit$se)))
      expect_equal(fit$se * scale, baseline$se, tolerance = .003)
      expect_equal(lapply(fit$diagnostics$starts, function(p) p * scale),
                   baseline$diagnostics$starts, tolerance = 1e-7)
      expect_equal(fit$diagnostics$starts[[1L]],
                   rep(.1 / scale, length(fit$rates)))
    }
  }
})

test_that("ARD scale recovery retains likelihood and ancestral probabilities", {
  set.seed(7)
  tree <- ape::rtree(48)
  x <- ape::rTraitDisc(tree, model = "ER", k = 3, rate = .2,
                      states = LETTERS[1:3])
  baseline <- suppressWarnings(fast_ace(x, tree, model = "ARD", progress = FALSE))
  for (scale in c(1e-9, 1000)) {
    scaled <- tree
    scaled$edge.length <- tree$edge.length * scale
    fit <- suppressWarnings(fast_ace(x, scaled, model = "ARD", progress = FALSE))
    expect_gte(fit$loglik, -19.15518)
    expect_equal(fit$loglik, baseline$loglik, tolerance = 1e-7)
    expect_equal(fit$lik.anc, baseline$lik.anc, tolerance = 1e-4)
    attempts <- fit$diagnostics$attempts
    # Different BLAS/optimizer paths need not exhaust the initial budget.
    # Diagnostics exist only when a scale candidate replaces the baseline.
    if (!is.null(attempts)) {
      expect_equal(attempts$objective[fit$diagnostics$selected_attempt],
                   -2 * fit$loglik)
    }
    expect_identical(fit$convergence, 0L)
  }
})

test_that("scalar ER refinement agrees with an independent analytic score", {
  tree <- ape::reorder.phylo(ape::read.tree(
    text = "((a:0.5,b:0.7):0.8,(c:0.4,d:0.6):0.9);"), "postorder")
  score <- function(q) {
    likelihood <- matrix(1, 7, 3)
    derivative <- matrix(0, 7, 3)
    likelihood[1:4, ] <- diag(3)[c(1, 1, 2, 3), ]
    for (i in seq_len(nrow(tree$edge))) {
      parent <- tree$edge[i, 1]; child <- tree$edge[i, 2]
      t <- tree$edge.length[i]; v <- exp(-3 * q * t)
      P <- matrix((1 - v) / 3, 3, 3) + diag(v, 3)
      dP <- matrix(t * v, 3, 3) - diag(3 * t * v, 3)
      a <- drop(P %*% likelihood[child, ])
      b <- drop(dP %*% likelihood[child, ] + P %*% derivative[child, ])
      derivative[parent, ] <- derivative[parent, ] * a + likelihood[parent, ] * b
      likelihood[parent, ] <- likelihood[parent, ] * a
    }
    sum(derivative[5, ]) / sum(likelihood[5, ])
  }
  optimum <- uniroot(score, c(.5, .8), tol = 1e-14)$root
  fit <- suppressWarnings(fast_ace(c(a = "A", b = "A", c = "B", d = "C"),
                                  tree, model = "ER", progress = FALSE))
  # State the rate accuracy in absolute units; the analytic score below
  # independently checks stationarity, rather than optimizer trajectory.
  expect_lt(abs(fit$rates - optimum), 1e-7)
  expect_lt(abs(score(fit$rates)), 1e-7)
  expect_equal(fit$loglik, -3.285635060172, tolerance = 1e-10)
})

test_that("unfinished scale candidates receive the extended reference budget", {
  calls <- list()
  testthat::local_mocked_bindings(nlminb = function(start, objective,
      lower, upper, control = list()) {
    extended <- identical(control, list(iter.max = 1500L, eval.max = 4000L))
    calls[[length(calls) + 1L]] <<- list(start = start, control = control)
    par <- rep(if (extended) .3 else .2, length(start))
    list(par = par, objective = objective(par),
      convergence = if (extended) 0L else 1L,
      iterations = if (extended) 200L else 150L,
      evaluations = c("function" = 250L),
      message = if (extended) "converged" else "iteration limit")
  }, .package = "stats")
  objective <- function(p) sum((p - .3)^2)
  first <- list(par = c(.1, .1), objective = objective(c(.1, .1)),
    convergence = 0L)
  fit <- .fast_ace_recheck_rate_scale(objective, first, first$par,
    list(edge_length = c(1, 1)))
  expect_length(calls, 6L)
  expect_true(all(vapply(calls[1:3], function(z) length(z$control) == 0L,
    logical(1))))
  for (i in seq_len(3L)) {
    expect_identical(calls[[i + 3L]]$start, calls[[i]]$start)
    expect_identical(calls[[i + 3L]]$control,
      list(iter.max = 1500L, eval.max = 4000L))
  }
  attempts <- fit$diagnostics$attempts
  expect_identical(attempts$convergence, c(0L, 1L, 1L, 1L, 0L, 0L, 0L))
  expect_equal(sum(grepl("_extended$", attempts$start)), 3L)
  expect_equal(attempts$objective[fit$diagnostics$selected_attempt], fit$objective)
  expect_equal(fit$par, c(.3, .3))
  expect_identical(fit$convergence, 0L)
})