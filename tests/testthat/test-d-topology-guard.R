test_that("D keeps ordinary binary controlled-state results finite", {
  skip_if_not_installed("ape")
  tree <- ape::rtree(6L)
  x <- stats::setNames(c(0, 0, 1, 1, 0, 1), tree$tip.label)
  nsim <- 3L
  random_states <- matrix(
    c(0, 1, 0, 1, 0, 1,
      1, 0, 1, 0, 1, 0,
      0, 0, 1, 1, 1, 0),
    nrow = 6L, ncol = nsim
  )
  brownian_states <- matrix(
    c(0, 0, 1, 0, 1, 1,
      1, 1, 0, 1, 0, 0,
      0, 1, 0, 1, 1, 0),
    nrow = 6L, ncol = nsim
  )

  fit <- fast_d(
    tree, x, nsim = nsim, random_states = random_states,
    brownian_states = brownian_states, return_sim = TRUE,
    verbose = FALSE, progress = FALSE
  )
  expect_true(all(is.finite(c(
    fit$DEstimate, fit$Pval1, fit$Pval0,
    fit$Parameters$Observed, fit$Parameters$MeanRandom,
    fit$Parameters$MeanBrownian
  ))))
})

test_that("D keeps rooted polytomy and non-ultrametric controlled results finite", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(
    text = "((a:1,b:2,c:1):3,d:2);"
  )
  x <- c(a = 0, b = 0, c = 1, d = 1)
  random_states <- matrix(
    c(0, 1, 0, 1,
      1, 0, 1, 0,
      0, 0, 1, 1),
    nrow = 4L, ncol = 3L
  )
  brownian_states <- matrix(
    c(0, 0, 1, 0,
      1, 1, 0, 1,
      0, 1, 0, 0),
    nrow = 4L, ncol = 3L
  )

  fit <- fast_d(
    tree, x, nsim = 3L, random_states = random_states,
    brownian_states = brownian_states, return_sim = TRUE,
    verbose = FALSE, progress = FALSE
  )
  expect_true(all(is.finite(c(
    fit$DEstimate, fit$Pval1, fit$Pval0,
    fit$Parameters$Observed, fit$Parameters$MeanRandom,
    fit$Parameters$MeanBrownian
  ))))
})

test_that("D rejects an internal single-child node at the R boundary", {
  skip_if_not_installed("ape")
  tree <- list(
    edge = matrix(c(3L, 4L, 4L, 1L, 3L, 2L), ncol = 2L, byrow = TRUE),
    edge.length = c(1, 1, 1),
    Nnode = 2L,
    tip.label = c("a", "b")
  )
  class(tree) <- "phylo"
  x <- c(a = 0, b = 1)

  expect_error(
    fast_d(tree, x, nsim = 2L, return_sim = FALSE,
           verbose = FALSE, progress = FALSE),
    "single-child"
  )
})

test_that("D canonicalizes a representation-only noncanonical structural root", {
  skip_if_not_installed("ape")
  tree <- list(
    edge = matrix(c(
      6L, 5L, 6L, 7L, 5L, 1L,
      5L, 2L, 7L, 3L, 7L, 4L
    ), ncol = 2L, byrow = TRUE),
    edge.length = rep(1, 6L),
    Nnode = 3L,
    tip.label = c("a", "b", "c", "d")
  )
  class(tree) <- "phylo"
  x <- c(a = 0, b = 0, c = 1, d = 1)

  random_states <- matrix(c(
    0, 1, 0, 1,
    1, 0, 1, 0
  ), nrow = 4L, ncol = 2L)
  brownian_states <- matrix(c(
    0, 0, 1, 1,
    1, 1, 0, 0
  ), nrow = 4L, ncol = 2L)
  fit <- fast_d(
    tree, x, nsim = 2L, return_sim = FALSE,
    random_states = random_states, brownian_states = brownian_states,
    verbose = FALSE, progress = FALSE
  )
  expect_s3_class(fit, "phylo.d")
  expect_true(isTRUE(attr(fit, "analysis_metadata")$tree_auto_normalized))
})
