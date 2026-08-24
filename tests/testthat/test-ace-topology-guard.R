test_that("fast_ace rejects mixed unary and polytomous topologies in R", {
  skip_if_not_installed("ape")

  # The edge count and Nnode value look binary for four tips (six edges and
  # three internal nodes), but node 6 is unary and node 7 is polytomous.
  # This used to reach the pairwise C++ pruning loop, which assumes every
  # internal node contributes exactly two consecutive edges.
  tree <- list(
    edge = matrix(c(
      5L, 1L,
      5L, 6L,
      6L, 7L,
      7L, 2L,
      7L, 3L,
      7L, 4L
    ), ncol = 2L, byrow = TRUE),
    edge.length = rep(1, 6L),
    Nnode = 3L,
    tip.label = c("a", "b", "c", "d")
  )
  class(tree) <- "phylo"
  x <- factor(c(a = "A", b = "A", c = "B", d = "B"))

  expect_error(
    fast_ace(x, tree, model = "ER", CI = FALSE, progress = FALSE),
    "rooted and fully dichotomous.*exactly two children"
  )
})

test_that("fast_ace topology guard accepts ordinary binary ER and ARD trees", {
  skip_if_not_installed("ape")
  set.seed(20260810L)
  tree <- ape::rtree(6L)
  x <- factor(sample(c("A", "B", "C"), 6L, replace = TRUE))
  names(x) <- tree$tip.label

  er <- suppressWarnings(fast_ace(
    x, tree, model = "ER", CI = FALSE, progress = FALSE
  ))
  ard <- suppressWarnings(fast_ace(
    x, tree, model = "ARD", CI = FALSE, progress = FALSE
  ))
  expect_s3_class(er, "ace")
  expect_s3_class(ard, "ace")
  expect_true(is.finite(er$loglik))
  expect_true(is.finite(ard$loglik))
})

test_that("fast_ace rejects zero branches and non-finite numeric states", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- c(a = 0, b = 0, c = 1, d = 1)

  zero <- tree
  zero$edge.length[[1L]] <- 0
  expect_error(
    fast_ace(x, zero, model = "ER", CI = FALSE, progress = FALSE),
    "strictly positive branch lengths"
  )
  # Raw-input positivity is checked before the kappa transform.  In
  # particular, kappa = 0 must not turn an input zero (0^0 == 1 in R) into a
  # silently accepted branch length.
  expect_error(
    fast_ace(x, zero, model = "ER", kappa = 0, CI = FALSE,
             progress = FALSE),
    "strictly positive branch lengths"
  )

  bad <- x
  bad[[2L]] <- Inf
  expect_error(
    fast_ace(bad, tree, model = "ER", CI = FALSE, progress = FALSE),
    "non-finite values"
  )
})
