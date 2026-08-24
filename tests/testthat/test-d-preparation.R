test_that("D preparation preserves raw/prepared controlled null results and legacy fields", {
  testthat::skip_if_not_installed("ape")
  set.seed(7401)
  tree <- ape::rtree(8L)
  x <- stats::setNames(rep(0:1, length.out = 8L), tree$tip.label)
  nsim <- 4L
  random_states <- matrix(rep(x, nsim), nrow = 8L)
  brownian_states <- matrix(rep(1 - x, nsim), nrow = 8L)

  raw <- fast_d(tree, x, nsim = nsim, test = TRUE, return_sim = TRUE,
                random_states = random_states,
                brownian_states = brownian_states,
                verbose = FALSE, progress = FALSE)
  ctx <- prepare_tree(tree)
  cached <- fast_d(ctx, x, nsim = nsim, test = TRUE, return_sim = TRUE,
                   random_states = random_states,
                   brownian_states = brownian_states,
                   verbose = FALSE, progress = FALSE)

  expect_equal(cached$DEstimate, raw$DEstimate, tolerance = 1e-12)
  expect_equal(cached$Pval1, raw$Pval1, tolerance = 1e-12)
  expect_equal(cached$Pval0, raw$Pval0, tolerance = 1e-12)
  expect_equal(cached$nsim_requested, raw$nsim_requested)
  expect_equal(cached$nsim_successful_random, raw$nsim_successful_random)
  expect_equal(cached$nsim_successful_brownian, raw$nsim_successful_brownian)
  expect_true(is.list(attr(raw, "analysis_metadata")))
  expect_true(is.list(attr(cached, "analysis_metadata")))
  expect_true(all(c("matching", "na_pattern_count", "na_pattern_keys") %in%
                    names(attr(raw, "analysis_metadata"))))
})

test_that("D canonicalizes a representation-only noncanonical root", {
  testthat::skip_if_not_installed("ape")
  tree <- list(
    edge = matrix(c(
      6L, 5L, 6L, 7L, 5L, 1L,
      5L, 2L, 7L, 3L, 7L, 4L
    ), ncol = 2L, byrow = TRUE),
    edge.length = rep(1, 6L), Nnode = 3L,
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
  fit <- fast_d(tree, x, nsim = 2L, random_states = random_states,
                brownian_states = brownian_states, return_sim = FALSE,
                verbose = FALSE, progress = FALSE)
  expect_s3_class(fit, "phylo.d")
  prep <- attr(fit, "analysis_metadata")
  expect_true(isTRUE(prep$tree_auto_normalized))
  expect_equal(prep$signal, "D")
})

test_that("D marks an unsafe retained NA subset instead of calling the kernel", {
  testthat::skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "(((a:1,b:1):1,c:1):1,d:1);")
  x <- data.frame(
    trait_a = c(0, 0, 1, NA),
    trait_b = c(0, NA, NA, NA),
    row.names = tree$tip.label
  )
  fit <- fast_d(tree, x, nsim = 2L, test = FALSE, return_sim = FALSE,
                verbose = FALSE, progress = FALSE)
  expect_equal(nrow(fit), 2L)
  expect_true(any(grepl("fewer than 2|not ready|check_tree",
                        fit$note, ignore.case = TRUE), na.rm = TRUE))
  prep <- attr(fit, "analysis_metadata")
  expect_true(is.list(prep))
  expect_true(length(prep$retained_tip_validation) >= 1L)
  expect_true(any(!prep$retained_tip_validation))
})
