test_that("D contrast sums are invariant to extreme branch scales", {
  tree <- ape::reorder.phylo(ape::read.tree(text =
    "(((a:0.3,b:0.8):0.6,c:0.4):0.7,(d:0.9,e:0.5):1.1);"), "postorder")
  states <- cbind(c(0, 1, 0, 1, 1), c(1, 1, 0, 0, 1))
  sums <- function(scale, threads = 1L) phylo_d_sums_cpp(
    states, tree$edge, tree$edge.length * scale, 5L, threads)
  expected <- sums(1)
  for (scale in c(1e200, 1e-200)) {
    expect_true(all(is.finite(sums(scale))))
    expect_equal(sums(scale), expected, tolerance = 1e-12)
    expect_identical(sums(scale, 2L), sums(scale))
  }
  stream <- function(scale) {
    set.seed(92301)
    result <- phylo_d_stream_cpp(states[, 1], tree$edge,
      tree$edge.length * scale, 5L, 37L, mean(states[, 1]), 32L, TRUE, 1L)
    list(result = result, seed = .Random.seed)
  }
  ordinary <- stream(1)
  for (scale in c(1e200, 1e-200)) {
    actual <- stream(scale)
    expect_identical(actual$seed, ordinary$seed)
    expect_equal(actual$result$observed, ordinary$result$observed, tolerance = 1e-12)
    expect_equal(actual$result$random, ordinary$result$random, tolerance = 1e-12)
    expect_equal(actual$result$brownian, ordinary$result$brownian, tolerance = 1e-12)
  }
})

test_that("K retains finite ratios when squared traits exceed double range", {
  tree <- ape::read.tree(text =
    "((a:0.3,b:0.8):0.6,(c:0.4,d:0.9):1.1);")
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 4L)
  x <- cbind(c(0, 2, -1, 4), c(3, 1, 0, -2))
  expected <- fast_k_tree_batch_cpp(core, x, 1L)
  for (scale in c(1e160, 1e-170)) {
    actual <- fast_k_tree_batch_cpp(core, x * scale, 1L)
    expect_true(all(is.finite(actual$K)))
    expect_equal(actual$K, expected$K, tolerance = 1e-12)
    expect_equal(actual$gls_mean / scale, expected$gls_mean, tolerance = 1e-12)
    expect_equal(fast_k_tree_batch_cpp(core, x * scale, 2L), actual,
      tolerance = 1e-12)
  }
  # 1'Qx exceeds double range, but the GLS mean and K remain finite.
  tiny_core <- compile_tree_cpp(tree$edge, tree$edge.length * 1e-200, 4L)
  joint <- fast_k_tree_batch_cpp(tiny_core, x * 1e160, 1L)
  expect_true(all(is.finite(joint$K)))
  expect_true(all(is.finite(joint$gls_mean)))
  expect_equal(joint$K, expected$K, tolerance = 1e-12)
  expect_equal(joint$gls_mean / 1e160, expected$gls_mean, tolerance = 1e-12)
  expect_equal(fast_k_tree_batch_cpp(tiny_core, x * 1e160, 2L), joint,
               tolerance = 1e-12)
  expect_true(all(is.nan(fast_k_tree_batch_cpp(core, matrix(1, 4, 1))$K)))
})

test_that("D falls back to finite weights for subnormal branch lengths", {
  tree <- ape::reorder.phylo(ape::read.tree(text =
    "(((a:0.3,b:0.8):0.6,c:0.4):0.7,(d:0.9,e:0.5):1.1);"), "postorder")
  states <- cbind(c(0, 1, 0, 1, 1), c(1, 1, 0, 0, 1))
  sums <- function(scale, threads = 1L) phylo_d_sums_cpp(
    states, tree$edge, tree$edge.length * scale, 5L, threads)
  for (scale in c(1e-308, 1e-309)) {
    expect_equal(sums(scale), sums(1), tolerance = 1e-12)
    expect_identical(sums(scale, 2L), sums(scale))
  }
  run <- function(scale) {
    scaled <- tree
    scaled$edge.length <- tree$edge.length * scale
    set.seed(42)
    fit <- fast_d(scaled, setNames(states[, 1], tree$tip.label),
      nsim = 19, return_sim = TRUE, verbose = FALSE, progress = FALSE)
    list(fit = fit, seed = .Random.seed)
  }
  ordinary <- run(1)
  tiny <- run(1e-308)
  expect_identical(tiny$seed, ordinary$seed)
  expect_equal(tiny$fit$estimate, ordinary$fit$estimate, tolerance = 1e-12)
  expect_identical(tiny$fit$status, ordinary$fit$status)
  expect_true(is.finite(tiny$fit$estimate))
  expect_equal(tiny$fit$nsim_successful_random, 19)
  expect_equal(tiny$fit$nsim_successful_brownian, 19)
})

test_that("K normalization stays finite when only the trace exceeds double range", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- setNames(c(0, 2, -1, 4), tree$tip.label)
  ordinary <- fast_k(tree, x, verbose = FALSE, progress = FALSE)
  tree$edge.length <- tree$edge.length * 5e307
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 4L)
  fitted <- fast_k_tree_batch_cpp(core, matrix(x, 4L), 1L)
  expect_true(is.finite(fitted$normalization))
  expect_equal(fitted$normalization / 5e307, 5/3, tolerance = 1e-12)
  actual <- fast_k(tree, x, verbose = FALSE, progress = FALSE)
  expect_equal(as.numeric(actual), as.numeric(ordinary), tolerance = 1e-12)
  set.seed(42)
  tested <- fast_k(tree, x, test = TRUE, nsim = 19,
    verbose = FALSE, progress = FALSE)
  expect_equal(as.numeric(actual), tested$K, tolerance = 1e-12)
})