test_that("K joint scale invariance does not require extended long double", {
  tree <- ape::read.tree(text = "((a:0.3,b:0.8):0.6,(c:0.4,d:0.9):1.1);")
  x <- cbind(c(0, 2, -1, 4), c(3, 1, 0, -2))
  perms <- rbind(1:4, 4:1, c(2:4, 1), c(3, 4, 1, 2))
  storage.mode(perms) <- "integer"
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 4L)
  reference <- fast_k_tree_batch_cpp(core, x, 1L)
  null_reference <- fast_k_tree_permutation_cpp(core, x, 4L, perms,
    1L, TRUE, TRUE, 1L, 2L)
  for (bs in c(1, 1e-200, 1e200)) for (xs in c(1, 1e160, 1e-170)) {
    scaled_core <- compile_tree_cpp(tree$edge, tree$edge.length * bs, 4L)
    value <- fast_k_tree_batch_cpp(scaled_core, x * xs, 1L)
    null <- fast_k_tree_permutation_cpp(scaled_core, x * xs, 4L, perms,
      1L, TRUE, TRUE, 1L, 2L)
    expect_equal(value$K, reference$K, tolerance = 1e-12)
    expect_equal(value$gls_mean / xs, reference$gls_mean, tolerance = 1e-12)
    expect_equal(value$normalization / bs, reference$normalization, tolerance = 1e-12)
    expect_equal(value$sum_inv * bs, reference$sum_inv, tolerance = 1e-12)
    expect_equal(null$K, null_reference$K, tolerance = 1e-12)
    expect_equal(null$sim_K, null_reference$sim_K, tolerance = 1e-12)
    expect_identical(null$P, null_reference$P)
    expect_equal(null$normalization / bs, null_reference$normalization, tolerance = 1e-12)
    expect_equal(null$sum_inv * bs, null_reference$sum_inv, tolerance = 1e-12)
  }
})

test_that("public K scaling preserves nulls and RNG across thread policies", {
  set.seed(711)
  tree <- ape::rtree(20)
  x <- matrix(rnorm(40), 20, 2, dimnames = list(tree$tip.label, c("a", "b")))
  run <- function(branch_scale, trait_scale, threads) {
    tr <- tree; tr$edge.length <- tr$edge.length * branch_scale
    set.seed(712)
    fit <- fast_k(tr, x * trait_scale, test = TRUE, nsim = 19,
      return_sim = TRUE, ncores = threads, verbose = FALSE, progress = FALSE)
    list(fit = fit, seed = .Random.seed)
  }
  reference <- run(1, 1, 1L)
  for (threads in c(1L, 2L)) for (bs in c(1e-200, 1e200)) {
    actual <- run(bs, 1e160, threads)
    expect_equal(actual$fit$K_fast, reference$fit$K_fast, tolerance = 1e-12)
    expect_equal(actual$fit$P_fast, reference$fit$P_fast, tolerance = 1e-12)
    expect_identical(actual$seed, reference$seed)
  }
})

test_that("K handles positive subnormal edges and restores finite normalization", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- matrix(c(0, 2, -1, 4), 4, 1)
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 4L)
  reference <- fast_k_tree_batch_cpp(core, x, 1L)
  for (bs in c(1e-309, 5e307)) {
    scaled <- compile_tree_cpp(tree$edge, tree$edge.length * bs, 4L)
    actual <- fast_k_tree_batch_cpp(scaled, x, 1L)
    expect_equal(actual$K, reference$K, tolerance = 1e-12)
    expect_equal(actual$normalization / bs, reference$normalization, tolerance = 1e-12)
  }
})

test_that("K BM null scale recovery consumes the same random stream", {
  tree <- ape::read.tree(text = "((a:0.3,b:0.8):0.6,(c:0.4,d:0.9):1.1);")
  run <- function(bs) {
    core <- compile_tree_cpp(tree$edge, tree$edge.length * bs, 4L)
    set.seed(713)
    result <- heterogeneity_bm_k_cpp(core, 19L, 4L)
    list(result = result, seed = .Random.seed)
  }
  ordinary <- run(1)
  for (bs in c(1e-200, 1e200)) {
    actual <- run(bs)
    expect_equal(actual$result, ordinary$result, tolerance = 1e-12)
    expect_identical(actual$seed, ordinary$seed)
  }
})

test_that("native permutation rejects parent ids before indexing", {
  tree <- ape::read.tree(text="((a:1,b:1):1,(c:1,d:1):1);")
  core <- compile_tree_cpp(tree$edge,tree$edge.length,4L)
  X <- matrix(c(0,2,-1,4),4L)
  for (bad in c(-2147483647L, -1L, 0L, 2147483647L)) {
    altered <- core
    altered$parent[1L] <- bad
    expect_error(fast_k_tree_permutation_cpp(altered,X,5L),"parent.*invalid")
  }
})