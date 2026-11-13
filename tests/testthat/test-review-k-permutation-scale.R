test_that("public K permutation results are invariant to extreme trait scales", {
  set.seed(1)
  tree <- ape::rtree(8)
  x <- setNames(rnorm(8), tree$tip.label)
  run <- function(tree, x, permutations = NULL, ncores = 1L) {
    set.seed(42)
    fit <- fast_k(tree, x, test = TRUE, nsim = 19,
      permutations = permutations, return_sim = TRUE, ncores = ncores,
      simulation_chunk = 7L, verbose = FALSE, progress = FALSE)
    list(fit = fit, seed = .Random.seed)
  }
  set.seed(71)
  fixed <- t(replicate(19, sample.int(8)))
  for (perms in list(NULL, fixed)) {
    ordinary <- run(tree, x, perms)
    for (branch_scale in c(1, 1e-200)) {
      scaled_tree <- tree
      scaled_tree$edge.length <- tree$edge.length * branch_scale
      for (trait_scale in c(1e160, 1e-170)) {
        scaled_x <- x * trait_scale
        observed <- as.numeric(fast_k(scaled_tree, scaled_x,
          verbose = FALSE, progress = FALSE))
        for (threads in c(1L, 2L)) {
          actual <- run(scaled_tree, scaled_x, perms, threads)
          expect_true(is.finite(actual$fit$K))
          expect_equal(actual$fit$K, observed, tolerance = 1e-11)
          expect_equal(actual$fit$K, ordinary$fit$K, tolerance = 1e-11)
          expect_length(actual$fit$sim.K, 19L)
          expect_true(all(is.finite(actual$fit$sim.K)))
          expect_equal(actual$fit$sim.K, ordinary$fit$sim.K, tolerance = 1e-11)
          expect_equal(actual$fit$P, ordinary$fit$P)
          expect_equal(actual$fit$nsim_successful, 19)
          expect_identical(actual$seed, ordinary$seed)
        }
      }
    }
  }
})

test_that("K permutation chunk paths retain extreme-scale energy ratios", {
  set.seed(812)
  tree <- ape::rtree(8)
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 8L)
  x <- matrix(rnorm(8 * 11), 8, 11)
  perms <- t(replicate(17, sample.int(8)))
  run <- function(values, chunk, threads) fast_k_tree_permutation_cpp(
    core, values, nsim = 17L, permutations = perms,
    trait_chunk = chunk, simulation_chunk = 5L,
    return_sim = TRUE, n_threads = threads)
  ordinary <- run(x, 1L, 1L)
  for (scale in c(1e160, 1e-170)) {
    for (chunk in c(1L, 4L, 64L)) for (threads in c(1L, 2L)) {
      actual <- run(x * scale, chunk, threads)
      expect_true(all(is.finite(actual$K)))
      expect_true(all(is.finite(actual$sim_K)))
      expect_equal(actual$K, ordinary$K, tolerance = 1e-11)
      expect_equal(actual$sim_K, ordinary$sim_K, tolerance = 1e-11)
      expect_equal(actual$P, ordinary$P)
      expect_equal(actual$nsim_successful, ordinary$nsim_successful)
    }
  }
})