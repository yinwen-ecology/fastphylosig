# Stage 4 production strategy contract and RNG invariance.

.skip_without_selector <- function() {
  skip_if_not(exists("k_strategy_cpp", asNamespace("fastphylosig")),
              "strategy selector is not present in this build")
}

.k_strategy <- function(n_tip, n_trait, nsim, budget_mb = 64,
                        trait_chunk = 0L, n_total = 0L, n_threads = 1L) {
  fastphylosig:::k_strategy_cpp(
    n_tip = as.integer(n_tip), n_trait = as.integer(n_trait),
    nsim = as.integer(nsim), memory_budget_mb = budget_mb,
    trait_chunk = as.integer(trait_chunk), n_total = as.integer(n_total),
    n_threads = as.integer(n_threads)
  )
}

.k_test_inputs <- function(n_tip = 120L, n_trait = 3L, seed = 20260924L) {
  set.seed(seed)
  tree <- ape::rtree(n_tip)
  tree$edge.length <- tree$edge.length / max(tree$edge.length)
  X <- matrix(stats::rnorm(n_tip * n_trait), nrow = n_tip, ncol = n_trait)
  colnames(X) <- sprintf("t%02d", seq_len(n_trait))
  rownames(X) <- tree$tip.label
  list(tree = tree, X = X)
}

test_that("the selector always keeps the hybrid path enabled", {
  .skip_without_selector()
  for (nsim in c(0L, 1L, 2L, 15L, 16L, 1000L)) {
    s <- .k_strategy(5000L, 1L, nsim)
    expect_true(s$hybrid, info = paste("nsim", nsim))
  }
})

test_that("fusion follows the conservative tree-size and total-trait window", {
  .skip_without_selector()
  cases <- data.frame(
    n = c(50L, 500L, 501L, 5000L, 10000L, 500L, 501L, 10000L, 10001L),
    p = c(10L, 10L, 5L, 5L, 5L, 11L, 6L, 6L, 1L),
    fuse = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
  )
  for (i in seq_len(nrow(cases))) {
    s <- .k_strategy(cases$n[[i]], cases$p[[i]], 1000L)
    expect_true(s$hybrid, info = paste("hybrid", i))
    expect_identical(isTRUE(s$fuse), cases$fuse[[i]],
                     info = paste("fusion window case", i))
    if (!cases$fuse[[i]]) expect_equal(s$block, 1L)
  }

  # The gate is based on total n_trait, even when the evaluator's working chunk
  # is narrower.  A one-trait chunk must not bypass the p=10 safety boundary.
  narrow_chunk <- .k_strategy(500L, 11L, 1000L, trait_chunk = 1L)
  expect_false(narrow_chunk$fuse)
  expect_equal(narrow_chunk$chunk, 1L)
})

test_that("fusion requires a one-permutation workspace and obeys block caps", {
  .skip_without_selector()

  too_small <- .k_strategy(500L, 1L, 1000L, budget_mb = 1e-6)
  expect_true(too_small$hybrid)
  expect_false(too_small$fuse)
  expect_equal(too_small$block, 1L)

  # A very large workspace exposes the explicit p-based cap rather than the
  # memory ladder's upper rung.
  expect_equal(.k_strategy(50L, 1L, 1000L, budget_mb = 1e6)$block, 32L)
  expect_equal(.k_strategy(50L, 2L, 1000L, budget_mb = 1e6)$block, 32L)
  expect_equal(.k_strategy(50L, 3L, 1000L, budget_mb = 1e6)$block, 4L)
  expect_equal(.k_strategy(50L, 10L, 1000L, budget_mb = 1e6)$block, 4L)
  expect_equal(.k_strategy(50L, 11L, 1000L, budget_mb = 1e6)$block, 1L)

  for (mb in c(0.5, 2, 8, 64)) {
    for (n in c(100L, 1000L, 10000L)) {
      for (p in c(1L, 2L, 5L, 10L)) {
        s <- .k_strategy(n, p, 1000L, budget_mb = mb)
        if (s$fuse) {
          expect_gt(s$block, 1L)
          expect_lte(s$workspace_mb, s$budget_mb * (1 + 1e-9))
        } else {
          expect_equal(s$block, 1L)
        }
      }
    }
  }
})

test_that("changing fusion block preserves K, nulls, tail counts, and RNG state", {
  .skip_without_selector()
  testthat::skip_if_not_installed("ape")
  inp <- .k_test_inputs(n_tip = 150L, n_trait = 2L)
  nsim <- 199L

  old <- getOption("fastphylosig.k_block_memory_mb")
  on.exit(options(fastphylosig.k_block_memory_mb = old), add = TRUE)

  options(fastphylosig.k_block_memory_mb = 64)
  s_big <- .k_strategy(150L, 2L, nsim, budget_mb = 64)
  set.seed(4242L)
  big <- fast_k(inp$tree, inp$X, test = TRUE, nsim = nsim,
                return_sim = TRUE, verbose = FALSE, progress = FALSE)
  seed_big <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)

  options(fastphylosig.k_block_memory_mb = 1e-6)
  s_small <- .k_strategy(150L, 2L, nsim, budget_mb = 1e-6)
  set.seed(4242L)
  small <- fast_k(inp$tree, inp$X, test = TRUE, nsim = nsim,
                  return_sim = TRUE, verbose = FALSE, progress = FALSE)
  seed_small <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)

  expect_true(s_big$fuse && s_big$block > 1L)
  expect_false(s_small$fuse)
  expect_equal(s_small$block, 1L)
  expect_identical(as.numeric(big$P_fast), as.numeric(small$P_fast))
  expect_identical(as.numeric(big$exceedance_count),
                   as.numeric(small$exceedance_count))
  expect_equal(as.numeric(big$K_fast), as.numeric(small$K_fast),
               tolerance = 1e-10)
  expect_equal(unlist(big[["sim.K_fast"]]), unlist(small[["sim.K_fast"]]),
               tolerance = 1e-10)
  expect_identical(seed_big, seed_small)
})

test_that("the RNG stream is independent of simulation chunk size", {
  .skip_without_selector()
  inp <- .k_test_inputs(n_tip = 100L, n_trait = 2L)
  nsim <- 99L
  set.seed(1234L)
  a <- fast_k(inp$tree, inp$X, test = TRUE, nsim = nsim, return_sim = TRUE,
              simulation_chunk = 16L, verbose = FALSE, progress = FALSE)
  seed_a <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  set.seed(1234L)
  b <- fast_k(inp$tree, inp$X, test = TRUE, nsim = nsim, return_sim = TRUE,
              simulation_chunk = 512L, verbose = FALSE, progress = FALSE)
  seed_b <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)

  expect_identical(seed_a, seed_b)
  expect_identical(as.numeric(a$P_fast), as.numeric(b$P_fast))
  expect_identical(as.numeric(a$exceedance_count),
                   as.numeric(b$exceedance_count))
  expect_equal(unlist(a[["sim.K_fast"]]), unlist(b[["sim.K_fast"]]),
               tolerance = 1e-10)
})

test_that("the selector remains internal to the package", {
  expect_false("k_strategy_cpp" %in% getNamespaceExports("fastphylosig"))
  expect_false(any(grepl("strategy", getNamespaceExports("fastphylosig"))))
})
