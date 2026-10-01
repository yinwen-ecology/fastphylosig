# Stage 5 release-tarball regression coverage for the public fast_k() path.
# The selector is inspected only to prove that the two public calls exercise
# the intended production branches; it remains an internal, unexported helper.

.stage5_tree <- function(n_tip = 24L, seed = 50501L) {
  set.seed(seed)
  tree <- ape::rtree(n_tip)
  tree$edge.length <- tree$edge.length / max(tree$edge.length)
  tree
}

.stage5_trait <- function(tree, seed = 50502L) {
  set.seed(seed)
  x <- matrix(stats::rnorm(length(tree$tip.label)), ncol = 1L,
              dimnames = list(tree$tip.label, "trait"))
  x
}

.stage5_strategy <- function(n_tip, nsim, memory_budget_mb) {
  fastphylosig:::k_strategy_cpp(
    n_tip = as.integer(n_tip), n_trait = 1L, nsim = as.integer(nsim),
    memory_budget_mb = memory_budget_mb, trait_chunk = 1L,
    n_total = 0L, n_threads = 1L
  )
}

.stage5_run_public <- function(tree, x, nsim, seed, memory_budget_mb,
                               permutations = NULL) {
  old_options <- options(fastphylosig.k_block_memory_mb = memory_budget_mb)
  on.exit(options(old_options), add = TRUE)
  set.seed(seed)
  fit <- suppressMessages(fastphylosig::fast_k(
    tree = tree, x = x, test = TRUE, nsim = nsim, return_sim = TRUE,
    permutations = permutations,
    verbose = FALSE, ncores = 1L, trait_chunk = 1L,
    simulation_chunk = 7L, progress = FALSE
  ))
  list(
    K = as.numeric(fit$K_fast[[1L]]),
    P = as.numeric(fit$P_fast[[1L]]),
    exceedance = as.numeric(fit$exceedance_count[[1L]]),
    nsim = as.integer(fit$nsim_successful[[1L]]),
    K_null = as.numeric(fit$sim.K_fast[[1L]]),
    rng_state = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  )
}

.stage5_tail <- function(values, observed) {
  values <- as.numeric(values)
  direct <- values >= observed
  scale <- pmax(1, abs(values), abs(observed))
  tie <- is.finite(values) & is.finite(observed) &
    (observed - values <= 8 * .Machine$double.eps * scale)
  direct | tie
}

.stage5_run_dense_reference <- function(tree, x, permutations) {
  suppressMessages(fastphylosig:::.fast_signal_with_engine(
    tree = tree, x = x, method = "K", test = TRUE,
    nsim = nrow(permutations), permutations = permutations,
    return_sim = TRUE, verbose = FALSE, ncores = 1L,
    engine = "dense", trait_chunk = 1L, simulation_chunk = 7L,
    keep_null = TRUE
  ))
}

test_that("the release tarball exercises fusion-on and fusion-off through fast_k", {
  nsim <- 37L
  tree <- .stage5_tree()
  x <- .stage5_trait(tree)

  # The 64 MiB default leaves ample room for a fused block at this tree size.
  fused_strategy <- .stage5_strategy(
    length(tree$tip.label), nsim, memory_budget_mb = 64
  )
  expect_true(fused_strategy$hybrid)
  expect_true(fused_strategy$fuse)
  expect_gt(fused_strategy$block, 1L)

  # The one-byte budget is below even the one-permutation workspace, forcing
  # the same production selector to retain hybrid evaluation and disable fusion.
  unfused_strategy <- .stage5_strategy(
    length(tree$tip.label), nsim, memory_budget_mb = 1e-6
  )
  expect_true(unfused_strategy$hybrid)
  expect_false(unfused_strategy$fuse)
  expect_identical(unfused_strategy$block, 1L)

  fused <- .stage5_run_public(
    tree, x, nsim, seed = 20260929L, memory_budget_mb = 64
  )
  unfused <- .stage5_run_public(
    tree, x, nsim, seed = 20260929L, memory_budget_mb = 1e-6
  )

  expect_true(all(is.finite(c(fused$K, fused$K_null))))
  expect_identical(length(fused$K_null), nsim)
  expect_identical(length(unfused$K_null), nsim)
  expect_lte(abs(fused$K - unfused$K), 1e-12)
  expect_equal(fused$K_null, unfused$K_null, tolerance = 1e-12)
  expect_identical(fused$P, unfused$P)
  expect_identical(fused$exceedance, unfused$exceedance)
  expect_identical(
    .stage5_tail(fused$K_null, fused$K),
    .stage5_tail(unfused$K_null, unfused$K)
  )
  expect_identical(fused$exceedance,
                   as.numeric(sum(.stage5_tail(fused$K_null, fused$K))))
  expect_identical(fused$P, fused$exceedance / nsim)
  expect_identical(fused$rng_state, unfused$rng_state)
})

test_that("fixed-seed star-tree K ties preserve the inclusive hybrid tail", {
  n_tip <- 12L
  nsim <- 41L
  tree_text <- paste0(
    "(", paste(sprintf("s%02d:1", seq_len(n_tip)), collapse = ","), ");"
  )
  tree <- ape::read.tree(text = tree_text)
  x <- matrix(seq_len(n_tip) - mean(seq_len(n_tip)), ncol = 1L,
              dimnames = list(tree$tip.label, "trait"))

  # Equal terminal branches make the star-tree K statistic permutation
  # invariant in exact arithmetic.  A fixed permutation seed therefore gives
  # a deterministic numerical boundary fixture: every null statistic is an
  # inclusive-tail tie with K_obs, including the identity row.
  set.seed(739201L)
  permutations <- t(replicate(nsim - 1L, sample.int(n_tip)))
  permutations <- rbind(seq_len(n_tip), permutations)

  hybrid_strategy <- .stage5_strategy(
    n_tip, nsim, memory_budget_mb = 1e-6
  )
  expect_true(hybrid_strategy$hybrid)
  expect_false(hybrid_strategy$fuse)

  candidate <- .stage5_run_public(
    tree, x, nsim, seed = 739202L, memory_budget_mb = 1e-6,
    permutations = permutations
  )
  reference <- .stage5_run_dense_reference(tree, x, permutations)
  reference_K <- as.numeric(reference$K_fast[[1L]])
  reference_null <- as.numeric(reference$sim.K_fast[[1L]])
  candidate_tail <- .stage5_tail(candidate$K_null, candidate$K)
  reference_tail <- .stage5_tail(reference_null, reference_K)

  # Stage 4's worst observed hybrid/reference difference was 6.03e-14;
  # 1e-12 keeps a small cross-platform rounding margin while remaining far
  # below any meaningful change to K or to the 8-epsilon inclusive-tail rule.
  expect_lte(abs(candidate$K - reference_K), 1e-12)
  expect_equal(candidate$K_null, reference_null, tolerance = 1e-12)
  expect_identical(candidate_tail, reference_tail)
  expect_true(all(candidate_tail))
  expect_true(all(reference_tail))
  expect_identical(candidate$exceedance, as.numeric(nsim))
  expect_identical(candidate$P, 1)
  expect_identical(as.numeric(reference$P_fast[[1L]]), 1)
})

test_that("Stage 5 strategy diagnostics do not add a public API", {
  expect_setequal(getNamespaceExports("fastphylosig"), c(
    "fast_signal", "fast_k", "fast_lambda", "fast_d", "fast_delta",
    "fast_ace", "prepare_tree", "cache_info", "match_tree_data",
    "match_phylo_data", "plot_signal", "check_tree", "resolve_tree"
  ))
  expect_identical(names(formals(fastphylosig::fast_k)), c(
    "tree", "x", "test", "nsim", "se", "permutations", "return_sim",
    "verbose", "ncores", "X", "prepared", "trait_chunk",
    "simulation_chunk", "keep_null", "progress", "data"
  ))
  expect_false("k_strategy_cpp" %in% getNamespaceExports("fastphylosig"))
})
