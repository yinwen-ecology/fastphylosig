test_that("prepare_tree reuses tree calculations without changing K", {
  set.seed(101)
  tree <- ape::rtree(18)
  X <- matrix(stats::rnorm(36), nrow = 18,
              dimnames = list(tree$tip.label, c("a", "b")))
  ctx <- prepare_tree(tree)
  initial <- cache_info(ctx)
  expect_equal(initial$n_numerical_entries, 0L)
  expect_equal(initial$bytes_used, 0)

  expect_s3_class(ctx, "fastphylosig_tree")
  expect_s3_class(ctx$inspection, "fastphylosig_tree_check")
  expect_named(ctx$inspection$ready_by_signal,
               c("K", "lambda", "D", "Delta"))
  raw_check <- check_tree(tree, signal = "K")
  cached_check <- check_tree(ctx, signal = "K")
  expect_identical(cached_check$ready_by_signal, raw_check$ready_by_signal)
  expect_identical(cached_check$issues$code, raw_check$issues$code)
  expect_true(cached_check$tree_summary$prepared)
  raw <- fast_signal(tree, X, method = "K", verbose = FALSE)
  cached <- fast_signal(ctx, X, method = "K", verbose = FALSE)
  after_k <- cache_info(ctx)
  expect_equal(after_k$n_numerical_entries, 0L)
  expect_false(any(after_k$entries$has_vcv))
  expect_equal(cached$K_fast, raw$K_fast, tolerance = 1e-12)
  expect_equal(cached$n_species, raw$n_species)
  expect_equal(cached$matched_species, raw$matched_species)

  ctx_bad <- ctx
  ctx_bad$tree$edge.length[[1L]] <- ctx_bad$tree$edge.length[[1L]] + 1
  expect_error(fast_signal(ctx_bad, X, method = "K", verbose = FALSE),
               "prepared tree was modified")

  ctx_bad_nnode <- ctx
  ctx_bad_nnode$tree$Nnode <- ctx_bad_nnode$tree$Nnode + 1L
  expect_error(check_tree(ctx_bad_nnode, signal = "K"),
               "prepared tree was modified")
  expect_error(fast_k(ctx_bad_nnode, X, test = FALSE, verbose = FALSE),
               "prepared tree was modified")
})

test_that("prepare_tree maps malformed structure to actionable diagnostics", {
  tree <- ape::rtree(5)
  tree$Nnode <- NA_real_
  expect_error(
    prepare_tree(tree),
    "Problem:|Nnode|full diagnosis"
  )
})

test_that("prepared contexts enforce selected zero-branch readiness", {
  testthat::skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:0,b:1):1,(c:1,d:1):1);")
  ctx <- prepare_tree(tree)

  continuous <- stats::setNames(stats::rnorm(4), tree$tip.label)
  expect_error(
    fast_k(ctx, continuous, progress = FALSE),
    "k_requires_positive_branches|strictly positive branch lengths"
  )
  expect_error(
    fast_lambda(ctx, continuous, progress = FALSE),
    "lambda_zero_terminal_branch|USER_ACTION_REQUIRED|positive branch"
  )

  binary <- stats::setNames(c(0, 0, 1, 1), tree$tip.label)
  expect_error(
    fast_d(ctx, binary, test = FALSE, return_sim = FALSE,
           verbose = FALSE, progress = FALSE),
    "d_requires_positive_branches|strictly positive branch lengths"
  )

  categorical <- stats::setNames(c("A", "A", "B", "B"), tree$tip.label)
  expect_error(
    fast_delta(ctx, categorical, mcmc_sim = 2, thin = 1, burn = 1,
               verbose = FALSE, progress = FALSE),
    "delta_requires_positive_branches|strictly positive branch lengths"
  )
})

test_that("lambda spectral cache matches dense likelihood and phytools", {
  set.seed(102)
  tree <- ape::rtree(16)
  x <- stats::rnorm(16)
  names(x) <- tree$tip.label
  ctx <- prepare_tree(tree)
  subset <- fastphylosig:::.prepared_tree_subset(ctx, need_lambda = TRUE)
  lambda <- c(0, 0.2, 0.7, 1)
  dense <- vapply(lambda, function(z) {
    fastphylosig:::lambda_loglik_cpp(z, subset$C, x)
  }, numeric(1))
  spectral <- vapply(lambda, function(z) {
    fastphylosig:::lambda_loglik_spectral_cpp(
      z, subset$lambda_spectral$values, subset$lambda_spectral$vectors,
      subset$lambda_spectral$inv_sqrt_diag,
      subset$lambda_spectral$log_diag, x
    )
  }, numeric(1))
  expect_equal(spectral, dense, tolerance = 1e-8)

  cached <- fast_signal(ctx, x, method = "lambda", test = TRUE,
                        verbose = FALSE)
  reference <- phytools::phylosig(tree, x, method = "lambda",
                                  test = TRUE, se = NULL)
  expect_equal(cached$lambda, unname(reference$lambda), tolerance = 5e-4)
  expect_equal(cached$logL, unname(reference$logL), tolerance = 1e-5)
  expect_equal(cached$P, unname(reference$P), tolerance = 1e-5)

  X <- cbind(a = x, b = -x)
  prof <- fast_signal(ctx, X, method = "lambda", lambda_profile = TRUE,
                      lambda_profile_points = 21, verbose = FALSE)
  expect_length(prof$lambda_profile_fast, 2)
  expect_true(all(vapply(prof$lambda_profile_fast, nrow, integer(1)) >= 21L))

  xc <- stats::setNames(rep(3, 16), tree$tip.label)
  fc <- suppressWarnings(
    fast_signal(tree, xc, method = "lambda", verbose = FALSE)
  )
  rc <- suppressWarnings(
    phytools::phylosig(tree, xc, method = "lambda", test = FALSE,
                       se = NULL)
  )
  # A constant trait has zero fitted variance, so lambda and its likelihood
  # are not statistically identified. The tree engine reports this explicitly
  # instead of preserving a platform-dependent dense-matrix round-off value.
  expect_true(is.na(fc$lambda))
  expect_identical(fc$logL, -Inf)
  expect_match(fc$note, "undefined|degenerate|variance", ignore.case = TRUE)
  expect_true(is.finite(unname(rc$lambda)) || is.na(unname(rc$lambda)))
})

test_that("D batch traversal is reproducible with a prepared tree", {
  set.seed(103)
  tree <- ape::rtree(14)
  X <- cbind(
    a = as.integer(stats::rbinom(14, 1, 0.45)),
    b = as.integer(stats::rbinom(14, 1, 0.55))
  )
  rownames(X) <- tree$tip.label
  nsim <- 7L
  random_states <- matrix(stats::runif(14 * nsim), 14, nsim)
  brownian_states <- matrix(stats::runif(14 * nsim), 14, nsim)

  raw <- fast_d(tree, X, nsim = nsim, test = TRUE,
                random_states = random_states,
                brownian_states = brownian_states,
                verbose = FALSE)
  ctx <- prepare_tree(tree)
  expect_equal(cache_info(ctx)$n_numerical_entries, 0L)
  cached <- fast_d(ctx, X, nsim = nsim, test = TRUE,
                   random_states = random_states,
                   brownian_states = brownian_states,
                   verbose = FALSE)
  expect_equal(cache_info(ctx)$n_numerical_entries, 0L)
  expect_equal(cached$D_fast, raw$D_fast, tolerance = 1e-12)
  expect_equal(cached$observed, raw$observed, tolerance = 1e-12)
  expect_equal(cached$mean_random, raw$mean_random, tolerance = 1e-12)
  expect_equal(cached$mean_brownian, raw$mean_brownian, tolerance = 1e-12)
  expect_equal(cached$random_fast, raw$random_fast, tolerance = 1e-12)
  expect_equal(cached$brownian_fast, raw$brownian_fast, tolerance = 1e-12)

  threaded <- fast_d(ctx, X, nsim = nsim, test = TRUE, ncores = 2,
                     random_states = random_states,
                     brownian_states = brownian_states,
                     verbose = FALSE)
  expect_equal(threaded$D_fast, raw$D_fast, tolerance = 1e-12)
  expect_equal(threaded$observed, raw$observed, tolerance = 1e-12)
})

test_that("tree preparation and ncores reject unsafe inputs", {
  set.seed(105)
  tree <- ape::rtree(6)
  ctx <- prepare_tree(tree)
  expect_null(fastphylosig:::.prepared_tree_subset(ctx)$lambda_spectral)
  ctx_small <- prepare_tree(tree, max_cached_subsets = 1L)
  fastphylosig:::.prepared_tree_subset(ctx_small, seq_len(5L))
  expect_length(ls(ctx_small$cache, all.names = TRUE), 1L)
  x <- stats::setNames(stats::rnorm(6), tree$tip.label)
  ctx_zero <- prepare_tree(tree, cache_budget = 0)
  fast_signal(ctx_zero, x, method = "K", verbose = FALSE)
  expect_equal(cache_info(ctx_zero)$n_numerical_entries, 0L)
  expect_equal(cache_info(ctx_zero)$bytes_used, 0)
  expect_error(fast_signal(ctx, x, method = "K", ncores = 1.5,
                           verbose = FALSE), "finite positive integer")
  two_tip <- ape::read.tree(text = "(a:1,b:1);")
  xx <- stats::setNames(c(0.2, -0.1), two_tip$tip.label)
  expect_s3_class(fast_signal(two_tip, xx, method = "lambda",
                              verbose = FALSE), "phylosig")
})

test_that("numerical cache obeys a byte-budget LRU", {
  set.seed(106)
  tree <- ape::rtree(12)
  x <- stats::setNames(stats::rnorm(12), tree$tip.label)

  probe <- prepare_tree(tree)
  fastphylosig:::.prepared_tree_subset(
    probe, need_lambda = TRUE, need_matrix = TRUE
  )
  one_entry_bytes <- cache_info(probe)$bytes_used
  expect_gt(one_entry_bytes, 0)

  ctx <- prepare_tree(tree, cache_budget = one_entry_bytes + 128)
  fastphylosig:::.prepared_tree_subset(
    ctx, need_lambda = TRUE, need_matrix = TRUE
  )
  x2 <- x[-1L]
  fastphylosig:::.prepared_tree_subset(
    ctx, keep = 2:12, need_lambda = TRUE, need_matrix = TRUE
  )
  info <- cache_info(ctx)
  expect_lte(info$bytes_used, info$budget)
  expect_equal(info$n_numerical_entries, 1L)
  expect_gte(info$evictions, 1L)
})

test_that("prepared contexts are accepted by fast_ace and fast_delta", {
  set.seed(104)
  tree <- ape::rtree(10)
  x <- sample(letters[1:3], ape::Ntip(tree), replace = TRUE)
  names(x) <- tree$tip.label
  ctx <- prepare_tree(tree)

  ace_raw <- suppressWarnings(
    fast_ace(x, tree, model = "ER", CI = FALSE)
  )
  ace_cached <- suppressWarnings(
    fast_ace(x, prepared = ctx, model = "ER", CI = FALSE)
  )
  expect_equal(cache_info(ctx)$n_numerical_entries, 0L)
  expect_equal(ace_cached$loglik, ace_raw$loglik, tolerance = 1e-8)
  expect_equal(ace_cached$rates, ace_raw$rates, tolerance = 1e-8)

  bad_tree <- tree
  bad_tree$edge.length[[1L]] <- Inf
  expect_error(fast_ace(x, bad_tree, model = "ER", CI = FALSE),
               "finite and non-negative")

  set.seed(105)
  delta_raw <- suppressWarnings(
    fast_delta(tree, x, test = FALSE, mcmc_sim = 40,
               thin = 5, burn = 10, verbose = FALSE)
  )
  set.seed(105)
  delta_cached <- suppressWarnings(
    fast_delta(ctx, x, test = FALSE, mcmc_sim = 40,
               thin = 5, burn = 10, verbose = FALSE)
  )
  expect_equal(cache_info(ctx)$n_numerical_entries, 0L)
  expect_s3_class(delta_cached, "phylo_delta")
  expect_equal(delta_cached$Delta, delta_raw$Delta, tolerance = 1e-12)
})

test_that("packed structural-cache keys support large trees", {
  key <- fastphylosig:::.tree_mask_key(seq_len(20000L), 20000L)
  expect_lt(nchar(key, type = "bytes"), 10000L)
  expect_identical(
    fastphylosig:::.tree_mask_key(c(9L, 2L, 1L), 20000L),
    fastphylosig:::.tree_mask_key(c(1L, 2L, 9L), 20000L)
  )

  set.seed(107)
  tree <- ape::rtree(5000L)
  ctx <- prepare_tree(tree)
  expect_s3_class(ctx, "fastphylosig_tree")
  expect_length(ls(ctx$structural_cache, all.names = TRUE), 1L)
})

test_that("large-tree analysis uses packed subset keys", {
  set.seed(108)
  tree <- ape::rtree(3000L)
  x <- stats::setNames(stats::rnorm(3000L), tree$tip.label)
  fit <- suppressMessages(suppressWarnings(fast_k(
    tree, x, test = FALSE, verbose = FALSE, progress = FALSE
  )))
  expect_s3_class(fit, "phylosig")
  expect_true(is.finite(as.numeric(fit)))
})
