test_that("native clade indexing and extraction preserve representations and crowns", {
  set.seed(110)
  trees <- list(ape::rtree(128), ape::stree(128, "balanced"),
                ape::stree(128, "star"), ape::stree(1024, "left"))
  for (tree in trees) {
    tree$edge.length <- rep(0.25, nrow(tree$edge))
    nt <- length(tree$tip.label)
    ids <- c(seq_len(nt), nt + rev(seq_len(tree$Nnode)))
    tree$edge[] <- ids[tree$edge]
    tree$edge <- tree$edge[nrow(tree$edge):1L, , drop = FALSE]
    core <- compile_tree_cpp(tree$edge, tree$edge.length, nt)
    index <- .hetero_index_reference(core)
    selected <- head(core$preorder[core$preorder > nt], 3L)
    plan <- heterogeneity_plan_cpp(core, tree$tip.label, selected, 3L, nt, 10000)
    expect_identical(plan$nodes, selected)
    expect_identical(plan$n_tip, index$n_tip[selected])
    expect_equal(plan$total, sum(index$n_tip[selected]))
    for (i in seq_along(selected)) {
      ref <- .hetero_subtree_reference(tree, core, index, selected[[i]])
      tip_ids <- attr(ref, "source_tip_index")
      attr(ref, "source_tip_index") <- NULL
      expect_identical(plan$trees[[i]], ref)
      expect_identical(plan$tips[[i]], tip_ids)
    }
  }
})

test_that("native scan shares exact NA patterns and matches the checked R kernels", {
  set.seed(144)
  tree <- ape::rtree(80)
  X <- matrix(rnorm(80 * 6), 80, 6,
    dimnames = list(tree$tip.label, paste0("trait", 1:6)))
  X[, 6] <- 1e12 + X[, 1]
  X[c(1:15, 30:45), 2:3] <- NA_real_
  X[c(16:29, 46:57), 4] <- NA_real_
  X[1, 5] <- NaN
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 80L)
  plan <- heterogeneity_plan_cpp(core, tree$tip.label, core$root, 3L, 80L, 10000)
  for (method in c("K", "lambda")) {
    fit <- heterogeneity_scan_cpp(plan, X, method, 3L, TRUE)
    ref <- .hetero_fit_reference(plan$trees[[1]], X, method)
    expect_equal(fit$estimate, ref$estimate, tolerance = 1e-12)
    expect_identical(fit$status, ref$status)
    expect_identical(fit$message, ref$message)
    expect_identical(fit$cores[[1]][[2]], fit$cores[[1]][[3]])
  }
})

test_that("native BM K preserves draws across RNG kinds and block widths", {
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  tree <- ape::stree(17, "star") # Odd edge count challenges Box-Muller caching.
  tree$edge.length <- seq(0.5, 2.1, length.out = 17)
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 17L)
  for (uniform in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
    for (normal in c("Inversion", "Box-Muller")) {
      RNGkind(uniform, normal)
      set.seed(125)
      ref <- fast_k_tree_batch_cpp(core, .hetero_bm_reference(core, 11L), 11L)$K
      seed <- .Random.seed
      for (chunk in c(1L, 3L, 11L)) {
        set.seed(125)
        fit <- heterogeneity_bm_k_cpp(core, 11L, chunk)
        expect_identical(.Random.seed, seed)
        expect_identical(fit, ref)
      }
    }
  }
})

test_that("memory-capped BM blocks retain all replicates and RNG order", {
  set.seed(555)
  tree <- ape::rtree(6000L)
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 6000L)
  # 401 samples exceed the native 64 MiB block budget on this tree.
  set.seed(77)
  small <- heterogeneity_bm_k_cpp(core, 401L, 7L)
  seed <- .Random.seed
  set.seed(77)
  large <- heterogeneity_bm_k_cpp(core, 401L, .Machine$integer.max)
  expect_identical(.Random.seed, seed)
  expect_identical(large, small)
  expect_length(large, 401L)
})

test_that("native calibration retains probabilities and reference envelopes", {
  tree <- ape::stree(16, "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 16L)
  set.seed(888)
  fit <- heterogeneity_calibrate_cpp(list(core, core), 0.5, 19L, 4L, TRUE)
  seed <- .Random.seed
  set.seed(888)
  draws <- lapply(1:2, function(i) fast_k_tree_batch_cpp(core, .hetero_bm_reference(core, 19L), 19L)$K)
  H <- apply(log(do.call(cbind, draws)), 1L, sd)
  expect_identical(.Random.seed, seed)
  expect_equal(fit$null, H, tolerance = 1e-12)
  expect_equal(fit$P_BM, (1 + sum(H >= 0.5)) / 20)
  expect_equal(fit$exceedance_count, sum(H >= 0.5))
  expect_equal(fit$clade_mean, vapply(draws, mean, numeric(1)), tolerance = 1e-12)
  expect_equal(fit$clade_lower, vapply(draws, quantile, numeric(1), probs = .025, names = FALSE))
  expect_equal(fit$clade_upper, vapply(draws, quantile, numeric(1), probs = .975, names = FALSE))
  set.seed(888)
  compact <- heterogeneity_calibrate_cpp(list(core, core), 0.5, 19L, 7L, FALSE)
  expect_null(compact$null)
  fit$null <- NULL
  compact$null <- NULL
  expect_equal(compact, fit, tolerance = 1e-12)
})
