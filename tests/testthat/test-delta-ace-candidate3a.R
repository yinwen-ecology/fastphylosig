test_that("Delta-only ACE buffer reuse preserves the production contract", {
  skip_if_not_installed("ape")
  set.seed(20260814)
  tree <- ape::reorder.phylo(ape::rtree(16), "postorder")
  tw <- fastphylosig:::.fast_ace_tree_workspace(tree)
  ws <- fastphylosig:::.fast_ace_workspace(tw, "ARD", 3L)
  y <- factor(sample(seq_len(3L), 16L, replace = TRUE), levels = seq_len(3L))
  names(y) <- ws$tip_label
  par <- rep(0.1, ws$rate$np)

  d0 <- fastphylosig:::fast_ace_discrete_deviance_cpp(
    ws$edge, ws$edge_length, as.integer(y), ws$rate$rate_index, par,
    reuse_buffers = FALSE
  )
  d1 <- fastphylosig:::fast_ace_discrete_deviance_cpp(
    ws$edge, ws$edge_length, as.integer(y), ws$rate$rate_index, par,
    reuse_buffers = TRUE
  )
  a0 <- fastphylosig:::fast_ace_discrete_liks_cpp(
    ws$edge, ws$edge_length, as.integer(y), ws$rate$rate_index, par,
    marginal = FALSE, reuse_buffers = FALSE
  )
  a1 <- fastphylosig:::fast_ace_discrete_liks_cpp(
    ws$edge, ws$edge_length, as.integer(y), ws$rate$rate_index, par,
    marginal = FALSE, reuse_buffers = TRUE
  )
  expect_identical(d0, d1)
  expect_identical(a0$deviance, a1$deviance)
  expect_identical(a0$lik.anc, a1$lik.anc)
  expect_false("reuse_buffers" %in% names(formals(fastphylosig::fast_ace)))
})

test_that("Delta private path uses Candidate 3A while public fast_ace keeps SE", {
  skip_if_not_installed("ape")
  set.seed(20260815)
  tree <- ape::reorder.phylo(ape::rtree(12), "postorder")
  y <- factor(sample(seq_len(3L), 12L, replace = TRUE), levels = seq_len(3L))
  names(y) <- tree$tip.label
  ws <- fastphylosig:::.fast_ace_workspace(
    fastphylosig:::.fast_ace_tree_workspace(tree), "ARD", 3L
  )
  set.seed(20260816)
  old <- suppressWarnings(fastphylosig:::.delta_one(
    tree, y, 0.1, 0.1, mcmc_sim = 8L, thin = 2L, burn = 2L,
    entropy_code = 2L, model = "ARD", ace_engine = "fast",
    keep_chains = TRUE, ace_workspace = ws, reuse_buffers = FALSE
  ))
  set.seed(20260816)
  new <- suppressWarnings(fastphylosig:::.delta_one(
    tree, y, 0.1, 0.1, mcmc_sim = 8L, thin = 2L, burn = 2L,
    entropy_code = 2L, model = "ARD", ace_engine = "fast",
    keep_chains = TRUE, ace_workspace = ws, reuse_buffers = TRUE
  ))
  expect_identical(old$delta, new$delta)
  expect_identical(old$alpha_chain, new$alpha_chain)
  expect_identical(old$beta_chain, new$beta_chain)
  expect_identical(old$MCSE_Delta, new$MCSE_Delta)
})

test_that("Candidate 3A fixed-rate gate covers state counts and branch scales", {
  skip_if_not_installed("ape")
  for (k in 2:4) {
    set.seed(20260820 + k)
    tree <- ape::reorder.phylo(ape::rtree(18), "postorder")
    tree$edge.length <- tree$edge.length * if (k == 2L) 0.01 else if (k == 3L) 1 else 10
    ws <- fastphylosig:::.fast_ace_workspace(
      fastphylosig:::.fast_ace_tree_workspace(tree), "ARD", k
    )
    y <- factor(sample(seq_len(k), 18L, replace = TRUE), levels = seq_len(k))
    names(y) <- ws$tip_label
    par <- rep(0.1, ws$rate$np)
    old <- fastphylosig:::fast_ace_discrete_liks_cpp(
      ws$edge, ws$edge_length, as.integer(y), ws$rate$rate_index, par,
      marginal = FALSE, reuse_buffers = FALSE
    )
    new <- fastphylosig:::fast_ace_discrete_liks_cpp(
      ws$edge, ws$edge_length, as.integer(y), ws$rate$rate_index, par,
      marginal = FALSE, reuse_buffers = TRUE
    )
    expect_identical(old$deviance, new$deviance)
    expect_identical(old$lik.anc, new$lik.anc)
  }
})
