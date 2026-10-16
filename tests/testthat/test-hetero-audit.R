test_that("ineligible K and lambda crowns do not receive local signal P values", {
  z <- .hetero_fixture()
  x <- z$X
  x[1:25, 1L] <- NA_real_
  for (method in c("K", "lambda")) {
    h <- fast_hetero(z$tree, x, method = method, nodes = z$nodes,
      min_tips = 20L, clade_test = TRUE, clade_nsim = 5L, progress = FALSE)
    expect_identical(h$clades$status, c("insufficient_data", "ok"))
    expect_true(is.na(h$clades$P_signal[1L]))
    expect_identical(h$clades$signal_nsim_successful[1L], 0L)
    expect_true(is.finite(h$clades$P_signal[2L]))
  }
})

test_that("lambda local LRT reuses the scan fit with specialist probabilities", {
  set.seed(441); tree <- ape::rtree(180)
  X <- matrix(rnorm(180*3), 180, 3, dimnames = list(tree$tip.label, c("a", "b", "c")))
  X[seq(5, 180, 9), 2] <- NA_real_
  X[,3] <- 1
  h <- fast_hetero(tree, X, method = "lambda", min_tips = 15, max_tips = 60,
    clade_test = TRUE, progress = FALSE)
  for (i in seq_along(h$selection$nodes)) {
    ids <- (i-1L)*3L+seq_len(3L)
    sub <- ape::extract.clade(tree, h$selection$nodes[i])
    ref <- fast_lambda(sub, X, test = TRUE, lambda_profile = FALSE,
      verbose = FALSE, progress = FALSE)
    eligible <- h$clades$status[ids] == "ok"
    expect_equal(h$clades$P_signal[ids[eligible]], ref$P_fast[eligible], tolerance = 1e-12)
    expect_true(all(is.na(h$clades$P_signal[ids[!eligible]])))
  }
})

test_that("cached Delta masks and local-only permutations preserve estimates and RNG", {
  set.seed(516); tree <- ape::rtree(150)
  X <- data.frame(a = sample(letters[1:4], 150, TRUE),
    b = sample(letters[1:3], 150, TRUE), c = sample(letters[1:4], 150, TRUE),
    row.names = tree$tip.label)
  X[c(3, 51, 78), 2L] <- NA
  control <- .heterogeneity_delta_control(list(model = "ER", mcmc_sim = 1200L,
    thin = 5L, burn = 100L))
  set.seed(628)
  h <- suppressWarnings(fast_hetero(tree, X, method = "Delta", min_tips = 15,
    max_tips = 45, clade_test = TRUE, clade_nsim = 3L,
    delta_control = control, progress = FALSE))
  seed <- .Random.seed
  set.seed(628)
  observed <- suppressWarnings(fast_hetero(tree, X, method = "Delta", min_tips = 15,
    max_tips = 45, delta_control = control, progress = FALSE))
  expected <- rep(NA_real_, nrow(h$clades))
  for (j in seq_len(ncol(X))) for (i in seq_along(h$selection$nodes)) {
    row <- (i-1L)*ncol(X)+j
    if (observed$clades$status[row] != "ok") next
    crown <- ape::extract.clade(tree, h$selection$nodes[i])
    keep <- crown$tip.label[!is.na(X[crown$tip.label,j])]
    ctx <- .heterogeneity_delta_context(crown, keep, tree$tip.label, control)
    draws <- replicate(3L, {
      full <- X[,j]
      full[ctx$rows] <- sample(full[ctx$rows], length(ctx$rows))
      .heterogeneity_delta_evaluate(ctx, full)$estimate
    })
    if (all(is.finite(draws))) expected[row] <- mean(draws >= observed$clades$estimate[row])
  }
  expect_equal(h$clades$estimate, observed$clades$estimate)
  expect_identical(h$clades$P_signal, expected)
  expect_identical(.Random.seed, seed)
  # Mutable ACE contexts belong to this call only and never enter saved results.
  expect_false("delta_contexts" %in% names(h))
})

test_that("all plus-one reference tails share complete-draw MCSE and summaries", {
  draws <- c(.2, .6, .8, 1.3, 1.5)
  q <- 3/5
  tail <- heterogeneity_bootstrap_summary_cpp(draws, .7, FALSE)
  expect_identical(tail$exceedance_count, 3L)
  expect_equal(tail$P_heterogeneity, 4/6)
  expect_equal(tail$MCSE_P_heterogeneity, sqrt(5*q*(1-q))/6)
  expect_equal(tail$null_mean, mean(draws))
  expect_equal(tail$null_upper, unname(quantile(draws, .95)))
  failed <- heterogeneity_bootstrap_summary_cpp(c(1, Inf, NA, 3), 2, TRUE)
  expect_identical(failed$nsim_successful, 2L)
  expect_true(is.na(failed$P_heterogeneity))
  expect_true(is.na(failed$exceedance_count))
  expect_identical(failed$null, c(1, Inf, NA_real_, 3))
})

test_that("Mk caches repeated lengths exactly and accommodates large trees", {
  tree <- ape::stree(32768L, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 32768L)
  # The former per-node cube was about 72 MiB and exceeded the 64 MiB cap.
  # All these transitions now share a single 12 x 12 matrix, without rounding.
  index <- matrix(1L, 12L, 12L); diag(index) <- 0L
  w <- heterogeneity_mk_prepare_cpp(core, index, .2)
  set.seed(534); x <- heterogeneity_mk_simulate_cpp(w, 2L); seed <- .Random.seed
  expect_identical(dim(x), c(32768L, 2L))
  expect_true(all(x >= 1L & x <= 12L))
  set.seed(534)
  expect_identical(cbind(heterogeneity_mk_simulate_cpp(w, 1L),
    heterogeneity_mk_simulate_cpp(w, 1L)), x)
  expect_identical(.Random.seed, seed)
  small <- ape::stree(16L, type = "balanced")
  small$edge.length <- rep(c(.1, .2, .3), length.out = nrow(small$edge))
  core <- compile_tree_cpp(small$edge, small$edge.length, 16L)
  rate <- .fast_ace_rate_index("ER", 4L)
  w <- heterogeneity_mk_prepare_cpp(core, rate$rate_index, .3)
  set.seed(13); native <- heterogeneity_mk_simulate_cpp(w, 11L); seed <- .Random.seed
  set.seed(13)
  expect_identical(native, .hetero_er_states(core, 4L, .3, 11L))
  expect_identical(.Random.seed, seed)
})
