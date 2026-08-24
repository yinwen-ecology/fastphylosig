test_that("streaming D kernel agrees on observed sums and keeps bounded output", {
  skip_if_not(exists("phylo_d_stream_cpp", envir = asNamespace("fastphylosig"),
                     inherits = FALSE))
  set.seed(7101)
  tree <- ape::reorder.phylo(ape::rtree(12), "pruningwise")
  edge <- matrix(as.integer(tree$edge), ncol = 2L)
  states <- sample(0:1, 12, replace = TRUE)
  if (length(unique(states)) < 2L) states[1:2] <- 0:1

  expected <- fastphylosig:::phylo_d_sums_cpp(
    matrix(states, ncol = 1L), edge = edge,
    edge_length = tree$edge.length, n_tip = 12L
  )[[1L]]
  streamed <- fastphylosig:::phylo_d_stream_cpp(
    observed = states, edge = edge, edge_length = tree$edge.length,
    n_tip = 12L, nsim = 64L, prop_state1 = mean(states == 0),
    chunk_size = 32L, return_sim = FALSE
  )

  expect_equal(streamed$observed, expected, tolerance = 1e-12)
  expect_false(any(c("random", "brownian") %in% names(streamed)))
  expect_equal(streamed$p_random, streamed$p_random[[1L]])
  expect_equal(streamed$p_brownian, streamed$p_brownian[[1L]])
})

test_that("streaming D retains only O(nsim) contrast sums when requested", {
  skip_if_not(exists("phylo_d_stream_cpp", envir = asNamespace("fastphylosig"),
                     inherits = FALSE))
  set.seed(7102)
  tree <- ape::reorder.phylo(ape::rtree(15), "pruningwise")
  edge <- matrix(as.integer(tree$edge), ncol = 2L)
  states <- rep(0:1, length.out = 15L)
  streamed <- fastphylosig:::phylo_d_stream_cpp(
    observed = states, edge = edge, edge_length = tree$edge.length,
    n_tip = 15L, nsim = 96L, prop_state1 = mean(states == 0),
    chunk_size = 64L, return_sim = TRUE
  )

  expect_length(streamed$random, 96L)
  expect_length(streamed$brownian, 96L)
  expect_equal(streamed$mean_random, mean(streamed$random), tolerance = 1e-12)
  expect_equal(streamed$mean_brownian, mean(streamed$brownian),
               tolerance = 1e-12)
  expect_equal(streamed$p_random, mean(streamed$random < streamed$observed),
               tolerance = 1e-12)
  expect_equal(streamed$p_brownian,
               mean(streamed$brownian > streamed$observed), tolerance = 1e-12)
})

test_that("streaming and controlled D null components have the same scale", {
  skip_if_not_installed("caper")
  skip_if_not(exists("phylo_d_stream_cpp", envir = asNamespace("fastphylosig"),
                     inherits = FALSE))
  set.seed(7103)
  tree <- ape::rtree(18)
  x <- sample(0:1, 18L, replace = TRUE)
  if (length(unique(x)) < 2L) x[1:2] <- 0:1
  names(x) <- tree$tip.label
  nsim <- 256L
  ds <- x[tree$tip.label]
  random_states <- replicate(nsim, sample(ds))
  C <- unclass(caper::VCV.array(tree))
  brownian_raw <- t(chol(C)) %*% matrix(stats::rnorm(18L * nsim), nrow = 18L)
  brownian_states <- fastphylosig:::brownian_threshold_cpp(
    brownian_raw, mean(ds == 0)
  )

  controlled <- fastphylosig::fast_d(
    tree, x, nsim = nsim, random_states = random_states,
    brownian_states = brownian_states, return_sim = TRUE, verbose = FALSE
  )
  set.seed(7104)
  streamed <- fastphylosig::fast_d(
    tree, x, nsim = nsim, return_sim = TRUE, chunk_size = 64L,
    verbose = FALSE
  )

  expect_equal(streamed$Parameters$Observed,
               controlled$Parameters$Observed, tolerance = 1e-12)
  expect_lt(abs(streamed$Parameters$MeanRandom -
                  controlled$Parameters$MeanRandom), 0.35)
  expect_lt(abs(streamed$Parameters$MeanBrownian -
                  controlled$Parameters$MeanBrownian), 0.35)
})

test_that("D streaming preserves branch-length validation", {
  set.seed(7105)
  tree <- ape::rtree(10)
  x <- rep(0:1, length.out = 10L)
  names(x) <- tree$tip.label

  short <- tree
  short$edge.length[1L] <- 1e-12
  fit <- fastphylosig::fast_d(
    short, x, nsim = 32L, return_sim = FALSE, verbose = FALSE
  )
  expect_true(is.finite(fit$DEstimate))

  zero <- tree
  zero$edge.length[1L] <- 0
  expect_error(
    fastphylosig::fast_d(zero, x, nsim = 32L,
                         return_sim = FALSE, verbose = FALSE),
    "zero"
  )
})
