test_that("prepared D factors match independent recursion over heterogeneous branches", {
  set.seed(9821)
  tree <- ape::rtree(32)
  tree$edge.length <- exp(seq(-12, 12, length.out = nrow(tree$edge)))
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 32L)
  ordered <- ape::reorder.phylo(tree, "postorder")
  states <- matrix(rbinom(32*7, 1, .5), 32, 7)
  reference <- apply(states, 2, function(x) .hetero_d_contrast(core, x))
  for (threads in c(1L, 2L)) {
    actual <- phylo_d_sums_cpp(states, ordered$edge, ordered$edge.length, 32L, threads)
    expect_equal(as.numeric(actual), reference, tolerance = 1e-12)
  }
  set.seed(9822)
  streamed <- phylo_d_stream_cpp(states[,1], ordered$edge, ordered$edge.length,
    32L, 33L, mean(states[,1]), 32L, TRUE, 1L)
  expect_equal(streamed$observed, reference[1], tolerance = 1e-12)
  before <- .Random.seed
  # The scaled heterogeneity path shares the prepared traversal while retaining
  # its extended accumulation convention. Checking an independent recursion
  # prevents both callers from agreeing through the same arithmetic defect.
  set.seed(9823)
  local <- heterogeneity_d_calibrate_cpp(list(core), list(states[,1]),
    19L, 0L, FALSE, FALSE)
  expect_equal(local$contrast, reference[1], tolerance = 1e-12)
  set.seed(9822)
  threaded <- phylo_d_stream_cpp(states[,1], ordered$edge, ordered$edge.length,
    32L, 33L, mean(states[,1]), 32L, TRUE, 2L)
  expect_identical(threaded, streamed)
  expect_identical(.Random.seed, before)
})
