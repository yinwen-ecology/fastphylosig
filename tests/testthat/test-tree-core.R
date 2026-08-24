test_that("compiled tree records binary topology, traversal, and distances", {
  edge <- matrix(
    c(
      5L, 1L,
      5L, 6L,
      6L, 2L,
      6L, 7L,
      7L, 3L,
      7L, 4L
    ),
    ncol = 2L,
    byrow = TRUE
  )
  edge_length <- c(1, 2, 3, 4, 5, 6)
  compiled <- fastphylosig:::compile_tree_cpp(edge, edge_length, 4L)

  expect_s3_class(compiled, "fastphylosig_compiled_tree")
  expect_equal(compiled$n_tip, 4L)
  expect_equal(compiled$n_node, 3L)
  expect_equal(compiled$n_total, 7L)
  expect_equal(compiled$root, 5L)
  expect_equal(compiled$parent, c(5L, 6L, 7L, 7L, 0L, 5L, 6L))
  expect_equal(compiled$child_ptr, c(1L, 1L, 1L, 1L, 1L, 3L, 5L, 7L))
  expect_equal(compiled$children, c(1L, 6L, 2L, 7L, 3L, 4L))
  expect_equal(
    compiled$branch_length_by_node,
    c(1, 3, 5, 6, 0, 2, 4)
  )
  expect_equal(compiled$preorder, c(5L, 1L, 6L, 2L, 7L, 3L, 4L))
  expect_equal(compiled$postorder, c(1L, 2L, 3L, 4L, 7L, 6L, 5L))
  expect_equal(compiled$root_distance, c(1, 5, 11, 12, 0, 2, 6))
  expect_equal(compiled$tip_order, c(1L, 2L, 3L, 4L))
  expect_equal(compiled$tip_index, c(1L, 2L, 3L, 4L, 0L, 0L, 0L))
  expect_true(is.finite(compiled$topology_bytes))
  expect_gt(compiled$topology_bytes, 0)
})

test_that("compiled tree handles a rooted polytomy", {
  edge <- matrix(
    c(6L, 1L, 6L, 2L, 6L, 3L, 6L, 4L, 6L, 5L),
    ncol = 2L,
    byrow = TRUE
  )
  compiled <- fastphylosig:::compile_tree_cpp(
    edge, c(0.1, 0.2, 0.3, 0.4, 0.5), 5L
  )

  expect_equal(compiled$root, 6L)
  expect_equal(compiled$preorder, c(6L, 1L, 2L, 3L, 4L, 5L))
  expect_equal(compiled$postorder, c(1L, 2L, 3L, 4L, 5L, 6L))
  expect_equal(compiled$child_ptr, c(1L, 1L, 1L, 1L, 1L, 1L, 6L))
  expect_equal(compiled$children, 1:5)
  expect_equal(compiled$root_distance, c(0.1, 0.2, 0.3, 0.4, 0.5, 0))
})

test_that("compiled tree validation rejects malformed rooted topologies", {
  cycle <- matrix(
    c(5L, 1L, 5L, 2L, 3L, 4L, 4L, 3L),
    ncol = 2L,
    byrow = TRUE
  )
  cycle_diag <- fastphylosig:::tree_core_validate_cpp(
    cycle, rep(1, nrow(cycle)), 2L
  )
  expect_false(cycle_diag$valid)
  expect_match(cycle_diag$message, "disconnected|cycle")
  expect_error(
    fastphylosig:::compile_tree_cpp(cycle, rep(1, nrow(cycle)), 2L),
    "disconnected|cycle"
  )

  duplicate_parent <- matrix(
    c(3L, 1L, 4L, 1L, 3L, 2L),
    ncol = 2L,
    byrow = TRUE
  )
  expect_error(
    fastphylosig:::compile_tree_cpp(
      duplicate_parent, rep(1, nrow(duplicate_parent)), 2L
    ),
    "at most one parent|duplicate child"
  )

  non_contiguous <- matrix(
    c(5L, 1L, 5L, 2L, 5L, 4L, 4L, 5L),
    ncol = 2L,
    byrow = TRUE
  )
  non_contiguous_diag <- fastphylosig:::tree_core_validate_cpp(
    non_contiguous, rep(1, nrow(non_contiguous)), 2L
  )
  expect_false(non_contiguous_diag$valid)
  expect_match(non_contiguous_diag$message, "contiguous")
  expect_error(
    fastphylosig:::compile_tree_cpp(
      non_contiguous, rep(1, nrow(non_contiguous)), 2L
    ),
    "contiguous"
  )
})

