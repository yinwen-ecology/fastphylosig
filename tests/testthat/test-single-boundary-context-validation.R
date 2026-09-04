.single_boundary_fixture <- function(n = 24L) {
  set.seed(20260904)
  tree <- ape::rtree(n)
  tree$node.label <- paste0("node", seq_len(tree$Nnode))
  trait <- stats::setNames(stats::rnorm(n), tree$tip.label)
  list(tree = tree, trait = trait)
}

.run_prepared_k <- function(ctx, trait) {
  fast_k(
    ctx,
    trait,
    test = FALSE,
    verbose = FALSE,
    progress = FALSE
  )
}

test_that("one prepared public call performs one complete fingerprint", {
  fixture <- .single_boundary_fixture()
  ctx <- prepare_tree(fixture$tree)
  original <- fastphylosig:::.tree_fingerprint
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  testthat::local_mocked_bindings(
    .tree_fingerprint = function(tree) {
      calls$n <- calls$n + 1L
      original(tree)
    },
    .package = "fastphylosig"
  )

  invisible(.run_prepared_k(ctx, fixture$trait))
  expect_identical(calls$n, 1L)
  expect_null(attr(ctx, "fastphylosig_validation_capability", exact = TRUE))

  invisible(.run_prepared_k(ctx, fixture$trait))
  expect_identical(calls$n, 2L)
  expect_null(attr(ctx, "fastphylosig_validation_capability", exact = TRUE))
})

test_that("each protected structural mutation is rejected on the next call", {
  fixture <- .single_boundary_fixture()
  mutations <- list(
    edge = function(ctx) {
      ctx$tree$edge[1L, 2L] <- ctx$tree$edge[2L, 2L]
      ctx
    },
    edge_length = function(ctx) {
      ctx$tree$edge.length[[1L]] <- ctx$tree$edge.length[[1L]] + 0.01
      ctx
    },
    tip_label = function(ctx) {
      ctx$tree$tip.label[[1L]] <- paste0(ctx$tree$tip.label[[1L]], "_changed")
      ctx
    },
    Nnode = function(ctx) {
      ctx$tree$Nnode <- as.integer(ctx$tree$Nnode + 1L)
      ctx
    }
  )

  for (name in names(mutations)) {
    ctx <- prepare_tree(fixture$tree)
    invisible(.run_prepared_k(ctx, fixture$trait))
    ctx <- mutations[[name]](ctx)
    expect_error(
      .run_prepared_k(ctx, fixture$trait),
      "prepared tree was modified after caching",
      info = name
    )
  }
})

test_that("non-computational node labels remain outside the fingerprint contract", {
  fixture <- .single_boundary_fixture()
  ctx <- prepare_tree(fixture$tree)
  before <- .run_prepared_k(ctx, fixture$trait)
  ctx$tree$node.label[[1L]] <- "metadata-only change"
  after <- .run_prepared_k(ctx, fixture$trait)

  expect_identical(unname(as.numeric(after)), unname(as.numeric(before)))
})
