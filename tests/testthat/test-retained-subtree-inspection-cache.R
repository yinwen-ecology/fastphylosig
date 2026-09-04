.subtree_cache_fixture <- function(n = 30L) {
  set.seed(20260904)
  tree <- ape::rtree(n)
  trait <- stats::setNames(stats::rnorm(n), tree$tip.label)
  list(tree = tree, trait = trait)
}

.count_internal_calls <- function(binding, code) {
  original <- getFromNamespace(binding, "fastphylosig")
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  replacement <- function(...) {
    calls$n <- calls$n + 1L
    original(...)
  }
  bindings <- stats::setNames(list(replacement), binding)
  testthat::local_mocked_bindings(
    !!!bindings,
    .package = "fastphylosig",
    .env = parent.frame()
  )
  force(code)
  calls$n
}

test_that("prepared single-trait analysis reuses one inspected structural entry", {
  fixture <- .subtree_cache_fixture()
  ctx <- prepare_tree(fixture$tree)

  inspect_calls <- .count_internal_calls(".inspect_tree_core", {
    invisible(fast_k(ctx, fixture$trait, test = FALSE,
                     verbose = FALSE, progress = FALSE))
  })
  expect_identical(inspect_calls, 0L)

  subset_calls <- .count_internal_calls(".prepared_tree_subset", {
    invisible(fast_k(ctx, fixture$trait, test = FALSE,
                     verbose = FALSE, progress = FALSE))
  })
  expect_identical(subset_calls, 1L)
})

test_that("retained entries store and reuse package-owned inspection evidence", {
  fixture <- .subtree_cache_fixture()
  ctx <- prepare_tree(fixture$tree)
  keep <- 2:25
  original <- fastphylosig:::.inspect_tree_core
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  testthat::local_mocked_bindings(
    .inspect_tree_core = function(...) {
      calls$n <- calls$n + 1L
      original(...)
    },
    .package = "fastphylosig"
  )

  first <- fastphylosig:::.prepared_tree_subset(ctx, keep,
                                                need_matrix = FALSE)
  expect_identical(calls$n, 1L)
  expect_s3_class(first$inspection, "fastphylosig_tree_check")
  expect_named(first$method_readiness, c("K", "lambda", "D", "Delta"))
  expect_true(is.list(first$topology_metadata))
  expect_identical(first$retained_parent_indices, as.integer(keep))
  expect_identical(first$parent_fingerprint, ctx$fingerprint)
  expect_identical(
    first$subtree_key,
    fastphylosig:::.tree_mask_key(keep, ctx$n_tip)
  )
  expect_identical(first$subtree_fingerprint,
                   fastphylosig:::.tree_fingerprint(first$tree))

  second <- fastphylosig:::.prepared_tree_subset(ctx, keep,
                                                 need_matrix = FALSE)
  expect_identical(calls$n, 1L)
  expect_identical(second$inspection, first$inspection)
})

test_that("different masks and stale evidence cannot reuse an inspection", {
  fixture <- .subtree_cache_fixture()
  ctx <- prepare_tree(fixture$tree)
  original <- fastphylosig:::.inspect_tree_core
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  testthat::local_mocked_bindings(
    .inspect_tree_core = function(...) {
      calls$n <- calls$n + 1L
      original(...)
    },
    .package = "fastphylosig"
  )

  keep_a <- 1:20
  keep_b <- 2:21
  entry_a <- fastphylosig:::.prepared_tree_subset(ctx, keep_a,
                                                  need_matrix = FALSE)
  entry_b <- fastphylosig:::.prepared_tree_subset(ctx, keep_b,
                                                  need_matrix = FALSE)
  expect_identical(calls$n, 2L)
  expect_false(identical(entry_a$subtree_key, entry_b$subtree_key))
  expect_false(identical(entry_a$tree$tip.label, entry_b$tree$tip.label))

  key_a <- fastphylosig:::.tree_mask_key(keep_a, ctx$n_tip)
  stale <- get(key_a, ctx$structural_cache, inherits = FALSE)
  stale$subtree_key <- "stale-key"
  assign(key_a, stale, ctx$structural_cache)
  rebuilt <- fastphylosig:::.prepared_tree_subset(ctx, keep_a,
                                                  need_matrix = FALSE)
  expect_identical(calls$n, 3L)
  expect_identical(rebuilt$subtree_key, key_a)
  expect_identical(rebuilt$tree$tip.label, fixture$tree$tip.label[keep_a])
})

test_that("parent mutation invalidates before any retained cache reuse", {
  fixture <- .subtree_cache_fixture()
  ctx <- prepare_tree(fixture$tree)
  keep <- 3:24
  invisible(fastphylosig:::.prepared_tree_subset(ctx, keep,
                                                 need_matrix = FALSE))
  ctx$tree$edge.length[[1L]] <- ctx$tree$edge.length[[1L]] + 0.5
  expect_error(
    fastphylosig:::.prepared_tree_subset(ctx, keep, need_matrix = FALSE),
    "prepared tree was modified after caching"
  )
})

test_that("bounded structural caches recreate uncached retained evidence", {
  fixture <- .subtree_cache_fixture()
  ctx <- prepare_tree(fixture$tree, max_cached_subsets = 1L)
  keep <- 1:22
  original <- fastphylosig:::.inspect_tree_core
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  testthat::local_mocked_bindings(
    .inspect_tree_core = function(...) {
      calls$n <- calls$n + 1L
      original(...)
    },
    .package = "fastphylosig"
  )

  first <- fastphylosig:::.prepared_tree_subset(ctx, keep,
                                                need_matrix = FALSE)
  second <- fastphylosig:::.prepared_tree_subset(ctx, keep,
                                                 need_matrix = FALSE)
  expect_identical(calls$n, 2L)
  expect_length(ls(ctx$structural_cache, all.names = TRUE), 1L)
  expect_identical(first$subtree_fingerprint, second$subtree_fingerprint)
})

test_that("one NA mask is inspected once and reused across public calls", {
  fixture <- .subtree_cache_fixture()
  fixture$trait[c(2L, 7L, 11L)] <- NA_real_
  ctx <- prepare_tree(fixture$tree)
  original <- fastphylosig:::.inspect_tree_core
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  testthat::local_mocked_bindings(
    .inspect_tree_core = function(...) {
      calls$n <- calls$n + 1L
      original(...)
    },
    .package = "fastphylosig"
  )

  first <- fast_k(ctx, fixture$trait, test = FALSE,
                  verbose = FALSE, progress = FALSE)
  expect_identical(calls$n, 1L)
  second <- fast_k(ctx, fixture$trait, test = FALSE,
                   verbose = FALSE, progress = FALSE)
  expect_identical(calls$n, 1L)
  expect_identical(unname(as.numeric(second)), unname(as.numeric(first)))
})
