test_that("resolve_tree normalizes representation and records provenance", {
  skip_if_not(exists("resolve_tree", asNamespace("fastphylosig")))
  tree <- ape::rtree(8)
  # A harmless edge-order permutation is representation-only.  Keep the
  # matching edge lengths with their rows so no path length changes.
  ord <- rev(seq_len(nrow(tree$edge)))
  scrambled <- tree
  scrambled$edge <- tree$edge[ord, , drop = FALSE]
  scrambled$edge.length <- tree$edge.length[ord]

  resolved <- fastphylosig::resolve_tree(scrambled, signal = "k")
  expect_s3_class(resolved, "phylo")
  expect_s3_class(resolved, "fastphylosig_resolved_tree")
  meta <- attr(resolved, "fastphylosig_resolution")
  expect_named(meta, c("signal", "original_status", "changes",
                       "final_status", "ready"))
  expect_identical(meta$signal, "K")
  expect_true(isTRUE(meta$ready))
  expect_identical(meta$final_status, "READY")
  expect_identical(resolved$tip.label, tree$tip.label)
  expect_equal(length(resolved$edge.length), length(tree$edge.length))
  expect_equal(sum(resolved$edge.length), sum(tree$edge.length))
  expect_equal(sort(resolved$edge.length), sort(tree$edge.length))
  expect_equal(ape::cophenetic.phylo(resolved), ape::cophenetic.phylo(tree))

  again <- fastphylosig::resolve_tree(resolved, signal = "K")
  expect_identical(again$tip.label, resolved$tip.label)
  expect_equal(ape::cophenetic.phylo(again), ape::cophenetic.phylo(resolved))
  expect_length(attr(again, "fastphylosig_resolution")$changes, 0L)
})

test_that("resolve_tree repairs equivalent internal node numbering", {
  skip_if_not(exists("resolve_tree", asNamespace("fastphylosig")))
  tree <- ape::rtree(7)
  n_tip <- ape::Ntip(tree)
  internal <- seq.int(n_tip + 1L, n_tip + tree$Nnode)
  # Reverse only internal IDs.  Parent-child relationships and the edge
  # lengths attached to each row remain unchanged, but the root is no longer
  # the canonical n_tip + 1 node required by D/Delta.
  renumber <- setNames(rev(internal), internal)
  equivalent <- tree
  equivalent$edge <- tree$edge
  internal_rows <- equivalent$edge > n_tip
  equivalent$edge[internal_rows] <- unname(
    renumber[as.character(equivalent$edge[internal_rows])]
  )

  before_dist <- ape::cophenetic.phylo(equivalent)
  before_lengths <- equivalent$edge.length
  resolved <- fastphylosig::resolve_tree(equivalent, signal = "D")
  meta <- attr(resolved, "fastphylosig_resolution")
  expect_true(isTRUE(meta$ready))
  expect_identical(meta$final_status, "READY")
  expect_identical(resolved$tip.label, equivalent$tip.label)
  expect_equal(ape::cophenetic.phylo(resolved), before_dist)
  expect_equal(sum(resolved$edge.length), sum(before_lengths))
  expect_equal(sort(resolved$edge.length), sort(before_lengths))
  expect_true(any(grepl("representation|node|reorder", meta$changes,
                         ignore.case = TRUE)))
})

test_that("canonical root normalization makes Delta ready without topology edits", {
  skip_if_not(exists("resolve_tree", asNamespace("fastphylosig")))
  tree <- ape::rtree(8)
  n_tip <- ape::Ntip(tree)
  internal <- seq.int(n_tip + 1L, n_tip + tree$Nnode)
  renumber <- setNames(rev(internal), internal)
  noncanonical <- tree
  noncanonical$edge[noncanonical$edge > n_tip] <- unname(
    renumber[as.character(noncanonical$edge[noncanonical$edge > n_tip])]
  )
  dist_before <- ape::cophenetic.phylo(noncanonical)

  resolved <- fastphylosig::resolve_tree(noncanonical, signal = "Delta")
  expect_true(isTRUE(attr(resolved, "fastphylosig_resolution")$ready))
  expect_equal(ape::cophenetic.phylo(resolved), dist_before)
  expect_identical(resolved$tip.label, noncanonical$tip.label)
})

test_that("resolve_tree rejects unsafe automatic repairs", {
  skip_if_not(exists("resolve_tree", asNamespace("fastphylosig")))
  tree <- ape::rtree(6)
  tree$edge.length[[1L]] <- 0
  expect_error(
    fastphylosig::resolve_tree(tree, signal = "K"),
    "USER_ACTION_REQUIRED|cannot produce a ready tree"
  )

  polytomy <- ape::stree(6, type = "star")
  polytomy$edge.length <- rep(1, nrow(polytomy$edge))
  expect_error(
    fastphylosig::resolve_tree(polytomy, signal = "Delta"),
    "polytom|USER_ACTION_REQUIRED|cannot produce a ready tree"
  )
})

test_that("resolve_tree validates signal and explicit outgroup", {
  skip_if_not(exists("resolve_tree", asNamespace("fastphylosig")))
  tree <- ape::rtree(5)
  expect_error(fastphylosig::resolve_tree(tree, signal = c("K", "D")),
               "exactly one")
  expect_error(fastphylosig::resolve_tree(tree, signal = "K",
                                           outgroup = "missing-tip"),
               "not present")
  expect_error(fastphylosig::resolve_tree(tree, signal = "K",
                                           outgroup = rep(tree$tip.label[[1L]], 2L)),
               "unique|outgroup")
})
