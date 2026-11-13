.renumber_internal_nodes <- function(tree, permutation) {
  n_tip <- length(tree$tip.label)
  old <- seq.int(n_tip + 1L, n_tip + tree$Nnode)
  stopifnot(length(permutation) == length(old))
  map <- stats::setNames(as.integer(permutation), as.character(old))
  out <- tree
  internal <- out$edge > n_tip
  out$edge[internal] <- unname(map[as.character(out$edge[internal])])
  if (!is.null(out$node.label)) {
    labels <- character(length(old))
    labels[as.integer(map) - n_tip] <- out$node.label
    out$node.label <- labels
  }
  out
}

.shuffle_edge_rows <- function(tree) {
  order <- rev(seq_len(nrow(tree$edge)))
  tree$edge <- tree$edge[order, , drop = FALSE]
  if (!is.null(tree$edge.length)) tree$edge.length <- tree$edge.length[order]
  tree
}

.make_pectinate_tree <- function(n) {
  stopifnot(n >= 2L)
  edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  if (n > 2L) {
    for (i in seq_len(n - 2L)) {
      parent <- n + i
      edge[2L * i - 1L, ] <- c(parent, i)
      edge[2L * i, ] <- c(parent, parent + 1L)
    }
  }
  last <- 2L * n - 1L
  edge[2L * n - 3L, ] <- c(last, n - 1L)
  edge[2L * n - 2L, ] <- c(last, n)
  structure(
    list(
      edge = edge,
      edge.length = rep(1, nrow(edge)),
      tip.label = paste0("tip", seq_len(n)),
      Nnode = n - 1L
    ),
    class = "phylo"
  )
}

test_that("V2 canonicalization preserves tree structure under representation changes", {
  set.seed(20260904)
  fixtures <- list(
    random = ape::rtree(17L),
    balanced = ape::compute.brlen(ape::stree(16L, type = "balanced")),
    pectinate = .make_pectinate_tree(40L),
    polytomy = ape::compute.brlen(ape::stree(12L, type = "star"))
  )
  fields <- c("edge", "edge.length", "tip.label", "Nnode", "node.label")
  for (name in names(fixtures)) {
    tree <- fixtures[[name]]
    tree$node.label <- paste0("node", seq_len(tree$Nnode))
    original <- serialize(tree, NULL)
    internals <- seq.int(length(tree$tip.label) + 1L,
                         length(tree$tip.label) + tree$Nnode)
    variants <- list(tree, .shuffle_edge_rows(tree),
                     .renumber_internal_nodes(tree, rev(internals)))
    canonical <- fastphylosig:::.safe_canonicalize_core(tree)
    fingerprint <- fastphylosig:::.v2_canonical_fingerprint_raw(tree)
    expect_identical(canonical$tip.label, tree$tip.label, info = name)
    expect_identical(sort(canonical$edge.length), sort(tree$edge.length), info = name)
    expect_equal(ape::cophenetic.phylo(canonical), ape::cophenetic.phylo(tree),
                 tolerance = 1e-12, info = name)
    expect_identical(serialize(tree, NULL), original, info = name)
    for (variant in variants) {
      normalized <- fastphylosig:::.safe_canonicalize_core(variant)
      expect_identical(normalized[fields], canonical[fields], info = name)
      expect_identical(fastphylosig:::.v2_canonical_fingerprint_raw(variant),
                       fingerprint, info = name)
      expect_identical(fastphylosig:::.safe_canonicalize_core(normalized)[fields],
                       canonical[fields], info = name)
    }
  }
})

test_that("V2 traversal does not depend on recursive expression depth", {
  expression_limit <- getOption("expressions")
  for (n in c(1000L, 5000L, 20000L)) {
    tree <- .make_pectinate_tree(n)
    graph <- fastphylosig:::.v2_build_graph(tree)
    expect_length(graph$canonical_preorder, 2L * n - 1L)
    expect_length(graph$canonical_postorder, 2L * n - 1L)
    expect_setequal(graph$canonical_preorder, seq_len(2L * n - 1L))
    expect_identical(graph$metrics$descendant_label_vectors_materialized, 0L)
    expect_identical(getOption("expressions"), expression_limit)
  }
})

test_that("raw lambda accepts a pectinate tree with 1000 tips", {
  tree <- .make_pectinate_tree(1000L)
  x <- stats::setNames(seq_len(1000L) / 1000, tree$tip.label)
  result <- fast_lambda(
    tree,
    x,
    test = FALSE,
    lambda_profile = FALSE,
    verbose = FALSE,
    progress = FALSE
  )
  expect_true(is.finite(result$lambda))
  expect_true(is.finite(result$logLik))
})
