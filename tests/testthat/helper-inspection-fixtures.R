# Fixtures for direct current tree-inspection contract tests.
.inspection_clone <- function(tree) {
  unserialize(serialize(tree, connection = NULL, version = 2L))
}

.inspection_make_pectinate <- function(n) {
  stopifnot(length(n) == 1L, n >= 2L)
  n <- as.integer(n)
  edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  if (n > 2L) {
    for (i in seq_len(n - 2L)) {
      parent <- n + i
      edge[2L * i - 1L, ] <- c(parent, i)
      edge[2L * i, ] <- c(parent, parent + 1L)
    }
  }
  root <- 2L * n - 1L
  edge[2L * n - 3L, ] <- c(root, n - 1L)
  edge[2L * n - 2L, ] <- c(root, n)
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

.inspection_shuffle_edges <- function(tree) {
  out <- .inspection_clone(tree)
  rows <- rev(seq_len(nrow(out$edge)))
  out$edge <- out$edge[rows, , drop = FALSE]
  if (!is.null(out$edge.length)) out$edge.length <- out$edge.length[rows]
  out
}

.inspection_renumber_internals <- function(tree, permutation) {
  out <- .inspection_clone(tree)
  n_tip <- length(out$tip.label)
  old <- seq.int(n_tip + 1L, n_tip + out$Nnode)
  stopifnot(length(permutation) == length(old),
            setequal(permutation, old))
  mapping <- stats::setNames(as.integer(permutation), as.character(old))
  internal <- out$edge > n_tip
  out$edge[internal] <- unname(mapping[as.character(out$edge[internal])])
  if (!is.null(out$node.label)) {
    labels <- character(length(old))
    labels[as.integer(mapping) - n_tip] <- out$node.label
    out$node.label <- labels
  }
  out
}

.inspection_unary_tree <- function() {
  structure(
    list(
      edge = matrix(c(3L, 4L, 4L, 1L, 3L, 2L), ncol = 2L, byrow = TRUE),
      edge.length = c(1, 1, 1),
      tip.label = c("a", "b"),
      Nnode = 2L
    ),
    class = "phylo"
  )
}

.inspection_disconnected_tree <- function() {
  structure(
    list(
      edge = matrix(c(5L, 1L, 5L, 2L, 6L, 3L, 6L, 4L),
                    ncol = 2L, byrow = TRUE),
      edge.length = rep(1, 4L),
      tip.label = paste0("tip", 1:4),
      Nnode = 2L
    ),
    class = "phylo"
  )
}

.inspection_fixture_set <- function() {
  random <- ape::rtree(14L)
  random$node.label <- paste0("node", seq_len(random$Nnode))
  balanced <- ape::compute.brlen(ape::stree(16L, type = "balanced"))
  pectinate <- .inspection_make_pectinate(40L)
  polytomy <- ape::compute.brlen(ape::stree(12L, type = "star"))
  shuffled <- .inspection_shuffle_edges(random)
  internal <- seq.int(length(random$tip.label) + 1L,
                      length(random$tip.label) + random$Nnode)
  renumbered <- .inspection_renumber_internals(random, rev(internal))
  root <- setdiff(unique(random$edge[, 1L]), unique(random$edge[, 2L]))
  other <- setdiff(internal, root)[[1L]]
  root_swap <- internal
  root_swap[match(root, internal)] <- other
  root_swap[match(other, internal)] <- root
  different_root_numbering <- .inspection_renumber_internals(random, root_swap)
  two_tip <- ape::read.tree(text = "(a:1,b:2);")
  large <- .inspection_make_pectinate(1000L)

  cycle <- .inspection_clone(random)
  cycle$edge[1L, 2L] <- cycle$edge[1L, 1L]
  disconnected <- .inspection_disconnected_tree()
  unary <- .inspection_unary_tree()
  malformed_edge <- .inspection_clone(random)
  malformed_edge$edge <- matrix(seq_len(9L), nrow = 3L, ncol = 3L)
  invalid_nnode <- .inspection_clone(random)
  invalid_nnode$Nnode <- NA_real_
  negative_branch <- .inspection_clone(random)
  negative_branch$edge.length[[1L]] <- -1
  zero_terminal <- .inspection_clone(random)
  terminal <- which(zero_terminal$edge[, 2L] <= length(zero_terminal$tip.label))[[1L]]
  zero_terminal$edge.length[[terminal]] <- 0
  invalid_root <- ape::unroot(.inspection_clone(random))

  list(
    balanced = balanced,
    random = random,
    pectinate = pectinate,
    polytomy = polytomy,
    shuffled_edge_order = shuffled,
    safe_internal_renumbering = renumbered,
    different_root_numbering = different_root_numbering,
    two_tip = two_tip,
    large_pectinate = large,
    cycle = cycle,
    disconnected = disconnected,
    unary = unary,
    malformed_edge = malformed_edge,
    invalid_nnode = invalid_nnode,
    negative_branch = negative_branch,
    zero_terminal_branch = zero_terminal,
    invalid_root_representation = invalid_root
  )
}

