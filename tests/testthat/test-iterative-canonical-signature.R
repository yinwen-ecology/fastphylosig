.legacy_canonical_tree_signature <- function(tree) {
  if (!inherits(tree, "phylo")) return(NULL)
  edge <- tree$edge
  if (!is.matrix(edge) || ncol(edge) != 2L || !nrow(edge)) return(NULL)
  tip <- as.character(tree$tip.label)
  n_tip <- length(tip)
  values <- suppressWarnings(as.numeric(edge))
  if (any(!is.finite(values)) || any(values != floor(values))) return(NULL)
  edge <- matrix(as.integer(values), ncol = 2L)
  ids <- sort(unique(c(edge[, 1L], edge[, 2L])))
  internal <- ids[ids > n_tip]
  root_candidates <- setdiff(unique(edge[, 1L]), unique(edge[, 2L]))
  root <- if (length(root_candidates) == 1L) root_candidates[[1L]] else NA_integer_
  descendants <- function(node, active = integer()) {
    if (node <= n_tip) return(tip[[node]])
    if (node %in% active) return(paste0("!cycle:", node))
    kids <- edge[edge[, 1L] == node, 2L]
    if (!length(kids)) return(paste0("!empty:", node))
    sort(unlist(lapply(kids, descendants, active = c(active, node)),
                use.names = FALSE))
  }
  root_key <- if (is.finite(root)) {
    paste(descendants(root), collapse = "\r")
  } else {
    NA_character_
  }
  outdegree <- tabulate(edge[, 1L], nbins = max(c(n_tip, ids)))
  polytomy <- if (length(internal)) sum(outdegree[internal] > 2L) else 0L
  list(
    tip_label = tip,
    root = root,
    root_descendants = root_key,
    root_degree = if (is.finite(root)) outdegree[[root]] else NA_integer_,
    polytomies = as.integer(polytomy),
    edge_length = if (!is.null(tree$edge.length) &&
                      length(tree$edge.length) == nrow(edge)) {
      as.numeric(tree$edge.length)
    } else NULL,
    edge = edge
  )
}

.legacy_internal_order <- function(tree) {
  before <- .legacy_canonical_tree_signature(tree)
  edge <- before$edge
  n_tip <- length(before$tip_label)
  ids <- sort(unique(c(edge)))
  internal <- ids[ids > n_tip]
  roots <- setdiff(unique(edge[, 1L]), unique(edge[, 2L]))
  root <- roots[[1L]]
  parent <- edge[, 1L]
  child <- edge[, 2L]
  children_of <- function(node) child[parent == node]
  memo <- new.env(parent = emptyenv())
  descendants <- function(node, active = integer()) {
    key <- as.character(node)
    if (exists(key, memo, inherits = FALSE)) {
      return(get(key, memo, inherits = FALSE))
    }
    if (node <= n_tip) {
      ans <- if (length(before$tip_label) >= node) {
        before$tip_label[[node]]
      } else {
        paste0("#", node)
      }
      assign(key, ans, memo)
      return(ans)
    }
    if (node %in% active) return(paste0("!cycle:", node))
    kids <- children_of(node)
    ans <- if (!length(kids)) {
      paste0("!empty:", node)
    } else {
      sort(unlist(lapply(kids, descendants, active = c(active, node)),
                  use.names = FALSE))
    }
    assign(key, ans, memo)
    ans
  }
  remaining <- setdiff(internal, root)
  if (length(remaining)) {
    keys <- vapply(remaining, function(node) {
      paste(descendants(node), collapse = "\r")
    }, character(1L))
    remaining <- remaining[order(keys, remaining)]
  }
  c(root, remaining)
}

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

test_that("iterative canonical signatures exactly match the recursive definition", {
  set.seed(20260904)
  random <- ape::rtree(17L)
  random$node.label <- paste0("node", seq_len(random$Nnode))
  balanced <- ape::compute.brlen(ape::stree(16L, type = "balanced"))
  pectinate <- .make_pectinate_tree(40L)
  polytomy <- ape::compute.brlen(ape::stree(12L, type = "star"))
  shuffled <- .shuffle_edge_rows(random)
  internals <- seq.int(length(random$tip.label) + 1L,
                       length(random$tip.label) + random$Nnode)
  renumbered <- .renumber_internal_nodes(random, rev(internals))
  root <- setdiff(unique(random$edge[, 1L]), unique(random$edge[, 2L]))
  other <- setdiff(internals, root)[[1L]]
  root_swapped <- .renumber_internal_nodes(
    random,
    replace(replace(internals, root - length(random$tip.label), other),
            other - length(random$tip.label), root)
  )
  fixtures <- list(
    balanced = balanced,
    pectinate = pectinate,
    random = random,
    polytomy = polytomy,
    shuffled_edge_order = shuffled,
    safe_node_renumbering = renumbered,
    different_root_numbering = root_swapped
  )

  for (name in names(fixtures)) {
    tree <- fixtures[[name]]
    expect_identical(
      fastphylosig:::.canonical_tree_signature(tree),
      .legacy_canonical_tree_signature(tree),
      info = name
    )
    normalized <- fastphylosig:::.safe_canonicalize_core(tree)
    mapping <- attr(normalized, "fastphylosig_canonicalization")$mapping
    expect_identical(
      unname(as.integer(mapping$old_internal_order)),
      unname(as.integer(.legacy_internal_order(tree))),
      info = paste(name, "canonical internal order")
    )
  }
})

test_that("canonical signatures do not depend on recursive expression depth", {
  expression_limit <- getOption("expressions")
  for (n in c(1000L, 5000L, 20000L)) {
    tree <- .make_pectinate_tree(n)
    signature <- fastphylosig:::.canonical_tree_signature(tree)
    expect_identical(length(signature$tip_label), n)
    expect_identical(
      signature$root_descendants,
      paste(sort(tree$tip.label), collapse = "\r"),
      info = paste("pectinate", n)
    )
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
