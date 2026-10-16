# Frozen pre-optimization traversal: test oracle only, never a production path.
.delta_opt_ref_compare <- function(left, right) {
  left <- as.raw(left)
  right <- as.raw(right)
  n <- min(length(left), length(right))
  if (n) {
    for (i in seq_len(n)) {
      a <- as.integer(left[[i]])
      b <- as.integer(right[[i]])
      if (a < b) return(-1L)
      if (a > b) return(1L)
    }
  }
  if (length(left) < length(right)) -1L else
    if (length(left) > length(right)) 1L else 0L
}

.delta_opt_ref_order <- function(keys) {
  n <- length(keys)
  if (n <= 1L) return(seq_len(n))
  index <- seq_len(n)
  work <- index
  width <- 1L
  while (width < n) {
    starts <- seq.int(1L, n, by = 2L * width)
    for (left in starts) {
      middle <- min(left + width - 1L, n)
      right <- min(left + 2L * width - 1L, n)
      if (middle >= right) next
      i <- left
      j <- middle + 1L
      k <- left
      while (i <= middle || j <= right) {
        take_left <- j > right ||
          (i <= middle && .delta_opt_ref_compare(
            keys[[index[[i]]]], keys[[index[[j]]]]
          ) <= 0L)
        if (take_left) {
          work[[k]] <- index[[i]]
          i <- i + 1L
        } else {
          work[[k]] <- index[[j]]
          j <- j + 1L
        }
        k <- k + 1L
      }
    }
    index <- work
    width <- width * 2L
  }
  index
}

.delta_opt_ref_graph <- function(tree) {
  fields <- .v2_source_fields(tree)
  edge <- fields$edge
  n_tip <- fields$n_tip
  n_node <- fields$n_node
  ids <- unique(as.numeric(edge))
  internal_ids <- ids[ids > n_tip]
  if (length(internal_ids) != n_node) {
    .v2_abort("edge endpoints do not define exactly Nnode internal IDs")
  }
  # Dense indices are temporary lookup keys only.  They never enter C2.
  node_ids <- c(seq_len(n_tip), internal_ids)
  index_by_id <- stats::setNames(seq_along(node_ids), as.character(node_ids))
  parent_idx <- unname(index_by_id[as.character(edge[, 1L])])
  child_idx <- unname(index_by_id[as.character(edge[, 2L])])
  if (anyNA(parent_idx) || anyNA(child_idx)) {
    .v2_abort("edge endpoints include an unrecognised node ID")
  }
  parent_idx <- as.integer(parent_idx)
  child_idx <- as.integer(child_idx)
  if (any(parent_idx == child_idx)) .v2_abort("tree contains a self-loop")

  total <- length(node_ids)
  degrees <- tabulate(parent_idx, nbins = total)
  children <- lapply(degrees, function(k) integer(k))
  child_cursor <- integer(total)
  for (i in seq_len(nrow(edge))) {
    parent <- parent_idx[[i]]
    child_cursor[[parent]] <- child_cursor[[parent]] + 1L
    children[[parent]][[child_cursor[[parent]]]] <- child_idx[[i]]
  }
  indegree <- tabulate(child_idx, nbins = total)
  tip_idx <- seq_len(n_tip)
  internal_idx <- n_tip + seq_len(n_node)
  if (any(degrees[tip_idx] != 0L)) {
    .v2_abort("tip nodes cannot be parents")
  }
  if (any(indegree[tip_idx] != 1L)) {
    .v2_abort("each tip must have exactly one parent")
  }
  if (any(indegree[internal_idx] > 1L)) {
    .v2_abort("internal nodes must not have multiple parents")
  }
  roots <- internal_idx[indegree[internal_idx] == 0L]
  if (length(roots) != 1L) .v2_abort("tree must have exactly one structural root")
  root_idx <- roots[[1L]]
  # Empty internal nodes are malformed; unary nodes remain valid here and are
  # left for the existing method-specific readiness checks.
  if (any(degrees[internal_idx] < 1L)) {
    .v2_abort("internal nodes must have at least one child")
  }

  # One iterative DFS performs reachability, cycle detection, and an input
  # adjacency postorder.  The source row order only affects this temporary
  # walk; canonical child ordering below is based on byte-derived ranks.
  colour <- integer(total)
  postorder <- integer(total)
  post_count <- 0L
  stack_node <- integer(max(4L, 2L * total + 4L))
  stack_exit <- logical(length(stack_node))
  top <- 1L
  stack_node[[top]] <- root_idx
  stack_peak <- top
  while (top > 0L) {
    node <- stack_node[[top]]
    is_exit <- stack_exit[[top]]
    top <- top - 1L
    if (is_exit) {
      colour[[node]] <- 2L
      post_count <- post_count + 1L
      postorder[[post_count]] <- node
      next
    }
    if (colour[[node]] == 1L) .v2_abort("tree contains a directed cycle")
    if (colour[[node]] == 2L) next
    colour[[node]] <- 1L
    top <- top + 1L
    stack_node[[top]] <- node
    stack_exit[[top]] <- TRUE
    kids <- children[[node]]
    if (length(kids)) {
      for (kid in rev(kids)) {
        if (colour[[kid]] == 1L) {
          .v2_abort("tree contains a directed cycle")
        }
        if (colour[[kid]] != 2L) {
          top <- top + 1L
          stack_node[[top]] <- kid
          stack_exit[[top]] <- FALSE
        }
      }
    }
    stack_peak <- max(stack_peak, top)
  }
  if (post_count != total || any(colour != 2L)) {
    .v2_abort("tree is disconnected from its structural root")
  }
  postorder <- postorder[seq_len(post_count)]

  labels <- fields$labels
  tip_order <- .delta_opt_ref_order(labels$bytes)
  tip_rank <- integer(n_tip)
  tip_rank[tip_order] <- seq_len(n_tip)
  min_tip_rank <- rep(NA_integer_, total)
  min_tip_rank[tip_idx] <- tip_rank
  ordered_children <- vector("list", total)
  for (node in postorder) {
    if (node <= n_tip) next
    kids <- children[[node]]
    child_mins <- min_tip_rank[kids]
    if (anyNA(child_mins) || anyDuplicated(child_mins)) {
      .v2_abort("child subtrees do not have distinct canonical minimum tip ranks")
    }
    child_order <- order(child_mins, method = "radix")
    ordered_children[[node]] <- kids[child_order]
    min_tip_rank[[node]] <- min(child_mins)
  }

  # Canonical pre-order assigns the root n_tip+1 and then assigns internal IDs
  # on first entry.  Tips retain their public IDs exactly.
  canonical_preorder <- integer(total)
  preorder_count <- 0L
  preorder_stack <- integer(max(2L, total + 1L))
  top <- 1L
  preorder_stack[[top]] <- root_idx
  while (top > 0L) {
    node <- preorder_stack[[top]]
    top <- top - 1L
    preorder_count <- preorder_count + 1L
    canonical_preorder[[preorder_count]] <- node
    kids <- ordered_children[[node]]
    if (length(kids)) {
      for (kid in rev(kids)) {
        top <- top + 1L
        preorder_stack[[top]] <- kid
      }
    }
  }
  canonical_preorder <- canonical_preorder[seq_len(preorder_count)]
  if (length(canonical_preorder) != total) {
    .v2_abort("canonical pre-order did not cover every node")
  }
  canonical_id <- integer(total)
  canonical_id[tip_idx] <- tip_idx
  internal_preorder <- canonical_preorder[canonical_preorder > n_tip]
  canonical_id[internal_preorder] <- seq.int(
    n_tip + 1L, length.out = length(internal_preorder)
  )

  # Canonical postorder emits child edges after each child subtree.  The
  # child-specific source row is retained for exact branch association.
  canonical_postorder <- integer(total)
  post_count <- 0L
  post_stack_node <- integer(max(4L, 2L * total + 4L))
  post_stack_exit <- logical(length(post_stack_node))
  top <- 1L
  post_stack_node[[top]] <- root_idx
  while (top > 0L) {
    node <- post_stack_node[[top]]
    is_exit <- post_stack_exit[[top]]
    top <- top - 1L
    if (is_exit) {
      post_count <- post_count + 1L
      canonical_postorder[[post_count]] <- node
      next
    }
    top <- top + 1L
    post_stack_node[[top]] <- node
    post_stack_exit[[top]] <- TRUE
    kids <- ordered_children[[node]]
    if (length(kids)) {
      for (kid in rev(kids)) {
        top <- top + 1L
        post_stack_node[[top]] <- kid
        post_stack_exit[[top]] <- FALSE
      }
    }
  }
  canonical_postorder <- canonical_postorder[seq_len(post_count)]
  if (length(canonical_postorder) != total) {
    .v2_abort("canonical post-order did not cover every node")
  }

  edge_len_by_child <- rep(NA_real_, total)
  edge_row_for_child <- integer(total)
  parent_for_child <- rep(NA_integer_, total)
  edge_len_by_child[child_idx] <- fields$edge_length
  edge_row_for_child[child_idx] <- seq_len(nrow(edge))
  parent_for_child[child_idx] <- parent_idx
  nonroot <- canonical_postorder[canonical_postorder != root_idx]
  if (anyNA(edge_len_by_child[nonroot]) || anyNA(parent_for_child[nonroot]) ||
      any(edge_row_for_child[nonroot] < 1L)) {
    .v2_abort("canonical edge traversal lost a source branch association")
  }

  # Descriptor IDs are independent from source dense indices.  Tip descriptors
  # use byte-derived ranks; internal descriptor IDs follow canonical IDs.
  descriptor_id <- canonical_id
  descriptor_id[tip_idx] <- tip_rank
  descriptor <- vector("list", total)
  for (node in canonical_postorder) {
    did <- descriptor_id[[node]]
    if (node <= n_tip) {
      descriptor[[did]] <- list(
        kind = "tip",
        label_rank = as.integer(tip_rank[[node]]),
        child_descriptor_ids = integer()
      )
    } else {
      kids <- ordered_children[[node]]
      descriptor[[did]] <- list(
        kind = "internal",
        child_count = as.integer(length(kids)),
        child_min_tip_ranks = as.integer(min_tip_rank[kids]),
        child_descriptor_ids = as.integer(descriptor_id[kids])
      )
    }
  }

  list(
    fields = fields,
    node_ids = node_ids,
    parent_idx = parent_idx,
    child_idx = child_idx,
    children = children,
    ordered_children = ordered_children,
    root_idx = root_idx,
    tip_idx = tip_idx,
    internal_idx = internal_idx,
    tip_order = tip_order,
    tip_rank = tip_rank,
    min_tip_rank = min_tip_rank,
    canonical_id = canonical_id,
    parent_for_child = parent_for_child,
    descriptor_id = descriptor_id,
    descriptor = descriptor,
    postorder = postorder,
    canonical_preorder = canonical_preorder,
    canonical_postorder = canonical_postorder,
    edge_len_by_child = edge_len_by_child,
    edge_row_for_child = edge_row_for_child,
    metrics = list(
      schema_version = .v2_canonical_contract_version,
      full_edge_scans = 1L,
      descendant_label_vectors_materialized = 0L,
      child_tuple_entries = as.integer(sum(lengths(
        ordered_children[internal_idx]
      ))),
      descriptor_entries = as.integer(length(descriptor)),
      child_sort_calls = as.integer(sum(vapply(
        ordered_children[internal_idx], length, integer(1L)
      ) >= 0L)),
      explicit_stack_peak = as.integer(stack_peak)
    )
  )
}

