# fastphylosig 0.2.0 Canonicalization Contract V2-A
#
# Audit-only prototype.  This file is deliberately independent of the package
# namespace and is not sourced by production code.  It defines a compact,
# exact representation for a rooted weighted labelled tree.  The representation
# uses UTF-8 bytes for label ordering, shared descriptor references rather than
# copied descendant-label vectors, and an exact protected snapshot for context
# validation.

.v2_fail <- function(message) {
  stop(paste0("V2 prototype: ", message), call. = FALSE)
}

.v2_concat <- function(parts) {
  if (!length(parts)) return(raw(0L))
  unlist(parts, recursive = TRUE, use.names = FALSE)
}

.v2_ascii_raw <- function(value) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value)) {
    .v2_fail("internal ASCII token is not scalar and non-missing")
  }
  charToRaw(value)
}

.v2_hex <- function(bytes) {
  if (!length(bytes)) return("")
  paste(sprintf("%02X", as.integer(bytes)), collapse = "")
}

.v2_label_bytes <- function(labels) {
  if (!is.character(labels)) {
    .v2_fail("tip.label must be a character vector")
  }
  if (!length(labels) || anyNA(labels) || any(!nzchar(labels))) {
    .v2_fail("tip.label must contain non-empty, non-missing labels")
  }
  stored_labels <- labels
  labels_utf8 <- tryCatch(enc2utf8(stored_labels), error = function(e) NULL)
  if (is.null(labels_utf8) || anyNA(labels_utf8)) {
    .v2_fail("tip.label contains a value that cannot be represented as UTF-8")
  }
  bytes <- lapply(labels_utf8, function(label) {
    tryCatch(charToRaw(label), error = function(e) NULL)
  })
  if (any(vapply(bytes, is.null, logical(1L))) ||
      any(!vapply(bytes, length, integer(1L))) ) {
    .v2_fail("tip.label contains an invalid or empty byte sequence")
  }
  # Hex is used only to detect exact duplicate byte sequences.  It is not an
  # ordering or probabilistic identity mechanism.
  byte_ids <- vapply(bytes, .v2_hex, character(1L))
  if (anyDuplicated(byte_ids)) {
    .v2_fail("tip.label byte sequences must be unique")
  }
  # `stored_labels` is the user value (including its encoding marker).  The
  # UTF-8 conversion is a comparison key only and is never written back.
  list(labels = stored_labels, utf8 = labels_utf8,
       bytes = bytes, byte_ids = byte_ids)
}

.v2_bytes_compare <- function(left, right) {
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

.v2_order_bytes <- function(keys) {
  n <- length(keys)
  if (n <= 1L) return(seq_len(n))
  # Bottom-up merge sort is locale-independent and keeps the comparator
  # explicit.  Ties retain input order, although duplicate labels are rejected.
  index <- seq_len(n)
  width <- 1L
  while (width < n) {
    next_index <- index
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
          (i <= middle && .v2_bytes_compare(
            keys[[index[[i]]]], keys[[index[[j]]]]
          ) <= 0L)
        if (take_left) {
          next_index[[k]] <- index[[i]]
          i <- i + 1L
        } else {
          next_index[[k]] <- index[[j]]
          j <- j + 1L
        }
        k <- k + 1L
      }
    }
    index <- next_index
    width <- width * 2L
  }
  index
}

.v2_float_bytes <- function(value) {
  value <- as.numeric(value)
  if (length(value) != 1L || !is.finite(value)) {
    .v2_fail("edge.length must contain only finite scalar values")
  }
  connection <- rawConnection(raw(0L), open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(value, connection, size = 8L, endian = "big")
  rawConnectionValue(connection)
}

.v2_u32_bytes <- function(value) {
  value <- as.integer(value)
  if (anyNA(value) || any(value < 0L)) {
    .v2_fail("unsigned integer encoding received an invalid value")
  }
  connection <- rawConnection(raw(0L), open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(value, connection, size = 4L, endian = "big")
  rawConnectionValue(connection)
}

.v2_count_payload <- function(n) {
  if (length(n) != 1L || is.na(n) || n < 0 || n != floor(n)) {
    .v2_fail("invalid vector length in exact encoding")
  }
  .v2_u32_bytes(n)
}

.v2_field <- function(tag, payload) {
  tag_bytes <- charToRaw(as.character(tag))
  payload <- as.raw(payload)
  # tag length, tag bytes, payload length, payload bytes: every boundary is
  # explicit, including labels containing separators or NUL bytes.
  .v2_concat(list(
    .v2_count_payload(length(tag_bytes)), tag_bytes,
    .v2_count_payload(length(payload)), payload
  ))
}

.v2_label_vector_payload <- function(label_bytes) {
  items <- lapply(label_bytes, function(bytes) {
    .v2_concat(list(.v2_count_payload(length(bytes)), bytes))
  })
  .v2_concat(c(list(.v2_count_payload(length(items))), items))
}

.v2_integer_vector_payload <- function(value) {
  value <- as.integer(value)
  .v2_concat(list(.v2_count_payload(length(value)),
                  .v2_u32_bytes(value)))
}

.v2_double_vector_payload <- function(value) {
  value <- as.numeric(value)
  if (any(!is.finite(value))) .v2_fail("non-finite value in exact encoding")
  .v2_concat(list(.v2_count_payload(length(value)),
                  if (length(value)) {
                    connection <- rawConnection(raw(0L), open = "wb")
                    on.exit(close(connection), add = TRUE)
                    writeBin(value, connection, size = 8L, endian = "big")
                    rawConnectionValue(connection)
                  } else raw(0L)))
}

.v2_graph <- function(tree) {
  if (!inherits(tree, "phylo") || !is.list(tree)) {
    .v2_fail("tree must be a phylo object")
  }
  labels <- .v2_label_bytes(tree$tip.label)
  n_tip <- length(labels$labels)
  if (n_tip < 2L) .v2_fail("tree must contain at least two tips")

  edge_input <- tree$edge
  if (!is.matrix(edge_input) || ncol(edge_input) != 2L ||
      !nrow(edge_input)) {
    .v2_fail("edge must be a non-empty two-column matrix")
  }
  edge_values <- suppressWarnings(as.numeric(edge_input))
  if (length(edge_values) != length(edge_input) || any(!is.finite(edge_values)) ||
      any(edge_values != floor(edge_values))) {
    .v2_fail("edge endpoints must be finite integer-valued numbers")
  }
  edge <- matrix(as.integer(edge_values), ncol = 2L)
  if (any(edge < 1L)) .v2_fail("edge endpoints must be positive")

  n_node <- tree$Nnode
  if (length(n_node) != 1L || !is.finite(n_node) || n_node < 1L ||
      n_node != floor(n_node)) {
    .v2_fail("Nnode must be a positive integer")
  }
  n_node <- as.integer(n_node)
  if (nrow(edge) != n_tip + n_node - 1L) {
    .v2_fail("edge row count is inconsistent with tip and internal-node counts")
  }

  edge_length <- tree$edge.length
  if (is.null(edge_length) || length(edge_length) != nrow(edge)) {
    .v2_fail("edge.length must be present with one value per edge")
  }
  edge_length <- suppressWarnings(as.numeric(edge_length))
  if (any(!is.finite(edge_length)) || any(edge_length < 0)) {
    .v2_fail("edge.length must be finite and non-negative")
  }

  ids <- sort(unique(as.integer(edge)))
  internal_ids <- ids[ids > n_tip]
  if (length(internal_ids) != n_node) {
    .v2_fail("edge endpoints do not define exactly Nnode internal IDs")
  }
  # Arbitrary positive internal IDs are supported; tip IDs retain the phylo
  # convention 1:n_tip.  Dense indices are internal implementation details.
  node_ids <- c(seq_len(n_tip), internal_ids)
  index_by_id <- stats::setNames(seq_along(node_ids), as.character(node_ids))
  parent_idx <- unname(index_by_id[as.character(edge[, 1L])])
  child_idx <- unname(index_by_id[as.character(edge[, 2L])])
  if (anyNA(parent_idx) || anyNA(child_idx)) {
    .v2_fail("edge endpoints include an unrecognised node ID")
  }
  parent_idx <- as.integer(parent_idx)
  child_idx <- as.integer(child_idx)

  degrees <- tabulate(parent_idx, nbins = length(node_ids))
  children <- lapply(degrees, integer)
  child_cursor <- integer(length(node_ids))
  for (i in seq_len(nrow(edge))) {
    p <- parent_idx[[i]]
    child_cursor[[p]] <- child_cursor[[p]] + 1L
    children[[p]][[child_cursor[[p]]]] <- child_idx[[i]]
  }
  indegree <- tabulate(child_idx, nbins = length(node_ids))
  tip_idx <- seq_len(n_tip)
  internal_idx <- n_tip + seq_len(n_node)
  if (any(indegree[tip_idx] != 1L)) {
    .v2_fail("each tip must have exactly one parent")
  }
  if (any(indegree[internal_idx] > 1L) || sum(indegree[internal_idx] == 0L) != 1L) {
    .v2_fail("internal nodes must have one structural root and no reticulation")
  }
  root_idx <- internal_idx[indegree[internal_idx] == 0L]
  if (length(root_idx) != 1L) .v2_fail("tree must have exactly one root")
  if (any(lengths(children[tip_idx]) != 0L)) {
    .v2_fail("tip nodes cannot be parents")
  }
  if (any(degrees[internal_idx] < 2L)) {
    .v2_fail("unary internal nodes are not valid rooted phylo trees")
  }

  edge_len_by_child <- rep(NA_real_, length(node_ids))
  edge_len_by_child[child_idx] <- edge_length
  parent_for_child <- rep(NA_integer_, length(node_ids))
  parent_for_child[child_idx] <- parent_idx

  # Reachability and cycle checks use an explicit stack.  No recursive call is
  # made at a depth proportional to the number of tips.
  colour <- integer(length(node_ids))
  stack_node <- integer(2L * length(node_ids) + 4L)
  stack_exit <- logical(length(stack_node))
  top <- 1L
  stack_node[[top]] <- root_idx
  max_stack <- top
  while (top > 0L) {
    node <- stack_node[[top]]
    exit <- stack_exit[[top]]
    top <- top - 1L
    if (exit) {
      colour[[node]] <- 2L
      next
    }
    if (colour[[node]] == 1L) .v2_fail("tree contains a directed cycle")
    if (colour[[node]] == 2L) next
    colour[[node]] <- 1L
    top <- top + 1L
    stack_node[[top]] <- node
    stack_exit[[top]] <- TRUE
    kids <- children[[node]]
    if (length(kids)) {
      for (kid in rev(kids)) {
        if (colour[[kid]] == 1L) .v2_fail("tree contains a directed cycle")
        top <- top + 1L
        stack_node[[top]] <- kid
        stack_exit[[top]] <- FALSE
      }
    }
    max_stack <- max(max_stack, top)
  }
  if (any(colour != 2L)) .v2_fail("tree is disconnected from its structural root")

  # Input-order postorder is used only to obtain bottom-up readiness.  The
  # subsequent child ordering is determined solely by UTF-8 tip ranks.
  postorder <- integer(length(node_ids))
  postorder_count <- 0L
  stack_node <- integer(2L * length(node_ids) + 4L)
  stack_exit <- logical(length(stack_node))
  top <- 1L
  stack_node[[top]] <- root_idx
  while (top > 0L) {
    node <- stack_node[[top]]
    exit <- stack_exit[[top]]
    top <- top - 1L
    if (exit) {
      postorder_count <- postorder_count + 1L
      postorder[[postorder_count]] <- node
      next
    }
    top <- top + 1L
    stack_node[[top]] <- node
    stack_exit[[top]] <- TRUE
    kids <- children[[node]]
    if (length(kids)) {
      for (kid in rev(kids)) {
        top <- top + 1L
        stack_node[[top]] <- kid
        stack_exit[[top]] <- FALSE
      }
    }
  }
  postorder <- postorder[seq_len(postorder_count)]

  tip_order <- .v2_order_bytes(labels$bytes)
  tip_rank <- integer(n_tip)
  tip_rank[tip_order] <- seq_len(n_tip)
  min_tip_rank <- rep(NA_integer_, length(node_ids))
  min_tip_rank[tip_idx] <- tip_rank
  ordered_children <- vector("list", length(node_ids))
  child_sort_calls <- 0L
  for (node in postorder) {
    if (node <= n_tip) next
    kids <- children[[node]]
    child_mins <- min_tip_rank[kids]
    if (anyNA(child_mins) || anyDuplicated(child_mins)) {
      .v2_fail("child subtrees do not have distinct canonical minimum tip ranks")
    }
    ordered_children[[node]] <- kids[order(child_mins, method = "radix")]
    min_tip_rank[[node]] <- min(child_mins)
    child_sort_calls <- child_sort_calls + 1L
  }

  # Canonical preorder assigns root n_tip+1 and then walks children by their
  # bytewise minimum tip rank.  This is independent of raw IDs and edge rows.
  canonical_preorder <- integer(length(node_ids))
  canonical_preorder_count <- 0L
  stack_node <- integer(length(node_ids) + 2L)
  top <- 1L
  stack_node[[top]] <- root_idx
  while (top > 0L) {
    node <- stack_node[[top]]
    top <- top - 1L
    canonical_preorder_count <- canonical_preorder_count + 1L
    canonical_preorder[[canonical_preorder_count]] <- node
    kids <- ordered_children[[node]]
    if (length(kids)) {
      for (kid in rev(kids)) {
        top <- top + 1L
        stack_node[[top]] <- kid
      }
    }
    max_stack <- max(max_stack, top)
  }
  canonical_preorder <- canonical_preorder[seq_len(canonical_preorder_count)]
  if (length(canonical_preorder) != length(node_ids)) {
    .v2_fail("canonical traversal did not cover every node")
  }

  canonical_id <- integer(length(node_ids))
  canonical_id[tip_idx] <- tip_idx
  internal_preorder <- canonical_preorder[canonical_preorder > n_tip]
  canonical_id[internal_preorder] <- seq.int(
    n_tip + 1L, length.out = length(internal_preorder)
  )
  descriptor_id <- canonical_id
  descriptor_id[tip_idx] <- tip_rank

  # Canonical postorder is also explicit.  It gives a deterministic edge-row
  # order while preserving each child edge's original length by child identity.
  canonical_postorder <- integer(length(node_ids))
  canonical_postorder_count <- 0L
  stack_node <- integer(2L * length(node_ids) + 4L)
  stack_exit <- logical(length(stack_node))
  top <- 1L
  stack_node[[top]] <- root_idx
  while (top > 0L) {
    node <- stack_node[[top]]
    exit <- stack_exit[[top]]
    top <- top - 1L
    if (exit) {
      canonical_postorder_count <- canonical_postorder_count + 1L
      canonical_postorder[[canonical_postorder_count]] <- node
      next
    }
    top <- top + 1L
    stack_node[[top]] <- node
    stack_exit[[top]] <- TRUE
    kids <- ordered_children[[node]]
    if (length(kids)) {
      for (kid in rev(kids)) {
        top <- top + 1L
        stack_node[[top]] <- kid
        stack_exit[[top]] <- FALSE
      }
    }
  }
  canonical_postorder <- canonical_postorder[seq_len(canonical_postorder_count)]

  # Exact descriptor DAG. Each node is stored once. Internal records contain
  # canonical child references and never copied descendant-label vectors.
  descriptor <- vector("list", length(node_ids))
  descriptor_builds <- 0L
  for (node in canonical_postorder) {
    cid <- descriptor_id[[node]]
    if (node <= n_tip) {
      descriptor[[cid]] <- list(
        kind = "tip",
        label_rank = as.integer(tip_rank[[node]]),
        label_bytes = labels$bytes[[node]],
        child_descriptor_ids = integer()
      )
    } else {
      kids <- ordered_children[[node]]
      descriptor[[cid]] <- list(
        kind = "internal",
        child_count = as.integer(length(kids)),
        child_min_tip_ranks = as.integer(min_tip_rank[kids]),
        child_descriptor_ids = as.integer(descriptor_id[kids])
      )
    }
    descriptor_builds <- descriptor_builds + 1L
  }

  total_tuple_entries <- sum(degrees[internal_idx])
  label_bytes_total <- sum(vapply(labels$bytes, length, integer(1L)))
  edge_payload_bytes <- 4L * total_tuple_entries + 8L * total_tuple_entries
  metrics <- list(
    schema_version = 2L,
    child_rank_tuple_entries = as.integer(total_tuple_entries),
    total_child_rank_tuple_entries = as.integer(total_tuple_entries),
    descriptor_entries = as.integer(length(descriptor)),
    total_stored_descriptor_entries = as.integer(length(descriptor)),
    peak_internal_object_proxy = as.integer(max(c(1L, degrees[internal_idx]))),
    descriptor_payload_bytes = as.numeric(label_bytes_total + edge_payload_bytes),
    descriptor_build_calls = as.integer(descriptor_builds),
    child_sort_calls = as.integer(child_sort_calls),
    explicit_stack_peak = as.integer(max_stack),
    full_edge_scans = 1L,
    descendant_label_vectors_materialized = 0L
  )

  list(
    labels = labels,
    n_tip = n_tip,
    n_node = n_node,
    edge = edge,
    edge_length = edge_length,
    node_ids = node_ids,
    children = children,
    parent_for_child = parent_for_child,
    edge_len_by_child = edge_len_by_child,
    root_idx = root_idx,
    degrees = degrees,
    tip_rank = tip_rank,
    tip_order = tip_order,
    min_tip_rank = min_tip_rank,
    ordered_children = ordered_children,
    canonical_id = canonical_id,
    descriptor_id = descriptor_id,
    canonical_preorder = canonical_preorder,
    canonical_postorder = canonical_postorder,
    descriptor = descriptor,
    metrics = metrics
  )
}

v2_canonicalize <- function(tree) {
  graph <- .v2_graph(tree)
  nonroot <- graph$canonical_postorder[
    graph$canonical_postorder != graph$root_idx
  ]
  canonical_edge <- cbind(
    as.integer(graph$canonical_id[graph$parent_for_child[nonroot]]),
    as.integer(graph$canonical_id[nonroot])
  )
  storage.mode(canonical_edge) <- "integer"
  canonical_lengths <- as.numeric(graph$edge_len_by_child[nonroot])

  # Keep only deterministic computational phylo fields.  node.label and other
  # metadata are deliberately outside the V2 protected representation.
  out <- list(
    edge = canonical_edge,
    tip.label = graph$labels$labels,
    Nnode = as.integer(graph$n_node),
    edge.length = canonical_lengths
  )
  class(out) <- "phylo"
  attr(out, "fastphylosig_v2_metrics") <- graph$metrics
  attr(out, "fastphylosig_v2_contract") <- list(
    schema_version = 2L,
    internal_mapping = seq.int(
      graph$n_tip + 1L, graph$n_tip + graph$n_node
    ),
    edge_order = seq_len(nrow(canonical_edge)),
    root = as.integer(graph$n_tip + 1L),
    root_id = as.integer(graph$n_tip + 1L),
    subtree_descriptors = graph$descriptor,
    internal_numbering = "canonical_preorder",
    edge_order_rule = "canonical_postorder",
    label_order = "UTF-8 unsigned-byte lexicographic"
  )
  out
}

v2_fingerprint <- function(tree) {
  canonical <- v2_canonicalize(tree)
  labels <- .v2_label_bytes(canonical$tip.label)
  encoded <- .v2_concat(list(
    .v2_field("domain", charToRaw("fastphylosig.canonical.tree")),
    .v2_field("schema_version", .v2_u32_bytes(2L)),
    .v2_field("edge_order", charToRaw("canonical_postorder")),
    .v2_field("tip.label", .v2_label_vector_payload(labels$bytes)),
    .v2_field("Nnode", .v2_integer_vector_payload(canonical$Nnode)),
    .v2_field("edge_dim", .v2_integer_vector_payload(dim(canonical$edge))),
    .v2_field("edge", .v2_integer_vector_payload(as.integer(t(canonical$edge)))),
    .v2_field("edge.length", .v2_double_vector_payload(canonical$edge.length))
  ))
  # The returned value is an exact encoded state, not a digest.  Hex is an
  # ASCII presentation of the unambiguous byte stream and cannot hide a hash
  # collision or delimiter ambiguity.
  .v2_hex(encoded)
}

v2_protected_snapshot <- function(tree) {
  graph <- .v2_graph(tree)
  encoded <- .v2_concat(list(
    .v2_field("domain", charToRaw("fastphylosig.protected.source.tree")),
    .v2_field("schema_version", .v2_u32_bytes(2L)),
    .v2_field("tip.label", .v2_label_vector_payload(graph$labels$bytes)),
    .v2_field("edge_dim", .v2_integer_vector_payload(dim(graph$edge))),
    .v2_field("edge", .v2_integer_vector_payload(as.integer(t(graph$edge)))),
    .v2_field("edge.length", .v2_double_vector_payload(graph$edge_length)),
    .v2_field("Nnode", .v2_integer_vector_payload(graph$n_node))
  ))
  list(
    schema_version = 2L,
    encoded = encoded,
    edge_dim = as.integer(dim(graph$edge)),
    edge = unname(graph$edge),
    edge_length = unname(graph$edge_length),
    tip_label_bytes = graph$labels$bytes,
    Nnode = as.integer(graph$n_node)
  )
}

.v2_snapshot_equal <- function(left, right) {
  is.list(left) && is.list(right) &&
    identical(left$schema_version, right$schema_version) &&
    identical(left$encoded, right$encoded) &&
    identical(left$edge_dim, right$edge_dim) &&
    identical(left$edge, right$edge) &&
    identical(left$edge_length, right$edge_length) &&
    identical(left$tip_label_bytes, right$tip_label_bytes) &&
    identical(left$Nnode, right$Nnode)
}

.v2_context_registry <- new.env(hash = TRUE, parent = emptyenv())
.v2_context_counter <- 0L

v2_prepare_context <- function(tree) {
  canonical <- v2_canonicalize(tree)
  snapshot <- v2_protected_snapshot(canonical)
  .v2_context_counter <<- .v2_context_counter + 1L
  token <- paste0("v2-context-", .v2_context_counter)
  assign(token, unserialize(serialize(snapshot, NULL, version = 2L)),
         envir = .v2_context_registry)
  out <- list(
    context_schema_version = 2L,
    canonical_contract_version = 2L,
    protected_snapshot_version = 2L,
    package_token = token,
    tree = canonical,
    protected_snapshot = snapshot,
    fingerprint = v2_fingerprint(canonical),
    descriptor_metrics = attr(canonical, "fastphylosig_v2_metrics", exact = TRUE)
  )
  class(out) <- c("fastphylosig_v2_context", "list")
  out
}

v2_validate_context <- function(ctx) {
  if (!is.list(ctx) || !inherits(ctx, "fastphylosig_v2_context")) {
    .v2_fail("context is not a V2 context; call v2_prepare_context()")
  }
  if (!identical(ctx$context_schema_version, 2L) ||
      !identical(ctx$canonical_contract_version, 2L) ||
      !identical(ctx$protected_snapshot_version, 2L)) {
    .v2_fail("context schema is unknown or not version 2; re-run v2_prepare_context()")
  }
  if (!is.list(ctx$protected_snapshot) ||
      !is.character(ctx$fingerprint) || length(ctx$fingerprint) != 1L) {
    .v2_fail("context is missing its V2 protected snapshot or fingerprint")
  }
  if (!inherits(ctx$tree, "phylo")) {
    .v2_fail("context tree is missing or is not a phylo object")
  }
  if (!is.character(ctx$package_token) || length(ctx$package_token) != 1L ||
      !exists(ctx$package_token, envir = .v2_context_registry,
              inherits = FALSE)) {
    .v2_fail("context token is missing or unknown; re-run v2_prepare_context()")
  }
  owned_snapshot <- get(ctx$package_token, envir = .v2_context_registry,
                        inherits = FALSE)
  if (!.v2_snapshot_equal(ctx$protected_snapshot, owned_snapshot)) {
    .v2_fail("stored protected snapshot changed; re-run v2_prepare_context()")
  }
  current_snapshot <- v2_protected_snapshot(ctx$tree)
  if (!.v2_snapshot_equal(current_snapshot, owned_snapshot)) {
    .v2_fail("protected tree fields changed; re-run v2_prepare_context()")
  }
  current_fingerprint <- v2_fingerprint(ctx$tree)
  if (!identical(current_fingerprint, ctx$fingerprint)) {
    .v2_fail("V2 fingerprint changed; re-run v2_prepare_context()")
  }
  TRUE
}
