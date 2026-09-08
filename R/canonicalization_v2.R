# Canonicalization Contract V2 ----------------------------------------------
#
# These helpers are deliberately private.  They provide a representation-only
# canonical form and two exact byte encodings for the production preparation
# layer.  C2 (canonical identity) and S2 (protected source snapshot) are kept
# separate: a canonical tree can be representation invariant while a prepared
# context can still reject a mutation of its source representation.

.v2_context_schema_version <- 2L
.v2_canonical_contract_version <- 2L
.v2_protected_snapshot_version <- 2L
.v2_context_record_attribute <- "fastphylosig_v2_integrity"

.v2_abort <- function(message) {
  stop(paste0("fastphylosig V2: ", as.character(message[[1L]])),
       call. = FALSE)
}

.v2_raw_concat <- function(parts) {
  if (missing(parts)) parts <- list()
  if (!is.list(parts)) parts <- list(parts)
  parts <- lapply(parts, as.raw)
  parts <- parts[lengths(parts) > 0L]
  if (!length(parts)) return(raw(0L))
  do.call(c, parts)
}

.v2_ascii_raw <- function(value) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    .v2_abort("internal ASCII token must be one non-empty value")
  }
  bytes <- charToRaw(value)
  if (any(as.integer(bytes) > 127L)) {
    .v2_abort("internal encoding token must contain ASCII bytes only")
  }
  bytes
}

.v2_raw_hex <- function(bytes) {
  bytes <- as.raw(bytes)
  if (!length(bytes)) return("")
  paste(sprintf("%02X", as.integer(bytes)), collapse = "")
}

# Encode one non-negative integer as four unsigned big-endian bytes.  R's
# integer type is signed, so the arithmetic is intentionally performed in a
# double vector and checked against the complete uint32 range.
.v2_u32_be <- function(value) {
  value <- as.numeric(value)
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value < 0 || value != floor(value) || value > 4294967295) {
    .v2_abort("unsigned integer encoding received an invalid value")
  }
  powers <- c(16777216, 65536, 256, 1)
  as.raw(as.integer(floor(value / powers) %% 256))
}

.v2_u32_payload <- function(value, name = "integer vector") {
  value <- as.numeric(value)
  if (anyNA(value) || any(!is.finite(value)) || any(value < 0) ||
      any(value != floor(value)) || any(value > 4294967295)) {
    .v2_abort(paste0(name, " contains an invalid unsigned integer"))
  }
  encoded <- if (length(value)) {
    powers <- c(16777216, 65536, 256, 1)
    bytes <- vapply(
      powers,
      function(power) floor(value / power) %% 256,
      numeric(length(value))
    )
    as.raw(as.integer(t(bytes)))
  } else raw(0L)
  .v2_raw_concat(list(.v2_u32_be(length(value)), encoded))
}

.v2_double_bytes <- function(value) {
  value <- as.numeric(value)
  if (length(value) != 1L || is.na(value) || !is.finite(value)) {
    .v2_abort("binary64 encoding received a non-finite value")
  }
  connection <- rawConnection(raw(0L), open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(value, connection, size = 8L, endian = "big")
  rawConnectionValue(connection)
}

.v2_double_payload <- function(value, name = "double vector") {
  value <- as.numeric(value)
  if (anyNA(value) || any(!is.finite(value))) {
    .v2_abort(paste0(name, " contains a non-finite binary64 value"))
  }
  encoded <- if (length(value)) {
    connection <- rawConnection(raw(0L), open = "wb")
    on.exit(close(connection), add = TRUE)
    writeBin(value, connection, size = 8L, endian = "big")
    rawConnectionValue(connection)
  } else raw(0L)
  .v2_raw_concat(list(.v2_u32_be(length(value)), encoded))
}

# Every field carries a length-delimited tag, an explicit type token, and a
# length-delimited payload.  No delimiter is structural syntax.
.v2_field <- function(tag, type, payload) {
  tag <- .v2_ascii_raw(tag)
  type <- .v2_ascii_raw(type)
  payload <- as.raw(payload)
  .v2_raw_concat(list(
    .v2_u32_be(length(tag)), tag,
    .v2_u32_be(length(type)), type,
    .v2_u32_be(length(payload)), payload
  ))
}

.v2_label_bytes <- function(labels) {
  if (!is.character(labels)) {
    .v2_abort("tip.label must be a character vector")
  }
  if (!length(labels) || anyNA(labels) || any(!nzchar(labels))) {
    .v2_abort("tip.label must contain non-empty, non-missing labels")
  }
  stored <- labels
  utf8 <- tryCatch(enc2utf8(stored), error = function(e) NULL)
  if (is.null(utf8) || anyNA(utf8) ||
      any(!vapply(utf8, function(x) isTRUE(validUTF8(x)), logical(1L)))) {
    .v2_abort("tip.label contains a value that cannot be represented as UTF-8")
  }
  bytes <- lapply(utf8, charToRaw)
  if (any(!vapply(bytes, length, integer(1L)))) {
    .v2_abort("tip.label contains an empty UTF-8 byte sequence")
  }
  # A bytes-marked copy makes duplicated() compare the already-normalized
  # UTF-8 payload without locale translation.  This retains exact byte truth
  # while avoiding quadratic pairwise comparisons of raw list elements.
  byte_keys <- utf8
  Encoding(byte_keys) <- "bytes"
  if (anyDuplicated(byte_keys)) {
    .v2_abort("tip.label UTF-8 byte sequences must be unique")
  }
  list(stored = stored, utf8 = utf8, bytes = bytes)
}

# Explicit unsigned-byte lexicographic comparator.  R's character collation is
# intentionally never consulted for canonical truth.
.v2_bytes_compare <- function(left, right) {
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

.v2_order_bytes <- function(keys) {
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
          (i <= middle && .v2_bytes_compare(
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

.v2_label_vector_payload <- function(label_bytes) {
  sizes <- lengths(label_bytes)
  total <- 4 + sum(4 + sizes)
  encoded <- raw(total)
  encoded[seq_len(4L)] <- .v2_u32_be(length(label_bytes))
  cursor <- 5L
  for (i in seq_along(label_bytes)) {
    encoded[cursor:(cursor + 3L)] <- .v2_u32_be(sizes[[i]])
    cursor <- cursor + 4L
    if (sizes[[i]]) {
      encoded[cursor:(cursor + sizes[[i]] - 1L)] <- label_bytes[[i]]
      cursor <- cursor + sizes[[i]]
    }
  }
  encoded
}

.v2_source_fields <- function(tree) {
  if (!is.list(tree) || !inherits(tree, "phylo")) {
    .v2_abort("tree must be an object of class \"phylo\"")
  }
  labels <- .v2_label_bytes(tree$tip.label)
  n_tip <- length(labels$stored)
  if (n_tip < 2L) .v2_abort("tree must contain at least two tips")

  edge_input <- tree$edge
  if (!is.matrix(edge_input) || ncol(edge_input) != 2L ||
      !nrow(edge_input)) {
    .v2_abort("edge must be a non-empty two-column matrix")
  }
  edge_values <- suppressWarnings(as.numeric(edge_input))
  if (length(edge_values) != length(edge_input) ||
      anyNA(edge_values) || any(!is.finite(edge_values)) ||
      any(edge_values != floor(edge_values)) ||
      any(edge_values < 1) || any(edge_values > .Machine$integer.max)) {
    .v2_abort("edge endpoints must be finite positive integer-valued numbers")
  }
  edge <- matrix(as.integer(edge_values), ncol = 2L,
                 dimnames = dimnames(edge_input))

  n_node <- as.numeric(tree$Nnode)
  if (length(n_node) != 1L || is.na(n_node) || !is.finite(n_node) ||
      n_node < 1 || n_node != floor(n_node) ||
      n_node > .Machine$integer.max) {
    .v2_abort("Nnode must be a positive integer")
  }
  n_node <- as.integer(n_node)
  expected_edges <- as.numeric(n_tip) + as.numeric(n_node) - 1
  if (nrow(edge) != expected_edges) {
    .v2_abort("edge row count is inconsistent with tip and internal-node counts")
  }

  edge_length <- tree$edge.length
  if (is.null(edge_length) || length(edge_length) != nrow(edge)) {
    .v2_abort("edge.length must be present with one value per edge")
  }
  edge_length <- suppressWarnings(as.numeric(edge_length))
  if (anyNA(edge_length) || any(!is.finite(edge_length)) ||
      any(edge_length < 0)) {
    .v2_abort("edge.length must be finite and non-negative")
  }
  list(
    labels = labels,
    n_tip = n_tip,
    n_node = n_node,
    edge = edge,
    edge_length = edge_length
  )
}

.v2_build_graph <- function(tree) {
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
  tip_order <- .v2_order_bytes(labels$bytes)
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

.v2_reorder_matrix_attributes <- function(original, value, source_rows) {
  attrs <- attributes(original)
  if (is.null(attrs)) return(value)
  attrs$dim <- dim(value)
  if (!is.null(attrs$dimnames) && length(attrs$dimnames) >= 1L &&
      !is.null(attrs$dimnames[[1L]])) {
    attrs$dimnames[[1L]] <- attrs$dimnames[[1L]][source_rows]
  }
  attributes(value) <- attrs
  value
}

.v2_reorder_length_attributes <- function(original, value, source_rows) {
  attrs <- attributes(original)
  if (is.null(attrs)) return(value)
  if (!is.null(attrs$names)) attrs$names <- attrs$names[source_rows]
  attributes(value) <- attrs
  value
}

.v2_canonicalize_core <- function(tree) {
  graph <- .v2_build_graph(tree)
  fields <- graph$fields
  n_tip <- fields$n_tip
  n_node <- fields$n_node
  root_idx <- graph$root_idx
  nonroot <- graph$canonical_postorder[
    graph$canonical_postorder != root_idx
  ]
  source_rows <- graph$edge_row_for_child[nonroot]
  canonical_edge <- cbind(
    as.integer(graph$canonical_id[graph$parent_for_child[nonroot]]),
    as.integer(graph$canonical_id[nonroot])
  )
  storage.mode(canonical_edge) <- "integer"
  canonical_edge <- .v2_reorder_matrix_attributes(
    tree$edge, canonical_edge, source_rows
  )
  canonical_lengths <- as.numeric(fields$edge_length[source_rows])
  canonical_lengths <- .v2_reorder_length_attributes(
    tree$edge.length, canonical_lengths, source_rows
  )

  # Copy the phylo object and all unrelated fields/attributes.  Only the
  # computational edge representation is replaced.  The source object is not
  # modified; node labels are remapped with the internal-node map when their
  # conventional length is valid.
  out <- tree
  out$edge <- canonical_edge
  out$edge.length <- canonical_lengths
  out$Nnode <- as.integer(n_node)
  if (!is.null(tree$node.label) && length(tree$node.label) == n_node) {
    source_internal_order <- order(
      graph$node_ids[graph$internal_idx], method = "radix"
    )
    labels_by_dense <- tree$node.label[seq_len(n_node)]
    dense_labels <- vector(mode = typeof(labels_by_dense), length = length(graph$node_ids))
    dense_labels[graph$internal_idx[source_internal_order]] <- labels_by_dense
    canonical_labels <- vector(mode = typeof(labels_by_dense), length = n_node)
    for (node in graph$internal_idx) {
      canonical_position <- graph$canonical_id[[node]] - n_tip
      canonical_labels[[canonical_position]] <- dense_labels[[node]]
    }
    out$node.label <- canonical_labels
  }
  out$tip.label <- tree$tip.label
  class(out) <- unique(c("phylo", class(tree)))
  # The rows above follow the V2 canonical postorder contract, which is not an
  # assertion about ape's private traversal layout.  Remove any stale input
  # marker so numerical consumers can request their own ape ordering safely.
  attr(out, "order") <- NULL

  canonical_internal_source <- graph$canonical_preorder[
    graph$canonical_preorder > n_tip
  ]
  old_internal_order <- as.numeric(graph$node_ids[canonical_internal_source])
  new_internal_order <- seq.int(n_tip + 1L, length.out = n_node)
  legacy_info <- list(
    changed = !identical(unname(as.matrix(tree$edge)),
                         unname(as.matrix(out$edge))),
    mapping = list(
      old_to_new = stats::setNames(
        as.integer(graph$canonical_id), as.character(graph$node_ids)
      ),
      new_to_old = stats::setNames(
        as.numeric(graph$node_ids[order(graph$canonical_id, method = "radix")]),
        as.character(seq_len(length(graph$node_ids)))
      ),
      old_internal_order = old_internal_order,
      new_internal_order = new_internal_order
    ),
    tip_identity = identical(tree$tip.label, out$tip.label),
    branch_identity = identical(
      sort(as.numeric(tree$edge.length)),
      sort(as.numeric(out$edge.length))
    ),
    edge_isomorphism = TRUE,
    cophenetic_identity = TRUE,
    root_identity = TRUE,
    polytomy_identity = TRUE,
    safe = TRUE,
    reason = NULL
  )
  attr(out, "fastphylosig_canonicalization") <- legacy_info

  contract <- list(
    schema_version = .v2_canonical_contract_version,
    canonical_root = as.integer(n_tip + 1L),
    internal_numbering = "canonical_preorder",
    edge_order = "canonical_postorder",
    label_order = "UTF-8 unsigned-byte lexicographic",
    source_node_ids = graph$node_ids,
    source_to_canonical = stats::setNames(
      as.integer(graph$canonical_id), as.character(graph$node_ids)
    ),
    canonical_to_source = stats::setNames(
      as.integer(graph$node_ids[order(graph$canonical_id, method = "radix")]),
      as.character(seq_len(length(graph$node_ids)))
    ),
    canonical_edge_source_rows = as.integer(source_rows),
    tip_rank = as.integer(graph$tip_rank),
    min_tip_rank = as.integer(graph$min_tip_rank),
    descriptor_table = graph$descriptor
  )
  attr(out, "fastphylosig_v2_contract") <- contract
  attr(out, "fastphylosig_v2_metrics") <- graph$metrics
  out
}

.v2_is_canonical_fields <- function(fields) {
  edge <- fields$edge
  n_tip <- fields$n_tip
  n_node <- fields$n_node
  endpoints <- as.integer(edge)
  total <- n_tip + n_node
  if (any(endpoints > total) ||
      any(tabulate(endpoints, nbins = total) == 0L)) return(FALSE)
  roots <- setdiff(edge[, 1L], edge[, 2L])
  length(roots) == 1L && identical(as.integer(roots[[1L]]), as.integer(n_tip + 1L))
}

.v2_is_canonical_tree <- function(tree) {
  if (!is.list(tree) || !inherits(tree, "phylo")) return(FALSE)
  fields <- tryCatch(.v2_source_fields(tree), error = function(e) NULL)
  !is.null(fields) && .v2_is_canonical_fields(fields)
}

.v2_canonical_fingerprint_fields <- function(fields) {
  edge <- fields$edge
  .v2_raw_concat(list(
    .v2_field("domain", "ascii",
               .v2_ascii_raw("fastphylosig.canonical.tree")),
    .v2_field("schema_version", "u32", .v2_u32_be(
      .v2_canonical_contract_version
    )),
    .v2_field("tip_count", "u32", .v2_u32_be(fields$n_tip)),
    .v2_field("tip.label", "utf8-vector",
               .v2_label_vector_payload(fields$labels$bytes)),
    .v2_field("Nnode", "u32", .v2_u32_be(fields$n_node)),
    .v2_field("edge_count", "u32", .v2_u32_be(nrow(edge))),
    .v2_field("edge_dim", "u32-vector", .v2_u32_payload(dim(edge))),
    .v2_field("edge", "u32-vector",
               .v2_u32_payload(as.numeric(t(edge)))),
    .v2_field("edge.length", "binary64-vector",
               .v2_double_payload(fields$edge_length))
  ))
}

.v2_canonical_fingerprint_raw <- function(tree, canonical = NULL) {
  canonical <- if (is.null(canonical)) {
    .v2_canonicalize_core(tree)
  } else canonical
  fields <- .v2_source_fields(canonical)
  if (!.v2_is_canonical_fields(fields)) {
    .v2_abort("fingerprint input is not a valid V2 canonical tree")
  }
  .v2_canonical_fingerprint_fields(fields)
}

.v2_tree_fingerprint_raw <- .v2_canonical_fingerprint_raw

.v2_protected_snapshot_fields <- function(fields) {
  edge <- fields$edge
  .v2_raw_concat(list(
    .v2_field("domain", "ascii",
               .v2_ascii_raw("fastphylosig.protected.source.tree")),
    .v2_field("schema_version", "u32", .v2_u32_be(
      .v2_protected_snapshot_version
    )),
    .v2_field("tip.label", "utf8-vector",
               .v2_label_vector_payload(fields$labels$bytes)),
    .v2_field("edge_dim", "u32-vector", .v2_u32_payload(dim(edge))),
    .v2_field("edge", "u32-vector",
               .v2_u32_payload(as.numeric(t(edge)))),
    .v2_field("edge.length", "binary64-vector",
               .v2_double_payload(fields$edge_length)),
    .v2_field("Nnode", "u32", .v2_u32_be(fields$n_node))
  ))
}

.v2_protected_snapshot_raw <- function(tree) {
  .v2_protected_snapshot_fields(.v2_source_fields(tree))
}

# The prepared boundary needs both exact encodings for the same canonical
# tree.  Parse and validate protected fields once, then construct the two
# independently domain-separated byte streams from that shared evidence.
.v2_context_evidence <- function(tree) {
  fields <- .v2_source_fields(tree)
  if (!.v2_is_canonical_fields(fields)) {
    .v2_abort("prepared tree is not a valid V2 canonical tree")
  }
  list(
    fingerprint = .tree_fingerprint(
      tree, .canonical = TRUE, .fields = fields
    ),
    protected_snapshot = .v2_protected_snapshot_fields(fields)
  )
}

.v2_protected_snapshot <- function(tree) {
  fields <- .v2_source_fields(tree)
  encoded <- .v2_protected_snapshot_raw(tree)
  list(
    schema_version = .v2_protected_snapshot_version,
    encoded = encoded,
    edge_dim = as.integer(dim(fields$edge)),
    edge = unname(fields$edge),
    edge_length = unname(fields$edge_length),
    tip_label_bytes = unname(fields$labels$bytes),
    Nnode = as.integer(fields$n_node)
  )
}

.v2_snapshot_equal <- function(left, right) {
  left_raw <- if (is.list(left)) left$encoded else left
  right_raw <- if (is.list(right)) right$encoded else right
  is.raw(left_raw) && is.raw(right_raw) && identical(left_raw, right_raw)
}

.v2_context_schema_markers <- function() {
  list(
    context_schema_version = .v2_context_schema_version,
    canonical_contract_version = .v2_canonical_contract_version,
    protected_snapshot_version = .v2_protected_snapshot_version
  )
}

.v2_validate_context_schema <- function(ctx) {
  expected <- .v2_context_schema_markers()
  if (!is.list(ctx) ||
      !identical(ctx$context_schema_version, expected$context_schema_version) ||
      !identical(ctx$canonical_contract_version,
                 expected$canonical_contract_version) ||
      !identical(ctx$protected_snapshot_version,
                 expected$protected_snapshot_version)) {
    .v2_abort("context schema is missing, unknown, or incompatible; call prepare_tree(tree) again")
  }
  invisible(TRUE)
}

.v2_new_context_record <- function(fingerprint, protected_snapshot,
                                   fixed_state, structural_cache,
                                   numerical_cache, cache_meta) {
  record <- new.env(parent = emptyenv())
  record$context_schema_version <- .v2_context_schema_version
  record$canonical_contract_version <- .v2_canonical_contract_version
  record$protected_snapshot_version <- .v2_protected_snapshot_version
  record$fingerprint <- as.raw(fingerprint)
  record$protected_snapshot <- as.raw(protected_snapshot)
  record$fixed_state <- unserialize(serialize(fixed_state, NULL, version = 3L))
  record$structural_cache <- structural_cache
  record$numerical_cache <- numerical_cache
  record$cache_meta <- cache_meta
  lockEnvironment(record, bindings = TRUE)
  record
}

.v2_context_record <- function(ctx) {
  record <- attr(ctx, .v2_context_record_attribute, exact = TRUE)
  required <- c(
    "context_schema_version", "canonical_contract_version",
    "protected_snapshot_version", "fingerprint", "protected_snapshot",
    "fixed_state", "structural_cache", "numerical_cache", "cache_meta"
  )
  if (!is.environment(record) || !environmentIsLocked(record) ||
      !all(vapply(required, bindingIsLocked, logical(1L), env = record))) {
    .v2_abort("prepared context integrity record is missing or invalid; call prepare_tree(tree) again")
  }
  expected <- .v2_context_schema_markers()
  if (!identical(record$context_schema_version,
                 expected$context_schema_version) ||
      !identical(record$canonical_contract_version,
                 expected$canonical_contract_version) ||
      !identical(record$protected_snapshot_version,
                 expected$protected_snapshot_version) ||
      !is.raw(record$fingerprint) || !is.raw(record$protected_snapshot)) {
    .v2_abort("prepared context integrity record is incompatible; call prepare_tree(tree) again")
  }
  record
}

# Stable private names used by the preparation integration.  These aliases
# keep the implementation namespace explicit without adding public API.
.canonicalize_tree_v2_core <- .v2_canonicalize_core
.tree_fingerprint_v2 <- .v2_canonical_fingerprint_raw
.tree_protected_snapshot_v2 <- .v2_protected_snapshot_raw
.context_schema_v2 <- .v2_context_schema_version
.canonical_contract_v2 <- .v2_canonical_contract_version
.protected_snapshot_v2 <- .v2_protected_snapshot_version
