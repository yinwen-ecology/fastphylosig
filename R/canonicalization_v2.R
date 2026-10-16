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
.v2_order_bytes <- function(keys) {
  canonical_order_bytes_cpp(keys)
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

# Source validation and attribute preservation stay in R; native traversal
# returns the identical V2 graph contract without descendant-label vectors.
.v2_build_graph <- function(tree) {
  fields <- .v2_source_fields(tree)
  graph <- canonical_tree_graph_cpp(
    fields$edge, fields$edge_length, fields$n_tip, fields$n_node,
    .v2_order_bytes(fields$labels$bytes)
  )
  c(list(fields = fields), graph)
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
