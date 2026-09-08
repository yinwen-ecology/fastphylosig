# Canonicalization Contract V2 production gates -----------------------------

# These tests exercise the private production boundary through the public
# preparation helpers.  The V2-A prototype remains an audit oracle only; no
# prototype code is sourced here.

.v2_prod_make_tree <- function(
    labels = c("a", "b", "c", "d"),
    edge_length = c(10.1, 20.2, 1.1, 2.2, 3.3, 4.4)) {
  stopifnot(length(labels) == 4L, length(edge_length) == 6L)
  tree <- list(
    edge = matrix(
      c(7L, 5L, 7L, 6L, 5L, 1L, 5L, 2L, 6L, 3L, 6L, 4L),
      ncol = 2L,
      byrow = TRUE
    ),
    edge.length = as.numeric(edge_length),
    Nnode = 3L,
    tip.label = as.character(labels)
  )
  class(tree) <- "phylo"
  tree
}

.v2_prod_renumber_internals <- function(tree) {
  # The map deliberately changes all internal IDs and keeps the root nonlocal.
  id_map <- stats::setNames(
    c(1L, 2L, 3L, 4L, 12L, 9L, 15L),
    as.character(seq_len(7L))
  )
  out <- tree
  values <- unname(id_map[as.character(tree$edge)])
  if (anyNA(values)) stop("test fixture contains an unmapped node", call. = FALSE)
  out$edge <- matrix(
    as.integer(values), nrow = nrow(tree$edge), ncol = 2L,
    dimnames = dimnames(tree$edge)
  )
  out
}

.v2_prod_shuffle_edges <- function(tree) {
  order_rows <- c(3L, 6L, 1L, 5L, 2L, 4L)
  out <- tree
  out$edge <- tree$edge[order_rows, , drop = FALSE]
  out$edge.length <- tree$edge.length[order_rows]
  out
}

.v2_prod_canonical_fields <- function(tree) {
  canonical <- fastphylosig:::.safe_canonicalize_core(tree)
  contract <- attr(canonical, "fastphylosig_v2_contract", exact = TRUE)
  list(
    edge = unname(canonical$edge),
    edge.length = unname(canonical$edge.length),
    tip.label = canonical$tip.label,
    Nnode = canonical$Nnode,
    fingerprint = fastphylosig:::.tree_fingerprint(
      canonical, .canonical = TRUE
    ),
    ordering = contract[c(
      "schema_version", "canonical_root", "internal_numbering",
      "edge_order", "label_order", "tip_rank", "descriptor_table"
    )],
    tree = canonical
  )
}

.v2_prod_raw_hex <- function(value) {
  bytes <- charToRaw(enc2utf8(as.character(value)))
  paste(sprintf("%02X", as.integer(bytes)), collapse = "")
}

.v2_prod_descendant_key <- function(tree, node) {
  n_tip <- length(tree$tip.label)
  if (node <= n_tip) return(.v2_prod_raw_hex(tree$tip.label[[node]]))
  rows <- which(tree$edge[, 1L] == node)
  child_keys <- vapply(
    rows,
    function(i) .v2_prod_descendant_key(tree, tree$edge[i, 2L]),
    character(1L)
  )
  paste(sort(child_keys), collapse = ";")
}

.v2_prod_edge_records <- function(tree) {
  keys <- vapply(
    tree$edge[, 2L],
    function(node) .v2_prod_descendant_key(tree, node),
    character(1L)
  )
  stats::setNames(as.numeric(tree$edge.length), keys)
}

.v2_prod_locale_fields <- function(tree, locale) {
  previous <- Sys.getlocale("LC_COLLATE")
  on.exit(try(Sys.setlocale("LC_COLLATE", previous), silent = TRUE),
          add = TRUE)
  selected <- tryCatch(
    Sys.setlocale("LC_COLLATE", locale),
    warning = function(e) NA_character_,
    error = function(e) NA_character_
  )
  if (length(selected) != 1L || is.na(selected) || !nzchar(selected)) {
    return(NULL)
  }
  .v2_prod_canonical_fields(tree)
}

.v2_prod_encoding_pair <- function() {
  candidates <- list(
    c("alpha", intToUtf8(0x4e2d), "omega", "beta"),
    c("alpha", intToUtf8(0xe9), "omega", "beta"),
    c("alpha", "beta", "omega", "delta")
  )
  for (labels in candidates) {
    utf8 <- enc2utf8(labels)
    native <- utf8
    Encoding(native) <- "unknown"
    if (identical(
      vapply(native, .v2_prod_raw_hex, character(1L)),
      vapply(utf8, .v2_prod_raw_hex, character(1L))
    ) && identical(enc2utf8(native), utf8)) {
      return(list(utf8 = utf8, native = native))
    }
    converted <- tryCatch(
      iconv(utf8, from = "UTF-8", to = "", sub = NA_character_),
      error = function(e) NULL
    )
    if (!is.null(converted) && !anyNA(converted) && identical(
      vapply(converted, .v2_prod_raw_hex, character(1L)),
      vapply(utf8, .v2_prod_raw_hex, character(1L))
    ) && identical(enc2utf8(converted), utf8)) {
      return(list(utf8 = utf8, native = converted))
    }
  }
  NULL
}

test_that("V2 canonical fields and fingerprint are representation invariant", {
  tree <- .v2_prod_make_tree(c("a", "b\rc", "a\rb", "c"))
  variants <- list(
    original = tree,
    renumbered = .v2_prod_renumber_internals(tree),
    shuffled = .v2_prod_shuffle_edges(tree),
    combined = .v2_prod_shuffle_edges(.v2_prod_renumber_internals(tree))
  )
  fields <- lapply(variants, .v2_prod_canonical_fields)
  reference <- fields[[1L]]
  for (candidate in fields[-1L]) {
    expect_identical(candidate$edge, reference$edge)
    expect_identical(candidate$edge.length, reference$edge.length)
    expect_identical(candidate$tip.label, reference$tip.label)
    expect_identical(candidate$Nnode, reference$Nnode)
    expect_identical(candidate$fingerprint, reference$fingerprint)
    expect_identical(candidate$ordering, reference$ordering)
  }
  expect_identical(reference$Nnode, 3L)
  expect_true(any(reference$edge[, 1L] == 5L))
  expect_identical(
    sort(unique(as.integer(reference$edge))),
    seq_len(length(tree$tip.label) + tree$Nnode)
  )
})

test_that("V2 preserves heterogeneous biological edge-length association", {
  tree <- .v2_prod_make_tree()
  canonical <- .v2_prod_canonical_fields(tree)$tree
  source_records <- .v2_prod_edge_records(tree)
  canonical_records <- .v2_prod_edge_records(canonical)
  source_records <- source_records[order(names(source_records))]
  canonical_records <- canonical_records[order(names(canonical_records))]
  expect_identical(canonical_records, source_records)
  expect_length(unique(canonical_records), nrow(tree$edge))
})

test_that("V2 canonical order and fingerprint do not depend on LC_COLLATE", {
  tree <- .v2_prod_make_tree(c("zeta", "A", "10", "punct!"))
  all_before <- Sys.getlocale()
  current <- Sys.getlocale("LC_COLLATE")
  current_fields <- .v2_prod_locale_fields(tree, current)
  c_fields <- .v2_prod_locale_fields(tree, "C")
  if (is.null(current_fields) || is.null(c_fields)) {
    skip("current locale or C LC_COLLATE is unavailable")
  }
  expect_identical(c_fields$edge, current_fields$edge)
  expect_identical(c_fields$edge.length, current_fields$edge.length)
  expect_identical(c_fields$tip.label, current_fields$tip.label)
  expect_identical(c_fields$Nnode, current_fields$Nnode)
  expect_identical(c_fields$fingerprint, current_fields$fingerprint)
  expect_identical(Sys.getlocale(), all_before)
})

test_that("V2 treats equivalent UTF-8 encodings identically and rejects byte duplicates", {
  pair <- .v2_prod_encoding_pair()
  if (is.null(pair)) skip("no equivalent non-lossy encoding available")

  tree_utf8 <- .v2_prod_make_tree(pair$utf8)
  tree_native <- .v2_prod_make_tree(pair$native)
  utf8_fields <- .v2_prod_canonical_fields(tree_utf8)
  native_fields <- .v2_prod_canonical_fields(tree_native)
  expect_identical(native_fields$edge, utf8_fields$edge)
  expect_identical(native_fields$edge.length, utf8_fields$edge.length)
  expect_identical(native_fields$Nnode, utf8_fields$Nnode)
  expect_identical(native_fields$fingerprint, utf8_fields$fingerprint)
  expect_identical(
    vapply(native_fields$tip.label, .v2_prod_raw_hex, character(1L)),
    vapply(utf8_fields$tip.label, .v2_prod_raw_hex, character(1L))
  )
  expect_identical(tree_utf8$tip.label, pair$utf8)
  expect_identical(tree_native$tip.label, pair$native)

  duplicate_labels <- pair$utf8
  duplicate_labels[[2L]] <- pair$native[[1L]]
  duplicate_tree <- .v2_prod_make_tree(duplicate_labels)
  expect_error(
    fastphylosig:::.safe_canonicalize_core(duplicate_tree),
    "UTF-8 byte sequences must be unique"
  )
})

test_that("V2 rejects missing, V1, and unknown prepared context schemas", {
  testthat::skip_if_not_installed("ape")
  tree <- .v2_prod_make_tree()
  trait <- stats::setNames(c(1, 2, 3, 4), tree$tip.label)
  context <- prepare_tree(tree)
  variants <- list(
    missing = context,
    v1 = context,
    unknown = context
  )
  variants$missing$context_schema_version <- NULL
  variants$missing$canonical_contract_version <- NULL
  variants$missing$protected_snapshot_version <- NULL
  variants$v1$context_schema_version <- 1L
  variants$v1$canonical_contract_version <- 1L
  variants$v1$protected_snapshot_version <- 1L
  variants$unknown$protected_snapshot_version <- 999L
  for (name in names(variants)) {
    expect_error(
      fast_k(variants[[name]], trait, test = FALSE,
             verbose = FALSE, progress = FALSE),
      "prepare_tree",
      info = name
    )
  }
})

test_that("V2 protected snapshot rejects mutations before cache or K estimator", {
  testthat::skip_if_not_installed("ape")
  tree <- .v2_prod_make_tree()
  trait <- stats::setNames(c(1, 2, 3, 4), tree$tip.label)
  context <- prepare_tree(tree)
  original_subset <- getFromNamespace(".prepared_tree_subset", "fastphylosig")
  original_kernel <- getFromNamespace("fast_k_tree_batch_cpp", "fastphylosig")
  calls <- new.env(parent = emptyenv())
  calls$subset <- 0L
  calls$kernel <- 0L
  testthat::local_mocked_bindings(
    .prepared_tree_subset = function(...) {
      calls$subset <- calls$subset + 1L
      original_subset(...)
    },
    fast_k_tree_batch_cpp = function(...) {
      calls$kernel <- calls$kernel + 1L
      original_kernel(...)
    },
    .package = "fastphylosig"
  )

  mutations <- list(
    delimiter = function(ctx) {
      ctx$tree$tip.label[[1L]] <- paste0(ctx$tree$tip.label[[1L]], "\rchanged")
      ctx
    },
    edge = function(ctx) {
      ctx$tree$edge[1L, 2L] <- ctx$tree$edge[1L, 2L] + 1L
      ctx
    },
    edge_length = function(ctx) {
      ctx$tree$edge.length[[1L]] <- ctx$tree$edge.length[[1L]] + 0.125
      ctx
    },
    Nnode = function(ctx) {
      ctx$tree$Nnode <- as.integer(ctx$tree$Nnode + 1L)
      ctx
    },
    combined = function(ctx) {
      ctx$tree$tip.label[[1L]] <- paste0(ctx$tree$tip.label[[1L]], "\rchanged")
      ctx$tree$edge[1L, 2L] <- ctx$tree$edge[1L, 2L] + 1L
      ctx$tree$edge.length[[1L]] <- ctx$tree$edge.length[[1L]] + 0.125
      ctx
    }
  )
  for (name in names(mutations)) {
    calls$subset <- 0L
    calls$kernel <- 0L
    mutated <- mutations[[name]](context)
    expect_error(
      fast_k(mutated, trait, test = FALSE,
             verbose = FALSE, progress = FALSE),
      "prepare_tree|modified|inconsistent",
      info = name
    )
    expect_identical(calls$subset, 0L)
    expect_identical(calls$kernel, 0L)
  }
})

test_that("the historical delimiter-collision mutation is rejected before reuse", {
  tree <- .v2_prod_make_tree(c("a", "b\rc", "a\rb", "c"))
  context <- prepare_tree(tree)
  trait <- stats::setNames(c(1, 2, 3, 4), tree$tip.label)
  calls <- new.env(parent = emptyenv())
  calls$subset <- 0L
  calls$kernel <- 0L
  original_subset <- getFromNamespace(".prepared_tree_subset", "fastphylosig")
  original_kernel <- getFromNamespace("fast_k_tree_batch_cpp", "fastphylosig")
  testthat::local_mocked_bindings(
    .prepared_tree_subset = function(...) {
      calls$subset <- calls$subset + 1L
      original_subset(...)
    },
    fast_k_tree_batch_cpp = function(...) {
      calls$kernel <- calls$kernel + 1L
      original_kernel(...)
    },
    .package = "fastphylosig"
  )

  context$tree$tip.label <- c("a\rb", "c", "a", "b\rc")
  expect_error(
    fast_k(context, trait, test = FALSE, verbose = FALSE, progress = FALSE),
    "prepared tree was modified|prepare_tree"
  )
  expect_identical(calls$subset, 0L)
  expect_identical(calls$kernel, 0L)
})

test_that("V2 package-owned integrity record cannot be replaced by public evidence", {
  tree <- .v2_prod_make_tree(c("a", "b\rc", "a\rb", "c"))
  context <- prepare_tree(tree)
  trait <- stats::setNames(c(1, 2, 3, 4), tree$tip.label)
  record <- attr(context, "fastphylosig_v2_integrity", exact = TRUE)
  expect_true(is.environment(record))
  expect_true(environmentIsLocked(record))
  expect_true(bindingIsLocked("protected_snapshot", record))
  expect_true(bindingIsLocked("fingerprint", record))

  synchronized <- context
  synchronized$tree$tip.label <- c("a\rb", "c", "a", "b\rc")
  synchronized$protected_snapshot <- getFromNamespace(
    ".tree_protected_snapshot_v2", "fastphylosig"
  )(synchronized$tree)
  synchronized$fingerprint <- getFromNamespace(
    ".tree_fingerprint", "fastphylosig"
  )(synchronized$tree, .canonical = TRUE)
  expect_error(
    fast_k(synchronized, trait, test = FALSE,
           verbose = FALSE, progress = FALSE),
    "integrity record|prepare_tree"
  )

  metadata <- context
  metadata$canonical_mapping <- rev(metadata$canonical_mapping)
  expect_error(
    fast_k(metadata, trait, test = FALSE,
           verbose = FALSE, progress = FALSE),
    "integrity record|prepare_tree"
  )

  replacement_cache <- context
  replacement_cache$structural_cache <- new.env(parent = emptyenv())
  replacement_cache$cache <- replacement_cache$structural_cache
  expect_error(
    fast_k(replacement_cache, trait, test = FALSE,
           verbose = FALSE, progress = FALSE),
    "integrity record|prepare_tree"
  )
})

test_that("V2 canonicalization is input immutable and metrics carry no descendant payload", {
  testthat::skip_if_not_installed("ape")
  tree <- .v2_prod_make_tree(c("a", "b\rc", "a\rb", "c"))
  before <- serialize(tree, NULL)
  canonical <- .v2_prod_canonical_fields(tree)$tree
  expect_identical(serialize(tree, NULL), before)
  context <- prepare_tree(tree)
  expect_identical(serialize(tree, NULL), before)

  metrics <- attr(canonical, "fastphylosig_v2_metrics", exact = TRUE)
  contract <- attr(canonical, "fastphylosig_v2_contract", exact = TRUE)
  expect_true(is.list(metrics))
  expect_identical(metrics$schema_version, 2L)
  expect_identical(metrics$descendant_label_vectors_materialized, 0L)
  expect_identical(metrics$full_edge_scans, 1L)
  expect_lte(metrics$descriptor_entries, length(tree$tip.label) + tree$Nnode)
  expect_false(any(c("descendant_labels", "descendant_keys") %in%
                     names(metrics)))
  expect_true(is.list(contract))
  expect_false(any(c("descendant_labels", "descendant_keys") %in%
                     names(contract)))
  expect_identical(context$context_schema_version, 2L)
  expect_identical(context$canonical_contract_version, 2L)
  expect_identical(context$protected_snapshot_version, 2L)
})
