# Shared tree preparation and cache ------------------------------------------

#' Prepare a phylogenetic tree for repeated signal calculations
#'
#' The returned object stores tree-only data and an explicit, invalidation-safe
#' cache. K, lambda, D, Delta, and fast_ace() accept this object anywhere a
#' phylo object is accepted. Dense numerical resources are prepared lazily, so
#' D/Delta/ACE-only workflows do not allocate a VCV matrix or factorization.
#'
#' @param tree A rooted phylo object with finite, non-negative branch lengths.
#' @param max_cached_subsets Maximum number of full or NA-pruned subset entries
#' to retain in the structural cache. Retained for compatibility; dense memory
#' is controlled independently by `cache_budget`.
#' @param cache_budget Maximum bytes retained in the numerical LRU cache.
#' @return An object of class `fastphylosig_tree`.
#' @export
prepare_tree <- function(tree, max_cached_subsets = 16L,
                         cache_budget = 512 * 1024^2) {
  if (inherits(tree, "fastphylosig_tree")) {
    .validate_prepared_context(tree)
    return(tree)
  }
  if (!inherits(tree, "phylo")) {
    stop("tree should be an object of class \"phylo\".", call. = FALSE)
  }
  max_cached_subsets <- .validate_cache_limit(max_cached_subsets)
  cache_budget <- .validate_cache_budget(cache_budget)

  # Cache representation diagnostics and method capabilities alongside the
  # structural context.  This pass is read-only and intentionally does not
  # allocate a VCV matrix, Cholesky factor, or eigendecomposition.
  generic_summary <- tryCatch(
    .inspect_tree_core(tree, signal = c("K", "lambda", "D", "Delta")),
    error = function(e) NULL
  )
  if (!is.list(generic_summary) ||
      !isTRUE(generic_summary$tree_summary$valid)) {
    generic_check <- generic_summary
    if (is.list(generic_check) && is.data.frame(generic_check$issues)) {
      generic_check$issues <- generic_check$issues[
        generic_check$issues$signal == "all", , drop = FALSE
      ]
    }
    stop(.format_actionable_condition(
      generic_check, signal = NULL,
      prefix = "prepare_tree() cannot compile this tree"
    ), call. = FALSE)
  }
  .validate_prepare_tree(tree)
  canonical_probe <- tryCatch(.safe_canonicalize_core(tree),
                              error = function(e) tree)
  canonical_info <- tryCatch(.canonicalization_info(canonical_probe),
                             error = function(e) list(changed = FALSE, safe = FALSE))
  structural_entry_validation <- list(
    valid = if (is.list(generic_summary)) isTRUE(generic_summary$tree_summary$valid) else FALSE,
    n_tip = ape::Ntip(tree),
    n_node = as.integer(tree$Nnode),
    edge_rows = nrow(tree$edge),
    branch_length_rows = length(tree$edge.length),
    canonicalizable = isTRUE(canonical_info$safe),
    issue = if (!is.null(canonical_info$reason)) canonical_info$reason else NULL
  )

  structural_cache <- new.env(parent = emptyenv())
  numerical_cache <- new.env(parent = emptyenv())
  cache_meta <- new.env(parent = emptyenv())
  cache_meta$clock <- 1
  cache_meta$hits <- 0L
  cache_meta$misses <- 0L
  cache_meta$evictions <- 0L
  cache_meta$bytes_used <- 0
  full <- .prepare_tree_subset(
    tree, need_lambda = FALSE, need_matrix = FALSE
  )
  full$.last_access <- cache_meta$clock
  assign(
    .tree_mask_key(seq_len(ape::Ntip(tree)), ape::Ntip(tree)),
    full, structural_cache
  )

  out <- list(
    tree = tree,
    tip.label = tree$tip.label,
    n_tip = ape::Ntip(tree),
    cache = structural_cache,
    structural_cache = structural_cache,
    numerical_cache = numerical_cache,
    cache_meta = cache_meta,
    cache_budget = cache_budget,
    max_cached_subsets = max_cached_subsets,
    fingerprint = .tree_fingerprint(tree),
    # These fields are metadata only; dense numerical resources remain lazy.
    inspection = generic_summary,
    canonical_mapping = canonical_info$mapping,
    canonical_summary = if (is.list(generic_summary)) generic_summary$tree_summary else NULL,
    generic_summary = if (is.list(generic_summary)) generic_summary$tree_summary else NULL,
    method_capability = if (is.list(generic_summary)) generic_summary$ready_by_signal else NULL,
    structural_entry_validation = structural_entry_validation
  )
  class(out) <- c("fastphylosig_tree", "list")
  out
}

.validate_prepare_tree <- function(tree) {
  if (is.null(tree$tip.label) || anyNA(tree$tip.label) ||
      anyDuplicated(tree$tip.label)) {
    stop("tree tip labels must be present, non-missing, and unique.",
         call. = FALSE)
  }
  if (is.null(tree$edge.length)) {
    stop("tree must contain branch lengths.", call. = FALSE)
  }
  edge_length <- as.numeric(tree$edge.length)
  if (length(edge_length) != nrow(tree$edge) ||
      any(!is.finite(edge_length)) || any(edge_length < 0)) {
    stop("tree branch lengths must be finite and non-negative.",
         call. = FALSE)
  }
  if (ape::Ntip(tree) < 2L) {
    stop("tree must contain at least two tips.", call. = FALSE)
  }
  invisible(tree)
}

.tree_fingerprint <- function(tree) {
  paste(
    paste(tree$tip.label, collapse = "\r"),
    paste(as.integer(tree$edge), collapse = ","),
    paste(formatC(as.numeric(tree$edge.length), digits = 17,
                  format = "fg"), collapse = ","),
    paste(as.integer(tree$Nnode), collapse = ","),
    sep = "|"
  )
}

.tree_mask_key <- function(keep, n_tip = NULL) {
  keep <- as.integer(keep)
  if (is.null(n_tip)) n_tip <- if (length(keep)) max(keep) else 0L
  n_tip <- as.integer(n_tip)
  if (length(n_tip) != 1L || is.na(n_tip) || n_tip < 1L ||
      anyNA(keep) || any(keep < 1L) || any(keep > n_tip)) {
    stop("keep and n_tip cannot form a valid tree mask.", call. = FALSE)
  }
  present <- matrix(FALSE, nrow = n_tip, ncol = 1L)
  present[keep, 1L] <- TRUE
  # Reuse the exact packed-mask encoder used to group trait NA patterns. This
  # keeps cache keys far below R's 10,000-byte environment-name limit for the
  # large trees supported by the package and is independent of keep order.
  group_na_masks_cpp(present)$key[[1L]]
}

.prepare_tree_subset <- function(tree, need_lambda = FALSE,
                                 need_matrix = TRUE) {
  need_matrix <- isTRUE(need_matrix) || isTRUE(need_lambda)
  C <- if (need_matrix) ape::vcv.phylo(tree) else NULL
  cholC <- if (need_matrix) tryCatch(chol(C), error = function(e) NULL) else NULL
  d_tree <- ape::reorder.phylo(tree, "pruningwise")
  edge <- matrix(as.integer(d_tree$edge), ncol = 2L)
  edge_length <- as.numeric(d_tree$edge.length)
  compiled_tree <- compile_tree_cpp(
    edge = edge, edge_length = edge_length, n_tip = ape::Ntip(tree)
  )
  out <- list(
    tree = tree,
    d_tree = d_tree,
    compiled_tree = compiled_tree,
    C = C,
    chol = cholC,
    traceC = if (need_matrix) sum(diag(C)) else NA_real_,
    edge = edge,
    edge_length = edge_length,
    n_tip = ape::Ntip(tree),
    lambda_spectral = NULL
  )
  if (isTRUE(need_lambda)) {
    out$lambda_spectral <- .lambda_spectral_cache(C)
  }
  out
}

.lambda_spectral_cache <- function(C) {
  d <- diag(C)
  if (any(!is.finite(d)) || any(d <= 0)) return(NULL)
  scale <- sqrt(outer(d, d))
  R <- C / scale
  eig <- tryCatch(eigen(R, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig) || any(!is.finite(eig$values)) ||
      any(!is.finite(eig$vectors))) {
    return(NULL)
  }
  list(
    values = as.numeric(eig$values),
    vectors = eig$vectors,
    inv_sqrt_diag = 1 / sqrt(d),
    log_diag = sum(log(d))
  )
}

.prepared_tree_subset <- function(ctx, keep = NULL, need_lambda = FALSE,
                                  need_matrix = TRUE) {
  if (!inherits(ctx, "fastphylosig_tree")) {
    ctx <- prepare_tree(ctx)
  } else {
    .validate_prepared_context(ctx)
  }
  if (is.null(keep)) keep <- seq_len(ctx$n_tip)
  keep <- as.integer(keep)
  if (length(keep) < 2L || any(keep < 1L) || any(keep > ctx$n_tip) ||
      anyDuplicated(keep)) {
    stop("keep must contain at least two unique valid tip indices.",
         call. = FALSE)
  }
  key <- .tree_mask_key(keep, ctx$n_tip)
  need_matrix <- isTRUE(need_matrix) || isTRUE(need_lambda)
  structural_cache <- .structural_cache(ctx)
  if (exists(key, structural_cache, inherits = FALSE)) {
    out <- get(key, structural_cache, inherits = FALSE)
    .validate_structural_entry(out, ctx)
    .cache_record(ctx, "hits")
    out$.last_access <- .cache_tick(ctx)
    assign(key, out, structural_cache)
  } else {
    .cache_record(ctx, "misses")
    drop <- setdiff(seq_len(ctx$n_tip), keep)
    subset_tree <- if (length(drop)) {
      ape::drop.tip(ctx$tree, ctx$tree$tip.label[drop])
    } else {
      ctx$tree
    }
    out <- .prepare_tree_subset(
      subset_tree, need_lambda = FALSE, need_matrix = FALSE
    )
    .validate_structural_entry(out, ctx)
    out$.last_access <- .cache_tick(ctx)
    max_entries <- if (is.null(ctx$max_cached_subsets)) 16L else
      ctx$max_cached_subsets
    if (length(ls(structural_cache, all.names = TRUE)) < max_entries) {
      assign(key, out, structural_cache)
    }
  }
  if (need_matrix) {
    payload <- .numerical_payload(ctx, key, out, need_lambda = need_lambda)
    out[names(payload)] <- payload
  }
  out
}

.validate_structural_entry <- function(entry, ctx = NULL) {
  valid <- is.list(entry) && inherits(entry$tree, "phylo") &&
    is.list(entry$compiled_tree) && !is.null(entry$edge) &&
    is.matrix(entry$edge) && ncol(entry$edge) == 2L &&
    !is.null(entry$edge_length) && length(entry$edge_length) == nrow(entry$edge)
  if (!isTRUE(valid)) {
    stop("invalid structural cache entry; call prepare_tree(tree) again.",
         call. = FALSE)
  }
  invisible(entry)
}

.structural_cache <- function(ctx) {
  if (!is.null(ctx$structural_cache)) ctx$structural_cache else ctx$cache
}

.numerical_payload <- function(ctx, key, structural, need_lambda = FALSE) {
  cache <- ctx$numerical_cache
  if (is.null(cache)) {
    payload <- .prepare_numerical_payload(structural$tree, need_lambda)
    return(payload)
  }
  if (exists(key, cache, inherits = FALSE)) {
    payload <- get(key, cache, inherits = FALSE)
    .cache_record(ctx, "hits")
    if (isTRUE(need_lambda) && is.null(payload$lambda_spectral)) {
      payload$lambda_spectral <- .lambda_spectral_cache(payload$C)
      payload <- .cache_store_payload(ctx, key, payload)
    } else {
      payload$.last_access <- .cache_tick(ctx)
      assign(key, payload, cache)
    }
    return(payload)
  }
  .cache_record(ctx, "misses")
  payload <- .prepare_numerical_payload(structural$tree, need_lambda)
  .cache_store_payload(ctx, key, payload)
}

.prepare_numerical_payload <- function(tree, need_lambda = FALSE) {
  C <- ape::vcv.phylo(tree)
  list(
    .n_tip = ape::Ntip(tree),
    C = C,
    chol = tryCatch(chol(C), error = function(e) NULL),
    traceC = sum(diag(C)),
    lambda_spectral = if (isTRUE(need_lambda)) {
      .lambda_spectral_cache(C)
    } else {
      NULL
    }
  )
}

.payload_bytes <- function(payload) {
  keep <- intersect(c("C", "chol", "lambda_spectral"), names(payload))
  as.numeric(utils::object.size(payload[keep]))
}

.cache_store_payload <- function(ctx, key, payload) {
  cache <- ctx$numerical_cache
  if (exists(key, cache, inherits = FALSE)) {
    old <- get(key, cache, inherits = FALSE)
    ctx$cache_meta$bytes_used <- max(
      0, ctx$cache_meta$bytes_used - as.numeric(old$.bytes)
    )
    rm(list = key, envir = cache)
  }
  payload$.bytes <- .payload_bytes(payload)
  payload$.last_access <- .cache_tick(ctx)
  budget <- ctx$cache_budget
  if (is.null(budget)) budget <- Inf
  if (payload$.bytes > budget) return(payload)

  while (ctx$cache_meta$bytes_used + payload$.bytes > budget) {
    keys <- ls(cache, all.names = TRUE)
    if (!length(keys)) break
    last <- vapply(keys, function(z) {
      as.numeric(get(z, cache, inherits = FALSE)$.last_access)
    }, numeric(1))
    victim <- keys[[which.min(last)]]
    old <- get(victim, cache, inherits = FALSE)
    ctx$cache_meta$bytes_used <- max(
      0, ctx$cache_meta$bytes_used - as.numeric(old$.bytes)
    )
    rm(list = victim, envir = cache)
    .cache_record(ctx, "evictions")
  }
  assign(key, payload, cache)
  ctx$cache_meta$bytes_used <- ctx$cache_meta$bytes_used + payload$.bytes
  payload
}

.cache_tick <- function(ctx) {
  if (is.null(ctx$cache_meta)) return(as.numeric(Sys.time()))
  ctx$cache_meta$clock <- ctx$cache_meta$clock + 1
  ctx$cache_meta$clock
}

.cache_record <- function(ctx, field) {
  if (!is.null(ctx$cache_meta)) {
    ctx$cache_meta[[field]] <- ctx$cache_meta[[field]] + 1L
  }
  invisible(NULL)
}

.validate_cache_limit <- function(value) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value != floor(value) ||
      value > .Machine$integer.max) {
    stop("max_cached_subsets must be a positive integer.", call. = FALSE)
  }
  as.integer(value)
}

.validate_cache_budget <- function(value) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0) {
    stop("cache_budget must be one finite non-negative byte value.",
         call. = FALSE)
  }
  as.numeric(value)
}

#' Inspect a prepared-tree cache
#'
#' @param ctx An object returned by `prepare_tree()`.
#' @return A list containing budget, current numerical bytes, cache counters,
#' and one row per structural entry.
#' @export
cache_info <- function(ctx) {
  .validate_prepared_context(ctx)
  structural <- .structural_cache(ctx)
  numerical <- ctx$numerical_cache
  skeys <- ls(structural, all.names = TRUE)
  nkeys <- if (is.null(numerical)) character() else
    ls(numerical, all.names = TRUE)
  keys <- unique(c(skeys, nkeys))
  entries <- if (!length(keys)) {
    data.frame(
      entry = character(), n_tip = integer(), bytes = numeric(),
      structural_bytes = numeric(), numerical_bytes = numeric(),
      has_vcv = logical(), has_chol = logical(), has_eigen = logical(),
      last_access = numeric(), stringsAsFactors = FALSE
    )
  } else {
    table <- do.call(rbind, lapply(keys, function(key) {
      s <- if (exists(key, structural, inherits = FALSE)) {
        get(key, structural, inherits = FALSE)
      } else NULL
      p <- if (!is.null(numerical) &&
               exists(key, numerical, inherits = FALSE)) {
        get(key, numerical, inherits = FALSE)
      } else NULL
      structural_bytes <- if (is.null(s)) 0 else
        as.numeric(utils::object.size(s$compiled_tree))
      numerical_bytes <- if (is.null(p)) 0 else as.numeric(p$.bytes)
      data.frame(
        entry = key,
        n_tip = if (!is.null(s)) as.integer(s$n_tip) else
          as.integer(p$.n_tip),
        bytes = structural_bytes + numerical_bytes,
        structural_bytes = structural_bytes,
        numerical_bytes = numerical_bytes,
        has_vcv = !is.null(p$C),
        has_chol = !is.null(p$chol),
        has_eigen = !is.null(p$lambda_spectral),
        last_access = max(
          c(if (!is.null(s)) as.numeric(s$.last_access) else NA_real_,
            if (!is.null(p)) as.numeric(p$.last_access) else NA_real_),
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }))
    rownames(table) <- NULL
    table
  }
  meta <- ctx$cache_meta
  list(
    budget = if (is.null(ctx$cache_budget)) Inf else ctx$cache_budget,
    bytes_used = if (is.null(meta)) NA_real_ else meta$bytes_used,
    n_structural_entries = length(skeys),
    n_numerical_entries = length(nkeys),
    hits = if (is.null(meta)) NA_integer_ else meta$hits,
    misses = if (is.null(meta)) NA_integer_ else meta$misses,
    evictions = if (is.null(meta)) NA_integer_ else meta$evictions,
    entries = entries
  )
}

.as_prepared_tree <- function(tree) {
  if (inherits(tree, "fastphylosig_tree")) {
    .validate_prepared_context(tree)
    return(tree)
  }
  prepare_tree(tree)
}

.validate_prepared_context <- function(ctx) {
  if (!is.list(ctx) || !inherits(ctx, "fastphylosig_tree") ||
      is.null(ctx$tree) || is.null(ctx$fingerprint)) {
    stop("invalid fastphylosig_tree context; call prepare_tree(tree) again.",
         call. = FALSE)
  }
  current <- .tree_fingerprint(ctx$tree)
  if (!identical(current, ctx$fingerprint)) {
    stop(
      "the prepared tree was modified after caching; call prepare_tree(tree) again.",
      call. = FALSE
    )
  }
  invisible(ctx)
}
