# Shared analysis preparation -------------------------------------------------
#
# This file contains the representation-only and data-matching layer shared by
# the public diagnostics and the numerical entry points.  It deliberately does
# not allocate VCV matrices, Cholesky factors, or eigendecompositions.

.canonical_signal_name <- function(signal) {
  if (missing(signal) || !is.character(signal) || length(signal) != 1L ||
      is.na(signal)) {
    stop("signal must be exactly one of \"K\", \"lambda\", \"D\", or \"Delta\".",
         call. = FALSE)
  }
  key <- tolower(trimws(as.character(signal)))
  out <- switch(
    key,
    k = "K",
    lambda = "lambda",
    pagel = "lambda",
    `pagel's lambda` = "lambda",
    d = "D",
    phylo.d = "D",
    phylo_d = "D",
    delta = "Delta",
    phylo.delta = "Delta",
    phylo_delta = "Delta",
    NULL
  )
  if (is.null(out)) {
    stop("signal must be exactly one of \"K\", \"lambda\", \"D\", or \"Delta\".",
         call. = FALSE)
  }
  out
}

.issue_table_for_action <- function(check, signal = NULL) {
  issues <- if (is.list(check)) check$issues else check
  if (!is.data.frame(issues) || !nrow(issues)) return(NULL)
  defaults <- list(
    code = rep("", nrow(issues)),
    severity = rep("ERROR", nrow(issues)),
    signal = rep("all", nrow(issues)),
    message = rep("tree is not ready", nrow(issues)),
    check = rep("tree", nrow(issues)),
    problem = rep("tree is not ready", nrow(issues)),
    action = rep("repair the tree, then re-run check_tree()", nrow(issues)),
    auto_fixable = rep(FALSE, nrow(issues))
  )
  for (nm in names(defaults)) {
    if (!nm %in% names(issues)) issues[[nm]] <- defaults[[nm]]
  }
  if (!is.null(signal)) {
    selected <- is.na(issues$signal) | !nzchar(as.character(issues$signal)) |
      tolower(as.character(issues$signal)) %in% c("all", tolower(signal))
    issues <- issues[selected, , drop = FALSE]
  }
  issues
}

# One formatter is used by preparation and tree resolution so callers receive
# the same actionable diagnostic regardless of which boundary found the issue.
.format_actionable_condition <- function(check = NULL, signal = NULL,
                                          prefix = "Analysis preparation failed",
                                          reason = NULL) {
  issues <- .issue_table_for_action(check, signal = signal)
  severe <- if (is.null(issues)) logical() else
    toupper(as.character(issues$severity)) %in% c("ERROR", "USER_ACTION_REQUIRED")
  if (!is.null(issues) && any(severe)) {
    # A representation issue is safe to repair automatically and should not be
    # reported as a blocking condition after canonicalization.
    severe <- severe & !(issues$auto_fixable %in% TRUE &
                           tolower(as.character(issues$check)) == "representation")
  }
  details <- character()
  actions <- character()
  checks <- character()
  codes <- character()
  if (!is.null(issues) && any(severe)) {
    details <- as.character(issues$problem[severe])
    details <- details[nzchar(details) & !is.na(details)]
    if (!length(details)) details <- as.character(issues$message[severe])
    actions <- as.character(issues$action[severe])
    actions <- actions[nzchar(actions) & !is.na(actions)]
    checks <- as.character(issues$check[severe])
    checks <- checks[nzchar(checks) & !is.na(checks)]
    codes <- as.character(issues$code[severe])
    codes <- codes[nzchar(codes) & !is.na(codes)]
  }
  if (!is.null(reason) && length(reason) && nzchar(as.character(reason[[1L]]))) {
    details <- c(as.character(reason[[1L]]), details)
  }
  details <- unique(details[nzchar(details)])
  if (!length(details)) details <- "the selected signal is not ready"
  actions <- unique(actions)
  if (!length(actions)) {
    actions <- "repair the reported input issue and run check_tree() again"
  }
  checks <- unique(checks)
  signal_name <- if (is.null(signal)) NULL else as.character(signal[[1L]])
  signal_label <- if (is.null(signal_name)) "analysis" else signal_name
  header <- if (identical(prefix, "Analysis preparation failed")) {
    paste0(if (is.null(signal_name)) "Analysis" else signal_name,
           " analysis cannot continue.")
  } else {
    paste0(prefix,
           if (is.null(signal_name)) "." else
             paste0(" for signal \"", signal_name, "\"."))
  }
  paste0(
    header,
    "\n\nProblem:\n- ", paste(details, collapse = "\n- "),
    "\n\nWhy it blocks the requested signal:\n",
    "The current ", signal_label, " requirements are not met",
    if (length(checks)) paste0(" (", paste(checks, collapse = ", "), ")") else "",
    ".",
    if (length(codes)) paste0("\nIssue code(s): ", paste(unique(codes), collapse = ", ")) else "",
    "\n\nRecommended action:\n- ", paste(actions, collapse = "\n- "),
    "\n\nNo automatic biological repair was applied. USER_ACTION_REQUIRED.",
    "\n\nFor a full diagnosis:\n",
    if (is.null(signal_name)) "check_tree(tree)" else
      paste0("check_tree(tree, signal = \"", signal_name, "\")")
  )
}

# Short alias retained for internal callers and downstream extensions.
.actionable_condition <- .format_actionable_condition

.canonical_tree_signature <- function(tree) {
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
  root_key <- if (is.finite(root)) paste(descendants(root), collapse = "\r") else NA_character_
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

.safe_canonicalize_core <- function(tree) {
  if (!inherits(tree, "phylo")) {
    stop("tree should be an object of class \"phylo\".", call. = FALSE)
  }
  original <- tree
  before <- .canonical_tree_signature(original)
  info <- list(
    changed = FALSE,
    mapping = NULL,
    tip_identity = FALSE,
    branch_identity = FALSE,
    edge_isomorphism = FALSE,
    cophenetic_identity = FALSE,
    root_identity = FALSE,
    polytomy_identity = FALSE,
    safe = FALSE,
    reason = NULL
  )
  if (is.null(before)) {
    attr(original, "fastphylosig_canonicalization") <- info
    return(original)
  }
  edge <- before$edge
  n_tip <- length(before$tip_label)
  n_node <- tree$Nnode
  if (length(n_node) != 1L || !is.finite(n_node) || n_node < 1L ||
      n_node != floor(n_node)) {
    info$reason <- "Nnode is not a finite positive integer"
    attr(original, "fastphylosig_canonicalization") <- info
    return(original)
  }
  internal <- sort(unique(c(edge[, 1L], edge[, 2L])[c(edge[, 1L], edge[, 2L]) > n_tip]))
  roots <- setdiff(unique(edge[, 1L]), unique(edge[, 2L]))
  if (!length(internal) || length(internal) != as.integer(n_node) ||
      length(roots) != 1L) {
    info$reason <- "edge topology does not have one structural root and contiguous internals"
    attr(original, "fastphylosig_canonicalization") <- info
    return(original)
  }
  root <- roots[[1L]]
  parent <- edge[, 1L]
  child <- edge[, 2L]
  children_of <- function(node) child[parent == node]
  memo <- new.env(parent = emptyenv())
  descendants <- function(node, active = integer()) {
    key <- as.character(node)
    if (exists(key, memo, inherits = FALSE)) return(get(key, memo, inherits = FALSE))
    if (node <= n_tip) {
      ans <- if (length(before$tip_label) >= node) before$tip_label[[node]] else paste0("#", node)
      assign(key, ans, memo)
      return(ans)
    }
    if (node %in% active) return(paste0("!cycle:", node))
    kids <- children_of(node)
    ans <- if (!length(kids)) paste0("!empty:", node) else
      sort(unlist(lapply(kids, descendants, active = c(active, node)), use.names = FALSE))
    assign(key, ans, memo)
    ans
  }
  remaining <- setdiff(internal, root)
  if (length(remaining)) {
    keys <- vapply(remaining, function(z) paste(descendants(z), collapse = "\r"), character(1L))
    remaining <- remaining[order(keys, remaining)]
  }
  old_order <- c(root, remaining)
  new_order <- seq.int(n_tip + 1L, length.out = length(old_order))
  map <- stats::setNames(new_order, as.character(old_order))
  mapped_edge <- cbind(
    unname(map[as.character(parent)]),
    ifelse(child <= n_tip, child, unname(map[as.character(child)]))
  )
  storage.mode(mapped_edge) <- "integer"
  out <- original
  out$edge <- mapped_edge
  if (!is.null(out$node.label) && length(out$node.label)) {
    old_labels <- as.character(out$node.label)
    labels <- rep(NA_character_, length(old_order))
    for (i in seq_along(old_order)) {
      j <- old_order[[i]] - n_tip
      if (j >= 1L && j <= length(old_labels)) labels[[i]] <- old_labels[[j]]
    }
    out$node.label <- labels
  }
  out <- tryCatch(ape::reorder.phylo(out, order = "postorder"),
                  error = function(e) out)
  after <- .canonical_tree_signature(out)
  info$mapping <- list(
    old_to_new = map,
    new_to_old = stats::setNames(old_order, as.character(new_order)),
    old_internal_order = old_order,
    new_internal_order = new_order
  )
  info$tip_identity <- !is.null(after) && identical(before$tip_label, after$tip_label)
  info$branch_identity <- !is.null(after) &&
    if (is.null(before$edge_length) || is.null(after$edge_length)) {
      is.null(before$edge_length) && is.null(after$edge_length)
    } else {
      isTRUE(all.equal(sort(before$edge_length), sort(after$edge_length),
                       tolerance = 0, check.attributes = FALSE))
    }
  info$root_identity <- !is.null(after) &&
    identical(before$root_descendants, after$root_descendants) &&
    identical(before$root_degree, after$root_degree)
  info$polytomy_identity <- !is.null(after) &&
    identical(before$polytomies, after$polytomies)
  candidate_changed <- !identical(unname(as.matrix(original$edge)),
                                  unname(as.matrix(out$edge))) ||
    !identical(as.character(original$node.label), as.character(out$node.label))
  edge_records <- function(edge_matrix, edge_length) {
    length_key <- if (is.null(edge_length)) {
      rep("<none>", nrow(edge_matrix))
    } else {
      sprintf("%.17g", as.numeric(edge_length))
    }
    sort(paste(edge_matrix[, 1L], edge_matrix[, 2L], length_key,
               sep = "\r"))
  }
  # Compare every mapped parent-child-length record after ape's reorder. This
  # is an exact weighted-graph isomorphism check and therefore proves all tip
  # patristic distances are preserved without constructing an O(n_tip^2)
  # cophenetic matrix or asking ape to parse the noncanonical input.
  info$edge_isomorphism <- !is.null(after) && identical(
    edge_records(mapped_edge, original$edge.length),
    edge_records(after$edge, out$edge.length)
  )
  info$cophenetic_identity <- isTRUE(info$edge_isomorphism)
  info$safe <- isTRUE(info$tip_identity) && isTRUE(info$branch_identity) &&
    isTRUE(info$cophenetic_identity) && isTRUE(info$root_identity) &&
    isTRUE(info$polytomy_identity)
  if (!isTRUE(info$safe)) {
    info$reason <- "canonicalization would alter tip identities, branch lengths, root, polytomies, or cophenetic distances"
    out <- original
  }
  info$changed <- isTRUE(info$safe) &&
    !identical(unname(as.matrix(original$edge)), unname(as.matrix(out$edge)))
  attr(out, "fastphylosig_canonicalization") <- info
  out
}

.canonicalization_info <- function(tree) {
  info <- attr(tree, "fastphylosig_canonicalization", exact = TRUE)
  if (is.list(info)) info else list(changed = FALSE, safe = FALSE)
}

.as_named_trait_table <- function(x, tree, verbose = TRUE,
                                  input_name = "x", vector_name = "x") {
  input_names <- if (is.null(dim(x))) names(x) else rownames(x)
  if (!is.null(input_names) && anyDuplicated(input_names)) {
    stop(input_name, " species names must be unique.", call. = FALSE)
  }

  if (is.data.frame(x)) {
    rn <- rownames(x)
    out <- as.data.frame(x, stringsAsFactors = FALSE)
    rownames(out) <- rn
  } else if (is.null(dim(x))) {
    out <- data.frame(value = I(x), stringsAsFactors = FALSE)
    names(out) <- vector_name
    rownames(out) <- names(x)
  } else {
    out <- as.data.frame(x, stringsAsFactors = FALSE)
  }

  rn <- rownames(out)
  default_rownames <- identical(rn, as.character(seq_len(nrow(out))))
  if ((is.null(rn) || default_rownames) &&
      nrow(out) == length(tree$tip.label)) {
    if (isTRUE(verbose)) {
      message(input_name, " has no species names; assuming tree tip order")
    }
    rownames(out) <- tree$tip.label
    rn <- rownames(out)
    default_rownames <- FALSE
  }
  if (is.null(rn) || default_rownames || anyNA(rn) || any(!nzchar(rn))) {
    stop(input_name, " must have non-empty species names as names/rownames.",
         call. = FALSE)
  }
  if (anyDuplicated(rn)) {
    stop(input_name, " species names must be unique.", call. = FALSE)
  }
  if (is.null(colnames(out))) {
    colnames(out) <- paste0("trait_", seq_len(ncol(out)))
  }
  out
}

.analysis_data_table <- function(data, tree, data_kind = NULL) {
  if (is.null(data)) stop("data must be supplied.", call. = FALSE)
  out <- .as_named_trait_table(
    data, tree, verbose = FALSE, input_name = "data", vector_name = "value"
  )
  if (!is.null(data_kind)) {
    key <- tolower(as.character(data_kind[[1L]]))
    if (!key %in% c("raw", "prepared", "continuous", "numeric", "trait",
                    "categorical", "discrete", "binary")) {
      stop("data_kind must be raw/prepared or a supported trait kind.", call. = FALSE)
    }
  }
  out
}

.analysis_na_groups <- function(present) {
  if (!is.matrix(present)) present <- as.matrix(present)
  if (exists("group_na_masks_cpp", mode = "function", inherits = TRUE)) {
    ans <- tryCatch(group_na_masks_cpp(present), error = function(e) NULL)
    if (is.list(ans) && !is.null(ans$keep) && !is.null(ans$columns)) return(ans)
  }
  if (!ncol(present)) return(list(n_group = 0L, columns = list(), keep = list(), key = character()))
  keys <- apply(present, 2L, function(z) paste(as.integer(z), collapse = ""))
  unique_keys <- unique(keys)
  list(
    n_group = length(unique_keys),
    columns = lapply(unique_keys, function(k) which(keys == k)),
    keep = lapply(unique_keys, function(k) which(present[, which(keys == k)[[1L]], drop = TRUE])),
    key = unique_keys
  )
}

.prepare_analysis <- function(tree, data, signal, data_kind,
                              verbose = FALSE) {
  signal <- .canonical_signal_name(signal)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }
  prepared_input <- inherits(tree, "fastphylosig_tree")
  tree_processing <- list(
    input = if (prepared_input) "prepared" else "raw",
    canonicalized = FALSE,
    canonical_mapping = NULL,
    original_fingerprint = NULL,
    final_fingerprint = NULL,
    check_before = NULL,
    check_after = NULL
  )
  if (prepared_input) {
    .validate_prepared_context(tree)
    ctx <- tree
    working_tree <- ctx$tree
    tree_processing$original_fingerprint <- ctx$fingerprint
    checked <- .inspect_tree_core(ctx, signal = signal)
    tree_processing$check_before <- checked
    tree_processing$check_after <- checked
    # A prepared context is generic and may legitimately cache a tree that is
    # only ready for a subset of methods (for example, zero internal branches
    # are a lambda boundary but block K/D/Delta).  Enforce the selected
    # method's readiness here, before matching, subset preparation, or a C++
    # kernel.  This keeps prepared and raw inputs on the same issue-registry
    # path without changing any numerical calculation.
    if (!isTRUE(checked$ready_by_signal[[signal]])) {
      stop(.format_actionable_condition(checked, signal), call. = FALSE)
    }
  } else {
    if (!inherits(tree, "phylo")) {
      stop(.format_actionable_condition(NULL, signal,
                                        reason = "tree must be a phylo object"),
           call. = FALSE)
    }
    tree_processing$original_fingerprint <- tryCatch(.tree_fingerprint(tree), error = function(e) NULL)
    before <- .inspect_tree_core(tree, signal = signal)
    tree_processing$check_before <- before
    working_tree <- .safe_canonicalize_core(tree)
    info <- .canonicalization_info(working_tree)
    tree_processing$canonicalized <- isTRUE(info$changed)
    tree_processing$canonical_mapping <- info$mapping
    # A failed invariant check returns an unchanged copy.  It is safe to
    # continue only when the selected method was already ready.
    if (!isTRUE(info$safe) && !isTRUE(before$ready_by_signal[[signal]])) {
      stop(.format_actionable_condition(before, signal), call. = FALSE)
    }
    checked <- .inspect_tree_core(working_tree, signal = signal)
    tree_processing$check_after <- checked
    if (!isTRUE(checked$ready_by_signal[[signal]])) {
      stop(.format_actionable_condition(checked, signal), call. = FALSE)
    }
    ctx <- prepare_tree(working_tree)
  }
  tree_processing$final_fingerprint <- ctx$fingerprint

  table <- .analysis_data_table(data, ctx$tree, data_kind = data_kind)
  matched <- .match_tree_data_core(ctx, data = table, prune = TRUE,
                                  verbose = verbose, allow_insufficient = TRUE)
  base_keep <- matched$base_keep
  retained <- length(base_keep)
  matching <- matched$report
  matching$retained_species <- retained
  matching$insufficient_retained <- retained < 2L
  matching$kernel_ready <- retained >= 2L
  base_group <- NULL
  if (retained >= 2L) {
    base_group <- .prepared_tree_subset(ctx, base_keep, need_matrix = FALSE)
  }

  matched_data <- matched$data
  present <- if (ncol(matched_data)) !is.na(as.matrix(matched_data)) else
    matrix(logical(), nrow = nrow(matched_data), ncol = 0L)
  na_patterns <- .analysis_na_groups(present)
  groups <- vector("list", na_patterns$n_group)
  checked_subsets <- new.env(parent = emptyenv())
  check_subset <- function(keep, packed_key = NULL) {
    keep <- sort(unique(as.integer(keep)))
    # Do not use a comma-separated tip-index list as an environment key: for
    # large trees it exceeds R's 10,000-byte symbol limit.  The packed key
    # generated by group_na_masks_cpp is compact, deterministic, and includes
    # the mask dimension.  Keep a short sentinel for an empty retained set.
    key <- if (!length(keep)) {
      "<none>"
    } else if (length(packed_key) == 1L && is.character(packed_key) &&
               nzchar(packed_key)) {
      packed_key
    } else {
      .tree_mask_key(keep, nrow(matched_data))
    }
    if (exists(key, checked_subsets, inherits = FALSE)) return(get(key, checked_subsets, inherits = FALSE))
    result <- if (length(keep) < 2L) {
      list(ready = FALSE, ready_by_signal = stats::setNames(FALSE, signal),
           issues = data.frame(code = "too_few_tips", severity = "ERROR",
                               signal = signal,
                               message = "fewer than two retained species",
                               check = "retained subset", problem = "fewer than two retained species",
                               action = "retain at least two species, then run check_tree()",
                               auto_fixable = FALSE, stringsAsFactors = FALSE))
    } else {
      .inspect_tree_core(.prepared_tree_subset(ctx, base_keep[keep], need_matrix = FALSE)$tree,
                         signal = signal)
    }
    assign(key, result, checked_subsets)
    result
  }
  if (na_patterns$n_group) {
    for (i in seq_len(na_patterns$n_group)) {
      cols <- as.integer(na_patterns$columns[[i]])
      keep <- as.integer(na_patterns$keep[[i]])
      packed_key <- if (length(na_patterns$key) >= i) {
        as.character(na_patterns$key[[i]])
      } else {
        NULL
      }
      check <- check_subset(keep, packed_key = packed_key)
      groups[[i]] <- list(
        index = i,
        columns = cols,
        keep = keep,
        retained = length(keep),
        data = if (length(cols)) matched_data[keep, cols, drop = FALSE] else matched_data[keep, , drop = FALSE],
        group = if (length(keep) >= 2L) .prepared_tree_subset(ctx, base_keep[keep], need_matrix = FALSE) else NULL,
        check = check,
        kernel_ready = length(keep) >= 2L && isTRUE(check$ready_by_signal[[signal]])
      )
    }
  }
  out <- list(
    ctx = ctx,
    matched_data = matched_data,
    base_group = base_group,
    groups = groups,
    matching = matching,
    matching_details = list(
      matched_species = matched$matched_species,
      tree_tips_removed = matched$tree_tips_removed,
      data_rows_removed = matched$data_rows_removed
    ),
    tree_processing = tree_processing,
    na_patterns = na_patterns,
    signal = signal,
    data_kind = data_kind,
    kernel_ready = retained >= 2L && isTRUE(checked$ready_by_signal[[signal]]),
    insufficient_retained = retained < 2L
  )
  class(out) <- c("fastphylosig_analysis_preparation", "list")
  out
}

.compact_analysis_metadata <- function(analysis) {
  if (!inherits(analysis, "fastphylosig_analysis_preparation")) return(NULL)
  patterns <- analysis$na_patterns
  list(
    signal = analysis$signal,
    input_tree_type = analysis$tree_processing$input,
    raw_or_prepared = analysis$tree_processing$input,
    tree_auto_normalized = isTRUE(analysis$tree_processing$canonicalized),
    canonical_changes = if (isTRUE(analysis$tree_processing$canonicalized)) {
      "internal node numbering and edge ordering"
    } else character(),
    matching = analysis$matching,
    matched_species = analysis$matching_details$matched_species,
    tree_tips_removed = analysis$matching_details$tree_tips_removed,
    data_rows_removed = analysis$matching_details$data_rows_removed,
    na_pattern_count = if (is.null(patterns$n_group)) 0L else
      as.integer(patterns$n_group),
    na_pattern_keys = if (is.null(patterns$key)) character() else
      as.character(patterns$key),
    retained_tip_validation = vapply(
      analysis$groups,
      function(z) isTRUE(z$kernel_ready),
      logical(1)
    )
  )
}
