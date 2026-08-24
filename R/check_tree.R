# Read-only tree preflight ----------------------------------------------------

.tree_check_select <- function(check, signal, prepared = FALSE) {
  if (!is.list(check) || is.null(check$ready_by_signal) ||
      is.null(check$tree_summary) || is.null(check$issues)) {
    stop("invalid cached tree inspection.", call. = FALSE)
  }
  signal <- unique(as.character(signal))
  if (!length(signal) || any(!signal %in% names(check$ready_by_signal))) {
    stop("cached tree inspection does not contain the requested signal.",
         call. = FALSE)
  }
  ready_by_signal <- check$ready_by_signal[signal]
  issues <- check$issues
  if (nrow(issues)) {
    issues <- issues[issues$signal %in% c("all", signal), , drop = FALSE]
    rownames(issues) <- NULL
  }
  tree_summary <- check$tree_summary
  tree_summary$prepared <- isTRUE(prepared)
  out <- list(
    ready = all(ready_by_signal),
    ready_by_signal = ready_by_signal,
    tree_summary = tree_summary,
    issues = issues
  )
  class(out) <- c("fastphylosig_tree_check", "list")
  out
}

#' Check a phylogenetic tree before signal calculations
#'
#' `check_tree()` inspects a `phylo` object (or a prepared
#' `fastphylosig_tree`) without reordering, pruning, compiling, or otherwise
#' changing it.  Readiness is reported separately for each selected signal.
#' K, lambda, and D use the tree-pruning representations and can retain
#' polytomies; Delta follows `fast_ace()` and therefore needs a rooted,
#' fully-dichotomous tree with strictly positive branch lengths.
#'
#' @param tree A `phylo` object or an object returned by [prepare_tree()].
#' @param signal Character vector selecting one or more signal methods.  The
#'   default checks all four public methods (`"K"`, `"lambda"`, `"D"`, and
#'   `"Delta"`).
#' @return An object of class `fastphylosig_tree_check`.  It contains
#'   `ready`, `ready_by_signal`, `tree_summary`, and an `issues` data frame
#'   with `code`, `severity`, `signal`, `message`, `check`, `problem`,
#'   `action`, and `auto_fixable` columns.  `issues` is
#'   diagnostic only; this function does not stop for an invalid tree.
#' @rdname check_tree
.inspect_tree_core <- function(tree, signal = c("K", "lambda", "D", "Delta")) {
  supported <- c("K", "lambda", "D", "Delta")
  signal_missing <- missing(signal)
  if (is.null(signal)) signal <- character()
  if (!is.character(signal) || anyNA(signal) || !length(signal)) {
    stop("signal must contain one or more of K, lambda, D, and Delta.",
         call. = FALSE)
  }

  # Accept the usual spelling variants while keeping the canonical names in
  # output.  Unknown names are an argument error; tree problems themselves are
  # returned as diagnostics instead of being raised as errors.
  signal_key <- tolower(trimws(signal))
  signal_key[signal_key == "k"] <- "k"
  signal_key[signal_key %in% c("lambda", "pagel", "pagel's lambda")] <-
    "lambda"
  signal_key[signal_key %in% c("d", "phylo.d", "phylo_d")] <- "d"
  signal_key[signal_key %in% c("delta", "phylo.delta", "phylo_delta")] <-
    "delta"
  canonical <- c(k = "K", lambda = "lambda", d = "D", delta = "Delta")
  if (any(!signal_key %in% names(canonical))) {
    stop("signal must contain only K, lambda, D, and Delta.", call. = FALSE)
  }
  selected <- unname(canonical[signal_key])
  selected <- selected[!duplicated(selected)]
  # Validate that every public method has a registered tree contract.  The
  # registry is descriptive; existing issue logic below remains unchanged.
  registry <- .tree_requirement_registry()
  if (any(!selected %in% names(registry))) {
    stop("internal tree requirement registry is incomplete.", call. = FALSE)
  }
  # `missing()` is deliberately retained as a small aid for debugging and to
  # make it explicit that the default is the complete method set.
  if (isTRUE(signal_missing)) selected <- supported

  # A prepared context already contains the complete read-only inspection.
  # Reuse it after fingerprint validation instead of rescanning the same edge
  # list on every analysis. Older serialized contexts safely fall through to
  # the ordinary scan when the cached field is absent.
  if (inherits(tree, "fastphylosig_tree") && is.list(tree) &&
      is.list(tree$inspection)) {
    .validate_prepared_context(tree)
    return(.tree_check_select(tree$inspection, selected, prepared = TRUE))
  }

  severity_levels <- c("INFO", "WARNING", "ERROR", "USER_ACTION_REQUIRED")
  blocking <- c("ERROR", "USER_ACTION_REQUIRED")
  issue_code <- character()
  issue_severity <- character()
  issue_signal <- character()
  issue_message <- character()
  issue_check <- character()
  issue_problem <- character()
  issue_action <- character()
  issue_auto_fixable <- logical()

  # Keep one issue vocabulary for every preparation entry point.  The
  # representation-only root-ID diagnostics are the only diagnostics that
  # can be repaired safely without changing the biological tree.
  issue_metadata <- function(code, severity, signal_name, message) {
    .tree_issue_metadata(code, severity, signal_name, message)
  }

  add_issue <- function(code, severity, signal_name, message) {
    severity <- as.character(severity[[1L]])
    signal_name <- as.character(signal_name[[1L]])
    if (!severity %in% severity_levels) {
      stop("internal check_tree error: unsupported issue severity.",
           call. = FALSE)
    }
    issue_code <<- c(issue_code, as.character(code[[1L]]))
    issue_severity <<- c(issue_severity, severity)
    issue_signal <<- c(issue_signal, signal_name)
    issue_message <<- c(issue_message, as.character(message[[1L]]))
    metadata <- issue_metadata(code, severity, signal_name, message)
    issue_check <<- c(issue_check, metadata$check)
    issue_problem <<- c(issue_problem, metadata$problem)
    issue_action <<- c(issue_action, metadata$action)
    issue_auto_fixable <<- c(issue_auto_fixable,
                             isTRUE(metadata$auto_fixable))
    invisible(NULL)
  }
  add_for <- function(code, severity, methods, message) {
    methods <- intersect(as.character(methods), selected)
    if (!length(methods)) return(invisible(NULL))
    for (method in methods) add_issue(code, severity, method, message)
    invisible(NULL)
  }
  add_global <- function(code, severity, message) {
    # A global diagnostic is intentionally represented once.  Readiness checks
    # below treat signal = "all" as applying to every requested method.
    add_issue(code, severity, "all", message)
    invisible(NULL)
  }

  # Work with the underlying phylo object only.  No helper here calls
  # reorder.phylo(), drop.tip(), prepare_tree(), or any mutating replacement.
  prepared <- inherits(tree, "fastphylosig_tree")
  phy <- if (prepared && is.list(tree)) tree$tree else tree
  input_is_phylo <- inherits(phy, "phylo")
  if (!input_is_phylo) {
    add_global(
      "invalid_input", "ERROR",
      "tree must be a phylo object or a valid fastphylosig_tree context."
    )
  }

  # Defaults keep the returned summary stable even for NULL, malformed, or
  # partially constructed objects.
  tip_label <- if (is.list(phy)) phy$tip.label else NULL
  n_tip <- if (is.null(tip_label)) 0L else length(tip_label)
  if (!is.numeric(n_tip) || length(n_tip) != 1L || !is.finite(n_tip)) {
    n_tip <- 0L
  }
  n_tip <- as.integer(n_tip)
  labels_present <- !is.null(tip_label) && length(tip_label) > 0L
  labels_character <- labels_present &&
    (is.character(tip_label) || is.factor(tip_label))
  labels <- if (labels_present) {
    tryCatch(as.character(tip_label), error = function(e) rep(NA_character_,
                                                                n_tip))
  } else {
    character()
  }
  n_missing_labels <- if (length(labels)) sum(is.na(labels)) else 0L
  n_empty_labels <- if (length(labels)) {
    sum(!is.na(labels) & !nzchar(trimws(labels)))
  } else {
    0L
  }
  n_duplicate_labels <- if (length(labels)) {
    z <- labels[!is.na(labels)]
    if (length(z)) sum(duplicated(z)) else 0L
  } else {
    0L
  }

  if (!labels_present || !labels_character || n_missing_labels > 0L) {
    add_global(
      "missing_tip_labels", "ERROR",
      "tree tip labels must be present, character, and non-missing."
    )
  }
  if (n_empty_labels > 0L) {
    add_global(
      "empty_tip_label", "ERROR",
      sprintf("tree contains %d empty or whitespace-only tip label(s).",
              n_empty_labels)
    )
  }
  if (n_duplicate_labels > 0L) {
    add_global(
      "duplicate_tip_labels", "ERROR",
      sprintf("tree contains %d duplicated tip label(s).", n_duplicate_labels)
    )
  }
  if (n_tip < 2L) {
    add_global(
      "too_few_tips", "ERROR",
      "tree must contain at least two tips for signal calculations."
    )
  }

  edge <- if (is.list(phy)) phy$edge else NULL
  edge_matrix <- is.matrix(edge) && length(dim(edge)) == 2L &&
    ncol(edge) == 2L && nrow(edge) >= 1L
  edge_numeric <- edge_matrix &&
    (is.numeric(edge) || is.integer(edge))
  edge_values <- if (edge_numeric) {
    suppressWarnings(as.numeric(edge))
  } else {
    numeric()
  }
  edge_integer <- edge_numeric && length(edge_values) > 0L &&
    all(is.finite(edge_values)) &&
    all(edge_values == floor(edge_values))
  if (!edge_matrix) {
    add_global("missing_edge", "ERROR",
               "tree must contain a two-column edge matrix.")
  } else if (!edge_numeric || !edge_integer) {
    add_global("invalid_edge", "ERROR",
               "tree edge node numbers must be finite integers.")
  }

  declared_nnode <- if (is.list(phy)) phy$Nnode else NULL
  declared_nnode_ok <- is.numeric(declared_nnode) &&
    length(declared_nnode) == 1L && is.finite(declared_nnode) &&
    declared_nnode == floor(declared_nnode) && declared_nnode >= 0
  if (!is.null(declared_nnode) && !declared_nnode_ok) {
    add_global(
      "invalid_nnode", "ERROR",
      "Nnode must be one finite non-negative integer consistent with edge."
    )
  }
  declared_nnode_value <- if (declared_nnode_ok) as.integer(declared_nnode)
    else NA_integer_

  n_total <- if (edge_integer && length(edge_values)) {
    max(c(n_tip, edge_values), na.rm = TRUE)
  } else if (declared_nnode_ok) {
    n_tip + declared_nnode_value
  } else {
    n_tip
  }
  if (!is.finite(n_total) || n_total < n_tip) n_total <- n_tip
  n_total <- as.integer(n_total)
  inferred_nnode <- max(0L, n_total - n_tip)
  if (declared_nnode_ok && declared_nnode_value != inferred_nnode) {
    add_global(
      "nnode_mismatch", "ERROR",
      sprintf("Nnode (%d) does not match the internal node IDs implied by edge (%d).",
              declared_nnode_value, inferred_nnode)
    )
  }

  parent <- child <- integer()
  edge_index <- integer()
  edge_range_ok <- FALSE
  topology_basic <- FALSE
  root <- NA_integer_
  root_candidates <- integer()
  indegree <- outdegree <- integer()
  connected <- FALSE
  tips_ok <- internals_ok <- FALSE
  n_polytomy <- n_single_child <- 0L
  root_degree <- NA_integer_
  if (edge_integer && length(edge_values) && n_total >= 1L) {
    edge_mat <- matrix(edge_values, nrow = nrow(edge), ncol = 2L)
    parent <- as.integer(edge_mat[, 1L])
    child <- as.integer(edge_mat[, 2L])
    edge_index <- seq_along(child)
    edge_range_ok <- all(parent >= 1L, parent <= n_total,
                         child >= 1L, child <= n_total)
    if (edge_range_ok) {
      indegree <- tabulate(child, nbins = n_total)
      outdegree <- tabulate(parent, nbins = n_total)
      tip_ids <- if (n_tip) seq_len(n_tip) else integer()
      internal_ids <- if (n_total > n_tip) seq.int(n_tip + 1L, n_total)
        else integer()
      tips_ok <- !length(tip_ids) || all(indegree[tip_ids] == 1L &
                                           outdegree[tip_ids] == 0L)
      internals_ok <- !length(internal_ids) ||
        all(indegree[internal_ids] <= 1L & outdegree[internal_ids] >= 1L)
      root_candidates <- internal_ids[indegree[internal_ids] == 0L]
      if (length(root_candidates) == 1L) root <- root_candidates[[1L]]
      if (length(internal_ids)) {
        n_polytomy <- sum(outdegree[internal_ids] > 2L)
        n_single_child <- sum(outdegree[internal_ids] == 1L)
      }
      root_degree <- if (is.finite(root)) outdegree[[root]] else NA_integer_

      # A valid rooted edge list is a connected acyclic graph with exactly one
      # source internal node.  This structural criterion deliberately accepts
      # root polytomies; conventional ape::is.rooted() is reported separately.
      if (is.finite(root) && length(parent) == n_total - 1L &&
          !anyDuplicated(child) && !any(parent == child) && tips_ok &&
          internals_ok) {
        children <- split(child, parent)
        seen <- rep(FALSE, n_total)
        stack <- root
        seen[[root]] <- TRUE
        while (length(stack)) {
          node <- stack[[length(stack)]]
          stack <- stack[-length(stack)]
          kids <- children[[as.character(node)]]
          if (is.null(kids)) next
          for (kid in kids) {
            kid <- as.integer(kid)
            if (seen[[kid]]) {
              # Repeated nodes/cycles are not connected rooted trees.
              stack <- integer()
              seen[] <- FALSE
              break
            }
            seen[[kid]] <- TRUE
            stack <- c(stack, kid)
          }
        }
        connected <- all(seen)
      }
      topology_basic <- edge_range_ok && tips_ok && internals_ok &&
        length(root_candidates) == 1L && connected &&
        length(parent) == n_total - 1L && !anyDuplicated(child) &&
        !any(parent == child)
    }
  }
  if (edge_matrix && edge_numeric && edge_integer && !edge_range_ok) {
    add_global("invalid_topology", "ERROR",
               "tree edge node numbers are outside the contiguous node range.")
  }
  if (edge_matrix && edge_numeric && edge_integer && !topology_basic) {
    add_global(
      "invalid_topology", "ERROR",
      "tree edge topology must be a connected rooted tree with one parent per node."
    )
  }

  # Branch-length diagnostics are kept independent of topology diagnostics so
  # callers can see all actionable problems in one pass.
  edge_length <- if (is.list(phy)) phy$edge.length else NULL
  branch_present <- !is.null(edge_length)
  branch_numeric <- branch_present &&
    (is.numeric(edge_length) || is.integer(edge_length))
  branch_values <- if (branch_numeric) suppressWarnings(as.numeric(edge_length))
    else numeric()
  branch_length_ok <- branch_numeric && edge_matrix &&
    length(branch_values) == nrow(edge)
  if (!branch_present) {
    add_global("missing_branch_lengths", "ERROR",
               "tree must contain one branch length for every edge.")
  } else if (!branch_numeric) {
    add_global("invalid_branch_lengths", "ERROR",
               "tree branch lengths must be numeric.")
  } else if (edge_matrix && length(branch_values) != nrow(edge)) {
    add_global(
      "branch_length_length_mismatch", "ERROR",
      "tree branch lengths must have one value for every edge."
    )
  }
  branch_finite <- branch_length_ok && all(is.finite(branch_values))
  if (branch_length_ok && !branch_finite) {
    add_global("nonfinite_branch_lengths", "ERROR",
               "tree branch lengths must be finite.")
  }
  branch_negative <- branch_length_ok && any(branch_values < 0, na.rm = TRUE)
  if (branch_negative) {
    add_global("negative_branch_lengths", "ERROR",
               "tree branch lengths must be non-negative.")
  }
  branch_valid_for_geometry <- branch_length_ok && branch_finite &&
    !branch_negative

  zero_terminal <- zero_internal <- near_zero <- 0L
  near_zero_threshold <- NA_real_
  if (branch_valid_for_geometry && edge_integer && length(child)) {
    terminal_edge <- child <= n_tip
    zero_terminal <- sum(branch_values[terminal_edge] == 0)
    zero_internal <- sum(branch_values[!terminal_edge] == 0)
    positive <- branch_values[branch_values > 0]
    if (length(positive)) {
      near_zero_threshold <- sqrt(.Machine$double.eps) *
        max(1, max(branch_values))
      near_zero <- sum(branch_values > 0 & branch_values <= near_zero_threshold)
    }
  }

  # Distances are computed from the edge list rather than ape::node.depth.* so
  # a malformed object remains diagnosable and no package helper can reorder or
  # mutate it.  The distances also expose zero-height terminal variances and
  # overflow before a production kernel is called.
  root_distance <- rep(NA_real_, n_total)
  root_distance_overflow <- FALSE
  if (topology_basic && branch_valid_for_geometry && length(edge_index)) {
    root_distance[] <- NA_real_
    root_distance[[root]] <- 0
    children <- split(seq_along(child), parent)
    stack <- root
    while (length(stack)) {
      node <- stack[[length(stack)]]
      stack <- stack[-length(stack)]
      edge_rows <- children[[as.character(node)]]
      if (is.null(edge_rows)) next
      for (row in edge_rows) {
        kid <- child[[row]]
        distance <- root_distance[[node]] + branch_values[[row]]
        if (!is.finite(distance)) {
          root_distance_overflow <- TRUE
          next
        }
        root_distance[[kid]] <- distance
        stack <- c(stack, kid)
      }
    }
  }
  if (root_distance_overflow) {
    add_global("root_distance_overflow", "ERROR",
               "root-to-tip distances overflow to a non-finite value.")
  }
  tip_distance <- if (length(root_distance) && n_tip >= 1L) {
    root_distance[seq_len(n_tip)]
  } else {
    numeric()
  }
  finite_tip_distance <- length(tip_distance) > 0L &&
    all(is.finite(tip_distance))
  root_height <- if (finite_tip_distance) max(tip_distance) else NA_real_
  ultrametric <- if (finite_tip_distance) {
    tol <- 1e-8 * max(1, abs(root_height))
    max(tip_distance) - min(tip_distance) <= tol
  } else {
    NA
  }
  structural_rooted <- isTRUE(topology_basic)
  # ape's rooted convention is represented by a degree-two structural root.
  # This keeps a rooted tree with a non-root polytomy rooted, while exposing a
  # root degree-three encoding (the usual unrooted representation) as
  # conventional_unrooted.  The structural root is still retained for the
  # K/lambda/D compatibility paths.
  conventional_rooted <- structural_rooted && n_tip >= 2L &&
    identical(root_degree, 2L)
  conventional_unrooted <- structural_rooted && !conventional_rooted

  # Conventional rootedness is useful to callers familiar with ape, but root
  # polytomies have the same edge representation as an unrooted encoding.  The
  # production K/lambda/D kernels use the structural source and can process a
  # rooted polytomy; Delta additionally requires a binary/rooted ACE input.
  if (!structural_rooted && edge_matrix && edge_numeric && edge_integer) {
    add_for("requires_rooted_tree", "ERROR", selected,
            "tree topology has no unique connected structural root.")
  } else if (conventional_unrooted) {
    add_for("conventional_unrooted", "WARNING",
            intersect(selected, c("K", "lambda", "D")),
            "ape's conventional rooted flag is false; a structural root is present, so K/lambda/D can proceed, but verify the intended rooting.")
    add_for("delta_requires_conventional_root", "ERROR", "Delta",
            "Delta/fast_ace requires a biologically justified conventional root with two descendants.")
  }

  if (n_polytomy > 0L) {
    add_for("polytomy_supported", "INFO",
            intersect(selected, c("K", "lambda", "D")),
            sprintf("tree contains %d polytomous internal node(s); this method uses its compatibility path.",
                    n_polytomy))
    add_for("d_polytomy_compatibility", "WARNING", "D",
            "D uses a compatibility split for polytomies; inspect the result and simulation diagnostics.")
    add_for("delta_requires_binary_tree", "ERROR", "Delta",
            "Delta/fast_ace requires a rooted fully-dichotomous tree; polytomies must be resolved first.")
  }
  if (n_single_child > 0L) {
    add_for("single_child_internal", "WARNING",
            intersect(selected, c("K", "lambda")),
            sprintf("tree contains %d single-child internal node(s); results use a degenerate compatibility path.",
                    n_single_child))
    add_for("d_requires_no_single_child", "ERROR", "D",
            "D's traversal does not retain a single-child branch; remove or collapse unary internal nodes first.")
    add_for("delta_requires_binary_tree", "ERROR", "Delta",
            "Delta/fast_ace does not support single-child internal nodes.")
  }
  if (branch_valid_for_geometry && (zero_terminal + zero_internal) > 0L) {
    add_for("k_requires_positive_branches", "ERROR", "K",
            sprintf("K's tree engine requires strictly positive branch lengths (%d zero branch(es) found).",
                    zero_terminal + zero_internal))
    if (zero_terminal > 0L) {
      add_for("lambda_zero_terminal_branch", "USER_ACTION_REQUIRED", "lambda",
              sprintf("lambda's upper-bound covariance is singular with %d zero terminal branch(es); verify or repair these branches before fitting.",
                      zero_terminal))
    }
    if (zero_internal > 0L) {
      add_for("lambda_zero_internal_branch", "WARNING", "lambda",
              sprintf("lambda can evaluate %d zero internal branch(es), but boundary likelihoods may be numerically delicate.",
                      zero_internal))
    }
    add_for("d_requires_positive_branches", "ERROR", "D",
            sprintf("D requires strictly positive branch lengths (%d zero branch(es) found).",
                    zero_terminal + zero_internal))
    add_for("delta_requires_positive_branches", "ERROR", "Delta",
            sprintf("Delta requires strictly positive branch lengths (%d zero branch(es) found).",
                    zero_terminal + zero_internal))
  }
  if (branch_valid_for_geometry && isTRUE(structural_rooted) &&
      length(tip_distance) && any(tip_distance <= 0, na.rm = TRUE)) {
    add_for("lambda_nonpositive_tip_height", "ERROR", "lambda",
            "lambda requires a positive root-to-tip variance for every tip.")
  }
  if (isTRUE(structural_rooted) && is.finite(root) && n_tip >= 2L &&
      root != n_tip + 1L) {
    # The D C++ traversal and ape's canonical pruning order use n_tip + 1 as
    # the root node.  K/lambda kernels use the explicit structural root.
    add_for("d_requires_canonical_root", "ERROR", "D",
            "D requires the canonical phylo root node n_tip + 1.")
    add_for("delta_requires_canonical_root", "ERROR", "Delta",
            "Delta/fast_ace requires the canonical phylo root node n_tip + 1.")
  }
  if (branch_valid_for_geometry && near_zero > 0L) {
    add_for("near_zero_branch", "WARNING", selected,
            sprintf("tree contains %d positive branch length(s) at or below the near-zero threshold (%.3g).",
                    near_zero, near_zero_threshold))
  }

  if (isTRUE(finite_tip_distance) && isFALSE(ultrametric)) {
    add_for("nonultrametric_tree", "INFO", selected,
            "tree is non-ultrametric; K, lambda, D, and Delta use branch lengths without requiring ultrametricity.")
  }

  # Delta's binary topology requirement is independent of its conventional
  # root requirement above. Polytomy and unary diagnostics already carry the
  # stable binary-topology issue code.
  internal_ids_for_shape <- if (n_total > n_tip) {
    seq.int(n_tip + 1L, n_total)
  } else {
    integer()
  }
  fully_dichotomous <- structural_rooted && conventional_rooted &&
    length(internal_ids_for_shape) > 0L &&
    all(outdegree[internal_ids_for_shape] == 2L)

  # Ensure Delta's positive branch contract is still reported if branch data
  # were malformed (the global branch diagnostic then carries the blocking
  # status, while this method-specific note tells callers why Delta cannot run).
  if (!branch_valid_for_geometry && "Delta" %in% selected) {
    add_issue("delta_requires_positive_branches", "ERROR", "Delta",
              "Delta requires finite, strictly positive branch lengths.")
  }
  if (!branch_valid_for_geometry && "D" %in% selected) {
    add_issue("d_requires_positive_branches", "ERROR", "D",
              "D requires finite, strictly positive branch lengths.")
  }
  if (!branch_valid_for_geometry && "K" %in% selected) {
    add_issue("k_requires_positive_branches", "ERROR", "K",
              "K's production tree engine requires finite, strictly positive branch lengths.")
  }
  if (!branch_valid_for_geometry && "lambda" %in% selected) {
    add_issue("lambda_requires_valid_branches", "ERROR", "lambda",
              "lambda requires finite, non-negative branch lengths and positive tip heights.")
  }

  # If malformed labels/topology/branch data prevented a method from reaching
  # its method-specific checks, the global diagnostics still block readiness.
  issues <- data.frame(
    code = issue_code,
    severity = issue_severity,
    signal = issue_signal,
    message = issue_message,
    check = issue_check,
    problem = issue_problem,
    action = issue_action,
    auto_fixable = issue_auto_fixable,
    stringsAsFactors = FALSE
  )
  ready_by_signal <- stats::setNames(vapply(selected, function(method) {
    if (!nrow(issues)) return(TRUE)
    hit <- issues$signal %in% c("all", method) &
      issues$severity %in% blocking &
      !(issues$auto_fixable %in% TRUE &
          tolower(as.character(issues$check)) == "representation")
    !any(hit)
  }, logical(1)), selected)

  branch_min <- if (branch_valid_for_geometry && length(branch_values)) {
    min(branch_values)
  } else {
    NA_real_
  }
  branch_max <- if (branch_valid_for_geometry && length(branch_values)) {
    max(branch_values)
  } else {
    NA_real_
  }
  tree_summary <- list(
    tips = n_tip,
    n_tip = n_tip,
    n_tips = n_tip,
    internal_nodes = inferred_nnode,
    n_node = inferred_nnode,
    n_internal = inferred_nnode,
    edges = if (edge_matrix) nrow(edge) else 0L,
    n_edge = if (edge_matrix) nrow(edge) else 0L,
    labels_present = labels_present,
    labels_character = labels_character,
    labels_unique = n_duplicate_labels == 0L && n_missing_labels == 0L,
    n_missing_labels = n_missing_labels,
    n_empty_labels = n_empty_labels,
    n_duplicate_labels = n_duplicate_labels,
    branch_lengths = branch_present && branch_length_ok,
    branch_lengths_present = branch_present,
    branch_lengths_finite = branch_finite,
    branch_lengths_nonnegative = branch_valid_for_geometry,
    branch_lengths_positive = branch_valid_for_geometry &&
      length(branch_values) > 0L && all(branch_values > 0),
    branch_min = branch_min,
    branch_max = branch_max,
    branch_length_count = if (branch_numeric) length(branch_values) else 0L,
    rooted = isTRUE(conventional_rooted),
    conventional_rooted = isTRUE(conventional_rooted),
    structural_rooted = isTRUE(structural_rooted),
    root = if (is.finite(root)) root else NA_integer_,
    root_degree = root_degree,
    topology_valid = isTRUE(topology_basic),
    connected = isTRUE(connected),
    polytomies = n_polytomy,
    n_polytomies = n_polytomy,
    has_polytomies = n_polytomy > 0L,
    single_child = n_single_child,
    n_single_child = n_single_child,
    has_single_child = n_single_child > 0L,
    zero_terminal = zero_terminal,
    n_zero_terminal = zero_terminal,
    zero_internal = zero_internal,
    n_zero_internal = zero_internal,
    near_zero = near_zero,
    n_near_zero = near_zero,
    near_zero_threshold = near_zero_threshold,
    root_height = root_height,
    ultrametric = ultrametric,
    fully_dichotomous = isTRUE(fully_dichotomous),
    prepared = prepared,
    valid = input_is_phylo && labels_present && labels_character &&
      n_missing_labels == 0L && n_empty_labels == 0L &&
      n_duplicate_labels == 0L && n_tip >= 2L &&
      (is.null(declared_nnode) || declared_nnode_ok) && topology_basic &&
      branch_valid_for_geometry
  )

  out <- list(
    ready = all(ready_by_signal),
    ready_by_signal = ready_by_signal,
    tree_summary = tree_summary,
    issues = issues
  )
  class(out) <- c("fastphylosig_tree_check", "list")
  out
}

#' @rdname check_tree
#' @export
check_tree <- function(tree, signal = c("K", "lambda", "D", "Delta")) {
  .inspect_tree_core(tree, signal = signal)
}

#' @export
print.fastphylosig_tree_check <- function(x, ...) {
  if (!inherits(x, "fastphylosig_tree_check")) {
    stop("x must be a fastphylosig_tree_check object.", call. = FALSE)
  }
  s <- x$tree_summary
  cat("Tree check\n\n")
  cat("Tips: ", if (length(s$n_tip)) as.character(s$n_tip[[1L]]) else "NA",
      "\n\n", sep = "")
  ready <- x$ready_by_signal
  if (length(ready)) {
    width <- max(nchar(names(ready))) + 2L
    for (i in seq_along(ready)) {
      cat(sprintf("%-*s%s\n", width, names(ready)[[i]],
                  if (isTRUE(ready[[i]])) "READY" else "NOT READY"))
    }
  }
  issues <- x$issues
  if (is.null(issues) || !nrow(issues)) {
    cat("\nIssues: none\n")
  } else {
    cat("\nIssues:\n")
    key <- paste(issues$severity, issues$message, sep = "\r")
    first <- which(!duplicated(key))
    max_show <- min(6L, length(first))
    for (j in seq_len(max_show)) {
      i <- first[[j]]
      same <- key == key[[i]]
      methods <- unique(as.character(issues$signal[same]))
      methods <- methods[!is.na(methods) & nzchar(methods) & methods != "all"]
      scope <- if (length(methods)) paste0("; ", paste(methods, collapse = "/")) else ""
      cat(sprintf("- [%s%s] %s\n", issues$severity[[i]], scope,
                  issues$message[[i]]))
    }
    if (length(first) > max_show) {
      cat(sprintf("- ... %d more issue(s)\n", length(first) - max_show))
    }
  }
  invisible(x)
}
