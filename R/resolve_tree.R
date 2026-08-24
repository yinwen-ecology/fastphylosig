# Conservative phylogenetic-tree resolution ---------------------------------

#' Resolve a phylogenetic tree for a signal calculation
#'
#' `resolve_tree()` performs only representation-preserving repairs that can
#' be made without changing the biological tree: edge ordering and internal
#' node numbering are canonicalised, and an unrooted tree is rooted only when
#' the selected signal requires a root and the caller supplies an explicit
#' outgroup.  It deliberately does not infer an outgroup, midpoint-root,
#' resolve polytomies, collapse unary nodes, alter branch lengths, add a small
#' value to zero branches, or drop tips.
#'
#' @param tree An object of class `"phylo"`.
#' @param signal One signal name: `"K"`, `"lambda"`, `"D"`, or `"Delta"`.
#'   Matching is case-insensitive; exactly one signal must be supplied.
#' @param outgroup Optional character vector of one or more existing tip
#'   labels.  It is used only when `check_tree()` reports that this signal
#'   requires rooting.
#' @return A `phylo` object.  The returned object has class
#'   `"fastphylosig_resolved_tree"` in addition to `"phylo"` and carries a
#'   `fastphylosig_resolution` attribute containing `signal`,
#'   `original_status`, `changes`, `final_status`, and `ready`.
#' @export
resolve_tree <- function(tree, signal, outgroup = NULL) {
  if (!inherits(tree, "phylo")) {
    stop("tree should be an object of class \"phylo\".", call. = FALSE)
  }
  signal <- .resolve_tree_signal(signal)

  # Validate the outgroup even when the current tree happens to be rooted.  A
  # misspelled explicit outgroup should never be silently ignored.
  if (!is.null(outgroup)) {
    outgroup <- .resolve_tree_outgroup(outgroup, tree$tip.label)
  }

  original <- check_tree(tree, signal = signal)
  original_ready <- .resolve_tree_ready(original, signal)
  original_status <- if (isTRUE(original_ready)) "READY" else "NOT READY"
  changes <- character()
  current <- tree

  # Reordering and internal-node relabelling are representation operations;
  # neither changes tip labels, topology, or branch-length values.
  normalized <- .safe_canonicalize_core(current)
  if (!.resolve_tree_representation_equal(current, normalized)) {
    current <- normalized
    changes <- c(changes, "representation normalized (edge order/node numbering)")
    checked <- check_tree(current, signal = signal)
  } else {
    checked <- original
  }

  # Root only in response to a method-specific blocking root issue and only
  # with an explicit, validated outgroup.  resolve.root = FALSE is important:
  # it avoids silently resolving an existing root polytomy.
  if (.resolve_tree_requires_root(checked, signal)) {
    if (is.null(outgroup)) {
      .resolve_tree_fail(checked, signal,
                         "an explicit outgroup is required to root this tree")
    }
    rooted <- tryCatch(
      ape::root(current, outgroup = outgroup, resolve.root = FALSE),
      error = function(e) e
    )
    if (inherits(rooted, "error")) {
      stop(
        paste0(
          "resolve_tree() could not root the tree with the supplied outgroup: ",
          conditionMessage(rooted),
          ". USER_ACTION_REQUIRED."
        ),
        call. = FALSE
      )
    }
    if (!inherits(rooted, "phylo")) {
      stop("resolve_tree() rooting did not return a phylo object.",
           call. = FALSE)
    }
    if (!.resolve_tree_patristic_equal(current, rooted)) {
      stop(
        paste0(
          "resolve_tree() refused a rooting operation that changed tip-to-tip ",
          "distances or branch-length totals. USER_ACTION_REQUIRED."
        ),
        call. = FALSE
      )
    }
    current <- .safe_canonicalize_core(rooted)
    changes <- c(changes, "rooted using explicit outgroup")
    checked <- check_tree(current, signal = signal)
  }

  final_ready <- .resolve_tree_ready(checked, signal)
  final_status <- if (isTRUE(final_ready)) "READY" else "NOT READY"
  if (!isTRUE(final_ready)) {
    .resolve_tree_fail(checked, signal,
                       "remaining blocking tree requirements cannot be repaired safely")
  }

  # Keep the metadata intentionally small and deterministic.  In particular,
  # do not embed the complete check object or tree in an attribute.
  changes <- unique(as.character(changes))
  attr(current, "fastphylosig_resolution") <- list(
    signal = signal,
    original_status = original_status,
    changes = changes,
    final_status = final_status,
    ready = isTRUE(final_ready)
  )
  class(current) <- unique(c("fastphylosig_resolved_tree", class(current)))
  current
}

.resolve_tree_signal <- function(signal) {
  if (missing(signal) || !is.character(signal) || length(signal) != 1L ||
      is.factor(signal) ||
      is.na(signal)) {
    stop(
      "signal must be exactly one of \"K\", \"lambda\", \"D\", or \"Delta\".",
      call. = FALSE
    )
  }
  value <- as.character(signal)
  key <- tolower(value)
  switch(
    key,
    k = "K",
    lambda = "lambda",
    d = "D",
    delta = "Delta",
    stop(
      "signal must be exactly one of \"K\", \"lambda\", \"D\", or \"Delta\".",
      call. = FALSE
    )
  )
}

.resolve_tree_outgroup <- function(outgroup, tip_label) {
  if (!is.character(outgroup) || !length(outgroup) || anyNA(outgroup) ||
      anyDuplicated(outgroup) || any(!nzchar(outgroup))) {
    stop(
      "outgroup must be one or more unique, non-missing existing tip labels.",
      call. = FALSE
    )
  }
  if (is.null(tip_label) || any(!outgroup %in% tip_label)) {
    missing_labels <- setdiff(outgroup, tip_label)
    stop(
      paste0(
        "outgroup contains tip label(s) not present in tree: ",
        paste(missing_labels, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }
  if (length(outgroup) >= length(tip_label)) {
    stop("outgroup must leave at least one non-outgroup tip.", call. = FALSE)
  }
  outgroup
}

.resolve_tree_ready <- function(check, signal) {
  by_signal <- check$ready_by_signal
  if (!is.null(by_signal) && length(by_signal)) {
    value <- unname(by_signal[[signal]])
    if (length(value)) return(isTRUE(value))
  }
  if (!is.null(check$ready) && length(check$ready) == 1L) {
    return(isTRUE(check$ready))
  }
  FALSE
}

.resolve_tree_issue_table <- function(check) {
  issues <- check$issues
  if (!is.data.frame(issues) || !nrow(issues)) return(NULL)
  if (!"severity" %in% names(issues)) issues$severity <- "ERROR"
  if (!"code" %in% names(issues)) issues$code <- ""
  if (!"message" %in% names(issues)) issues$message <- ""
  if (!"check" %in% names(issues)) issues$check <- "tree"
  if (!"problem" %in% names(issues)) issues$problem <- issues$message
  if (!"action" %in% names(issues)) issues$action <- "repair the tree and re-run check_tree()"
  if (!"auto_fixable" %in% names(issues)) issues$auto_fixable <- FALSE
  issues
}

.resolve_tree_requires_root <- function(check, signal) {
  issues <- .resolve_tree_issue_table(check)
  if (is.null(issues)) return(FALSE)
  severe <- toupper(as.character(issues$severity)) %in%
    c("ERROR", "USER_ACTION_REQUIRED")
  text <- paste(as.character(issues$code), as.character(issues$message))
  root_issue <- severe & grepl("root|unroot|rooted", text, ignore.case = TRUE)
  if (!any(root_issue)) return(FALSE)
  # A root issue for another signal should not trigger a mutation.  The
  # check_tree contract normally filters issues by selected signal, but this
  # guard also handles callers returning a multi-signal issue table.
  if ("signal" %in% names(issues)) {
    issue_signal <- as.character(issues$signal)
    issue_signal <- issue_signal[!is.na(issue_signal) & nzchar(issue_signal)]
    if (length(issue_signal) && !any(tolower(issue_signal) == tolower(signal))) {
      return(FALSE)
    }
  }
  TRUE
}

.resolve_tree_fail <- function(check, signal, reason) {
  stop(.format_actionable_condition(
    check = check, signal = signal, prefix = "resolve_tree() cannot produce a ready tree",
    reason = reason
  ), call. = FALSE)
}

.resolve_tree_representation_equal <- function(x, y) {
  if (!inherits(y, "phylo")) return(FALSE)
  same_edge <- identical(unname(as.matrix(x$edge)), unname(as.matrix(y$edge)))
  same_length <- if (is.null(x$edge.length) || is.null(y$edge.length)) {
    is.null(x$edge.length) && is.null(y$edge.length)
  } else {
    isTRUE(all.equal(as.numeric(x$edge.length), as.numeric(y$edge.length),
                     check.attributes = FALSE))
  }
  same_tips <- identical(as.character(x$tip.label), as.character(y$tip.label))
  same_nodes <- identical(as.character(x$node.label), as.character(y$node.label))
  same_edge && same_length && same_tips && same_nodes
}

.resolve_tree_patristic_equal <- function(x, y, tolerance = 1e-8) {
  if (is.null(x$edge.length) || is.null(y$edge.length)) return(TRUE)
  # ape::cophenetic.phylo is the relevant invariant for rerooting: it compares
  # all tip-to-tip paths while allowing edge rows and node ids to change.
  dx <- tryCatch(ape::cophenetic.phylo(x), error = function(e) NULL)
  dy <- tryCatch(ape::cophenetic.phylo(y), error = function(e) NULL)
  if (is.null(dx) || is.null(dy) || !identical(sort(rownames(dx)), sort(rownames(dy)))) {
    return(FALSE)
  }
  dy <- dy[rownames(dx), colnames(dx), drop = FALSE]
  ok_dist <- isTRUE(all.equal(unname(dx), unname(dy),
                              tolerance = tolerance, check.attributes = FALSE))
  total_x <- sum(as.numeric(x$edge.length)) +
    if (is.null(x$root.edge)) 0 else sum(as.numeric(x$root.edge))
  total_y <- sum(as.numeric(y$edge.length)) +
    if (is.null(y$root.edge)) 0 else sum(as.numeric(y$root.edge))
  ok_total <- isTRUE(all.equal(total_x, total_y,
                               tolerance = tolerance, check.attributes = FALSE))
  ok_dist && ok_total
}

#' @export
print.fastphylosig_resolved_tree <- function(x, ...) {
  info <- attr(x, "fastphylosig_resolution", exact = TRUE)
  signal <- if (is.list(info) && length(info$signal)) info$signal else "unknown"
  cat("Tree resolved for ", signal, "\n\n", sep = "")
  if (is.list(info) && !is.null(info$signal)) {
    status <- if (!is.null(info$final_status) && length(info$final_status)) {
      info$final_status
    } else {
      "UNKNOWN"
    }
    cat("Changes:\n")
    if (length(info$changes)) {
      for (change in info$changes) cat("- ", change, "\n", sep = "")
    } else {
      cat("- none\n")
    }
    cat("\nStatus: ", status, "\n", sep = "")
  }
  invisible(x)
}
