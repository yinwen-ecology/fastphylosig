# Independent R orchestration and simulation oracles, kept out of the namespace.
# These deliberately do not call the native clade planner/scan/BM simulator.

# Analytic ER transitions, independent of the ACE eigensystem and simulator.
.hetero_er_states <- function(core, k, q, nsim) {
  out <- matrix(0L, core$n_tip, nsim)
  for (b in seq_len(nsim)) {
    states <- integer(core$n_total)
    states[core$root] <- as.integer(floor(runif(1, 0, k)))
    for (node in core$preorder[-1]) {
      t <- core$branch_length_by_node[node]
      dec <- if(t==0)1 else exp(-k*q*t)
      probabilities <- rep((1-dec)/k, k)
      p <- states[core$parent[node]] + 1L
      probabilities[p] <- probabilities[p] + dec
      states[node] <- which(runif(1) < cumsum(probabilities))[1L] - 1L
    }
    out[,b] <- states[seq_len(core$n_tip)] + 1L
  }
  out
}

.hetero_fixture <- function() {
  tree <- ape::stree(64, type = "balanced")
  tree$edge.length <- rep(1, nrow(tree$edge))
  nodes <- tree$edge[tree$edge[, 1L] == 65L, 2L]
  list(tree = tree, nodes = nodes,
       X = matrix(sin(seq_len(64)), ncol = 1L,
                  dimnames = list(tree$tip.label, "trait")))
}

.hetero_index_reference <- function(core) {
  count <- integer(core$n_total)
  count[seq_len(core$n_tip)] <- 1L
  size <- rep(1L, core$n_total)
  for (node in core$postorder) {
    parent <- core$parent[[node]]
    if (parent) {
      count[[parent]] <- count[[parent]] + count[[node]]
      size[[parent]] <- size[[parent]] + size[[node]]
    }
  }
  start <- integer(core$n_total)
  start[core$preorder] <- seq_len(core$n_total)
  list(n_tip = count, start = start, end = start + size - 1L)
}

.hetero_subtree_reference <- function(tree, core, index, node) {
  ids <- core$preorder[seq.int(index$start[[node]], index$end[[node]])]
  tips <- sort(ids[ids <= core$n_tip])
  internal <- ids[ids > core$n_tip]
  # Matching only local IDs avoids a whole-tree lookup allocation per clade.
  new_id <- c(seq_along(tips), length(tips) + seq_along(internal))
  old_id <- c(tips, internal)
  child <- ids[-1L]
  edge <- cbind(new_id[match(core$parent[child], old_id)],
                new_id[match(child, old_id)])
  storage.mode(edge) <- "integer"
  out <- list(edge = edge, edge.length = core$branch_length_by_node[child],
              tip.label = tree$tip.label[tips], Nnode = length(internal))
  if (!is.null(tree$node.label)) {
    out$node.label <- tree$node.label[internal - core$n_tip]
  }
  class(out) <- "phylo"
  attr(out, "source_tip_index") <- tips
  out
}

.hetero_fit_reference <- function(tree, X, method) {
  # Check the complete crown before NA pruning, as the specialist APIs do.
  checked <- .inspect_tree_core(tree, signal = method)
  if (!isTRUE(checked$ready_by_signal[[method]])) {
    stop(.format_actionable_condition(checked, method), call. = FALSE)
  }
  p <- ncol(X)
  estimate <- rep(NA_real_, p)
  status <- rep("ok", p)
  message <- rep(NA_character_, p)
  cores <- vector("list", p)
  invalid <- colSums(is.infinite(X) | is.nan(X)) > 0L
  status[invalid] <- "invalid_trait"
  message[invalid] <- "trait contains non-finite values other than NA"
  groups <- .analysis_na_groups(!is.na(X))
  for (g in seq_len(groups$n_group)) {
    columns <- groups$columns[[g]]
    columns <- columns[!invalid[columns]]
    if (!length(columns)) next
    keep <- groups$keep[[g]]
    local_tree <- if (length(keep) == nrow(X)) tree else
      ape::keep.tip(tree, tree$tip.label[keep])
    compiled <- compile_tree_cpp(local_tree$edge, local_tree$edge.length,
                                 length(local_tree$tip.label))
    cores[columns] <- rep(list(compiled), length(columns))
    local <- X[match(local_tree$tip.label, rownames(X)), columns, drop = FALSE]
    if (method == "K") {
      fitted <- fast_k_tree_batch_cpp(compiled, local, trait_chunk = 64L)
      estimate[columns] <- fitted$K
      undefined <- !is.finite(fitted$K)
      message[columns[undefined]] <- "K is undefined for a constant or numerically degenerate trait"
    } else {
      fitted <- fast_lambda_tree_optimize_cpp(compiled, local,
        max_lambda = .max_lambda(local_tree), test = FALSE, profile = FALSE,
        profile_points = 101L, trait_chunk = 64L, n_threads = 1L)
      estimate[columns] <- fitted$lambda
      undefined <- !fitted$valid
      message[columns[undefined]] <- fitted$status[undefined]
    }
    status[columns[undefined]] <- "undefined"
  }
  list(estimate = estimate, status = status, message = message, cores = cores)
}

.hetero_bm_reference <- function(core, nsim) {
  child <- core$preorder[-1L]
  increments <- matrix(stats::rnorm(length(child) * nsim), nrow = length(child)) *
    sqrt(core$branch_length_by_node[child])
  values <- matrix(0, core$n_total, nsim)
  for (i in seq_along(child)) {
    node <- child[[i]]
    values[node, ] <- values[core$parent[[node]], ] + increments[i, ]
  }
  values[seq_len(core$n_tip), , drop = FALSE]
}

# Shared independent R D recursion; loaded before the test files.
.hetero_d_contrast <- function(core, tips) {
  nv <- c(tips, rep(0, core$n_total-core$n_tip))
  lengths <- core$branch_length_by_node
  total <- 0
  for (node in core$postorder) if (node > core$n_tip) {
    child <- core$children[seq.int(core$child_ptr[node], core$child_ptr[node+1]-1)]
    v <- lengths[child]; inv <- 1/v
    nv[node] <- sum(nv[child]*inv)/sum(inv)
    total <- total+sum(abs(nv[child]-nv[node]))
    lengths[node] <- lengths[node]+prod(v)/sum(v)
  }
  total
}
