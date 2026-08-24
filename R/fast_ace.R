# Fast discrete ancestral-state reconstruction -------------------------------

fast_ace <- function(x, phy = NULL, type = "discrete", method = "ML", CI = TRUE,
                     model = c("ER", "ARD"), kappa = 1, ip = 0.1,
                     marginal = FALSE, prepared = NULL,
                     progress = interactive()) {
  runtime <- .runtime_begin("ACE", progress = progress, verbose = TRUE)
  on.exit(.runtime_on_exit(runtime), add = TRUE)
  .runtime_stage(runtime, "Checking tree...")
  if (!is.null(prepared)) phy <- prepared
  if (is.null(phy)) {
    stop("phy or prepared must be supplied.", call. = FALSE)
  }
  if (inherits(phy, "fastphylosig_tree")) {
    .validate_prepared_context(phy)
    phy <- phy$tree
  } else if (!inherits(phy, "phylo")) {
    stop("phy should be an object of class \"phylo\".", call. = FALSE)
  } else {
    .validate_prepare_tree(phy)
  }
  if (is.null(phy$edge.length)) {
    stop("tree has no branch lengths.", call. = FALSE)
  }
  # ACE's public contract requires the supplied tree itself to have strictly
  # positive branch lengths.  Check before applying kappa: with kappa = 0,
  # R's 0^0 convention would otherwise turn an input zero into one and let an
  # unsafe representation pass the post-transform check.
  stop_ace_branch_error <- function(stage) {
    stop(paste0(
      "fast_ace cannot continue.\n\n",
      "Problem:\n", stage, " contains a non-finite or non-positive branch length.\n\n",
      "Requirement:\nThe ACE production path requires strictly positive branch lengths.\n\n",
      "Issue code: ace_requires_positive_branches\n\n",
      "Recommended action:\nProvide scientifically justified positive branch lengths; ",
      "no automatic branch repair was applied. USER_ACTION_REQUIRED."
    ), call. = FALSE)
  }
  raw_edge_length <- suppressWarnings(as.numeric(phy$edge.length))
  raw_edge_rows <- if (is.matrix(phy$edge) && length(dim(phy$edge)) == 2L) {
    nrow(phy$edge)
  } else {
    NA_integer_
  }
  if (!is.finite(raw_edge_rows) || length(raw_edge_length) != raw_edge_rows ||
      any(!is.finite(raw_edge_length)) || any(raw_edge_length <= 0)) {
    stop_ace_branch_error("The input tree")
  }
  if (!is.logical(CI) || length(CI) != 1L || is.na(CI) ||
      !is.logical(marginal) || length(marginal) != 1L || is.na(marginal)) {
    stop("CI and marginal must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa) ||
      !is.numeric(ip) || length(ip) != 1L || !is.finite(ip) || ip < 0) {
    stop("kappa must be finite and ip must be finite and non-negative.",
         call. = FALSE)
  }
  type <- match.arg(type, c("discrete"))
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("method must be a single, non-missing character string.",
         call. = FALSE)
  }
  if (!identical(method, "ML")) {
    stop("fast_ace currently supports only method = \"ML\".", call. = FALSE)
  }
  model <- match.arg(model)

  nb_tip <- length(phy$tip.label)
  nb_node <- .fast_ace_topology_guard(phy)
  if (length(x) != nb_tip) {
    stop("length of trait data and tree tips do not match.", call. = FALSE)
  }
  if (!is.null(names(x))) {
    if (anyDuplicated(names(x))) {
      stop("names of x must be unique.", call. = FALSE)
    }
    if (!setequal(names(x), phy$tip.label)) {
      stop("names of x must match phy$tip.label exactly.", call. = FALSE)
    }
    x <- x[phy$tip.label]
  }
  if (is.numeric(x) && any(is.infinite(x) | is.nan(x))) {
    stop("fast_ace cannot continue: x contains non-finite values other than NA.",
         call. = FALSE)
  }
  if (kappa != 1) {
    phy$edge.length <- phy$edge.length^kappa
  }
  if (any(!is.finite(phy$edge.length)) || any(phy$edge.length <= 0)) {
    stop_ace_branch_error("The transformed tree")
  }

  .runtime_stage(runtime, "Preparing data...")
  x <- droplevels(factor(x))
  lvls <- levels(x)
  nl <- nlevels(x)
  if (nl < 2L) {
    stop("x must contain at least two observed states.", call. = FALSE)
  }
  tip_state <- as.integer(x)
  .runtime_stage(runtime, "Estimating ancestral states...")
  tree_workspace <- .fast_ace_tree_workspace(phy, nb_node = nb_node)
  workspace <- .fast_ace_workspace(tree_workspace, model = model, nl = nl)
  # CI = FALSE is the explicit likelihood-only path; avoid a numerical
  # Hessian/SE pass that the caller did not request.
  obj <- .fast_ace_fit_workspace(
    tip_state = tip_state, lvls = lvls, workspace = workspace, ip = ip,
    CI = CI, marginal = marginal, estimate_se = isTRUE(CI)
  )

  obj$call <- match.call()
  obj$engine <- "fastphylosig::fast_ace"
  class(obj) <- "ace"
  timing <- .runtime_close(runtime, success = TRUE)
  .runtime_attach(obj, timing)
}

.fast_ace_tree_workspace <- function(phy, nb_node = NULL,
                                     include_fingerprint = FALSE) {
  # Delta calls this once per retained NA-mask group.  Direct fast_ace() calls
  # have already run the same preflight, so they supply the validated Nnode.
  if (is.null(nb_node)) {
    if (is.null(phy$edge.length)) {
      stop("tree has no branch lengths.", call. = FALSE)
    }
    raw_edge_length <- suppressWarnings(as.numeric(phy$edge.length))
    raw_edge_rows <- if (is.matrix(phy$edge) && length(dim(phy$edge)) == 2L) {
      nrow(phy$edge)
    } else {
      NA_integer_
    }
    if (!is.finite(raw_edge_rows) ||
        length(raw_edge_length) != raw_edge_rows ||
        any(!is.finite(raw_edge_length)) || any(raw_edge_length <= 0)) {
      stop(paste0(
        "fast_ace cannot continue.\n\n",
        "Problem:\nThe input tree contains a non-finite or non-positive branch length.\n\n",
        "Requirement:\nThe ACE production path requires strictly positive branch lengths.\n\n",
        "Issue code: ace_requires_positive_branches\n\n",
        "Recommended action:\nProvide scientifically justified positive branch lengths; ",
        "no automatic branch repair was applied. USER_ACTION_REQUIRED."
      ), call. = FALSE)
    }
    nb_node <- .fast_ace_topology_guard(phy)
  }
  postorder <- ape::reorder.phylo(phy, "postorder")
  edge <- matrix(as.integer(postorder$edge), ncol = 2L)
  edge_length <- as.numeric(postorder$edge.length)
  workspace <- list(
    tip_label = postorder$tip.label,
    nb_tip = length(postorder$tip.label),
    nb_node = as.integer(nb_node),
    edge = edge,
    edge_length = edge_length
  )
  if (isTRUE(include_fingerprint)) {
    workspace$fingerprint <- paste(
      c(postorder$tip.label, as.integer(edge), sprintf("%.17g", edge_length)),
      collapse = "|"
    )
  }
  workspace
}

.fast_ace_workspace <- function(tree_workspace, model, nl) {
  rate <- .fast_ace_rate_index(model, nl)
  c(tree_workspace, list(rate = rate, model = model, nl = as.integer(nl)))
}

.fast_ace_fit_workspace <- function(tip_state, lvls, workspace, ip,
                                    CI = TRUE, marginal = FALSE,
                                    estimate_se = TRUE,
                                    reuse_buffers = !isTRUE(estimate_se)) {
  rate <- workspace$rate
  par0 <- rep(ip, length.out = rate$np)
  dev <- function(p) {
    fast_ace_discrete_deviance_cpp(
      edge = workspace$edge, edge_length = workspace$edge_length,
      tip_state = tip_state, rate_index = rate$rate_index, par = p,
      reuse_buffers = isTRUE(reuse_buffers)
    )
  }

  fit <- stats::nlminb(
    start = par0, objective = dev,
    lower = rep(0, rate$np), upper = rep(1e50, rate$np)
  )

  rate_se <- rep(NaN, rate$np)
  if (isTRUE(estimate_se)) {
    hessian_fit <- try(
      suppressWarnings(stats::nlm(
        f = dev, p = fit$par, iterlim = 1, stepmax = 0, hessian = TRUE
      )),
      silent = TRUE
    )
    rate_se <- if (inherits(hessian_fit, "try-error") ||
        any(diag(hessian_fit$hessian) == 0)) {
      rep(NaN, rate$np)
    } else {
      solved <- try(solve(hessian_fit$hessian), silent = TRUE)
      if (inherits(solved, "try-error")) {
        rep(NaN, rate$np)
      } else {
        sqrt(diag(solved))
      }
    }
  }

  obj <- list(
    loglik = -0.5 * fit$objective,
    rates = fit$par,
    se = rate_se,
    index.matrix = rate$index_matrix,
    convergence = fit$convergence,
    message = fit$message
  )
  if (fit$convergence != 0L) {
    warning("transition-rate optimization did not converge: ", fit$message,
            call. = FALSE)
  }
  if (isTRUE(CI)) {
    lik <- fast_ace_discrete_liks_cpp(
      edge = workspace$edge, edge_length = workspace$edge_length,
      tip_state = tip_state, rate_index = rate$rate_index, par = fit$par,
      marginal = marginal, reuse_buffers = isTRUE(reuse_buffers)
    )
    lik_anc <- lik$lik.anc
    rownames(lik_anc) <- workspace$nb_tip + seq_len(workspace$nb_node)
    colnames(lik_anc) <- lvls
    obj$lik.anc <- lik_anc
  }
  obj
}

.fast_ace_lik_anc_workspace <- function(y, workspace, ip = 0.1,
                                        reuse_buffers = TRUE) {
  y <- y[workspace$tip_label]
  if (!is.factor(y) || nlevels(y) != workspace$nl) {
    stop("ACE workspace does not match the categorical trait.", call. = FALSE)
  }
  .fast_ace_fit_workspace(
    tip_state = as.integer(y), lvls = levels(y), workspace = workspace,
    ip = ip, CI = TRUE, marginal = FALSE, estimate_se = FALSE,
    reuse_buffers = reuse_buffers
  )$lik.anc
}

.fast_ace_topology_guard <- function(phy) {
  # The ACE C++ kernel consumes a postorder binary edge matrix and indexes
  # every internal node as n_tip + 1, ..., n_tip + Nnode.  Validate the full
  # rooted-tree contract here so malformed edge lists fail at the R boundary
  # with an actionable diagnostic instead of being interpreted as a binary
  # tree by the pairwise pruning loop.
  fail <- function(detail) {
    stop(
      paste0("\"phy\" must be rooted and fully dichotomous: ", detail),
      call. = FALSE
    )
  }

  tip_label <- phy$tip.label
  n_tip <- if (is.null(tip_label)) 0L else length(tip_label)
  if (length(n_tip) != 1L || !is.finite(n_tip) || n_tip < 2L) {
    fail("the tree must contain at least two tips.")
  }
  n_tip <- as.integer(n_tip)

  declared_nnode <- phy$Nnode
  declared_ok <- (is.numeric(declared_nnode) || is.integer(declared_nnode)) &&
    length(declared_nnode) == 1L && is.finite(declared_nnode) &&
    declared_nnode >= 1 && declared_nnode <= .Machine$integer.max &&
    declared_nnode == floor(declared_nnode)
  if (!declared_ok) {
    fail("Nnode must be one positive integer matching the edge topology.")
  }
  n_node <- as.integer(declared_nnode)
  if (n_node != n_tip - 1L) {
    fail(sprintf("Nnode (%d) must equal n_tip - 1 (%d).", n_node,
                 n_tip - 1L))
  }

  edge <- phy$edge
  edge_ok <- is.matrix(edge) && length(dim(edge)) == 2L &&
    ncol(edge) == 2L && nrow(edge) > 0L &&
    (is.numeric(edge) || is.integer(edge))
  if (!edge_ok) {
    fail("edge must be a two-column numeric/integer matrix.")
  }
  edge_values <- suppressWarnings(as.numeric(edge))
  if (nrow(edge) != 2L * n_node || any(!is.finite(edge_values)) ||
      any(edge_values != floor(edge_values))) {
    fail(sprintf("edge must contain exactly %d finite integer rows for a binary tree.",
                 2L * n_node))
  }
  edge <- matrix(as.integer(edge_values), nrow = nrow(edge), ncol = 2L)

  n_total <- n_tip + n_node
  parent <- edge[, 1L]
  child <- edge[, 2L]
  if (any(parent < 1L | parent > n_total |
          child < 1L | child > n_total)) {
    fail(sprintf("edge node IDs must be canonical integers in 1:%d.",
                 n_total))
  }

  indegree <- tabulate(child, nbins = n_total)
  outdegree <- tabulate(parent, nbins = n_total)
  tip_ids <- seq_len(n_tip)
  internal_ids <- seq.int(n_tip + 1L, n_total)

  if (any(parent %in% tip_ids) || any(outdegree[tip_ids] != 0L)) {
    fail("tips must never occur as parents.")
  }
  if (any(indegree[tip_ids] != 1L)) {
    fail("each tip must occur exactly once as a child.")
  }
  if (any(parent == child)) {
    fail("edge topology contains a self-loop.")
  }
  if (any(indegree[internal_ids] > 1L)) {
    fail("each non-root internal node must have indegree one.")
  }

  root_candidates <- internal_ids[indegree[internal_ids] == 0L]
  if (length(root_candidates) != 1L) {
    fail("topology must contain exactly one internal structural root.")
  }
  root <- root_candidates[[1L]]
  if (root != n_tip + 1L) {
    fail(sprintf("the structural root must use canonical node ID n_tip + 1 (%d).",
                 n_tip + 1L))
  }
  nonroot_internal <- internal_ids[internal_ids != root]
  if (length(nonroot_internal) && any(indegree[nonroot_internal] != 1L)) {
    fail("each non-root internal node must have indegree one.")
  }
  if (any(outdegree[internal_ids] != 2L)) {
    bad <- internal_ids[outdegree[internal_ids] != 2L]
    fail(sprintf("every internal node must have exactly two children (invalid node IDs: %s).",
                 paste(bad, collapse = ", ")))
  }

  # With one parent per non-root node, a DFS from the unique root detects both
  # disconnected components and cycles/repeated visits.  Keep this check
  # explicit because edge count and degree checks alone do not establish a
  # connected acyclic rooted tree for arbitrary malformed inputs.
  children <- split(child, parent)
  seen <- rep(FALSE, n_total)
  stack <- root
  seen[[root]] <- TRUE
  while (length(stack)) {
    node <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    kids <- children[[as.character(node)]]
    if (is.null(kids)) next
    for (kid in as.integer(kids)) {
      if (seen[[kid]]) {
        fail("edge topology must be connected and acyclic.")
      }
      seen[[kid]] <- TRUE
      stack <- c(stack, kid)
    }
  }
  if (!all(seen)) {
    fail("edge topology must be connected and acyclic.")
  }

  n_node
}

.fast_ace_rate_index <- function(model, nl) {
  rate_index <- matrix(0L, nl, nl)
  if (model == "ER") {
    np <- 1L
    rate_index[col(rate_index) != row(rate_index)] <- 1L
  } else if (model == "ARD") {
    np <- nl * (nl - 1L)
    rate_index[col(rate_index) != row(rate_index)] <- seq_len(np)
  } else {
    stop("fast_ace currently supports only model = \"ER\" or \"ARD\".",
         call. = FALSE)
  }

  index_matrix <- rate_index
  diag(index_matrix) <- NA_integer_
  list(rate_index = rate_index, index_matrix = index_matrix, np = np)
}
