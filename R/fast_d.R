fast_d <- function(tree, x = NULL, test = TRUE, nsim = 1000,
                   se = NULL, return_sim = test, verbose = TRUE,
                   rnd.bias = NULL, permutations = NULL,
                   random_states = NULL, brownian_states = NULL,
                   ncores = 1, X = NULL, permut = NULL) {
  if (is.null(x)) x <- X
  if (is.null(x)) {
    stop("x must be a named binary vector, matrix, or data.frame.",
         call. = FALSE)
  }
  ncores <- .normalize_ncores(ncores)
  if (!is.logical(test) || length(test) != 1L || is.na(test)) {
    stop("test must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(return_sim) || length(return_sim) != 1L ||
      is.na(return_sim)) {
    stop("return_sim must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(permut)) nsim <- permut
  if (!is.null(se)) {
    stop("fast_d does not use se; provide binary states in x.",
         call. = FALSE)
  }
  if (!inherits(tree, "phylo")) {
    stop("tree should be an object of class \"phylo\".", call. = FALSE)
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim < 1 || nsim != floor(nsim)) {
    stop("nsim/permut must be a positive integer.", call. = FALSE)
  }
  nsim <- as.integer(nsim)

  vector_input <- is.null(dim(x))
  X <- .as_binary_table(x, tree)

  tree_tips <- tree$tip.label
  matched <- tree_tips[tree_tips %in% rownames(X)]
  removed_tree <- setdiff(tree_tips, rownames(X))
  removed_data <- setdiff(rownames(X), tree_tips)
  if (length(matched) < 2) {
    stop("Fewer than two species are shared by tree and x.", call. = FALSE)
  }

  base_tree <- if (length(removed_tree)) {
    ape::drop.tip(tree, removed_tree)
  } else {
    tree
  }
  X0 <- X[base_tree$tip.label, , drop = FALSE]

  if (isTRUE(verbose)) {
    message(sprintf(
      paste(
        "Matched %d species; removed %d tree tips and %d data rows before",
        "trait-wise NA pruning."
      ),
      length(matched), length(removed_tree), length(removed_data)
    ))
  }

  bias0 <- NULL
  if (!is.null(rnd.bias)) {
    if (is.null(names(rnd.bias))) {
      if (length(rnd.bias) != length(tree_tips)) {
        stop("rnd.bias must be named or have one value per tree tip.",
             call. = FALSE)
      }
      names(rnd.bias) <- tree_tips
    }
    bias0 <- as.numeric(rnd.bias[base_tree$tip.label])
    names(bias0) <- base_tree$tip.label
    if (anyNA(bias0)) {
      stop("rnd.bias is missing values for matched species.", call. = FALSE)
    }
  }

  D_fast <- Pval1_fast <- Pval0_fast <- observed <- mean_random <-
    mean_brownian <- rep(NA_real_, ncol(X0))
  n_species <- n_removed_na <- rep(NA_integer_, ncol(X0))
  notes <- rep(NA_character_, ncol(X0))
  state_tables <- vector("list", ncol(X0))
  random_list <- brownian_list <- vector("list", ncol(X0))

  for (j in seq_len(ncol(X0))) {
    keep <- !is.na(X0[, j])
    n_species[j] <- sum(keep)
    n_removed_na[j] <- sum(!keep)
    if (n_species[j] < 2) {
      notes[j] <- "fewer than 2 non-NA matched species"
      next
    }

    group_tree <- if (all(keep)) {
      base_tree
    } else {
      ape::drop.tip(base_tree, base_tree$tip.label[!keep])
    }
    d <- .binary_state(stats::setNames(X0[keep, j], rownames(X0)[keep]))
    if (!is.null(d$error)) {
      notes[j] <- d$error
      next
    }

    phy <- ape::reorder.phylo(group_tree, "pruningwise")
    edge <- matrix(as.integer(phy$edge), ncol = 2)
    edge_length <- as.numeric(phy$edge.length)
    n_tip <- length(phy$tip.label)
    .check_phylo_d_branch_lengths(phy)

    ds <- d$values[phy$tip.label]
    bias <- if (is.null(bias0)) NULL else as.numeric(bias0[phy$tip.label])

    rand <- .phylo_d_random_states(
      ds = ds, nsim = nsim, rnd.bias = bias, permutations = permutations,
      random_states = random_states
    )
    brown <- .phylo_d_brownian_states(
      tree = group_tree, nsim = nsim, prop_state1 = d$prop_state1,
      brownian_states = brownian_states
    )

    ran_sums <- phylo_d_sums_cpp(
      cbind(Obs = ds, rand), edge = edge, edge_length = edge_length,
      n_tip = n_tip, n_threads = ncores
    )
    phy_sums <- phylo_d_sums_cpp(
      cbind(Obs = ds, brown), edge = edge, edge_length = edge_length,
      n_tip = n_tip, n_threads = ncores
    )

    if (round(ran_sums[[1]], 6) != round(phy_sums[[1]], 6)) {
      stop("Problem with character change calculation in fast_d.",
           call. = FALSE)
    }

    rans <- ran_sums[-1]
    phys <- phy_sums[-1]
    observed[j] <- ran_sums[[1]]
    mean_random[j] <- mean(rans)
    mean_brownian[j] <- mean(phys)
    D_fast[j] <- (observed[j] - mean_brownian[j]) /
      (mean_random[j] - mean_brownian[j])
    Pval1_fast[j] <- sum(rans < observed[j]) / nsim
    Pval0_fast[j] <- sum(phys > observed[j]) / nsim
    state_tables[[j]] <- d$states_table
    if (return_sim) {
      random_list[[j]] <- rans
      brownian_list[[j]] <- phys
    }
  }

  out <- data.frame(
    trait = colnames(X0),
    D_fast = D_fast,
    Pval1_fast = if (test) Pval1_fast else NA_real_,
    Pval0_fast = if (test) Pval0_fast else NA_real_,
    observed = observed,
    mean_random = mean_random,
    mean_brownian = mean_brownian,
    n_species = n_species,
    n_removed_na = n_removed_na,
    matched_species = length(matched),
    removed_tree_tips = length(removed_tree),
    removed_data_rows = length(removed_data),
    note = notes,
    stringsAsFactors = FALSE
  )
  if (return_sim) {
    out$random_fast <- I(random_list)
    out$brownian_fast <- I(brownian_list)
  }

  attr(out, "match_report") <- list(
    original_tree_tips = ape::Ntip(tree),
    input_rows = nrow(X),
    matched_species = length(matched),
    removed_tree_tips = length(removed_tree),
    removed_data_rows = length(removed_data),
    tree_tips_removed = removed_tree,
    data_rows_removed = removed_data
  )

  if (!vector_input) {
    return(out)
  }

  object <- list(
    DEstimate = out$D_fast[[1]],
    Pval1 = out$Pval1_fast[[1]],
    Pval0 = out$Pval0_fast[[1]],
    Parameters = list(
      Observed = out$observed[[1]],
      MeanRandom = out$mean_random[[1]],
      MeanBrownian = out$mean_brownian[[1]]
    ),
    StatesTable = state_tables[[1]],
    binvar = out$trait[[1]],
    data = list(data.name = deparse(substitute(x)),
                phy.name = deparse(substitute(tree))),
    nPermut = nsim,
    rnd.bias = rnd.bias
  )
  if (return_sim) {
    object$Permutations <- list(
      random = out$random_fast[[1]],
      brownian = out$brownian_fast[[1]]
    )
  }
  class(object) <- "phylo.d"
  object
}

.as_binary_table <- function(x, tree) {
  if (is.data.frame(x)) {
    rn <- rownames(x)
    out <- as.data.frame(x, stringsAsFactors = FALSE)
    rownames(out) <- rn
  } else if (is.null(dim(x))) {
    out <- data.frame(x = I(x), stringsAsFactors = FALSE)
    rownames(out) <- names(x)
  } else {
    out <- as.data.frame(x, stringsAsFactors = FALSE)
  }

  default_rownames <- identical(rownames(out), as.character(seq_len(nrow(out))))
  if ((is.null(rownames(out)) || default_rownames) &&
      nrow(out) == ape::Ntip(tree)) {
    message("x has no names; assuming x is in tree$tip.label order")
    rownames(out) <- tree$tip.label
    default_rownames <- FALSE
  }
  if (is.null(rownames(out)) || default_rownames) {
    stop("x must have species names as names/rownames.", call. = FALSE)
  }
  if (anyDuplicated(rownames(out))) {
    stop("x species names must be unique.", call. = FALSE)
  }
  if (is.null(colnames(out))) {
    colnames(out) <- paste0("trait_", seq_len(ncol(out)))
  }
  out
}

.binary_state <- function(x) {
  x <- x[!is.na(x)]
  u <- unique(x)
  if (length(u) > 2) {
    return(list(error = "contains more than two states"))
  }
  if (length(u) < 2) {
    return(list(error = "contains a single state"))
  }
  # D must be invariant to labels such as 0/1, 1/2, or 10/20. Recode every
  # valid binary trait to equally spaced states before comparing it with the
  # simulated 0/1 Brownian null.
  x_factor <- droplevels(factor(x))
  states_table <- unclass(table(x_factor))
  prop_state1 <- states_table[[1]] / sum(states_table)
  values <- as.numeric(x_factor) - 1
  names(values) <- names(x)
  list(
    values = values,
    states_table = states_table,
    prop_state1 = prop_state1,
    error = NULL
  )
}

.check_phylo_d_branch_lengths <- function(phy) {
  el <- phy$edge.length
  if (is.null(el)) {
    stop("fast_d requires a tree with branch lengths.", call. = FALSE)
  }
  if (any(!is.finite(el)) || any(el < 0)) {
    stop("fast_d requires finite non-negative branch lengths.",
         call. = FALSE)
  }
  tip_edge <- phy$edge[, 2] <= length(phy$tip.label)
  if (any(el[tip_edge] == 0)) {
    stop("Phylogeny contains pairs of tips on zero branch lengths.",
         call. = FALSE)
  }
  if (any(el[!tip_edge] == 0)) {
    stop("Phylogeny contains zero length internal branches. Use di2multi.",
         call. = FALSE)
  }
}

.phylo_d_random_states <- function(ds, nsim, rnd.bias, permutations,
                                   random_states) {
  n <- length(ds)
  if (!is.null(random_states)) {
    random_states <- as.matrix(random_states)
    storage.mode(random_states) <- "double"
    if (nrow(random_states) != n || ncol(random_states) != nsim) {
      stop("random_states must be an n_species x nsim matrix.",
           call. = FALSE)
    }
    return(random_states)
  }
  if (!is.null(permutations)) {
    permutations <- .permutation_matrix(permutations, n = n, nsim = nsim)
    return(matrix(ds[t(permutations)], nrow = n, ncol = nsim))
  }
  replicate(nsim, sample(ds, prob = rnd.bias))
}

.phylo_d_brownian_states <- function(tree, nsim, prop_state1,
                                     brownian_states) {
  n <- length(tree$tip.label)
  if (!is.null(brownian_states)) {
    brownian_states <- as.matrix(brownian_states)
    storage.mode(brownian_states) <- "double"
    if (nrow(brownian_states) != n || ncol(brownian_states) != nsim) {
      stop("brownian_states must be an n_species x nsim matrix.",
           call. = FALSE)
    }
    return(brownian_states)
  }
  phy <- ape::reorder.phylo(tree, "cladewise")
  edge <- matrix(as.integer(phy$edge), ncol = 2)
  brownian_tree_threshold_cpp(
    edge = edge,
    edge_length = as.numeric(phy$edge.length),
    n_tip = n,
    nsim = nsim,
    prop_state1 = prop_state1
  )
}
