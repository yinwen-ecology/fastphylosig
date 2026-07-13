# Delta statistic for categorical traits --------------------------------------

fast_delta <- function(tree, x = NULL, test = FALSE, nsim = 1000,
                       se = NULL, mcmc_sim = 10000, thin = 10, burn = 100,
                       lambda0 = 0.1, proposal_sd = 0.5,
                       entropy = c("LSE", "SE", "GINI"), model = "ARD",
                       ace_engine = c("fast", "ape"),
                       permutations = NULL, return_sim = test,
                       verbose = TRUE, ncores = 1, X = NULL) {
  if (is.null(x)) x <- X
  if (is.null(x)) {
    stop("x must be a named categorical vector, matrix, or data.frame.",
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
  if (!is.null(se)) {
    stop("fast_delta does not use se; use proposal_sd for the MCMC proposal.",
         call. = FALSE)
  }
  if (!inherits(tree, "phylo")) {
    stop("tree should be an object of class \"phylo\".", call. = FALSE)
  }
  .check_delta_tree(tree)
  if (test && (!is.numeric(nsim) || length(nsim) != 1L ||
      !is.finite(nsim) || nsim < 1 || nsim != floor(nsim))) {
    stop("nsim must be a positive integer when test = TRUE.", call. = FALSE)
  }
  nsim <- as.integer(nsim)
  mcmc_args <- c(mcmc_sim = mcmc_sim, thin = thin, burn = burn)
  if (!is.numeric(mcmc_args) || any(lengths(list(mcmc_sim, thin, burn)) != 1L) ||
      any(!is.finite(mcmc_args)) || any(mcmc_args != floor(mcmc_args))) {
    stop("mcmc_sim, thin, and burn must be finite integers.", call. = FALSE)
  }
  mcmc_sim <- as.integer(mcmc_sim)
  thin <- as.integer(thin)
  burn <- as.integer(burn)
  if (mcmc_sim < 1L || thin < 1L || burn < 1L || burn > mcmc_sim) {
    stop("mcmc_sim, thin, and burn must satisfy 1 <= burn <= mcmc_sim.",
         call. = FALSE)
  }
  if (!is.numeric(lambda0) || length(lambda0) != 1L ||
      !is.finite(lambda0) || !is.numeric(proposal_sd) ||
      length(proposal_sd) != 1L || !is.finite(proposal_sd) ||
      lambda0 <= 0 || proposal_sd <= 0) {
    stop("lambda0 and proposal_sd must be positive.", call. = FALSE)
  }
  entropy <- match.arg(entropy)
  ace_engine <- match.arg(ace_engine)
  model <- match.arg(model, c("ER", "ARD"))
  entropy_code <- match(entropy, c("LSE", "SE", "GINI"))

  vector_input <- is.null(dim(x))
  X <- .as_categorical_table(x, tree)

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

  Delta_fast <- P_fast <- alpha_mean <- beta_mean <- rep(NA_real_, ncol(X0))
  n_saved <- rep(NA_real_, ncol(X0))
  n_species <- n_removed_na <- n_states <- rep(NA_integer_, ncol(X0))
  n_failed_sim <- integer(ncol(X0))
  notes <- rep(NA_character_, ncol(X0))
  sim_list <- vector("list", ncol(X0))
  delta_cluster <- NULL

  for (j in seq_len(ncol(X0))) {
    keep <- !is.na(X0[[j]])
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
    y <- X0[[j]][keep]
    names(y) <- rownames(X0)[keep]
    y <- factor(y[group_tree$tip.label])
    n_states[j] <- length(levels(y))
    if (n_states[j] < 2) {
      notes[j] <- "contains a single state"
      next
    }

    calc <- tryCatch(
      .delta_one(
        tree = group_tree, y = y, lambda0 = lambda0,
        proposal_sd = proposal_sd, mcmc_sim = mcmc_sim, thin = thin,
        burn = burn, entropy_code = entropy_code, model = model,
        ace_engine = ace_engine
      ),
      error = function(e) e
    )
    if (inherits(calc, "error")) {
      notes[j] <- conditionMessage(calc)
      next
    }
    Delta_fast[j] <- calc$delta
    alpha_mean[j] <- calc$alpha_mean
    beta_mean[j] <- calc$beta_mean
    n_saved[j] <- calc$n_saved

    if (test) {
      perms <- .permutation_matrix(
        permutations, n = length(y), nsim = nsim,
        include_observed = FALSE
      )
      if (is.null(delta_cluster)) {
        if (ncores > 1L) {
          delta_cluster <- parallel::makeCluster(min(ncores, nsim))
          on.exit(parallel::stopCluster(delta_cluster), add = TRUE)
          parallel::clusterEvalQ(delta_cluster, library(fastphylosig))
          parallel::clusterSetRNGStream(
            delta_cluster, sample.int(.Machine$integer.max, 1L)
          )
        }
      }
      if (is.null(delta_cluster)) {
        sim_delta <- vapply(
          seq_len(nrow(perms)), .delta_permutation_worker, numeric(1),
          y = y, perms = perms, tree = group_tree, lambda0 = lambda0,
          proposal_sd = proposal_sd, mcmc_sim = mcmc_sim, thin = thin,
          burn = burn, entropy_code = entropy_code, model = model,
          ace_engine = ace_engine
        )
      } else {
        sim_delta <- unlist(
          parallel::parLapply(
            delta_cluster, seq_len(nrow(perms)), .delta_permutation_worker,
            y = y, perms = perms, tree = group_tree, lambda0 = lambda0,
            proposal_sd = proposal_sd, mcmc_sim = mcmc_sim, thin = thin,
            burn = burn, entropy_code = entropy_code, model = model,
            ace_engine = ace_engine
          ),
          use.names = FALSE
        )
      }
      finite_sim <- is.finite(sim_delta)
      n_failed_sim[j] <- sum(!finite_sim)
      if (!any(finite_sim)) {
        notes[j] <- "all Delta permutations failed"
      } else {
        P_fast[j] <- mean(sim_delta[finite_sim] >= Delta_fast[j])
        if (n_failed_sim[j] > 0L) {
          notes[j] <- sprintf("%d Delta permutations failed", n_failed_sim[j])
        }
      }
      if (return_sim) {
        sim_list[[j]] <- sim_delta
      }
    }
  }

  out <- data.frame(
    trait = colnames(X0),
    Delta_fast = Delta_fast,
    alpha_mean = alpha_mean,
    beta_mean = beta_mean,
    n_saved = n_saved,
    n_states = n_states,
    n_species = n_species,
    n_removed_na = n_removed_na,
    matched_species = length(matched),
    removed_tree_tips = length(removed_tree),
    removed_data_rows = length(removed_data),
    note = notes,
    stringsAsFactors = FALSE
  )
  if (test) {
    out$P_fast <- P_fast
    out$n_failed_sim <- n_failed_sim
    if (return_sim) {
      out$sim.Delta_fast <- I(sim_list)
    }
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
  attr(out, "delta_control") <- list(
    lambda0 = lambda0,
    proposal_sd = proposal_sd,
    mcmc_sim = mcmc_sim,
    thin = thin,
    burn = burn,
    entropy = entropy,
    model = model,
    ace_engine = ace_engine,
    ncores = ncores
  )

  if (!vector_input) {
    return(out)
  }

  object <- list(
    delta = out$Delta_fast[[1]],
    alpha_mean = out$alpha_mean[[1]],
    beta_mean = out$beta_mean[[1]],
    n_saved = out$n_saved[[1]],
    n_states = out$n_states[[1]],
    n_species = out$n_species[[1]],
    parameters = attr(out, "delta_control")
  )
  if (test) {
    object$P <- out$P_fast[[1]]
    object$n_failed_sim <- out$n_failed_sim[[1]]
    if (return_sim) {
      object$sim.delta <- out$sim.Delta_fast[[1]]
    }
  }
  class(object) <- "phylo_delta"
  object
}

.delta_permutation_worker <- function(i, y, perms, tree, lambda0, proposal_sd,
                                      mcmc_sim, thin, burn, entropy_code,
                                      model, ace_engine) {
  yp <- y[perms[i, ]]
  names(yp) <- names(y)
  sim_calc <- tryCatch(
    .delta_one(
      tree = tree, y = yp, lambda0 = lambda0, proposal_sd = proposal_sd,
      mcmc_sim = mcmc_sim, thin = thin, burn = burn,
      entropy_code = entropy_code, model = model, ace_engine = ace_engine
    ),
    error = function(e) e
  )
  if (inherits(sim_calc, "error")) {
    return(NA_real_)
  }
  sim_calc$delta
}

.delta_one <- function(tree, y, lambda0, proposal_sd, mcmc_sim, thin, burn,
                       entropy_code, model, ace_engine) {
  y <- y[tree$tip.label]
  ar <- if (identical(ace_engine, "fast")) {
    fast_ace(
      y, tree, type = "discrete", method = "ML", model = model,
      CI = TRUE, marginal = FALSE
    )$lik.anc
  } else {
    ape::ace(
      y, tree, type = "discrete", method = "ML", model = model
    )$lik.anc
  }
  if (is.complex(ar)) {
    ar <- Re(ar)
  }
  ar <- as.matrix(ar)
  storage.mode(ar) <- "double"
  delta_mcmc_cpp(
    probabilities = ar, lambda0 = lambda0, proposal_sd = proposal_sd,
    sim = mcmc_sim, thin = thin, burn = burn, entropy_type = entropy_code
  )
}

.as_categorical_table <- function(x, tree) {
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

.check_delta_tree <- function(tree) {
  if (is.null(tree$edge.length)) {
    stop("Delta requires a tree with branch lengths.", call. = FALSE)
  }
  if (anyNA(tree$edge.length) || any(tree$edge.length <= 0)) {
    stop("Delta requires strictly positive branch lengths.", call. = FALSE)
  }
}
