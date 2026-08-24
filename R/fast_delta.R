# Delta statistic for categorical traits --------------------------------------

fast_delta <- function(tree, x = NULL, test = FALSE, nsim = 1000,
                       se = NULL, mcmc_sim = 10000, thin = 10, burn = 100,
                       lambda0 = 0.1, proposal_sd = 0.5,
                       entropy = c("LSE", "SE", "GINI"), model = "ARD",
                       permutations = NULL, return_sim = test,
                       verbose = TRUE, ncores = 1, X = NULL,
                       prepared = NULL, progress = interactive()) {
  if (!missing(x) && !missing(X)) {
    stop("Supply only one of x or X, not both.", call. = FALSE)
  }
  runtime <- .runtime_begin("Delta", progress = progress, verbose = verbose)
  on.exit(.runtime_on_exit(runtime), add = TRUE)
  verbose <- isTRUE(verbose) && isTRUE(progress)
  .runtime_stage(runtime, "Checking tree...")
  .runtime_stage(runtime, "Preparing data...")
  out <- .fast_delta_with_engine(
    tree = tree, x = x, test = test, nsim = nsim, se = se,
    mcmc_sim = mcmc_sim, thin = thin, burn = burn, lambda0 = lambda0,
    proposal_sd = proposal_sd, entropy = entropy, model = model,
    ace_engine = "fast", permutations = permutations,
    return_sim = return_sim, verbose = verbose, ncores = ncores, X = X,
    prepared = prepared, .runtime = runtime
  )
  timing <- .runtime_close(runtime, success = TRUE)
  .runtime_attach(out, timing)
}

.fast_delta_with_engine <- function(
    tree, x = NULL, test = FALSE, nsim = 1000, se = NULL,
    mcmc_sim = 10000, thin = 10, burn = 100, lambda0 = 0.1,
    proposal_sd = 0.5, entropy = c("LSE", "SE", "GINI"), model = "ARD",
    ace_engine = c("fast", "ape"), permutations = NULL,
    return_sim = test, verbose = TRUE, ncores = 1, X = NULL,
    prepared = NULL, .runtime = NULL) {
  if (!is.null(prepared)) tree <- prepared
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
  if (test && (!is.numeric(nsim) || length(nsim) != 1L ||
      !is.finite(nsim) || nsim < 1 || nsim != floor(nsim))) {
    stop("nsim must be a positive integer when test = TRUE.", call. = FALSE)
  }
  if (test) nsim <- as.integer(nsim)
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

  # Capture only non-sensitive stochastic metadata.  The raw .Random.seed is
  # intentionally never copied into a result object.
  rng_kind_start <- tryCatch(RNGkind(), error = function(e) character())
  seed_metadata_start <- .delta_seed_metadata()

  vector_input <- is.null(dim(x))
  # The shared preparation layer performs Delta's method-specific tree
  # preflight, safe representation canonicalisation, matching, and packed NA
  # grouping.  Keep the categorical conversion here for the historical input
  # contract and to retain the exact input-row metadata in the result.
  X <- .as_named_trait_table(
    x, if (inherits(tree, "fastphylosig_tree")) tree$tree else tree,
    verbose = FALSE, input_name = "x"
  )
  analysis <- .prepare_analysis(
    tree = tree, data = X, signal = "Delta", data_kind = "categorical",
    verbose = verbose
  )
  tree_check <- analysis$tree_processing$check_after
  if (!isTRUE(tree_check$ready_by_signal[["Delta"]])) {
    stop(.format_actionable_condition(tree_check, signal = "Delta"),
         call. = FALSE)
  }
  ctx <- analysis$ctx
  tree <- ctx$tree
  tree_tips <- tree$tip.label
  X0 <- analysis$matched_data
  matched <- rownames(X0)
  matching_details <- analysis$matching_details
  removed_tree <- if (is.list(matching_details)) {
    matching_details$tree_tips_removed
  } else {
    setdiff(tree_tips, matched)
  }
  removed_data <- if (is.list(matching_details)) {
    matching_details$data_rows_removed
  } else {
    setdiff(rownames(X), tree_tips)
  }
  if (length(matched) < 2) {
    stop("Fewer than two species are shared by tree and x.", call. = FALSE)
  }

  base_group <- analysis$base_group
  if (is.null(base_group)) {
    stop("Fewer than two species are shared by tree and x.", call. = FALSE)
  }
  base_tree <- base_group$tree
  X0 <- X0[base_tree$tip.label, , drop = FALSE]

  .runtime_stage(.runtime, "Estimating ancestral states...")
  .runtime_stage(.runtime, "Running MCMC...")
  if (isTRUE(test)) {
    .runtime_stage(.runtime, "Running permutations...")
  }

  if (isTRUE(verbose) && (length(removed_tree) || length(removed_data))) {
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
  alpha_sd <- beta_sd <- ESS_alpha <- ESS_beta <- rep(NA_real_, ncol(X0))
  split_Rhat_alpha <- split_Rhat_beta <- rep(NA_real_, ncol(X0))
  alpha_beta_cov <- MCSE_Delta <- rep(NA_real_, ncol(X0))
  diagnostics_available <- rep(FALSE, ncol(X0))
  diagnostics_note <- rep(NA_character_, ncol(X0))
  saved_per_chain_requested <- .delta_saved_iterations(mcmc_sim, thin, burn)
  n_saved_requested <- rep(
    as.numeric(2 * saved_per_chain_requested), ncol(X0)
  )
  n_saved_successful <- rep(0, ncol(X0))
  requested_iterations <- rep(as.numeric(2 * mcmc_sim), ncol(X0))
  successful_iterations <- n_saved_successful
  n_species <- n_removed_na <- n_states <- rep(NA_integer_, ncol(X0))
  n_failed_sim <- integer(ncol(X0))
  requested_simulations <- rep(if (test) as.integer(nsim) else NA_integer_,
                               ncol(X0))
  successful_simulations <- rep(NA_integer_, ncol(X0))
  requested_permutations <- requested_simulations
  successful_permutations <- successful_simulations
  notes <- rep(NA_character_, ncol(X0))
  nonfinite_trait <- vapply(X0, function(z) {
    is.numeric(z) && any(is.infinite(z) | is.nan(z))
  }, logical(1))
  if (any(nonfinite_trait)) {
    notes[nonfinite_trait] <-
      "trait contains non-finite values other than NA"
  }
  sim_list <- vector("list", ncol(X0))
  delta_cluster <- NULL
  # These structures are invariant to categorical-state permutations.  Cache
  # them locally so no state escapes into a result object or later call.
  mask_groups <- analysis$na_patterns
  ace_tree_workspaces <- vector("list", mask_groups$n_group)
  ace_workspaces <- new.env(parent = emptyenv())

  # Reuse the packed NA masks and validated structural groups produced by
  # preparation; no second mask grouping or tree matching pass is performed.
  group_id <- mask_groups$group_id
  if (is.null(group_id)) {
    group_id <- integer(ncol(X0))
    for (i in seq_len(mask_groups$n_group)) {
      group_id[mask_groups$columns[[i]]] <- i
    }
  }

  for (j in seq_len(ncol(X0))) {
    if (isTRUE(nonfinite_trait[[j]])) {
      n_species[j] <- sum(!is.na(X0[[j]]))
      n_removed_na[j] <- nrow(X0) - n_species[j]
      next
    }
    group_index <- group_id[[j]]
    keep <- mask_groups$keep[[group_index]]
    n_species[j] <- length(keep)
    n_removed_na[j] <- nrow(X0) - length(keep)
    if (n_species[j] < 2) {
      notes[j] <- "fewer than 2 non-NA matched species"
      next
    }

    prepared_group <- analysis$groups[[group_index]]
    if (is.null(prepared_group) || !isTRUE(prepared_group$kernel_ready)) {
      group_check <- if (is.list(prepared_group)) {
        prepared_group$check
      } else {
        NULL
      }
      notes[j] <- if (exists(".format_actionable_condition",
                             mode = "function", inherits = TRUE) &&
                     !is.null(group_check)) {
        .format_actionable_condition(
          group_check, signal = "Delta",
          prefix = "Delta retained subset is not ready"
        )
      } else {
        "Delta retained subset is not ready; run check_tree()"
      }
      next
    }
    group_tree <- prepared_group$group$tree
    y <- X0[[j]][keep]
    if (is.numeric(y) && any(!is.finite(y[!is.na(y)]))) {
      notes[j] <- "trait contains non-finite values other than NA"
      next
    }
    names(y) <- rownames(X0)[keep]
    y <- factor(y[group_tree$tip.label])
    n_states[j] <- length(levels(y))
    if (n_states[j] < 2) {
      notes[j] <- "contains a single state"
      next
    }

    ace_workspace <- NULL
    if (identical(ace_engine, "fast")) {
      tree_workspace <- ace_tree_workspaces[[group_index]]
      if (is.null(tree_workspace)) {
        tree_workspace <- .fast_ace_tree_workspace(
          group_tree, include_fingerprint = TRUE
        )
        ace_tree_workspaces[[group_index]] <- tree_workspace
      }
      # The cache is local to this call and group_index identifies the
      # retained subtree.  Keep the full fingerprint in the workspace for
      # auditability, but do not use it as an environment key: large trees
      # can exceed R's 10,000-byte variable-name limit.
      workspace_key <- paste(group_index, model, n_states[j], sep = ":")
      ace_workspace <- get0(workspace_key, envir = ace_workspaces,
                            inherits = FALSE)
      if (is.null(ace_workspace)) {
        ace_workspace <- .fast_ace_workspace(
          tree_workspace, model = model, nl = n_states[j]
        )
        assign(workspace_key, ace_workspace, envir = ace_workspaces)
      }
    }

    calc <- tryCatch(
      .delta_one(
        tree = group_tree, y = y, lambda0 = lambda0,
        proposal_sd = proposal_sd, mcmc_sim = mcmc_sim, thin = thin,
        burn = burn, entropy_code = entropy_code, model = model,
        ace_engine = ace_engine, keep_chains = TRUE,
        ace_workspace = ace_workspace
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
    n_saved_successful[j] <- if (is.finite(calc$n_saved)) {
      as.numeric(calc$n_saved)
    } else {
      0
    }
    successful_iterations[j] <- as.numeric(2 * mcmc_sim)

    diag <- .delta_mcmc_diagnostics(calc)
    alpha_sd[j] <- diag$alpha_sd
    beta_sd[j] <- diag$beta_sd
    ESS_alpha[j] <- diag$ESS_alpha
    ESS_beta[j] <- diag$ESS_beta
    split_Rhat_alpha[j] <- diag$split_Rhat_alpha
    split_Rhat_beta[j] <- diag$split_Rhat_beta
    alpha_beta_cov[j] <- diag$alpha_beta_cov
    MCSE_Delta[j] <- diag$MCSE_Delta
    diagnostics_available[j] <- isTRUE(diag$available)
    diagnostics_note[j] <- diag$note

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
          ace_engine = ace_engine, ace_workspace = ace_workspace
        )
      } else {
        sim_delta <- unlist(
          parallel::parLapply(
            delta_cluster, seq_len(nrow(perms)), .delta_permutation_worker,
            y = y, perms = perms, tree = group_tree, lambda0 = lambda0,
            proposal_sd = proposal_sd, mcmc_sim = mcmc_sim, thin = thin,
            burn = burn, entropy_code = entropy_code, model = model,
            ace_engine = ace_engine, ace_workspace = ace_workspace
          ),
          use.names = FALSE
        )
      }
      finite_sim <- is.finite(sim_delta)
      n_failed_sim[j] <- sum(!finite_sim)
      successful_simulations[j] <- sum(finite_sim)
      successful_permutations[j] <- successful_simulations[j]
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

  .runtime_stage(.runtime, "Finalizing diagnostics...")
  status <- rep("failed", ncol(X0))
  status[is.finite(Delta_fast)] <- "ok"
  status[grepl("fewer than 2", notes, fixed = TRUE) %in% TRUE] <-
    "insufficient_data"
  status[grepl("single state|non-finite", notes, ignore.case = TRUE) %in% TRUE] <-
    "invalid_trait"
  if (isTRUE(test)) {
    partial <- is.finite(Delta_fast) &
      (!is.finite(successful_simulations) |
         successful_simulations < requested_simulations)
    status[partial %in% TRUE] <- "partial"
  }

  small_retained <- which(
    is.finite(n_species) & n_species >= 2L & n_species < 20L
  )
  sample_size_status <- ifelse(
    n_species < 2L, "error",
    ifelse(n_species < 20L, "warning", "ok")
  )
  delta_retained_species <- list(
    warning_threshold = 20L,
    warning_emitted = length(small_retained) > 0L,
    by_trait = data.frame(
      trait = colnames(X0),
      n_species = as.integer(n_species),
      n_removed_na = as.integer(n_removed_na),
      sample_size_status = sample_size_status,
      stringsAsFactors = FALSE
    ),
    affected_traits = data.frame(
      trait = colnames(X0)[small_retained],
      n_species = as.integer(n_species[small_retained]),
      n_removed_na = as.integer(n_removed_na[small_retained]),
      stringsAsFactors = FALSE
    )
  )
  if (length(small_retained)) {
    retained_detail <- paste(
      sprintf(
        "%s (n = %d)", colnames(X0)[small_retained],
        as.integer(n_species[small_retained])
      ),
      collapse = ", "
    )
    warning(
      sprintf(
        paste(
          "Delta retained 2 to 19 species after matching and trait-specific",
          "NA pruning for %d trait%s: %s. Estimates are returned, but",
          "should be interpreted cautiously."
        ),
        length(small_retained), if (length(small_retained) == 1L) "" else "s",
        retained_detail
      ),
      call. = FALSE
    )
  }

  # A finite estimate is still returned when its MCMC summary is unreliable.
  # This is deliberately post-sampling: warnings must never trigger a hidden
  # rerun or change the caller's requested iteration count.
  delta_mcmc_diagnostics <- .delta_mcmc_warning_summary(
    trait = colnames(X0), delta = Delta_fast, ess_alpha = ESS_alpha,
    ess_beta = ESS_beta, rhat_alpha = split_Rhat_alpha,
    rhat_beta = split_Rhat_beta, mcse_delta = MCSE_Delta
  )
  affected_mcmc <- which(delta_mcmc_diagnostics$by_trait$warning)
  if (length(affected_mcmc)) {
    mcmc_detail <- paste(
      sprintf(
        "%s (%s)",
        delta_mcmc_diagnostics$by_trait$trait[affected_mcmc],
        delta_mcmc_diagnostics$by_trait$reason[affected_mcmc]
      ),
      collapse = ", "
    )
    warning(
      sprintf(
        paste(
          "Delta MCMC diagnostics require inspection for %d trait%s: %s.",
          "Estimates are returned unchanged; no MCMC chains were rerun.",
          "Inspect ESS, split R-hat, and MCSE_Delta before interpretation."
        ),
        length(affected_mcmc),
        if (length(affected_mcmc) == 1L) "" else "s",
        mcmc_detail
      ),
      call. = FALSE
    )
  }

  out <- data.frame(
    trait = colnames(X0),
    Delta_fast = Delta_fast,
    alpha_mean = alpha_mean,
    beta_mean = beta_mean,
    n_saved = n_saved,
    n_saved_requested = n_saved_requested,
    n_saved_successful = n_saved_successful,
    alpha_sd = alpha_sd,
    beta_sd = beta_sd,
    ESS_alpha = ESS_alpha,
    ESS_beta = ESS_beta,
    split_Rhat_alpha = split_Rhat_alpha,
    split_Rhat_beta = split_Rhat_beta,
    Rhat_alpha = split_Rhat_alpha,
    Rhat_beta = split_Rhat_beta,
    alpha_beta_cov = alpha_beta_cov,
    MCSE_Delta = MCSE_Delta,
    diagnostics_available = diagnostics_available,
    diagnostics_note = diagnostics_note,
    requested_simulations = requested_simulations,
    successful_simulations = successful_simulations,
    requested_iterations = requested_iterations,
    successful_iterations = successful_iterations,
    requested_permutations = requested_permutations,
    successful_permutations = successful_permutations,
    n_states = n_states,
    n_species = n_species,
    n_removed_na = n_removed_na,
    matched_species = length(matched),
    removed_tree_tips = length(removed_tree),
    removed_data_rows = length(removed_data),
    status = status,
    note = notes,
    message = notes,
    stringsAsFactors = FALSE
  )
  if (test) {
    out$P_fast <- P_fast
    out$P_MCSE <- vapply(seq_len(nrow(out)), function(j) {
      if (!is.finite(P_fast[j])) return(NA_real_)
      .delta_p_mcse(P_fast[j], successful_simulations[j])
    }, numeric(1))
    # Lower-case alias follows the existing lower-case R result conventions.
    out$P_mcse <- out$P_MCSE
    out$MCSE_P <- out$P_MCSE
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
  analysis_metadata <- if (exists(".compact_analysis_metadata",
                                  mode = "function", inherits = TRUE)) {
    .compact_analysis_metadata(analysis)
  } else {
    NULL
  }
  if (is.list(analysis_metadata)) {
    analysis_metadata$delta_retained_species <- delta_retained_species
    analysis_metadata$delta_mcmc_diagnostics <- delta_mcmc_diagnostics
    attr(out, "analysis_metadata") <- analysis_metadata
    attr(out, "match_report")$analysis_metadata <- analysis_metadata
  }
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

  rng_kind_end <- tryCatch(RNGkind(), error = function(e) character())
  seed_metadata_end <- .delta_seed_metadata()
  successful_saved_total <- sum(n_saved_successful[is.finite(n_saved_successful)])
  successful_permutation_total <- if (test) {
    sum(successful_simulations[is.finite(successful_simulations)])
  } else {
    0L
  }
  seed_summary <- seed_metadata_end
  seed_summary$start <- seed_metadata_start
  seed_summary$end <- seed_metadata_end
  seed_summary$raw_exposed <- FALSE
  attr(out, "delta_stochastic") <- list(
    rng_kind = rng_kind_end,
    rng_kind_start = rng_kind_start,
    rng_kind_end = rng_kind_end,
    seed = seed_summary,
    seed_metadata = list(
      start = seed_metadata_start,
      end = seed_metadata_end,
      raw_exposed = FALSE
    ),
    mcmc_sim_requested = as.integer(mcmc_sim),
    chains_requested = 2L,
    saved_iterations_requested = as.numeric(2 * saved_per_chain_requested),
    saved_iterations_successful = successful_saved_total,
    requested_iterations = as.numeric(2 * mcmc_sim),
    successful_iterations = as.numeric(
      2 * mcmc_sim * sum(is.finite(Delta_fast))
    ),
    requested_simulations = if (test) as.integer(nsim) else NA_integer_,
    successful_simulations = successful_permutation_total,
    requested_permutations = if (test) as.integer(nsim) else NA_integer_,
    successful_permutations = successful_permutation_total,
    diagnostics_available = any(diagnostics_available)
  )
  attr(out, "stochastic") <- attr(out, "delta_stochastic")
  attr(out, "rng_kind") <- rng_kind_end
  attr(out, "seed_metadata") <- list(
    start = seed_metadata_start,
    end = seed_metadata_end,
    raw_exposed = FALSE
  )

  if (!vector_input) {
    object <- .decorate_fastphylosig_result(
      out, method = "Delta", vector_input = FALSE
    )
    if (is.list(analysis_metadata)) {
      attr(object, "analysis_metadata") <- analysis_metadata
    }
    timing <- .runtime_close(.runtime, success = TRUE)
    return(.attach_fast_signal_workflow(
      .runtime_attach(object, timing), "Delta", "fast_delta"
    ))
  }

  object <- list(
    delta = out$Delta_fast[[1]],
    alpha_mean = out$alpha_mean[[1]],
    beta_mean = out$beta_mean[[1]],
    n_saved = out$n_saved[[1]],
    n_saved_requested = out$n_saved_requested[[1]],
    n_saved_successful = out$n_saved_successful[[1]],
    alpha_sd = out$alpha_sd[[1]],
    beta_sd = out$beta_sd[[1]],
    ESS_alpha = out$ESS_alpha[[1]],
    ESS_beta = out$ESS_beta[[1]],
    split_Rhat_alpha = out$split_Rhat_alpha[[1]],
    split_Rhat_beta = out$split_Rhat_beta[[1]],
    Rhat_alpha = out$Rhat_alpha[[1]],
    Rhat_beta = out$Rhat_beta[[1]],
    alpha_beta_cov = out$alpha_beta_cov[[1]],
    MCSE_Delta = out$MCSE_Delta[[1]],
    diagnostics_available = out$diagnostics_available[[1]],
    diagnostics_note = out$diagnostics_note[[1]],
    requested_simulations = out$requested_simulations[[1]],
    successful_simulations = out$successful_simulations[[1]],
    requested_iterations = out$requested_iterations[[1]],
    successful_iterations = out$successful_iterations[[1]],
    requested_permutations = out$requested_permutations[[1]],
    successful_permutations = out$successful_permutations[[1]],
    n_states = out$n_states[[1]],
    n_species = out$n_species[[1]],
    status = out$status[[1]],
    note = out$note[[1]],
    message = out$message[[1]],
    parameters = attr(out, "delta_control")
  )
  object$stochastic <- attr(out, "delta_stochastic")
  object$rng_kind <- attr(out, "rng_kind")
  object$seed_metadata <- attr(out, "seed_metadata")
  if (test) {
    object$P <- out$P_fast[[1]]
    object$P_MCSE <- out$P_MCSE[[1]]
    object$P_mcse <- out$P_mcse[[1]]
    object$MCSE_P <- out$MCSE_P[[1]]
    object$n_failed_sim <- out$n_failed_sim[[1]]
    if (return_sim) {
      object$sim.delta <- out$sim.Delta_fast[[1]]
    }
  }
  class(object) <- "phylo_delta"
  attr(object, "match_report") <- attr(out, "match_report")
  if (is.list(analysis_metadata)) {
    attr(object, "analysis_metadata") <- analysis_metadata
  }
  object <- .decorate_fastphylosig_result(
    object, method = "Delta", vector_input = TRUE
  )
  timing <- .runtime_close(.runtime, success = TRUE)
  .attach_fast_signal_workflow(
    .runtime_attach(object, timing), "Delta", "fast_delta"
  )
}

# Number of retained samples from one MCMC chain under the sampler's inclusive
# burn-in rule.  Controls are validated by fast_delta() before this helper is
# called, but keeping the checks local makes it safe for internal callers too.
.delta_saved_iterations <- function(sim, thin, burn) {
  if (!is.numeric(sim) || length(sim) != 1L || !is.finite(sim) ||
      !is.numeric(thin) || length(thin) != 1L || !is.finite(thin) ||
      !is.numeric(burn) || length(burn) != 1L || !is.finite(burn) ||
      sim < 1 || thin < 1 || burn < 1 || burn > sim) {
    return(NA_real_)
  }
  as.numeric(1 + floor((sim - burn) / thin))
}

.delta_seed_metadata <- function() {
  present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  rng_kind <- tryCatch(RNGkind(), error = function(e) character())
  seed_length <- 0L
  seed_type <- NA_character_
  checksum <- NA_real_
  if (isTRUE(present)) {
    seed <- tryCatch(
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
      error = function(e) NULL
    )
    if (!is.null(seed)) {
      seed_length <- length(seed)
      seed_type <- typeof(seed)
      vals <- suppressWarnings(as.numeric(seed))
      vals <- vals[is.finite(vals)]
      if (length(vals)) {
        # A compact numeric fingerprint is enough to identify accidental RNG
        # state changes without exposing any raw seed entries.
        idx <- seq_along(vals)
        checksum <- sum((vals %% 2147483627) *
                          (idx %% 2147483627)) %% 2147483627
      }
    }
  }
  list(
    present = isTRUE(present),
    captured = FALSE,
    raw_exposed = FALSE,
    length = as.integer(seed_length),
    type = seed_type,
    rng_kind = rng_kind,
    rng_kind_code = paste(rng_kind, collapse = "|"),
    checksum = checksum,
    fingerprint = checksum,
    note = "Only seed presence is recorded; raw .Random.seed is not stored."
  )
}

.delta_p_mcse <- function(p, n_successful) {
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) ||
      !is.numeric(n_successful) || length(n_successful) != 1L ||
      !is.finite(n_successful) || n_successful < 1) {
    return(NA_real_)
  }
  p <- min(1, max(0, p))
  sqrt(p * (1 - p) / n_successful)
}

.delta_mcmc_warning_summary <- function(
    trait, delta, ess_alpha, ess_beta, rhat_alpha, rhat_beta, mcse_delta,
    min_ess = 20, max_rhat = 1.1) {
  fields <- list(
    trait = trait, delta = delta, ess_alpha = ess_alpha, ess_beta = ess_beta,
    rhat_alpha = rhat_alpha, rhat_beta = rhat_beta, mcse_delta = mcse_delta
  )
  n <- length(delta)
  if (any(lengths(fields) != n)) {
    stop("Delta diagnostic vectors must have a common length.",
         call. = FALSE)
  }
  if (!is.numeric(min_ess) || length(min_ess) != 1L ||
      !is.finite(min_ess) || min_ess <= 0 ||
      !is.numeric(max_rhat) || length(max_rhat) != 1L ||
      !is.finite(max_rhat) || max_rhat <= 1) {
    stop("Delta diagnostic thresholds must be finite and positive.",
         call. = FALSE)
  }

  finite_delta <- is.finite(delta)
  reason <- vapply(seq_len(n), function(i) {
    if (!finite_delta[i]) return(NA_character_)
    issues <- character()
    if (!is.finite(ess_alpha[i])) {
      issues <- c(issues, "ESS_alpha is non-finite")
    } else if (ess_alpha[i] < min_ess) {
      issues <- c(issues, sprintf("ESS_alpha < %g", min_ess))
    }
    if (!is.finite(ess_beta[i])) {
      issues <- c(issues, "ESS_beta is non-finite")
    } else if (ess_beta[i] < min_ess) {
      issues <- c(issues, sprintf("ESS_beta < %g", min_ess))
    }
    if (!is.finite(rhat_alpha[i])) {
      issues <- c(issues, "split_Rhat_alpha is non-finite")
    } else if (rhat_alpha[i] > max_rhat) {
      issues <- c(issues, sprintf("split_Rhat_alpha > %.1f", max_rhat))
    }
    if (!is.finite(rhat_beta[i])) {
      issues <- c(issues, "split_Rhat_beta is non-finite")
    } else if (rhat_beta[i] > max_rhat) {
      issues <- c(issues, sprintf("split_Rhat_beta > %.1f", max_rhat))
    }
    if (!is.finite(mcse_delta[i])) {
      issues <- c(issues, "MCSE_Delta is non-finite")
    }
    if (length(issues)) paste(issues, collapse = "; ") else NA_character_
  }, character(1))

  by_trait <- data.frame(
    trait = as.character(trait),
    finite_delta = finite_delta,
    ESS_alpha = as.numeric(ess_alpha),
    ESS_beta = as.numeric(ess_beta),
    split_Rhat_alpha = as.numeric(rhat_alpha),
    split_Rhat_beta = as.numeric(rhat_beta),
    MCSE_Delta = as.numeric(mcse_delta),
    warning = !is.na(reason),
    reason = reason,
    stringsAsFactors = FALSE
  )
  list(
    ess_minimum = as.numeric(min_ess),
    split_rhat_maximum = as.numeric(max_rhat),
    warning_emitted = any(by_trait$warning),
    by_trait = by_trait,
    affected_traits = by_trait[by_trait$warning, , drop = FALSE]
  )
}

.delta_chain_matrix <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x) && !is.data.frame(x)) {
    if (!length(x) || any(!vapply(x, is.numeric, logical(1)))) return(NULL)
    lens <- lengths(x)
    if (!length(lens) || any(lens != lens[[1L]]) || lens[[1L]] < 1L) {
      return(NULL)
    }
    return(do.call(cbind, lapply(x, as.numeric)))
  }
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) && !is.array(x)) return(NULL)
  x <- as.matrix(x)
  if (!is.numeric(x) || !nrow(x) || !ncol(x)) return(NULL)
  # The package kernel returns rows = saved iterations, columns = chains.
  # Accept the common transposed form for direct diagnostic callers when it
  # is unambiguous (two/few rows and many columns).
  if (nrow(x) <= 4L && ncol(x) > nrow(x)) x <- t(x)
  storage.mode(x) <- "double"
  x
}

.delta_chain_ess <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n <= 1L) return(as.numeric(n))
  center <- x - mean(x)
  denominator <- sum(center * center)
  if (!is.finite(denominator) || denominator <= 0) return(as.numeric(n))
  max_lag <- min(n - 1L, 1000L)
  tau <- 1
  if (max_lag >= 1L) {
    for (lag in seq.int(1L, max_lag, by = 2L)) {
      rho1 <- sum(center[seq_len(n - lag)] *
                    center[(lag + 1L):n]) / denominator
      pair <- rho1
      if (lag + 1L <= max_lag) {
        rho2 <- sum(center[seq_len(n - lag - 1L)] *
                      center[(lag + 2L):n]) / denominator
        pair <- pair + rho2
      }
      if (!is.finite(pair) || pair <= 0) break
      tau <- tau + 2 * pair
    }
  }
  if (!is.finite(tau) || tau < 1) tau <- 1
  max(1, min(n, n / tau))
}

.delta_split_rhat <- function(x) {
  x <- as.matrix(x)
  if (ncol(x) < 2L || nrow(x) < 4L) return(NA_real_)
  half <- floor(nrow(x) / 2L)
  if (half < 2L) return(NA_real_)
  blocks <- vector("list", 2L * ncol(x))
  k <- 0L
  for (j in seq_len(ncol(x))) {
    k <- k + 1L
    blocks[[k]] <- x[seq_len(half), j]
    k <- k + 1L
    blocks[[k]] <- x[seq.int(nrow(x) - half + 1L, nrow(x)), j]
  }
  means <- vapply(blocks, mean, numeric(1))
  W <- mean(vapply(blocks, function(z) {
    if (length(z) < 2L) return(0)
    sum((z - mean(z))^2) / (length(z) - 1L)
  }, numeric(1)))
  B <- half * stats::var(means)
  if (!is.finite(W) || !is.finite(B) || W < 0 || B < 0) return(NA_real_)
  if (W == 0) return(if (B == 0) 1 else NA_real_)
  sqrt(max(0, ((half - 1) / half * W + B / half) / W))
}

.delta_chain_mean_covariance <- function(alpha, beta) {
  n <- length(alpha)
  vals <- cbind(alpha, beta)
  if (n < 2L) return(matrix(0, nrow = 2L, ncol = 2L))
  batch <- max(1L, floor(sqrt(n)))
  n_batch <- floor(n / batch)
  if (n_batch < 2L) {
    return(stats::cov(vals) / n)
  }
  means <- do.call(rbind, lapply(seq_len(n_batch), function(i) {
    idx <- ((i - 1L) * batch + 1L):(i * batch)
    colMeans(vals[idx, , drop = FALSE])
  }))
  stats::cov(means) * batch / (n_batch * batch)
}

.delta_diagnostics_fallback <- function(alpha_chain, beta_chain) {
  alpha_chain <- .delta_chain_matrix(alpha_chain)
  beta_chain <- .delta_chain_matrix(beta_chain)
  if (is.null(alpha_chain) || is.null(beta_chain) ||
      !identical(dim(alpha_chain), dim(beta_chain)) ||
      any(!is.finite(alpha_chain)) || any(!is.finite(beta_chain))) {
    return(list(
      available = FALSE,
      note = paste(
        "MCMC traces are unavailable; delta_mcmc_cpp must return",
        "alpha_chain and beta_chain for ESS, R-hat, and MCSE diagnostics."
      )
    ))
  }

  n <- nrow(alpha_chain)
  chains <- ncol(alpha_chain)
  aa <- as.numeric(alpha_chain)
  bb <- as.numeric(beta_chain)
  alpha_mean <- mean(aa)
  beta_mean <- mean(bb)
  mean_cov <- matrix(0, nrow = 2L, ncol = 2L)
  ess_alpha <- ess_beta <- 0
  chain_alpha_mean <- chain_beta_mean <- numeric(chains)
  chain_alpha_sd <- chain_beta_sd <- numeric(chains)
  for (j in seq_len(chains)) {
    a <- alpha_chain[, j]
    b <- beta_chain[, j]
    ess_alpha <- ess_alpha + .delta_chain_ess(a)
    ess_beta <- ess_beta + .delta_chain_ess(b)
    chain_alpha_mean[j] <- mean(a)
    chain_beta_mean[j] <- mean(b)
    chain_alpha_sd[j] <- if (length(a) > 1L) stats::sd(a) else 0
    chain_beta_sd[j] <- if (length(b) > 1L) stats::sd(b) else 0
    mean_cov <- mean_cov + .delta_chain_mean_covariance(a, b) / chains^2
  }
  g <- c(-beta_mean / alpha_mean^2, 1 / alpha_mean)
  mcse <- if (is.finite(alpha_mean) && alpha_mean != 0 &&
              all(is.finite(mean_cov))) {
    variance_terms <- c(
      g[1L]^2 * mean_cov[1L, 1L],
      2 * g[1L] * g[2L] * mean_cov[1L, 2L],
      g[2L]^2 * mean_cov[2L, 2L]
    )
    var_delta <- sum(variance_terms)
    roundoff_tol <- 64 * .Machine$double.eps *
      max(1, sum(abs(variance_terms)))
    if (is.finite(var_delta) && var_delta >= -roundoff_tol) {
      sqrt(max(0, var_delta))
    } else {
      NA_real_
    }
  } else {
    NA_real_
  }
  list(
    available = TRUE,
    note = "Diagnostics computed from saved alpha/beta MCMC traces.",
    alpha_mean = alpha_mean,
    beta_mean = beta_mean,
    alpha_sd = if (length(aa) > 1L) stats::sd(aa) else 0,
    beta_sd = if (length(bb) > 1L) stats::sd(bb) else 0,
    chain_alpha_mean = chain_alpha_mean,
    chain_beta_mean = chain_beta_mean,
    chain_alpha_sd = chain_alpha_sd,
    chain_beta_sd = chain_beta_sd,
    alpha_chain_mean = chain_alpha_mean,
    beta_chain_mean = chain_beta_mean,
    alpha_chain_sd = chain_alpha_sd,
    beta_chain_sd = chain_beta_sd,
    ESS_alpha = max(1, min(length(aa), ess_alpha)),
    ESS_beta = max(1, min(length(bb), ess_beta)),
    split_Rhat_alpha = .delta_split_rhat(alpha_chain),
    split_Rhat_beta = .delta_split_rhat(beta_chain),
    Rhat_alpha = .delta_split_rhat(alpha_chain),
    Rhat_beta = .delta_split_rhat(beta_chain),
    alpha_beta_cov = mean_cov[1L, 2L],
    mean_covariance = mean_cov,
    MCSE_Delta = mcse,
    n_saved = length(aa),
    n_chains = chains,
    n_iter = n
  )
}

.delta_mcmc_diagnostics <- function(calc) {
  if (!is.list(calc)) {
    return(.delta_diagnostics_fallback(NULL, NULL))
  }
  alpha_chain <- if (!is.null(calc$alpha_chain)) {
    calc$alpha_chain
  } else {
    calc$alpha_chains
  }
  beta_chain <- if (!is.null(calc$beta_chain)) {
    calc$beta_chain
  } else {
    calc$beta_chains
  }

  # A generated Rcpp wrapper is present after compileAttributes().  During
  # source-only testing it may be absent, in which case the same formulas are
  # evaluated by the R fallback below.
  cpp_fun <- get0("delta_diagnostics_cpp", mode = "function", inherits = TRUE)
  if (is.function(cpp_fun) && !is.null(alpha_chain) && !is.null(beta_chain)) {
    cpp <- tryCatch(
      cpp_fun(.delta_chain_matrix(alpha_chain), .delta_chain_matrix(beta_chain)),
      error = function(e) NULL
    )
    if (is.list(cpp)) {
      cpp$available <- TRUE
      cpp$note <- "Diagnostics computed from saved alpha/beta MCMC traces."
      return(cpp)
    }
  }

  diag <- .delta_diagnostics_fallback(alpha_chain, beta_chain)
  if (isFALSE(diag$available)) {
    diag$alpha_mean <- if (!is.null(calc$alpha_mean)) {
      calc$alpha_mean
    } else {
      NA_real_
    }
    diag$beta_mean <- if (!is.null(calc$beta_mean)) {
      calc$beta_mean
    } else {
      NA_real_
    }
    diag$alpha_sd <- NA_real_
    diag$beta_sd <- NA_real_
    diag$ESS_alpha <- NA_real_
    diag$ESS_beta <- NA_real_
    diag$split_Rhat_alpha <- NA_real_
    diag$split_Rhat_beta <- NA_real_
    diag$alpha_beta_cov <- NA_real_
    diag$MCSE_Delta <- NA_real_
  }
  diag
}

.delta_permutation_worker <- function(i, y, perms, tree, lambda0, proposal_sd,
                                       mcmc_sim, thin, burn, entropy_code,
                                       model, ace_engine, ace_workspace = NULL,
                                       reuse_buffers = TRUE) {
  yp <- y[perms[i, ]]
  names(yp) <- names(y)
  sim_calc <- tryCatch(
    .delta_one(
      tree = tree, y = yp, lambda0 = lambda0, proposal_sd = proposal_sd,
      mcmc_sim = mcmc_sim, thin = thin, burn = burn,
      entropy_code = entropy_code, model = model, ace_engine = ace_engine,
       keep_chains = FALSE, ace_workspace = ace_workspace,
       reuse_buffers = reuse_buffers
    ),
    error = function(e) e
  )
  if (inherits(sim_calc, "error")) {
    return(NA_real_)
  }
  sim_calc$delta
}

.delta_one <- function(tree, y, lambda0, proposal_sd, mcmc_sim, thin, burn,
                       entropy_code, model, ace_engine, keep_chains = TRUE,
                       ace_workspace = NULL, reuse_buffers = TRUE) {
  y <- y[tree$tip.label]
  ar <- if (identical(ace_engine, "fast")) {
    if (is.null(ace_workspace)) {
      fast_ace(
        y, tree, type = "discrete", method = "ML", model = model,
        CI = TRUE, marginal = FALSE, progress = FALSE
      )$lik.anc
    } else {
      .fast_ace_lik_anc_workspace(
        y, ace_workspace, reuse_buffers = reuse_buffers
      )
    }
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
    sim = mcmc_sim, thin = thin, burn = burn, entropy_type = entropy_code,
    return_chains = isTRUE(keep_chains)
  )
}
