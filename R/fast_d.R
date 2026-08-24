fast_d <- function(tree, x = NULL, test = TRUE, nsim = 1000,
                   se = NULL, return_sim = test, verbose = TRUE,
                   rnd.bias = NULL, permutations = NULL,
                   random_states = NULL, brownian_states = NULL,
                   ncores = 1, X = NULL, permut = NULL, prepared = NULL,
                   chunk_size = 128L, keep_null = return_sim,
                   progress = interactive()) {
  if (!missing(x) && !missing(X)) {
    stop("Supply only one of x or X, not both.", call. = FALSE)
  }
  runtime <- .runtime_begin("D", progress = progress, verbose = verbose)
  on.exit(.runtime_on_exit(runtime), add = TRUE)
  # Keep the legacy verbose argument as a status-output master switch.  This
  # also suppresses the historical matching notices when progress is FALSE.
  verbose <- isTRUE(verbose) && isTRUE(progress)
  .runtime_stage(runtime, "Checking tree...")
  .runtime_stage(runtime, "Preparing data...")
  if (!is.null(prepared)) tree <- prepared
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
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      !is.finite(chunk_size) || chunk_size != floor(chunk_size) ||
      chunk_size < 32L || chunk_size > 512L) {
    stop("chunk_size must be an integer between 32 and 512.",
         call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)
  if (!is.logical(keep_null) || length(keep_null) != 1L ||
      is.na(keep_null)) {
    stop("keep_null must be TRUE or FALSE.", call. = FALSE)
  }
  store_null <- isTRUE(return_sim) && isTRUE(keep_null)
  if (!is.null(permut)) nsim <- permut
  if (!is.null(se)) {
    stop("fast_d does not use se; provide binary states in x.",
         call. = FALSE)
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim < 1 || nsim != floor(nsim)) {
    stop("nsim/permut must be a positive integer.", call. = FALSE)
  }
  nsim <- as.integer(nsim)

  vector_input <- is.null(dim(x))
  # The shared preparation layer performs the method-specific D preflight,
  # safe representation canonicalisation, matching, and packed NA grouping.
  # Keep the binary-table conversion here for the historical input contract;
  # preparation itself is read-only with respect to the caller's data.
  X <- .as_named_trait_table(
    x, if (inherits(tree, "fastphylosig_tree")) tree$tree else tree,
    verbose = FALSE, input_name = "x"
  )
  analysis <- .prepare_analysis(
    tree = tree, data = X, signal = "D", data_kind = "binary",
    verbose = verbose
  )
  ctx <- analysis$ctx
  tree <- ctx$tree
  X0 <- analysis$matched_data
  tree_tips <- tree$tip.label
  matched <- rownames(X0)
  removed_tree <- setdiff(tree_tips, matched)
  removed_data <- setdiff(rownames(X), tree_tips)
  if (length(matched) < 2L) {
    stop("Fewer than two species are shared by tree and x.", call. = FALSE)
  }

  base_keep <- analysis$matching$base_keep
  if (is.null(base_keep)) base_keep <- match(matched, tree_tips)
  base_group <- analysis$base_group
  if (is.null(base_group)) {
    stop("Fewer than two species are shared by tree and x.", call. = FALSE)
  }
  base_tree <- base_group$tree
  base_idx <- match(base_tree$tip.label, tree_tips)
  # `matched_data` is already in tree-tip order.  Keep this explicit subset to
  # preserve the exact row ordering expected by the legacy output contract.
  X0 <- X0[base_tree$tip.label, , drop = FALSE]

  # Keep result provenance compact: prepared contexts contain environments and
  # compiled structural entries, so never attach the full analysis object.
  analysis_metadata <- .compact_analysis_metadata(analysis)

  .runtime_stage(runtime, "Calculating observed D...")
  if (isTRUE(test)) {
    .runtime_stage(runtime, "Simulating random and Brownian nulls...")
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
  MCSE_P_random <- MCSE_P_Brownian <- rep(NA_real_, ncol(X0))
  # nsim_requested is the user-facing request.  The two successful counts
  # are tracked separately because a future null generator can fail for one
  # null while the other remains usable.  A failed count is NA until a trait
  # actually reaches the simulation stage (invalid traits do not silently
  # pretend that nsim draws were attempted).
  nsim_requested <- rep.int(nsim, ncol(X0))
  nsim_successful_random <- nsim_successful_brownian <-
    rep.int(0L, ncol(X0))
  nsim_failed_random <- nsim_failed_brownian <-
    rep(NA_integer_, ncol(X0))
  n_species <- n_removed_na <- rep(NA_integer_, ncol(X0))
  notes <- rep(NA_character_, ncol(X0))
  nonfinite_trait <- vapply(X0, function(z) {
    is.numeric(z) && any(is.infinite(z) | is.nan(z))
  }, logical(1))
  if (any(nonfinite_trait)) {
    notes[nonfinite_trait] <-
      "trait contains non-finite values other than NA"
  }
  state_tables <- vector("list", ncol(X0))
  random_list <- brownian_list <- vector("list", ncol(X0))

  # Reuse the packed NA masks and structural groups produced by preparation;
  # no second mask grouping or tree matching pass is performed here.
  mask_groups <- analysis$na_patterns
  for (group_index in seq_len(mask_groups$n_group)) {
    idx <- mask_groups$columns[[group_index]]
    keep <- mask_groups$keep[[group_index]]
    n_species[idx] <- length(keep)
    n_removed_na[idx] <- nrow(X0) - length(keep)
    if (length(keep) < 2L) {
      notes[idx] <- "fewer than 2 non-NA matched species"
      next
    }

    prepared_group <- analysis$groups[[group_index]]
    group <- if (is.list(prepared_group)) prepared_group$group else NULL
    if (is.null(group)) {
      notes[idx] <- "retained species subset is not ready for D; run check_tree()"
      next
    }
    group_tree <- group$tree
    # Dropping unmatched/NA tips can expose a unary internal node.  The shared
    # preparation check records that method-specific issue once per unique
    # retained subset; do not silently send an unsafe group to C++.
    group_check <- if (is.list(prepared_group)) prepared_group$check else NULL
    if (!is.null(group_check) &&
        !isTRUE(group_check$ready_by_signal[["D"]])) {
      issue_text <- if (exists(".format_actionable_condition",
                              mode = "function", inherits = TRUE)) {
        .format_actionable_condition(
          group_check, signal = "D", prefix = "D retained subset is not ready"
        )
      } else {
        "D retained subset is not ready; run check_tree()"
      }
      notes[idx] <- issue_text
      next
    }
    phy <- group$d_tree
    edge <- group$edge
    edge_length <- group$edge_length
    n_tip <- group$n_tip

    for (j in idx) {
      if (isTRUE(nonfinite_trait[[j]])) next
      trait_values <- X0[keep, j]
      if (is.numeric(trait_values) &&
          any(!is.finite(trait_values[!is.na(trait_values)]))) {
        notes[j] <- "trait contains non-finite values other than NA"
        next
      }
      d <- .binary_state(stats::setNames(X0[keep, j], rownames(X0)[keep]))
      if (!is.null(d$error)) {
        notes[j] <- d$error
        next
      }

      ds <- d$values[phy$tip.label]
      bias <- if (is.null(bias0)) NULL else as.numeric(bias0[phy$tip.label])
      stream <- NULL
      stream_candidate <- is.null(bias) && is.null(permutations) &&
        is.null(random_states) && is.null(brownian_states) &&
        .phylo_d_binary_edge(edge, n_tip) &&
        exists("phylo_d_stream_cpp", mode = "function", inherits = TRUE)
      if (stream_candidate) {
        # The streaming kernel fuses random/Brownian generation with the
        # contrast traversal and only retains O(nsim) sums when requested.
        stream <- phylo_d_stream_cpp(
          observed = as.numeric(ds), edge = edge,
          edge_length = edge_length, n_tip = n_tip, nsim = nsim,
          prop_state1 = d$prop_state1, chunk_size = chunk_size,
          return_sim = store_null,
          n_threads = ncores
        )
        observed_j <- stream$observed
        # With keep_null = FALSE the streaming C++ kernel only returns online
        # summaries.  Its input and branch-length validation guarantees that
        # each generated draw is finite; expose that contract explicitly as
        # nsim successful draws without allocating the null vectors.
        if (store_null) {
          random_summary <- .d_null_summary(
            stream$random, observed_j, direction = "less"
          )
          brownian_summary <- .d_null_summary(
            stream$brownian, observed_j, direction = "greater"
          )
        } else {
          random_summary <- .d_stream_summary(
            stream$mean_random, stream$p_random, nsim,
            label = "random"
          )
          brownian_summary <- .d_stream_summary(
            stream$mean_brownian, stream$p_brownian, nsim,
            label = "Brownian"
          )
        }
        mean_random[j] <- random_summary$mean
        mean_brownian[j] <- brownian_summary$mean
        nsim_successful_random[j] <- random_summary$n_successful
        nsim_successful_brownian[j] <- brownian_summary$n_successful
        nsim_failed_random[j] <- random_summary$n_failed
        nsim_failed_brownian[j] <- brownian_summary$n_failed
        Pval1_fast[j] <- random_summary$p
        Pval0_fast[j] <- brownian_summary$p
        MCSE_P_random[j] <- random_summary$mcse
        MCSE_P_Brownian[j] <- brownian_summary$mcse
        if (length(random_summary$note)) {
          notes[j] <- .append_d_note(notes[j], random_summary$note)
        }
        if (length(brownian_summary$note)) {
          notes[j] <- .append_d_note(notes[j], brownian_summary$note)
        }
        denom <- mean_random[j] - mean_brownian[j]
        D_fast[j] <- if (is.finite(denom) && denom != 0) {
          (observed_j - mean_brownian[j]) / denom
        } else {
          denominator_note <- if (!is.finite(denom)) {
            "random or Brownian null mean is unavailable"
          } else {
            "random and Brownian null means are identical"
          }
          notes[j] <- .append_d_note(
            notes[j], denominator_note
          )
          NA_real_
        }
        if (store_null) {
          random_list[[j]] <- as.numeric(stream$random)
          brownian_list[[j]] <- as.numeric(stream$brownian)
        }
      } else {
        rand <- .phylo_d_random_states(
          ds = ds, nsim = nsim, rnd.bias = bias, permutations = permutations,
          random_states = random_states
        )
        brown <- .phylo_d_brownian_states(
          tree = group_tree, nsim = nsim, prop_state1 = d$prop_state1,
          brownian_states = brownian_states
        )

        # One C++ batch traversal handles observed, random, and Brownian states.
        all_sums <- phylo_d_sums_cpp(
          cbind(Obs = ds, rand, brown), edge = edge,
          edge_length = edge_length, n_tip = n_tip, n_threads = ncores
        )
        observed_j <- all_sums[[1L]]
        rans <- all_sums[seq.int(2L, nsim + 1L)]
        phys <- all_sums[seq.int(nsim + 2L, 2L * nsim + 1L)]

        random_summary <- .d_null_summary(
          rans, observed_j, direction = "less"
        )
        brownian_summary <- .d_null_summary(
          phys, observed_j, direction = "greater"
        )
        mean_random[j] <- random_summary$mean
        mean_brownian[j] <- brownian_summary$mean
        nsim_successful_random[j] <- random_summary$n_successful
        nsim_successful_brownian[j] <- brownian_summary$n_successful
        nsim_failed_random[j] <- random_summary$n_failed
        nsim_failed_brownian[j] <- brownian_summary$n_failed
        Pval1_fast[j] <- random_summary$p
        Pval0_fast[j] <- brownian_summary$p
        MCSE_P_random[j] <- random_summary$mcse
        MCSE_P_Brownian[j] <- brownian_summary$mcse
        if (length(random_summary$note)) {
          notes[j] <- .append_d_note(notes[j], random_summary$note)
        }
        if (length(brownian_summary$note)) {
          notes[j] <- .append_d_note(notes[j], brownian_summary$note)
        }
        denom <- mean_random[j] - mean_brownian[j]
        D_fast[j] <- if (is.finite(denom) && denom != 0) {
          (observed_j - mean_brownian[j]) / denom
        } else {
          denominator_note <- if (!is.finite(denom)) {
            "random or Brownian null mean is unavailable"
          } else {
            "random and Brownian null means are identical"
          }
          notes[j] <- .append_d_note(
            notes[j], denominator_note
          )
          NA_real_
        }
        if (store_null) {
          random_list[[j]] <- rans
          brownian_list[[j]] <- phys
        }
      }

      observed[j] <- observed_j
      state_tables[[j]] <- d$states_table
    }
  }

  status <- rep("failed", ncol(X0))
  status[is.finite(D_fast)] <- "ok"
  status[grepl("fewer than 2", notes, fixed = TRUE) %in% TRUE] <-
    "insufficient_data"
  status[grepl("non-finite|two states|single state|binary", notes,
               ignore.case = TRUE) %in% TRUE] <- "invalid_trait"
  partial <- is.finite(D_fast) & isTRUE(test) &
    (nsim_successful_random < nsim_requested |
       nsim_successful_brownian < nsim_requested)
  status[partial %in% TRUE] <- "partial"

  out <- data.frame(
    trait = colnames(X0),
    D_fast = D_fast,
    Pval1_fast = if (test) Pval1_fast else NA_real_,
    Pval0_fast = if (test) Pval0_fast else NA_real_,
    # Clear null-specific aliases are additive; the legacy Pval1/Pval0 names
    # above remain unchanged for caper-compatible callers.
    P_random = if (test) Pval1_fast else NA_real_,
    P_Brownian = if (test) Pval0_fast else NA_real_,
    MCSE_P_random = if (test) MCSE_P_random else NA_real_,
    MCSE_P_Brownian = if (test) MCSE_P_Brownian else NA_real_,
    nsim_requested = nsim_requested,
    nsim_successful_random = nsim_successful_random,
    nsim_successful_brownian = nsim_successful_brownian,
    nsim_failed_random = nsim_failed_random,
    nsim_failed_brownian = nsim_failed_brownian,
    observed = observed,
    mean_random = mean_random,
    mean_brownian = mean_brownian,
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
  if (store_null) {
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
    data_rows_removed = removed_data,
    analysis_metadata = analysis_metadata
  )
  attr(out, "analysis_metadata") <- analysis_metadata

  if (!vector_input) {
    object <- .decorate_fastphylosig_result(
      out, method = "D", vector_input = FALSE
    )
    attr(object, "analysis_metadata") <- analysis_metadata
    timing <- .runtime_close(runtime, success = TRUE)
    return(.attach_fast_signal_workflow(
      .runtime_attach(object, timing), "D", "fast_d"
    ))
  }

  object <- list(
    DEstimate = out$D_fast[[1]],
    Pval1 = out$Pval1_fast[[1]],
    Pval0 = out$Pval0_fast[[1]],
    P_random = out$P_random[[1]],
    P_Brownian = out$P_Brownian[[1]],
    MCSE_P_random = out$MCSE_P_random[[1]],
    MCSE_P_Brownian = out$MCSE_P_Brownian[[1]],
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
    nsim_requested = out$nsim_requested[[1]],
    nsim_successful_random = out$nsim_successful_random[[1]],
    nsim_successful_brownian = out$nsim_successful_brownian[[1]],
    nsim_failed_random = out$nsim_failed_random[[1]],
    nsim_failed_brownian = out$nsim_failed_brownian[[1]],
    status = out$status[[1]],
    note = out$note[[1]],
    message = out$message[[1]],
    rnd.bias = rnd.bias
  )
  if (store_null) {
    object$Permutations <- list(
      random = out$random_fast[[1]],
      brownian = out$brownian_fast[[1]]
    )
  }
  class(object) <- "phylo.d"
  attr(object, "match_report") <- attr(out, "match_report")
  attr(object, "analysis_metadata") <- analysis_metadata
  object <- .decorate_fastphylosig_result(
    object, method = "D", vector_input = TRUE
  )
  timing <- .runtime_close(runtime, success = TRUE)
  .attach_fast_signal_workflow(
    .runtime_attach(object, timing), "D", "fast_d"
  )
}

.binary_state <- function(x) {
  if (any(is.infinite(x), na.rm = TRUE)) {
    return(list(error = "contains infinite values"))
  }
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

.phylo_d_binary_edge <- function(edge, n_tip) {
  if (is.null(dim(edge)) || ncol(edge) != 2L || nrow(edge) < 2L) {
    return(FALSE)
  }
  parent_count <- table(edge[, 1L])
  internal <- as.integer(names(parent_count))
  # The streaming path is intentionally limited to ordinary rooted binary
  # trees.  Any polytomy (or malformed unary internal node) stays on the
  # legacy oracle path, which retains the existing caper-compatible split.
  if (length(internal) == 0L || any(internal <= n_tip)) return(FALSE)
  all(as.integer(parent_count) == 2L)
}

.phylo_d_random_states <- function(ds, nsim, rnd.bias, permutations,
                                   random_states) {
  n <- length(ds)
  if (!is.null(random_states)) {
    random_states <- as.matrix(random_states)
    storage.mode(random_states) <- "double"
    if (any(!is.finite(random_states))) {
      stop("random_states must contain only finite values.", call. = FALSE)
    }
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
    if (any(!is.finite(brownian_states))) {
      stop("brownian_states must contain only finite values.", call. = FALSE)
    }
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

# Summarise one null distribution without changing the tail definition used
# by caper::phylo.d(): random uses values < observed, while Brownian uses
# values > observed.  Non-finite simulation statistics are not silently
# counted as non-extremes; they are excluded from the denominator and exposed
# through the successful/failed counts and note.  The resulting P estimate is
# therefore explicitly conditional on successful finite draws.
.d_null_summary <- function(values, observed, direction = c("less", "greater")) {
  direction <- match.arg(direction)
  values <- as.numeric(values)
  finite <- is.finite(values)
  n_successful <- sum(finite)
  n_failed <- length(values) - n_successful
  if (n_successful < 1L || !is.finite(observed)) {
    note <- if (n_failed > 0L) {
      sprintf("%s null has %d non-finite simulation(s); no successful draws",
              if (direction == "less") "random" else "Brownian", n_failed)
    } else {
      sprintf("%s null has no successful finite simulations",
              if (direction == "less") "random" else "Brownian")
    }
    return(list(mean = NA_real_, p = NA_real_, mcse = NA_real_,
                n_successful = n_successful, n_failed = n_failed,
                note = note))
  }

  finite_values <- values[finite]
  mean_value <- mean(finite_values)
  extreme <- if (direction == "less") {
    sum(finite_values < observed)
  } else {
    sum(finite_values > observed)
  }
  p_value <- extreme / n_successful
  # The plug-in Bernoulli Monte Carlo standard error is defined for the
  # boundary cases too: P=0 and P=1 both have MCSE exactly zero.
  mcse <- sqrt(p_value * (1 - p_value) / n_successful)
  note <- if (n_failed > 0L) {
    sprintf("%s null excluded %d non-finite simulation(s) from P and MCSE",
            if (direction == "less") "random" else "Brownian", n_failed)
  } else {
    character()
  }
  list(mean = mean_value, p = p_value, mcse = mcse,
       n_successful = n_successful, n_failed = n_failed, note = note)
}

# Summary path for the streaming kernel when null vectors are deliberately not
# retained.  The kernel reports online mean/tail summaries and guarantees
# finite generated draws after validating the tree and observed states.  If a
# future kernel reports a non-finite summary, mark all draws unsuccessful
# rather than inventing a P value.
.d_stream_summary <- function(mean_value, p_value, nsim, label) {
  ok <- is.finite(mean_value) && is.finite(p_value) &&
    p_value >= 0 && p_value <= 1
  if (!ok) {
    return(list(mean = NA_real_, p = NA_real_, mcse = NA_real_,
                n_successful = 0L, n_failed = as.integer(nsim),
                note = sprintf(
                  "%s streaming null returned a non-finite or invalid summary",
                  label)))
  }
  mcse <- sqrt(p_value * (1 - p_value) / nsim)
  list(mean = as.numeric(mean_value), p = as.numeric(p_value),
       mcse = mcse, n_successful = as.integer(nsim), n_failed = 0L,
       note = character())
}

.append_d_note <- function(existing, new_note) {
  if (length(new_note) == 0L || is.na(new_note) || !nzchar(new_note)) {
    return(existing)
  }
  if (is.na(existing) || !nzchar(existing)) new_note else {
    paste(existing, new_note, sep = "; ")
  }
}
