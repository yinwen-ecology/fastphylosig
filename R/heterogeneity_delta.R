# Thin orchestration of the existing ACE workspace and Delta MCMC engines.
# No alternative Delta estimator or entropy formula is implemented here.
.heterogeneity_delta_control <- function(control) {
  defaults <- list(mcmc_sim = 10000L, thin = 10L, burn = 100L,
    lambda0 = .1, proposal_sd = .5, entropy = "LSE", model = "ARD")
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyNA(names(control)) || any(!nzchar(names(control))) || anyDuplicated(names(control)))))
    stop("delta_control must be a uniquely named list.", call. = FALSE)
  unknown <- setdiff(names(control), c(names(defaults), "reference_fit", "reference_unavailable",
    "reference_test", "state_space", "nuisance_control", "statistic"))
  if (length(unknown)) stop("Unknown delta_control: ", paste(unknown, collapse = ", "), call. = FALSE)
  z <- utils::modifyList(defaults, control, keep.null = TRUE)
  for (nm in c("mcmc_sim", "thin", "burn")) z[[nm]] <- .heterogeneity_integer(z[[nm]], nm, 1L)
  if (z$burn > z$mcmc_sim) stop("burn must not exceed mcmc_sim.", call. = FALSE)
  for (nm in c("lambda0", "proposal_sd"))
    if (!is.numeric(z[[nm]]) || length(z[[nm]]) != 1L || !is.finite(z[[nm]]) || z[[nm]] <= 0)
      stop(nm, " must be a positive finite number.", call. = FALSE)
  z$entropy <- match.arg(z$entropy, c("LSE", "SE", "GINI"))
  z$model <- match.arg(z$model, c("ARD", "ER"))
  if ("reference_fit" %in% names(z)) {
    z$reference_fit <- match.arg(z$reference_fit, c("nlminb", "profile", "multistart"))
    if (z$reference_fit == "profile" && z$model != "ER")
      stop("delta_control$reference_fit='profile' currently requires model='ER'.", call. = FALSE)
    if (z$reference_fit == "multistart" && z$model != "ARD")
      stop("delta_control$reference_fit='multistart' requires model='ARD'.", call. = FALSE)
  }
  if ("reference_unavailable" %in% names(z))
    z$reference_unavailable <- match.arg(z$reference_unavailable, c("strict", "conservative"))
  if ("reference_test" %in% names(z))
    z$reference_test <- match.arg(z$reference_test, c("plugin", "nuisance"))
  if (identical(z$reference_test, "nuisance")) {
    if("statistic"%in%names(z))z$statistic <- match.arg(z$statistic,c("local_global","between_clades"))
    if (is.null(z$state_space)) stop("Nuisance testing requires an explicit state_space.", call. = FALSE)
    if ("reference_fit" %in% names(z) || "reference_unavailable" %in% names(z))
      stop("reference_fit/reference_unavailable apply to plugin tests; nuisance mode retains single-state upper tails and actual failures.", call. = FALSE)
    z$nuisance_control <- .heterogeneity_nuisance_control(if (is.null(z$nuisance_control)) list() else z$nuisance_control)
    if(z$nuisance_control$proposal_family=="ARD" && z$model!="ARD")
      stop("proposal_family='ARD' requires the three-state ARD model.",call.=FALSE)
  } else if (any(c("state_space","nuisance_control","statistic")%in%names(z)))
    stop("state_space/nuisance_control/statistic require reference_test='nuisance'.", call. = FALSE)
  z
}

.heterogeneity_delta_batch <- function(trees, data, columns, control, n, rows = NULL) {
  out <- list(estimate = rep(NA_real_, n), status = rep("undefined", n),
    message = rep(NA_character_, n), details = vector("list", n))
  contexts <- lapply(trees, function(tree) .heterogeneity_delta_context(tree,
    tree$tip.label, tree$tip.label, control, canonicalized = TRUE))
  group <- position <- integer(n)
  for (i in seq_along(columns)) {
    cols <- columns[[i]] + 1L
    group[cols] <- i; position[cols] <- seq_along(cols)
  }
  # Reuse NA/model workspaces while retaining original trait-column RNG order.
  for (j in seq_len(n)) {
    i <- group[j]
    if (!i) next
    fit <- .heterogeneity_delta_evaluate(contexts[[i]], data[[i]][,position[j]])
    out$estimate[j] <- fit$estimate; out$status[j] <- fit$status
    out$message[j] <- fit$message; out$details[[j]] <- fit$details
  }
  if (!is.null(rows)) {
    # Only the tests need these contexts. Remap their local scan rows to the
    # original full-tree matrix after all observed traits have been evaluated.
    for (i in seq_along(contexts)) contexts[[i]]$rows <- rows[[i]]
    out$contexts <- lapply(group, function(i) if (i) contexts[[i]] else NULL)
  }
  out
}

.heterogeneity_categorical_matrix <- function(table, tips, method = "Delta") {
  X <- matrix(NA_real_, nrow(table), ncol(table), dimnames = dimnames(table))
  labels <- vector("list", ncol(table))
  matched <- rownames(table) %in% tips
  for (j in seq_len(ncol(table))) {
    z <- table[[j]]
    if (is.numeric(z) && any((is.infinite(z) | is.nan(z)) & matched))
      stop(method, " trait '", names(table)[j], "' contains non-finite values other than NA.", call. = FALSE)
    observed <- !is.na(z) & matched
    y <- droplevels(factor(z[observed]))
    labels[[j]] <- levels(y)
    if (method == "D" && nlevels(y) > 2L)
      stop("D trait '", names(table)[j], "' contains more than two states.", call. = FALSE)
    if (any(observed)) X[observed,j] <- as.integer(y) - if (method == "D") 1L else 0L
  }
  attr(X, "state_labels") <- labels
  X
}

# Fixed tree/NA contexts are prepared once per trait, not once per replicate.
# Local alphabets shrink in the same way as public fast_delta(). The global
# alphabet remains fixed for fitting and simulation of the common reference.
.heterogeneity_delta_context <- function(tree, keep, full_tips, control, canonicalized = FALSE) {
  # Match public Delta's canonicalize-then-prune order, including entropy rows
  # and endpoint-jitter draw order. Reordering ancestors after pruning would
  # change finite MCMC results even though the tree is isomorphic.
  if (!canonicalized) tree <- .safe_canonicalize_core(tree)
  if (length(keep) != length(tree$tip.label)) tree <- ape::keep.tip(tree, keep)
  tree_workspace <- .fast_ace_tree_workspace(tree)
  list(tree = tree, rows = match(tree$tip.label, full_tips),
    tree_workspace = tree_workspace, workspaces = new.env(parent = emptyenv()),
    control = control)
}

.heterogeneity_delta_evaluate <- function(context, x) {
  .heterogeneity_delta_values(context, x[context$rows])
}

.heterogeneity_delta_values <- function(context, values, diagnostics = TRUE) {
  y <- factor(values); names(y) <- context$tree$tip.label
  details <- list(n_saved = NA_real_, ESS_alpha = NA_real_, ESS_beta = NA_real_,
    ESS_Delta = NA_real_, split_Rhat_Delta = NA_real_,
    split_Rhat_alpha = NA_real_, split_Rhat_beta = NA_real_, MCSE_Delta = NA_real_,
    MCMC_warning = FALSE)
  details[.delta_rank_fields()] <- rep(list(NA_real_), length(.delta_rank_fields()))
  details <- c(details, .delta_calculation_flags())
  out <- list(estimate = NA_real_, warning = FALSE, status = "undefined",
    message = NA_character_, fit_warning = FALSE, single_state = FALSE, details = details)
  if (nlevels(y) < 2L) {
    out$single_state <- TRUE
    out$message <- "Delta is undefined for a single-state crown"; return(out)
  }
  key <- as.character(nlevels(y))
  if (is.null(context$workspaces[[key]])) context$workspaces[[key]] <-
    .fast_ace_workspace(context$tree_workspace, context$control$model, nlevels(y),
      native_pool=attr(context$control,"ace_pool"))
  notes <- character()
  calc <- tryCatch(withCallingHandlers(do.call(.delta_one,
    c(list(tree = context$tree, y = y,
      entropy_code = match(context$control$entropy, c("LSE", "SE", "GINI")),
      ace_engine = "fast", keep_chains = diagnostics,
      ace_workspace = context$workspaces[[key]]),
      context$control[setdiff(names(context$control), c("entropy", "reference_fit", "reference_unavailable",
        "reference_test", "state_space", "nuisance_control", "statistic"))])),
    warning = function(w) { notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e)
  if (inherits(calc, "error")) {
    out$message <- conditionMessage(calc); out$warning <- TRUE
    flags <- calc$calculation_diagnostics
    if (is.null(flags)) {
      flags <- .delta_calculation_flags()
      flags$calculation_failed <- TRUE; flags$failure_stage <- "calculation"
    }
    out$details[names(flags)] <- flags
    out$fit_warning <- flags$optimizer_warning; return(out)
  }
  flags <- calc$calculation_diagnostics
  out$details[names(flags)] <- flags
  out$fit_warning <- flags$optimizer_warning
  out$details$n_saved <- calc$n_saved
  if (diagnostics) {
    d <- .delta_mcmc_diagnostics(calc)
    native_flag <- get0("delta_warning_reason_cpp",mode="function",inherits=TRUE)
    # R's paste() respects OutDec/scipen whereas sprintf() does not. Retain the
    # existing formatter for custom options rather than changing warning text.
    native_format <- identical(getOption("OutDec"),".") && isTRUE(getOption("scipen")==0)
    reason <- if (is.function(native_flag) && native_format) native_flag(calc$delta,d) else
      .delta_mcmc_warning_summary("reference", calc$delta,
      d$ESS_alpha, d$ESS_beta, d$split_Rhat_alpha, d$split_Rhat_beta,
      d$MCSE_Delta, ess_delta = d$ESS_Delta, rhat_delta = d$split_Rhat_Delta,
      rank_diagnostics = d[.delta_rank_fields()])$by_trait$reason
    for (nm in setdiff(names(details), c("n_saved", "MCMC_warning", names(flags)))) out$details[[nm]] <- d[[nm]]
    out$details$MCMC_warning <- !is.na(reason)
    if (out$details$MCMC_warning) notes <- c(notes, paste0("Delta MCMC: ", reason))
  }
  out$warning <- out$details$MCMC_warning || length(notes) > 0L
  if (length(notes)) out$message <- paste(unique(notes), collapse = "; ")
  if (is.finite(calc$delta) && calc$delta > 0) {
    out$estimate <- calc$delta; out$status <- "ok"
  } else {
    out$message <- "Delta is non-positive or non-finite"
    out$details$calculation_failed <- TRUE; out$details$failure_stage <- "estimate"
  }
  out
}

.heterogeneity_delta_tests <- function(tree, core, local_contexts, X, nodes, clades,
    overall, control, min_tips, test, nsim, keep_null, clade_test, clade_nsim) {
  active_model <- NULL
  on.exit(if(!is.null(active_model))heterogeneity_mk_release_cpp(active_model),add=TRUE)
  summary <- data.frame(trait = colnames(X), status = "insufficient_clades",
    statistic = NA_real_, test_statistic = "RMS(log(Delta_local / Delta_overall))",
    P_heterogeneity = NA_real_, MCSE_P_heterogeneity = NA_real_,
    null_mean = NA_real_, null_upper = NA_real_, nsim_successful = 0L,
    exceedance_count = NA_integer_, nsim_MCMC_warning = 0L,
    nsim_fit_warning = 0L,
    nsim_optimizer_warning = 0L, nsim_numerical_recovery = 0L,
    nsim_entropy_recovery = 0L, nsim_other_warning = 0L,
    nsim_calculation_failure = 0L, nsim_unavailable = 0L,
    nsim_single_state = 0L, finite_exceedance_count = NA_integer_,
    P_heterogeneity_lower = NA_real_, P_heterogeneity_upper = NA_real_,
    P_heterogeneity_kind = "unavailable",
    reference_unavailable = if (identical(control$reference_unavailable, "conservative")) "conservative" else "strict",
    reference_fit_optimizer_warning = FALSE, reference_fit_numerical_recovery = FALSE,
    reference_fit_attempts = NA_integer_, reference_fit_deviance_improvement = NA_real_,
    message = "", stringsAsFactors = FALSE)
  null <- stats::setNames(vector("list", ncol(X)), colnames(X))
  models <- stats::setNames(vector("list", ncol(X)), colnames(X))
  summary$reference_fit_diagnostics <- I(vector("list", ncol(X)))
  nuisance <- identical(control$reference_test, "nuisance")
  between <- identical(control$statistic,"between_clades")
  if(between)summary$test_statistic <- "SD(log(Delta_local))"
  if (nuisance) {
    summary$P_nuisance_lower <- summary$P_nuisance_upper <- summary$nuisance_gap <- NA_real_
    summary$nuisance_converged <- FALSE
    summary$nuisance_cells <- summary$nuisance_evaluations <- 0L
  }
  canonical_tree <- if (test) tryCatch(.safe_canonicalize_core(tree),
    error = function(e) NULL) else NULL
  common_tree <- if (!is.null(canonical_tree)) tryCatch(.fast_ace_tree_workspace(canonical_tree),
    error = function(e) NULL) else NULL
  for (j in seq_len(ncol(X))) {
    ids <- (seq_along(nodes) - 1L) * ncol(X) + j
    eligible <- which(clades$n_species[ids] >= min_tips)
    if (!length(eligible)) next
    contexts <- lapply(eligible, function(i) {
      z <- local_contexts[[i]]
      if (is.null(z)) NULL else z[[j]]
    })
    observed_x <- X[,j]
    missing <- is.na(observed_x)
    if (clade_test) for (k in seq_along(eligible)) {
      row <- ids[eligible[k]]; ctx <- contexts[[k]]
      if (is.null(ctx) || clades$status[row] != "ok") next
      values <- rep(NA_real_, clade_nsim)
      local_x <- observed_x[ctx$rows]
      for (b in seq_len(clade_nsim)) {
        # Public Delta permutation workers also skip trace diagnostics. Only
        # the crown vector is shuffled; no full-tree vector is copied here.
        values[b] <- .heterogeneity_delta_values(ctx,
          sample(local_x, length(local_x)), diagnostics = FALSE)$estimate
      }
      clades$signal_nsim_successful[row] <- sum(is.finite(values))
      # Match the public Delta signal-test tail; failed draws are explicit and
      # make this probability unavailable rather than changing its reference.
      if (all(is.finite(values))) clades$P_signal[row] <- mean(values >= clades$estimate[row])
    }
    if (!test || length(eligible) < 2L) next
    if (any(clades$status[ids[eligible]] != "ok") || any(vapply(contexts, is.null, logical(1)))) {
      summary$status[j] <- "unavailable_clades"
      summary$message[j] <- "All size-eligible crowns must be estimable; state-based exclusions are not silently frozen for calibration"
      next
    }
    global <- overall$estimate[j]
    if (!is.finite(global) || global <= 0) {
      summary$status[j] <- "unavailable_overall"; next
    }
    observed <- log(clades$estimate[ids[eligible]]) - log(global)
    summary$statistic[j] <- if(between)heterogeneity_delta_between_cpp(clades$estimate[ids[eligible]]) else
      sqrt(mean(observed^2))
    global_ctx <- if (!is.null(canonical_tree)) tryCatch(.heterogeneity_delta_context(canonical_tree,
      tree$tip.label[!missing], tree$tip.label, control, canonicalized = TRUE),
      error = function(e) NULL) else NULL
    if (is.null(global_ctx) || is.null(common_tree)) { summary$status[j] <- "reference_fit_failure"; next }
    if (nuisance) {
      bound <- .heterogeneity_delta_nuisance(tree, core, global_ctx, contexts,
        observed_x, summary$statistic[j], attr(X, "state_labels")[[j]], control, nsim, keep_null)
      for (nm in names(bound$summary)) {
        if (is.null(summary[[nm]])) summary[[nm]] <- 0L
        summary[[nm]][j] <- bound$summary[[nm]]
      }
      models[[j]] <- bound$model
      if (keep_null) null[[j]] <- bound$null
      next
    }
    # Fit Q on the original root/geometry with unobserved tips integrated out.
    # Pruning can move the root; for ARD a uniform mixture at that new root
    # would be a different likelihood from the full-tree simulation reference.
    y <- factor(X[match(common_tree$tip_label, tree$tip.label), j])
    workspace <- .fast_ace_workspace(common_tree, control$model, nlevels(y))
    notes <- character()
    summary$reference_fit_optimizer_warning[j] <- FALSE
    common <- tryCatch(withCallingHandlers(if (identical(control$reference_fit, "profile"))
      .heterogeneity_delta_er_reference(as.integer(y), workspace) else .fast_ace_fit_workspace(
        as.integer(y), levels(y), workspace, ip = .1, CI = FALSE, estimate_se = FALSE,
        fit_strategy = if (identical(control$reference_fit, "multistart")) "multistart" else "default"),
      warning = function(w) {
        if (inherits(w, "fastphylosig_optimizer_warning")) summary$reference_fit_optimizer_warning[j] <<- TRUE
        if (inherits(w, "fastphylosig_numerical_recovery_warning")) summary$reference_fit_numerical_recovery[j] <<- TRUE
        notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning")
      }),
      error = identity)
    if (!inherits(common, "error")) {
      summary$reference_fit_diagnostics[j] <- list(common$diagnostics)
      summary$reference_fit_attempts[j] <- if (is.null(common$diagnostics$attempts)) 1L else nrow(common$diagnostics$attempts)
      summary$reference_fit_deviance_improvement[j] <- if (is.null(common$diagnostics$deviance_improvement)) 0 else common$diagnostics$deviance_improvement
    }
    if (inherits(common, "error") || common$convergence != 0L || !is.finite(common$loglik)) {
      summary$status[j] <- "reference_fit_failure"
      if (inherits(common, "error")) notes <- c(notes, conditionMessage(common))
      summary$message[j] <- paste(notes, collapse = "; "); next
    }
    model <- tryCatch(heterogeneity_mk_prepare_cpp(core, workspace$rate$rate_index,
      common$rates), error = function(e) e)
    if (inherits(model, "error")) {
      summary$status[j] <- "reference_fit_failure"
      summary$message[j] <- conditionMessage(model); next
    }
    active_model <- model
    states <- attr(X, "state_labels")[[j]][as.integer(levels(y))]
    dimnames(common$index.matrix) <- list(states, states)
    models[[j]] <- list(model = control$model, rates = common$rates,
      index_matrix = common$index.matrix, states = states, root = "uniform",
      loglik = common$loglik, convergence = common$convergence)
    if (!is.null(common$diagnostics)) models[[j]]$fit_diagnostics <- common$diagnostics
    if (!is.null(common$numerical_recovery)) {
      models[[j]]$numerical_recovery <- common$numerical_recovery
      summary$message[j] <- paste(unique(notes),collapse="; ")
    }
    draws <- rep(NA_real_, nsim)
    count_fields <- c("nsim_MCMC_warning","nsim_optimizer_warning",
      "nsim_numerical_recovery","nsim_entropy_recovery","nsim_other_warning",
      "nsim_calculation_failure","nsim_unavailable","nsim_single_state")
    counts <- stats::setNames(integer(length(count_fields)),count_fields)
    retries <- c(nsim_optimizer_retry=0L,nsim_optimizer_recovered=0L)
    # Interleave each trait simulation and its MCMC: chunk settings cannot
    # change the serial R RNG stream. Native transitions and scratch are reused.
    for (b in seq_len(nsim)) {
        x <- heterogeneity_mk_simulate_cpp(model, 1L)[,1L]
        x[missing] <- NA_integer_
        full <- .heterogeneity_delta_evaluate(global_ctx, x)
        local <- lapply(contexts, .heterogeneity_delta_evaluate, x = x)
        estimates <- vapply(local, `[[`, numeric(1), "estimate")
        flags <- heterogeneity_delta_flags_cpp(full, local)
        retry <- attr(flags,"optimizer_retry",exact=TRUE)
        if(!is.null(retry)) {
          retries <- retries+retry
          attr(flags,"optimizer_retry") <- NULL
        }
        counts <- counts + flags
        if (is.finite(full$estimate) && all(is.finite(estimates)))
          draws[b] <- sqrt(mean((log(estimates) - log(full$estimate))^2))
    }
    heterogeneity_mk_release_cpp(model);active_model <- NULL
    for (nm in count_fields) summary[[nm]][j] <- counts[[nm]]
    if(retries[["nsim_optimizer_retry"]]>0L) {
      for(nm in names(retries)) {
        if(is.null(summary[[nm]]))summary[[nm]] <- 0L
        summary[[nm]][j] <- retries[[nm]]
      }
    }
    summary$nsim_fit_warning[j] <- summary$nsim_optimizer_warning[j]
    tail <- heterogeneity_bootstrap_summary_cpp(draws, summary$statistic[j], FALSE)
    for (nm in setdiff(names(tail), "null")) summary[[nm]][j] <- tail[[nm]]
    summary$status[j] <- if (tail$nsim_successful == nsim) "ok" else "simulation_failure"
    bounds <- heterogeneity_bootstrap_bounds_cpp(draws, summary$statistic[j])
    for (nm in names(bounds)) summary[[nm]][j] <- bounds[[nm]]
    if (tail$nsim_successful == nsim) summary$P_heterogeneity_kind[j] <- "complete_reference"
    if (summary$status[j] == "simulation_failure") summary$message[j] <-
      "At least one full/local Delta replicate was unavailable; no failed draws were discarded"
    if (summary$reference_unavailable[j] == "conservative" &&
        counts[["nsim_unavailable"]] > 0L && counts[["nsim_calculation_failure"]] == 0L &&
        counts[["nsim_single_state"]] == counts[["nsim_unavailable"]]) {
      summary$P_heterogeneity[j] <- bounds$P_heterogeneity_upper
      q <- (bounds$finite_exceedance_count + counts[["nsim_unavailable"]]) / nsim
      summary$MCSE_P_heterogeneity[j] <- sqrt(nsim * q * (1 - q)) / (nsim + 1)
      summary$status[j] <- "conservative_reference"
      summary$P_heterogeneity_kind[j] <- "single_state_upper_tail"
      summary$message[j] <- "Single-state reference discrepancies counted in the upper tail; Delta remains undefined, no draws discarded or replaced"
    }
    if (keep_null) null[[j]] <- draws
  }
  # Store fitted reference information separately from signal summaries.
  summary$reference_model <- I(models)
  list(clades = clades, summary = summary, null = null)
}

.heterogeneity_delta_er_reference <- function(states, workspace) {
  pointer <- fast_ace_fit_prepare_cpp(workspace$edge, workspace$edge_length,
    states, workspace$rate$rate_index, 1L)
  on.exit(fast_ace_fit_release_cpp(pointer))
  fit <- fast_ace_er_reference_cpp(pointer)
  fit$index.matrix <- workspace$rate$index_matrix
  fit
}
