# Continuous-trait signal API -------------------------------------------------

# The public API always uses the validated tree-pruning production engines.
# Alternative dense/spectral paths remain internal numerical oracles for tests.
# Private continuous-trait production wrapper.  The public entry points are
# `fast_k()` and `fast_lambda()` below; keeping the numerical orchestration in
# one function makes it impossible for the specialist APIs to drift apart.
.fast_signal_continuous <- function(tree, x = NULL,
                        method = c("K", "lambda"),
                        test = FALSE, nsim = 1000, se = NULL,
                        start = NULL, control = list(),
                        permutations = NULL, return_sim = FALSE,
                        verbose = TRUE, ncores = 1, X = NULL,
                        lambda_profile = NULL,
                        lambda_profile_points = 101, prepared = NULL,
                        trait_chunk = 64L,
                        simulation_chunk = 128L,
                        keep_null = return_sim,
                        progress = interactive()) {
  method <- match.arg(method)
  runtime <- .runtime_begin(method, progress = progress, verbose = verbose)
  on.exit(.runtime_on_exit(runtime), add = TRUE)
  out <- .fast_signal_with_engine(
    tree = tree, x = x, method = method, test = test, nsim = nsim,
    se = se, start = start, control = control,
    permutations = permutations, return_sim = return_sim,
    # `progress` is the master status switch; suppress the historical
    # matching notices whenever status output is disabled.
    verbose = isTRUE(verbose) && isTRUE(progress), ncores = ncores, X = X,
    lambda_profile = lambda_profile,
    lambda_profile_points = lambda_profile_points, prepared = prepared,
    engine = "tree", trait_chunk = trait_chunk,
    validate_tolerance = 1e-8, simulation_chunk = simulation_chunk,
    keep_null = keep_null, .runtime = runtime
  )
  timing <- .runtime_close(runtime, success = TRUE)
  .runtime_attach(out, timing)
}

.fast_signal_with_engine <- function(
    tree, x = NULL, method = c("K", "lambda"), test = FALSE, nsim = 1000,
    se = NULL, start = NULL, control = list(), permutations = NULL,
    return_sim = FALSE, verbose = TRUE, ncores = 1, X = NULL,
    lambda_profile = NULL, lambda_profile_points = 101, prepared = NULL,
    engine = NULL, trait_chunk = 64L, validate_tolerance = 1e-8,
    simulation_chunk = 128L, keep_null = return_sim, .runtime = NULL) {
  method <- match.arg(method)
  if (is.null(engine)) {
    engine <- "tree"
  }
  engine <- match.arg(
    engine,
    if (method == "K") c("dense", "tree", "validate") else
      c("spectral", "tree", "validate")
  )
  if (!is.numeric(trait_chunk) || length(trait_chunk) != 1L ||
      !is.finite(trait_chunk) || trait_chunk < 1 ||
      trait_chunk != floor(trait_chunk) ||
      trait_chunk > .Machine$integer.max) {
    stop("trait_chunk must be a finite positive integer.", call. = FALSE)
  }
  trait_chunk <- as.integer(trait_chunk)
  if (!is.numeric(simulation_chunk) || length(simulation_chunk) != 1L ||
      !is.finite(simulation_chunk) || simulation_chunk < 1 ||
      simulation_chunk != floor(simulation_chunk) ||
      simulation_chunk > .Machine$integer.max) {
    stop("simulation_chunk must be a finite positive integer.",
         call. = FALSE)
  }
  simulation_chunk <- as.integer(simulation_chunk)
  if (!is.numeric(validate_tolerance) || length(validate_tolerance) != 1L ||
      !is.finite(validate_tolerance) || validate_tolerance < 0) {
    stop("validate_tolerance must be one finite non-negative value.",
         call. = FALSE)
  }
  if (!is.null(prepared)) tree <- prepared
  ncores <- .normalize_ncores(ncores)
  if (!is.logical(test) || length(test) != 1L || is.na(test)) {
    stop("test must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(return_sim) || length(return_sim) != 1L ||
      is.na(return_sim)) {
    stop("return_sim must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(keep_null) || length(keep_null) != 1L ||
      is.na(keep_null)) {
    stop("keep_null must be TRUE or FALSE.", call. = FALSE)
  }
  store_null <- isTRUE(return_sim) || isTRUE(keep_null)
  if (is.null(x)) {
    x <- X
  }
  if (is.null(x)) {
    stop("x must be a named vector, matrix, or data.frame.", call. = FALSE)
  }
  if (!is.null(se)) {
    stop("fast_signal currently supports se = NULL only.", call. = FALSE)
  }

  vector_input <- is.null(dim(x))
  if (is.null(lambda_profile)) {
    lambda_profile <- method == "lambda" && vector_input
  } else {
    lambda_profile <- isTRUE(lambda_profile)
  }
  if (!is.numeric(lambda_profile_points) ||
      length(lambda_profile_points) != 1L ||
      !is.finite(lambda_profile_points) ||
      lambda_profile_points != floor(lambda_profile_points) ||
      lambda_profile_points < 5 ||
      lambda_profile_points > .Machine$integer.max) {
    stop("lambda_profile_points must be a finite integer >= 5.",
         call. = FALSE)
  }
  lambda_profile_points <- as.integer(lambda_profile_points)

  out <- .fast_signal_batch(
    tree = tree, X = x, method = method, test = test, nsim = nsim,
    permutations = permutations, return_sim = store_null, verbose = verbose,
    ncores = ncores, lambda_profile = lambda_profile,
    lambda_profile_points = lambda_profile_points, engine = engine,
    trait_chunk = trait_chunk, validate_tolerance = validate_tolerance,
    simulation_chunk = simulation_chunk, .runtime = .runtime
  )
  analysis_metadata <- attr(out, "analysis_metadata", exact = TRUE)

  if (!vector_input) {
    result <- .decorate_fastphylosig_result(
      out, method = method, vector_input = FALSE
    )
    attr(result, "analysis_metadata") <- analysis_metadata
    return(result)
  }

  if (method == "K") {
    object <- if (!test) {
      out$K_fast[[1]]
    } else {
      z <- list(K = out$K_fast[[1]], P = out$P_fast[[1]])
      z$status <- out$status[[1]]
      z$message <- out$message[[1]]
      z$nsim <- out$nsim_requested[[1]]
      z$nsim_requested <- out$nsim_requested[[1]]
      z$nsim_successful <- out$nsim_successful[[1]]
      z$MCSE_P <- out$MCSE_P[[1]]
      z$exceedance_count <- out$exceedance_count[[1]]
      if (store_null && "sim.K_fast" %in% names(out)) {
        z$sim.K <- out$sim.K_fast[[1]]
      }
      z
    }
  } else {
    object <- list(
      lambda = out$lambda_fast[[1]],
      logL = out$logL_fast[[1]],
      gls_mean = out$gls_mean_fast[[1]],
      sig2 = out$sig2_fast[[1]],
      note = out$note[[1]],
      status = out$status[[1]],
      message = out$message[[1]]
    )
    if (test) {
      object$logL0 <- out$logL0_fast[[1]]
      object$LR <- out$LR_fast[[1]]
      object$P <- out$P_fast[[1]]
    }
    if (lambda_profile && "lambda_profile_fast" %in% names(out)) {
      object$lambda_profile <- out$lambda_profile_fast[[1]]
      object$lambda_CI <- c(
        lower = out$lambda_CI_lower_fast[[1]],
        upper = out$lambda_CI_upper_fast[[1]]
      )
      object$lambda_CI_level <- 0.95
      object$lambda_CI_cutoff <- out$lambda_CI_cutoff_fast[[1]]
    }
  }

  attr(object, "class") <- "phylosig"
  attr(object, "method") <- method
  attr(object, "test") <- test
  attr(object, "se") <- FALSE
  attr(object, "engine") <- engine
  if (method == "K" && !test) {
    attr(object, "status") <- out$status[[1]]
    attr(object, "message") <- out$message[[1]]
  }
  if (method == "K" && engine == "validate") {
    attr(object, "validation") <- out[1L, c(
      "K_tree", "K_dense", "absolute_error", "relative_error",
      "validation_pass"
    ), drop = FALSE]
  }
  if (method == "lambda" && engine == "validate") {
    attr(object, "validation") <- out[1L, c(
      "lambda_tree", "lambda_spectral", "lambda_absolute_error",
      "logL_tree", "logL_spectral", "logL_absolute_error",
      "LR_tree", "LR_spectral", "LR_absolute_error",
      "validation_pass"
    ), drop = FALSE]
  }
  result <- .decorate_fastphylosig_result(
    object, method = method, vector_input = TRUE
  )
  attr(result, "analysis_metadata") <- analysis_metadata
  result
}


# Public specialist entry points and the unified dispatcher ------------------

# `fast_signal()` deliberately does not infer a method from the input.  This
# helper accepts the small, documented alias set and rejects both missing and
# vector-valued methods before any tree/data preparation starts.
.normalize_fast_signal_method <- function(method) {
  allowed <- "K/k, lambda/Lambda/\u03bb, D/d, or Delta/delta"
  if (is.null(method) || length(method) != 1L || is.na(method)) {
    stop(
      sprintf("method must be exactly one of %s.", allowed),
      call. = FALSE
    )
  }
  key <- tryCatch(as.character(method), error = function(e) character())
  if (length(key) != 1L || is.na(key) || !nzchar(trimws(key))) {
    stop(
      sprintf("method must be exactly one of %s.", allowed),
      call. = FALSE
    )
  }
  key <- trimws(key)
  lower <- tolower(key)
  if (identical(key, "\u03bb") || identical(lower, "lambda")) {
    return("lambda")
  }
  if (identical(lower, "k")) return("K")
  if (identical(lower, "d")) return("D")
  if (identical(lower, "delta")) return("Delta")
  stop(
    sprintf("Unknown method '%s'; method must be exactly one of %s.",
            key, allowed),
    call. = FALSE
  )
}

.resolve_fast_signal_data <- function(x = NULL, X = NULL, data = NULL) {
  supplied <- !vapply(list(x, X, data), is.null, logical(1))
  if (sum(supplied) > 1L) {
    stop("Supply only one of x, X, or data.", call. = FALSE)
  }
  if (!is.null(data)) return(data)
  if (!is.null(x)) return(x)
  X
}

.attach_fast_signal_workflow <- function(x, method, production_function) {
  # Keep the numerical result and its legacy classes untouched.  Workflow
  # metadata lives in attributes so scalar K results remain atomic numerics and
  # the D/Delta caper-compatible objects remain their original classes.
  attr(x, "method") <- method
  attr(x, "production_function") <- production_function
  workflow <- attr(x, "workflow", exact = TRUE)
  if (!is.list(workflow)) workflow <- list()
  workflow$method <- method
  workflow$production_function <- production_function
  analysis_metadata <- attr(x, "analysis_metadata", exact = TRUE)
  if (is.list(analysis_metadata)) {
    workflow$input_tree_type <- analysis_metadata$input_tree_type
    workflow$raw_or_prepared <- analysis_metadata$raw_or_prepared
    workflow$tree_auto_normalized <-
      isTRUE(analysis_metadata$tree_auto_normalized)
    workflow$matching_performed <- TRUE
    workflow$tree_tips_removed <-
      length(analysis_metadata$tree_tips_removed)
    workflow$data_rows_removed <-
      length(analysis_metadata$data_rows_removed)
    workflow$NA_patterns <- analysis_metadata$na_pattern_count
    attr(x, "matching") <- analysis_metadata$matching
    attr(x, "tree_processing") <- list(
      input_tree_type = analysis_metadata$input_tree_type,
      tree_auto_normalized = analysis_metadata$tree_auto_normalized,
      canonical_changes = analysis_metadata$canonical_changes,
      retained_tip_validation = analysis_metadata$retained_tip_validation
    )
  }
  timing <- attr(x, "timing", exact = TRUE)
  workflow$total_elapsed <- if (is.list(timing)) {
    as.numeric(timing$total_elapsed)
  } else NA_real_
  attr(x, "workflow") <- workflow
  x
}

.dispatch_fast_signal_dots <- function(dots, target, method, x = NULL,
                                       allowed = NULL) {
  if (length(dots)) {
    dot_names <- names(dots)
    if (is.null(dot_names)) dot_names <- rep("", length(dots))
    dot_names[is.na(dot_names)] <- ""
  } else {
    dot_names <- character()
  }

  # `data` is the generic dispatcher spelling; the production functions keep
  # their historical x/X aliases.  Resolve it before validating method-
  # specific arguments so duplicate aliases fail deterministically.
  data_idx <- which(dot_names == "data")
  if (length(data_idx) > 1L) {
    stop("Supply data only once.", call. = FALSE)
  }
  data <- if (length(data_idx)) dots[[data_idx[[1L]]]] else NULL
  if (length(data_idx)) {
    dots <- dots[-data_idx[[1L]]]
    dot_names <- dot_names[-data_idx[[1L]]]
  }

  X_idx <- which(dot_names == "X")
  if (length(X_idx) > 1L) {
    stop("Supply X only once.", call. = FALSE)
  }
  X <- if (length(X_idx)) dots[[X_idx[[1L]]]] else NULL
  if (length(X_idx)) {
    dots <- dots[-X_idx[[1L]]]
    dot_names <- dot_names[-X_idx[[1L]]]
  }
  input <- .resolve_fast_signal_data(x = x, X = X, data = data)

  if (is.null(allowed)) {
    target_names <- names(formals(target))
    allowed <- setdiff(target_names, c("tree", "x", "method"))
  }
  named <- dot_names != ""
  used <- dot_names[named]
  unknown_named <- setdiff(used, allowed)
  if (length(unknown_named)) {
    stop(
      sprintf(
        "Arguments not applicable to method '%s': %s.", method,
        paste(unique(unknown_named), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Preserve the legacy positional API (test, nsim, ...) while still keeping
  # strict dots.  Unnamed values are assigned to the first target formals not
  # already supplied by name; all other names are rejected above.
  unnamed <- which(!named)
  if (length(unnamed)) {
    available <- setdiff(allowed, used)
    if (length(unnamed) > length(available)) {
      stop(
        sprintf("Too many unnamed arguments for method '%s'.", method),
        call. = FALSE
      )
    }
    dot_names[unnamed] <- available[seq_along(unnamed)]
  }
  if (anyDuplicated(dot_names)) {
    duplicated_names <- unique(dot_names[duplicated(dot_names)])
    stop(
      sprintf("Arguments supplied more than once: %s.",
              paste(duplicated_names, collapse = ", ")),
      call. = FALSE
    )
  }
  names(dots) <- dot_names
  list(x = input, args = dots)
}

# The specialist signatures expose only arguments that affect their method.
# `data` is a named alias for x; all numerical work remains in
# `.fast_signal_continuous()`.
fast_k <- function(tree, x = NULL, test = FALSE, nsim = 1000, se = NULL,
                   permutations = NULL, return_sim = FALSE, verbose = TRUE,
                   ncores = 1, X = NULL, prepared = NULL, trait_chunk = 64L,
                   simulation_chunk = 128L, keep_null = return_sim,
                   progress = interactive(), data = NULL) {
  x <- .resolve_fast_signal_data(x = x, X = X, data = data)
  out <- .fast_signal_continuous(
    tree = tree, x = x, method = "K", test = test, nsim = nsim, se = se,
    start = NULL, control = list(), permutations = permutations,
    return_sim = return_sim, verbose = verbose, ncores = ncores, X = NULL,
    lambda_profile = NULL, lambda_profile_points = 101L, prepared = prepared,
    trait_chunk = trait_chunk, simulation_chunk = simulation_chunk,
    keep_null = keep_null, progress = progress
  )
  .attach_fast_signal_workflow(out, "K", "fast_k")
}

fast_lambda <- function(tree, x = NULL, test = FALSE, se = NULL,
                        verbose = TRUE, ncores = 1, X = NULL,
                        lambda_profile = NULL,
                        lambda_profile_points = 101, prepared = NULL,
                        trait_chunk = 64L, progress = interactive(),
                        data = NULL) {
  x <- .resolve_fast_signal_data(x = x, X = X, data = data)
  out <- .fast_signal_continuous(
    tree = tree, x = x, method = "lambda", test = test, nsim = 1000L,
    se = se, start = NULL, control = list(), permutations = NULL,
    return_sim = FALSE, verbose = verbose, ncores = ncores, X = NULL,
    lambda_profile = lambda_profile,
    lambda_profile_points = lambda_profile_points, prepared = prepared,
    trait_chunk = trait_chunk, simulation_chunk = 128L,
    keep_null = FALSE, progress = progress
  )
  .attach_fast_signal_workflow(out, "lambda", "fast_lambda")
}

# Unified public entry point.  It performs dispatch only; progress/runtime
# contexts are created exactly once by the selected production target.
fast_signal <- function(tree, data = NULL, method = NULL, ..., x = NULL,
                        X = NULL, prepared = NULL, trait_chunk = 64L) {
  method <- .normalize_fast_signal_method(method)
  if (!is.null(data) && (!is.null(x) || !is.null(X))) {
    stop("Supply only one of data, x, or X.", call. = FALSE)
  }
  input <- .resolve_fast_signal_data(x = x, X = X, data = data)
  target <- switch(method,
    K = .fast_signal_continuous,
    lambda = .fast_signal_continuous,
    D = fast_d,
    Delta = fast_delta
  )
  dots <- list(...)
  allowed <- switch(method,
    K = c("test", "nsim", "se", "permutations", "return_sim",
          "verbose", "ncores", "prepared", "trait_chunk",
          "simulation_chunk", "keep_null", "progress"),
    lambda = c("test", "se", "verbose", "ncores", "lambda_profile",
               "lambda_profile_points", "prepared", "trait_chunk",
               "progress"),
    D = setdiff(names(formals(fast_d)), c("tree", "x", "X")),
    Delta = setdiff(names(formals(fast_delta)), c("tree", "x", "X"))
  )
  # Keep the long-standing prepared/trait_chunk formals visible for callers
  # that inspect the API, while forwarding them only when explicitly supplied.
  # (The default trait chunk is owned by the selected target.)
  if (!missing(prepared)) dots$prepared <- prepared
  if (!missing(trait_chunk)) dots$trait_chunk <- trait_chunk
  captured <- .dispatch_fast_signal_dots(
    dots, target = target, method = method, x = input, allowed = allowed
  )
  call <- c(
    list(tree = tree, x = captured$x),
    if (method %in% c("K", "lambda")) list(method = method) else list(),
    captured$args
  )
  out <- do.call(target, call)
  .attach_fast_signal_workflow(
    out, method = method,
    production_function = switch(method,
      K = "fast_k", lambda = "fast_lambda", D = "fast_d", Delta = "fast_delta"
    )
  )
}


# One internal batch engine ----------------------------------------------------

.fast_signal_batch <- function(tree, X, method, test, nsim, permutations,
                               return_sim, verbose, ncores,
                               lambda_profile, lambda_profile_points, engine,
                               trait_chunk, validate_tolerance,
                               simulation_chunk, .runtime = NULL,
                               .analysis = NULL) {
  .runtime_stage(.runtime, "Checking tree...")
  .runtime_stage(.runtime, "Preparing data...")
  if (is.null(.analysis)) {
    # Shared preparation owns representation checks, safe canonicalization,
    # species matching, and NA-mask grouping.  Keep it in one place so the
    # K/lambda production path cannot silently diverge from the public
    # workflow diagnostics.
    .analysis <- .prepare_analysis(
      tree = tree, data = X, signal = method, data_kind = "continuous",
      verbose = FALSE
    )
  } else {
    if (!inherits(.analysis, "fastphylosig_analysis_preparation") ||
        !identical(.analysis$signal, method)) {
      stop(
        "analysis must be a fastphylosig_analysis_preparation for the selected method.",
        call. = FALSE
      )
    }
  }
  ctx <- .analysis$ctx
  tree <- ctx$tree
  if (method == "K" && test) {
    if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
        nsim < 1 || nsim != floor(nsim)) {
      stop("nsim must be a positive integer when test = TRUE.",
           call. = FALSE)
    }
    nsim <- as.integer(nsim)
  }

  if (isTRUE(.analysis$insufficient_retained)) {
    stop(
      "Fewer than two species are shared by tree and x; retain at least two matching species.",
      call. = FALSE
    )
  }
  # `.prepare_analysis()` keeps the matched table in tree-tip order.  The
  # continuous validator is still applied here to preserve the historical
  # numeric/finite-value contract without repeating matching or NA grouping.
  X_raw <- .as_trait_matrix(.analysis$matched_data, tree, verbose = FALSE)
  tree_tips <- tree$tip.label
  matching_details <- .analysis$matching_details
  matched <- as.character(matching_details$matched_species)
  removed_tree <- as.character(matching_details$tree_tips_removed)
  removed_data <- as.character(matching_details$data_rows_removed)

  base_keep <- match(matched, tree_tips)
  need_dense <- !identical(engine, "tree")
  need_lambda_cache <- identical(method, "lambda") &&
    !identical(engine, "tree")
  base_group <- .prepared_tree_subset(
    ctx, base_keep, need_lambda = need_lambda_cache,
    need_matrix = need_dense
  )
  base_tree <- base_group$tree
  base_idx <- match(base_tree$tip.label, tree_tips)
  X0 <- X_raw[base_tree$tip.label, , drop = FALSE]
  n0 <- nrow(X0)
  traits <- colnames(X0)
  p <- ncol(X0)

  .runtime_stage(.runtime, if (identical(method, "K")) {
    "Calculating K..."
  } else {
    "Optimizing lambda..."
  })
  if (identical(method, "K") && isTRUE(test)) {
    .runtime_stage(.runtime, sprintf("Running %d randomizations...", nsim))
  }
  if (identical(method, "lambda") && isTRUE(test)) {
    .runtime_stage(.runtime, "Testing lambda = 0...")
  }

  present <- !is.na(X0)
  n_species <- as.integer(colSums(present))
  n_removed_na <- n0 - n_species
  # Reuse the shared NA-mask grouping exactly; this also preserves the
  # historical per-trait <2-retained-species note/NA behavior.
  mask_groups <- .analysis$na_patterns

  if (isTRUE(verbose) && (length(removed_tree) || length(removed_data))) {
    message(sprintf(
      paste(
        "Matched %d species; removed %d tree tips and %d data rows before",
        "trait-wise NA pruning."
      ),
      length(matched), length(removed_tree), length(removed_data)
    ))
  }

  K_fast <- P_fast <- lambda_fast <- logL_fast <- logL0_fast <-
    rep(NA_real_, p)
  lambda_gls_mean_fast <- lambda_sig2_fast <- lambda_LR_fast <-
    rep(NA_real_, p)
  K_MCSE_P <- rep(NA_real_, p)
  K_exceedance_count <- rep(NA_real_, p)
  K_nsim_requested <- rep(if (method == "K" && test) nsim else NA_integer_, p)
  K_nsim_successful <- rep(NA_integer_, p)
  K_tree <- K_dense <- K_absolute_error <- K_relative_error <-
    rep(NA_real_, p)
  K_validation_pass <- rep(NA, p)
  lambda_tree <- lambda_spectral <- lambda_absolute_error <-
    rep(NA_real_, p)
  logL_tree <- logL_spectral <- logL_absolute_error <- rep(NA_real_, p)
  logL0_tree <- P_tree <- lambda_gls_mean_tree <- lambda_sig2_tree <-
    rep(NA_real_, p)
  LR_tree <- LR_spectral <- LR_absolute_error <- rep(NA_real_, p)
  lambda_validation_pass <- rep(NA, p)
  lambda_CI_lower_fast <- lambda_CI_upper_fast <- lambda_CI_cutoff_fast <-
    rep(NA_real_, p)
  sim_list <- vector("list", p)
  lambda_profile_list <- vector("list", p)
  notes <- rep(NA_character_, p)
  status <- rep("ok", p)
  nonfinite_trait <- vapply(
    seq_len(p),
    function(j) any(is.infinite(X0[, j]) | is.nan(X0[, j])),
    logical(1)
  )
  if (any(nonfinite_trait)) {
    status[nonfinite_trait] <- "invalid_trait"
    notes[nonfinite_trait] <- "trait contains non-finite values other than NA"
  }
  for (group_index in seq_len(mask_groups$n_group)) {
    analysis_group <- if (length(.analysis$groups) >= group_index) {
      .analysis$groups[[group_index]]
    } else NULL
    idx <- if (is.list(analysis_group) && !is.null(analysis_group$columns)) {
      as.integer(analysis_group$columns)
    } else {
      mask_groups$columns[[group_index]]
    }
    idx <- idx[!nonfinite_trait[idx]]
    if (!length(idx)) next
    keep <- if (is.list(analysis_group) && !is.null(analysis_group$keep)) {
      as.integer(analysis_group$keep)
    } else {
      mask_groups$keep[[group_index]]
    }
    if (length(keep) < 2) {
      notes[idx] <- "fewer than 2 non-NA matched species"
      status[idx] <- "insufficient_data"
      next
    }

    group <- .prepared_tree_subset(
      ctx, base_idx[keep], need_lambda = need_lambda_cache,
      need_matrix = need_dense
    )
    group_tree <- group$tree
    C <- group$C
    Xg <- X0[group_tree$tip.label, idx, drop = FALSE]

    if (method == "K") {
      tree_result <- NULL
      if (engine %in% c("tree", "validate")) {
        if (!test) {
          tree_result <- fast_k_tree_batch_cpp(
            group$compiled_tree, Xg, trait_chunk = trait_chunk
          )
        } else {
          tree_permutations <- permutations
          if (engine == "validate" && is.null(tree_permutations)) {
            # Development validation deliberately materializes one controlled
            # sequence so tree and dense oracles receive identical draws.
            tree_permutations <- .permutation_matrix(
              NULL, n = nrow(Xg), nsim = nsim, include_observed = TRUE
            )
          } else if (!is.null(tree_permutations)) {
            tree_permutations <- .permutation_matrix(
              tree_permutations, n = nrow(Xg), nsim = nsim,
              include_observed = FALSE
            )
          }
          tree_result <- fast_k_tree_permutation_cpp(
            compiled_tree = group$compiled_tree, X = Xg, nsim = nsim,
            permutations = tree_permutations, trait_chunk = trait_chunk,
            simulation_chunk = simulation_chunk,
            return_sim = return_sim,
            include_observed = is.null(tree_permutations),
            n_threads = ncores
          )
          P_fast[idx] <- as.numeric(tree_result$P)
          K_MCSE_P[idx] <- as.numeric(tree_result$MCSE_P)
          K_exceedance_count[idx] <- as.numeric(
            tree_result$exceedance_count
          )
          K_nsim_successful[idx] <- as.integer(
            tree_result$nsim_successful
          )
          if (return_sim && !is.null(tree_result$sim_K)) {
            for (j in seq_along(idx)) {
              sim_list[[idx[[j]]]] <- tree_result$sim_K[, j]
            }
          }
        }
        K_tree[idx] <- as.numeric(tree_result$K)
        invalid_tree <- !is.finite(K_tree[idx])
        if (any(invalid_tree)) {
          notes[idx[invalid_tree]] <-
            "K is undefined for a constant or numerically degenerate trait"
          status[idx[invalid_tree]] <- "undefined"
        }
      }
      if (engine == "tree") {
        K_fast[idx] <- K_tree[idx]
        next
      }
      cholC <- group$chol
      if (is.null(cholC)) {
        stop("VCV matrix is not positive definite for K calculation.",
             call. = FALSE)
      }
      traceC <- group$traceC

      if (!test) {
        K_fast[idx] <- fast_k_chol_batch_cpp(Xg, cholC, traceC)
      } else {
        n <- nrow(C)
        dense_permutations <- if (engine == "validate") {
          tree_permutations
        } else {
          .permutation_matrix(
            permutations, n = n, nsim = nsim, include_observed = TRUE
          )
        }
        kres <- if (return_sim) {
          fast_k_chol_permutation_cpp(
            Xg, cholC, traceC, dense_permutations, n_threads = ncores
          )
        } else {
          fast_k_chol_permutation_p_cpp(
            Xg, cholC, traceC, dense_permutations, n_threads = ncores
          )
        }
        K_fast[idx] <- as.numeric(kres$K)
        dense_P <- as.numeric(kres$P)
        if (engine != "validate") P_fast[idx] <- dense_P
        dense_success <- ifelse(
          is.finite(K_fast[idx]) & is.finite(dense_P), nsim, 0L
        )
        dense_count <- dense_P * nsim
        dense_mcse <- .k_permutation_mcse(
          exceedance = dense_count, nsim = nsim,
          include_observed = is.null(permutations)
        )
        if (engine != "validate") {
          K_exceedance_count[idx] <- dense_count
          K_nsim_successful[idx] <- as.integer(dense_success)
          K_MCSE_P[idx] <- dense_mcse
        }
        if (return_sim) {
          if (engine != "validate") {
            for (j in seq_along(idx)) {
              sim_list[[idx[[j]]]] <- kres$sim_K[, j]
            }
          }
        }
      }
      K_dense[idx] <- K_fast[idx]
      if (engine == "validate") {
        K_absolute_error[idx] <- abs(K_tree[idx] - K_dense[idx])
        K_relative_error[idx] <- K_absolute_error[idx] /
          pmax(abs(K_dense[idx]), .Machine$double.eps)
        K_validation_pass[idx] <- K_absolute_error[idx] <=
          validate_tolerance
        K_fast[idx] <- K_tree[idx]
      }
    } else {
      maxlam <- .max_lambda(group_tree)
      if (engine %in% c("tree", "validate")) {
        tree_fit <- fast_lambda_tree_optimize_cpp(
          compiled_tree = group$compiled_tree,
          X = Xg,
          max_lambda = maxlam,
          test = test,
          profile = lambda_profile,
          profile_points = lambda_profile_points,
          trait_chunk = trait_chunk,
          n_threads = ncores
        )
        lambda_tree[idx] <- as.numeric(tree_fit$lambda)
        logL_tree[idx] <- as.numeric(tree_fit$logLik)
        logL0_tree[idx] <- as.numeric(tree_fit$logLik0)
        LR_tree[idx] <- as.numeric(tree_fit$LR)
        P_tree[idx] <- as.numeric(tree_fit$P)
        lambda_gls_mean_tree[idx] <- as.numeric(tree_fit$gls_mean)
        lambda_sig2_tree[idx] <- as.numeric(tree_fit$sigma2)

        invalid <- !as.logical(tree_fit$valid)
        if (any(invalid)) {
          notes[idx[invalid]] <- as.character(tree_fit$status[invalid])
          status[idx[invalid]] <- "undefined"
        }
        if (lambda_profile && !is.null(tree_fit$profile_lambda) &&
            !is.null(tree_fit$profile_logLik)) {
          profile_lambda <- as.numeric(tree_fit$profile_lambda)
          profile_logL <- as.matrix(tree_fit$profile_logLik)
          for (j in seq_along(idx)) {
            col <- idx[[j]]
            profile <- data.frame(
              lambda = profile_lambda,
              logL = profile_logL[, j],
              stringsAsFactors = FALSE
            )
            ci <- .lambda_profile_ci(
              profile = profile,
              lambda_hat = lambda_tree[[col]],
              logL_hat = logL_tree[[col]]
            )
            lambda_profile_list[[col]] <- profile
            lambda_CI_lower_fast[col] <- ci[["lower"]]
            lambda_CI_upper_fast[col] <- ci[["upper"]]
            lambda_CI_cutoff_fast[col] <- ci[["cutoff"]]
          }
        }
        if (engine == "tree") {
          lambda_fast[idx] <- lambda_tree[idx]
          logL_fast[idx] <- logL_tree[idx]
          logL0_fast[idx] <- logL0_tree[idx]
          lambda_LR_fast[idx] <- LR_tree[idx]
          P_fast[idx] <- P_tree[idx]
          lambda_gls_mean_fast[idx] <- lambda_gls_mean_tree[idx]
          lambda_sig2_fast[idx] <- lambda_sig2_tree[idx]
          next
        }
      }

      spectral <- group$lambda_spectral
      lik_batch <- NULL
      wproj <- zproj <- NULL
      spectral_ok <- rep(FALSE, ncol(Xg))
      if (!is.null(spectral)) {
        # Project every trait into the cached eigenbasis once per NA pattern.
        # Subsequent optimizer evaluations only touch O(n) vectors.
        wproj <- as.numeric(crossprod(
          spectral$vectors, spectral$inv_sqrt_diag
        ))
        scaled_Xg <- sweep(
          Xg, 1L, spectral$inv_sqrt_diag, FUN = "*"
        )
        zproj <- crossprod(spectral$vectors, scaled_Xg)
        centered_Xg <- sweep(Xg, 2L, colMeans(Xg), FUN = "-")
        relative_sd <- sqrt(colSums(centered_Xg * centered_Xg) /
                              nrow(Xg)) /
          pmax(1, apply(abs(Xg), 2L, max))
        # Near-constant traits suffer cancellation after the GLS mean is
        # projected into the eigenbasis. Dense likelihood is the same
        # phytools formula and is numerically safer for those columns.
        spectral_ok <- is.finite(relative_sd) & relative_sd > 1e-8
        if (any(!spectral_ok)) {
          notes[idx[!spectral_ok]] <-
            "dense lambda likelihood fallback for near-constant trait"
        }
        lik_batch <- function(lambda) {
          out_batch <- lambda_loglik_spectral_projected_batch_cpp(
            lambda = lambda,
            eigval = spectral$values,
            wproj = wproj,
            zproj = zproj,
            log_diag = spectral$log_diag,
            n_threads = ncores
          )
          bad <- which(!spectral_ok)
          if (length(bad)) {
            for (j_bad in bad) {
              out_batch[, j_bad] <- vapply(
                lambda, .lambda_loglik_reference, numeric(1), C = C,
                y = Xg[, j_bad]
              )
            }
          }
          out_batch
        }
      }
      lik_one <- function(lambda, y, trait_index = NULL) {
        use_spectral <- !is.null(spectral) && !is.null(trait_index) &&
          isTRUE(spectral_ok[[trait_index]])
        if (!use_spectral) {
          .lambda_loglik_reference(lambda, C, y)
        } else {
          lambda_loglik_spectral_projected_cpp(
            lambda = lambda,
            eigval = spectral$values,
            wproj = wproj,
            zproj = zproj[, trait_index, drop = TRUE],
            log_diag = spectral$log_diag
          )
        }
      }
      for (j in seq_along(idx)) {
        col <- idx[[j]]
        y <- Xg[, j]
        lik <- function(lambda) lik_one(lambda, y, j)
        opt <- stats::optimize(lik, c(0, maxlam), maximum = TRUE)
        if (!is.finite(opt$objective)) {
          notes[col] <- "lambda likelihood is undefined on this trait/tree"
          status[col] <- "undefined"
          next
        }
        lambda_fast[col] <- opt$maximum
        logL_fast[col] <- opt$objective
        components <- .lambda_components_reference(
          lambda = opt$maximum, C = C, y = y
        )
        lambda_gls_mean_fast[col] <- components$gls_mean
        lambda_sig2_fast[col] <- components$sigma2
        if (test) {
          logL0_fast[col] <- lik(0)
          if (is.finite(logL0_fast[col])) {
            lambda_LR_fast[col] <- max(
              0, 2 * (opt$objective - logL0_fast[col])
            )
            P_fast[col] <- stats::pchisq(
              lambda_LR_fast[col],
              df = 1, lower.tail = FALSE
            )
          }
        }
      }
      if (lambda_profile && engine != "validate") {
        hats <- lambda_fast[idx]
        fits <- logL_fast[idx]
        grid_profiles <- .lambda_profile_data_batch(
          max_lambda = maxlam, lambda_hat = hats, logL_hat = fits,
          n_points = lambda_profile_points, batch_lik = lik_batch,
          scalar_lik = function(lambda, y) lik_one(lambda, y), Y = Xg
        )
        for (j in seq_along(idx)) {
          col <- idx[[j]]
          profile <- grid_profiles[[j]]
          ci <- .lambda_profile_ci(
            profile = profile, lambda_hat = lambda_fast[[col]],
            logL_hat = logL_fast[[col]]
          )
          lambda_profile_list[[col]] <- profile
          lambda_CI_lower_fast[col] <- ci[["lower"]]
          lambda_CI_upper_fast[col] <- ci[["upper"]]
          lambda_CI_cutoff_fast[col] <- ci[["cutoff"]]
        }
      }
      if (engine == "validate") {
        lambda_spectral[idx] <- lambda_fast[idx]
        logL_spectral[idx] <- logL_fast[idx]
        LR_spectral[idx] <- lambda_LR_fast[idx]
        lambda_absolute_error[idx] <- abs(
          lambda_tree[idx] - lambda_spectral[idx]
        )
        logL_absolute_error[idx] <- abs(logL_tree[idx] - logL_spectral[idx])
        LR_absolute_error[idx] <- abs(LR_tree[idx] - LR_spectral[idx])
        lambda_validation_pass[idx] <-
          lambda_absolute_error[idx] <= validate_tolerance &
          logL_absolute_error[idx] <= validate_tolerance &
          (is.na(LR_absolute_error[idx]) |
             LR_absolute_error[idx] <= validate_tolerance)
        lambda_fast[idx] <- lambda_tree[idx]
        logL_fast[idx] <- logL_tree[idx]
        logL0_fast[idx] <- logL0_tree[idx]
        lambda_LR_fast[idx] <- LR_tree[idx]
        P_fast[idx] <- P_tree[idx]
        lambda_gls_mean_fast[idx] <- lambda_gls_mean_tree[idx]
        lambda_sig2_fast[idx] <- lambda_sig2_tree[idx]
      }
    }
  }

  common <- data.frame(
    trait = traits,
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
  out <- if (method == "K") {
    data.frame(trait = traits, K_fast = K_fast, common[-1],
               stringsAsFactors = FALSE)
  } else {
    data.frame(trait = traits, lambda_fast = lambda_fast,
               logL_fast = logL_fast,
               gls_mean_fast = lambda_gls_mean_fast,
               sig2_fast = lambda_sig2_fast,
               common[-1], stringsAsFactors = FALSE)
  }

  if (method == "K" && test) {
    out$P_fast <- P_fast
    out$MCSE_P <- K_MCSE_P
    out$nsim <- K_nsim_requested
    out$nsim_requested <- K_nsim_requested
    out$nsim_successful <- K_nsim_successful
    out$exceedance_count <- K_exceedance_count
    if (return_sim) out$sim.K_fast <- I(sim_list)
  }
  if (method == "K" && engine == "validate") {
    out$K_tree <- K_tree
    out$K_dense <- K_dense
    out$absolute_error <- K_absolute_error
    out$relative_error <- K_relative_error
    out$validation_pass <- K_validation_pass
  }
  if (method == "lambda" && test) {
    out$logL0_fast <- logL0_fast
    out$LR_fast <- lambda_LR_fast
    out$P_fast <- P_fast
  }
  if (method == "lambda" && lambda_profile) {
    out$lambda_CI_lower_fast <- lambda_CI_lower_fast
    out$lambda_CI_upper_fast <- lambda_CI_upper_fast
    out$lambda_CI_cutoff_fast <- lambda_CI_cutoff_fast
    out$lambda_profile_fast <- I(lambda_profile_list)
  }
  if (method == "lambda" && engine == "validate") {
    out$lambda_tree <- lambda_tree
    out$lambda_spectral <- lambda_spectral
    out$lambda_absolute_error <- lambda_absolute_error
    out$logL_tree <- logL_tree
    out$logL_spectral <- logL_spectral
    out$logL_absolute_error <- logL_absolute_error
    out$LR_tree <- LR_tree
    out$LR_spectral <- LR_spectral
    out$LR_absolute_error <- LR_absolute_error
    out$validation_pass <- lambda_validation_pass
  }

  attr(out, "match_report") <- list(
    original_tree_tips = ape::Ntip(tree),
    input_rows = as.integer(.analysis$matching$input_rows[[1L]]),
    matched_species = length(matched),
    removed_tree_tips = length(removed_tree),
    removed_data_rows = length(removed_data),
    tree_tips_removed = removed_tree,
    data_rows_removed = removed_data
  )
  attr(out, "engine") <- engine
  # Keep only compact provenance on the numerical table; never expose the
  # prepared context/cache or matched data through a result attribute.
  attr(out, "analysis_metadata") <- .compact_analysis_metadata(.analysis)
  out
}

.max_lambda <- function(tree) {
  if (!ape::is.ultrametric(tree)) {
    return(1)
  }
  depth <- ape::node.depth.edgelength(tree)
  parent_height <- depth[tree$edge[, 1]]
  child_height <- depth[tree$edge[, 2]]
  denom <- max(parent_height)
  if (!is.finite(denom) || denom <= 0) return(1)
  out <- max(child_height) / denom
  if (!is.finite(out) || out <= 0) 1 else out
}

# Dense fallback for numerically near-constant traits. This mirrors the
# reference likelihood in phytools::phylosig() and is used only when the
# eigenbasis projection would lose meaningful residual precision.
.lambda_components_reference <- function(lambda, C, y) {
  invalid <- list(gls_mean = NA_real_, sigma2 = NA_real_, logL = -Inf)
  n <- nrow(C)
  Cl <- lambda * C
  diag(Cl) <- diag(C)
  invCl <- tryCatch(solve(Cl), error = function(e) NULL)
  if (is.null(invCl)) return(invalid)
  sum_inv <- sum(invCl)
  if (!is.finite(sum_inv) || sum_inv == 0) return(invalid)
  a <- sum(invCl %*% y) / sum_inv
  centered <- y - a
  quad <- drop(t(centered) %*% invCl %*% centered)
  sig2 <- quad / n
  if (!is.finite(sig2) || sig2 <= 0) {
    invalid$gls_mean <- as.numeric(a)
    return(invalid)
  }
  logL <- tryCatch(
    -drop(t(centered) %*% ((1 / sig2) * invCl) %*% centered) / 2 -
      n * log(2 * pi) / 2 -
      determinant(sig2 * Cl, logarithm = TRUE)$modulus / 2,
    error = function(e) NA_real_
  )
  list(
    gls_mean = as.numeric(a),
    sigma2 = as.numeric(sig2),
    logL = if (!is.finite(logL)) -Inf else as.numeric(logL)
  )
}

.lambda_loglik_reference <- function(lambda, C, y) {
  .lambda_components_reference(lambda, C, y)$logL
}

.lambda_profile_data_batch <- function(max_lambda, lambda_hat, logL_hat,
                                       n_points, batch_lik = NULL,
                                       scalar_lik, Y) {
  grid <- seq(0, max_lambda, length.out = n_points)
  extra <- c(0, lambda_hat)
  if (max_lambda >= 1) extra <- c(extra, 1)
  lambda <- sort(unique(c(grid, extra)))
  logL <- if (!is.null(batch_lik)) {
    batch_lik(lambda)
  } else {
    vapply(lambda, function(z) {
      vapply(seq_len(ncol(Y)), function(j) scalar_lik(z, Y[, j]),
             numeric(1))
    }, numeric(ncol(Y)))
  }
  if (is.null(dim(logL))) logL <- matrix(logL, ncol = ncol(Y))
  lapply(seq_len(ncol(Y)), function(j) {
    profile_logL <- logL[, j]
    hit <- which.min(abs(lambda - lambda_hat[[j]]))
    if (length(hit) == 1L && is.finite(logL_hat[[j]])) {
      profile_logL[hit] <- logL_hat[[j]]
    }
    data.frame(lambda = lambda, logL = profile_logL,
               stringsAsFactors = FALSE)
  })
}

.lambda_profile_ci <- function(profile, lambda_hat, logL_hat, level = 0.95) {
  profile_max <- max(c(profile$logL, logL_hat), na.rm = TRUE)
  if (!is.finite(profile_max)) {
    profile_max <- logL_hat
  }
  cutoff <- profile_max - 0.5 * stats::qchisq(level, df = 1)
  lambda <- profile$lambda
  delta <- profile$logL - cutoff
  finite <- is.finite(lambda) & is.finite(delta)
  lambda <- lambda[finite]
  delta <- delta[finite]
  if (length(lambda) < 2L || !is.finite(lambda_hat)) {
    return(c(lower = NA_real_, upper = NA_real_, cutoff = cutoff))
  }
  left <- which(lambda <= lambda_hat)
  right <- which(lambda >= lambda_hat)
  lower <- .lambda_profile_bound(lambda[left], delta[left],
                                 side = "left")
  upper <- .lambda_profile_bound(lambda[right], delta[right],
                                 side = "right")
  c(lower = lower, upper = upper, cutoff = cutoff)
}

.lambda_profile_bound <- function(lambda, delta, side = c("left", "right")) {
  side <- match.arg(side)
  if (length(lambda) < 1L) return(NA_real_)
  if (length(lambda) == 1L) return(lambda[[1L]])

  if (side == "left") {
    if (!any(delta < 0) || !any(delta >= 0)) {
      return(min(lambda))
    }
    idx <- which(delta[-length(delta)] < 0 & delta[-1L] >= 0)
    if (!length(idx)) {
      return(min(lambda))
    }
    i <- idx[[length(idx)]]
  } else {
    if (!any(delta < 0) || !any(delta >= 0)) {
      return(max(lambda))
    }
    idx <- which(delta[-length(delta)] >= 0 & delta[-1L] < 0)
    if (!length(idx)) {
      return(max(lambda))
    }
    i <- idx[[1L]]
  }

  d1 <- delta[[i]]
  d2 <- delta[[i + 1L]]
  if (d1 == d2) {
    return(mean(lambda[c(i, i + 1L)]))
  }
  lambda[[i]] + (0 - d1) * (lambda[[i + 1L]] - lambda[[i]]) / (d2 - d1)
}


# Small utilities --------------------------------------------------------------

.as_trait_matrix <- function(X, tree, verbose = TRUE) {
  if (is.data.frame(X)) {
    if (!all(vapply(X, is.numeric, logical(1)))) {
      stop("All columns of x must be numeric.", call. = FALSE)
    }
    rn <- rownames(X)
    X <- data.matrix(X)
    rownames(X) <- rn
  } else if (is.null(dim(X))) {
    nm <- names(X)
    X <- matrix(as.numeric(X), ncol = 1L)
    rownames(X) <- nm
    colnames(X) <- "x"
  } else {
    if (!is.numeric(X)) {
      stop("x must be numeric.", call. = FALSE)
    }
    X <- as.matrix(X)
    storage.mode(X) <- "double"
  }

  default_rownames <- identical(rownames(X), as.character(seq_len(nrow(X))))
  if ((is.null(rownames(X)) || default_rownames) &&
      nrow(X) == ape::Ntip(tree)) {
    if (isTRUE(verbose)) {
      message("x has no names; assuming x is in tree$tip.label order")
    }
    rownames(X) <- tree$tip.label
    default_rownames <- FALSE
  }
  if (is.null(rownames(X)) || default_rownames) {
    stop("x must have species names as names/rownames.", call. = FALSE)
  }
  if (anyDuplicated(rownames(X))) {
    stop("x species names must be unique.", call. = FALSE)
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("trait_", seq_len(ncol(X)))
  }
  X
}

.permutation_matrix <- function(permutations, n, nsim,
                                include_observed = FALSE) {
  if (is.null(permutations)) {
    out <- matrix(NA_integer_, nrow = nsim, ncol = n)
    first_random <- 1L
    if (isTRUE(include_observed)) {
      out[1L, ] <- seq_len(n)
      first_random <- 2L
    }
    if (first_random <= nsim) {
      for (i in seq.int(first_random, nsim)) {
        out[i, ] <- sample.int(n)
      }
    }
    return(out)
  }
  permutations <- as.matrix(permutations)
  if (!is.numeric(permutations) || any(!is.finite(permutations)) ||
      any(permutations != floor(permutations))) {
    stop("permutations must contain finite integer indices.", call. = FALSE)
  }
  storage.mode(permutations) <- "integer"
  if (nrow(permutations) != nsim) {
    stop("permutations must have exactly nsim rows.", call. = FALSE)
  }
  if (ncol(permutations) != n) {
    stop("permutations must have one column per retained species.",
         call. = FALSE)
  }
  if (anyNA(permutations) || any(permutations < 1L) ||
      any(permutations > n)) {
    stop("permutations must contain 1-based indices in 1:n.",
         call. = FALSE)
  }
  duplicated_row <- vapply(
    seq_len(nrow(permutations)),
    function(i) anyDuplicated(permutations[i, ]) > 0L,
    logical(1)
  )
  if (any(duplicated_row)) {
    stop("each row of permutations must contain every index in 1:n once.",
         call. = FALSE)
  }
  permutations
}

.k_permutation_mcse <- function(exceedance, nsim,
                                include_observed = FALSE) {
  out <- rep(NA_real_, length(exceedance))
  ok <- is.finite(exceedance) & is.finite(nsim) & nsim > 0
  if (!any(ok)) return(out)

  if (isTRUE(include_observed)) {
    n_random <- nsim - 1
    if (n_random <= 0) return(out)
    p_random <- (exceedance[ok] - 1) / n_random
    p_random <- pmin(1, pmax(0, p_random))
    out[ok] <- sqrt(n_random * p_random * (1 - p_random)) / nsim
  } else {
    p <- exceedance[ok] / nsim
    p <- pmin(1, pmax(0, p))
    out[ok] <- sqrt(p * (1 - p) / nsim)
  }
  out
}

.normalize_ncores <- function(ncores) {
  if (!is.numeric(ncores) || length(ncores) != 1L || is.na(ncores) ||
      !is.finite(ncores) || ncores < 1 || ncores != floor(ncores) ||
      ncores > .Machine$integer.max) {
    stop("ncores must be a finite positive integer.", call. = FALSE)
  }
  detected <- parallel::detectCores(logical = TRUE)
  if (is.finite(detected) && detected >= 1) {
    ncores <- min(ncores, detected)
  }
  as.integer(floor(ncores))
}
