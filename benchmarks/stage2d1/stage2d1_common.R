# Shared helpers for the Stage 2D1 private K permutation audit.
#
# This file is an audit harness dependency only.  It never changes the
# package namespace, source tree, tests, or public API.

stage2d1_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  smoke <- any(args %in% c("--smoke", "quick"))
  formal <- "--formal" %in% args
  args <- args[!args %in% c("--smoke", "quick", "--formal")]
  repo <- if (length(args)) args[[1L]] else getwd()
  out <- if (length(args) >= 2L) args[[2L]] else
    file.path(repo, "benchmarks", "stage2d1", "results")
  list(
    smoke = smoke,
    formal = formal,
    repo = normalizePath(repo, winslash = "/", mustWork = FALSE),
    out = normalizePath(out, winslash = "/", mustWork = FALSE)
  )
}

stage2d1_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

stage2d1_ints <- function(value, name, minimum = 1L) {
  out <- suppressWarnings(as.integer(trimws(strsplit(
    as.character(value), ",", fixed = TRUE
  )[[1L]])))
  if (!length(out) || anyNA(out) || any(out < minimum)) {
    stop(name, " must contain comma-separated integers >= ", minimum, ".",
         call. = FALSE)
  }
  unique(out)
}

stage2d1_grid <- function(name, formal, smoke, minimum = 1L) {
  value <- stage2d1_env(name)
  if (nzchar(value)) stage2d1_ints(value, name, minimum) else
    if (isTRUE(smoke)) {
      smoke
    } else formal
}

stage2d1_clock <- function() unname(as.numeric(proc.time()[["elapsed"]]))

stage2d1_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "NA")
  invisible(path)
}

stage2d1_append_rows <- function(rows) {
  if (!length(rows)) return(data.frame(stringsAsFactors = FALSE))
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(cols, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[cols]
  })
  do.call(rbind, rows)
}

stage2d1_restore_seed <- function(seed, fun) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(fun)()
}

stage2d1_capture <- function(fun) {
  warnings <- character()
  status <- "ok"
  value <- NULL
  error_message <- ""
  error_class <- ""
  started <- stage2d1_clock()
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      status <<- "error"
      error_message <<- conditionMessage(e)
      error_class <<- paste(class(e), collapse = "/")
      NULL
    }
  )
  list(
    value = value,
    elapsed_s = max(0, stage2d1_clock() - started),
    status = status,
    warnings = unique(warnings),
    warning_text = paste(unique(warnings), collapse = " | "),
    error_message = error_message,
    error_class = error_class
  )
}

stage2d1_make_balanced <- function(n) {
  n <- as.integer(n)
  root <- n + 1L
  internal <- n - 1L
  edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  q_node <- q_lo <- q_hi <- integer(internal)
  q_node[[1L]] <- root
  q_lo[[1L]] <- 1L
  q_hi[[1L]] <- n
  head <- 1L
  tail <- 1L
  next_internal <- root + 1L
  edge_i <- 0L
  while (head <= tail) {
    node <- q_node[[head]]
    lo <- q_lo[[head]]
    hi <- q_hi[[head]]
    head <- head + 1L
    mid <- floor((lo + hi) / 2L)
    child_lo <- c(lo, mid + 1L)
    child_hi <- c(mid, hi)
    for (j in 1:2) {
      child <- if (child_lo[[j]] == child_hi[[j]]) child_lo[[j]] else {
        new_node <- next_internal
        next_internal <- next_internal + 1L
        tail <- tail + 1L
        q_node[[tail]] <- new_node
        q_lo[[tail]] <- child_lo[[j]]
        q_hi[[tail]] <- child_hi[[j]]
        new_node
      }
      edge_i <- edge_i + 1L
      edge[edge_i, ] <- c(node, child)
    }
  }
  edge <- edge[seq_len(edge_i), , drop = FALSE]
  tree <- list(
    edge = edge,
    tip.label = paste0("sp", seq_len(n)),
    edge.length = 0.5 + seq_len(nrow(edge)) / nrow(edge),
    Nnode = internal
  )
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

stage2d1_make_pectinate <- function(n) {
  n <- as.integer(n)
  root <- n + 1L
  if (n == 2L) {
    edge <- matrix(c(root, 1L, root, 2L), ncol = 2L, byrow = TRUE)
  } else {
    edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
    current <- root
    next_internal <- root + 1L
    edge_i <- 0L
    for (tip in seq_len(n - 2L)) {
      edge_i <- edge_i + 1L
      edge[edge_i, ] <- c(current, tip)
      edge_i <- edge_i + 1L
      edge[edge_i, ] <- c(current, next_internal)
      current <- next_internal
      next_internal <- next_internal + 1L
    }
    edge_i <- edge_i + 1L
    edge[edge_i, ] <- c(current, n - 1L)
    edge_i <- edge_i + 1L
    edge[edge_i, ] <- c(current, n)
    edge <- edge[seq_len(edge_i), , drop = FALSE]
  }
  tree <- list(
    edge = edge,
    tip.label = paste0("sp", seq_len(n)),
    edge.length = 0.5 + seq_len(nrow(edge)) / nrow(edge),
    Nnode = n - 1L
  )
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

stage2d1_make_random <- function(n) {
  stage2d1_restore_seed(20260908L + as.integer(n), function() {
    tree <- ape::reorder.phylo(ape::rtree(n), order = "postorder")
    tree$tip.label <- paste0("sp", seq_len(n))
    tree$edge.length <- 0.5 + seq_len(nrow(tree$edge)) / nrow(tree$edge)
    tree
  })
}

stage2d1_make_tree <- function(n, shape = "balanced") {
  switch(shape,
         balanced = stage2d1_make_balanced(n),
         pectinate = stage2d1_make_pectinate(n),
         random = stage2d1_make_random(n),
         stop("unknown tree shape: ", shape, call. = FALSE))
}

stage2d1_make_trait <- function(tree, offset = 0L) {
  n <- length(tree$tip.label)
  i <- seq_len(n) + as.integer(offset)
  stats::setNames(sin(i) + cos(i / 7) + i / (n + 11), tree$tip.label)
}

stage2d1_make_matrix <- function(tree, traits = 1L) {
  n <- length(tree$tip.label)
  i <- seq_len(n)
  j <- seq_len(traits)
  x <- outer(i, j, function(a, b) sin(a + b / 11) +
               cos(a / 7 + b) + (a + b) / (n + traits + 17))
  dimnames(x) <- list(tree$tip.label, paste0("trait", j))
  x
}

stage2d1_make_permutations <- function(n, nsim, seed = 1L) {
  stage2d1_restore_seed(seed, function() {
    out <- matrix(NA_integer_, nrow = nsim, ncol = n)
    out[1L, ] <- seq_len(n)
    if (nsim >= 2L) out[2L, ] <- rev(seq_len(n))
    if (nsim >= 3L) {
      for (i in 3:nsim) out[i, ] <- sample.int(n)
    }
    out
  })
}

stage2d1_make_na_matrix <- function(tree, kind = "none", traits = 4L) {
  x <- stage2d1_make_matrix(tree, traits)
  n <- nrow(x)
  if (kind == "shared") {
    x[seq.int(2L, n, by = max(2L, floor(n / 10L))), ] <- NA_real_
  } else if (kind == "repeated") {
    for (j in seq_len(ncol(x))) {
      x[seq.int(1L + (j %% 3L), n, by = max(3L, floor(n / 12L))), j] <- NA_real_
    }
  } else if (kind == "near_unique") {
    for (j in seq_len(ncol(x))) x[((j - 1L) %% n) + 1L, j] <- NA_real_
  } else if (kind != "none") {
    stop("unknown NA mask kind: ", kind, call. = FALSE)
  }
  x
}

stage2d1_get_private <- function(name, ns) {
  if (!exists(name, envir = ns, inherits = FALSE)) return(NULL)
  get(name, envir = ns, inherits = FALSE)
}

stage2d1_load_package <- function(repo, library_dir = stage2d1_env(
  "FASTPHYLOSIG_STAGE2D1_LIBRARY"
)) {
  if (nzchar(library_dir)) {
    library_dir <- normalizePath(library_dir, winslash = "/", mustWork = TRUE)
    .libPaths(c(library_dir, .libPaths()))
    suppressPackageStartupMessages(library("fastphylosig", lib.loc = library_dir,
                                           character.only = TRUE))
  } else {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("pkgload is required when no installed library is supplied.",
           call. = FALSE)
    }
    pkgload::load_all(repo, quiet = TRUE)
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("ape is required for the fixed fixtures.", call. = FALSE)
  }
  asNamespace("fastphylosig")
}

stage2d1_find_function <- function(name, envs) {
  for (env in envs) {
    if (!is.null(env) && exists(name, envir = env, inherits = FALSE)) {
      value <- get(name, envir = env, inherits = FALSE)
      if (is.function(value)) return(value)
    }
  }
  NULL
}

stage2d1_load_hooks <- function(repo) {
  hook_env <- new.env(parent = globalenv())
  path <- stage2d1_env("FASTPHYLOSIG_STAGE2D1_HOOK_FILE")
  cpp <- stage2d1_env("FASTPHYLOSIG_STAGE2D1_PROTO_CPP")
  if (!nzchar(cpp)) {
    sibling_cpp <- file.path(repo, "benchmarks", "stage2d1", "prototype",
                              "k_permutation_prototype.cpp")
    sibling_exists <- tryCatch(file.exists(sibling_cpp), error = function(e) FALSE)
    if (sibling_exists) cpp <- sibling_cpp
  }
  if (nzchar(path)) {
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    sys.source(path, envir = hook_env)
  }
  if (nzchar(cpp)) {
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      stop("Rcpp is required for FASTPHYLOSIG_STAGE2D1_PROTO_CPP.",
           call. = FALSE)
    }
    Rcpp::sourceCpp(normalizePath(cpp, winslash = "/", mustWork = TRUE),
                    env = hook_env)
  }
  envs <- list(hook_env, .GlobalEnv)
  if (exists("stage2d1_hooks", envir = hook_env, inherits = FALSE)) {
    value <- get("stage2d1_hooks", envir = hook_env, inherits = FALSE)
    if (is.list(value)) {
      envs <- c(list(value), envs)
    }
  }
  controlled_name <- stage2d1_env(
    "FASTPHYLOSIG_STAGE2D1_PROTO_CONTROLLED",
    "stage2d1_proto_controlled"
  )
  rng_name <- stage2d1_env(
    "FASTPHYLOSIG_STAGE2D1_PROTO_RNG",
    "stage2d1_proto_rng"
  )
  batch_name <- stage2d1_env(
    "FASTPHYLOSIG_STAGE2D1_PROTO_BATCH",
    "stage2d1_proto_batch"
  )
  controlled <- stage2d1_find_function(controlled_name, envs)
  rng <- stage2d1_find_function(rng_name, envs)
  # The checked-in private prototype exposes one mode-selectable entry point
  # rather than three R wrappers.  Auto-adapt that entry point so the same
  # harness can be used with it or with a separately named sourceCpp hook.
  prototype <- stage2d1_find_function("stage2d1_k_permutation_prototype", envs)
  prototype_mode <- stage2d1_env("FASTPHYLOSIG_STAGE2D1_PROTO_MODE", "A")
  if (!is.null(prototype) && (is.null(controlled) || is.null(rng))) {
    make_wrapper <- function(target, selected_mode) {
      force(target)
      force(selected_mode)
      function(...) {
        args <- list(...)
        args$return_ordered <- identical(
          tolower(stage2d1_env("FASTPHYLOSIG_STAGE2D1_RETURN_ORDERED", "true")),
          "true"
        )
        args$mode <- selected_mode
        if ("collect_counters" %in% names(formals(target))) {
          args$collect_counters <- identical(
            tolower(stage2d1_env("FASTPHYLOSIG_STAGE2D1_COLLECT_COUNTERS", "false")),
            "true"
          )
        }
        do.call(target, args)
      }
    }
    if (is.null(controlled)) controlled <- make_wrapper(prototype, prototype_mode)
    if (is.null(rng)) rng <- make_wrapper(prototype, prototype_mode)
  }
  list(
    env = hook_env,
    controlled = controlled,
    rng = rng,
    batch = stage2d1_find_function(batch_name, envs),
    prototype = prototype,
    prototype_mode = prototype_mode,
    controlled_name = controlled_name,
    rng_name = rng_name,
    batch_name = batch_name,
    source_file = path,
    source_cpp = cpp
  )
}

stage2d1_production_cpp <- function(ns) {
  fun <- stage2d1_get_private("fast_k_tree_permutation_cpp", ns)
  if (!is.function(fun)) {
    stop("fast_k_tree_permutation_cpp is unavailable in package namespace.",
         call. = FALSE)
  }
  fun
}

stage2d1_compile_tree <- function(tree, ns) {
  compiler <- stage2d1_get_private("compile_tree_cpp", ns)
  if (!is.function(compiler)) {
    stop("compile_tree_cpp is unavailable in package namespace.", call. = FALSE)
  }
  compiler(tree$edge, tree$edge.length, length(tree$tip.label))
}

stage2d1_engine_args <- function(compiled_tree, X, nsim, permutations,
                                 trait_chunk, return_sim, include_observed,
                                 ncores, simulation_chunk) {
  list(
    compiled_tree = compiled_tree,
    X = X,
    nsim = as.integer(nsim),
    permutations = permutations,
    trait_chunk = as.integer(trait_chunk),
    return_sim = isTRUE(return_sim),
    include_observed = isTRUE(include_observed),
    n_threads = as.integer(ncores),
    simulation_chunk = as.integer(simulation_chunk)
  )
}

stage2d1_call_with_seed <- function(fun, args, seed = NULL) {
  if (is.null(seed)) return(do.call(fun, args))
  fml <- tryCatch(names(formals(fun)), error = function(e) character())
  seed_arg <- stage2d1_env("FASTPHYLOSIG_STAGE2D1_SEED_ARG", "seed")
  if (seed_arg %in% fml) {
    return(do.call(fun, c(args, stats::setNames(list(as.integer(seed)), seed_arg))))
  }
  stage2d1_restore_seed(seed, function() do.call(fun, args))
}

stage2d1_call_engine <- function(fun, mode, compiled_tree, X, nsim,
                                  permutations = NULL, trait_chunk = 64L,
                                  return_sim = TRUE, include_observed = TRUE,
                                  ncores = 1L, simulation_chunk = 128L,
                                  seed = NULL) {
  args <- stage2d1_engine_args(
    compiled_tree = compiled_tree, X = X, nsim = nsim,
    permutations = permutations, trait_chunk = trait_chunk,
    return_sim = return_sim, include_observed = include_observed,
    ncores = ncores, simulation_chunk = simulation_chunk
  )
  if (identical(mode, "controlled")) {
    return(do.call(fun, args))
  }
  if (!identical(mode, "rng")) stop("unknown engine mode: ", mode, call. = FALSE)
  stage2d1_call_with_seed(fun, args, seed)
}

stage2d1_pick <- function(x, names, default = NULL) {
  if (is.null(x) || !is.list(x)) return(default)
  for (nm in names) if (!is.null(x[[nm]])) return(x[[nm]])
  default
}

stage2d1_as_sim_matrix <- function(value, nsim, p) {
  if (is.null(value)) return(NULL)
  value <- as.matrix(value)
  if (identical(dim(value), c(as.integer(nsim), as.integer(p)))) return(value)
  if (identical(dim(value), c(as.integer(p), as.integer(nsim)))) return(t(value))
  as.numeric(value)
}

stage2d1_normalize_result <- function(value, nsim, p) {
  if (is.null(value)) {
    return(list(K = NULL, P = NULL, MCSE_P = NULL, exceedance = NULL,
                requested = NULL, successful = NULL, failed = NULL,
                success = NULL, sim_K = NULL, names = character()))
  }
  if (is.atomic(value) && is.null(dim(value))) {
    value <- list(K = value)
  }
  if (!is.list(value)) value <- list(K = value)
  sim <- stage2d1_pick(value, c(
    "sim_K", "sim.K", "sim_K_fast", "ordered_null_K", "ordered_null",
    "simulation", "null_K", "null_k"
  ))
  success <- stage2d1_pick(value, c(
    "success", "successful", "valid", "simulation_success", "replicate_success"
  ))
  if (is.null(success) && !is.null(sim)) success <- rep(TRUE, nsim)
  requested <- stage2d1_pick(value, c(
    "nsim_requested", "requested", "n_sim_requested", "nPerm", "nsim"
  ), nsim)
  successful <- stage2d1_pick(value, c(
    "nsim_successful", "successful_simulations", "n_sim_successful",
    "successful_count"
  ))
  failed <- stage2d1_pick(value, c(
    "nsim_failed", "failed_simulations", "n_sim_failed", "failed_count"
  ))
  if (is.null(successful) && length(requested) == 1L && !is.null(failed)) {
    successful <- as.numeric(requested) - as.numeric(failed)
  }
  if (is.null(failed) && length(requested) == 1L && !is.null(successful)) {
    failed <- as.numeric(requested) - as.numeric(successful)
  }
  list(
    K = stage2d1_pick(value, c("K", "K_fast", "observed", "observed_K")),
    P = stage2d1_pick(value, c("P", "P_fast", "p_value", "p")),
    MCSE_P = stage2d1_pick(value, c("MCSE_P", "P_MCSE", "p_mcse", "mcse")),
    exceedance = stage2d1_pick(value, c(
      "exceedance_count", "exceedance", "n_exceed", "tail_count"
    )),
    requested = requested,
    successful = successful,
    failed = failed,
    success = success,
    sim_K = stage2d1_as_sim_matrix(sim, nsim, p),
    names = names(value)
  )
}

stage2d1_num_equal <- function(x, y, tolerance = 0) {
  if (is.null(x) || is.null(y)) return(is.null(x) && is.null(y))
  if (identical(x, y)) return(TRUE)
  isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = TRUE))
}

stage2d1_max_abs_diff <- function(x, y) {
  if (is.null(x) || is.null(y)) return(NA_real_)
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) return(Inf)
  if (!length(x)) return(0)
  z <- abs(x - y)
  if (all(is.na(z))) 0 else max(z, na.rm = TRUE)
}

stage2d1_compare_results <- function(left, right, nsim, p,
                                     numeric_tolerance = 0) {
  a <- stage2d1_normalize_result(left, nsim, p)
  b <- stage2d1_normalize_result(right, nsim, p)
  fields <- c("K", "P", "MCSE_P", "exceedance", "requested", "successful",
              "failed", "success", "sim_K")
  pass <- vapply(fields, function(nm) {
    stage2d1_num_equal(a[[nm]], b[[nm]], numeric_tolerance)
  }, logical(1))
  list(
    pass = all(pass), fields = pass,
    max_abs = c(
      K = stage2d1_max_abs_diff(a$K, b$K),
      P = stage2d1_max_abs_diff(a$P, b$P),
      MCSE_P = stage2d1_max_abs_diff(a$MCSE_P, b$MCSE_P),
      exceedance = stage2d1_max_abs_diff(a$exceedance, b$exceedance),
      sim_K = stage2d1_max_abs_diff(a$sim_K, b$sim_K)
    ),
    left = a,
    right = b
  )
}

stage2d1_compare_captures <- function(left, right, nsim, p,
                                      numeric_tolerance = 0) {
  detail <- if (left$status == "ok" && right$status == "ok") {
    stage2d1_compare_results(left$value, right$value, nsim, p,
                             numeric_tolerance)
  } else list(pass = identical(left$status, right$status), fields = logical(),
              max_abs = numeric())
  warning_pass <- identical(left$warning_text, right$warning_text)
  status_pass <- identical(left$status, right$status)
  list(
    pass = isTRUE(detail$pass) && status_pass && warning_pass,
    result = detail,
    status_pass = status_pass,
    warning_pass = warning_pass,
    left_status = left$status,
    right_status = right$status,
    left_warning = left$warning_text,
    right_warning = right$warning_text,
    left_error = left$error_message,
    right_error = right$error_message
  )
}

stage2d1_call_description <- function(mode, candidate, n, nsim, shape,
                                      seed, ncores, chunk, traits = 1L) {
  paste(mode, candidate, shape, paste0("n=", n), paste0("nsim=", nsim),
        paste0("seed=", seed), paste0("ncores=", ncores),
        paste0("chunk=", chunk), paste0("traits=", traits), sep = "/")
}

stage2d1_workspace_proxy <- function(n, nsim, p = 1L, trait_chunk = 64L,
                                     simulation_chunk = 128L,
                                     return_sim = TRUE, candidate = FALSE,
                                     reported = NA_real_) {
  n <- as.double(n)
  nsim <- as.double(nsim)
  p <- as.double(p)
  chunk <- min(nsim, as.double(simulation_chunk))
  base <- 4 * n
  out <- 8 * nsim * p
  trait <- 8 * (2 * (n + max(1, n - 1)) * min(p, trait_chunk) +
                  2 * min(p, trait_chunk))
  index_chunk <- if (isTRUE(candidate)) base else base * chunk
  data.frame(
    candidate = isTRUE(candidate), n = n, nsim = nsim, traits = p,
    per_replicate_index_bytes = base,
    active_index_chunk_bytes = index_chunk,
    trait_kernel_workspace_bytes = trait,
    requested_sim_output_bytes = if (isTRUE(return_sim)) out else 0,
    controlled_input_matrix_bytes = if (isTRUE(candidate)) NA_real_ else 4 * n * nsim,
    reported_workspace_bytes = reported,
    stringsAsFactors = FALSE
  )
}

stage2d1_frozen_contract <- function(repo) {
  lines <- c(
    "# K RNG Contract: Stage 2D1",
    "",
    "Status: extracted from the current source and user-facing documentation;",
    "this document adds no new public RNG promise.",
    "",
    "## Contract used by this audit",
    "",
    "- A replay comparison holds the seed, arguments, retained species order,",
    "  simulation chunk, and `ncores` fixed. Changing thread count is not used",
    "  as evidence of a new worker-count-invariant stochastic contract.",
    "- A supplied permutation matrix is the controlled-input path. It has one",
    "  1-based permutation row per requested replicate and one column per",
    "  retained species; each row must contain 1:n exactly once.",
    "- Internal RNG mode uses the existing C++ Fisher-Yates implementation",
    "  backed by R's RNG (`R::runif`). The first row is identity only when the",
    "  existing `include_observed` behavior is enabled; subsequent rows are",
    "  random permutations.",
    "- The current source comments state that internal-RNG permutation draws",
    "  are generated serially before OpenMP evaluation and are independent of",
    "  thread count and `simulation_chunk`. The harness records this as an",
    "  observed implementation property; it does not broaden the public API",
    "  contract beyond the existing tests and documentation.",
    "- README reproducibility guidance requires the seed, explicit stochastic",
    "  inputs when supplied, and `ncores` to be held fixed. Delta's process",
    "  count is explicitly not draw-for-draw invariant; that statement is not",
    "  transferred to K.",
    "- `P` and `MCSE_P` are compared from the existing inclusive upper-tail",
    "  accounting. Requested, successful, and failed counts are compared",
    "  without changing the estimator or its tolerance.",
    "",
    "## Source anchors",
    "",
    "- `src/k_permutation.cpp`: `fisher_yates`, the controlled branch, and the",
    "  serial RNG generation before chunked OpenMP evaluation.",
    "- `R/fast_signal.R`: `.permutation_matrix`, `simulation_chunk`, and the",
    "  public K arguments.",
    "- `README.md`: Reproducibility and performance section, especially the",
    "  fixed seed/inputs/`ncores` guidance.",
    "",
    "## Non-contract claims",
    "",
    "The audit does not promise identical stochastic trajectories after changing",
    "`ncores`, changing RNG kind, changing compiler/platform, or changing the",
    "retained-tip order. Such observations are reported separately.",
    ""
  )
  path <- file.path(repo, "benchmarks", "stage2d1", "K_RNG_CONTRACT_STAGE2D1.md")
  # A sibling audit may already own this contract file.  Do not overwrite an
  # existing evidence artifact merely because this harness was sourced.
  if (!file.exists(path)) writeLines(lines, path, useBytes = TRUE)
  path
}

stage2d1_provenance <- function(repo, stage, args, hooks) {
  cores <- tryCatch(parallel::detectCores(logical = TRUE),
                    error = function(e) NA_integer_)
  source_commit <- stage2d1_env("FASTPHYLOSIG_STAGE2D1_SOURCE_COMMIT")
  if (!nzchar(source_commit)) {
    source_commit <- tryCatch(trimws(system2(
      "git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE,
      stderr = TRUE
    )), error = function(e) "<unavailable>")
  }
  data.frame(
    stage = stage,
    mode = if (isTRUE(args$smoke)) "smoke" else "formal",
    source_commit = source_commit,
    package_version = as.character(utils::packageVersion("fastphylosig")),
    R_version = R.version.string,
    platform = R.version$platform,
    os = Sys.info()[["sysname"]],
    machine = Sys.info()[["machine"]],
    logical_cores = cores,
    BLAS = tryCatch(extSoftVersion()[["BLAS"]], error = function(e) NA_character_),
    LAPACK = tryCatch(extSoftVersion()[["LAPACK"]], error = function(e) NA_character_),
    locale_collate = Sys.getlocale("LC_COLLATE"),
    prototype_hook_file = if (nzchar(hooks$source_file)) hooks$source_file else "",
    prototype_cpp = if (nzchar(hooks$source_cpp)) hooks$source_cpp else "",
    controlled_hook = hooks$controlled_name,
    rng_hook = hooks$rng_name,
    batch_hook = hooks$batch_name,
    prototype_mode = if (!is.null(hooks$prototype_mode)) hooks$prototype_mode else "",
    production_code_changed = "NO",
    stringsAsFactors = FALSE
  )
}
