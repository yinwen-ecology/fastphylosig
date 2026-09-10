# Shared helpers for the Stage 2D2 exactness audit.
#
# This file is audit infrastructure only.  It never changes the package
# namespace, production source, tests, public API, or Stage 2D1 evidence.
# The private hook contract is intentionally explicit: an oracle and the
# Candidate A/B evaluators must be supplied by one sourceCpp translation unit
# (or by an R hook that declares the same-translation-unit property).

stage2d2_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

stage2d2_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  smoke <- "--smoke" %in% args
  correctness <- "--correctness" %in% args
  if (!smoke && !correctness) {
    stop("Stage 2D2 correctness is opt-in: use --smoke or --correctness.",
         call. = FALSE)
  }
  args <- args[!args %in% c("--smoke", "--correctness")]
  if (length(args) > 2L) {
    stop("usage: run_stage2d2_correctness.R [--smoke|--correctness] [repo] [out].",
         call. = FALSE)
  }
  repo <- if (length(args)) args[[1L]] else getwd()
  out <- if (length(args) >= 2L) args[[2L]] else
    file.path(repo, "benchmarks", "stage2d2", "results", "correctness")
  list(
    level = if (isTRUE(correctness)) "targeted" else "smoke",
    smoke = isTRUE(smoke),
    correctness = isTRUE(correctness),
    repo = normalizePath(repo, winslash = "/", mustWork = FALSE),
    out = normalizePath(out, winslash = "/", mustWork = FALSE)
  )
}

stage2d2_ints <- function(value, name, minimum = 1L) {
  out <- suppressWarnings(as.integer(trimws(strsplit(
    as.character(value), ",", fixed = TRUE
  )[[1L]])))
  if (!length(out) || anyNA(out) || any(out < minimum)) {
    stop(name, " must contain comma-separated integers >= ", minimum, ".",
         call. = FALSE)
  }
  unique(out)
}

stage2d2_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "NA")
  invisible(path)
}

stage2d2_rbind <- function(rows) {
  if (!length(rows)) return(data.frame(stringsAsFactors = FALSE))
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(cols, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[cols]
  })
  do.call(rbind, rows)
}

stage2d2_clock <- function() unname(as.numeric(proc.time()[["elapsed"]]))

stage2d2_capture <- function(fun) {
  warnings <- list()
  status <- "ok"
  value <- NULL
  error_message <- ""
  error_class <- character()
  started <- stage2d2_clock()
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- list(
          text = conditionMessage(w), classes = class(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      status <<- "error"
      error_message <<- conditionMessage(e)
      error_class <<- class(e)
      NULL
    }
  )
  warning_text <- vapply(warnings, `[[`, character(1), "text")
  warning_classes <- unique(unlist(lapply(warnings, `[[`, "classes"),
                                 use.names = FALSE))
  list(
    value = value,
    elapsed_s = max(0, stage2d2_clock() - started),
    status = status,
    warnings = warnings,
    warning_text = paste(unique(warning_text), collapse = " | "),
    warning_classes = paste(warning_classes, collapse = "/"),
    error_message = error_message,
    error_class = paste(error_class, collapse = "/")
  )
}

stage2d2_restore_seed <- function(seed, fun) {
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

stage2d2_get_private <- function(name, ns) {
  if (!exists(name, envir = ns, inherits = FALSE)) return(NULL)
  value <- get(name, envir = ns, inherits = FALSE)
  if (is.function(value)) value else NULL
}

stage2d2_load_package <- function(repo, library_dir = stage2d2_env(
  "FASTPHYLOSIG_STAGE2D2_LIBRARY"
)) {
  if (nzchar(library_dir)) {
    library_dir <- normalizePath(library_dir, winslash = "/", mustWork = TRUE)
    .libPaths(c(library_dir, .libPaths()))
    suppressPackageStartupMessages(library(
      "fastphylosig", lib.loc = library_dir, character.only = TRUE
    ))
  } else {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("pkgload is required when FASTPHYLOSIG_STAGE2D2_LIBRARY is unset.",
           call. = FALSE)
    }
    pkgload::load_all(repo, quiet = TRUE)
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("ape is required for Stage 2D2 fixed tree fixtures.", call. = FALSE)
  }
  asNamespace("fastphylosig")
}

stage2d2_compile_tree <- function(tree, ns) {
  compiler <- stage2d2_get_private("compile_tree_cpp", ns)
  if (!is.function(compiler)) {
    stop("compile_tree_cpp is unavailable in fastphylosig namespace.",
         call. = FALSE)
  }
  compiler(tree$edge, tree$edge.length, length(tree$tip.label))
}

# A deterministic set of valid permutations.  The first rows are fixed and
# later rows use a separately restored R stream; controlled audits consume no
# RNG inside either evaluator.
stage2d2_make_permutations <- function(n, nsim, seed = 1L) {
  stage2d2_restore_seed(seed, function() {
    out <- matrix(NA_integer_, nrow = nsim, ncol = n)
    out[1L, ] <- seq_len(n)
    if (nsim >= 2L) out[2L, ] <- rev(seq_len(n))
    if (nsim >= 3L) {
      for (i in 3:nsim) out[i, ] <- sample.int(n)
    }
    out
  })
}

stage2d2_make_balanced <- function(n) {
  n <- as.integer(n)
  if (n < 2L) stop("n must be at least two.", call. = FALSE)
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
  tree <- list(
    edge = edge[seq_len(edge_i), , drop = FALSE],
    tip.label = paste0("sp", seq_len(n)),
    edge.length = 0.5 + seq_len(edge_i) / edge_i,
    Nnode = internal
  )
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

stage2d2_make_pectinate <- function(n) {
  n <- as.integer(n)
  if (n < 2L) stop("n must be at least two.", call. = FALSE)
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

stage2d2_make_random <- function(n) {
  stage2d2_restore_seed(20260908L + as.integer(n), function() {
    tree <- ape::reorder.phylo(ape::rtree(n), order = "postorder")
    tree$tip.label <- paste0("sp", seq_len(n))
    tree$edge.length <- 0.5 + seq_len(nrow(tree$edge)) / nrow(tree$edge)
    tree
  })
}

stage2d2_make_tree <- function(n, shape = "balanced") {
  switch(shape,
         balanced = stage2d2_make_balanced(n),
         random = stage2d2_make_random(n),
         pectinate = stage2d2_make_pectinate(n),
         stop("unknown tree shape: ", shape, call. = FALSE))
}

stage2d2_make_matrix <- function(tree, traits = 1L, kind = "normal",
                                  seed = 20260909L) {
  n <- length(tree$tip.label)
  p <- as.integer(traits)
  if (p < 1L) stop("traits must be positive.", call. = FALSE)
  z <- stage2d2_restore_seed(seed + p, function() {
    matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  })
  if (kind == "large_offset") {
    z <- z + 1e12
  } else if (kind == "near_constant") {
    z <- 3 + z * 1e-8
  } else if (kind == "constant") {
    z[,] <- 3
  } else if (kind != "normal") {
    stop("unknown trait kind: ", kind, call. = FALSE)
  }
  dimnames(z) <- list(tree$tip.label, paste0("trait", seq_len(p)))
  z
}

stage2d2_masked_cases <- function(tree, traits = 8L, kind = "shared") {
  X <- stage2d2_make_matrix(tree, traits = traits, kind = "normal",
                            seed = 20261001L)
  n <- nrow(X)
  p <- ncol(X)
  missing <- matrix(FALSE, nrow = n, ncol = p)
  if (kind == "shared") {
    missing[seq.int(2L, n, by = max(2L, floor(n / 4L))), ] <- TRUE
  } else if (kind == "multiple") {
    for (j in seq_len(p)) {
      start <- 1L + ((j - 1L) %% max(1L, min(3L, n - 1L)))
      missing[seq.int(start, n, by = max(3L, floor(n / 5L))), j] <- TRUE
    }
  } else {
    stop("unknown mask kind: ", kind, call. = FALSE)
  }
  X[missing] <- NA_real_
  pattern <- apply(is.na(X), 2L, paste, collapse = ",")
  groups <- split(seq_len(p), pattern)
  out <- list()
  for (cols in groups) {
    keep <- !is.na(X[, cols[[1L]]])
    if (sum(keep) < 2L) next
    retained <- tree$tip.label[keep]
    sub_tree <- ape::keep.tip(tree, retained)
    sub_tree <- ape::reorder.phylo(sub_tree, order = "postorder")
    sub_X <- X[keep, cols, drop = FALSE]
    sub_X <- sub_X[sub_tree$tip.label, , drop = FALSE]
    out[[length(out) + 1L]] <- list(
      tree = sub_tree, X = sub_X, trait_columns = cols,
      mask_kind = kind, retained = sub_tree$tip.label
    )
  }
  out
}

stage2d2_find_function <- function(name, envs) {
  for (env in envs) {
    if (is.null(env)) next
    if (is.environment(env) && exists(name, envir = env, inherits = FALSE)) {
      value <- get(name, envir = env, inherits = FALSE)
    } else if (is.list(env) && !is.null(env[[name]])) {
      value <- env[[name]]
    } else {
      next
    }
    if (is.function(value)) return(value)
  }
  NULL
}

stage2d2_make_mode_wrapper <- function(target, mode) {
  force(target)
  force(mode)
  function(...) {
    args <- list(...)
    args$mode <- mode
    do.call(target, stage2d2_filter_call_args(target, args))
  }
}

stage2d2_load_hooks <- function(repo) {
  hook_env <- new.env(parent = globalenv())
  hook_file <- stage2d2_env("FASTPHYLOSIG_STAGE2D2_HOOK_FILE")
  cpp_file <- stage2d2_env("FASTPHYLOSIG_STAGE2D2_PROTO_CPP")
  if (!nzchar(hook_file) && !nzchar(cpp_file)) {
    return(list(ready = FALSE, reason = paste0(
      "No opt-in Stage 2D2 hook configured. Set ",
      "FASTPHYLOSIG_STAGE2D2_HOOK_FILE or FASTPHYLOSIG_STAGE2D2_PROTO_CPP."
    )))
  }
  if (nzchar(hook_file)) {
    hook_file <- normalizePath(hook_file, winslash = "/", mustWork = TRUE)
    sys.source(hook_file, envir = hook_env)
  }
  names_before_cpp <- ls(envir = hook_env, all.names = TRUE)
  marker_false_before_cpp <- any(vapply(
    c("stage2d2_oracle_same_translation_unit",
      "oracle_same_translation_unit"),
    function(nm) exists(nm, envir = hook_env, inherits = FALSE) &&
      identical(get(nm, envir = hook_env, inherits = FALSE), FALSE),
    logical(1)
  ))
  cpp_export_names <- character()
  if (nzchar(cpp_file)) {
    cpp_file <- normalizePath(cpp_file, winslash = "/", mustWork = TRUE)
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      stop("Rcpp is required for FASTPHYLOSIG_STAGE2D2_PROTO_CPP.",
           call. = FALSE)
    }
    Rcpp::sourceCpp(cpp_file, env = hook_env, rebuild = TRUE,
                    showOutput = FALSE, verbose = FALSE)
    cpp_export_names <- setdiff(ls(envir = hook_env, all.names = TRUE),
                                names_before_cpp)
  }
  envs <- list(hook_env, .GlobalEnv)
  hook_list <- if (exists("stage2d2_hooks", envir = hook_env, inherits = FALSE))
    get("stage2d2_hooks", envir = hook_env, inherits = FALSE) else NULL
  if (is.list(hook_list)) envs <- c(list(hook_list), envs)

  oracle <- stage2d2_find_function("stage2d2_k_permutation_oracle", envs)
  if (is.null(oracle)) oracle <- stage2d2_find_function(
    "stage2d2_k_null_oracle", envs
  )
  if (is.null(oracle)) oracle <- stage2d2_find_function("oracle", envs)
  candidate_a <- stage2d2_find_function(
    "stage2d2_k_permutation_candidate_a", envs
  )
  if (is.null(candidate_a)) candidate_a <- stage2d2_find_function(
    "stage2d2_k_null_candidate_a", envs
  )
  if (is.null(candidate_a)) candidate_a <- stage2d2_find_function(
    "candidate_a", envs
  )
  candidate_b <- stage2d2_find_function(
    "stage2d2_k_permutation_candidate_b", envs
  )
  if (is.null(candidate_b)) candidate_b <- stage2d2_find_function(
    "stage2d2_k_null_candidate_b", envs
  )
  if (is.null(candidate_b)) candidate_b <- stage2d2_find_function(
    "candidate_b", envs
  )
  candidate_c <- stage2d2_find_function(
    "stage2d2_k_permutation_candidate_c", envs
  )
  if (is.null(candidate_c)) candidate_c <- stage2d2_find_function(
    "stage2d2_k_null_candidate_c", envs
  )
  generic <- stage2d2_find_function(
    "stage2d2_k_permutation_prototype", envs
  )
  if (is.null(generic)) generic <- stage2d2_find_function(
    "stage2d2_k_null_prototype", envs
  )
  oracle_from_generic <- FALSE
  candidate_a_from_generic <- FALSE
  candidate_b_from_generic <- FALSE
  candidate_c_from_generic <- FALSE
  if (!is.null(generic)) {
    if (is.null(oracle)) {
      oracle <- stage2d2_make_mode_wrapper(generic, "oracle")
      oracle_from_generic <- TRUE
    }
    if (is.null(candidate_a)) {
      candidate_a <- stage2d2_make_mode_wrapper(generic, "A")
      candidate_a_from_generic <- TRUE
    }
    if (is.null(candidate_b)) {
      candidate_b <- stage2d2_make_mode_wrapper(generic, "B")
      candidate_b_from_generic <- TRUE
    }
    if (is.null(candidate_c)) {
      candidate_c <- stage2d2_make_mode_wrapper(generic, "C")
      candidate_c_from_generic <- TRUE
    }
  }

  b_gate <- toupper(trimws(stage2d2_env(
    "FASTPHYLOSIG_STAGE2D2_CANDIDATE_B_GATE", "NOT_AUTHORIZED"
  )))
  b_authorized <- b_gate %in% c("GO", "AUTHORIZED", "RUN", "PASS")
  candidate_mode <- toupper(trimws(stage2d2_env(
    "FASTPHYLOSIG_STAGE2D2_CANDIDATE_MODE", "A"
  )))
  if (!candidate_mode %in% c("A", "B", "C")) {
    stop("FASTPHYLOSIG_STAGE2D2_CANDIDATE_MODE must be A, B, or C.", call. = FALSE)
  }
  if (candidate_mode == "B" && !b_authorized) {
    stop("Candidate B-only audit requires an authorized Candidate B gate.",
         call. = FALSE)
  }
  b_status <- if (is.null(candidate_b)) "NOT_AUTHORIZED" else
    if (b_authorized) "AUTHORIZED" else "NOT_AUTHORIZED"

  # A function is sourceCpp-bound only when the actual function object is one
  # of the bindings created by sourceCpp in this exact environment.  Merely
  # setting an environment variable or placing an R wrapper in hook_env is
  # not evidence of a same-translation-unit oracle.
  cpp_bound <- function(fun) {
    if (!nzchar(cpp_file) || is.null(fun) || !length(cpp_export_names)) {
      return(FALSE)
    }
    any(vapply(cpp_export_names, function(nm) {
      exists(nm, envir = hook_env, inherits = FALSE) &&
        identical(get(nm, envir = hook_env, inherits = FALSE), fun)
    }, logical(1)))
  }
  generic_bound <- cpp_bound(generic)
  oracle_bound <- cpp_bound(oracle) || (oracle_from_generic && generic_bound)
  candidate_a_bound <- cpp_bound(candidate_a) ||
    (candidate_a_from_generic && generic_bound)
  candidate_b_bound <- cpp_bound(candidate_b) ||
    (candidate_b_from_generic && generic_bound)
  candidate_c_bound <- cpp_bound(candidate_c) ||
    (candidate_c_from_generic && generic_bound)
  marker_false_after_cpp <- any(vapply(
    c("stage2d2_oracle_same_translation_unit",
      "oracle_same_translation_unit"),
    function(nm) exists(nm, envir = hook_env, inherits = FALSE) &&
      identical(get(nm, envir = hook_env, inherits = FALSE), FALSE),
    logical(1)
  ))
  required_candidate_bound <- switch(
    candidate_mode, A = candidate_a_bound, B = candidate_b_bound,
    C = candidate_c_bound
  )
  declared <- oracle_bound && required_candidate_bound &&
    (candidate_mode == "B" || !b_authorized || candidate_b_bound) &&
    !marker_false_before_cpp && !marker_false_after_cpp
  required_candidate <- switch(
    candidate_mode, A = candidate_a, B = candidate_b, C = candidate_c
  )
  ready <- is.function(oracle) && is.function(required_candidate) &&
    isTRUE(declared) && (candidate_mode == "B" || !b_authorized ||
      is.function(candidate_b))
  reason <- if (ready) "" else paste0(
    "Stage 2D2 requires a sourceCpp-bound oracle and Candidate A; Candidate B ",
    "is required only when its phase gate is authorized. ",
    "same-translation-unit oracle. Found oracle=", !is.null(oracle),
    ", candidate_a=", !is.null(candidate_a), ", candidate_b=",
    !is.null(candidate_b), ", B_gate=", b_gate,
    ", same_translation_unit=", declared, "."
  )
  list(
    ready = ready, reason = reason, env = hook_env,
    oracle = oracle, candidate_a = candidate_a, candidate_b = candidate_b,
    candidate_c = candidate_c,
    hook_file = hook_file, cpp_file = cpp_file,
    source_cpp_env = hook_env, source_cpp_exports = cpp_export_names,
    same_translation_unit = declared, candidate_b_gate = b_gate,
    candidate_b_authorized = b_authorized, candidate_b_status = b_status,
    candidate_mode = candidate_mode
  )
}

stage2d2_filter_call_args <- function(fun, args) {
  fml <- tryCatch(names(formals(fun)), error = function(e) character())
  if (!length(fml) || "..." %in% fml) return(args)
  args[intersect(names(args), fml)]
}

stage2d2_call_with_seed <- function(fun, args, seed = NULL) {
  invoke <- function() do.call(fun, stage2d2_filter_call_args(fun, args))
  if (is.null(seed)) return(invoke())
  stage2d2_restore_seed(seed, invoke)
}

stage2d2_call_engine <- function(fun, compiled_tree, X, nsim,
                                 permutations = NULL, trait_chunk = 64L,
                                 return_sim = TRUE, include_observed = TRUE,
                                 n_threads = 1L, simulation_chunk = 128L,
                                 seed = NULL) {
  fml <- tryCatch(names(formals(fun)), error = function(e) character())
  args <- list(
    compiled_tree = compiled_tree,
    X = X,
    nsim = as.integer(nsim),
    permutations = permutations,
    trait_chunk = as.integer(trait_chunk),
    return_sim = isTRUE(return_sim),
    include_observed = isTRUE(include_observed),
    n_threads = as.integer(n_threads),
    simulation_chunk = as.integer(simulation_chunk),
    return_ordered = TRUE
  )
  if (!("n_threads" %in% fml) && "ncores" %in% fml) {
    args$ncores <- args$n_threads
    args$n_threads <- NULL
  }
  stage2d2_call_with_seed(fun, args, seed)
}

stage2d2_pick <- function(x, names, default = NULL) {
  if (is.null(x) || !is.list(x)) return(default)
  for (nm in names) if (!is.null(x[[nm]])) return(x[[nm]])
  default
}

stage2d2_as_sim_matrix <- function(value, nsim, p) {
  if (is.null(value)) return(NULL)
  if (is.matrix(value) && identical(dim(value), c(as.integer(nsim),
                                                   as.integer(p)))) {
    return(value)
  }
  if (is.matrix(value) && identical(dim(value), c(as.integer(p),
                                                   as.integer(nsim)))) return(value)
  if (length(value) != nsim * p) return(value)
  dim(value) <- c(as.integer(nsim), as.integer(p))
  value
}

stage2d2_normalize_result <- function(value, nsim, p) {
  if (is.null(value)) return(NULL)
  if (is.atomic(value) && is.null(dim(value))) value <- list(K = value)
  if (!is.list(value)) value <- list(K = value)
  sim <- stage2d2_pick(value, c(
    "sim_K", "sim.K", "ordered_null_K", "ordered_null",
    "simulation", "null_K", "null_k"
  ))
  success <- stage2d2_pick(value, c(
    "replicate_status", "status_by_replicate", "success", "valid",
    "simulation_success", "replicate_success"
  ))
  requested <- stage2d2_pick(value, c(
    "nsim_requested", "requested", "n_sim_requested", "nPerm", "nsim"
  ))
  successful <- stage2d2_pick(value, c(
    "nsim_successful", "successful_simulations", "n_sim_successful",
    "successful_count"
  ))
  failed <- stage2d2_pick(value, c(
    "nsim_failed", "failed_simulations", "n_sim_failed", "failed_count"
  ))
  failure_reason <- stage2d2_pick(value, c(
    "failure_reason", "failure_reasons", "failure_reason_by_replicate",
    "failure_message", "failure_messages", "failure", "reason"
  ))
  list(
    K = stage2d2_pick(value, c("K", "K_fast", "observed", "observed_K")),
    P = stage2d2_pick(value, c("P", "P_fast", "p_value", "p")),
    MCSE_P = stage2d2_pick(value, c("MCSE_P", "P_MCSE", "p_mcse", "mcse")),
    exceedance = stage2d2_pick(value, c(
      "exceedance_count", "exceedance", "n_exceed", "tail_count"
    )),
    requested = requested,
    successful = successful,
    failed = failed,
    success = success,
    sim_K = stage2d2_as_sim_matrix(sim, nsim, p),
    status = stage2d2_pick(value, c("status", "result_status", "overall_status")),
    failure_reason = failure_reason,
    n_randomizations = stage2d2_pick(value, c(
      "n_randomizations", "n_random", "randomizations"
    )),
    include_observed = stage2d2_pick(value, c("include_observed", "identity_first")),
    permutation_mode = stage2d2_pick(value, c("permutation_mode")),
    names = names(value)
  )
}

stage2d2_exact_equal <- function(x, y) {
  if (is.null(x) || is.null(y)) return(is.null(x) && is.null(y))
  if (is.numeric(x) && is.numeric(y)) {
    if (!identical(dim(x), dim(y)) || length(x) != length(y)) return(FALSE)
    if (length(x) == 0L) return(TRUE)
    na_x <- is.na(x)
    na_y <- is.na(y)
    if (!identical(na_x, na_y)) return(FALSE)
    nan_x <- is.nan(x)
    nan_y <- is.nan(y)
    if (!identical(nan_x, nan_y)) return(FALSE)
    finite_x <- !na_x
    if (any(finite_x)) {
      block <- 65536L
      for (start in seq.int(1L, length(x), by = block)) {
        finish <- min(length(x), start + block - 1L)
        idx <- start:finish
        idx <- idx[finite_x[idx]]
        if (length(idx) && !identical(as.numeric(x[idx]),
                                      as.numeric(y[idx]))) return(FALSE)
      }
    }
    return(TRUE)
  }
  if (is.logical(x) && is.logical(y)) return(identical(x, y))
  if (is.character(x) && is.character(y)) return(identical(x, y))
  identical(x, y)
}

stage2d2_max_abs <- function(x, y) {
  if (is.null(x) || is.null(y)) return(NA_real_)
  if (length(x) != length(y)) return(Inf)
  if (!length(x)) return(0)
  block <- 65536L
  maximum <- 0
  saw_finite <- FALSE
  for (start in seq.int(1L, length(x), by = block)) {
    finish <- min(length(x), start + block - 1L)
    idx <- start:finish
    z <- abs(as.numeric(x[idx]) - as.numeric(y[idx]))
    finite <- !is.na(z)
    if (any(finite)) {
      saw_finite <- TRUE
      maximum <- max(maximum, z[finite])
    }
  }
  if (saw_finite) maximum else 0
}

stage2d2_compare_results <- function(left, right, nsim, p) {
  a <- stage2d2_normalize_result(left, nsim, p)
  b <- stage2d2_normalize_result(right, nsim, p)
  if (is.null(a) || is.null(b)) {
    return(list(pass = is.null(a) && is.null(b), fields = logical(),
                required = FALSE, left = a, right = b,
                max_abs = c(sim_K = NA_real_)))
  }
  fields <- c("K", "P", "MCSE_P", "exceedance", "requested",
              "successful", "failed", "success", "sim_K", "status",
              "failure_reason", "n_randomizations", "include_observed",
              "permutation_mode")
  field_pass <- vapply(fields, function(nm) {
    stage2d2_exact_equal(a[[nm]], b[[nm]])
  }, logical(1))
  required <- all(vapply(c(
    "K", "P", "MCSE_P", "exceedance", "requested", "successful",
    "failed", "success", "sim_K", "status", "failure_reason",
    "n_randomizations", "include_observed"
  ),
                         function(nm) !is.null(a[[nm]]) && !is.null(b[[nm]]),
                         logical(1)))
  list(
    pass = all(field_pass) && required,
    fields = field_pass,
    required = required,
    left = a,
    right = b,
    max_abs = c(
      K = stage2d2_max_abs(a$K, b$K),
      P = stage2d2_max_abs(a$P, b$P),
      MCSE_P = stage2d2_max_abs(a$MCSE_P, b$MCSE_P),
      exceedance = stage2d2_max_abs(a$exceedance, b$exceedance),
      sim_K = stage2d2_max_abs(a$sim_K, b$sim_K)
    )
  )
}

stage2d2_compare_captures <- function(left, right, nsim, p) {
  detail <- if (left$status == "ok" && right$status == "ok")
    stage2d2_compare_results(left$value, right$value, nsim, p) else
    list(pass = identical(left$status, right$status) &&
           identical(left$error_class, right$error_class) &&
           identical(left$error_message, right$error_message),
         fields = logical(), required = TRUE,
         left = NULL, right = NULL, max_abs = numeric())
  warning_pass <- identical(left$warning_classes, right$warning_classes) &&
    identical(left$warning_text, right$warning_text)
  status_pass <- identical(left$status, right$status)
  error_pass <- if (left$status == "error" && right$status == "error")
    identical(left$error_class, right$error_class) &&
      identical(left$error_message, right$error_message) else TRUE
  list(
    pass = isTRUE(detail$pass) && status_pass && warning_pass && error_pass,
    result = detail,
    status_pass = status_pass,
    warning_pass = warning_pass,
    error_pass = error_pass,
    left_status = left$status,
    right_status = right$status,
    left_warning = left$warning_text,
    right_warning = right$warning_text,
    left_warning_classes = left$warning_classes,
    right_warning_classes = right$warning_classes,
    left_error = left$error_message,
    right_error = right$error_message
  )
}

stage2d2_validate_payload <- function(capture, nsim, p) {
  if (capture$status != "ok") return(FALSE)
  result <- stage2d2_normalize_result(capture$value, nsim, p)
  if (is.null(result) || is.null(result$sim_K)) return(FALSE)
  if (!identical(dim(result$sim_K), c(as.integer(nsim), as.integer(p)))) {
    return(FALSE)
  }
  if (is.null(result$success) || length(result$success) != nsim) {
    return(FALSE)
  }
  required <- c("K", "P", "MCSE_P", "exceedance", "requested",
                "successful", "failed", "status", "failure_reason",
                "n_randomizations", "include_observed")
  fields_present <- all(vapply(required, function(nm) !is.null(result[[nm]]),
                               logical(1)))
  if (!fields_present) return(FALSE)
  if (length(result$requested) != 1L || length(result$n_randomizations) != 1L) {
    return(FALSE)
  }
  if (!is.logical(result$include_observed) ||
      length(result$include_observed) != 1L ||
      is.na(result$include_observed)) return(FALSE)
  all(vapply(c("successful", "failed"), function(nm) {
    length(result[[nm]]) %in% c(1L, p)
  }, logical(1))) &&
    length(result$K) == p && length(result$P) == p &&
    length(result$MCSE_P) == p && length(result$exceedance) == p
}

stage2d2_replicate_rows <- function(comparison, config) {
  if (is.null(comparison$result$left) || is.null(comparison$result$right)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  a <- comparison$result$left
  b <- comparison$result$right
  if (is.null(a$sim_K) || is.null(b$sim_K)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  nsim <- nrow(a$sim_K)
  p <- ncol(a$sim_K)
  success_a <- a$success
  success_b <- b$success
  if (is.null(success_a) || is.null(success_b)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  max_rows <- suppressWarnings(as.integer(stage2d2_env(
    "FASTPHYLOSIG_STAGE2D2_REPLICATE_LOG_MAX", "100000"
  )))
  if (is.na(max_rows) || max_rows < 1L) max_rows <- 100000L
  row_limit <- min(as.double(nsim) * as.double(p), as.double(max_rows))
  out <- list()
  written <- 0L
  for (i in seq_len(nsim)) for (j in seq_len(p)) {
    if (written >= row_limit) break
    out[[length(out) + 1L]] <- data.frame(
      check = config$check, candidate = config$candidate,
      shape = config$shape, n = config$n, traits = p,
      nsim = nsim, seed = if (is.null(config$seed)) NA_integer_ else config$seed,
      n_threads = config$n_threads, simulation_chunk = config$chunk,
      trait_kind = config$trait_kind, mask_kind = config$mask_kind,
      replicate = i, trait = j,
      oracle_null_K = a$sim_K[i, j], candidate_null_K = b$sim_K[i, j],
      null_exact = stage2d2_exact_equal(a$sim_K[i, j], b$sim_K[i, j]),
      oracle_success = if (length(success_a) >= i) success_a[[i]] else NA,
      candidate_success = if (length(success_b) >= i) success_b[[i]] else NA,
      success_exact = if (length(success_a) >= i && length(success_b) >= i)
        stage2d2_exact_equal(success_a[[i]], success_b[[i]]) else FALSE,
      stringsAsFactors = FALSE
    )
    written <- written + 1L
  }
  stage2d2_rbind(out)
}

stage2d2_memory_contract <- function() {
  max_rows <- suppressWarnings(as.integer(stage2d2_env(
    "FASTPHYLOSIG_STAGE2D2_REPLICATE_LOG_MAX", "100000"
  )))
  if (is.na(max_rows) || max_rows < 1L) max_rows <- 100000L
  data.frame(
    ordered_null_policy = "one returned payload per side; no retained copy",
    comparison_block_elements = 65536L,
    replicate_log_max_rows = max_rows,
    formal_grid_in_this_runner = "NOT_RUN",
    public_test_false_batch_guard = "INDEPENDENT_NOT_RUN_HERE",
    stringsAsFactors = FALSE
  )
}

stage2d2_provenance <- function(repo, args, hooks) {
  commit <- tryCatch(trimws(system2(
    "git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE
  )), error = function(e) "<unavailable>")
  data.frame(
    stage = "Stage 2D2 exactness correctness",
    level = args$level,
    source_commit = commit,
    package_version = tryCatch(as.character(utils::packageVersion(
      "fastphylosig")), error = function(e) NA_character_),
    R_version = R.version.string, platform = R.version$platform,
    os = Sys.info()[["sysname"]], machine = Sys.info()[["machine"]],
    omp = stage2d2_env("OMP_NUM_THREADS", ""),
    locale_collate = Sys.getlocale("LC_COLLATE"),
    hook_file = hooks$hook_file %||% "",
    prototype_cpp = hooks$cpp_file %||% "",
    same_translation_unit = isTRUE(hooks$same_translation_unit),
    candidate_b_gate = hooks$candidate_b_gate %||% "NOT_AUTHORIZED",
    candidate_b_status = hooks$candidate_b_status %||% "NOT_AUTHORIZED",
    memory_ordered_null_policy = "one payload per side; no retained copy",
    memory_comparison_block_elements = 65536L,
    memory_replicate_log_max_rows = suppressWarnings(as.integer(stage2d2_env(
      "FASTPHYLOSIG_STAGE2D2_REPLICATE_LOG_MAX", "100000"
    ))),
    public_test_false_batch_guard = "INDEPENDENT_NOT_RUN_HERE",
    production_code_changed = "NO",
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
