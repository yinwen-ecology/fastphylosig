# Shared helpers for the fastphylosig Stage 2D2B correctness gate.
#
# This file is audit infrastructure only.  It does not modify the package
# namespace, production source, tests, public API, or the frozen Stage 2D2
# evidence.  The oracle and Candidate B must be sourceCpp exports from one
# translation unit so candidate provenance cannot be satisfied by a fallback
# call to the production evaluator.

# Reuse the frozen Stage 2D2 fixture and comparison helpers.  `local = TRUE`
# keeps the imported symbols private to the runner environment.
stage2d2b_common_file <- tryCatch(
  normalizePath(sys.frame(1L)$ofile, winslash = "/", mustWork = FALSE),
  error = function(e) ""
)
stage2d2b_repo_from_common <- function() {
  if (nzchar(stage2d2b_common_file)) {
    return(normalizePath(file.path(dirname(stage2d2b_common_file), "..",
                                   "stage2d2"),
                        winslash = "/", mustWork = FALSE))
  }
  file.path(getwd(), "benchmarks", "stage2d2")
}
source(file.path(stage2d2b_repo_from_common(), "stage2d2_common.R"),
       local = TRUE)

stage2d2b_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

stage2d2b_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  smoke <- "--smoke" %in% args
  correctness <- "--correctness" %in% args
  if (smoke == correctness) {
    stop("Use exactly one of --smoke or --correctness.", call. = FALSE)
  }
  positional <- args[!args %in% c("--smoke", "--correctness")]
  if (length(positional) > 2L) {
    stop("usage: run_stage2d2b_correctness.R [--smoke|--correctness] " ,
         "[repo] [out].", call. = FALSE)
  }
  repo <- if (length(positional)) positional[[1L]] else getwd()
  out <- if (length(positional) >= 2L) positional[[2L]] else
    file.path(repo, "benchmarks", "stage2d2b", "results", "correctness")
  list(
    mode = if (isTRUE(correctness)) "formal" else "smoke",
    smoke = isTRUE(smoke), correctness = isTRUE(correctness),
    repo = normalizePath(repo, winslash = "/", mustWork = FALSE),
    out = normalizePath(out, winslash = "/", mustWork = FALSE)
  )
}

stage2d2b_ints <- function(value, name, minimum = 1L) {
  out <- suppressWarnings(as.integer(trimws(strsplit(
    as.character(value), ",", fixed = TRUE
  )[[1L]])))
  if (!length(out) || anyNA(out) || any(out < minimum)) {
    stop(name, " must contain comma-separated integers >= ", minimum,
         ".", call. = FALSE)
  }
  unique(out)
}

stage2d2b_strings <- function(value, name, allowed = NULL) {
  out <- unique(trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1L]]))
  if (!length(out) || any(!nzchar(out)) ||
      (!is.null(allowed) && any(!out %in% allowed))) {
    suffix <- if (is.null(allowed)) "" else paste0(
      " Allowed: ", paste(allowed, collapse = ", "), "."
    )
    stop(name, " must contain comma-separated values.", suffix, call. = FALSE)
  }
  out
}

stage2d2b_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "NA")
  invisible(path)
}

stage2d2b_rbind <- function(rows) {
  if (!length(rows)) return(data.frame(stringsAsFactors = FALSE))
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(cols, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[cols]
  })
  do.call(rbind, rows)
}

stage2d2b_capture <- function(fun) {
  # Keep warning classes/text as part of the failure contract.  The warning is
  # muffled only after it has been recorded, so the runner remains readable.
  warnings <- list()
  status <- "ok"
  value <- NULL
  error_message <- ""
  error_class <- character()
  started <- unname(as.numeric(proc.time()[["elapsed"]]))
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
  elapsed <- max(0, unname(as.numeric(proc.time()[["elapsed"]])) - started)
  warning_text <- vapply(warnings, `[[`, character(1), "text")
  warning_classes <- unique(unlist(lapply(warnings, `[[`, "classes"),
                                  use.names = FALSE))
  list(
    value = value, elapsed_s = elapsed, status = status,
    warnings = warnings, warning_text = paste(unique(warning_text), collapse = " | "),
    warning_classes = paste(warning_classes, collapse = "/"),
    error_message = error_message, error_class = paste(error_class, collapse = "/")
  )
}

stage2d2b_call_with_seed <- function(fun, args, seed = NULL) {
  invoke <- function() do.call(fun, stage2d2_filter_call_args(fun, args))
  if (is.null(seed)) return(invoke())
  stage2d2_restore_seed(seed, invoke)
}

stage2d2b_call_engine <- function(fun, compiled_tree, X, nsim,
                                  permutations = NULL, block_size = 2L,
                                  n_threads = 1L, trait_chunk = 64L,
                                  simulation_chunk = 128L,
                                  include_observed = FALSE, seed = NULL,
                                  return_sim = TRUE) {
  fml <- tryCatch(names(formals(fun)), error = function(e) character())
  args <- list(
    compiled_tree = compiled_tree, X = X, nsim = as.integer(nsim),
    permutations = permutations, trait_chunk = as.integer(trait_chunk),
    return_sim = isTRUE(return_sim), include_observed = isTRUE(include_observed),
    n_threads = as.integer(n_threads), simulation_chunk = as.integer(simulation_chunk),
    return_ordered = TRUE, block_size = as.integer(block_size),
    collect_counters = TRUE
  )
  if (!("n_threads" %in% fml) && "ncores" %in% fml) {
    args$ncores <- args$n_threads
    args$n_threads <- NULL
  }
  stage2d2b_call_with_seed(fun, args, seed)
}

stage2d2b_find_export <- function(env, candidates) {
  for (name in candidates) {
    if (exists(name, envir = env, inherits = FALSE)) {
      value <- get(name, envir = env, inherits = FALSE)
      if (is.function(value)) return(list(name = name, fun = value))
    }
  }
  NULL
}

stage2d2b_mode_wrapper <- function(target, mode) {
  force(target)
  force(mode)
  function(...) {
    args <- list(...)
    args$mode <- mode
    do.call(target, stage2d2_filter_call_args(target, args))
  }
}

stage2d2b_load_hooks <- function(repo) {
  cpp_file <- stage2d2b_env("FASTPHYLOSIG_STAGE2D2B_PROTO_CPP")
  if (!nzchar(cpp_file)) {
    stop("Set FASTPHYLOSIG_STAGE2D2B_PROTO_CPP to the one Candidate B " ,
         "audit translation unit.", call. = FALSE)
  }
  cpp_file <- normalizePath(cpp_file, winslash = "/", mustWork = TRUE)
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Rcpp is required for FASTPHYLOSIG_STAGE2D2B_PROTO_CPP.",
         call. = FALSE)
  }
  hook_env <- new.env(parent = globalenv())
  before <- ls(envir = hook_env, all.names = TRUE)
  Rcpp::sourceCpp(cpp_file, env = hook_env, rebuild = TRUE,
                  showOutput = FALSE, verbose = FALSE)
  exports <- setdiff(ls(envir = hook_env, all.names = TRUE), before)
  cpp_bound <- function(entry) {
    !is.null(entry) && entry$name %in% exports &&
      exists(entry$name, envir = hook_env, inherits = FALSE) &&
      identical(get(entry$name, envir = hook_env, inherits = FALSE), entry$fun)
  }

  generic <- stage2d2b_find_export(hook_env, c(
    "stage2d2b_k_permutation_prototype", "stage2d2b_k_null_prototype",
    "stage2d2_k_permutation_prototype", "stage2d2_k_null_prototype"
  ))
  oracle_raw <- stage2d2b_find_export(hook_env, c(
    "stage2d2b_k_permutation_oracle", "stage2d2b_k_null_oracle",
    "stage2d2_k_permutation_oracle", "stage2d2_k_null_oracle", "oracle"
  ))
  candidate_raw <- stage2d2b_find_export(hook_env, c(
    "stage2d2b_k_permutation_candidate_b", "stage2d2b_k_null_candidate_b",
    "stage2d2_k_permutation_candidate_b", "stage2d2_k_null_candidate_b",
    "candidate_b"
  ))
  if (!is.null(generic)) {
    if (!cpp_bound(generic)) stop("Generic export is not sourceCpp-bound.",
                                 call. = FALSE)
    oracle <- stage2d2b_mode_wrapper(generic$fun, "oracle")
    candidate <- stage2d2b_mode_wrapper(generic$fun, "B")
    oracle_name <- paste0(generic$name, "(mode=oracle)")
    candidate_name <- paste0(generic$name, "(mode=B)")
    same_tu <- cpp_bound(generic)
    layout <- "generic_mode"
  } else {
    if (is.null(oracle_raw) || is.null(candidate_raw) ||
        !cpp_bound(oracle_raw) || !cpp_bound(candidate_raw)) {
      stop("Candidate B requires sourceCpp-bound oracle and candidate_b " ,
           "exports from the same translation unit.", call. = FALSE)
    }
    oracle <- oracle_raw$fun
    candidate <- candidate_raw$fun
    oracle_name <- oracle_raw$name
    candidate_name <- candidate_raw$name
    same_tu <- cpp_bound(oracle_raw) && cpp_bound(candidate_raw)
    layout <- "separate_exports"
  }
  if (!isTRUE(same_tu)) stop("The oracle/Candidate B pair is not same-TU bound.",
                             call. = FALSE)
  gate <- toupper(trimws(stage2d2b_env(
    "FASTPHYLOSIG_STAGE2D2B_GATE", "NOT_AUTHORIZED"
  )))
  authorized <- gate %in% c("GO", "AUTHORIZED", "RUN", "PASS")
  if (!authorized) {
    stop("Stage 2D2B Candidate B requires an authorized gate: " ,
         "FASTPHYLOSIG_STAGE2D2B_GATE=AUTHORIZED.", call. = FALSE)
  }
  list(
    cpp_file = cpp_file, env = hook_env, exports = exports,
    mode = "B", layout = layout, oracle = oracle, candidate = candidate,
    oracle_name = oracle_name, candidate_name = candidate_name,
    same_translation_unit = same_tu, gate = gate
  )
}

stage2d2b_load_package <- function(repo) {
  library_dir <- stage2d2b_env("FASTPHYLOSIG_STAGE2D2B_LIBRARY")
  if (nzchar(library_dir)) {
    library_dir <- normalizePath(library_dir, winslash = "/", mustWork = TRUE)
    .libPaths(c(library_dir, .libPaths()))
    suppressPackageStartupMessages(library(
      "fastphylosig", lib.loc = library_dir, character.only = TRUE
    ))
  } else {
    stage2d2d_library <- stage2d2_env("FASTPHYLOSIG_STAGE2D2_LIBRARY")
    if (nzchar(stage2d2d_library)) {
      library_dir <- normalizePath(stage2d2d_library, winslash = "/",
                                   mustWork = TRUE)
      .libPaths(c(library_dir, .libPaths()))
      suppressPackageStartupMessages(library(
        "fastphylosig", lib.loc = library_dir, character.only = TRUE
      ))
    } else {
      if (!requireNamespace("pkgload", quietly = TRUE)) {
        stop("pkgload is required when no audited library is configured.",
             call. = FALSE)
      }
      pkgload::load_all(repo, quiet = TRUE)
    }
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("ape is required for Stage 2D2B fixtures.", call. = FALSE)
  }
  asNamespace("fastphylosig")
}

stage2d2b_extract_named <- function(value, keys, depth = 0L) {
  if (!is.list(value) || depth > 3L) return(NULL)
  hit <- intersect(keys, names(value))
  if (length(hit)) return(value[[hit[[1L]]]])
  # Only descend through metadata containers.  Never walk large simulation
  # payloads while looking for audit counters.
  containers <- intersect(c("candidate_metadata", "counters", "audit",
                            "diagnostics", "metadata", "provenance",
                            "engine_metadata", "memory"), names(value))
  for (name in containers) {
    out <- stage2d2b_extract_named(value[[name]], keys, depth + 1L)
    if (!is.null(out)) return(out)
  }
  NULL
}

stage2d2b_counter <- function(value, keys) {
  raw <- stage2d2b_extract_named(value, keys)
  if (is.null(raw) || !length(raw)) {
    return(list(present = FALSE, value = NA_real_))
  }
  numeric <- suppressWarnings(as.numeric(raw[[1L]]))
  if (!length(numeric) || is.na(numeric)) {
    return(list(present = FALSE, value = NA_real_))
  }
  list(present = TRUE, value = numeric[[1L]])
}

stage2d2b_counter_row <- function(captured, cfg) {
  value <- captured$value
  old <- stage2d2b_counter(value, c(
    "candidate_old_compute_one_calls", "old_compute_one_calls",
    "candidate_compute_one_calls"
  ))
  oracle <- stage2d2b_counter(value, c(
    "oracle_call_count", "candidate_oracle_call_count"
  ))
  compute <- stage2d2b_counter(value, c(
    "candidate_compute_count", "candidate_compute_calls",
    "compute_block_count", "blocks_processed"
  ))
  block <- stage2d2b_counter(value, c("block_size"))
  generated <- stage2d2b_counter(value, c(
    "generated_permutations", "permutations_generated"
  ))
  old_equivalent <- stage2d2b_counter(value, c("compute_one_equivalents"))
  present <- old$present && oracle$present && compute$present
  pass <- captured$status == "ok" && present && old$value == 0 &&
    oracle$value == 0 && compute$value > 0
  data.frame(
    check = cfg$check, implementation = cfg$implementation,
    shape = cfg$shape, n = as.integer(cfg$n), traits = as.integer(cfg$traits),
    nsim = as.integer(cfg$nsim), pattern = cfg$pattern,
    block_size = as.integer(cfg$block_size), n_threads = as.integer(cfg$n_threads),
    seed = if (is.null(cfg$seed)) NA_integer_ else as.integer(cfg$seed),
    old_compute_one_calls = old$value, old_counter_present = old$present,
    oracle_call_count = oracle$value, oracle_counter_present = oracle$present,
    candidate_compute_count = compute$value,
    candidate_counter_present = compute$present, reported_block_size = block$value,
    generated_permutations = generated$value,
    compute_one_equivalents = old_equivalent$value,
    counter_gate_pass = pass, status = captured$status,
    error_message = captured$error_message, stringsAsFactors = FALSE
  )
}

stage2d2b_permutation_patterns <- c(
  "identity", "single_swap", "reverse", "cyclic_shift", "fixed_random",
  "adversarial"
)

stage2d2b_make_permutations <- function(n, nsim, pattern, seed = 1L) {
  n <- as.integer(n)
  nsim <- as.integer(nsim)
  if (n < 2L || nsim < 1L || !pattern %in% stage2d2b_permutation_patterns) {
    stop("invalid controlled permutation fixture.", call. = FALSE)
  }
  stage2d2_restore <- function(fun) stage2d2_restore_seed(seed, fun)
  stage2d2_restore(function() {
    identity <- seq_len(n)
    out <- matrix(identity, nrow = nsim, ncol = n, byrow = TRUE)
    if (nsim < 2L) return(out)
    for (i in 2:nsim) {
      row <- identity
      if (pattern == "single_swap") {
        left <- ((i - 2L) %% (n - 1L)) + 1L
        row[c(left, left + 1L)] <- row[c(left + 1L, left)]
      } else if (pattern == "reverse") {
        row <- rev(identity)
      } else if (pattern == "cyclic_shift") {
        shift <- (i - 1L) %% n
        if (shift > 0L) row <- c(identity[(shift + 1L):n], identity[seq_len(shift)])
      } else if (pattern == "fixed_random") {
        row <- sample.int(n, n, replace = FALSE)
      } else if (pattern == "adversarial") {
        mode <- (i - 2L) %% 4L
        if (mode == 0L) {
          row[c(1L, n)] <- row[c(n, 1L)]
        } else if (mode == 1L) {
          shift <- max(1L, floor(n / 2L))
          row <- c(identity[(shift + 1L):n], identity[seq_len(shift)])
        } else if (mode == 2L) {
          row <- rev(identity)
          if (n >= 4L) row[c(2L, n - 1L)] <- row[c(n - 1L, 2L)]
        } else {
          even <- seq.int(2L, n, by = 2L)
          odd <- seq.int(1L, n, by = 2L)
          row <- c(even, odd)
        }
      }
      out[i, ] <- row
    }
    out
  })
}

stage2d2b_normalize <- function(value, nsim, p) {
  stage2d2_normalize_result(value, nsim, p)
}

stage2d2b_numeric_equal <- function(x, y) stage2d2_exact_equal(x, y)

stage2d2b_compare <- function(old, new, nsim, p) {
  compared <- stage2d2_compare_captures(old, new, nsim, p)
  # This gate is deliberately explicit for the production contract fields;
  # a candidate may not pass merely because a converter omitted P or MCSE.
  a <- stage2d2b_normalize(old$value, nsim, p)
  b <- stage2d2b_normalize(new$value, nsim, p)
  required <- c("K", "P", "MCSE_P", "exceedance", "requested",
                "successful", "failed", "success", "sim_K", "status",
                "failure_reason", "n_randomizations", "include_observed")
  fields_present <- old$status == "ok" && new$status == "ok" &&
    !is.null(a) && !is.null(b) && all(vapply(required, function(name) {
      !is.null(a[[name]]) && !is.null(b[[name]])
    }, logical(1)))
  compared$contract_fields_present <- fields_present
  compared$pass <- isTRUE(compared$pass) && isTRUE(fields_present)
  compared
}

stage2d2b_accounting <- function(captured, nsim, p, include_observed,
                                  controlled = FALSE, require_identity = FALSE) {
  if (captured$status != "ok") return(list(pass = FALSE, detail = "error"))
  value <- stage2d2b_normalize(captured$value, nsim, p)
  if (is.null(value)) return(list(pass = FALSE, detail = "missing result"))
  scalar <- function(x) if (is.null(x) || length(x) != 1L) NA_real_ else
    suppressWarnings(as.numeric(x[[1L]]))
  requested <- scalar(value$requested)
  randomizations <- scalar(value$n_randomizations)
  expected_randomizations <- if (controlled) nsim else
    if (include_observed) max(0, nsim - 1L) else nsim
  count_ok <- is.finite(requested) && requested == nsim &&
    is.finite(randomizations) && randomizations == expected_randomizations
  successful <- suppressWarnings(as.numeric(value$successful))
  failed <- suppressWarnings(as.numeric(value$failed))
  count_ok <- count_ok && length(successful) %in% c(1L, p) &&
    length(failed) %in% c(1L, p)
  if (count_ok) {
    count_ok <- all((rep(successful, length.out = max(length(successful),
                                                     length(failed))) +
                     rep(failed, length.out = max(length(successful),
                                                   length(failed)))) == nsim)
  }
  success_ok <- !is.null(value$success) && length(value$success) == nsim
  include_ok <- !is.null(value$include_observed) &&
    length(value$include_observed) == 1L &&
    identical(as.logical(value$include_observed), isTRUE(include_observed))
  identity_ok <- TRUE
  if (require_identity && !is.null(value$sim_K) && !is.null(value$K) &&
      nrow(value$sim_K) >= 1L) {
    identity_ok <- stage2d2b_numeric_equal(value$sim_K[1L, ], value$K)
  }
  p_ok <- TRUE
  if (!is.null(value$P) && !is.null(value$exceedance) &&
      is.finite(requested) && requested > 0 &&
      is.finite(randomizations) && randomizations >= 0) {
    # fast_k's frozen K contract reports exceedances over the requested
    # ordered rows, including the identity row when it is present.  The
    # denominator is therefore nsim, while n_randomizations is retained as
    # metadata for the MCSE contract.
    expected_p <- as.numeric(value$exceedance) /
      as.numeric(value$requested[[1L]])
    finite <- is.finite(as.numeric(value$P)) & is.finite(expected_p)
    if (any(finite)) p_ok <- stage2d2b_numeric_equal(
      as.numeric(value$P)[finite], expected_p[finite]
    )
  }
  list(
    pass = isTRUE(count_ok) && isTRUE(success_ok) && isTRUE(include_ok) &&
      isTRUE(identity_ok) && isTRUE(p_ok),
    requested_ok = isTRUE(is.finite(requested) && requested == nsim),
    randomizations_ok = isTRUE(is.finite(randomizations) &&
                               randomizations == expected_randomizations),
    counts_ok = isTRUE(count_ok), success_ok = isTRUE(success_ok),
    include_observed_ok = isTRUE(include_ok), identity_first_ok = isTRUE(identity_ok),
    p_formula_ok = isTRUE(p_ok), detail = paste(
      "requested", requested, "randomizations", randomizations,
      "identity", identity_ok, "P_formula", p_ok
    )
  )
}

stage2d2b_replicate_rows <- function(comparison, cfg) {
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
  out <- list()
  for (i in seq_len(nsim)) for (j in seq_len(p)) {
    oracle <- a$sim_K[i, j]
    candidate <- b$sim_K[i, j]
    out[[length(out) + 1L]] <- data.frame(
      check = cfg$check, shape = cfg$shape, n = as.integer(cfg$n),
      traits = p, nsim = nsim, pattern = cfg$pattern,
      block_size = as.integer(cfg$block_size), n_threads = as.integer(cfg$n_threads),
      replicate = i, trait = j, oracle_null_K = oracle,
      candidate_null_K = candidate,
      bitwise_exact = stage2d2b_numeric_equal(oracle, candidate),
      absolute_difference = stage2d2_max_abs(oracle, candidate),
      oracle_success = if (length(a$success) >= i) a$success[[i]] else NA,
      candidate_success = if (length(b$success) >= i) b$success[[i]] else NA,
      success_exact = if (length(a$success) >= i && length(b$success) >= i)
        identical(a$success[[i]], b$success[[i]]) else FALSE,
      stringsAsFactors = FALSE
    )
  }
  stage2d2b_rbind(out)
}

stage2d2b_source_sha256 <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256"))
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    return(paste(as.character(openssl::sha256(file(path))), collapse = ""))
  }
  NA_character_
}

stage2d2b_git_value <- function(repo) {
  value <- tryCatch(system2("git", c("-C", repo, "rev-parse", "HEAD"),
                            stdout = TRUE, stderr = FALSE),
                    error = function(e) character())
  if (length(value)) paste(value, collapse = "\n") else "UNAVAILABLE"
}
