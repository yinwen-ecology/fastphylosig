#!/usr/bin/env Rscript

# fastphylosig 0.2.0 Stage 2C method profiling harness.
#
# This file is evidence-only.  It loads an installed package, measures the
# lambda, D, and Delta public routes, and writes CSV evidence below the
# requested output directory.  It never edits package source or tests.  Full
# runs are deliberately serial: one R process owns the benchmark at a time.

options(stringsAsFactors = FALSE)

raw_args <- commandArgs(trailingOnly = TRUE)
quick_mode <- "quick" %in% raw_args
help_mode <- any(raw_args %in% c("--help", "-h"))
args <- raw_args[!raw_args %in% c("quick", "--help", "-h")]

usage <- paste(
  "usage: run_stage2c_methods.R LIBRARY OUTPUT_DIR [quick]",
  "",
  "The installed library must contain fastphylosig and ape.",
  "quick runs a tiny smoke-sized grid; omit it for the Stage 2C grid.",
  sep = "\n"
)
if (help_mode) {
  cat(usage, "\n")
  quit(save = "no", status = 0L)
}
if (length(args) < 2L) {
  cat(usage, "\n", file = stderr())
  quit(save = "no", status = 2L)
}

library_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Keep benchmark ordering and string handling stable.  The original locale is
# retained in provenance before the benchmark process is made deterministic.
locale_before <- tryCatch(Sys.getlocale(), error = function(e) "<unavailable>")
Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")
try(Sys.setlocale("LC_COLLATE", "C"), silent = TRUE)

base_libs <- .libPaths()
dependency_user_libs <- base_libs[grepl(
  "win-library|R-library|site-library", base_libs
)]
Sys.setenv(R_LIBS_USER = paste(
  c(library_dir, dependency_user_libs), collapse = .Platform$path.sep
))
.libPaths(unique(c(library_dir, base_libs)))

suppressPackageStartupMessages(library(fastphylosig))
suppressPackageStartupMessages(library(ape))

if (!quick_mode && !requireNamespace("callr", quietly = TRUE)) {
  # callr is not needed by the harness itself, but a missing package is not a
  # reason to abort an otherwise valid serial run.  Record it below.
  callr_available <- FALSE
} else {
  callr_available <- requireNamespace("callr", quietly = TRUE)
}

ns <- asNamespace("fastphylosig")
pkg_version <- as.character(utils::packageVersion("fastphylosig"))

env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

parse_int_vector <- function(name, default, min_value = 1L) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  out <- suppressWarnings(as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])))
  if (!length(out) || anyNA(out) || any(out < min_value)) {
    stop(name, " must contain comma-separated integers >= ", min_value,
         call. = FALSE)
  }
  unique(out)
}

positive_bytes_option <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = "")))
  if (length(value) != 1L || is.na(value) || !is.finite(value) || value <= 0) {
    return(as.numeric(default))
  }
  value
}

# phytools::phylosig(lambda) uses a dense n-by-n covariance representation.
# Keep the reference path bounded and record skipped cells explicitly instead
# of warming up or timing an allocation that is known to exceed the policy.
reference_dense_cap_bytes <- positive_bytes_option(
  "FASTPHYLOSIG_STAGE2C_REFERENCE_DENSE_CAP_BYTES", 512 * 1024^2
)
reference_cubic_work_cap <- positive_bytes_option(
  "FASTPHYLOSIG_STAGE2C_LAMBDA_REFERENCE_CUBIC_CAP", 1e9
)

reference_precheck <- function(n) {
  dense_bytes <- as.double(n) * as.double(n) * 8
  # Lambda optimization evaluates the dense likelihood repeatedly.  Twenty
  # evaluations is a conservative preflight proxy, not a measured FLOP count.
  estimated_cubic_work <- 20 * as.double(n)^3
  reason <- if (as.integer(n) >= 10000L) {
    sprintf(
      "single phytools lambda reference is not run at n >= 10000; estimated dense covariance %.0f bytes",
      dense_bytes
    )
  } else if (dense_bytes > reference_dense_cap_bytes) {
    sprintf(
      "estimated dense covariance %.0f bytes exceeds cap %.0f bytes",
      dense_bytes, reference_dense_cap_bytes
    )
  } else if (estimated_cubic_work > reference_cubic_work_cap) {
    sprintf(
      paste0(
        "estimated repeated dense cubic work %.0f exceeds lambda ",
        "reference cap %.0f"
      ),
      estimated_cubic_work, reference_cubic_work_cap
    )
  } else {
    ""
  }
  list(
    limited = nzchar(reason), estimated_dense_bytes = dense_bytes,
    estimated_cubic_work = estimated_cubic_work,
    reason = reason
  )
}

parse_shape_vector <- function(env_name, default) {
  value <- Sys.getenv(env_name, unset = "")
  if (!nzchar(value)) return(default)
  out <- unique(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
  out <- out[nzchar(out)]
  if (!length(out) || any(!out %in% c("balanced", "random", "pectinate"))) {
    stop(env_name, " must contain balanced, random, and/or pectinate.",
         call. = FALSE)
  }
  out
}

method_text <- Sys.getenv(
  "FASTPHYLOSIG_STAGE2C_METHODS", unset = "lambda,D,Delta"
)
method_tokens <- unique(tolower(trimws(strsplit(method_text, ",", fixed = TRUE)[[1L]])))
method_tokens <- method_tokens[nzchar(method_tokens)]
method_alias <- c(lambda = "lambda", d = "D", delta = "Delta")
if (!length(method_tokens) || any(!method_tokens %in% names(method_alias))) {
  stop("FASTPHYLOSIG_STAGE2C_METHODS must contain lambda, d, and/or delta.",
       call. = FALSE)
}
selected_methods <- unname(method_alias[method_tokens])
run_lambda <- "lambda" %in% selected_methods
run_d <- "D" %in% selected_methods
run_delta <- "Delta" %in% selected_methods

repeat_override <- suppressWarnings(as.integer(Sys.getenv(
  "FASTPHYLOSIG_STAGE2C_REPEATS", unset = ""
)))
repeats <- if (length(repeat_override) && is.finite(repeat_override) &&
               repeat_override > 0L) repeat_override else if (quick_mode) 1L else 10L

repeats_for <- function(method, workload, n, ncores = 1L) {
  if (length(repeat_override) && is.finite(repeat_override) &&
      repeat_override > 0L) {
    return(as.integer(repeat_override))
  }
  if (quick_mode) return(1L)
  if (identical(method, "Delta")) return(3L)
  if (as.integer(n) >= 5000L) return(5L)
  10L
}

timeout_seconds <- suppressWarnings(as.numeric(Sys.getenv(
  "FASTPHYLOSIG_STAGE2C_TIMEOUT", unset = ""
)))
if (!length(timeout_seconds) || is.na(timeout_seconds) ||
    !is.finite(timeout_seconds) || timeout_seconds <= 0) {
  timeout_seconds <- if (quick_mode) 30 else 300
}

lambda_n <- parse_int_vector(
  "FASTPHYLOSIG_STAGE2C_LAMBDA_N",
  if (quick_mode) c(100L) else c(100L, 500L, 2000L, 5000L, 10000L),
  min_value = 2L
)
d_n <- parse_int_vector(
  "FASTPHYLOSIG_STAGE2C_D_N",
  if (quick_mode) c(100L) else c(100L, 500L, 2000L),
  min_value = 2L
)
delta_n <- parse_int_vector(
  "FASTPHYLOSIG_STAGE2C_DELTA_N",
  if (quick_mode) c(50L) else c(50L, 100L, 500L),
  min_value = 2L
)
lambda_shapes <- parse_shape_vector(
  "FASTPHYLOSIG_STAGE2C_LAMBDA_SHAPES",
  c("balanced", "random", "pectinate")
)
d_shapes <- parse_shape_vector(
  "FASTPHYLOSIG_STAGE2C_D_SHAPES", c("random")
)
delta_shapes <- parse_shape_vector(
  "FASTPHYLOSIG_STAGE2C_DELTA_SHAPES", c("random")
)

d_nsim <- parse_int_vector(
  "FASTPHYLOSIG_STAGE2C_D_NSIM",
  if (quick_mode) c(199L) else c(199L, 999L, 9999L),
  min_value = 1L
)
delta_nsim <- parse_int_vector(
  "FASTPHYLOSIG_STAGE2C_DELTA_NSIM",
  if (quick_mode) c(19L) else c(19L, 199L),
  min_value = 1L
)

# These are the production Delta controls at the time of Stage 2C.  They are
# intentionally not reduced in quick mode: quick changes only the grid and
# repeat count, never the production MCMC definition.
delta_mcmc_sim <- 10000L
delta_thin <- 10L
delta_burn <- 100L

lambda_ncores <- parse_int_vector(
  "FASTPHYLOSIG_STAGE2C_LAMBDA_NCORES",
  c(1L), min_value = 1L
)
d_ncores <- parse_int_vector(
  "FASTPHYLOSIG_STAGE2C_D_NCORES",
  if (quick_mode) c(1L) else c(1L, 2L), min_value = 1L
)
delta_ncores <- parse_int_vector(
  "FASTPHYLOSIG_STAGE2C_DELTA_NCORES",
  if (quick_mode) c(1L) else c(1L, 2L), min_value = 1L
)

now <- function() unname(as.numeric(proc.time()[["elapsed"]]))

scalar_numeric <- function(value) {
  out <- tryCatch(as.numeric(value), error = function(e) numeric())
  if (length(out)) unname(out[[1L]]) else NA_real_
}

status_of <- function(value) {
  if (is.null(value)) return(NA_character_)
  if (is.data.frame(value) && "status" %in% names(value)) {
    return(as.character(value$status[[1L]]))
  }
  if (is.list(value) && !is.null(value$status)) {
    return(as.character(value$status[[1L]]))
  }
  NA_character_
}

estimate_of <- function(value, method) {
  candidates <- switch(
    method,
    lambda = c("lambda", "lambda_fast", "estimate"),
    D = c("D_fast", "D", "estimate"),
    Delta = c("Delta_fast", "delta", "Delta", "estimate")
  )
  if (is.data.frame(value)) {
    for (nm in candidates) if (nm %in% names(value)) {
      return(scalar_numeric(value[[nm]]))
    }
  }
  if (is.list(value)) {
    for (nm in candidates) if (!is.null(value[[nm]])) {
      return(scalar_numeric(value[[nm]]))
    }
  }
  scalar_numeric(value)
}

field_of <- function(value, candidates) {
  if (is.data.frame(value)) {
    for (nm in candidates) if (nm %in% names(value)) {
      return(scalar_numeric(value[[nm]]))
    }
  }
  if (is.list(value)) {
    for (nm in candidates) if (!is.null(value[[nm]])) {
      return(scalar_numeric(value[[nm]]))
    }
  }
  NA_real_
}

warning_text <- function(warnings) {
  if (!length(warnings)) "" else paste(unique(warnings), collapse = " || ")
}

run_timed <- function(fun, seed = NULL, timeout = timeout_seconds) {
  if (!is.null(seed)) set.seed(seed)
  gc(FALSE)
  warnings <- character()
  status <- "ok"
  error_message <- ""
  error_class <- ""
  value <- NULL
  started <- now()
  if (is.finite(timeout) && timeout > 0) {
    setTimeLimit(cpu = Inf, elapsed = timeout, transient = TRUE)
  }
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      status <<- if (inherits(e, "elapsedTimeLimit")) {
        "RESOURCE_LIMIT"
      } else {
        "error"
      }
      error_message <<- conditionMessage(e)
      error_class <<- paste(class(e), collapse = "/")
      NULL
    }
  )
  elapsed <- max(0, now() - started)
  if (is.finite(timeout) && timeout > 0) {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  }
  list(
    elapsed_ms = 1000 * elapsed,
    value = value,
    status = status,
    error_message = error_message,
    error_class = error_class,
    warnings = warnings
  )
}

tip_labels <- function(n) paste0("sp", seq_len(n))

make_balanced <- function(n) {
  n <- as.integer(n)
  n_internal <- n - 1L
  edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  edge_i <- 0L
  next_internal <- n + 1L
  root <- next_internal
  next_internal <- next_internal + 1L
  queue_node <- integer(n_internal)
  queue_lo <- integer(n_internal)
  queue_hi <- integer(n_internal)
  head <- 1L
  tail <- 1L
  queue_node[[1L]] <- root
  queue_lo[[1L]] <- 1L
  queue_hi[[1L]] <- n
  while (head <= tail) {
    node <- queue_node[[head]]
    lo <- queue_lo[[head]]
    hi <- queue_hi[[head]]
    head <- head + 1L
    mid <- floor((lo + hi) / 2)
    child_lo <- c(lo, mid + 1L)
    child_hi <- c(mid, hi)
    for (j in 1:2) {
      child <- if (child_lo[[j]] == child_hi[[j]]) {
        child_lo[[j]]
      } else {
        next_node <- next_internal
        next_internal <- next_internal + 1L
        tail <- tail + 1L
        queue_node[[tail]] <- next_node
        queue_lo[[tail]] <- child_lo[[j]]
        queue_hi[[tail]] <- child_hi[[j]]
        next_node
      }
      edge_i <- edge_i + 1L
      edge[edge_i, ] <- c(node, child)
    }
  }
  edge <- edge[seq_len(edge_i), , drop = FALSE]
  tree <- list(
    edge = edge,
    tip.label = tip_labels(n),
    edge.length = 0.5 + seq_len(nrow(edge)) / nrow(edge),
    Nnode = n_internal
  )
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

make_pectinate <- function(n) {
  n <- as.integer(n)
  root <- n + 1L
  edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  edge_i <- 0L
  current <- root
  next_internal <- root + 1L
  if (n == 2L) {
    edge <- matrix(c(root, 1L, root, 2L), ncol = 2L, byrow = TRUE)
  } else {
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
    tip.label = tip_labels(n),
    edge.length = 0.5 + seq_len(nrow(edge)) / nrow(edge),
    Nnode = n - 1L
  )
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

make_random <- function(n) {
  set.seed(as.integer(20260908L + n))
  tree <- ape::reorder.phylo(ape::rtree(n), order = "postorder")
  tree$tip.label <- tip_labels(n)
  tree$edge.length <- 0.5 + seq_len(nrow(tree$edge)) / nrow(tree$edge)
  tree
}

make_fixture <- function(shape, n) {
  switch(
    shape,
    balanced = make_balanced(n),
    random = make_random(n),
    pectinate = make_pectinate(n),
    stop("unknown tree shape: ", shape, call. = FALSE)
  )
}

make_trait <- function(tree) {
  n <- length(tree$tip.label)
  stats::setNames(
    sin(seq_len(n)) + cos(seq_len(n) / 7) + seq_len(n) / (100 * n),
    tree$tip.label
  )
}

make_binary <- function(tree) {
  n <- length(tree$tip.label)
  x <- (seq_len(n) %% 2L == 0L) * 1
  stats::setNames(as.numeric(x), tree$tip.label)
}

make_categorical <- function(tree) {
  n <- length(tree$tip.label)
  x <- ifelse(seq_len(n) %% 2L == 0L, "a", "b")
  stats::setNames(x, tree$tip.label)
}

clone_fixture <- function(value) unserialize(serialize(value, NULL))

call_lambda <- function(path, tree, ctx, trait, test = FALSE, ncores = 1L) {
  target <- if (identical(path, "raw")) tree else ctx
  fastphylosig::fast_lambda(
    target, trait, test = test, verbose = FALSE, progress = FALSE,
    ncores = ncores, lambda_profile = FALSE
  )
}

call_d <- function(path, tree, ctx, trait, test, nsim, ncores,
                   random_states = NULL, brownian_states = NULL) {
  target <- if (identical(path, "raw")) tree else ctx
  fastphylosig::fast_d(
    target, trait, test = test, nsim = nsim, return_sim = FALSE,
    keep_null = FALSE, verbose = FALSE, progress = FALSE, ncores = ncores,
    random_states = random_states, brownian_states = brownian_states
  )
}

call_delta <- function(path, tree, ctx, trait, test, nsim, ncores) {
  target <- if (identical(path, "raw")) tree else ctx
  fastphylosig::fast_delta(
    target, trait, test = test, nsim = nsim,
    mcmc_sim = delta_mcmc_sim, thin = delta_thin, burn = delta_burn,
    return_sim = FALSE, verbose = FALSE, progress = FALSE, ncores = ncores
  )
}

# The phase probes are intentionally separate from formal wall-time rows.
# Trace events are inclusive and therefore never presented as additive timing
# budgets.  Formal timings remain untraced and are the authoritative totals.
trace_state <- new.env(parent = emptyenv())
trace_state$stack <- list()
trace_state$events <- list()
trace_state$counter <- 0L

stage2c_trace_enter <- function(label) {
  trace_state$stack[[length(trace_state$stack) + 1L]] <- list(
    label = label, start = now()
  )
  invisible(NULL)
}

stage2c_trace_exit <- function(label) {
  if (!length(trace_state$stack)) return(invisible(NULL))
  frame <- trace_state$stack[[length(trace_state$stack)]]
  trace_state$stack <- if (length(trace_state$stack) <= 1L) {
    list()
  } else {
    trace_state$stack[-length(trace_state$stack)]
  }
  trace_state$counter <- trace_state$counter + 1L
  trace_state$events[[trace_state$counter]] <- data.frame(
    event = trace_state$counter, label = as.character(label),
    elapsed_ms = 1000 * max(0, now() - frame$start),
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

trace_specs <- list(
  lambda = list(
    preparation = c(".prepare_analysis"),
    context_validation = c(".validate_prepared_context"),
    retained_subtree = c(".prepared_tree_subset"),
    optimizer_kernel = c("fast_lambda_tree_optimize_cpp"),
    result_packaging = c(".decorate_fastphylosig_result")
  ),
  D = list(
    preparation = c(".prepare_analysis"),
    context_validation = c(".validate_prepared_context"),
    retained_subtree = c(".prepared_tree_subset"),
    random_null = c(".phylo_d_random_states"),
    brownian_null = c(".phylo_d_brownian_states"),
    observed_and_kernel = c("phylo_d_stream_cpp", "phylo_d_sums_cpp"),
    result_packaging = c(".decorate_fastphylosig_result")
  ),
  Delta = list(
    preparation = c(".prepare_analysis"),
    context_validation = c(".validate_prepared_context"),
    retained_subtree = c(".prepared_tree_subset"),
    observed_ace = c(".fast_ace_lik_anc_workspace"),
    mcmc = c("delta_mcmc_cpp"),
    diagnostics = c(".delta_mcmc_diagnostics"),
    permutation_ace = c(".delta_permutation_worker"),
    result_packaging = c(".decorate_fastphylosig_result")
  )
)

install_traces <- function(method) {
  installed <- list()
  specs <- trace_specs[[method]]
  for (label in names(specs)) {
    for (name in specs[[label]]) {
      if (!exists(name, envir = ns, inherits = FALSE)) next
      ok <- tryCatch({
        trace(
          name,
          tracer = substitute(stage2c_trace_enter(LABEL),
                              list(LABEL = label)),
          exit = substitute(stage2c_trace_exit(LABEL),
                            list(LABEL = label)),
          print = FALSE, where = ns
        )
        TRUE
      }, error = function(e) FALSE)
      if (ok) installed[[length(installed) + 1L]] <- list(
        name = name, where = ns
      )
    }
  }
  installed
}

remove_traces <- function(installed) {
  for (entry in rev(installed)) {
    try(untrace(entry$name, where = entry$where), silent = TRUE)
  }
  invisible(NULL)
}

phase_probe <- function(method, call_fun, seed, shape, n, path, workload,
                        nsim = NA_integer_, ncores = 1L) {
  trace_state$stack <- list()
  trace_state$events <- list()
  trace_state$counter <- 0L
  baseline <- run_timed(call_fun, seed = seed, timeout = timeout_seconds)
  installed <- install_traces(method)
  traced <- tryCatch(
    run_timed(call_fun, seed = seed, timeout = timeout_seconds),
    error = function(e) list(
      elapsed_ms = NA_real_, value = NULL, status = "trace_error",
      error_message = conditionMessage(e), error_class = paste(class(e), collapse = "/"),
      warnings = character()
    )
  )
  remove_traces(installed)
  events <- if (length(trace_state$events)) {
    do.call(rbind, trace_state$events)
  } else {
    data.frame(event = integer(), label = character(), elapsed_ms = numeric(),
               stringsAsFactors = FALSE)
  }
  if (!nrow(events)) {
    return(data.frame(
      method = method, shape = shape, n = as.integer(n), workload = workload,
      path = path, nsim = as.integer(nsim), ncores = as.integer(ncores),
      phase = "<no_trace_event>", phase_elapsed_ms = NA_real_,
      formal_elapsed_ms = baseline$elapsed_ms,
      traced_elapsed_ms = traced$elapsed_ms,
      trace_overhead_pct = if (is.finite(baseline$elapsed_ms) &&
                               baseline$elapsed_ms > 0) {
        100 * (traced$elapsed_ms - baseline$elapsed_ms) / baseline$elapsed_ms
      } else NA_real_,
      installed_traces = length(installed), phase_authoritative = FALSE,
      baseline_status = baseline$status, traced_status = traced$status,
      note = "No trace event was observed; phase accounting unavailable.",
      stringsAsFactors = FALSE
    ))
  }
  overhead <- if (is.finite(baseline$elapsed_ms) && baseline$elapsed_ms > 0) {
    100 * (traced$elapsed_ms - baseline$elapsed_ms) / baseline$elapsed_ms
  } else NA_real_
  out <- do.call(rbind, lapply(split(events, events$label, drop = TRUE), function(z) {
    data.frame(
      method = method, shape = shape, n = as.integer(n), workload = workload,
      path = path, nsim = as.integer(nsim), ncores = as.integer(ncores),
      phase = z$label[[1L]], phase_elapsed_ms = sum(z$elapsed_ms),
      formal_elapsed_ms = baseline$elapsed_ms,
      traced_elapsed_ms = traced$elapsed_ms, trace_overhead_pct = overhead,
      installed_traces = length(installed),
      phase_authoritative = is.finite(overhead) && abs(overhead) <= 5,
      baseline_status = baseline$status, traced_status = traced$status,
      note = "Inclusive trace duration; phases may overlap and are not additive.",
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

timing_rows <- list()
timing_i <- 0L
add_timing <- function(method, shape, n, workload, path, pair, position,
                       nsim, ncores, measured, value = NULL,
                       repeats_used = repeats) {
  timing_i <<- timing_i + 1L
  timing_rows[[timing_i]] <<- data.frame(
    method = method, shape = shape, n = as.integer(n), workload = workload,
    path = path, pair = as.integer(pair), position = as.integer(position),
    repeats = as.integer(repeats_used), nsim = as.integer(nsim),
    ncores = as.integer(ncores), elapsed_ms = as.numeric(measured$elapsed_ms),
    status = as.character(measured$status),
    error_message = as.character(measured$error_message),
    error_class = as.character(measured$error_class),
    warning_count = length(measured$warnings),
    warning_text = warning_text(measured$warnings),
    estimate = estimate_of(value, method),
    stringsAsFactors = FALSE
  )
}

run_pair <- function(method, shape, n, workload, raw_fun, prepared_fun,
                     nsim = NA_integer_, ncores = 1L, include_reference = FALSE,
                     reference_fun = NULL, reference_limit = NULL,
                     repeat_count = NULL,
                     seed_base = 200000L) {
  repeat_count <- if (is.null(repeat_count)) {
    repeats_for(method, workload, n, ncores)
  } else {
    as.integer(repeat_count)
  }
  # Equal warmups are outside formal observations.
  invisible(run_timed(raw_fun, seed = seed_base, timeout = timeout_seconds))
  invisible(run_timed(prepared_fun, seed = seed_base, timeout = timeout_seconds))
  for (pair in seq_len(repeat_count)) {
    order_for_pair <- if (pair %% 2L) c("raw", "prepared") else c("prepared", "raw")
    for (position in seq_along(order_for_pair)) {
      path <- order_for_pair[[position]]
      fun <- if (identical(path, "raw")) raw_fun else prepared_fun
      measured <- run_timed(fun, seed = seed_base + pair,
                            timeout = timeout_seconds)
      add_timing(method, shape, n, workload, path, pair, position,
                 nsim, ncores, measured, measured$value,
                 repeats_used = repeat_count)
    }
    if (include_reference &&
        (is.function(reference_fun) || is.list(reference_limit))) {
      measured <- if (is.list(reference_limit) &&
                      isTRUE(reference_limit$limited)) {
        list(
          elapsed_ms = NA_real_, value = NULL, status = "RESOURCE_LIMIT",
          error_message = reference_limit$reason,
          error_class = "resource_limit", warnings = character()
        )
      } else {
        run_timed(reference_fun, seed = seed_base + pair,
                  timeout = timeout_seconds)
      }
      add_timing(method, shape, n, workload, "reference", pair, NA_integer_,
                 nsim, ncores, measured, measured$value,
                 repeats_used = repeat_count)
    }
  }
}

# Provenance is written before timing so an interrupted run still records the
# execution environment.  The exact commit is supplied by the caller when the
# script is used on a staging library; the fallback is explicit, never guessed.
git_commit <- env_or("FASTPHYLOSIG_STAGE2C_COMMIT", "<unspecified>")
if (identical(git_commit, "<unspecified>")) {
  git_probe <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE,
                                stderr = FALSE), error = function(e) character())
  if (length(git_probe) == 1L && grepl("^[0-9a-f]{40}$", git_probe)) {
    git_commit <- git_probe
  }
}

openmp_capability <- tryCatch({
  val <- capabilities("OpenMP")
  if (length(val) && !is.na(val)) as.character(val) else "<unavailable>"
}, error = function(e) "<unavailable>")
blas <- tryCatch(as.character(extSoftVersion()[["BLAS"]]),
                 error = function(e) "<unavailable>")
compiler <- if (length(R.version$compiler)) as.character(R.version$compiler[[1L]]) else
  "<unavailable>"
hardware <- paste(
  c(
    paste(names(Sys.info()), as.character(Sys.info()), sep = "="),
    paste("PROCESSOR_IDENTIFIER", Sys.getenv("PROCESSOR_IDENTIFIER", unset = "<unset>"), sep = "="),
    paste("NUMBER_OF_PROCESSORS", Sys.getenv("NUMBER_OF_PROCESSORS", unset = "<unset>"), sep = "=")
  ), collapse = " || "
)
provenance <- data.frame(
  stage = "Stage 2C methods",
  package = "fastphylosig",
  package_version = pkg_version,
  source_commit = git_commit,
  R_version = R.version.string,
  R_platform = R.version$platform,
  R_arch = R.version$arch,
  OS = paste(Sys.info()[c("sysname", "release", "version")], collapse = " | "),
  compiler = compiler,
  OpenMP_capability = openmp_capability,
  BLAS = blas,
  LC_ALL = Sys.getenv("LC_ALL", unset = "<unset>"),
  LC_COLLATE = tryCatch(Sys.getlocale("LC_COLLATE"), error = function(e) "<unavailable>"),
  locale_before_benchmark = locale_before,
  hardware = hardware,
  callr_available = callr_available,
  benchmark_serial = TRUE,
  quick_mode = quick_mode,
  selected_methods = paste(selected_methods, collapse = ","),
  lambda_shapes = paste(lambda_shapes, collapse = ","),
  d_shapes = paste(d_shapes, collapse = ","),
  delta_shapes = paste(delta_shapes, collapse = ","),
  timeout_seconds = timeout_seconds,
  repeats = repeats,
  lambda_reference_dense_cap_bytes = reference_dense_cap_bytes,
  lambda_reference_cubic_work_cap = reference_cubic_work_cap,
  lambda_reference_policy = paste(
    "skip before allocation when dense bytes or repeated cubic-work proxy exceeds cap"
  ),
  delta_mcmc_sim = delta_mcmc_sim,
  delta_thin = delta_thin,
  delta_burn = delta_burn,
  stringsAsFactors = FALSE
)
utils::write.csv(provenance, file.path(output_dir, "stage2c_methods_provenance.csv"),
                 row.names = FALSE)

message("[stage2c] serial selected-method profiling starting: ",
        paste(selected_methods, collapse = ","))

# 1. Lambda: raw, prepared compute, and phytools reference.  Preparation is
# separately timed, while prepared end-to-end is reported as prepare + compute
# by the downstream summary consumer; prepared rows themselves are compute only.
if (run_lambda) {
  for (shape in lambda_shapes) {
    for (n in lambda_n) {
      message("[stage2c] lambda shape=", shape, " n=", n)
      tree <- make_fixture(shape, n)
      trait <- make_trait(tree)
      ctx <- fastphylosig::prepare_tree(clone_fixture(tree))
      prep <- run_timed(function() fastphylosig::prepare_tree(clone_fixture(tree)),
                        seed = 10000L + n, timeout = timeout_seconds)
      add_timing("lambda", shape, n, "test_false", "prepare_tree", 0L, NA_integer_,
                 NA_integer_, 1L, prep, prep$value)
      raw_fun <- function() call_lambda("raw", clone_fixture(tree), NULL, trait,
                                        test = FALSE, ncores = lambda_ncores[[1L]])
      prepared_fun <- function() call_lambda("prepared", NULL, ctx, trait,
                                             test = FALSE, ncores = lambda_ncores[[1L]])
      reference_fun <- if (requireNamespace("phytools", quietly = TRUE)) {
        function() phytools::phylosig(tree, trait, method = "lambda", test = FALSE)
      } else NULL
      reference_limit <- reference_precheck(n)
      # The formal lambda table is intentionally limited to raw/prepared and
      # the reference implementation.  test=TRUE remains a guard-only check.
      run_pair("lambda", shape, n, "test_false", raw_fun, prepared_fun,
               nsim = NA_integer_, ncores = lambda_ncores[[1L]],
               include_reference = !is.null(reference_fun),
               reference_fun = reference_fun,
               reference_limit = reference_limit,
               seed_base = 210000L + n)
    }
  }
}

# 2. D: raw/prepared, test=FALSE and random/Brownian null workloads.  D's
# ncores is an OpenMP/kernel control, not a worker-process pool; startup is
# therefore recorded as NA in the phase table with an explicit note.
if (run_d) {
  for (shape in d_shapes) {
    for (n in d_n) {
    message("[stage2c] D shape=", shape, " n=", n)
    tree <- make_fixture(shape, n)
    trait <- make_binary(tree)
    ctx <- fastphylosig::prepare_tree(clone_fixture(tree))
    prep <- run_timed(function() fastphylosig::prepare_tree(clone_fixture(tree)),
                      seed = 30000L + n, timeout = timeout_seconds)
    add_timing("D", shape, n, "test_false", "prepare_tree", 0L, NA_integer_,
               NA_integer_, 1L, prep, prep$value)
      for (ncores in d_ncores) {
        raw_fun <- function() call_d("raw", clone_fixture(tree), NULL, trait,
                                     test = FALSE, nsim = 1000L, ncores = ncores)
        prepared_fun <- function() call_d("prepared", NULL, ctx, trait,
                                          test = FALSE, nsim = 1000L, ncores = ncores)
        run_pair("D", shape, n, "test_false_nsim_1000", raw_fun, prepared_fun,
                 nsim = 1000L, ncores = ncores, seed_base = 310000L + n + ncores)
        for (nsim in d_nsim) {
          raw_test <- function() call_d("raw", clone_fixture(tree), NULL, trait,
                                        test = TRUE, nsim = nsim, ncores = ncores)
          prepared_test <- function() call_d("prepared", NULL, ctx, trait,
                                             test = TRUE, nsim = nsim, ncores = ncores)
          run_pair("D", shape, n, paste0("test_true_nsim_", nsim), raw_test,
                   prepared_test, nsim = nsim, ncores = ncores,
                   seed_base = 320000L + n + nsim + ncores)
        }
      }
    }
  }
}

# 3. Delta: production MCMC controls stay fixed.  test=FALSE is measured only
# serially; test=TRUE covers the requested serial and representative worker
# configurations.  A warning is data in the evidence, not a reason to
# suppress or rerun the call.
if (run_delta) {
  for (shape in delta_shapes) {
    for (n in delta_n) {
      message("[stage2c] Delta shape=", shape, " n=", n)
      tree <- make_fixture(shape, n)
      trait <- make_categorical(tree)
      ctx <- fastphylosig::prepare_tree(clone_fixture(tree))
      prep <- run_timed(function() fastphylosig::prepare_tree(clone_fixture(tree)),
                        seed = 40000L + n, timeout = timeout_seconds)
      add_timing("Delta", shape, n, "test_false", "prepare_tree", 0L, NA_integer_,
                 NA_integer_, 1L, prep, prep$value)

      raw_false <- function() call_delta("raw", clone_fixture(tree), NULL, trait,
                                         test = FALSE, nsim = 19L, ncores = 1L)
      prepared_false <- function() call_delta("prepared", NULL, ctx, trait,
                                              test = FALSE, nsim = 19L, ncores = 1L)
      run_pair("Delta", shape, n, "test_false", raw_false, prepared_false,
               nsim = NA_integer_, ncores = 1L, seed_base = 410000L + n)

      for (ncores in delta_ncores) {
        for (nsim in delta_nsim) {
          raw_test <- function() call_delta("raw", clone_fixture(tree), NULL, trait,
                                            test = TRUE, nsim = nsim, ncores = ncores)
          prepared_test <- function() call_delta("prepared", NULL, ctx, trait,
                                                 test = TRUE, nsim = nsim, ncores = ncores)
          run_pair("Delta", shape, n, paste0("test_true_nsim_", nsim), raw_test,
                   prepared_test, nsim = nsim, ncores = ncores,
                   seed_base = 420000L + n + nsim + ncores)
        }
      }
    }
  }
}

timings <- if (length(timing_rows)) do.call(rbind, timing_rows) else data.frame()
utils::write.csv(timings, file.path(output_dir, "stage2c_methods_timings.csv"),
                 row.names = FALSE)

make_summary <- function(timings) {
  if (!nrow(timings)) return(data.frame())
  key_fields <- timings[c(
    "method", "shape", "n", "workload", "path", "nsim", "ncores"
  )]
  key_fields[] <- lapply(key_fields, function(value) {
    value <- as.character(value)
    value[is.na(value)] <- "<NA>"
    value
  })
  key <- do.call(paste, c(key_fields, sep = "\r"))
  rows <- lapply(split(seq_len(nrow(timings)), key), function(ii) {
    z <- timings[ii, , drop = FALSE]
    finite <- z$elapsed_ms[is.finite(z$elapsed_ms)]
    med <- if (length(finite)) stats::median(finite) else NA_real_
    data.frame(
      method = z$method[[1L]], shape = z$shape[[1L]], n = z$n[[1L]],
      workload = z$workload[[1L]], path = z$path[[1L]], nsim = z$nsim[[1L]],
      ncores = z$ncores[[1L]], repeats_requested = z$repeats[[1L]],
      observations = nrow(z), ok_observations = sum(z$status == "ok"),
      resource_limit_observations = sum(z$status == "RESOURCE_LIMIT"),
      error_observations = sum(z$status == "error"),
      median_ms = med, IQR_ms = if (length(finite)) stats::IQR(finite) else NA_real_,
      traits_per_sec = if (is.finite(med) && med > 0) 1000 / med else NA_real_,
      per_trait_ms = med,
      warning_calls = sum(z$warning_count),
      statuses = paste(sort(unique(z$status)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

summary <- make_summary(timings)
utils::write.csv(summary, file.path(output_dir, "stage2c_methods_summary.csv"),
                 row.names = FALSE)

# Low-overhead phase probes are done after the formal totals.  They use one
# representative n per method and every requested shape, so trace overhead is
# visible but never contaminates the authoritative wall-time table.
phase_rows <- list()
phase_i <- 0L
add_phase <- function(...) {
  phase_i <<- phase_i + 1L
  phase_rows[[phase_i]] <<- list(...)
}

phase_n_lambda <- lambda_n[[which.min(abs(lambda_n - if (quick_mode) 100L else 5000L))]]
phase_n_d <- d_n[[which.min(abs(d_n - if (quick_mode) 100L else 500L))]]
phase_n_delta <- delta_n[[which.min(abs(delta_n - if (quick_mode) 50L else 100L))]]

if (run_lambda) {
  for (shape in lambda_shapes) {
    tree <- make_fixture(shape, phase_n_lambda)
    trait <- make_trait(tree)
    ctx <- fastphylosig::prepare_tree(clone_fixture(tree))
    for (path in c("raw", "prepared")) {
      fun <- if (path == "raw") {
        function() call_lambda("raw", clone_fixture(tree), NULL, trait)
      } else {
        function() call_lambda("prepared", NULL, ctx, trait)
      }
      phase_i <- phase_i + 1L
      phase_rows[[phase_i]] <- phase_probe(
        "lambda", fun, 510000L + phase_i, shape, phase_n_lambda, path,
        "test_false", ncores = lambda_ncores[[1L]]
      )
    }
  }
}

if (run_d) {
  for (shape in d_shapes) {
    tree <- make_fixture(shape, phase_n_d)
    trait <- make_binary(tree)
    ctx <- fastphylosig::prepare_tree(clone_fixture(tree))
    for (ncores in d_ncores) {
      for (path in c("raw", "prepared")) {
        fun <- if (path == "raw") {
          function() call_d("raw", clone_fixture(tree), NULL, trait,
                            test = TRUE, nsim = d_nsim[[1L]], ncores = ncores)
        } else {
          function() call_d("prepared", NULL, ctx, trait,
                            test = TRUE, nsim = d_nsim[[1L]], ncores = ncores)
        }
        phase_i <- phase_i + 1L
        phase_rows[[phase_i]] <- phase_probe(
          "D", fun, 520000L + phase_i, shape, phase_n_d, path,
          paste0("test_true_nsim_", d_nsim[[1L]]),
          nsim = d_nsim[[1L]], ncores = ncores
        )
      }
    }
  }
}

if (run_delta) {
  for (shape in delta_shapes) {
    tree <- make_fixture(shape, phase_n_delta)
    trait <- make_categorical(tree)
    ctx <- fastphylosig::prepare_tree(clone_fixture(tree))
    for (ncores in delta_ncores) {
      for (path in c("raw", "prepared")) {
        fun <- if (path == "raw") {
          function() call_delta("raw", clone_fixture(tree), NULL, trait,
                                test = TRUE, nsim = delta_nsim[[1L]], ncores = ncores)
        } else {
          function() call_delta("prepared", NULL, ctx, trait,
                                test = TRUE, nsim = delta_nsim[[1L]], ncores = ncores)
        }
        phase_i <- phase_i + 1L
        phase_rows[[phase_i]] <- phase_probe(
          "Delta", fun, 530000L + phase_i, shape, phase_n_delta, path,
          paste0("test_true_nsim_", delta_nsim[[1L]]),
          nsim = delta_nsim[[1L]], ncores = ncores
        )
      }
    }
  }
}

phase <- if (length(phase_rows)) do.call(rbind, phase_rows) else data.frame()
utils::write.csv(phase, file.path(output_dir, "stage2c_methods_phase.csv"),
                 row.names = FALSE)

# Direct startup probe is deliberately outside public timing.  D uses OpenMP
# threads and has no PSOCK worker startup; Delta's 2-worker configuration does.
startup_rows <- list()
startup_i <- 0L
startup_probe <- function(method, workers, seed) {
  started <- now()
  status <- "ok"
  error_message <- ""
  if (identical(method, "D") && workers > 1L) {
    return(data.frame(
      method = method, workers = workers, startup_ms = NA_real_, status = "NOT_APPLICABLE",
      error_message = "",
      note = "D ncores is an OpenMP/kernel control; no PSOCK worker startup is used.",
      stringsAsFactors = FALSE
    ))
  }
  if (identical(method, "Delta") && workers > 1L) {
    set.seed(seed)
    cl <- NULL
    on.exit({
      if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
    }, add = TRUE)
    cl <- tryCatch(parallel::makeCluster(workers, type = "PSOCK"),
                   error = function(e) {
                     status <<- "error"
                     error_message <<- conditionMessage(e)
                     NULL
                   })
    if (!is.null(cl)) {
      try(parallel::clusterEvalQ(cl, NULL), silent = TRUE)
      try(parallel::clusterSetRNGStream(cl, seed), silent = TRUE)
    }
  } else {
    # A serial route has no worker startup; retain an explicit zero rather than
    # implying that an unmeasured process launch was included.
    return(data.frame(
      method = method, workers = workers, startup_ms = 0,
      status = "NOT_APPLICABLE",
      error_message = "",
      note = "Serial route has no worker startup.", stringsAsFactors = FALSE
    ))
  }
  data.frame(
    method = method, workers = workers, startup_ms = 1000 * max(0, now() - started),
    status = status, error_message = error_message,
    note = "Direct PSOCK makeCluster/EvalQ/RNG-stream startup probe; not included in public total.",
    stringsAsFactors = FALSE
  )
}
if (run_d) {
  for (workers in d_ncores) {
    startup_i <- startup_i + 1L
    startup_rows[[startup_i]] <- startup_probe("D", workers, 610000L + workers)
  }
}
if (run_delta) {
  for (workers in delta_ncores) {
    startup_i <- startup_i + 1L
    startup_rows[[startup_i]] <- startup_probe("Delta", workers, 620000L + workers)
  }
}
startup <- if (length(startup_rows)) {
  do.call(rbind, startup_rows)
} else {
  data.frame()
}
utils::write.csv(startup, file.path(output_dir, "stage2c_methods_startup.csv"),
                 row.names = FALSE)

# Correctness guards use small fixed fixtures and controlled inputs.  They are
# not a replacement for the package test suite; they stop interpretation of a
# timing run if raw/prepared parity or warning/status contracts fail.
guard_rows <- list()
guard_i <- 0L
add_guard <- function(method, shape, n, workload, raw_run, prepared_run,
                      seed, fields, controlled = FALSE) {
  raw <- run_timed(raw_run, seed = seed, timeout = timeout_seconds)
  prepared <- run_timed(prepared_run, seed = seed, timeout = timeout_seconds)
  raw_values <- vapply(fields, function(x) field_of(raw$value, x), numeric(1))
  prepared_values <- vapply(fields, function(x) field_of(prepared$value, x), numeric(1))
  finite <- is.finite(raw_values) & is.finite(prepared_values)
  diffs <- abs(raw_values - prepared_values)
  max_abs <- if (any(finite)) max(diffs[finite]) else NA_real_
  guard_i <<- guard_i + 1L
  guard_rows[[guard_i]] <<- data.frame(
    method = method, shape = shape, n = n, workload = workload,
    controlled_inputs = controlled,
    raw_status = raw$status, prepared_status = prepared$status,
    status_equal = identical(raw$status, prepared$status),
    raw_result_status = status_of(raw$value),
    prepared_result_status = status_of(prepared$value),
    result_status_equal = identical(status_of(raw$value), status_of(prepared$value)),
    raw_warning_text = warning_text(raw$warnings),
    prepared_warning_text = warning_text(prepared$warnings),
    warning_equal = identical(raw$warnings, prepared$warnings),
    max_abs_difference = max_abs,
    numerical_equal_1e_10 = is.finite(max_abs) && max_abs <= 1e-10,
    guard_pass = identical(raw$status, prepared$status) &&
      identical(status_of(raw$value), status_of(prepared$value)) &&
      identical(raw$warnings, prepared$warnings) &&
      is.finite(max_abs) && max_abs <= 1e-10,
    stringsAsFactors = FALSE
  )
}

guard_shape <- "random"
# Use the same representative n=50 fixture in quick and formal guard runs.
# With n=20 Delta can legitimately enter its small-retained-sample failure
# path, which would test a diagnostic boundary rather than raw/prepared parity.
guard_tree <- make_fixture(guard_shape, 50L)
guard_n <- length(guard_tree$tip.label)
guard_ctx <- fastphylosig::prepare_tree(clone_fixture(guard_tree))
guard_trait <- make_trait(guard_tree)
if (run_lambda) {
  add_guard(
    "lambda", guard_shape, guard_n, "test_true",
    function() call_lambda("raw", clone_fixture(guard_tree), NULL, guard_trait, TRUE),
    function() call_lambda("prepared", NULL, guard_ctx, guard_trait, TRUE),
    seed = 710001L, fields = list(c("lambda", "lambda_fast"), c("logL", "logLik"),
                                  c("LR"), c("P"))
  )
}

guard_binary <- make_binary(guard_tree)
guard_rand <- matrix(rep(seq_len(19L) %% 2L, each = guard_n),
                     nrow = guard_n, ncol = 19L)
guard_brown <- matrix(rep(seq(0.25, 0.75, length.out = guard_n), 19L),
                      nrow = guard_n, ncol = 19L)
if (run_d) {
  add_guard(
    "D", guard_shape, guard_n, "controlled_null_nsim_19",
    function() call_d("raw", clone_fixture(guard_tree), NULL, guard_binary,
                      TRUE, 19L, 1L, guard_rand, guard_brown),
    function() call_d("prepared", NULL, guard_ctx, guard_binary,
                      TRUE, 19L, 1L, guard_rand, guard_brown),
    seed = 710002L, controlled = TRUE,
    # Vector fast_d() returns the caper-compatible DEstimate/Pval1/Pval0 names;
    # the data-frame route uses D_fast/Pval1_fast/Pval0_fast.
    fields = list(c("DEstimate", "D_fast", "D"),
                  c("Pval1", "Pval1_fast", "P_random"),
                  c("Pval0", "Pval0_fast", "P_Brownian"))
  )
}

# Use a fixed-seed random two-state composition for the Delta parity guard.
# A deterministic alternating composition can be a legitimate difficult ACE
# likelihood surface and is not the representative production fixture used by
# the Stage 2A Delta evidence.
if (run_delta) {
  set.seed(710004L)
  guard_cat <- stats::setNames(
    sample(c("a", "b"), guard_n, replace = TRUE), guard_tree$tip.label
  )
  if (length(unique(guard_cat)) < 2L) guard_cat[[1L]] <- "b"
  add_guard(
    "Delta", guard_shape, guard_n, "fixed_seed_nsim_19",
    function() call_delta("raw", clone_fixture(guard_tree), NULL, guard_cat,
                          TRUE, 19L, 1L),
    function() call_delta("prepared", NULL, guard_ctx, guard_cat,
                          TRUE, 19L, 1L),
    seed = 710003L,
    fields = list(c("Delta_fast", "delta"), c("P_fast", "P"),
                  c("alpha_mean"), c("beta_mean")), controlled = TRUE
  )
}

guards <- if (length(guard_rows)) do.call(rbind, guard_rows) else data.frame()
utils::write.csv(guards, file.path(output_dir, "stage2c_methods_correctness.csv"),
                 row.names = FALSE)

# A concise machine-readable run status helps the parent review distinguish an
# interrupted/resource-limited grid from a completed one without implying that
# a small Delta chain is suitable for biological inference.
run_status <- data.frame(
  stage = "Stage 2C methods",
  quick_mode = quick_mode,
  formal_rows = nrow(timings),
  formal_ok_rows = if (nrow(timings)) sum(timings$status == "ok") else 0L,
  resource_limit_rows = if (nrow(timings)) sum(timings$status == "RESOURCE_LIMIT") else 0L,
  correctness_guard_rows = nrow(guards),
  correctness_guard_pass = if (nrow(guards)) sum(guards$guard_pass) else 0L,
  correctness_all_pass = nrow(guards) > 0L && all(guards$guard_pass),
  production_code_changed = FALSE,
  benchmark_processes_concurrent = FALSE,
  note = "Phase trace rows are inclusive and non-authoritative when overhead exceeds 5 percent.",
  stringsAsFactors = FALSE
)
utils::write.csv(run_status, file.path(output_dir, "stage2c_methods_status.csv"),
                 row.names = FALSE)

message("[stage2c] complete: ", output_dir)
