#!/usr/bin/env Rscript

# fastphylosig 0.2.0 Stage 2B2B Candidate 1 paired benchmark.
#
# This harness is read-only with respect to the package source. It compares
# the frozen Stage 2B2B oracle with the current inspector on fixed tree
# fixtures, and measures both direct inspection and prepare_tree(). The old
# inspector is installed in the package namespace before each timer starts;
# namespace setup and restoration are outside the timed interval.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(index, default) {
  if (length(args) >= index && nzchar(args[[index]])) args[[index]] else default
}

repo <- normalizePath(arg_value(1L, getwd()), winslash = "/", mustWork = TRUE)
out_dir <- arg_value(2L, file.path(tempdir(), "fastphylosig-stage2b2b1"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")
options(stringsAsFactors = FALSE)

prefix <- "FASTPHYLOSIG_STAGE2B2B1_"
quick_mode <- identical(Sys.getenv(paste0(prefix, "QUICK"), "0"), "1")
shape_names <- c("balanced", "random", "pectinate")
n_grid <- if (quick_mode) c(20L, 50L, 100L) else
  c(500L, 1000L, 2000L, 5000L, 10000L, 20000L)

parse_csv <- function(text, what) {
  values <- trimws(strsplit(text, ",", fixed = TRUE)[[1L]])
  values <- values[nzchar(values)]
  if (!length(values)) stop(sprintf("%s must not be empty.", what), call. = FALSE)
  values
}

shape_override <- Sys.getenv(paste0(prefix, "SHAPES"), "")
if (nzchar(shape_override)) shape_names <- unique(parse_csv(
  shape_override, paste0(prefix, "SHAPES")
))
if (!length(shape_names) || any(!shape_names %in% c("balanced", "random", "pectinate"))) {
  stop("SHAPES must contain only balanced, random, and/or pectinate.", call. = FALSE)
}

n_override <- Sys.getenv(paste0(prefix, "N_GRID"), "")
if (nzchar(n_override)) {
  n_grid <- suppressWarnings(as.integer(parse_csv(
    n_override, paste0(prefix, "N_GRID")
  )))
  if (anyNA(n_grid) || any(n_grid < 2L)) {
    stop("N_GRID must contain integers >= 2.", call. = FALSE)
  }
}
n_grid <- sort(unique(n_grid))

parse_positive_integer <- function(text, what) {
  value <- suppressWarnings(as.integer(text))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(sprintf("%s must be a positive integer.", what), call. = FALSE)
  }
  value
}

repeats <- parse_positive_integer(
  Sys.getenv(paste0(prefix, "REPEATS"), if (quick_mode) "2" else "10"),
  paste0(prefix, "REPEATS")
)
repeat_schedule_text <- Sys.getenv(paste0(prefix, "REPEATS_BY_N"), "")
repeat_schedule <- integer()
if (nzchar(repeat_schedule_text)) {
  schedule_parts <- parse_csv(
    repeat_schedule_text, paste0(prefix, "REPEATS_BY_N")
  )
  schedule_pairs <- strsplit(schedule_parts, ":", fixed = TRUE)
  if (any(lengths(schedule_pairs) != 2L)) {
    stop("REPEATS_BY_N must use n:repeats pairs.", call. = FALSE)
  }
  schedule_n <- suppressWarnings(vapply(schedule_pairs, function(x) {
    as.integer(x[[1L]])
  }, integer(1)))
  schedule_repeats <- suppressWarnings(vapply(schedule_pairs, function(x) {
    as.integer(x[[2L]])
  }, integer(1)))
  if (anyNA(schedule_n) || any(schedule_n < 2L) ||
      anyNA(schedule_repeats) || any(schedule_repeats < 1L)) {
    stop("REPEATS_BY_N contains invalid n:repeats values.", call. = FALSE)
  }
  repeat_schedule <- schedule_repeats
  names(repeat_schedule) <- as.character(schedule_n)
}

repeats_for_n <- function(n) {
  key <- as.character(as.integer(n))
  value <- if (length(repeat_schedule) && key %in% names(repeat_schedule)) {
    repeat_schedule[[key]]
  } else repeats
  required <- if (as.integer(n) <= 5000L) 10L else 5L
  if (!quick_mode && value < required) {
    stop(sprintf(
      "n=%d requires at least %d paired repeats (set QUICK=1 for smoke runs).",
      n, required
    ), call. = FALSE)
  }
  as.integer(value)
}

timing_timeout <- suppressWarnings(as.numeric(
  Sys.getenv(paste0(prefix, "TIMING_TIMEOUT"), "0")
))
if (!is.finite(timing_timeout) || timing_timeout <= 0) timing_timeout <- Inf

partition_label <- Sys.getenv(paste0(prefix, "PARTITION"), "full")

# Load the candidate source from the requested checkout, or an explicitly
# supplied installed library in an ASCII staging copy.
installed_lib <- Sys.getenv(paste0(prefix, "LIB"), "")
if (nzchar(installed_lib) && dir.exists(file.path(installed_lib, "fastphylosig"))) {
  suppressPackageStartupMessages(
    library("fastphylosig", lib.loc = installed_lib, character.only = TRUE)
  )
  package_load <- paste0("installed library: ", normalizePath(
    installed_lib, winslash = "/", mustWork = TRUE
  ))
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  suppressPackageStartupMessages(pkgload::load_all(repo, quiet = TRUE))
  package_load <- "pkgload::load_all(source checkout)"
} else {
  stop("pkgload is required unless FASTPHYLOSIG_STAGE2B2B1_LIB is set.",
       call. = FALSE)
}
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("ape is required for fixed tree fixtures.", call. = FALSE)
}
suppressPackageStartupMessages(library(ape))

ns <- asNamespace("fastphylosig")
get_internal <- function(name) get(name, envir = ns, inherits = FALSE)
current_inspector <- get_internal(".inspect_tree_core")
signals <- c("K", "lambda", "D", "Delta")

oracle_file <- file.path(
  repo, "tests", "testthat", "fixtures", "stage2b2b-old-inspect-tree.R"
)
if (!file.exists(oracle_file)) {
  stop("frozen Stage 2B2B oracle fixture is missing.", call. = FALSE)
}
oracle_env <- new.env(parent = ns)
sys.source(oracle_file, envir = oracle_env)
old_inspector <- get(".stage2b2b_old_inspect_tree_core", envir = oracle_env)

elapsed_seconds <- function() unname(proc.time()[["elapsed"]])

capture_call <- function(fun) {
  warning_count <- 0L
  error_message <- ""
  error_class <- ""
  status <- "ok"
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(condition) {
        warning_count <<- warning_count + 1L
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      status <<- if (inherits(condition, "elapsedTimeLimit")) {
        "censored"
      } else {
        "error"
      }
      error_message <<- conditionMessage(condition)
      error_class <<- paste(class(condition), collapse = "/")
      NULL
    }
  )
  list(value = value, status = status, warning_count = warning_count,
       error_message = error_message, error_class = error_class)
}

time_call <- function(fun, timeout_seconds = Inf) {
  limited <- is.finite(timeout_seconds) && timeout_seconds > 0
  if (limited) {
    setTimeLimit(cpu = Inf, elapsed = timeout_seconds, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
  }
  start <- elapsed_seconds()
  captured <- capture_call(fun)
  elapsed <- max(0, elapsed_seconds() - start)
  if (limited && identical(captured$status, "ok") &&
      elapsed > timeout_seconds) {
    captured$status <- "censored_overrun"
    captured$error_message <- sprintf(
      "elapsed time %.3fs exceeded %.3fs timing limit.",
      elapsed, timeout_seconds
    )
    captured$error_class <- "timing_overrun"
  }
  captured$elapsed_ms <- 1000 * elapsed
  captured
}

clone_fixture <- function(tree) unserialize(serialize(tree, NULL, version = 2L))

# Install one inspector before timing and restore its original binding after
# timing. The helper is idempotent so error and timeout paths are covered.
install_inspector <- function(implementation) {
  replacement <- switch(
    implementation,
    old = old_inspector,
    new = current_inspector,
    stop("unknown inspector implementation", call. = FALSE)
  )
  was_locked <- bindingIsLocked(".inspect_tree_core", ns)
  previous <- get(".inspect_tree_core", envir = ns, inherits = FALSE)
  if (was_locked) unlockBinding(".inspect_tree_core", ns)
  assign(".inspect_tree_core", replacement, envir = ns)
  restored <- FALSE
  function() {
    if (isTRUE(restored)) return(invisible(NULL))
    if (bindingIsLocked(".inspect_tree_core", ns)) {
      unlockBinding(".inspect_tree_core", ns)
    }
    assign(".inspect_tree_core", previous, envir = ns)
    if (was_locked) lockBinding(".inspect_tree_core", ns)
    restored <<- TRUE
    invisible(NULL)
  }
}

run_route <- function(tree, implementation, workload, timed = FALSE) {
  restore <- install_inspector(implementation)
  on.exit(restore(), add = TRUE)
  call <- function() {
    inspector <- get(".inspect_tree_core", envir = ns, inherits = FALSE)
    if (identical(workload, "direct_inspection")) {
      inspector(tree, signal = signals)
    } else if (identical(workload, "prepare_tree")) {
      fastphylosig::prepare_tree(tree)
    } else {
      stop("unknown benchmark workload", call. = FALSE)
    }
  }
  result <- if (timed) time_call(call, timing_timeout) else capture_call(call)
  restore()
  result
}

timed_record <- function(tree, implementation, workload) {
  # Clone before namespace setup and before the timer starts. The returned
  # record deliberately drops the potentially large prepared context.
  result <- run_route(clone_fixture(tree), implementation, workload, timed = TRUE)
  record <- list(
    status = result$status,
    elapsed_ms = result$elapsed_ms,
    warning_count = result$warning_count,
    error_class = result$error_class,
    error_message = result$error_message
  )
  result$value <- NULL
  rm(result)
  record
}

compare_inspection <- function(old, new, label) {
  if (!identical(old$status, "ok") || !identical(new$status, "ok")) {
    stop(sprintf("%s warmup failed: old=%s new=%s.", label,
                 old$status, new$status), call. = FALSE)
  }
  if (!identical(old$value, new$value)) {
    stop(sprintf("%s old/new inspection result differs.", label), call. = FALSE)
  }
  invisible(NULL)
}

compare_prepared <- function(old, new, label) {
  if (!identical(old$status, "ok") || !identical(new$status, "ok")) {
    stop(sprintf("%s warmup failed: old=%s new=%s.", label,
                 old$status, new$status), call. = FALSE)
  }
  fields <- c(
    "tree", "tip.label", "n_tip", "cache_budget", "max_cached_subsets",
    "fingerprint", "inspection", "canonical_mapping", "canonical_summary",
    "generic_summary", "method_capability", "structural_entry_validation"
  )
  for (field in fields) {
    if (!identical(old$value[[field]], new$value[[field]])) {
      stop(sprintf("%s old/new prepared field differs: %s.", label, field),
           call. = FALSE)
    }
  }
  invisible(NULL)
}

# -------------------------------------------------------------------------
# Deterministic fixed fixtures.

tip_labels <- function(n) paste0("sp", seq_len(n))

make_balanced <- function(n) {
  n <- as.integer(n)
  if (n < 2L) stop("n must be at least two.", call. = FALSE)
  n_internal <- n - 1L
  max_edges <- 2L * n_internal
  edges <- matrix(0L, nrow = max_edges, ncol = 2L)
  queue_node <- integer(n_internal)
  queue_lo <- integer(n_internal)
  queue_hi <- integer(n_internal)
  root <- n + 1L
  next_internal <- root + 1L
  head <- 1L
  tail <- 1L
  queue_node[[1L]] <- root
  queue_lo[[1L]] <- 1L
  queue_hi[[1L]] <- n
  edge_i <- 0L
  while (head <= tail) {
    node <- queue_node[[head]]
    lo <- queue_lo[[head]]
    hi <- queue_hi[[head]]
    head <- head + 1L
    mid <- floor((lo + hi) / 2L)
    child_lo <- c(lo, mid + 1L)
    child_hi <- c(mid, hi)
    for (j in 1:2) {
      child <- if (child_lo[[j]] == child_hi[[j]]) {
        child_lo[[j]]
      } else {
        child <- next_internal
        next_internal <- next_internal + 1L
        tail <- tail + 1L
        queue_node[[tail]] <- child
        queue_lo[[tail]] <- child_lo[[j]]
        queue_hi[[tail]] <- child_hi[[j]]
        child
      }
      edge_i <- edge_i + 1L
      edges[edge_i, ] <- c(node, child)
    }
  }
  tree <- list(
    edge = edges[seq_len(edge_i), , drop = FALSE],
    tip.label = tip_labels(n),
    edge.length = rep(1, edge_i),
    Nnode = n_internal,
    root = root
  )
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

make_pectinate <- function(n) {
  n <- as.integer(n)
  if (n < 2L) stop("n must be at least two.", call. = FALSE)
  root <- n + 1L
  edges <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  if (n == 2L) {
    edges[1L, ] <- c(root, 1L)
    edges[2L, ] <- c(root, 2L)
  } else {
    next_internal <- root + 1L
    current <- next_internal
    row <- 1L
    edges[row, ] <- c(root, 1L)
    row <- row + 1L
    edges[row, ] <- c(root, current)
    row <- row + 1L
    next_internal <- next_internal + 1L
    if (n > 3L) {
      for (tip in 2:(n - 2L)) {
        edges[row, ] <- c(current, tip)
        row <- row + 1L
        edges[row, ] <- c(current, next_internal)
        row <- row + 1L
        current <- next_internal
        next_internal <- next_internal + 1L
      }
    }
    edges[row, ] <- c(current, n - 1L)
    row <- row + 1L
    edges[row, ] <- c(current, n)
  }
  tree <- list(
    edge = edges, tip.label = tip_labels(n),
    edge.length = rep(1, nrow(edges)), Nnode = n - 1L, root = root
  )
  class(tree) <- "phylo"
  tree
}

make_random <- function(n) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(2072000L + as.integer(n))
  tree <- ape::rtree(n, tip.label = tip_labels(n))
  tree <- ape::reorder.phylo(tree, order = "postorder")
  tree$edge.length <- rep(1, nrow(tree$edge))
  tree
}

make_fixture <- function(shape, n) {
  switch(shape,
         balanced = make_balanced(n),
         random = make_random(n),
         pectinate = make_pectinate(n),
         stop("unknown fixture shape", call. = FALSE))
}

validate_fixture <- function(tree, shape, n) {
  n_total <- n + n - 1L
  if (!inherits(tree, "phylo") || length(tree$tip.label) != n ||
      as.integer(tree$Nnode) != n - 1L || nrow(tree$edge) != 2L * (n - 1L) ||
      length(tree$edge.length) != nrow(tree$edge) ||
      any(!is.finite(tree$edge.length)) || any(tree$edge.length <= 0) ||
      anyNA(tree$tip.label) || anyDuplicated(tree$tip.label)) {
    stop(sprintf("invalid %s n=%d fixture.", shape, n), call. = FALSE)
  }
  if (max(tree$edge) > n_total || min(tree$edge) < 1L) {
    stop(sprintf("out-of-range %s n=%d fixture.", shape, n), call. = FALSE)
  }
  check <- current_inspector(tree, signal = signals)
  if (!isTRUE(check$tree_summary$valid)) {
    stop(sprintf("%s n=%d fixture failed candidate inspection.", shape, n),
         call. = FALSE)
  }
  invisible(tree)
}

root_of <- function(tree) {
  edge <- tree$edge
  roots <- setdiff(unique(edge[, 1L]), unique(edge[, 2L]))
  if (length(roots) == 1L) as.integer(roots[[1L]]) else NA_integer_
}

fixture_checksum <- function(tree) {
  edge <- as.double(tree$edge)
  # This checksum is a compact fixture identity aid; the construction code
  # and deterministic random seed remain the reproducibility source.
  sum(edge * seq_along(edge))
}

# -------------------------------------------------------------------------
# Source-aligned operation counts. These are deterministic replays of the
# valid fixed fixtures and are never collected inside a timing interval.

old_stack_metrics <- function(tree) {
  edge <- tree$edge
  n_tip <- length(tree$tip.label)
  n_total <- n_tip + as.integer(tree$Nnode)
  root <- root_of(tree)
  children <- split(as.integer(edge[, 2L]), as.integer(edge[, 1L]))
  stack <- integer(n_total)
  old_character_lookup_count <- 0
  old_character_key_materialization_count <- 0
  old_stack_pop_calls <- 0
  old_stack_append_calls <- 0
  old_stack_pop_copied_slots <- 0
  old_stack_append_copied_slots <- 0
  for (pass in 1:2) {
    top <- 1L
    stack[[top]] <- root
    while (top > 0L) {
      stack_length <- top
      node <- stack[[top]]
      top <- top - 1L
      old_stack_pop_calls <- old_stack_pop_calls + 1
      old_stack_pop_copied_slots <- old_stack_pop_copied_slots +
        max(0, stack_length - 1L)
      old_character_lookup_count <- old_character_lookup_count + 1
      old_character_key_materialization_count <-
        old_character_key_materialization_count + 1
      kids <- children[[as.character(node)]]
      if (is.null(kids)) next
      for (kid in kids) {
        old_stack_append_calls <- old_stack_append_calls + 1
        # c(stack, kid) copies the active pre-append stack slots. The new
        # element itself is counted separately as an append call.
        old_stack_append_copied_slots <- old_stack_append_copied_slots + top
        top <- top + 1L
        stack[[top]] <- as.integer(kid)
      }
    }
  }
  data.frame(
    old_character_lookup_count = old_character_lookup_count,
    old_character_key_materialization_count = old_character_key_materialization_count,
    old_stack_pop_calls = old_stack_pop_calls,
    old_stack_append_calls = old_stack_append_calls,
    old_stack_pop_copied_slots = old_stack_pop_copied_slots,
    old_stack_append_copied_slots = old_stack_append_copied_slots,
    old_stack_copied_slots_total = old_stack_pop_copied_slots +
      old_stack_append_copied_slots,
    stringsAsFactors = FALSE
  )
}

new_integer_metrics <- function(tree) {
  edge <- tree$edge
  parent <- as.integer(edge[, 1L])
  child <- as.integer(edge[, 2L])
  n_tip <- length(tree$tip.label)
  n_total <- n_tip + as.integer(tree$Nnode)
  n_edge <- length(child)
  child_counts <- tabulate(parent, nbins = n_total)
  child_offsets <- integer(n_total + 1L)
  child_offsets[-1L] <- cumsum(child_counts)
  adjacency_children <- integer(n_edge)
  adjacency_edge_rows <- integer(n_edge)
  adjacency_cursor <- child_offsets[seq_len(n_total)] + 1L
  fill_cursor_reads <- 0
  fill_cursor_updates <- 0
  fill_writes <- 0
  for (row in seq_along(child)) {
    parent_node <- parent[[row]]
    position <- adjacency_cursor[[parent_node]]
    fill_cursor_reads <- fill_cursor_reads + 1
    adjacency_children[[position]] <- child[[row]]
    adjacency_edge_rows[[position]] <- row
    fill_writes <- fill_writes + 2
    adjacency_cursor[[parent_node]] <- position + 1L
    fill_cursor_updates <- fill_cursor_updates + 1
  }

  stack <- integer(n_total)
  top <- 0L
  offset_reads <- 0
  child_reads <- 0
  edge_row_reads <- 0
  stack_pop_calls <- 0
  stack_push_writes <- 0
  root <- root_of(tree)
  for (pass in 1:2) {
    top <- 1L
    stack[[top]] <- root
    stack_push_writes <- stack_push_writes + 1
    while (top > 0L) {
      node <- stack[[top]]
      top <- top - 1L
      stack_pop_calls <- stack_pop_calls + 1
      start <- child_offsets[[node]] + 1L
      end <- child_offsets[[node + 1L]]
      offset_reads <- offset_reads + 2
      if (start > end) next
      for (position in seq.int(start, end)) {
        child_reads <- child_reads + 1
        if (pass == 2L) edge_row_reads <- edge_row_reads + 1
        top <- top + 1L
        stack[[top]] <- adjacency_children[[position]]
        stack_push_writes <- stack_push_writes + 1
      }
    }
  }
  data.frame(
    new_integer_adjacency_fill_writes = fill_writes,
    new_integer_cursor_reads = fill_cursor_reads,
    new_integer_cursor_updates = fill_cursor_updates,
    new_integer_adjacency_child_reads = child_reads,
    new_integer_adjacency_edge_row_reads = edge_row_reads,
    new_integer_offset_reads = offset_reads,
    new_stack_pop_calls = stack_pop_calls,
    new_stack_push_writes = stack_push_writes,
    new_stack_resize_copied_slots = 0,
    stringsAsFactors = FALSE
  )
}

operation_counts <- function(tree, shape, n) {
  old <- old_stack_metrics(tree)
  new <- new_integer_metrics(tree)
  base <- data.frame(
    shape = shape, n = as.integer(n),
    n_tip = length(tree$tip.label),
    n_internal = as.integer(tree$Nnode),
    n_total = length(tree$tip.label) + as.integer(tree$Nnode),
    n_edge = nrow(tree$edge),
    stringsAsFactors = FALSE
  )
  cbind(base, old, new,
        operation_count_method = "source-aligned deterministic replay",
        stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------------
# Build fixtures, warm up, check parity, and run serial alternating pairs.

stage_message <- function(text) message("[stage2b2b1] ", text)
stage_message(sprintf(
  "starting; partition=%s; n=%s; shapes=%s; quick=%s",
  partition_label, paste(n_grid, collapse = ","),
  paste(shape_names, collapse = ","), quick_mode
))

fixture_cache <- new.env(parent = emptyenv())
fixture_rows <- list()
operation_rows <- list()
fixture_id <- 0L
for (shape in shape_names) {
  for (n in n_grid) {
    key <- paste(shape, n, sep = "::")
    stage_message(sprintf("constructing %s n=%d", shape, n))
    tree <- make_fixture(shape, n)
    validate_fixture(tree, shape, n)
    fixture_cache[[key]] <- tree
    fixture_id <- fixture_id + 1L
    fixture_rows[[fixture_id]] <- data.frame(
      shape = shape, n = as.integer(n), root = root_of(tree),
      n_tip = length(tree$tip.label), n_internal = as.integer(tree$Nnode),
      n_edge = nrow(tree$edge), edge_checksum = fixture_checksum(tree),
      fixed_branch_length = TRUE, valid = TRUE, stringsAsFactors = FALSE
    )
    operation_rows[[fixture_id]] <- operation_counts(tree, shape, n)
  }
}

timing_rows <- list()
timing_id <- 0L
workloads <- c("direct_inspection", "prepare_tree")
for (shape in shape_names) {
  for (n in n_grid) {
    tree <- fixture_cache[[paste(shape, n, sep = "::")]]
    reps <- repeats_for_n(n)
    stage_message(sprintf("warming %s n=%d; pairs=%d", shape, n, reps))

    warm_old_inspection <- run_route(
      clone_fixture(tree), "old", "direct_inspection", timed = FALSE
    )
    warm_new_inspection <- run_route(
      clone_fixture(tree), "new", "direct_inspection", timed = FALSE
    )
    compare_inspection(
      warm_old_inspection, warm_new_inspection,
      sprintf("%s n=%d direct inspection", shape, n)
    )
    warm_old_prepare <- run_route(
      clone_fixture(tree), "old", "prepare_tree", timed = FALSE
    )
    warm_new_prepare <- run_route(
      clone_fixture(tree), "new", "prepare_tree", timed = FALSE
    )
    compare_prepared(
      warm_old_prepare, warm_new_prepare,
      sprintf("%s n=%d prepare_tree", shape, n)
    )
    warm_old_inspection$value <- NULL
    warm_new_inspection$value <- NULL
    warm_old_prepare$value <- NULL
    warm_new_prepare$value <- NULL
    rm(warm_old_inspection, warm_new_inspection,
       warm_old_prepare, warm_new_prepare)
    gc(FALSE)

    for (pair in seq_len(reps)) {
      order <- if (pair %% 2L) c("old", "new") else c("new", "old")
      for (workload in workloads) {
        records <- vector("list", 2L)
        names(records) <- order
        for (implementation in order) {
          records[[implementation]] <- timed_record(
            tree, implementation, workload
          )
        }
        old_record <- records[["old"]]
        new_record <- records[["new"]]
        old_ok <- identical(old_record$status, "ok") &&
          is.finite(old_record$elapsed_ms)
        new_ok <- identical(new_record$status, "ok") &&
          is.finite(new_record$elapsed_ms)
        pair_ok <- old_ok && new_ok
        ratio <- if (pair_ok && old_record$elapsed_ms > 0) {
          new_record$elapsed_ms / old_record$elapsed_ms
        } else NA_real_
        slowdown <- if (is.finite(ratio)) 100 * (ratio - 1) else NA_real_
        timing_id <- timing_id + 1L
        timing_rows[[timing_id]] <- data.frame(
          shape = shape, n = as.integer(n), pair = as.integer(pair),
          order = paste(order, collapse = ">"), workload = workload,
          old_status = old_record$status, new_status = new_record$status,
          old_elapsed_ms = old_record$elapsed_ms,
          new_elapsed_ms = new_record$elapsed_ms,
          old_warning_count = old_record$warning_count,
          new_warning_count = new_record$warning_count,
          old_error_class = old_record$error_class,
          new_error_class = new_record$error_class,
          old_error_message = old_record$error_message,
          new_error_message = new_record$error_message,
          pair_ok = pair_ok, new_over_old = ratio,
          slowdown_pct = slowdown, stringsAsFactors = FALSE
        )
        rm(records, old_record, new_record)
        gc(FALSE)
      }
    }
  }
}

bind_rows <- function(rows, empty) {
  if (length(rows)) do.call(rbind, rows) else empty[FALSE, , drop = FALSE]
}

timings <- bind_rows(timing_rows, data.frame(
  shape = character(), n = integer(), pair = integer(), order = character(),
  workload = character(), old_status = character(), new_status = character(),
  old_elapsed_ms = numeric(), new_elapsed_ms = numeric(),
  old_warning_count = integer(), new_warning_count = integer(),
  old_error_class = character(), new_error_class = character(),
  old_error_message = character(), new_error_message = character(),
  pair_ok = logical(), new_over_old = numeric(), slowdown_pct = numeric(),
  stringsAsFactors = FALSE
))
fixtures <- bind_rows(fixture_rows, data.frame(
  shape = character(), n = integer(), root = integer(), n_tip = integer(),
  n_internal = integer(), n_edge = integer(), edge_checksum = numeric(),
  fixed_branch_length = logical(), valid = logical(), stringsAsFactors = FALSE
))
operations <- bind_rows(operation_rows, data.frame(
  shape = character(), n = integer(), n_tip = integer(), n_internal = integer(),
  n_total = integer(), n_edge = integer(),
  old_character_lookup_count = numeric(),
  old_character_key_materialization_count = numeric(),
  old_stack_pop_calls = numeric(), old_stack_append_calls = numeric(),
  old_stack_pop_copied_slots = numeric(), old_stack_append_copied_slots = numeric(),
  old_stack_copied_slots_total = numeric(),
  new_integer_adjacency_fill_writes = numeric(),
  new_integer_cursor_reads = numeric(), new_integer_cursor_updates = numeric(),
  new_integer_adjacency_child_reads = numeric(),
  new_integer_adjacency_edge_row_reads = numeric(),
  new_integer_offset_reads = numeric(), new_stack_pop_calls = numeric(),
  new_stack_push_writes = numeric(), new_stack_resize_copied_slots = numeric(),
  operation_count_method = character(), stringsAsFactors = FALSE
))

summarise_timings <- function(data) {
  rows <- list()
  row_id <- 0L
  for (shape in shape_names) {
    for (n in n_grid) {
      for (workload in workloads) {
        z <- data[data$shape == shape & data$n == n &
                    data$workload == workload, , drop = FALSE]
        ok <- z$pair_ok %in% TRUE & is.finite(z$old_elapsed_ms) &
          is.finite(z$new_elapsed_ms)
        old <- z$old_elapsed_ms[ok]
        new <- z$new_elapsed_ms[ok]
        row_id <- row_id + 1L
        before_median <- if (length(old)) stats::median(old) else NA_real_
        after_median <- if (length(new)) stats::median(new) else NA_real_
        ratio <- if (is.finite(before_median) && before_median > 0 &&
                     is.finite(after_median)) after_median / before_median else NA_real_
        rows[[row_id]] <- data.frame(
          shape = shape, n = as.integer(n), workload = workload,
          repeats = nrow(z), successful_pairs = length(old),
          failed_pairs = nrow(z) - length(old),
          before_median_ms = before_median,
          before_iqr_ms = if (length(old)) stats::IQR(old) else NA_real_,
          after_median_ms = after_median,
          after_iqr_ms = if (length(new)) stats::IQR(new) else NA_real_,
          speedup_before_over_after = if (is.finite(ratio) && ratio > 0) 1 / ratio else NA_real_,
          new_over_old = ratio,
          slowdown_pct = if (is.finite(ratio)) 100 * (ratio - 1) else NA_real_,
          median_pair_slowdown_pct = if (length(old)) {
            stats::median(z$slowdown_pct[ok], na.rm = TRUE)
          } else NA_real_,
          paired_new_faster_fraction = if (length(old)) mean(new < old) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

summary <- summarise_timings(timings)

empirical_p <- function(data) {
  rows <- list()
  row_id <- 0L
  for (shape in shape_names) {
    for (workload in workloads) {
      z <- data[data$shape == shape & data$workload == workload &
                  data$successful_pairs > 0, , drop = FALSE]
      z <- z[order(z$n), , drop = FALSE]
      if (nrow(z) < 2L) next
      for (i in seq_len(nrow(z) - 1L)) {
        if (!is.finite(z$before_median_ms[[i]]) ||
            !is.finite(z$before_median_ms[[i + 1L]]) ||
            !is.finite(z$after_median_ms[[i]]) ||
            !is.finite(z$after_median_ms[[i + 1L]])) next
        n_from <- z$n[[i]]
        n_to <- z$n[[i + 1L]]
        p_before <- log(max(z$before_median_ms[[i + 1L]], .Machine$double.xmin) /
                        max(z$before_median_ms[[i]], .Machine$double.xmin)) /
          log(n_to / n_from)
        p_after <- log(max(z$after_median_ms[[i + 1L]], .Machine$double.xmin) /
                       max(z$after_median_ms[[i]], .Machine$double.xmin)) /
          log(n_to / n_from)
        row_id <- row_id + 1L
        rows[[row_id]] <- data.frame(
          shape = shape, workload = workload,
          n_from = as.integer(n_from), n_to = as.integer(n_to),
          before_median_ms = z$before_median_ms[[i]],
          before_median_ms_to = z$before_median_ms[[i + 1L]],
          after_median_ms = z$after_median_ms[[i]],
          after_median_ms_to = z$after_median_ms[[i + 1L]],
          empirical_p_before = p_before, empirical_p_after = p_after,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows)) return(do.call(rbind, rows))
  data.frame(
    shape = character(), workload = character(), n_from = integer(),
    n_to = integer(), before_median_ms = numeric(),
    before_median_ms_to = numeric(), after_median_ms = numeric(),
    after_median_ms_to = numeric(), empirical_p_before = numeric(),
    empirical_p_after = numeric(), stringsAsFactors = FALSE
  )
}

empirical <- empirical_p(summary)

representative_n <- if (5000L %in% n_grid) 5000L else {
  lower <- n_grid[n_grid <= 5000L]
  if (length(lower)) max(lower) else min(n_grid)
}
gate_rows <- list()
gate_id <- 0L
for (workload in workloads) {
  z <- summary[summary$workload == workload & summary$n == representative_n,
               , drop = FALSE]
  for (shape in shape_names) {
    cell <- z[z$shape == shape, , drop = FALSE]
    gate_id <- gate_id + 1L
    slowdown <- if (nrow(cell)) cell$slowdown_pct[[1L]] else NA_real_
    status <- if (!is.finite(slowdown)) "INSUFFICIENT" else
      if (slowdown > 5) "FAIL" else "PASS"
    gate_rows[[gate_id]] <- data.frame(
      scope = "cell", shape = shape, workload = workload,
      representative_n = representative_n,
      before_median_ms = if (nrow(cell)) cell$before_median_ms[[1L]] else NA_real_,
      after_median_ms = if (nrow(cell)) cell$after_median_ms[[1L]] else NA_real_,
      slowdown_pct = slowdown, threshold_pct = 5,
      status = status, stringsAsFactors = FALSE
    )
  }
  finite <- z$slowdown_pct[is.finite(z$slowdown_pct)]
  gate_id <- gate_id + 1L
  aggregate_slowdown <- if (length(finite)) stats::median(finite) else NA_real_
  aggregate_status <- if (!is.finite(aggregate_slowdown)) "INSUFFICIENT" else
    if (aggregate_slowdown > 5) "FAIL" else "PASS"
  gate_rows[[gate_id]] <- data.frame(
    scope = "all_shapes_median", shape = "all_shapes", workload = workload,
    representative_n = representative_n,
    before_median_ms = if (nrow(z)) stats::median(z$before_median_ms, na.rm = TRUE) else NA_real_,
    after_median_ms = if (nrow(z)) stats::median(z$after_median_ms, na.rm = TRUE) else NA_real_,
    slowdown_pct = aggregate_slowdown, threshold_pct = 5,
    status = aggregate_status, stringsAsFactors = FALSE
  )
}
gate <- do.call(rbind, gate_rows)

# -------------------------------------------------------------------------
# Exact provenance and output files.

git_capture <- function(git_args) {
  git_marker <- file.path(repo, ".git")
  if (!file.exists(git_marker) && !dir.exists(git_marker)) return(NA_character_)
  value <- tryCatch(
    system2("git", c("-C", repo, git_args), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  value <- trimws(paste(value, collapse = "\n"))
  if (nzchar(value)) value else NA_character_
}

source_commit <- Sys.getenv(paste0(prefix, "SOURCE_COMMIT"), "")
if (!nzchar(source_commit)) source_commit <- git_capture(c("rev-parse", "HEAD"))
baseline_commit <- Sys.getenv(paste0(prefix, "BASELINE_COMMIT"), "")
if (!nzchar(baseline_commit)) {
  baseline_commit <- git_capture(c("rev-parse", "4dc4fa3^{commit}"))
}
if (is.na(baseline_commit) || !nzchar(baseline_commit)) {
  baseline_commit <- "4dc4fa3b4eb6a0815789bdbd5276ec3063dcbda4"
}
git_status <- git_capture(c("status", "--short", "--untracked-files=no"))
if (is.na(git_status)) {
  git_status <- "not available in ASCII staging; source commit supplied"
}
md5_or_na <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(as.character(tools::md5sum(path)))
}
ascii_only <- function(value) !grepl("[^ -~]", value)
scalar_or_na <- function(value) {
  if (!length(value) || is.na(value[[1L]]) || !nzchar(as.character(value[[1L]]))) {
    return(NA_character_)
  }
  as.character(value[[1L]])
}

provenance <- data.frame(
  key = c(
    "source_commit", "baseline_oracle_commit", "git_status_tracked",
    "check_tree_md5", "oracle_fixture_md5", "R_version", "R_platform",
    "R_arch", "R_home", "R_compiler", "CC", "CXX", "CXX14", "CXX17",
    "package_version", "ape_version", "package_load", "repo_path_ascii",
    "output_path_ascii", "protocol", "operation_count_method",
    "partition", "quick_mode", "repeats", "repeat_schedule", "n_grid",
    "shapes", "timing_timeout_seconds", "representative_n",
    "slowdown_gate_threshold_pct"
  ),
  value = c(
    source_commit, baseline_commit, git_status,
    md5_or_na(file.path(repo, "R", "check_tree.R")), md5_or_na(oracle_file),
    R.version$version.string, R.version$platform, R.version$arch, R.home(),
    scalar_or_na(R.version$compiler), Sys.getenv("CC", "<unset>"),
    Sys.getenv("CXX", "<unset>"), Sys.getenv("CXX14", "<unset>"),
    Sys.getenv("CXX17", "<unset>"),
    as.character(utils::packageVersion("fastphylosig")),
    as.character(utils::packageVersion("ape")), package_load,
    as.character(ascii_only(repo)), as.character(ascii_only(out_dir)),
    "fixed balanced/random/pectinate; positive unit branches; warmup; parity check; fresh serialized clones; alternating old/new order; serial timing; median/IQR",
    "source-aligned deterministic replay; old character lookups and R stack copied slots versus integer adjacency/cursor writes",
    partition_label, as.character(quick_mode), as.character(repeats),
    if (nzchar(repeat_schedule_text)) repeat_schedule_text else "<fixed>",
    paste(n_grid, collapse = ","), paste(shape_names, collapse = ","),
    if (is.finite(timing_timeout)) as.character(timing_timeout) else "Inf",
    as.character(representative_n), "5"
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(fixtures, file.path(out_dir, "fixture_manifest.csv"), row.names = FALSE)
utils::write.csv(operations, file.path(out_dir, "operation_counts.csv"), row.names = FALSE)
utils::write.csv(timings, file.path(out_dir, "paired_timings.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "paired_summary.csv"), row.names = FALSE)
utils::write.csv(empirical, file.path(out_dir, "empirical_p_before_after.csv"), row.names = FALSE)
utils::write.csv(gate, file.path(out_dir, "representative_slowdown_gate.csv"), row.names = FALSE)
utils::write.csv(provenance, file.path(out_dir, "provenance.csv"), row.names = FALSE)

manifest <- data.frame(
  file = sort(list.files(out_dir, full.names = FALSE)),
  stringsAsFactors = FALSE
)
utils::write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)

stage_message(sprintf("completed; output=%s", out_dir))
print(summary)
print(gate)
