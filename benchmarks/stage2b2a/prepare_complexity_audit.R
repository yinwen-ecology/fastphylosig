#!/usr/bin/env Rscript

# fastphylosig 0.2.0 Stage 2B2A: preparation complexity root-cause audit.
#
# This is a read-only benchmark. It loads the current package source, builds
# deterministic tree-shape fixtures, and writes evidence below the requested
# output directory. It does not edit package source or test files.

args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(
  if (length(args) >= 1L) args[[1L]] else getwd(),
  winslash = "/", mustWork = TRUE
)
out_dir <- if (length(args) >= 2L && nzchar(args[[2L]])) {
  normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
} else {
  file.path(tempdir(), "fastphylosig-stage2b2a")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")
options(stringsAsFactors = FALSE)

quick_mode <- identical(Sys.getenv("FASTPHYLOSIG_STAGE2B2A_QUICK", "0"), "1")
if (length(args) >= 3L && identical(args[[3L]], "quick")) quick_mode <- TRUE
shape_names <- c("balanced", "random", "pectinate")
n_grid <- if (quick_mode) c(50L, 100L, 500L, 1000L, 2000L) else
  c(500L, 1000L, 2000L, 5000L, 10000L, 20000L)
shape_override <- Sys.getenv("FASTPHYLOSIG_STAGE2B2A_SHAPES", "")
if (nzchar(shape_override)) {
  shape_names <- unique(trimws(strsplit(shape_override, ",", fixed = TRUE)[[1L]]))
  shape_names <- shape_names[nzchar(shape_names)]
  if (!length(shape_names) || any(!shape_names %in% c("balanced", "random", "pectinate"))) {
    stop("FASTPHYLOSIG_STAGE2B2A_SHAPES must contain balanced, random, and/or pectinate.",
         call. = FALSE)
  }
}
n_override <- Sys.getenv("FASTPHYLOSIG_STAGE2B2A_N_GRID", "")
if (nzchar(n_override)) {
  n_grid <- suppressWarnings(as.integer(trimws(strsplit(n_override, ",", fixed = TRUE)[[1L]])))
  if (!length(n_grid) || anyNA(n_grid) || any(n_grid < 2L)) {
    stop("FASTPHYLOSIG_STAGE2B2A_N_GRID must contain integers >= 2.", call. = FALSE)
  }
  n_grid <- unique(n_grid)
}
repeats <- as.integer(Sys.getenv(
  "FASTPHYLOSIG_STAGE2B2A_REPEATS", if (quick_mode) "3" else "10"
))
if (!is.finite(repeats) || repeats < 1L) repeats <- if (quick_mode) 3L else 10L
repeats <- max(1L, repeats)
repeat_schedule_text <- Sys.getenv("FASTPHYLOSIG_STAGE2B2A_REPEATS_BY_N", "")
repeat_schedule <- numeric()
if (nzchar(repeat_schedule_text)) {
  schedule_parts <- trimws(strsplit(repeat_schedule_text, ",", fixed = TRUE)[[1L]])
  schedule_pairs <- strsplit(schedule_parts, ":", fixed = TRUE)
  if (any(lengths(schedule_pairs) != 2L)) {
    stop("FASTPHYLOSIG_STAGE2B2A_REPEATS_BY_N must use n:repeats pairs.",
         call. = FALSE)
  }
  repeat_schedule <- suppressWarnings(vapply(schedule_pairs, function(x) {
    as.numeric(x[[2L]])
  }, numeric(1)))
  schedule_n <- suppressWarnings(vapply(schedule_pairs, function(x) {
    as.numeric(x[[1L]])
  }, numeric(1)))
  if (any(!is.finite(schedule_n)) || any(schedule_n < 2) ||
      any(!is.finite(repeat_schedule)) || any(repeat_schedule < 1) ||
      any(repeat_schedule != floor(repeat_schedule))) {
    stop("FASTPHYLOSIG_STAGE2B2A_REPEATS_BY_N contains invalid n:repeats values.",
         call. = FALSE)
  }
  names(repeat_schedule) <- as.character(as.integer(schedule_n))
  repeats <- max(repeats, max(repeat_schedule))
}
repeats_for_n <- function(n) {
  value <- if (length(repeat_schedule) && as.character(n) %in% names(repeat_schedule)) {
    repeat_schedule[[as.character(n)]]
  } else repeats
  as.integer(value)
}
timing_timeout <- as.numeric(Sys.getenv("FASTPHYLOSIG_STAGE2B2A_TIMING_TIMEOUT", "0"))
if (!is.finite(timing_timeout) || timing_timeout <= 0) timing_timeout <- Inf
probe_helper <- Sys.getenv("FASTPHYLOSIG_STAGE2B2A_PROBE_HELPER", "")
probe_shape <- Sys.getenv("FASTPHYLOSIG_STAGE2B2A_PROBE_SHAPE", "pectinate")
probe_n <- as.integer(Sys.getenv("FASTPHYLOSIG_STAGE2B2A_PROBE_N", "20000"))
probe_timeout <- as.numeric(Sys.getenv("FASTPHYLOSIG_STAGE2B2A_PROBE_TIMEOUT", "60"))
probe_mode <- nzchar(probe_helper)
if (probe_mode) {
  if (!probe_shape %in% c("balanced", "random", "pectinate")) {
    stop("FASTPHYLOSIG_STAGE2B2A_PROBE_SHAPE is not supported.", call. = FALSE)
  }
  if (!is.finite(probe_n) || probe_n < 2L) {
    stop("FASTPHYLOSIG_STAGE2B2A_PROBE_N must be at least two.", call. = FALSE)
  }
  if (!is.finite(probe_timeout) || probe_timeout <= 0) probe_timeout <- 60
  shape_names <- probe_shape
  n_grid <- as.integer(probe_n)
  repeats <- 1L
}

elapsed_seconds <- function() unname(proc.time()[["elapsed"]])
time_call <- function(fun, timeout_seconds = Inf) {
  limited <- is.finite(timeout_seconds) && timeout_seconds > 0
  if (limited) setTimeLimit(cpu = Inf, elapsed = timeout_seconds, transient = TRUE)
  start <- elapsed_seconds()
  status <- "ok"
  error_message <- ""
  error_class <- ""
  value <- tryCatch(
    fun(),
    error = function(e) {
      status <<- if (inherits(e, "elapsedTimeLimit")) "censored" else "error"
      error_message <<- conditionMessage(e)
      error_class <<- paste(class(e), collapse = "/")
      NULL
    }
  )
  elapsed <- max(0, elapsed_seconds() - start)
  if (limited && identical(status, "ok") && elapsed > timeout_seconds) {
    status <- "censored_overrun"
    error_message <- sprintf("elapsed time %.3fs exceeded %.3fs timing limit.",
                             elapsed, timeout_seconds)
    error_class <- "timing_overrun"
  }
  if (limited) setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  list(value = value, elapsed = elapsed, status = status,
       error_message = error_message, error_class = error_class)
}
time_only <- function(fun) {
  start <- elapsed_seconds()
  invisible(fun())
  max(0, elapsed_seconds() - start)
}
clone_fixture <- function(tree) unserialize(serialize(tree, NULL))

# Load the current source commit without writing into the repository. The
# installed-library route is useful for an ASCII staging copy where pkgload is
# unavailable; source loading remains the default for an audit checkout.
installed_lib <- Sys.getenv("FASTPHYLOSIG_STAGE2B2A_LIB", "")
if (nzchar(installed_lib) && dir.exists(file.path(installed_lib, "fastphylosig"))) {
  suppressPackageStartupMessages(
    library("fastphylosig", lib.loc = installed_lib, character.only = TRUE)
  )
  package_load <- paste0("installed library: ", normalizePath(installed_lib,
                                                               winslash = "/"))
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  suppressPackageStartupMessages(pkgload::load_all(repo, quiet = TRUE))
  package_load <- "pkgload::load_all(source commit)"
} else {
  stop("pkgload is required unless FASTPHYLOSIG_STAGE2B2A_LIB is set.",
       call. = FALSE)
}
suppressPackageStartupMessages(library(ape))
ns <- asNamespace("fastphylosig")
get_internal <- function(name) get(name, envir = ns, inherits = FALSE)

git_commit <- Sys.getenv("FASTPHYLOSIG_SOURCE_COMMIT", "")
if (!nzchar(git_commit)) {
  git_commit <- tryCatch(
    system2("git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE,
            stderr = FALSE)[[1L]],
    error = function(e) NA_character_
  )
}
git_commit <- if (length(git_commit) && nzchar(git_commit)) git_commit else NA_character_

compiler_text <- paste(
  paste0("R.version.compiler=", R.version$compiler),
  paste0("CXX=", Sys.getenv("CXX", unset = "<unset>")),
  paste0("CXX14=", Sys.getenv("CXX14", unset = "<unset>")),
  sep = "; "
)
provenance <- data.frame(
  key = c("source_commit", "source_version", "R_version", "R_platform",
          "R_arch", "R_home", "compiler", "OpenMP", "package_load",
          "protocol", "quick_mode", "repeats", "repeat_schedule",
          "timing_timeout_seconds", "n_grid", "shapes"),
  value = c(
    git_commit,
    as.character(utils::packageVersion("fastphylosig")),
    R.version.string, R.version$platform, R.version$arch, R.home(),
    compiler_text,
    if (isTRUE(capabilities("long.double"))) "R capabilities recorded; build OpenMP not inferred" else "R capabilities recorded",
    package_load,
    "fixed balanced/random/pectinate topology; fixed positive branch lengths; warmup; alternating copies; serial timing; median/IQR",
    as.character(quick_mode), as.character(repeats),
    if (nzchar(repeat_schedule_text)) repeat_schedule_text else "<fixed>",
    if (is.finite(timing_timeout)) as.character(timing_timeout) else "Inf",
    paste(n_grid, collapse = ","), paste(shape_names, collapse = ",")),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(out_dir, "provenance.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# Deterministic tree fixtures.

tip_labels <- function(n) paste0("sp", seq_len(n))

make_balanced <- function(n) {
  n <- as.integer(n)
  if (n < 2L) stop("n must be at least two.", call. = FALSE)
  edges <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  edge_i <- 0L
  n_internal <- n - 1L
  next_internal <- n + 1L
  root <- next_internal
  next_internal <- next_internal + 1L
  # Use an explicit interval queue so even the largest fixture has no
  # recursion-depth dependency.
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
    mid <- floor((lo + hi) / 2L)
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
      edges[edge_i, ] <- c(node, child)
    }
  }
  edges <- edges[seq_len(edge_i), , drop = FALSE]
  tree <- list(edge = matrix(as.integer(edges), ncol = 2L),
               tip.label = tip_labels(n),
               edge.length = rep(1, nrow(edges)), Nnode = n - 1L)
  class(tree) <- "phylo"
  tree$root <- root
  ape::reorder.phylo(tree, order = "postorder")
}

make_pectinate <- function(n) {
  n <- as.integer(n)
  if (n < 2L) stop("n must be at least two.", call. = FALSE)
  root <- n + 1L
  next_internal <- root + 1L
  if (n == 2L) {
    edges <- matrix(c(root, 1L, root, 2L), ncol = 2L, byrow = TRUE)
  } else {
    edges <- matrix(c(root, 1L, root, next_internal), ncol = 2L,
                    byrow = TRUE)
    current <- next_internal
    next_internal <- next_internal + 1L
    if (n > 3L) {
      for (tip in 2:(n - 2L)) {
        edges <- rbind(edges, c(current, tip),
                       c(current, next_internal))
        current <- next_internal
        next_internal <- next_internal + 1L
      }
    }
    edges <- rbind(edges, c(current, n - 1L), c(current, n))
  }
  tree <- list(edge = matrix(as.integer(edges), ncol = 2L),
               tip.label = tip_labels(n),
               edge.length = rep(1, nrow(edges)), Nnode = n - 1L)
  class(tree) <- "phylo"
  tree$root <- root
  ape::reorder.phylo(tree, order = "postorder")
}

make_random <- function(n) {
  set.seed(2072000L + as.integer(n))
  tree <- ape::rtree(n, tip.label = tip_labels(n))
  tree <- ape::reorder.phylo(tree, order = "postorder")
  # Identical positive edge lengths isolate shape effects from branch-length
  # generation and avoid introducing a numerical workload into this audit.
  tree$edge.length <- rep(1, nrow(tree$edge))
  tree
}

make_fixture <- function(shape, n) {
  switch(shape,
         balanced = make_balanced(n),
         random = make_random(n),
         pectinate = make_pectinate(n),
         stop("unknown shape", call. = FALSE))
}

validate_fixture <- function(tree, n, shape) {
  stopifnot(inherits(tree, "phylo"), ape::Ntip(tree) == n,
            ape::Nnode(tree) == n - 1L, nrow(tree$edge) == 2L * n - 2L,
            length(tree$edge.length) == nrow(tree$edge),
            all(is.finite(tree$edge.length)), all(tree$edge.length > 0),
            !anyDuplicated(tree$tip.label))
  check <- get_internal(".inspect_tree_core")(tree, signal = c("K", "lambda", "D", "Delta"))
  if (!isTRUE(check$tree_summary$valid)) {
    stop(sprintf("%s n=%d fixture did not pass generic inspection.", shape, n),
         call. = FALSE)
  }
  invisible(tree)
}

# -------------------------------------------------------------------------
# Shape-independent primitive phases. These call the actual private helpers
# without tracing, so their wall times are direct helper timings.

root_of <- function(tree) {
  edge <- tree$edge
  roots <- setdiff(unique(edge[, 1L]), unique(edge[, 2L]))
  if (length(roots) == 1L) as.integer(roots[[1L]]) else NA_integer_
}

primitive_calls <- function(tree) {
  edge <- tree$edge
  tip <- as.character(tree$tip.label)
  n_tip <- length(tip)
  root <- root_of(tree)
  children <- get_internal(".edge_children_index")(edge)
  internals <- sort(unique(c(edge[, 1L], edge[, 2L])))
  internals <- internals[internals > n_tip]
  list(
    inspect_tree_core = function() get_internal(".inspect_tree_core")(
      tree, signal = c("K", "lambda", "D", "Delta")
    ),
    safe_canonicalize_core = function() {
      get_internal(".safe_canonicalize_core")(tree)
    },
    canonical_tree_signature = function() {
      get_internal(".canonical_tree_signature")(tree)
    },
    edge_normalization = function() {
      values <- suppressWarnings(as.numeric(edge))
      matrix(as.integer(values), ncol = 2L)
    },
    root_identification = function() root_of(tree),
    parent_child_construction = function() get_internal(".edge_children_index")(edge),
    preorder_connectivity = function() {
      if (!is.finite(root)) return(FALSE)
      seen <- rep(FALSE, n_tip + as.integer(tree$Nnode))
      stack <- as.integer(root)
      seen[[root]] <- TRUE
      while (length(stack)) {
        node <- stack[[length(stack)]]
        stack <- stack[-length(stack)]
        kids <- children[[as.character(node)]]
        if (!is.null(kids)) for (kid in kids) {
          kid <- as.integer(kid)
          if (seen[[kid]]) next
          seen[[kid]] <- TRUE
          stack <- c(stack, kid)
        }
      }
      all(seen)
    },
    postorder_reordering = function() ape::reorder.phylo(tree, order = "postorder"),
    root_distance = function() {
      n_total <- n_tip + as.integer(tree$Nnode)
      distance <- rep(NA_real_, n_total)
      if (!is.finite(root)) return(distance)
      distance[[root]] <- 0
      rows_by_parent <- split(seq_len(nrow(edge)), edge[, 1L])
      stack <- as.integer(root)
      while (length(stack)) {
        node <- stack[[length(stack)]]
        stack <- stack[-length(stack)]
        rows <- rows_by_parent[[as.character(node)]]
        if (is.null(rows)) next
        for (row in rows) {
          kid <- as.integer(edge[row, 2L])
          distance[[kid]] <- distance[[node]] + tree$edge.length[[row]]
          stack <- c(stack, kid)
        }
      }
      distance
    },
    tip_mapping = function() match(tip, tip),
    signature_descendant_traversal = function() {
      if (!is.finite(root)) return(character())
      get_internal(".descendant_labels_iterative")(root, children, tip, n_tip)
    },
    signature_descendant_keys = function() {
      if (!length(internals) || !is.finite(root)) return(character())
      get_internal(".descendant_keys_iterative")(
        setdiff(internals, root), children, tip, n_tip
      )
    }
  )
}

time_primitives <- function(tree, timeout_seconds = Inf) {
  calls <- primitive_calls(tree)
  out <- lapply(names(calls), function(name) {
    timed <- time_call(calls[[name]], timeout_seconds = timeout_seconds)
    data.frame(phase = name, elapsed_ms = 1000 * timed$elapsed,
               status = timed$status, error_class = timed$error_class,
               error_message = timed$error_message,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# -------------------------------------------------------------------------
# Source-derived operation counts. These counts are intentionally separate
# from timing; they expose repeated descendant aggregation without claiming a
# mathematical complexity proof.

descendant_tip_counts <- function(tree) {
  edge <- tree$edge
  n_tip <- length(tree$tip.label)
  nodes <- sort(unique(c(edge[, 1L], edge[, 2L])))
  children_named <- get_internal(".edge_children_index")(edge)
  # Convert names to integer slots once. The counting pass itself must not
  # repeat children[[as.character(node)]] or counts[as.character(kids)],
  # otherwise the audit's own instrumentation would create a second
  # name-lookup scaling problem on pectinate fixtures.
  children <- vector("list", max(nodes))
  child_keys <- suppressWarnings(as.integer(names(children_named)))
  children[child_keys] <- unname(children_named)
  counts <- numeric(max(nodes))
  counts[seq_len(n_tip)] <- 1
  order_nodes <- rev(nodes[nodes > n_tip])
  for (node in order_nodes) {
    kids <- children[[node]]
    if (length(kids)) counts[[node]] <- sum(counts[kids])
  }
  counts
}

source_operation_counts <- function(tree) {
  edge <- tree$edge
  n_tip <- length(tree$tip.label)
  n_edge <- nrow(edge)
  n_total <- n_tip + as.integer(tree$Nnode)
  internals <- sort(unique(c(edge[, 1L], edge[, 2L])))
  internals <- internals[internals > n_tip]
  root <- root_of(tree)
  tips_by_node <- descendant_tip_counts(tree)
  descendant_key_tip_visits <- if (length(internals) && is.finite(root)) {
    sum(tips_by_node[setdiff(internals, root)])
  } else 0
  n_internal <- length(internals)
  n_parent_groups <- length(unique(edge[, 1L]))
  n_total_nodes <- n_total
  # For valid binary fixtures, inspect_tree_core indexes its named adjacency
  # list once for every visited node in each of its two traversals. The
  # canonicalization path indexes the list once per internal node in each
  # root signature and twice for each non-root internal in the memoized key
  # pass. These are dynamic lookup counts, not a claim that R implements
  # named lookup as a linear scan.
  adjacency_lookup_inspection <- 2L * n_total_nodes
  adjacency_lookup_signature <- 2L * n_internal
  adjacency_lookup_canonicalization <- 2L * n_internal +
    2L * max(0L, n_internal - 1L)
  adjacency_lookup_total <- adjacency_lookup_inspection +
    adjacency_lookup_canonicalization
  linear_name_comparison_upper <- as.double(adjacency_lookup_total) *
    as.double(n_parent_groups)
  # The production canonicalizer calls the iterative signature before and
  # after remapping, and its all-internal key pass once in between. Each
  # signature root traversal emits one label per tip.
  data.frame(
    full_edge_rows = n_edge,
    edge_scan_passes_inspection = 8L,
    edge_scan_passes_signature = 5L,
    edge_scan_passes_canonicalization = 7L,
    per_node_full_edge_searches_inspection = 0L,
    per_node_full_edge_searches_canonicalization = 0L,
    named_adjacency_lookup_inspection = adjacency_lookup_inspection,
    named_adjacency_lookup_signature = adjacency_lookup_signature,
    named_adjacency_lookup_canonicalization = adjacency_lookup_canonicalization,
    named_adjacency_lookup_total = adjacency_lookup_total,
    named_adjacency_parent_groups = n_parent_groups,
    estimated_name_comparisons_linear_upper = linear_name_comparison_upper,
    descendant_traversal_calls_inspection = 0L,
    descendant_traversal_calls_signature = if (is.finite(root)) 1L else 0L,
    descendant_traversal_calls_canonicalization = if (is.finite(root)) 2L else 0L,
    descendant_key_passes_canonicalization = if (length(setdiff(internals, root))) 1L else 0L,
    descendant_tip_visits_signature = if (is.finite(root)) n_tip else 0L,
    descendant_tip_visits_canonicalization = if (is.finite(root)) {
      2 * n_tip + descendant_key_tip_visits
    } else 0,
    descendant_key_tip_visits = descendant_key_tip_visits,
    descendant_key_sort_calls = max(0L, n_internal - 1L),
    descendant_key_collapse_calls = max(0L, n_internal - 1L),
    complete_parent_child_reconstructions = 5L,
    order_calls_canonicalization = if (length(setdiff(internals, root))) 1L else 0L,
    sort_calls_canonicalization = if (is.finite(root)) n_internal + 6L else 0L,
    n_tip = n_tip, n_internal = n_internal, n_total = n_total,
    stringsAsFactors = FALSE
  )
}

source_pattern_audit <- function() {
  bodies <- list(
    inspect_tree_core = get_internal(".inspect_tree_core"),
    safe_canonicalize_core = get_internal(".safe_canonicalize_core"),
    canonical_tree_signature = get_internal(".canonical_tree_signature")
  )
  patterns <- c(
    "which(", "match(", "%in%", "edge[edge", "edge[, 1L] ==",
    "edge[, 2L] ==", "descendant", "order(", "sort(", "duplicated(",
    "anyDuplicated(", "split(", "tabulate("
  )
  rows <- list()
  id <- 0L
  for (fun_name in names(bodies)) {
    text <- paste(deparse(body(bodies[[fun_name]])), collapse = "\n")
    for (pattern in patterns) {
      id <- id + 1L
      rows[[id]] <- data.frame(
        function_name = fun_name, pattern = pattern,
        occurrences = lengths(regmatches(text, gregexpr(pattern, text, fixed = TRUE))),
        per_node_loop = grepl("for \\(.*\\).*", text) && grepl(pattern, text, fixed = TRUE),
        interpretation = if (grepl(pattern, text, fixed = TRUE)) "source occurrence; inspect loop context" else "not present",
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# -------------------------------------------------------------------------
# One-shot helper call-count audit. Bindings are restored immediately and no
# count run contributes to the timing tables.

with_helper_counts <- function(fun) {
  names_to_count <- c(
    ".inspect_tree_core", ".safe_canonicalize_core",
    ".canonical_tree_signature", ".edge_children_index",
    ".descendant_labels_iterative", ".descendant_keys_iterative",
    ".prepare_tree_core", ".prepare_tree_subset", ".attach_structural_evidence"
  )
  counts <- new.env(parent = emptyenv())
  for (name in names_to_count) counts[[name]] <- 0L
  originals <- lapply(names_to_count, get_internal)
  names(originals) <- names_to_count
  locked <- vapply(names_to_count, bindingIsLocked, logical(1), env = ns)
  for (name in names_to_count) {
    if (locked[[name]]) unlockBinding(name, ns)
    original <- originals[[name]]
    assign(name, local({
      nm <- name
      fn <- original
      function(...) {
        counts[[nm]] <- counts[[nm]] + 1L
        fn(...)
      }
    }), ns)
  }
  on.exit({
    for (name in names_to_count) {
      if (bindingIsLocked(name, ns)) unlockBinding(name, ns)
      assign(name, originals[[name]], ns)
      if (locked[[name]]) lockBinding(name, ns)
    }
  }, add = TRUE)
  value <- fun()
  list(value = value, counts = setNames(
    vapply(names_to_count, function(name) as.integer(counts[[name]]), integer(1)),
    names_to_count
  ))
}

# -------------------------------------------------------------------------
# Build fixtures and run the serial timing protocol.

stage_message <- function(text) message("[stage2b2a] ", text)
stage_message(sprintf("starting; commit=%s; n=%s; repeats=%d; quick=%s",
                      git_commit, paste(n_grid, collapse = ","), repeats,
                      quick_mode))

fixture_cache <- new.env(parent = emptyenv())
for (shape in shape_names) {
  for (n in n_grid) {
    key <- paste(shape, n, sep = "::")
    stage_message(sprintf("constructing %s n=%d", shape, n))
    tree <- make_fixture(shape, n)
    validate_fixture(tree, n, shape)
    fixture_cache[[key]] <- tree
  }
}

fixture_only <- identical(Sys.getenv("FASTPHYLOSIG_STAGE2B2A_FIXTURE_ONLY", "0"), "1")
if (fixture_only) {
  fixture_rows <- list()
  fixture_id <- 0L
  for (shape in shape_names) {
    for (n in n_grid) {
      tree <- fixture_cache[[paste(shape, n, sep = "::")]]
      fixture_id <- fixture_id + 1L
      fixture_rows[[fixture_id]] <- data.frame(
        shape = shape, n = n, root = root_of(tree),
        n_tip = ape::Ntip(tree), n_internal = ape::Nnode(tree),
        n_edge = nrow(tree$edge), valid = TRUE, stringsAsFactors = FALSE
      )
    }
  }
  write.csv(do.call(rbind, fixture_rows),
            file.path(out_dir, "fixture_validation.csv"), row.names = FALSE)
  stage_message("fixture-only validation completed")
  quit(save = "no", status = 0L)
}

if (probe_mode) {
  tree <- fixture_cache[[paste(probe_shape, probe_n, sep = "::")]]
  edge <- tree$edge
  tip <- as.character(tree$tip.label)
  n_tip <- length(tip)
  root <- root_of(tree)
  children <- get_internal(".edge_children_index")(edge)
  internal_nodes <- sort(unique(c(edge[, 1L], edge[, 2L])))
  internal_nodes <- internal_nodes[internal_nodes > n_tip]
  probe_fun <- switch(
    probe_helper,
    prepare_tree = function() prepare_tree(tree),
    inspect_tree_core = function() get_internal(".inspect_tree_core")(
      tree, signal = c("K", "lambda", "D", "Delta")
    ),
    safe_canonicalize_core = function() get_internal(".safe_canonicalize_core")(tree),
    canonical_tree_signature = function() get_internal(".canonical_tree_signature")(tree),
    edge_children_index = function() get_internal(".edge_children_index")(edge),
    descendant_labels_iterative = function() get_internal(
      ".descendant_labels_iterative")(root, children, tip, n_tip),
    descendant_keys_iterative = function() get_internal(
      ".descendant_keys_iterative")(setdiff(internal_nodes, root), children, tip, n_tip),
    stop("unknown FASTPHYLOSIG_STAGE2B2A_PROBE_HELPER.", call. = FALSE)
  )
  setTimeLimit(cpu = Inf, elapsed = probe_timeout, transient = TRUE)
  started <- elapsed_seconds()
  status <- "ok"
  error_message <- ""
  error_class <- ""
  result <- tryCatch(
    {
      value <- probe_fun()
      invisible(value)
      TRUE
    },
    error = function(e) {
      status <<- if (inherits(e, "elapsedTimeLimit")) "timeout" else "error"
      error_message <<- conditionMessage(e)
      error_class <<- paste(class(e), collapse = "/")
      FALSE
    }
  )
  elapsed <- elapsed_seconds() - started
  if (identical(status, "ok") && elapsed > probe_timeout) {
    status <- "censored_overrun"
    error_message <- sprintf("elapsed time %.3fs exceeded %.3fs probe limit.",
                             elapsed, probe_timeout)
    error_class <- "timing_overrun"
  }
  setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  write.csv(data.frame(
    shape = probe_shape, n = probe_n, helper = probe_helper,
    timeout_seconds = probe_timeout, elapsed_seconds = elapsed,
    status = status, error_class = error_class,
    error_message = error_message, stringsAsFactors = FALSE
  ), file.path(out_dir, "single_helper_probe.csv"), row.names = FALSE)
  stage_message(sprintf("probe %s %s n=%d: %s in %.3fs",
                        probe_shape, probe_helper, probe_n, status, elapsed))
  quit(save = "no", status = if (isTRUE(result)) 0L else 2L)
}

timing_rows <- list()
primitive_rows <- list()
count_rows <- list()
operation_rows <- list()
calibration_rows <- list()
id_timing <- id_primitive <- id_count <- id_operation <- id_cal <- 0L

for (shape in shape_names) {
  for (n in n_grid) {
    key <- paste(shape, n, sep = "::")
    tree_a <- fixture_cache[[key]]
    tree_b <- clone_fixture(tree_a)
    n_repeats <- repeats_for_n(n)
    stage_message(sprintf("timing %s n=%d (%d repeats, serial%s)", shape, n,
                          n_repeats,
                          if (is.finite(timing_timeout)) {
                            paste0(", timeout=", timing_timeout, "s")
                          } else ""))

    # Warmup each direct helper and public preparation path.
    invisible(time_call(function() prepare_tree(tree_a), timeout_seconds = timing_timeout))
    invisible(time_primitives(tree_a, timeout_seconds = timing_timeout))

    for (rep in seq_len(n_repeats)) {
      fixture_key <- if (rep %% 2L) "a" else "b"
      tree <- if (identical(fixture_key, "a")) tree_a else tree_b
      timed <- time_call(function() prepare_tree(tree),
                         timeout_seconds = timing_timeout)
      id_timing <- id_timing + 1L
      timing_rows[[id_timing]] <- data.frame(
        shape = shape, n = n, `repeat` = rep, position = 1L,
        order = fixture_key, elapsed_ms = 1000 * timed$elapsed,
        status = timed$status, error_class = timed$error_class,
        error_message = timed$error_message,
        stringsAsFactors = FALSE
      )
      phases <- time_primitives(tree, timeout_seconds = timing_timeout)
      for (j in seq_len(nrow(phases))) {
        id_primitive <- id_primitive + 1L
        primitive_rows[[id_primitive]] <- data.frame(
          shape = shape, n = n, `repeat` = rep, phase = phases$phase[[j]],
          elapsed_ms = phases$elapsed_ms[[j]], status = phases$status[[j]],
          error_class = phases$error_class[[j]],
          error_message = phases$error_message[[j]], stringsAsFactors = FALSE
        )
      }
    }

    # Helper calls and operation counts are deliberately not timed. They are
    # collected once on a clone after the formal timing block.
    counted_timed <- time_call(function() with_helper_counts(function() {
      prepare_tree(clone_fixture(tree_a))
    }), timeout_seconds = timing_timeout)
    counted <- counted_timed$value
    counted_names <- if (is.list(counted) && !is.null(counted$counts)) {
      names(counted$counts)
    } else {
      c(".inspect_tree_core", ".safe_canonicalize_core",
        ".canonical_tree_signature", ".edge_children_index",
        ".descendant_labels_iterative", ".descendant_keys_iterative",
        ".prepare_tree_core", ".prepare_tree_subset", ".attach_structural_evidence")
    }
    for (helper in counted_names) {
      id_count <- id_count + 1L
      count_rows[[id_count]] <- data.frame(
        shape = shape, n = n, route = "prepare_tree", helper = helper,
        calls = if (is.list(counted) && !is.null(counted$counts)) {
          counted$counts[[helper]]
        } else NA_integer_,
        status = counted_timed$status, error_class = counted_timed$error_class,
        error_message = counted_timed$error_message, stringsAsFactors = FALSE
      )
    }
    ops <- source_operation_counts(tree_a)
    ops$shape <- shape
    ops$n <- n
    operation_rows[[length(operation_rows) + 1L]] <- ops
    # Repeat the wrapper calibration in alternating order. The wrapper is
    # installed only for a single non-timed count call and is restored before
    # the next run. Its timing is diagnostic and never replaces direct helper
    # timings above.
    calibration_repeats <- max(3L, min(5L, n_repeats))
    for (cal_rep in seq_len(calibration_repeats)) {
      cal_order <- if (cal_rep %% 2L) c("plain", "counted") else
        c("counted", "plain")
      cal_values <- list()
      for (route in cal_order) {
        if (identical(route, "counted")) {
          counted <- with_helper_counts(function() {
            time_only(function() get_internal(".canonical_tree_signature")(tree_a))
          })
          cal_values$counted <- 1000 * counted$value
        } else {
          cal_values$plain <- 1000 * time_only(
            function() get_internal(".canonical_tree_signature")(tree_a)
          )
        }
      }
      plain_ms <- cal_values$plain
      counted_ms <- cal_values$counted
      id_cal <- id_cal + 1L
      calibration_rows[[id_cal]] <- data.frame(
        shape = shape, n = n, `repeat` = cal_rep,
        order = paste(cal_order, collapse = ">"),
        plain_signature_ms = plain_ms,
        counted_signature_ms = counted_ms,
        overhead_pct = if (plain_ms > 0) {
          100 * (counted_ms - plain_ms) / plain_ms
        } else NA_real_,
        stringsAsFactors = FALSE
      )
    }
    rm(tree_a, tree_b)
    gc(FALSE)
  }
}

timings <- do.call(rbind, timing_rows)
primitive_timings <- do.call(rbind, primitive_rows)
helper_counts <- do.call(rbind, count_rows)
operation_counts <- do.call(rbind, operation_rows)
calibration <- do.call(rbind, calibration_rows)
write.csv(timings, file.path(out_dir, "prepare_tree_timings.csv"), row.names = FALSE)
write.csv(primitive_timings, file.path(out_dir, "helper_phase_timings.csv"), row.names = FALSE)
write.csv(helper_counts, file.path(out_dir, "helper_call_counts.csv"), row.names = FALSE)
write.csv(operation_counts, file.path(out_dir, "operation_counts.csv"), row.names = FALSE)
write.csv(calibration, file.path(out_dir, "instrumentation_calibration.csv"), row.names = FALSE)

calibration_summary <- do.call(rbind, lapply(split(
  calibration, interaction(calibration$shape, calibration$n,
                           drop = TRUE, lex.order = TRUE), drop = TRUE
), function(z) {
  finite_overhead <- z$overhead_pct[is.finite(z$overhead_pct)]
  data.frame(
    shape = as.character(z$shape[[1L]]), n = as.integer(z$n[[1L]]),
    repeats = nrow(z),
    median_plain_signature_ms = stats::median(z$plain_signature_ms),
    median_counted_signature_ms = stats::median(z$counted_signature_ms),
    median_overhead_pct = if (length(finite_overhead)) {
      stats::median(finite_overhead)
    } else NA_real_,
    max_abs_overhead_pct = if (length(finite_overhead)) {
      max(abs(finite_overhead))
    } else NA_real_,
    stringsAsFactors = FALSE
  )
}))
write.csv(calibration_summary,
          file.path(out_dir, "instrumentation_calibration_summary.csv"),
          row.names = FALSE)

summarise <- function(data, groups) {
  split_data <- split(data, interaction(data[groups], drop = TRUE, lex.order = TRUE), drop = TRUE)
  do.call(rbind, lapply(split_data, function(z) {
    out <- as.list(z[1L, groups, drop = FALSE])
    out$repeats <- nrow(z)
    ok <- if ("status" %in% names(z)) z$status == "ok" else rep(TRUE, nrow(z))
    elapsed <- z$elapsed_ms[ok & is.finite(z$elapsed_ms)]
    out$successful_repeats <- length(elapsed)
    out$censored_or_failed <- nrow(z) - length(elapsed)
    out$median_ms <- if (length(elapsed)) stats::median(elapsed) else NA_real_
    out$IQR_ms <- if (length(elapsed)) stats::IQR(elapsed) else NA_real_
    as.data.frame(out, stringsAsFactors = FALSE)
  }))
}
timing_summary <- summarise(timings, c("shape", "n"))
primitive_summary <- summarise(primitive_timings, c("shape", "n", "phase"))
write.csv(timing_summary, file.path(out_dir, "prepare_tree_summary.csv"), row.names = FALSE)
write.csv(primitive_summary, file.path(out_dir, "helper_phase_summary.csv"), row.names = FALSE)

empirical_exponent <- function(summary, time_col = "median_ms") {
  groups <- intersect(c("shape", "phase"), names(summary))
  split_key <- interaction(summary[groups], drop = TRUE, lex.order = TRUE)
  rows <- list()
  id <- 0L
  for (key in unique(split_key)) {
    z <- summary[split_key == key, , drop = FALSE]
    z <- z[order(z$n), , drop = FALSE]
    z <- z[is.finite(z[[time_col]]) & z$n > 0, , drop = FALSE]
    if (nrow(z) < 2L) next
    for (i in seq_len(nrow(z) - 1L)) {
      id <- id + 1L
      t1 <- max(as.numeric(z[[time_col]][[i]]), .Machine$double.xmin)
      t2 <- max(as.numeric(z[[time_col]][[i + 1L]]), .Machine$double.xmin)
      rows[[id]] <- cbind(z[i, groups, drop = FALSE], data.frame(
        n_from = z$n[[i]], n_to = z$n[[i + 1L]],
        elapsed_from_ms = z[[time_col]][[i]], elapsed_to_ms = z[[time_col]][[i + 1L]],
        empirical_p = log(t2 / t1) / log(z$n[[i + 1L]] / z$n[[i]]),
        stringsAsFactors = FALSE
      ))
    }
  }
  do.call(rbind, rows)
}
write.csv(empirical_exponent(timing_summary), file.path(out_dir, "prepare_tree_empirical_exponent.csv"), row.names = FALSE)
write.csv(empirical_exponent(primitive_summary), file.path(out_dir, "helper_phase_empirical_exponent.csv"), row.names = FALSE)
write.csv(source_pattern_audit(), file.path(out_dir, "source_operation_pattern_audit.csv"), row.names = FALSE)

manifest <- data.frame(
  file = sort(list.files(out_dir, full.names = FALSE)),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)

stage_message("completed")
print(timing_summary)
print(primitive_summary)
print(empirical_exponent(timing_summary))
