#!/usr/bin/env Rscript

# fastphylosig 0.2.0 Stage 2B2B Candidate 2.
# Candidate-1 post-freeze direct baseline. This script is read-only with
# respect to package source and writes only to the requested output directory.

args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(
  if (length(args) >= 1L) args[[1L]] else getwd(),
  winslash = "/", mustWork = TRUE
)
out_dir <- if (length(args) >= 2L && nzchar(args[[2L]])) {
  normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
} else {
  file.path(tempdir(), "fastphylosig-stage2b2b2-baseline")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")
options(stringsAsFactors = FALSE)

quick_mode <- identical(Sys.getenv("FASTPHYLOSIG_STAGE2B2B2_QUICK", "0"), "1")
if (length(args) >= 3L && identical(args[[3L]], "quick")) quick_mode <- TRUE

shape_names <- c("balanced", "random", "pectinate")
n_grid <- c(2000L, 5000L, 10000L, 20000L)
if (quick_mode) n_grid <- c(50L, 100L, 500L)

shape_override <- Sys.getenv("FASTPHYLOSIG_STAGE2B2B2_SHAPES", "")
if (nzchar(shape_override)) {
  shape_names <- unique(trimws(strsplit(shape_override, ",", fixed = TRUE)[[1L]]))
  shape_names <- shape_names[nzchar(shape_names)]
  if (!length(shape_names) || any(!shape_names %in% c("balanced", "random", "pectinate"))) {
    stop("FASTPHYLOSIG_STAGE2B2B2_SHAPES must contain balanced, random, and/or pectinate.",
         call. = FALSE)
  }
}

n_override <- Sys.getenv("FASTPHYLOSIG_STAGE2B2B2_N_GRID", "")
if (nzchar(n_override)) {
  n_grid <- suppressWarnings(as.integer(trimws(strsplit(n_override, ",", fixed = TRUE)[[1L]])))
  if (!length(n_grid) || anyNA(n_grid) || any(n_grid < 2L)) {
    stop("FASTPHYLOSIG_STAGE2B2B2_N_GRID must contain integers >= 2.", call. = FALSE)
  }
  n_grid <- unique(n_grid)
}

default_repeats <- if (quick_mode) 2L else 5L
repeats <- as.integer(Sys.getenv("FASTPHYLOSIG_STAGE2B2B2_REPEATS", default_repeats))
if (!is.finite(repeats) || repeats < 1L) repeats <- default_repeats
repeats <- max(1L, repeats)
repeat_schedule_text <- Sys.getenv("FASTPHYLOSIG_STAGE2B2B2_REPEATS_BY_N", "")
repeat_schedule <- numeric()
if (nzchar(repeat_schedule_text)) {
  pieces <- trimws(strsplit(repeat_schedule_text, ",", fixed = TRUE)[[1L]])
  pairs <- strsplit(pieces, ":", fixed = TRUE)
  if (any(lengths(pairs) != 2L)) {
    stop("FASTPHYLOSIG_STAGE2B2B2_REPEATS_BY_N must use n:repeats pairs.",
         call. = FALSE)
  }
  schedule_n <- suppressWarnings(vapply(pairs, function(x) as.numeric(x[[1L]]), numeric(1)))
  schedule_r <- suppressWarnings(vapply(pairs, function(x) as.numeric(x[[2L]]), numeric(1)))
  if (any(!is.finite(schedule_n)) || any(schedule_n < 2) ||
      any(!is.finite(schedule_r)) || any(schedule_r < 1) ||
      any(schedule_r != floor(schedule_r))) {
    stop("FASTPHYLOSIG_STAGE2B2B2_REPEATS_BY_N contains invalid n:repeats values.",
         call. = FALSE)
  }
  names(schedule_r) <- as.character(as.integer(schedule_n))
  repeat_schedule <- schedule_r
}
repeats_for_n <- function(n) {
  if (length(repeat_schedule) && as.character(n) %in% names(repeat_schedule)) {
    as.integer(repeat_schedule[[as.character(n)]])
  } else {
    as.integer(repeats)
  }
}

elapsed_seconds <- function() unname(proc.time()[["elapsed"]])
time_call <- function(fun) {
  started <- elapsed_seconds()
  status <- "ok"
  error_class <- ""
  error_message <- ""
  tryCatch(
    invisible(fun()),
    error = function(e) {
      status <<- "error"
      error_class <<- paste(class(e), collapse = "/")
      error_message <<- conditionMessage(e)
      invisible(NULL)
    }
  )
  elapsed <- max(0, elapsed_seconds() - started)
  gc(FALSE)
  list(elapsed_ms = 1000 * elapsed, status = status,
       error_class = error_class, error_message = error_message)
}
clone_fixture <- function(tree) unserialize(serialize(tree, NULL))

# Load source without modifying the repository. An installed-library override
# is supported for ASCII staging environments where pkgload is unavailable.
installed_lib <- Sys.getenv("FASTPHYLOSIG_STAGE2B2B2_LIB", "")
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
  stop("pkgload is required unless FASTPHYLOSIG_STAGE2B2_LIB is set.", call. = FALSE)
}
suppressPackageStartupMessages(library(ape))
ns <- asNamespace("fastphylosig")
get_internal <- function(name) get(name, envir = ns, inherits = FALSE)

source_commit <- Sys.getenv(
  "FASTPHYLOSIG_STAGE2B2B2_SOURCE_COMMIT",
  "1685d8f58dea6dfe0a3e3d519d58581b3bbae5b9"
)
if (!nzchar(source_commit)) source_commit <- NA_character_
compiler_text <- paste(
  paste0("R.version.compiler=", R.version$compiler),
  paste0("CC=", Sys.getenv("CC", unset = "<unset>")),
  paste0("CXX=", Sys.getenv("CXX", unset = "<unset>")),
  paste0("CXX14=", Sys.getenv("CXX14", unset = "<unset>")),
  sep = "; "
)
openmp_capability <- if ("OpenMP" %in% names(capabilities())) {
  as.character(capabilities()[["OpenMP"]])
} else {
  "unavailable"
}
provenance <- data.frame(
  key = c(
    "source_commit", "source_commit_role", "source_version", "R_version",
    "R_platform", "R_arch", "R_home", "compiler", "OpenMP_capability",
    "package_load", "protocol", "quick_mode", "repeats_default",
    "repeat_schedule", "n_grid", "shapes", "working_directory"
  ),
  value = c(
    source_commit, "Candidate-1 frozen production baseline",
    as.character(utils::packageVersion("fastphylosig")), R.version.string,
    R.version$platform, R.version$arch, R.home(), compiler_text,
    openmp_capability, package_load,
    "direct private-helper and prepare_tree timing; fixed fixtures; warmup; alternating serialized copies; median/IQR",
    as.character(quick_mode), as.character(repeats),
    if (nzchar(repeat_schedule_text)) repeat_schedule_text else "<fixed>",
    paste(n_grid, collapse = ","), paste(shape_names, collapse = ","), getwd()
  ),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(out_dir, "provenance.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# Deterministic topology fixtures. These are the same construction rules used
# for the Stage 2A shape audit; branch lengths are fixed at one.

tip_labels <- function(n) paste0("sp", seq_len(n))

make_balanced <- function(n) {
  n <- as.integer(n)
  if (n < 2L) stop("n must be at least two.", call. = FALSE)
  edges <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  edge_i <- 0L
  root <- n + 1L
  next_internal <- root + 1L
  queue_node <- integer(n - 1L)
  queue_lo <- integer(n - 1L)
  queue_hi <- integer(n - 1L)
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
        edges <- rbind(edges, c(current, tip), c(current, next_internal))
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
  ape::reorder.phylo(tree, order = "postorder")
}

make_random <- function(n) {
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
         stop("unknown shape", call. = FALSE))
}

root_of <- function(tree) {
  edge <- tree$edge
  roots <- setdiff(unique(edge[, 1L]), unique(edge[, 2L]))
  if (length(roots) == 1L) as.integer(roots[[1L]]) else NA_integer_
}

helper_inputs <- function(tree) {
  edge <- tree$edge
  n_tip <- length(tree$tip.label)
  root <- root_of(tree)
  children <- get_internal(".edge_children_index")(edge)
  nodes <- sort(unique(c(edge[, 1L], edge[, 2L])))
  nodes <- nodes[nodes > n_tip]
  nodes <- nodes[is.finite(root) & nodes != root]
  list(
    descendant_keys_iterative = function() get_internal(
      ".descendant_keys_iterative")(nodes, children,
                                     as.character(tree$tip.label), n_tip),
    safe_canonicalize_core = function() get_internal(
      ".safe_canonicalize_core")(tree),
    prepare_tree = function() prepare_tree(tree)
  )
}

validate_fixture <- function(tree, shape, n) {
  stopifnot(inherits(tree, "phylo"), ape::Ntip(tree) == n,
            ape::Nnode(tree) == n - 1L, nrow(tree$edge) == 2L * n - 2L,
            length(tree$edge.length) == nrow(tree$edge),
            all(is.finite(tree$edge.length)), all(tree$edge.length > 0),
            !anyDuplicated(tree$tip.label))
  check <- get_internal(".inspect_tree_core")(
    tree, signal = c("K", "lambda", "D", "Delta")
  )
  if (!isTRUE(check$tree_summary$valid)) {
    stop(sprintf("%s n=%d fixture did not pass generic inspection.", shape, n),
         call. = FALSE)
  }
  invisible(tree)
}

stage_message <- function(text) message("[stage2b2b2 baseline] ", text)
stage_message(sprintf("starting; commit=%s; n=%s; shapes=%s; repeats=%d; quick=%s",
                      source_commit, paste(n_grid, collapse = ","),
                      paste(shape_names, collapse = ","), repeats, quick_mode))

fixture_rows <- list()
timing_rows <- list()
fixture_id <- 0L
timing_id <- 0L

for (shape in shape_names) {
  for (n in n_grid) {
    stage_message(sprintf("constructing %s n=%d", shape, n))
    tree_a <- make_fixture(shape, n)
    validate_fixture(tree_a, shape, n)
    tree_b <- clone_fixture(tree_a)
    inputs_a <- helper_inputs(tree_a)
    inputs_b <- helper_inputs(tree_b)
    fixture_id <- fixture_id + 1L
    fixture_rows[[fixture_id]] <- data.frame(
      shape = shape, n = n, root = root_of(tree_a),
      n_tip = ape::Ntip(tree_a), n_internal = ape::Nnode(tree_a),
      n_edge = nrow(tree_a$edge), valid = TRUE, stringsAsFactors = FALSE
    )

    # One untimed warmup per direct helper and fixture copy. The warmup result
    # is discarded before formal timing and cannot enter the evidence tables.
    stage_message(sprintf("warming %s n=%d", shape, n))
    for (helper in names(inputs_a)) invisible(time_call(inputs_a[[helper]]))
    for (helper in names(inputs_b)) invisible(time_call(inputs_b[[helper]]))

    n_repeats <- repeats_for_n(n)
    stage_message(sprintf("timing %s n=%d (%d repeats, serial)",
                          shape, n, n_repeats))
    for (rep in seq_len(n_repeats)) {
      copy_name <- if (rep %% 2L) "a" else "b"
      inputs <- if (identical(copy_name, "a")) inputs_a else inputs_b
      helper_order <- if (rep %% 2L) {
        c("safe_canonicalize_core", "descendant_keys_iterative", "prepare_tree")
      } else {
        c("prepare_tree", "descendant_keys_iterative", "safe_canonicalize_core")
      }
      for (position in seq_along(helper_order)) {
        helper <- helper_order[[position]]
        timed <- time_call(inputs[[helper]])
        timing_id <- timing_id + 1L
        timing_rows[[timing_id]] <- data.frame(
          shape = shape, n = n, `repeat` = rep, position = position,
          fixture_copy = copy_name, helper = helper,
          elapsed_ms = timed$elapsed_ms, status = timed$status,
          error_class = timed$error_class, error_message = timed$error_message,
          stringsAsFactors = FALSE
        )
        if (!identical(timed$status, "ok")) {
          stop(sprintf("timed %s %s n=%d failed: %s", shape, helper, n,
                       timed$error_message), call. = FALSE)
        }
      }
    }
    rm(tree_a, tree_b, inputs_a, inputs_b)
    gc(FALSE)
  }
}

timings <- do.call(rbind, timing_rows)
fixtures <- do.call(rbind, fixture_rows)
write.csv(timings, file.path(out_dir, "raw_timings.csv"), row.names = FALSE)
write.csv(fixtures, file.path(out_dir, "fixture_manifest.csv"), row.names = FALSE)

summarise <- function(data, groups) {
  key <- interaction(data[groups], drop = TRUE, lex.order = TRUE)
  pieces <- split(data, key, drop = TRUE)
  out <- lapply(pieces, function(z) {
    ok <- identical(z$status, "ok") | z$status == "ok"
    values <- z$elapsed_ms[ok & is.finite(z$elapsed_ms)]
    row <- as.list(z[1L, groups, drop = FALSE])
    row$repeats <- nrow(z)
    row$successful_repeats <- length(values)
    row$failed_repeats <- nrow(z) - length(values)
    row$median_ms <- if (length(values)) stats::median(values) else NA_real_
    row$IQR_ms <- if (length(values)) stats::IQR(values) else NA_real_
    row$min_ms <- if (length(values)) min(values) else NA_real_
    row$max_ms <- if (length(values)) max(values) else NA_real_
    as.data.frame(row, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
summary <- summarise(timings, c("shape", "n", "helper"))
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)

empirical_exponent <- function(summary_data) {
  rows <- list()
  id <- 0L
  for (shape in unique(summary_data$shape)) {
    for (helper in unique(summary_data$helper)) {
      z <- summary_data[summary_data$shape == shape &
                        summary_data$helper == helper, , drop = FALSE]
      z <- z[order(z$n), , drop = FALSE]
      z <- z[is.finite(z$median_ms) & z$median_ms > 0, , drop = FALSE]
      if (nrow(z) < 2L) next
      for (i in seq_len(nrow(z) - 1L)) {
        id <- id + 1L
        rows[[id]] <- data.frame(
          shape = shape, helper = helper,
          n_from = z$n[[i]], n_to = z$n[[i + 1L]],
          elapsed_from_ms = z$median_ms[[i]],
          elapsed_to_ms = z$median_ms[[i + 1L]],
          empirical_p = log(z$median_ms[[i + 1L]] / z$median_ms[[i]]) /
            log(z$n[[i + 1L]] / z$n[[i]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}
write.csv(empirical_exponent(summary), file.path(out_dir, "empirical_exponent.csv"),
          row.names = FALSE)

manifest <- data.frame(
  file = sort(list.files(out_dir, full.names = FALSE)),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)

stage_message("completed")
print(summary)
