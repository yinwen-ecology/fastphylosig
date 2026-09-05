#!/usr/bin/env Rscript

# fastphylosig 0.2.0 Stage 2B2B Candidate 1 small-resolution supplement.
#
# This read-only harness repeats the n=500 direct inspection enough times per
# timed block to make the Windows wall-clock resolution negligible. Namespace
# binding changes, fixture cloning, and old/new equality checks are outside
# the timer. The run is intentionally one serial R process.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(index, default) {
  if (length(args) >= index && nzchar(args[[index]])) args[[index]] else default
}

repo <- normalizePath(arg_value(1L, getwd()), winslash = "/", mustWork = TRUE)
out_dir <- arg_value(2L, file.path(tempdir(), "fastphylosig-stage2b2b1-small"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")
options(stringsAsFactors = FALSE)

prefix <- "FASTPHYLOSIG_STAGE2B2B1_SMALL_"
positive_integer <- function(text, what) {
  value <- suppressWarnings(as.integer(text))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(sprintf("%s must be a positive integer.", what), call. = FALSE)
  }
  value
}

batch_calls <- positive_integer(
  Sys.getenv(paste0(prefix, "BATCH_CALLS"), "50"),
  paste0(prefix, "BATCH_CALLS")
)
paired_blocks <- positive_integer(
  Sys.getenv(paste0(prefix, "BLOCKS"), "24"),
  paste0(prefix, "BLOCKS")
)
if (paired_blocks < 20L) {
  stop("BLOCKS must be at least 20 for the supplemental benchmark.",
       call. = FALSE)
}
n <- 500L
shape_names <- c("balanced", "random", "pectinate")

# Load the candidate source from the requested ASCII staging checkout, or
# from an explicitly supplied installed library.
installed_lib <- Sys.getenv("FASTPHYLOSIG_STAGE2B2B1_LIB", "")
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
current_inspector <- get(".inspect_tree_core", envir = ns, inherits = FALSE)
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

clone_fixture <- function(tree) unserialize(serialize(tree, NULL, version = 2L))
elapsed_seconds <- function() unname(proc.time()[["elapsed"]])

capture_call <- function(fun) {
  warning_records <- list()
  status <- "ok"
  error_message <- ""
  error_class <- ""
  value <- tryCatch(
    withCallingHandlers(
      fun(),
      warning = function(condition) {
        warning_records[[length(warning_records) + 1L]] <<- list(
          message = conditionMessage(condition), class = class(condition)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      status <<- "error"
      error_message <<- conditionMessage(condition)
      error_class <<- paste(class(condition), collapse = "/")
      NULL
    }
  )
  list(value = value, status = status, warnings = warning_records,
       error_message = error_message, error_class = error_class)
}

# Install one inspector before a call and restore the exact original binding
# afterwards. This is deliberately outside every timed batch interval.
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

invoke_once <- function(tree, implementation) {
  # The clone is made before binding setup. No timer surrounds either step.
  tree_copy <- clone_fixture(tree)
  restore <- install_inspector(implementation)
  on.exit(restore(), add = TRUE)
  inspector <- get(".inspect_tree_core", envir = ns, inherits = FALSE)
  result <- capture_call(function() inspector(tree_copy, signal = signals))
  restore()
  result
}

check_exact_equality <- function(tree, shape, block) {
  old <- invoke_once(tree, "old")
  new <- invoke_once(tree, "new")
  same_conditions <- identical(old$status, new$status) &&
    identical(old$warnings, new$warnings) &&
    identical(old$error_message, new$error_message) &&
    identical(old$error_class, new$error_class)
  if (!same_conditions || !identical(old$value, new$value) ||
      !identical(old$status, "ok")) {
    stop(sprintf(
      "exact old/new equality failed before %s block %d.", shape, block
    ), call. = FALSE)
  }
  old$value <- NULL
  new$value <- NULL
  rm(old, new)
  invisible(TRUE)
}

timed_batch <- function(tree, implementation, batch_size) {
  # Clone and install before the timer. The inspector closure is captured
  # after installation, then exactly batch_size calls are timed.
  tree_copy <- clone_fixture(tree)
  restore <- install_inspector(implementation)
  on.exit(restore(), add = TRUE)
  inspector <- get(".inspect_tree_core", envir = ns, inherits = FALSE)
  warning_count <- 0L
  status <- "ok"
  error_message <- ""
  error_class <- ""
  start <- elapsed_seconds()
  tryCatch(
    withCallingHandlers(
      for (call_id in seq_len(batch_size)) {
        invisible(inspector(tree_copy, signal = signals))
      },
      warning = function(condition) {
        warning_count <<- warning_count + 1L
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      status <<- "error"
      error_message <<- conditionMessage(condition)
      error_class <<- paste(class(condition), collapse = "/")
    }
  )
  elapsed_ms <- 1000 * max(0, elapsed_seconds() - start)
  restore()
  list(status = status, elapsed_ms = elapsed_ms,
       per_call_ms = elapsed_ms / batch_size,
       warning_count = warning_count, error_message = error_message,
       error_class = error_class)
}

# -------------------------------------------------------------------------
# Deterministic fixed n=500 fixtures. These definitions match the formal
# Stage 2B2B1 harness: unit positive branch lengths and fixed random seed.

tip_labels <- function(size) paste0("sp", seq_len(size))

make_balanced <- function(size) {
  size <- as.integer(size)
  n_internal <- size - 1L
  edges <- matrix(0L, nrow = 2L * n_internal, ncol = 2L)
  queue_node <- integer(n_internal)
  queue_lo <- integer(n_internal)
  queue_hi <- integer(n_internal)
  root <- size + 1L
  next_internal <- root + 1L
  head <- 1L
  tail <- 1L
  queue_node[[1L]] <- root
  queue_lo[[1L]] <- 1L
  queue_hi[[1L]] <- size
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
  tree <- list(edge = edges[seq_len(edge_i), , drop = FALSE],
               tip.label = tip_labels(size), edge.length = rep(1, edge_i),
               Nnode = n_internal, root = root)
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

make_pectinate <- function(size) {
  size <- as.integer(size)
  root <- size + 1L
  edges <- matrix(0L, nrow = 2L * (size - 1L), ncol = 2L)
  if (size == 2L) {
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
    if (size > 3L) {
      for (tip in 2:(size - 2L)) {
        edges[row, ] <- c(current, tip)
        row <- row + 1L
        edges[row, ] <- c(current, next_internal)
        row <- row + 1L
        current <- next_internal
        next_internal <- next_internal + 1L
      }
    }
    edges[row, ] <- c(current, size - 1L)
    row <- row + 1L
    edges[row, ] <- c(current, size)
  }
  tree <- list(edge = edges, tip.label = tip_labels(size),
               edge.length = rep(1, nrow(edges)), Nnode = size - 1L,
               root = root)
  class(tree) <- "phylo"
  tree
}

make_random <- function(size) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(2072000L + as.integer(size))
  tree <- ape::rtree(size, tip.label = tip_labels(size))
  tree <- ape::reorder.phylo(tree, order = "postorder")
  tree$edge.length <- rep(1, nrow(tree$edge))
  tree
}

make_fixture <- function(shape) {
  switch(shape,
         balanced = make_balanced(n),
         random = make_random(n),
         pectinate = make_pectinate(n),
         stop("unknown fixture shape", call. = FALSE))
}

validate_fixture <- function(tree, shape) {
  if (!inherits(tree, "phylo") || length(tree$tip.label) != n ||
      as.integer(tree$Nnode) != n - 1L || nrow(tree$edge) != 2L * (n - 1L) ||
      length(tree$edge.length) != nrow(tree$edge) ||
      any(!is.finite(tree$edge.length)) || any(tree$edge.length <= 0) ||
      anyNA(tree$tip.label) || anyDuplicated(tree$tip.label)) {
    stop(sprintf("invalid %s n=%d fixture.", shape, n), call. = FALSE)
  }
  check <- current_inspector(tree, signal = signals)
  if (!isTRUE(check$tree_summary$valid)) {
    stop(sprintf("%s n=%d fixture failed candidate inspection.", shape, n),
         call. = FALSE)
  }
  invisible(tree)
}

stage_message <- function(text) message("[stage2b2b1-small] ", text)
stage_message(sprintf("starting; n=%d; batch_calls=%d; blocks=%d",
                      n, batch_calls, paired_blocks))

fixture_cache <- new.env(parent = emptyenv())
for (shape in shape_names) {
  stage_message(sprintf("constructing %s n=%d", shape, n))
  tree <- make_fixture(shape)
  validate_fixture(tree, shape)
  fixture_cache[[shape]] <- tree
}

raw_rows <- list()
raw_id <- 0L
for (shape in shape_names) {
  tree <- fixture_cache[[shape]]
  for (block in seq_len(paired_blocks)) {
    # Equality is checked before each block and is not included in either
    # timed batch. Any mismatch aborts rather than producing ambiguous data.
    equality_ok <- check_exact_equality(tree, shape, block)
    order <- if (block %% 2L) c("old", "new") else c("new", "old")
    for (implementation in order) {
      result <- timed_batch(tree, implementation, batch_calls)
      raw_id <- raw_id + 1L
      raw_rows[[raw_id]] <- data.frame(
        shape = shape, n = n, block = as.integer(block),
        order = paste(order, collapse = ">"),
        implementation = implementation, equality_ok = equality_ok,
        batch_calls = batch_calls, elapsed_ms = result$elapsed_ms,
        per_call_ms = result$per_call_ms, status = result$status,
        warning_count = result$warning_count,
        error_class = result$error_class, error_message = result$error_message,
        stringsAsFactors = FALSE
      )
      if (!identical(result$status, "ok")) {
        stop(sprintf("timed %s %s block %d failed: %s.",
                     shape, implementation, block, result$error_message),
             call. = FALSE)
      }
      rm(result)
    }
    gc(FALSE)
  }
}

raw <- do.call(rbind, raw_rows)

summary_rows <- list()
for (shape in shape_names) {
  z <- raw[raw$shape == shape, , drop = FALSE]
  old <- z[z$implementation == "old", , drop = FALSE]
  new <- z[z$implementation == "new", , drop = FALSE]
  old <- old[order(old$block), , drop = FALSE]
  new <- new[order(new$block), , drop = FALSE]
  if (nrow(old) != paired_blocks || nrow(new) != paired_blocks ||
      !identical(old$block, new$block) || any(!raw$equality_ok)) {
    stop(sprintf("incomplete or non-equal paired data for %s.", shape),
         call. = FALSE)
  }
  old_total <- old$elapsed_ms
  new_total <- new$elapsed_ms
  old_median <- stats::median(old_total)
  new_median <- stats::median(new_total)
  ratio <- if (old_median > 0) new_median / old_median else NA_real_
  slowdown <- if (is.finite(ratio)) 100 * (ratio - 1) else NA_real_
  pair_ratio <- new_total / old_total
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    shape = shape, n = n, paired_blocks = paired_blocks,
    batch_calls = batch_calls,
    old_median_block_ms = old_median,
    old_iqr_block_ms = stats::IQR(old_total),
    new_median_block_ms = new_median,
    new_iqr_block_ms = stats::IQR(new_total),
    old_median_per_call_ms = stats::median(old$per_call_ms),
    old_iqr_per_call_ms = stats::IQR(old$per_call_ms),
    new_median_per_call_ms = stats::median(new$per_call_ms),
    new_iqr_per_call_ms = stats::IQR(new$per_call_ms),
    old_over_new_speedup = if (is.finite(ratio) && ratio > 0) 1 / ratio else NA_real_,
    new_over_old = ratio, median_slowdown_pct = slowdown,
    median_pair_slowdown_pct = stats::median(100 * (pair_ratio - 1)),
    new_median_slowdown_gt5 = is.finite(slowdown) && slowdown > 5,
    gate_status = if (!is.finite(slowdown)) "INSUFFICIENT" else
      if (slowdown > 5) "FAIL" else "PASS",
    stringsAsFactors = FALSE
  )
}
summary <- do.call(rbind, summary_rows)

gate <- data.frame(
  n = n, batch_calls = batch_calls, paired_blocks = paired_blocks,
  threshold_pct = 5, any_shape_new_median_slowdown_gt5 = any(
    summary$new_median_slowdown_gt5
  ), status = if (any(summary$new_median_slowdown_gt5)) "FAIL" else "PASS",
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------------------
# Provenance. R.version$compiler is not populated on this Windows build, so
# also retain the installed Makeconf compiler identity and C++ command.

git_capture <- function(git_args) {
  value <- tryCatch(
    system2("git", c("-C", repo, git_args), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  value <- trimws(paste(value, collapse = "\n"))
  if (nzchar(value)) value else NA_character_
}
source_commit <- Sys.getenv("FASTPHYLOSIG_STAGE2B2B1_SOURCE_COMMIT", "")
if (!nzchar(source_commit)) source_commit <- git_capture(c("rev-parse", "HEAD"))
baseline_commit <- Sys.getenv("FASTPHYLOSIG_STAGE2B2B1_BASELINE_COMMIT", "")
if (!nzchar(baseline_commit)) baseline_commit <- git_capture(
  c("rev-parse", "4dc4fa3^{commit}")
)
if (is.na(baseline_commit) || !nzchar(baseline_commit)) {
  baseline_commit <- "4dc4fa3b4eb6a0815789bdbd5276ec3063dcbda4"
}
makeconf_candidates <- c(
  file.path(R.home("etc"), "x64", "Makeconf"),
  file.path(R.home("etc"), "Makeconf")
)
makeconf_file <- makeconf_candidates[file.exists(makeconf_candidates)][1L]
makeconf_lines <- if (length(makeconf_file) && !is.na(makeconf_file)) {
  readLines(makeconf_file, warn = FALSE)
} else character()
makeconf_value <- function(key) {
  hit <- grep(paste0("^", key, "\\s*="), makeconf_lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub(paste0("^", key, "\\s*=\\s*"), "", hit[[1L]]))
}
scalar_or_na <- function(value) {
  if (!length(value) || is.na(value[[1L]]) || !nzchar(as.character(value[[1L]]))) {
    return(NA_character_)
  }
  as.character(value[[1L]])
}
ascii_only <- function(value) !grepl("[^ -~]", value)
md5_or_na <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(as.character(tools::md5sum(path)))
}

provenance <- data.frame(
  key = c(
    "source_commit", "baseline_oracle_commit", "check_tree_md5",
    "oracle_fixture_md5", "R_version", "R_platform", "R_arch", "R_home",
    "R_compiler", "R_compiled_by", "R_makeconf_CC",
    "R_makeconf_CXX_default", "R_makeconf_CXX20STD",
    "package_declared_system_requirements",
    "package_version", "ape_version", "package_load", "repo_path_ascii",
    "output_path_ascii", "protocol", "n", "batch_calls", "paired_blocks"
  ),
  value = c(
    source_commit, baseline_commit,
    md5_or_na(file.path(repo, "R", "check_tree.R")), md5_or_na(oracle_file),
    R.version$version.string, R.version$platform, R.version$arch, R.home(),
    scalar_or_na(R.version$compiler), makeconf_value("COMPILED_BY"),
    makeconf_value("CC"), makeconf_value("CXX"), makeconf_value("CXX20STD"),
    scalar_or_na(read.dcf(file.path(repo, "DESCRIPTION"))[
      1L, "SystemRequirements"
    ]),
    as.character(utils::packageVersion("fastphylosig")),
    as.character(utils::packageVersion("ape")), package_load,
    as.character(ascii_only(repo)), as.character(ascii_only(out_dir)),
    "fixed balanced/random/pectinate n=500; unit positive branches; exact equality before every block; fresh clone/setup outside timer; batch calls per timed block; alternating order; serial",
    as.character(n), as.character(batch_calls), as.character(paired_blocks)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(raw, file.path(out_dir, "raw_blocks.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
utils::write.csv(gate, file.path(out_dir, "slowdown_gate.csv"), row.names = FALSE)
utils::write.csv(provenance, file.path(out_dir, "provenance.csv"), row.names = FALSE)
manifest <- data.frame(
  file = sort(list.files(out_dir, full.names = FALSE)), stringsAsFactors = FALSE
)
utils::write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)

stage_message(sprintf("completed; output=%s", out_dir))
print(summary)
print(gate)
