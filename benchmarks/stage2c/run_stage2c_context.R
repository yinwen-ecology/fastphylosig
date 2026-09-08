#!/usr/bin/env Rscript

# fastphylosig 0.2.0 Stage 2C context memory/persistence harness.
#
# This file is an audit harness. It does not modify package source, tests, or
# package state. The default mode is a small smoke run. The formal context
# grid is opt-in with --formal because it includes a 20,000-tip object.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_candidates <- c(
  Sys.getenv("FASTPHYLOSIG_STAGE2C_SCRIPT", unset = ""),
  if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else "",
  file.path(getwd(), "benchmarks", "stage2c", "run_stage2c_context.R")
)
script_path <- NA_character_
for (candidate in script_candidates) {
  candidate <- sub('^"|"$', "", as.character(candidate[[1L]]))
  if (!nzchar(candidate) || nchar(candidate, type = "bytes") > 4096L) next
  exists <- tryCatch(file.exists(candidate), error = function(e) FALSE)
  if (isTRUE(exists)) {
    script_path <- tryCatch(
      normalizePath(candidate, winslash = "/", mustWork = TRUE),
      error = function(e) candidate
    )
    break
  }
}

cli <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript run_stage2c_context.R LIBRARY OUTPUT_DIR [--smoke]",
      "Rscript run_stage2c_context.R LIBRARY OUTPUT_DIR --formal",
      "",
      "The default is --smoke. --formal enables n=500,1000,5000,10000,20000.",
      sep = "\n"
    ),
    "\n",
    sep = ""
  )
}

if ("--help" %in% cli || "-h" %in% cli) {
  usage()
  quit(status = 0L)
}

# The child is deliberately implemented in this same file. That keeps the
# fresh-process protocol self-contained and avoids a second helper script.
if ("--child" %in% cli) {
  child_library <- Sys.getenv("FASTPHYLOSIG_STAGE2C_CHILD_LIBRARY", unset = "")
  child_rds <- Sys.getenv("FASTPHYLOSIG_STAGE2C_CHILD_RDS", unset = "")
  child_result <- Sys.getenv("FASTPHYLOSIG_STAGE2C_CHILD_RESULT", unset = "")
  if (!nzchar(child_library) || !nzchar(child_rds) ||
      !nzchar(child_result)) {
    stop(
      paste0(
        "child mode requires FASTPHYLOSIG_STAGE2C_CHILD_LIBRARY, ",
        "FASTPHYLOSIG_STAGE2C_CHILD_RDS, and ",
        "FASTPHYLOSIG_STAGE2C_CHILD_RESULT."
      ),
      call. = FALSE
    )
  }
  child_timeout <- suppressWarnings(as.numeric(Sys.getenv(
    "FASTPHYLOSIG_STAGE2C_CHILD_TIMEOUT_SECONDS", unset = "120"
  )))
  if (length(child_timeout) != 1L || is.na(child_timeout) ||
      !is.finite(child_timeout) || child_timeout <= 0) {
    stop("FASTPHYLOSIG_STAGE2C_CHILD_TIMEOUT_SECONDS must be positive.",
         call. = FALSE)
  }
  # The context child is pure R read/validation work. A transient R-level
  # elapsed limit guarantees that a malformed or unexpectedly slow child
  # returns a recorded result instead of leaving a live process behind.
  setTimeLimit(cpu = Inf, elapsed = child_timeout, transient = TRUE)
  result <- tryCatch({
    .libPaths(unique(c(child_library, .libPaths())))
    suppressPackageStartupMessages(
      library("fastphylosig", lib.loc = child_library, character.only = TRUE)
    )
    started <- unname(proc.time()[["elapsed"]])
    ctx <- readRDS(child_rds)
    after_read <- unname(proc.time()[["elapsed"]])
    info <- fastphylosig::cache_info(ctx)
    after_validation <- unname(proc.time()[["elapsed"]])
    if (!is.list(info) || !is.data.frame(info$entries)) {
      stop("cache_info() did not return the expected context report.")
    }
    list(
      status = "ok",
      read_rds_ms = 1000 * (after_read - started),
      first_validation_ms = 1000 * (after_validation - after_read),
      read_plus_first_validation_ms = 1000 * (after_validation - started),
      pid = Sys.getpid(),
      context_class = paste(class(ctx), collapse = "/"),
      n_structural_entries = as.integer(info$n_structural_entries),
      n_numerical_entries = as.integer(info$n_numerical_entries),
      timeout_seconds = child_timeout,
      validation_message = "cache_info() completed before any estimator call",
      error_message = ""
    )
  }, error = function(e) {
    list(
      status = if (inherits(e, "elapsedTimeLimit")) "timeout" else "error",
      read_rds_ms = NA_real_,
      first_validation_ms = NA_real_,
      read_plus_first_validation_ms = NA_real_,
      pid = Sys.getpid(),
      context_class = "",
      n_structural_entries = NA_integer_,
      n_numerical_entries = NA_integer_,
      timeout_seconds = child_timeout,
      validation_message = "",
      error_message = conditionMessage(e)
    )
  })
  # Reset before writing the status file so timeout handling itself cannot be
  # interrupted. The child then exits synchronously after saveRDS().
  setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  dir.create(dirname(child_result), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, child_result, version = 3L, compress = FALSE)
  if (identical(result$status, "ok")) {
    cat("STAGE2C_CHILD=PASS\n")
    quit(status = 0L)
  }
  cat("STAGE2C_CHILD=FAIL\n")
  quit(status = 1L)
}

formal_mode <- "--formal" %in% cli
positional <- cli[!cli %in% c("--formal", "--smoke")]
if (length(positional) < 2L) {
  usage()
  stop("LIBRARY and OUTPUT_DIR are required.", call. = FALSE)
}

package_library <- normalizePath(positional[[1L]], winslash = "/",
                                  mustWork = TRUE)
output_dir <- normalizePath(positional[[2L]], winslash = "/",
                             mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(file.path(package_library, "fastphylosig", "DESCRIPTION"))) {
  stop("LIBRARY must contain the installed fastphylosig package.", call. = FALSE)
}

options(stringsAsFactors = FALSE)
.libPaths(unique(c(package_library, .libPaths())))
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("ape is required for the deterministic tree fixture.", call. = FALSE)
}
suppressPackageStartupMessages(
  library("fastphylosig", lib.loc = package_library, character.only = TRUE)
)

mode <- if (formal_mode) "formal" else "smoke"
n_grid <- if (formal_mode) {
  c(500L, 1000L, 5000L, 10000L, 20000L)
} else {
  # Two small sizes exercise the complete route without accidentally turning
  # a smoke invocation into the formal grid.
  c(8L, 32L)
}
n_override <- Sys.getenv("FASTPHYLOSIG_STAGE2C_N_GRID", unset = "")
if (nzchar(n_override)) {
  if (!formal_mode) {
    stop(
      "FASTPHYLOSIG_STAGE2C_N_GRID is accepted only together with --formal.",
      call. = FALSE
    )
  }
  parsed <- suppressWarnings(as.integer(
    trimws(strsplit(n_override, ",", fixed = TRUE)[[1L]])
  ))
  if (!length(parsed) || anyNA(parsed) || any(parsed < 2L)) {
    stop("FASTPHYLOSIG_STAGE2C_N_GRID must contain integers >= 2.",
         call. = FALSE)
  }
  n_grid <- unique(parsed)
}

as_flag <- function(value, default = FALSE) {
  if (!nzchar(value)) return(default)
  key <- tolower(trimws(value))
  if (key %in% c("1", "true", "yes", "y")) return(TRUE)
  if (key %in% c("0", "false", "no", "n")) return(FALSE)
  stop("expected a boolean environment value, got: ", value, call. = FALSE)
}

keep_rds <- as_flag(
  Sys.getenv("FASTPHYLOSIG_STAGE2C_KEEP_RDS", unset = ""),
  default = FALSE
)
child_timeout_seconds <- suppressWarnings(as.numeric(Sys.getenv(
  "FASTPHYLOSIG_STAGE2C_CHILD_TIMEOUT_SECONDS", unset = "120"
)))
if (length(child_timeout_seconds) != 1L || is.na(child_timeout_seconds) ||
    !is.finite(child_timeout_seconds) || child_timeout_seconds <= 0) {
  stop(
    "FASTPHYLOSIG_STAGE2C_CHILD_TIMEOUT_SECONDS must be positive.",
    call. = FALSE
  )
}

make_balanced <- function(n) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("n must be an integer >= 2.", call. = FALSE)
  }
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
      edge[edge_i, ] <- c(node, child)
    }
  }
  edge <- edge[seq_len(edge_i), , drop = FALSE]
  tree <- list(
    edge = edge,
    tip.label = paste0("sp", seq_len(n)),
    edge.length = 0.5 + seq_len(nrow(edge)) / nrow(edge),
    Nnode = n_internal
  )
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

elapsed_ms <- function(start, end) 1000 * (end - start)

safe_scalar <- function(value, fallback = NA_real_) {
  tryCatch({
    out <- as.numeric(value)
    if (length(out) == 1L && is.finite(out)) out else fallback
  }, error = function(e) fallback)
}

cache_payload_bytes <- function(cache) {
  if (!is.environment(cache)) return(NA_real_)
  keys <- ls(cache, all.names = TRUE)
  if (!length(keys)) return(0)
  sum(vapply(keys, function(key) {
    as.numeric(utils::object.size(get(key, cache, inherits = FALSE)))
  }, numeric(1)))
}

cache_metric <- function(ctx) {
  info <- fastphylosig::cache_info(ctx)
  entries <- info$entries
  structural_reported <- if (nrow(entries)) {
    sum(entries$structural_bytes, na.rm = TRUE)
  } else 0
  numerical_reported <- if (nrow(entries)) {
    sum(entries$numerical_bytes, na.rm = TRUE)
  } else 0
  structural_payload <- cache_payload_bytes(ctx$structural_cache)
  numerical_payload <- cache_payload_bytes(ctx$numerical_cache)
  list(
    info = info,
    cache_info_bytes_used = safe_scalar(info$bytes_used),
    cache_reported_entry_bytes = structural_reported + numerical_reported,
    cache_reported_structural_bytes = structural_reported,
    cache_reported_numerical_bytes = numerical_reported,
    structural_cache_payload_bytes = structural_payload,
    numerical_cache_payload_bytes = numerical_payload,
    cache_payload_bytes = structural_payload + numerical_payload,
    cache_entries = as.integer(info$n_structural_entries +
                                info$n_numerical_entries)
  )
}

short_text <- function(value, limit = 500L) {
  value <- paste(as.character(value), collapse = " | ")
  if (nchar(value, type = "bytes") > limit) {
    paste0(substr(value, 1L, limit - 3L), "...")
  } else value
}

rscript_path <- file.path(R.home("bin"), "Rscript")
if (.Platform$OS.type == "windows" && !file.exists(rscript_path)) {
  rscript_path <- paste0(rscript_path, ".exe")
}
if (!file.exists(rscript_path)) {
  stop("could not locate Rscript for the fresh-process check.", call. = FALSE)
}

run_fresh_validation <- function(rds_path, n) {
  script_available <- !is.na(script_path) && length(script_path) == 1L &&
    nzchar(script_path) && tryCatch(file.exists(script_path),
                                     error = function(e) FALSE)
  if (!isTRUE(script_available)) {
    return(list(
      status = "not_run_script_path",
      read_rds_ms = NA_real_, first_validation_ms = NA_real_,
      read_plus_first_validation_ms = NA_real_, pid = NA_integer_,
      context_class = "", n_structural_entries = NA_integer_,
      n_numerical_entries = NA_integer_, validation_message = "",
      timeout_seconds = child_timeout_seconds, child_output = "",
      error_message = "the harness script path was not available"
    ))
  }
  result_path <- tempfile(
    sprintf("fastphylosig-stage2c-child-n%s-", n), fileext = ".rds"
  )
  # On Windows, system2(env=) can replace rather than extend the inherited
  # environment. Preserve it explicitly so R_HOME, PATH, and DLL lookup stay
  # available to the child, then override only the child protocol variables.
  child_env <- Sys.getenv()
  child_env["FASTPHYLOSIG_STAGE2C_CHILD_LIBRARY"] <- package_library
  child_env["FASTPHYLOSIG_STAGE2C_CHILD_RDS"] <- rds_path
  child_env["FASTPHYLOSIG_STAGE2C_CHILD_RESULT"] <- result_path
  child_env["FASTPHYLOSIG_STAGE2C_CHILD_TIMEOUT_SECONDS"] <-
    as.character(child_timeout_seconds)
  env <- paste(names(child_env), child_env, sep = "=")
  output <- tryCatch(
    system2(
      rscript_path,
      args = c("--vanilla", script_path, "--child"),
      stdout = TRUE,
      stderr = TRUE,
      env = env
    ),
    error = function(e) structure(
      paste0("system2 error: ", conditionMessage(e)), status = 1L
    )
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  child <- if (file.exists(result_path)) {
    tryCatch(readRDS(result_path), error = function(e) NULL)
  } else NULL
  unlink(result_path, force = TRUE)
  if (!is.list(child)) {
    return(list(
      status = if (as.integer(status) == 0L) "missing_result" else "process_error",
      read_rds_ms = NA_real_, first_validation_ms = NA_real_,
      read_plus_first_validation_ms = NA_real_, pid = NA_integer_,
      context_class = "", n_structural_entries = NA_integer_,
      n_numerical_entries = NA_integer_, validation_message = "",
      timeout_seconds = child_timeout_seconds,
      child_output = short_text(output),
      error_message = short_text(output)
    ))
  }
  child$process_status <- as.integer(status)
  child$child_output <- short_text(output)
  child
}

safe_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else "<unset>"
}

safe_version <- function(value) {
  if (length(value) && !is.null(value) && nzchar(as.character(value[[1L]]))) {
    as.character(value[[1L]])
  } else "<unavailable>"
}

resolve_commit <- function() {
  supplied <- Sys.getenv("FASTPHYLOSIG_STAGE2C_COMMIT", unset = "")
  if (nzchar(supplied)) {
    return(list(value = supplied, source = "FASTPHYLOSIG_STAGE2C_COMMIT"))
  }
  source_dir <- Sys.getenv(
    "FASTPHYLOSIG_STAGE2C_SOURCE_DIR", unset = getwd()
  )
  value <- tryCatch(
    system2("git", c("-C", source_dir, "rev-parse", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  value <- trimws(paste(value, collapse = ""))
  if (length(value) && grepl("^[0-9a-fA-F]{40}$", value)) {
    return(list(value = tolower(value), source = "git rev-parse HEAD"))
  }
  list(value = "<unavailable>", source = "not resolved")
}

commit <- resolve_commit()
session_lines <- capture.output(sessionInfo())
writeLines(session_lines, file.path(output_dir, "stage2c_sessionInfo.txt"),
           useBytes = TRUE)

soft_versions <- tryCatch(extSoftVersion(), error = function(e) character())
blas <- if (length(soft_versions) && "BLAS" %in% names(soft_versions)) {
  safe_version(soft_versions[["BLAS"]])
} else "<unavailable>"
lapack <- if (length(soft_versions) && "LAPACK" %in% names(soft_versions)) {
  safe_version(soft_versions[["LAPACK"]])
} else "<unavailable>"
capability_values <- tryCatch(capabilities(), error = function(e) logical())
openmp_capability <- if (length(capability_values) &&
                         "OpenMP" %in% names(capability_values)) {
  as.character(unname(capability_values[["OpenMP"]]))
} else "not_exposed_by_R"
openmp_flags <- paste(
  Filter(nzchar, Sys.getenv(
    c("SHLIB_OPENMP_CFLAGS", "SHLIB_OPENMP_CXXFLAGS", "OPENMP_CXXFLAGS"),
    unset = ""
  )),
  collapse = ";"
)
if (!nzchar(openmp_flags)) openmp_flags <- "<unset>"

sys <- Sys.info()
logical_cores <- tryCatch(parallel::detectCores(logical = TRUE),
                          error = function(e) NA_integer_)
physical_cores <- tryCatch(parallel::detectCores(logical = FALSE),
                           error = function(e) NA_integer_)
compiler <- if (length(R.version$compiler)) {
  paste(as.character(R.version$compiler), collapse = ";")
} else "<unavailable>"

provenance <- data.frame(
  run_mode = mode,
  run_started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  full_commit = commit$value,
  commit_resolution = commit$source,
  package_version = as.character(packageVersion("fastphylosig")),
  R_version = R.version.string,
  R_platform = R.version$platform,
  R_arch = R.version$arch,
  OS = safe_version(sys[["sysname"]]),
  OS_release = safe_version(sys[["release"]]),
  machine = safe_version(sys[["machine"]]),
  cpu_cores_logical = as.integer(logical_cores),
  cpu_cores_physical = as.integer(physical_cores),
  processor_identifier = safe_env("PROCESSOR_IDENTIFIER"),
  compiler = compiler,
  CC = safe_env("CC"),
  CXX = safe_env("CXX"),
  CXX17 = safe_env("CXX17"),
  BLAS = blas,
  LAPACK = lapack,
  OpenMP_status = "not exercised by context-only path",
  OpenMP_capability = openmp_capability,
  OpenMP_flags = openmp_flags,
  OMP_NUM_THREADS = safe_env("OMP_NUM_THREADS"),
  OMP_DYNAMIC = safe_env("OMP_DYNAMIC"),
  OMP_PROC_BIND = safe_env("OMP_PROC_BIND"),
  LC_ALL = safe_env("LC_ALL"),
  LC_COLLATE = safe_env("LC_COLLATE"),
  LANG = safe_env("LANG"),
  locale_all = safe_version(Sys.getlocale()),
  saveRDS_version = 3L,
  saveRDS_compress = FALSE,
  formal_grid = paste(n_grid, collapse = ","),
  stringsAsFactors = FALSE
)
utils::write.csv(provenance,
                 file.path(output_dir, "stage2c_provenance.csv"),
                 row.names = FALSE)

rds_dir <- if (keep_rds) {
  path <- file.path(output_dir, "context_rds")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
} else {
  tempdir()
}

rows <- vector("list", length(n_grid))
for (i in seq_along(n_grid)) {
  n <- n_grid[[i]]
  message("[stage2c-context] mode=", mode, " n=", n)
  tree <- make_balanced(n)
  prepare_started <- unname(proc.time()[["elapsed"]])
  context <- tryCatch(fastphylosig::prepare_tree(tree), error = identity)
  prepare_finished <- unname(proc.time()[["elapsed"]])
  if (inherits(context, "error")) {
    rows[[i]] <- data.frame(
      run_mode = mode, shape = "balanced", n = n,
      prepare_tree_ms = elapsed_ms(prepare_started, prepare_finished),
      prepare_status = "error", error_message = conditionMessage(context),
      context_object_bytes = NA_real_, fingerprint_bytes = NA_real_,
      protected_snapshot_bytes = NA_real_, serialized_context_bytes = NA_real_,
      cache_info_bytes_used = NA_real_, cache_reported_entry_bytes = NA_real_,
      cache_reported_structural_bytes = NA_real_,
      cache_reported_numerical_bytes = NA_real_,
      structural_cache_payload_bytes = NA_real_,
      numerical_cache_payload_bytes = NA_real_, cache_payload_bytes = NA_real_,
      cache_entries = NA_integer_, rds_path = "", saveRDS_ms = NA_real_,
      saveRDS_bytes = NA_real_, fresh_process_status = "not_run",
      readRDS_ms = NA_real_, first_validation_ms = NA_real_,
      read_plus_first_validation_ms = NA_real_, fresh_process_pid = NA_integer_,
      fresh_process_timeout_seconds = child_timeout_seconds,
      fresh_process_error = "", fresh_process_child_output = "",
      stringsAsFactors = FALSE
    )
    next
  }
  metrics <- tryCatch(cache_metric(context), error = identity)
  if (inherits(metrics, "error")) {
    stop("cache metric collection failed at n=", n, ": ",
         conditionMessage(metrics), call. = FALSE)
  }
  serialized_bytes <- tryCatch(
    as.numeric(length(serialize(context, NULL, version = 3L))),
    error = function(e) NA_real_
  )
  rds_path <- file.path(rds_dir, sprintf("stage2c_context_n%s.rds", n))
  save_started <- unname(proc.time()[["elapsed"]])
  save_error <- tryCatch({
    saveRDS(context, rds_path, version = 3L, compress = FALSE)
    NULL
  }, error = identity)
  save_finished <- unname(proc.time()[["elapsed"]])
  save_bytes <- if (is.null(save_error) && file.exists(rds_path)) {
    as.numeric(file.info(rds_path)$size)
  } else NA_real_
  fresh <- if (is.null(save_error)) {
    run_fresh_validation(rds_path, n)
  } else {
    list(
      status = "save_error", read_rds_ms = NA_real_,
      first_validation_ms = NA_real_,
      read_plus_first_validation_ms = NA_real_, pid = NA_integer_,
      timeout_seconds = child_timeout_seconds,
      error_message = conditionMessage(save_error)
    )
  }
  rows[[i]] <- data.frame(
    run_mode = mode, shape = "balanced", n = n,
    prepare_tree_ms = elapsed_ms(prepare_started, prepare_finished),
    prepare_status = "ok", error_message = "",
    context_object_bytes = as.numeric(utils::object.size(context)),
    fingerprint_bytes = as.numeric(length(context$fingerprint)),
    protected_snapshot_bytes = as.numeric(length(context$protected_snapshot)),
    serialized_context_bytes = serialized_bytes,
    cache_info_bytes_used = metrics$cache_info_bytes_used,
    cache_reported_entry_bytes = metrics$cache_reported_entry_bytes,
    cache_reported_structural_bytes = metrics$cache_reported_structural_bytes,
    cache_reported_numerical_bytes = metrics$cache_reported_numerical_bytes,
    structural_cache_payload_bytes = metrics$structural_cache_payload_bytes,
    numerical_cache_payload_bytes = metrics$numerical_cache_payload_bytes,
    cache_payload_bytes = metrics$cache_payload_bytes,
    cache_entries = metrics$cache_entries,
    rds_path = if (keep_rds) rds_path else "temporary; removed after run",
    saveRDS_ms = if (is.null(save_error)) {
      elapsed_ms(save_started, save_finished)
    } else NA_real_,
    saveRDS_bytes = save_bytes,
    fresh_process_status = as.character(fresh$status),
    readRDS_ms = safe_scalar(fresh$read_rds_ms),
    first_validation_ms = safe_scalar(fresh$first_validation_ms),
    read_plus_first_validation_ms = safe_scalar(
      fresh$read_plus_first_validation_ms
    ),
    fresh_process_pid = suppressWarnings(as.integer(fresh$pid)),
    fresh_process_timeout_seconds = safe_scalar(
      fresh$timeout_seconds, child_timeout_seconds
    ),
    fresh_process_error = if (identical(fresh$status, "ok")) "" else
      short_text(fresh$error_message),
    fresh_process_child_output = short_text(fresh$child_output),
    stringsAsFactors = FALSE
  )
  if (!keep_rds && file.exists(rds_path)) unlink(rds_path, force = TRUE)
  rm(tree, context, metrics, fresh)
  invisible(gc(FALSE))
}

memory <- do.call(rbind, rows)
utils::write.csv(memory,
                 file.path(output_dir, "stage2c_context_memory.csv"),
                 row.names = FALSE)

cat("STAGE2C_CONTEXT_HARNESS=PASS\n")
cat("MODE=", mode, "\n", sep = "")
cat("N_GRID=", paste(n_grid, collapse = ","), "\n", sep = "")
cat("OUTPUT_DIR=", output_dir, "\n", sep = "")
cat("FORMAL_GRID_EXECUTED=", if (formal_mode) "YES" else "NO", "\n", sep = "")
