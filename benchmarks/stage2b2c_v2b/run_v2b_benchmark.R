#!/usr/bin/env Rscript

# fastphylosig 0.2.0 Stage 2B2 V2 integration benchmark.
#
# This is an audit harness only.  It loads one or two installed package
# libraries in isolated callr processes, times the same fixed fixtures, and
# writes evidence below the requested output directory.  It never edits
# package source, tests, documentation, or the persistence decision.

raw_args <- commandArgs(trailingOnly = TRUE)
single_mode <- "--single" %in% raw_args
quick_mode <- "quick" %in% raw_args
args <- raw_args[!raw_args %in% c("--single", "quick")]

if (single_mode) {
  if (length(args) < 2L) {
    stop(
      "usage: run_v2b_benchmark.R --single LIBRARY OUTPUT_DIR [quick]",
      call. = FALSE
    )
  }
  before_lib <- NULL
  after_lib <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
  out_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
  single_label <- Sys.getenv("FASTPHYLOSIG_V2B_LABEL", unset = "current")
  if (!nzchar(single_label)) single_label <- "current"
} else {
  if (length(args) < 3L) {
    stop(
      paste(
        "usage: run_v2b_benchmark.R BEFORE_LIB AFTER_LIB OUTPUT_DIR [quick]",
        "or: run_v2b_benchmark.R --single LIBRARY OUTPUT_DIR [quick]"
      ),
      call. = FALSE
    )
  }
  before_lib <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
  after_lib <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
  out_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
  single_label <- NULL
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(stringsAsFactors = FALSE)

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("ape is required for the benchmark fixtures.", call. = FALSE)
}
if (!requireNamespace("callr", quietly = TRUE)) {
  stop("callr is required for isolated sequential benchmark processes.",
       call. = FALSE)
}

parse_integer_list <- function(value, name) {
  out <- suppressWarnings(as.integer(
    trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  ))
  if (!length(out) || anyNA(out) || any(out < 2L)) {
    stop(name, " must contain comma-separated integers >= 2.", call. = FALSE)
  }
  unique(out)
}

n_grid_default <- if (quick_mode) c(500L, 1000L) else
  c(500L, 1000L, 2000L, 5000L, 10000L, 20000L)
n_grid_text <- Sys.getenv("FASTPHYLOSIG_V2B_N_GRID", unset = "")
n_grid <- if (nzchar(n_grid_text)) {
  parse_integer_list(n_grid_text, "FASTPHYLOSIG_V2B_N_GRID")
} else {
  n_grid_default
}

shape_names <- c("balanced", "random", "pectinate")
shape_text <- Sys.getenv("FASTPHYLOSIG_V2B_SHAPES", unset = "")
if (nzchar(shape_text)) {
  shape_names <- unique(trimws(strsplit(shape_text, ",", fixed = TRUE)[[1L]]))
  shape_names <- shape_names[nzchar(shape_names)]
  if (!length(shape_names) ||
      any(!shape_names %in% c("balanced", "random", "pectinate"))) {
    stop(
      "FASTPHYLOSIG_V2B_SHAPES must contain balanced, random, and/or pectinate.",
      call. = FALSE
    )
  }
}

repeat_override <- suppressWarnings(as.integer(
  Sys.getenv("FASTPHYLOSIG_V2B_REPEATS", unset = "")
))
if (length(repeat_override) && !is.na(repeat_override) &&
    repeat_override > 0L) {
  repeats_for <- function(n) repeat_override
} else if (quick_mode) {
  repeats_for <- function(n) 2L
} else {
  # Large fixtures use fewer paired observations; the schedule is recorded in
  # provenance so a result is never mistaken for a full small-n repeat count.
  repeats_for <- function(n) if (n >= 10000L) 5L else 10L
}

timeout_seconds <- suppressWarnings(as.numeric(
  Sys.getenv("FASTPHYLOSIG_V2B_TIMEOUT", unset = "")
))
if (!length(timeout_seconds) || is.na(timeout_seconds) ||
    !is.finite(timeout_seconds) || timeout_seconds <= 0) {
  timeout_seconds <- if (quick_mode) Inf else 30
}

complexity_n <- suppressWarnings(as.integer(Sys.getenv(
  "FASTPHYLOSIG_V2B_COMPLEXITY_N", unset = if (quick_mode) "500" else "20000"
)))
if (length(complexity_n) != 1L || is.na(complexity_n) || complexity_n < 2L) {
  stop("FASTPHYLOSIG_V2B_COMPLEXITY_N must be an integer >= 2.",
       call. = FALSE)
}

route_names <- c(
  "canonicalization", "prepare_tree", "raw_fast_k", "prepared_fast_k",
  "fingerprint", "context_validation"
)

tip_labels <- function(n) paste0("sp", seq_len(n))

make_balanced <- function(n) {
  n <- as.integer(n)
  if (n < 2L) stop("n must be at least two.", call. = FALSE)
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
    tip.label = tip_labels(n),
    edge.length = 0.5 + seq_len(nrow(edge)) / nrow(edge),
    Nnode = n_internal
  )
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

make_pectinate <- function(n) {
  n <- as.integer(n)
  if (n < 2L) stop("n must be at least two.", call. = FALSE)
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
  stats::setNames(sin(seq_len(n)) + cos(seq_len(n) / 7), tree$tip.label)
}

clone_fixture <- function(value) unserialize(serialize(value, NULL))

package_version_at <- function(lib) {
  desc <- tryCatch(
    utils::packageDescription("fastphylosig", lib.loc = lib),
    error = function(e) NULL
  )
  if (is.null(desc) || is.null(desc$Version)) "<unknown>" else
    as.character(desc$Version)
}

commit_for <- function(version) {
  env_name <- if (identical(version, "candidate1")) {
    "FASTPHYLOSIG_V2B_BEFORE_COMMIT"
  } else if (identical(version, "v2")) {
    "FASTPHYLOSIG_V2B_AFTER_COMMIT"
  } else {
    "FASTPHYLOSIG_V2B_COMMIT"
  }
  value <- Sys.getenv(env_name, unset = "")
  if (nzchar(value)) return(value)
  switch(
    version,
    candidate1 = "7d0d6e38756f271e48a74a4ac890f64b186674f3",
    v2 = "34ae5de",
    "<unspecified>"
  )
}

libraries <- if (single_mode) {
  stats::setNames(list(after_lib), single_label)
} else {
  list(candidate1 = before_lib, v2 = after_lib)
}
versions <- names(libraries)

time_route <- function(session, route, tree, trait, timeout_seconds = Inf) {
  n <- length(tree$tip.label)
  inner_iterations <- if (n <= 1000L) {
    if (route %in% c("fingerprint", "context_validation", "prepared_fast_k")) 10L else 3L
  } else if (n <= 5000L) {
    if (route %in% c("fingerprint", "context_validation", "prepared_fast_k")) 3L else 1L
  } else {
    1L
  }
  session$run(
    function(route, tree, trait, timeout_seconds, inner_iterations) {
      ns <- asNamespace("fastphylosig")
      elapsed_seconds <- function() unname(proc.time()[["elapsed"]])
      error_result <- function(status, message, class = "") {
        list(
          elapsed_ms = NA_real_, statistic = NA_real_, status = status,
          error_message = as.character(message), error_class = as.character(class)
        )
      }
      safe_statistic <- function(value) {
        tryCatch({
          numeric_value <- as.numeric(value)
          if (length(numeric_value)) numeric_value[[1L]] else NA_real_
        }, error = function(e) NA_real_)
      }
      get_private <- function(name) {
        if (!exists(name, envir = ns, inherits = FALSE)) return(NULL)
        get(name, envir = ns, inherits = FALSE)
      }
      setup <- tryCatch({
        run <- NULL
        if (identical(route, "canonicalization")) {
          canonicalize <- get_private(".safe_canonicalize_core")
          if (!is.function(canonicalize)) {
            stop(".safe_canonicalize_core is unavailable")
          }
          run <- function() canonicalize(tree)
        } else if (identical(route, "prepare_tree")) {
          run <- function() fastphylosig::prepare_tree(tree)
        } else if (identical(route, "raw_fast_k")) {
          run <- function() fastphylosig::fast_k(
            tree, trait, test = FALSE, verbose = FALSE, progress = FALSE
          )
        } else if (identical(route, "prepared_fast_k")) {
          context <- fastphylosig::prepare_tree(tree)
          run <- function() fastphylosig::fast_k(
            context, trait, test = FALSE, verbose = FALSE, progress = FALSE
          )
        } else if (identical(route, "fingerprint")) {
          canonicalize <- get_private(".safe_canonicalize_core")
          fingerprint <- get_private(".tree_fingerprint")
          if (!is.function(canonicalize) || !is.function(fingerprint)) {
            stop("canonicalization/fingerprint helper is unavailable")
          }
          canonical <- canonicalize(tree)
          canonical_argument <- ".canonical" %in% names(formals(fingerprint))
          run <- if (canonical_argument) {
            function() fingerprint(canonical, .canonical = TRUE)
          } else {
            function() fingerprint(canonical)
          }
        } else if (identical(route, "context_validation")) {
          validate <- get_private(".validate_prepared_context")
          if (!is.function(validate)) {
            stop(".validate_prepared_context is unavailable")
          }
          context <- fastphylosig::prepare_tree(tree)
          run <- function() validate(context)
        } else {
          stop("unknown timing route: ", route)
        }
        invisible(run())
        run
      }, error = function(e) e)
      if (inherits(setup, "error")) {
        return(error_result("warmup_error", conditionMessage(setup),
                            paste(class(setup), collapse = "/")))
      }
      run <- setup
      gc(FALSE)
      limited <- is.finite(timeout_seconds) && timeout_seconds > 0
      if (limited) {
        setTimeLimit(cpu = Inf, elapsed = timeout_seconds, transient = TRUE)
      }
      status <- "ok"
      error_message <- ""
      error_class <- ""
      start <- elapsed_seconds()
      value <- tryCatch(
        {
          result <- NULL
          for (i in seq_len(inner_iterations)) result <- run()
          result
        },
        error = function(e) {
          status <<- if (inherits(e, "elapsedTimeLimit")) "censored" else "error"
          error_message <<- conditionMessage(e)
          error_class <<- paste(class(e), collapse = "/")
          NULL
        }
      )
      elapsed <- max(0, elapsed_seconds() - start)
      if (limited) setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
      list(
        elapsed_ms = 1000 * elapsed / inner_iterations,
        inner_iterations = inner_iterations,
        statistic = safe_statistic(value),
        status = status,
        error_message = error_message,
        error_class = error_class
      )
    },
    args = list(
      route = route, tree = tree, trait = trait,
      timeout_seconds = timeout_seconds,
      inner_iterations = inner_iterations
    )
  )
}

# Keep one isolated R session per installed package version.  Formal timings
# remain serialized and alternating; persistence removes thousands of process
# launches without allowing the two namespaces or DLLs to share a process.
sessions <- lapply(libraries, function(lib) {
  session <- callr::r_session$new(options = callr::r_session_options(
    user_profile = FALSE,
    system_profile = FALSE
  ))
  session$run(function(lib) {
    .libPaths(c(lib, .libPaths()))
    suppressPackageStartupMessages(
      library("fastphylosig", lib.loc = lib, character.only = TRUE)
    )
    invisible(TRUE)
  }, args = list(lib = lib))
  session
})
on.exit({
  for (session in sessions) {
    try(session$close(), silent = TRUE)
  }
}, add = TRUE)

timing_rows <- list()
timing_id <- 0L
checkpoint_file <- file.path(out_dir, "v2b_benchmark_timings_checkpoint.csv")
message("[stage2b2c-v2b] sequential paired benchmark starting")
for (shape in shape_names) {
  for (n in n_grid) {
    message("[stage2b2c-v2b] shape=", shape, " n=", n)
    fixture <- make_fixture(shape, n)
    trait <- make_trait(fixture)
    repeats <- repeats_for(n)
    for (route in route_names) {
      # Route-specific warmups are deliberately outside formal observations.
      warmups <- list()
      for (version in versions) {
        warmups[[version]] <- time_route(
          sessions[[version]], route, clone_fixture(fixture),
          clone_fixture(trait), timeout_seconds
        )
      }
      for (pair in seq_len(repeats)) {
        order_for_pair <- if (length(versions) == 2L && pair %% 2L == 0L) {
          rev(versions)
        } else {
          versions
        }
        for (position in seq_along(order_for_pair)) {
          version <- order_for_pair[[position]]
          measured <- if (!identical(warmups[[version]]$status, "ok")) {
            list(
              elapsed_ms = NA_real_, statistic = NA_real_,
              status = paste0("not_run_after_warmup_", warmups[[version]]$status),
              error_message = warmups[[version]]$error_message,
              error_class = warmups[[version]]$error_class,
              inner_iterations = warmups[[version]]$inner_iterations
            )
          } else tryCatch(
            time_route(
              sessions[[version]], route, clone_fixture(fixture),
              clone_fixture(trait), timeout_seconds
            ),
            error = function(e) list(
              elapsed_ms = NA_real_, statistic = NA_real_, status = "process_error",
              error_message = conditionMessage(e),
              error_class = paste(class(e), collapse = "/")
            )
          )
          timing_id <- timing_id + 1L
          timing_rows[[timing_id]] <- data.frame(
            version = version,
            shape = shape,
            n = as.integer(n),
            route = route,
            pair = as.integer(pair),
            position = as.integer(position),
            repeats_for_n = as.integer(repeats),
            elapsed_ms = as.numeric(measured$elapsed_ms),
            inner_iterations = if (!is.null(measured$inner_iterations)) {
              as.integer(measured$inner_iterations)
            } else NA_integer_,
            statistic = as.numeric(measured$statistic),
            status = as.character(measured$status),
            error_message = as.character(measured$error_message),
            error_class = as.character(measured$error_class),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    utils::write.csv(
      if (length(timing_rows)) do.call(rbind, timing_rows) else data.frame(),
      checkpoint_file, row.names = FALSE
    )
    rm(fixture, trait)
    gc(FALSE)
  }
}

timings <- if (length(timing_rows)) do.call(rbind, timing_rows) else
  data.frame()
timing_file <- if (single_mode) {
  file.path(out_dir, paste0("v2b_benchmark_timings_", single_label, ".csv"))
} else {
  file.path(out_dir, "v2b_benchmark_timings.csv")
}
utils::write.csv(timings, timing_file, row.names = FALSE)

summary_rows <- list()
summary_id <- 0L
if (nrow(timings)) {
  key <- interaction(
    timings$version, timings$shape, timings$n, timings$route,
    drop = TRUE, sep = "\r"
  )
  for (indices in split(seq_len(nrow(timings)), key)) {
    z <- timings[indices, , drop = FALSE]
    finite <- z$elapsed_ms[is.finite(z$elapsed_ms)]
    summary_id <- summary_id + 1L
    summary_rows[[summary_id]] <- data.frame(
      version = z$version[[1L]],
      shape = z$shape[[1L]],
      n = z$n[[1L]],
      route = z$route[[1L]],
      repeats_for_n = z$repeats_for_n[[1L]],
      observations = nrow(z),
      ok_observations = sum(z$status == "ok"),
      error_observations = sum(z$status == "error" | z$status == "process_error"),
      censored_observations = sum(z$status == "censored"),
      median_ms = if (length(finite)) stats::median(finite) else NA_real_,
      IQR_ms = if (length(finite)) stats::IQR(finite) else NA_real_,
      min_ms = if (length(finite)) min(finite) else NA_real_,
      max_ms = if (length(finite)) max(finite) else NA_real_,
      statuses = paste(sort(unique(z$status)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}
summary <- if (length(summary_rows)) do.call(rbind, summary_rows) else
  data.frame()
summary_file <- if (single_mode) {
  file.path(out_dir, paste0("v2b_benchmark_summary_", single_label, ".csv"))
} else {
  file.path(out_dir, "v2b_benchmark_summary.csv")
}
utils::write.csv(summary, summary_file, row.names = FALSE)

if (length(versions) == 2L && nrow(summary)) {
  paired_rows <- list()
  paired_id <- 0L
  paired_key <- interaction(summary$shape, summary$n, summary$route,
                            drop = TRUE, sep = "\r")
  for (indices in split(seq_len(nrow(summary)), paired_key)) {
    z <- summary[indices, , drop = FALSE]
    before <- z[z$version == versions[[1L]], , drop = FALSE]
    after <- z[z$version == versions[[2L]], , drop = FALSE]
    if (!nrow(before) || !nrow(after)) next
    paired_id <- paired_id + 1L
    before_median <- before$median_ms[[1L]]
    after_median <- after$median_ms[[1L]]
    paired_rows[[paired_id]] <- data.frame(
      shape = z$shape[[1L]],
      n = z$n[[1L]],
      route = z$route[[1L]],
      before_version = versions[[1L]],
      after_version = versions[[2L]],
      before_median_ms = before_median,
      before_IQR_ms = before$IQR_ms[[1L]],
      after_median_ms = after_median,
      after_IQR_ms = after$IQR_ms[[1L]],
      speedup_before_over_after = if (is.finite(before_median) &&
                                      is.finite(after_median) &&
                                      after_median > 0) {
        before_median / after_median
      } else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  utils::write.csv(
    if (length(paired_rows)) do.call(rbind, paired_rows) else data.frame(),
    file.path(out_dir, "v2b_paired_summary.csv"), row.names = FALSE
  )
}

byte_length <- function(value) {
  tryCatch({
    if (is.raw(value)) return(as.numeric(length(value)))
    as.numeric(length(charToRaw(enc2utf8(paste(as.character(value), collapse = "")))))
  }, error = function(e) NA_real_)
}

complexity_for <- function(lib, version, tree, commit) {
  callr::r(
    function(lib, version, tree, commit) {
      suppressPackageStartupMessages(
        library("fastphylosig", lib.loc = lib, character.only = TRUE)
      )
      ns <- asNamespace("fastphylosig")
      get_private <- function(name) {
        if (!exists(name, envir = ns, inherits = FALSE)) return(NULL)
        get(name, envir = ns, inherits = FALSE)
      }
      byte_length <- function(value) {
        tryCatch({
          if (is.raw(value)) return(as.numeric(length(value)))
          as.numeric(length(charToRaw(enc2utf8(
            paste(as.character(value), collapse = "")
          ))))
        }, error = function(e) NA_real_)
      }
      missing_row <- function(status, message = "", class = "") {
        data.frame(
          version = version, shape = "pectinate", n = length(tree$tip.label),
          edge_count = nrow(tree$edge),
          child_tuple_entries = NA_real_, descriptor_entries = NA_real_,
          fingerprint_bytes = NA_real_, snapshot_bytes = NA_real_,
          context_object_bytes = NA_real_, descendant_payload_count = NA_real_,
          context_schema_version = NA_real_, status = status,
          error_message = message, error_class = class,
          source_commit = commit, R_version = R.version.string,
          R_platform = R.version$platform, R_arch = R.version$arch,
          compiler = if (length(R.version$compiler)) {
            as.character(R.version$compiler[[1L]])
          } else "<unknown>",
          package_version = as.character(utils::packageVersion(
            "fastphylosig", lib.loc = lib
          )),
          stringsAsFactors = FALSE
        )
      }
      result <- tryCatch({
        canonicalize <- get_private(".safe_canonicalize_core")
        fingerprint <- get_private(".tree_fingerprint")
        if (!is.function(canonicalize) || !is.function(fingerprint)) {
          stop("canonicalization/fingerprint helper is unavailable")
        }
        canonical <- canonicalize(tree)
        canonical_argument <- ".canonical" %in% names(formals(fingerprint))
        fingerprint_value <- if (canonical_argument) {
          fingerprint(canonical, .canonical = TRUE)
        } else {
          fingerprint(canonical)
        }
        metrics <- attr(canonical, "fastphylosig_v2_metrics", exact = TRUE)
        contract <- attr(canonical, "fastphylosig_v2_contract", exact = TRUE)
        descriptor_entries <- if (is.list(metrics) &&
                                  !is.null(metrics$descriptor_entries)) {
          as.numeric(metrics$descriptor_entries)
        } else if (is.list(contract) && is.list(contract$descriptor_table)) {
          as.numeric(length(contract$descriptor_table))
        } else NA_real_
        child_tuple_entries <- if (is.list(contract) &&
                                   is.list(contract$descriptor_table)) {
          as.numeric(sum(vapply(
            contract$descriptor_table,
            function(x) length(x$child_descriptor_ids), integer(1L)
          )))
        } else NA_real_
        descendant_payload_count <- if (is.list(metrics) &&
                                        !is.null(metrics$descendant_label_vectors_materialized)) {
          as.numeric(metrics$descendant_label_vectors_materialized)
        } else NA_real_
        snapshot_fun <- get_private(".tree_protected_snapshot_v2")
        context <- fastphylosig::prepare_tree(tree)
        snapshot_value <- if (is.function(snapshot_fun)) {
          snapshot_fun(context$tree)
        } else NULL
        data.frame(
          version = version, shape = "pectinate", n = length(tree$tip.label),
          edge_count = nrow(tree$edge),
          child_tuple_entries = child_tuple_entries,
          descriptor_entries = descriptor_entries,
          fingerprint_bytes = byte_length(fingerprint_value),
          snapshot_bytes = byte_length(snapshot_value),
          context_object_bytes = as.numeric(object.size(context)),
          descendant_payload_count = descendant_payload_count,
          context_schema_version = if (!is.null(context$context_schema_version)) {
            as.numeric(context$context_schema_version)
          } else NA_real_,
          status = "ok", error_message = "", error_class = "",
          source_commit = commit, R_version = R.version.string,
          R_platform = R.version$platform, R_arch = R.version$arch,
          compiler = if (length(R.version$compiler)) {
            as.character(R.version$compiler[[1L]])
          } else "<unknown>",
          package_version = as.character(utils::packageVersion(
            "fastphylosig", lib.loc = lib
          )),
          stringsAsFactors = FALSE
        )
      }, error = function(e) missing_row(
        "error", conditionMessage(e), paste(class(e), collapse = "/")
      ))
      result
    },
    args = list(lib = lib, version = version, tree = tree, commit = commit),
    libpath = c(lib, .libPaths()),
    spinner = FALSE
  )
}

complexity_tree <- make_pectinate(complexity_n)
complexity_rows <- list()
for (version in versions) {
  complexity_rows[[length(complexity_rows) + 1L]] <- tryCatch(
    complexity_for(
      libraries[[version]], version, clone_fixture(complexity_tree),
      commit_for(version)
    ),
    error = function(e) data.frame(
      version = version, shape = "pectinate", n = complexity_n,
      edge_count = nrow(complexity_tree$edge),
      child_tuple_entries = NA_real_, descriptor_entries = NA_real_,
      fingerprint_bytes = NA_real_, snapshot_bytes = NA_real_,
      context_object_bytes = NA_real_, descendant_payload_count = NA_real_,
      context_schema_version = NA_real_, status = "process_error",
      error_message = conditionMessage(e),
      error_class = paste(class(e), collapse = "/"),
      source_commit = commit_for(version), R_version = R.version.string,
      R_platform = R.version$platform, R_arch = R.version$arch,
      compiler = if (length(R.version$compiler)) {
        as.character(R.version$compiler[[1L]])
      } else "<unknown>",
      package_version = package_version_at(libraries[[version]]),
      stringsAsFactors = FALSE
    )
  )
}
utils::write.csv(
  do.call(rbind, complexity_rows),
  file.path(out_dir, "v2b_pectinate_20000_complexity.csv"), row.names = FALSE
)

provenance_rows <- list()
for (version in versions) {
  provenance_rows[[length(provenance_rows) + 1L]] <- data.frame(
    version = version,
    library = libraries[[version]],
    source_commit = commit_for(version),
    package_version = package_version_at(libraries[[version]]),
    R_version = R.version.string,
    R_platform = R.version$platform,
    R_arch = R.version$arch,
    R_home = R.home(),
    compiler = if (length(R.version$compiler)) {
      as.character(R.version$compiler[[1L]])
    } else "<unknown>",
    OpenMP = "not inferred by this R harness",
    shapes = paste(shape_names, collapse = ","),
    n_grid = paste(n_grid, collapse = ","),
    complexity_n = as.integer(complexity_n),
    repeats = if (length(repeat_override) && !is.na(repeat_override)) {
      as.character(repeat_override)
    } else "n<10000:10; n>=10000:5; quick:2",
    warmup = "one route-specific warmup per version/fixture/route",
    alternation = if (length(versions) == 2L) {
      "odd pairs candidate1->v2; even pairs v2->candidate1"
    } else "single-library run; no paired alternation",
    timer_scope = "child timer excludes process startup, serialization, and library load",
    benchmark_parallelism = "none; callr children run sequentially",
    timeout_seconds = if (is.finite(timeout_seconds)) as.character(timeout_seconds) else "Inf",
    source_change_scope = "audit harness only; no production estimator change",
    stringsAsFactors = FALSE
  )
}
utils::write.csv(
  do.call(rbind, provenance_rows),
  file.path(out_dir, "v2b_provenance.csv"), row.names = FALSE
)
capture.output(sessionInfo(), file = file.path(out_dir, "v2b_sessionInfo.txt"))
message("[stage2b2c-v2b] sequential paired benchmark finished")
