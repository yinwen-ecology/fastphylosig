#!/usr/bin/env Rscript

# fastphylosig 0.2.0 Stage 2C post-V2 performance reprofile.
#
# Audit harness only.  This script writes evidence below a caller-supplied
# output directory.  It never edits package source, tests, documentation, or
# frozen V2 production bindings.  Formal runs are deliberately serialized in
# this one R process; private helper recorders are installed only for separate
# low-overhead diagnostic probes and are restored immediately.

raw_args <- commandArgs(trailingOnly = TRUE)
smoke <- any(raw_args %in% c("--smoke", "quick")) ||
  identical(Sys.getenv("FASTPHYLOSIG_STAGE2C_QUICK", unset = ""), "1")
raw_args <- raw_args[!raw_args %in% c("--smoke", "quick")]

repo_arg <- if (length(raw_args)) raw_args[[1L]] else getwd()
out_arg <- if (length(raw_args) >= 2L) raw_args[[2L]] else
  file.path(repo_arg, "benchmarks", "stage2c", "results")
repo <- normalizePath(repo_arg, winslash = "/", mustWork = FALSE)
out_dir <- normalizePath(out_arg, winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE)
Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")

env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

parse_ints <- function(value, name, minimum = 1L) {
  answer <- suppressWarnings(as.integer(
    trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1L]])
  ))
  if (!length(answer) || anyNA(answer) || any(answer < minimum)) {
    stop(name, " must contain comma-separated integers >= ", minimum, ".",
         call. = FALSE)
  }
  unique(answer)
}

parse_shapes <- function(value, name) {
  answer <- unique(trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1L]]))
  answer <- answer[nzchar(answer)]
  if (!length(answer) || any(!answer %in% c("balanced", "random", "pectinate"))) {
    stop(name, " must contain balanced, random, and/or pectinate.",
         call. = FALSE)
  }
  answer
}

n_grid <- function(env_name, formal, smoke_value, minimum = 2L) {
  value <- Sys.getenv(env_name, unset = "")
  if (nzchar(value)) parse_ints(value, env_name, minimum = minimum) else
    if (smoke) smoke_value else formal
}

repeats_for <- local({
  override <- suppressWarnings(as.integer(Sys.getenv(
    "FASTPHYLOSIG_STAGE2C_REPEATS", unset = ""
  )))
  function(n) {
    if (length(override) == 1L && is.finite(override) && override > 0L) {
      return(override)
    }
    if (smoke) return(1L)
    if (n >= 10000L) 5L else 10L
  }
})

timeout_seconds <- suppressWarnings(as.numeric(Sys.getenv(
  "FASTPHYLOSIG_STAGE2C_TIMEOUT", unset = ""
)))
if (length(timeout_seconds) != 1L || is.na(timeout_seconds) ||
    !is.finite(timeout_seconds) || timeout_seconds <= 0) {
  timeout_seconds <- if (smoke) Inf else 60
}

positive_bytes_option <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = "")))
  if (length(value) != 1L || is.na(value) || !is.finite(value) || value <= 0) {
    return(as.numeric(default))
  }
  value
}

# phytools::phylosig(K) constructs a dense n x n covariance matrix.  These
# deterministic preflight limits prevent a formal reference warmup from
# allocating an unbounded dense object.  A skipped reference is recorded as
# resource_limit and never converted into a speedup estimate.
reference_dense_cap_bytes <- positive_bytes_option(
  "FASTPHYLOSIG_STAGE2C_REFERENCE_DENSE_CAP_BYTES", 512 * 1024^2
)
reference_batch_work_cap_bytes <- positive_bytes_option(
  "FASTPHYLOSIG_STAGE2C_REFERENCE_BATCH_WORK_CAP_BYTES", 2 * 1024^3
)
reference_single_cubic_cap <- positive_bytes_option(
  "FASTPHYLOSIG_STAGE2C_REFERENCE_SINGLE_CUBIC_CAP", 1e8
)
reference_batch_cubic_cap <- positive_bytes_option(
  "FASTPHYLOSIG_STAGE2C_REFERENCE_BATCH_CUBIC_CAP", 1e8
)

mcmc_sim <- suppressWarnings(as.integer(Sys.getenv(
  "FASTPHYLOSIG_STAGE2C_MCMC_SIM", unset = "1000"
)))
thin <- suppressWarnings(as.integer(Sys.getenv(
  "FASTPHYLOSIG_STAGE2C_THIN", unset = "10"
)))
burn <- suppressWarnings(as.integer(Sys.getenv(
  "FASTPHYLOSIG_STAGE2C_BURN", unset = "100"
)))
if (anyNA(c(mcmc_sim, thin, burn)) || mcmc_sim < 1L || thin < 1L ||
    burn < 1L || burn > mcmc_sim) {
  stop("MCMC controls must satisfy 1 <= burn <= mcmc_sim and thin >= 1.",
       call. = FALSE)
}

parse_ncores <- function() {
  value <- Sys.getenv("FASTPHYLOSIG_STAGE2C_NCORES", unset = "")
  if (nzchar(value)) return(parse_ints(value, "FASTPHYLOSIG_STAGE2C_NCORES"))
  detected <- tryCatch(parallel::detectCores(logical = TRUE),
                       error = function(e) 1L)
  if (smoke || !is.finite(detected) || detected < 2L) 1L else c(1L, 2L)
}
ncores_grid <- parse_ncores()

single_n <- n_grid(
  "FASTPHYLOSIG_STAGE2C_K_SINGLE_N",
  c(50L, 100L, 300L, 500L, 1000L, 2000L, 5000L, 10000L, 20000L),
  c(50L, 100L)
)
batch_n <- n_grid("FASTPHYLOSIG_STAGE2C_K_BATCH_N", c(500L, 2000L, 5000L),
                  c(500L))
batch_traits <- n_grid("FASTPHYLOSIG_STAGE2C_K_BATCH_TRAITS",
                       c(1L, 8L, 32L, 100L), c(1L, 8L), minimum = 1L)
perm_n <- n_grid("FASTPHYLOSIG_STAGE2C_K_PERM_N",
                 c(500L, 2000L, 5000L, 10000L), c(500L))
perm_nsim <- n_grid("FASTPHYLOSIG_STAGE2C_K_PERM_NSIM",
                    c(199L, 999L, 9999L), c(19L))
lambda_n <- n_grid("FASTPHYLOSIG_STAGE2C_LAMBDA_N",
                   c(100L, 500L, 2000L, 5000L, 10000L), c(100L))
d_n <- n_grid("FASTPHYLOSIG_STAGE2C_D_N", c(100L, 500L, 2000L), c(100L))
d_nsim <- n_grid("FASTPHYLOSIG_STAGE2C_D_NSIM",
                  c(199L, 999L, 9999L), c(19L))
delta_n <- n_grid("FASTPHYLOSIG_STAGE2C_DELTA_N", c(50L, 100L, 500L), c(50L))
delta_nsim <- n_grid("FASTPHYLOSIG_STAGE2C_DELTA_NSIM", c(19L, 199L), c(19L))
memory_n <- n_grid("FASTPHYLOSIG_STAGE2C_MEMORY_N",
                   c(500L, 1000L, 5000L, 10000L, 20000L), c(500L))

single_shapes <- if (nzchar(Sys.getenv("FASTPHYLOSIG_STAGE2C_SINGLE_SHAPES", ""))) {
  parse_shapes(Sys.getenv("FASTPHYLOSIG_STAGE2C_SINGLE_SHAPES"),
               "FASTPHYLOSIG_STAGE2C_SINGLE_SHAPES")
} else "random"
batch_shapes <- if (nzchar(Sys.getenv("FASTPHYLOSIG_STAGE2C_BATCH_SHAPES", ""))) {
  parse_shapes(Sys.getenv("FASTPHYLOSIG_STAGE2C_BATCH_SHAPES"),
               "FASTPHYLOSIG_STAGE2C_BATCH_SHAPES")
} else "random"
perm_shapes <- if (nzchar(Sys.getenv("FASTPHYLOSIG_STAGE2C_PERM_SHAPES", ""))) {
  parse_shapes(Sys.getenv("FASTPHYLOSIG_STAGE2C_PERM_SHAPES"),
               "FASTPHYLOSIG_STAGE2C_PERM_SHAPES")
} else "random"
lambda_shapes <- if (nzchar(Sys.getenv("FASTPHYLOSIG_STAGE2C_LAMBDA_SHAPES", ""))) {
  parse_shapes(Sys.getenv("FASTPHYLOSIG_STAGE2C_LAMBDA_SHAPES"),
               "FASTPHYLOSIG_STAGE2C_LAMBDA_SHAPES")
} else if (smoke) "balanced" else c("balanced", "random", "pectinate")
d_shapes <- if (nzchar(Sys.getenv("FASTPHYLOSIG_STAGE2C_D_SHAPES", ""))) {
  parse_shapes(Sys.getenv("FASTPHYLOSIG_STAGE2C_D_SHAPES"),
               "FASTPHYLOSIG_STAGE2C_D_SHAPES")
} else "random"
delta_shapes <- if (nzchar(Sys.getenv("FASTPHYLOSIG_STAGE2C_DELTA_SHAPES", ""))) {
  parse_shapes(Sys.getenv("FASTPHYLOSIG_STAGE2C_DELTA_SHAPES"),
               "FASTPHYLOSIG_STAGE2C_DELTA_SHAPES")
} else "random"

if (!file.exists(file.path(repo, "DESCRIPTION")) &&
    !nzchar(Sys.getenv("FASTPHYLOSIG_STAGE2C_LIBRARY", unset = ""))) {
  stop("repo must contain DESCRIPTION or FASTPHYLOSIG_STAGE2C_LIBRARY must be set.",
       call. = FALSE)
}

library_dir <- Sys.getenv("FASTPHYLOSIG_STAGE2C_LIBRARY", unset = "")
if (nzchar(library_dir)) {
  library_dir <- normalizePath(library_dir, winslash = "/", mustWork = TRUE)
  .libPaths(c(library_dir, .libPaths()))
  suppressPackageStartupMessages(library("fastphylosig", character.only = TRUE,
                                         lib.loc = library_dir))
} else {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required when loading fastphylosig from source.", call. = FALSE)
  }
  pkgload::load_all(repo, quiet = TRUE)
}
suppressPackageStartupMessages(requireNamespace("ape"))
has_phytools <- requireNamespace("phytools", quietly = TRUE)
has_callr <- requireNamespace("callr", quietly = TRUE)

ns <- asNamespace("fastphylosig")

clock <- function() unname(as.numeric(proc.time()[["elapsed"]]))

clone_object <- function(x) unserialize(serialize(x, NULL, version = 3L))

safe_scalar <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  if (is.data.frame(x)) return(NA_real_)
  y <- suppressWarnings(as.numeric(x[[1L]]))
  if (length(y)) y[[1L]] else NA_real_
}

estimate_of <- function(value, method = NULL) {
  if (is.data.frame(value)) {
    candidates <- switch(method,
      K = c("K_fast", "K"), lambda = c("lambda_fast", "lambda"),
      D = c("D_fast", "DEstimate"), Delta = c("Delta_fast", "delta"),
      c("K_fast", "lambda_fast", "D_fast", "Delta_fast", "DEstimate")
    )
    for (nm in candidates) if (nm %in% names(value)) return(safe_scalar(value[[nm]]))
  }
  if (is.list(value)) {
    candidates <- switch(method,
      K = c("K", "K_fast"), lambda = c("lambda", "lambda_fast"),
      D = c("DEstimate", "D_fast"), Delta = c("Delta", "Delta_fast"),
      c("K", "lambda", "DEstimate", "Delta", "K_fast", "lambda_fast")
    )
    for (nm in candidates) if (!is.null(value[[nm]])) return(safe_scalar(value[[nm]]))
  }
  safe_scalar(value)
}

warnings_text <- function(warnings) {
  if (!length(warnings)) "" else paste(unique(warnings), collapse = " | ")
}

timed_call <- function(fun, timeout = timeout_seconds) {
  warnings <- character()
  status <- "ok"
  error_message <- ""
  error_class <- ""
  value <- NULL
  started <- clock()
  limited <- is.finite(timeout) && timeout > 0
  if (limited) {
    setTimeLimit(cpu = Inf, elapsed = timeout, transient = TRUE)
    # Always clear the transient limit, including non-local exits and errors.
    on.exit(
      setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE),
      add = TRUE
    )
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
        "resource_limit"
      } else {
        "error"
      }
      error_message <<- conditionMessage(e)
      error_class <<- paste(class(e), collapse = "/")
      NULL
    }
  )
  elapsed <- max(0, clock() - started)
  list(value = value, elapsed_s = elapsed, status = status,
       warnings = warnings, warning_text = warnings_text(warnings),
       error_message = error_message, error_class = error_class)
}

reference_precheck <- function(scope, n, traits = 1L) {
  n <- as.double(n); traits <- as.double(traits)
  dense_bytes <- n * n * 8
  workload_bytes <- dense_bytes * max(1, traits)
  cubic_work <- n * n * n * max(1, traits)
  reason <- character()
  if (identical(scope, "single_reference") && n >= 10000) {
    reason <- sprintf(
      "single phytools reference is not run at n >= 10000; estimated dense covariance %.0f bytes",
      dense_bytes
    )
  } else if (dense_bytes > reference_dense_cap_bytes) {
    reason <- sprintf(
      "estimated dense covariance %.0f bytes exceeds cap %.0f bytes",
      dense_bytes, reference_dense_cap_bytes
    )
  } else if (identical(scope, "batch_reference") &&
             workload_bytes > reference_batch_work_cap_bytes) {
    reason <- sprintf(
      "estimated reference dense/workload bytes %.0f (n=%d, traits=%d) exceeds cap %.0f bytes",
      workload_bytes, as.integer(n), as.integer(traits), reference_batch_work_cap_bytes
    )
  } else if (identical(scope, "single_reference") &&
             cubic_work > reference_single_cubic_cap) {
    reason <- sprintf(
      "estimated dense cubic work %.0f exceeds single-reference cap %.0f",
      cubic_work, reference_single_cubic_cap
    )
  } else if (identical(scope, "batch_reference") &&
             cubic_work > reference_batch_cubic_cap) {
    reason <- sprintf(
      paste0(
        "estimated repeated dense cubic work %.0f (n=%d, traits=%d) ",
        "exceeds batch-reference cap %.0f"
      ),
      cubic_work, as.integer(n), as.integer(traits),
      reference_batch_cubic_cap
    )
  }
  list(
    limited = length(reason) > 0L,
    estimated_dense_bytes = dense_bytes,
    estimated_workload_bytes = workload_bytes,
    estimated_cubic_work = cubic_work,
    reason = if (length(reason)) reason else ""
  )
}

resource_limited_timing <- function(resource) {
  list(
    value = NULL, elapsed_s = NA_real_, status = "resource_limit",
    warnings = character(), warning_text = "",
    error_message = resource$reason, error_class = "resource_limit"
  )
}

empty_timing_row <- function(...) {
  data.frame(..., elapsed_s = NA_real_, status = "not_run",
             error_message = "", error_class = "", warning_text = "",
             statistic = NA_real_, stringsAsFactors = FALSE)
}

summarize_timings <- function(dat, keys, value = "elapsed_s") {
  if (!nrow(dat)) return(data.frame())
  key <- do.call(paste, c(dat[keys], sep = "\r"))
  groups <- split(seq_len(nrow(dat)), key, drop = TRUE)
  out <- lapply(groups, function(ii) {
    z <- dat[ii, , drop = FALSE]
    finite <- z[[value]][is.finite(z[[value]])]
    first <- z[1L, keys, drop = FALSE]
    cbind(
      first,
      data.frame(
        observations = nrow(z), ok_observations = sum(z$status == "ok"),
        error_observations = sum(z$status == "error"),
        resource_limit_observations = sum(z$status == "resource_limit"),
        warning_observations = sum(nzchar(z$warning_text)),
        median_s = if (length(finite)) stats::median(finite) else NA_real_,
        IQR_s = if (length(finite)) stats::IQR(finite) else NA_real_,
        min_s = if (length(finite)) min(finite) else NA_real_,
        max_s = if (length(finite)) max(finite) else NA_real_,
        statuses = paste(sort(unique(z$status)), collapse = ";"),
        warnings = paste(unique(z$warning_text[nzchar(z$warning_text)]),
                         collapse = " | "),
        stringsAsFactors = FALSE
      )
    )
  })
  do.call(rbind, out)
}

write_csv <- function(dat, name) {
  utils::write.csv(dat, file.path(out_dir, name), row.names = FALSE,
                   na = "NA")
  invisible(dat)
}

append_csv <- function(dat, name) {
  if (!nrow(dat)) return(invisible(NULL))
  path <- file.path(out_dir, name)
  if (file.exists(path)) {
    old <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE),
                    error = function(e) NULL)
    if (!is.null(old) && ncol(old) == ncol(dat) &&
        identical(names(old), names(dat))) dat <- rbind(old, dat)
  }
  write_csv(dat, name)
}

make_balanced <- function(n) {
  n <- as.integer(n)
  root <- n + 1L
  n_internal <- n - 1L
  edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  queue_node <- queue_lo <- queue_hi <- integer(n_internal)
  queue_node[[1L]] <- root; queue_lo[[1L]] <- 1L; queue_hi[[1L]] <- n
  head <- 1L; tail <- 1L; edge_i <- 0L; next_internal <- root + 1L
  while (head <= tail) {
    node <- queue_node[[head]]; lo <- queue_lo[[head]]; hi <- queue_hi[[head]]
    head <- head + 1L; mid <- floor((lo + hi) / 2)
    child_lo <- c(lo, mid + 1L); child_hi <- c(mid, hi)
    for (j in 1:2) {
      child <- if (child_lo[[j]] == child_hi[[j]]) child_lo[[j]] else {
        next_node <- next_internal; next_internal <- next_internal + 1L
        tail <- tail + 1L; queue_node[[tail]] <- next_node
        queue_lo[[tail]] <- child_lo[[j]]; queue_hi[[tail]] <- child_hi[[j]]
        next_node
      }
      edge_i <- edge_i + 1L; edge[edge_i, ] <- c(node, child)
    }
  }
  edge <- edge[seq_len(edge_i), , drop = FALSE]
  tree <- list(edge = edge, tip.label = paste0("sp", seq_len(n)),
               edge.length = 0.5 + seq_len(nrow(edge)) / nrow(edge),
               Nnode = n_internal)
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

make_pectinate <- function(n) {
  n <- as.integer(n); root <- n + 1L
  if (n == 2L) {
    edge <- matrix(c(root, 1L, root, 2L), ncol = 2L, byrow = TRUE)
  } else {
    edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
    current <- root; next_internal <- root + 1L; edge_i <- 0L
    for (tip in seq_len(n - 2L)) {
      edge_i <- edge_i + 1L; edge[edge_i, ] <- c(current, tip)
      edge_i <- edge_i + 1L; edge[edge_i, ] <- c(current, next_internal)
      current <- next_internal; next_internal <- next_internal + 1L
    }
    edge_i <- edge_i + 1L; edge[edge_i, ] <- c(current, n - 1L)
    edge_i <- edge_i + 1L; edge[edge_i, ] <- c(current, n)
    edge <- edge[seq_len(edge_i), , drop = FALSE]
  }
  tree <- list(edge = edge, tip.label = paste0("sp", seq_len(n)),
               edge.length = 0.5 + seq_len(nrow(edge)) / nrow(edge),
               Nnode = n - 1L)
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

make_random <- function(n) {
  set.seed(20260908L + as.integer(n))
  tree <- ape::reorder.phylo(ape::rtree(n), order = "postorder")
  tree$tip.label <- paste0("sp", seq_len(n))
  tree$edge.length <- 0.5 + seq_len(nrow(tree$edge)) / nrow(tree$edge)
  tree
}

make_fixture <- function(shape, n) switch(
  shape, balanced = make_balanced(n), random = make_random(n),
  pectinate = make_pectinate(n), stop("unknown tree shape: ", shape, call. = FALSE)
)

make_continuous_trait <- function(tree, offset = 0L) {
  n <- length(tree$tip.label); i <- seq_len(n) + offset
  stats::setNames(sin(i) + cos(i / 7) + i / (n + 11), tree$tip.label)
}

make_continuous_matrix <- function(tree, p) {
  n <- length(tree$tip.label); i <- seq_len(n); j <- seq_len(p)
  values <- outer(i, j, function(a, b) sin(a + b / 11) +
                    cos(a / 7 + b) + (a + b) / (n + p + 17))
  dimnames(values) <- list(tree$tip.label, paste0("trait", j))
  values
}

make_binary_trait <- function(tree, seed = 1L) {
  set.seed(seed); n <- length(tree$tip.label)
  value <- sample(c(0, 1), n, replace = TRUE)
  if (length(unique(value)) < 2L) value[[1L]] <- 1 - value[[1L]]
  stats::setNames(value, tree$tip.label)
}

make_categorical_trait <- function(tree, seed = 1L) {
  set.seed(seed); n <- length(tree$tip.label)
  value <- sample(c("a", "b", "c"), n, replace = TRUE)
  if (length(unique(value)) < 2L) value[[1L]] <- if (value[[1L]] == "a") "b" else "a"
  stats::setNames(value, tree$tip.label)
}

private_get <- function(name, where = ns) {
  if (!exists(name, envir = where, inherits = FALSE)) return(NULL)
  get(name, envir = where, inherits = FALSE)
}

# Install low-overhead cumulative recorders in an in-memory namespace only.
# The wrapper closure does no allocations apart from the two scalar updates;
# all authoritative wall-times are collected without these recorders.
record_bindings <- function(specs) {
  if (!length(specs)) return(list(state = new.env(parent = emptyenv()), old = list()))
  state <- new.env(parent = emptyenv())
  state$elapsed <- setNames(rep(0, length(specs)), vapply(specs, `[[`, "", "label"))
  state$calls <- setNames(integer(length(specs)), vapply(specs, `[[`, "", "label"))
  old <- list()
  installed <- character()
  for (spec in specs) {
    name <- spec$name; label <- spec$label; where <- spec$where
    if (!exists(name, envir = where, inherits = FALSE)) next
    old[[paste(environmentName(where), name, sep = "::")]] <-
      list(name = name, where = where, value = get(name, envir = where, inherits = FALSE),
           locked = bindingIsLocked(name, where))
    ok <- tryCatch({
      old_value <- old[[length(old)]]$value
      wrapped <- local({
        target <- old_value; target_label <- label; target_state <- state
        function(...) {
          started <- clock()
          on.exit({
            target_state$elapsed[[target_label]] <-
              target_state$elapsed[[target_label]] + max(0, clock() - started)
            target_state$calls[[target_label]] <-
              target_state$calls[[target_label]] + 1L
          }, add = TRUE)
          target(...)
        }
      })
      if (bindingIsLocked(name, where)) unlockBinding(name, where)
      assign(name, wrapped, envir = where)
      if (isTRUE(old[[length(old)]]$locked)) lockBinding(name, where)
      installed <<- c(installed, paste(environmentName(where), name, sep = "::"))
      TRUE
    }, error = function(e) FALSE)
    if (!ok) {
      item <- old[[paste(environmentName(where), name, sep = "::")]]
      if (!is.null(item)) {
        try({
          if (bindingIsLocked(name, where)) unlockBinding(name, where)
          assign(name, item$value, envir = where)
          if (isTRUE(item$locked)) lockBinding(name, where)
        }, silent = TRUE)
        old[[paste(environmentName(where), name, sep = "::")]] <- NULL
      }
    }
  }
  list(state = state, old = old, installed = installed)
}

restore_bindings <- function(recorder) {
  if (is.null(recorder$old) || !length(recorder$old)) return(invisible(NULL))
  for (item in rev(recorder$old)) {
    name <- item$name; where <- item$where
    try({
      if (bindingIsLocked(name, where)) unlockBinding(name, where)
      assign(name, item$value, envir = where)
      if (isTRUE(item$locked)) lockBinding(name, where)
    }, silent = TRUE)
  }
  invisible(NULL)
}

parallel_record_specs <- function() {
  if (!exists("parallel", mode = "namespace")) return(list())
  pns <- asNamespace("parallel")
  lapply(c("makeCluster", "clusterEvalQ", "clusterSetRNGStream",
           "parLapply", "stopCluster"), function(name) {
    list(name = name, label = paste0("parallel_", name), where = pns)
  })
}

phase_probe <- function(fun, labels, timeout = timeout_seconds,
                        include_parallel = FALSE) {
  specs <- lapply(labels, function(name) {
    list(name = name, label = name, where = ns)
  })
  parallel_labels <- character()
  if (isTRUE(include_parallel)) {
    ps <- parallel_record_specs()
    specs <- c(specs, ps)
    parallel_labels <- vapply(ps, `[[`, "", "label")
  }
  recorder <- record_bindings(specs)
  on.exit(restore_bindings(recorder), add = TRUE)
  result <- timed_call(fun, timeout = timeout)
  phase <- data.frame(
    label = names(recorder$state$elapsed),
    elapsed_s = as.numeric(recorder$state$elapsed),
    calls = as.integer(recorder$state$calls),
    stringsAsFactors = FALSE
  )
  result$phase <- phase
  result$instrumented_bindings <- length(recorder$installed)
  result$parallel_labels <- parallel_labels
  result
}

call_fast_k <- function(route, tree, context, trait, test = FALSE,
                        nsim = 199L, permutations = NULL, ncores = 1L) {
  target <- if (identical(route, "prepared")) context else tree
  fastphylosig::fast_k(
    target, trait, test = test, nsim = nsim, permutations = permutations,
    return_sim = FALSE, keep_null = FALSE, verbose = FALSE, progress = FALSE,
    ncores = ncores
  )
}

call_fast_signal_k <- function(tree, trait) {
  fastphylosig::fast_signal(
    tree, x = trait, method = "K", test = FALSE, verbose = FALSE,
    progress = FALSE, ncores = 1L
  )
}

call_fast_lambda <- function(route, tree, context, trait, test = FALSE) {
  target <- if (identical(route, "prepared")) context else tree
  fastphylosig::fast_lambda(
    target, trait, test = test, verbose = FALSE, progress = FALSE,
    ncores = 1L, lambda_profile = FALSE
  )
}

call_fast_d <- function(route, tree, context, trait, test = FALSE,
                        nsim = 199L, ncores = 1L,
                        random_states = NULL, brownian_states = NULL) {
  target <- if (identical(route, "prepared")) context else tree
  fastphylosig::fast_d(
    target, trait, test = test, nsim = nsim, return_sim = FALSE,
    keep_null = FALSE, verbose = FALSE, progress = FALSE, ncores = ncores,
    random_states = random_states, brownian_states = brownian_states
  )
}

call_fast_delta <- function(route, tree, context, trait, test = FALSE,
                            nsim = 19L, ncores = 1L,
                            permutations = NULL) {
  target <- if (identical(route, "prepared")) context else tree
  fastphylosig::fast_delta(
    target, trait, test = test, nsim = nsim, mcmc_sim = mcmc_sim,
    thin = thin, burn = burn, permutations = permutations,
    return_sim = FALSE, verbose = FALSE, progress = FALSE, ncores = ncores
  )
}

call_reference_k <- function(tree, trait) {
  if (!has_phytools) stop("phytools unavailable", call. = FALSE)
  phytools::phylosig(tree, trait, method = "K", test = FALSE)
}

call_reference_lambda <- function(tree, trait) {
  if (!has_phytools) stop("phytools unavailable", call. = FALSE)
  phytools::phylosig(tree, trait, method = "lambda", test = FALSE)
}

result_status <- function(value) {
  if (is.data.frame(value) && "status" %in% names(value)) {
    as.character(value$status[[1L]])
  } else if (is.list(value) && !is.null(value$status)) {
    as.character(value$status[[1L]])
  } else "ok"
}

result_p <- function(value, method) {
  field <- switch(method, K = "P", lambda = "P", D = "P_random", Delta = "P_fast")
  if (is.data.frame(value) && field %in% names(value)) return(safe_scalar(value[[field]]))
  if (is.list(value) && !is.null(value[[field]])) return(safe_scalar(value[[field]]))
  NA_real_
}

result_components <- function(value, method) {
  names <- switch(method,
    K = c("K", "P", "MCSE_P", "status"),
    lambda = c("lambda", "logL", "LR", "P", "status"),
    D = c("DEstimate", "Pval1", "Pval0", "status"),
    Delta = c("Delta", "P", "P_fast", "alpha_mean", "beta_mean", "status"),
    character()
  )
  answer <- setNames(as.list(rep(NA_real_, length(names))), names)
  for (nm in names) {
    if (identical(nm, "status")) {
      answer[[nm]] <- result_status(value)
    } else if (is.data.frame(value) && nm %in% names(value)) {
      answer[[nm]] <- safe_scalar(value[[nm]])
    } else if (is.list(value) && !is.null(value[[nm]])) {
      answer[[nm]] <- safe_scalar(value[[nm]])
    }
  }
  answer
}

timing_row <- function(method, workload, shape, n, route, repeat_id,
                       position, nsim = NA_integer_, ncores = 1L,
                       timed, value = NULL, extra = list(), resource = NULL) {
  components <- result_components(value, method)
  if (is.null(resource)) {
    resource <- list(limited = FALSE, estimated_dense_bytes = NA_real_,
                     estimated_workload_bytes = NA_real_,
                     estimated_cubic_work = NA_real_, reason = "")
  }
  row <- data.frame(
    method = method, workload = workload, shape = shape, n = as.integer(n),
    route = route, repeat_id = as.integer(repeat_id), position = as.integer(position),
    nsim = as.integer(nsim), ncores = as.integer(ncores),
    elapsed_s = as.numeric(timed$elapsed_s), status = timed$status,
    error_message = timed$error_message, error_class = timed$error_class,
    warning_text = timed$warning_text, statistic = estimate_of(value, method),
    resource_limit = isTRUE(resource$limited),
    estimated_dense_bytes = as.numeric(resource$estimated_dense_bytes),
    estimated_workload_bytes = as.numeric(resource$estimated_workload_bytes),
    estimated_cubic_work = as.numeric(resource$estimated_cubic_work),
    resource_reason = as.character(resource$reason),
    P = as.numeric(components$P %||% components$P_fast %||% components$Pval1),
    stringsAsFactors = FALSE
  )
  if (length(extra)) row <- cbind(row, as.data.frame(extra, stringsAsFactors = FALSE))
  row
}

`%||%` <- function(x, y) if (is.null(x)) y else x

phase_rows <- function(probe, method, workload, shape, n, route,
                       nsim = NA_integer_, ncores = 1L, mode = "diagnostic",
                       total_s = probe$elapsed_s, note = "") {
  if (is.null(probe$phase) || !nrow(probe$phase)) return(data.frame())
  z <- probe$phase
  known <- z$elapsed_s
  share <- if (is.finite(total_s) && total_s > 0) known / total_s else NA_real_
  data.frame(
    method = method, workload = workload, shape = shape, n = as.integer(n),
    route = route, nsim = as.integer(nsim), ncores = as.integer(ncores),
    mode = mode, phase = z$label, elapsed_s = z$elapsed_s,
    share_total = share, calls = z$calls,
    measurement_scope = "diagnostic; inclusive helpers may overlap",
    note = note, stringsAsFactors = FALSE
  )
}

run_k_single <- function() {
  rows <- list(); at <- 0L
  for (shape in single_shapes) for (n in single_n) {
    message("[stage2c] K single shape=", shape, " n=", n)
    tree <- make_fixture(shape, n); trait <- make_continuous_trait(tree)
    context_timed <- timed_call(function() fastphylosig::prepare_tree(tree))
    context <- context_timed$value
    if (is.null(context)) next
    routes <- c("A_reference", "B_fast_k_raw", "C_fast_k_prepared",
                "D_fast_signal_raw")
    reference_resource <- reference_precheck("single_reference", n, traits = 1L)
    route_fun <- list(
      A_reference = function() call_reference_k(clone_object(tree), trait),
      B_fast_k_raw = function() call_fast_k("raw", clone_object(tree), context, trait),
      C_fast_k_prepared = function() call_fast_k("prepared", tree, context, trait),
      D_fast_signal_raw = function() call_fast_signal_k(clone_object(tree), trait)
    )
    warmup <- lapply(routes, function(route) {
      if (identical(route, "A_reference") && !has_phytools) {
        return(list(status = "unavailable"))
      }
      if (identical(route, "A_reference") && reference_resource$limited) {
        return(resource_limited_timing(reference_resource))
      }
      set.seed(810000L + n); timed_call(route_fun[[route]], timeout = timeout_seconds)
    })
    reps <- repeats_for(n)
    for (rep in seq_len(reps)) {
      order <- if (rep %% 2L == 0L) rev(routes) else routes
      for (pos in seq_along(order)) {
        route <- order[[pos]]
        timed <- if (identical(route, "A_reference") && !has_phytools) {
          list(value = NULL, elapsed_s = NA_real_, status = "unavailable",
               error_message = "phytools unavailable", error_class = "",
               warning_text = "", warnings = character())
        } else if (identical(route, "A_reference") && reference_resource$limited) {
          resource_limited_timing(reference_resource)
        } else {
          set.seed(810000L + n + rep)
          timed_call(route_fun[[route]], timeout = timeout_seconds)
        }
        at <- at + 1L
        rows[[at]] <- timing_row(
          "K", "single_trait_crossover", shape, n, route, rep, pos,
          nsim = NA_integer_, ncores = 1L, timed = timed, value = timed$value,
          extra = list(warmup_status = if (is.null(warmup[[route]])) "unavailable" else
                         warmup[[route]]$status),
          resource = if (identical(route, "A_reference")) reference_resource else NULL
        )
      }
      prep_timed <- timed_call(function() fastphylosig::prepare_tree(tree),
                                timeout = timeout_seconds)
      at <- at + 1L
      rows[[at]] <- timing_row(
        "K", "single_trait_crossover", shape, n, "prepare_tree", rep, 1L,
        timed = prep_timed, value = prep_timed$value,
        extra = list(warmup_status = "n/a")
      )
    }
    rm(tree, trait, context); gc(FALSE)
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  write_csv(out, "stage2c_k_single_timings.csv")
  summary <- summarize_timings(
    out, c("method", "workload", "shape", "n", "route", "nsim", "ncores")
  )
  if (nrow(summary)) {
    prepare <- summary[summary$route == "prepare_tree", , drop = FALSE]
    compute <- summary[summary$route == "C_fast_k_prepared", , drop = FALSE]
    if (nrow(prepare) && nrow(compute)) {
      end_to_end <- merge(prepare[, c("shape", "n", "median_s", "IQR_s")],
                          compute[, c("shape", "n", "median_s", "IQR_s")],
                          by = c("shape", "n"), suffixes = c("_prepare", "_compute"))
      end_to_end$method <- "K"; end_to_end$workload <- "single_trait_crossover"
      end_to_end$route <- "C_prepared_end_to_end_derived"
      end_to_end$prepared_end_to_end_s <- end_to_end$median_s_prepare + end_to_end$median_s_compute
      write_csv(summary, "stage2c_k_single_summary.csv")
      write_csv(end_to_end, "stage2c_k_single_prepared_end_to_end.csv")
    } else write_csv(summary, "stage2c_k_single_summary.csv")
  } else write_csv(summary, "stage2c_k_single_summary.csv")
  out
}

run_k_batch <- function() {
  rows <- list(); at <- 0L
  for (shape in batch_shapes) for (n in batch_n) {
    tree <- make_fixture(shape, n); context <- timed_call(function() fastphylosig::prepare_tree(tree))$value
    if (is.null(context)) next
    X <- make_continuous_matrix(tree, max(batch_traits))
    for (p in batch_traits) {
      Xp <- X[, seq_len(p), drop = FALSE]
      reference_resource <- reference_precheck("batch_reference", n, traits = p)
      loops <- list(
        raw_matrix = function() call_fast_k("raw", tree, context, Xp),
        prepared_matrix = function() call_fast_k("prepared", tree, context, Xp),
        raw_column_loop = function() for (j in seq_len(p)) call_fast_k("raw", tree, context, Xp[, j]),
        prepared_column_loop = function() for (j in seq_len(p)) call_fast_k("prepared", tree, context, Xp[, j]),
        reference_column_loop = function() {
          if (!has_phytools) stop("phytools unavailable", call. = FALSE)
          for (j in seq_len(p)) phytools::phylosig(tree, Xp[, j], method = "K", test = FALSE)
        }
      )
      routes <- names(loops); reps <- repeats_for(n)
      for (route in routes) {
        if (identical(route, "reference_column_loop") && !has_phytools) next
        if (identical(route, "reference_column_loop") && reference_resource$limited) next
        invisible(timed_call(loops[[route]], timeout = timeout_seconds))
      }
      for (rep in seq_len(reps)) {
        order <- if (rep %% 2L == 0L) rev(routes) else routes
        for (pos in seq_along(order)) {
          route <- order[[pos]]
          if (identical(route, "reference_column_loop") && !has_phytools) {
            timed <- list(value = NULL, elapsed_s = NA_real_, status = "unavailable",
                          error_message = "phytools unavailable", error_class = "",
                          warning_text = "", warnings = character())
          } else if (identical(route, "reference_column_loop") && reference_resource$limited) {
            timed <- resource_limited_timing(reference_resource)
          } else {
            set.seed(820000L + n + p + rep); timed <- timed_call(loops[[route]], timeout_seconds)
          }
          at <- at + 1L
          rows[[at]] <- timing_row(
            "K", "batch", shape, n, route, rep, pos, ncores = 1L,
            timed = timed, value = timed$value,
            extra = list(traits = p, traits_per_sec = if (is.finite(timed$elapsed_s) && timed$elapsed_s > 0) p / timed$elapsed_s else NA_real_,
                         per_trait_s = if (is.finite(timed$elapsed_s) && timed$elapsed_s > 0) timed$elapsed_s / p else NA_real_),
            resource = if (identical(route, "reference_column_loop")) reference_resource else NULL
          )
        }
      }
    }
    rm(tree, context, X); gc(FALSE)
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  write_csv(out, "stage2c_k_batch_timings.csv")
  write_csv(summarize_timings(out, c("method", "workload", "shape", "n", "route", "traits", "ncores")),
            "stage2c_k_batch_summary.csv")
  out
}

make_permutation_matrix <- function(n, nsim) {
  if (as.double(n) * as.double(nsim) > 5e6) return(NULL)
  private <- private_get(".permutation_matrix")
  if (!is.function(private)) return(NULL)
  private(NULL, n = n, nsim = nsim, include_observed = FALSE)
}

run_k_permutation <- function() {
  rows <- list(); phase_rows_all <- list(); at <- 0L; pt <- 0L
  for (shape in perm_shapes) for (n in perm_n) {
    message("[stage2c] K permutation shape=", shape, " n=", n)
    tree <- make_fixture(shape, n); trait <- make_continuous_trait(tree)
    context <- timed_call(function() fastphylosig::prepare_tree(tree), timeout_seconds)$value
    if (is.null(context)) next
    for (nsim in perm_nsim) {
      funcs <- list(
        raw = function() call_fast_k("raw", tree, context, trait, test = TRUE, nsim = nsim),
        prepared = function() call_fast_k("prepared", tree, context, trait, test = TRUE, nsim = nsim)
      )
      routes <- names(funcs); reps <- repeats_for(n)
      for (route in routes) invisible(timed_call(funcs[[route]], timeout_seconds))
      for (rep in seq_len(reps)) {
        order <- if (rep %% 2L == 0L) rev(routes) else routes
        for (pos in seq_along(order)) {
          route <- order[[pos]]; set.seed(830000L + n + nsim + rep)
          timed <- timed_call(funcs[[route]], timeout_seconds); at <- at + 1L
          rows[[at]] <- timing_row(
            "K", "permutation", shape, n, route, rep, pos, nsim = nsim,
            timed = timed, value = timed$value,
            extra = list(permutation_mode = "internal_rng")
          )
        }
      }
      # One low-overhead phase probe per route/cell.  It is intentionally
      # separate from authoritative timings and may overlap inclusive helpers.
      for (route in routes) {
        set.seed(835000L + n + nsim)
        probe <- phase_probe(
          funcs[[route]],
          c(".prepare_analysis", ".validated_context", ".validate_prepared_context",
            ".v2_context_evidence", ".tree_fingerprint", ".v2_protected_snapshot",
            ".prepared_tree_subset", ".analysis_na_groups", ".as_trait_matrix",
            ".permutation_matrix", "fast_k_tree_batch_cpp",
            "fast_k_tree_permutation_cpp", ".k_permutation_mcse")
        )
        known <- phase_rows(probe, "K", "permutation", shape, n, route,
                            nsim = nsim, mode = "internal_rng")
        if (nrow(known)) {
          known$probe_status <- probe$status; known$probe_warning <- probe$warning_text
          known$probe_error <- probe$error_message; known$instrumented_bindings <- probe$instrumented_bindings
          pt <- pt + 1L; phase_rows_all[[pt]] <- known
          remainder <- data.frame(method = "K", workload = "permutation", shape = shape,
            n = n, route = route, nsim = nsim, ncores = 1L, mode = "internal_rng",
            phase = "result_construction_and_unattributed_remainder",
            elapsed_s = if (is.finite(probe$elapsed_s)) max(0, probe$elapsed_s - sum(probe$phase$elapsed_s)) else NA_real_,
            share_total = NA_real_, calls = NA_integer_,
            measurement_scope = "derived residual; phase helpers are inclusive",
            note = "C++ permutation call includes internal RNG, trait gather, null K and exceedance reduction",
            probe_status = probe$status, probe_warning = probe$warning_text,
            probe_error = probe$error_message, instrumented_bindings = probe$instrumented_bindings,
            stringsAsFactors = FALSE)
          pt <- pt + 1L; phase_rows_all[[pt]] <- remainder
        }
        # A controlled matrix is used where its bounded memory is explicit;
        # this separates R permutation generation from the C++ null kernel.
        if (as.double(n) * as.double(nsim) <= 5e6) {
          set.seed(836000L + n + nsim)
          controlled <- make_permutation_matrix(n, nsim)
          if (!is.null(controlled)) {
            controlled_fun <- function() call_fast_k(
              route, tree, context, trait, test = TRUE, nsim = nsim,
              permutations = controlled
            )
            cp <- phase_probe(controlled_fun, c(
              ".prepare_analysis", ".validated_context", ".validate_prepared_context",
              ".v2_context_evidence", ".tree_fingerprint", ".v2_protected_snapshot",
              ".prepared_tree_subset", ".analysis_na_groups", ".as_trait_matrix",
              ".permutation_matrix", "fast_k_tree_permutation_cpp", ".k_permutation_mcse"
            ))
            cr <- phase_rows(cp, "K", "permutation", shape, n, route, nsim,
                             mode = "controlled_permutations",
                             note = "controlled matrix generated outside public timed call")
            if (nrow(cr)) {
              cr$probe_status <- cp$status
              cr$probe_warning <- cp$warning_text
              cr$probe_error <- cp$error_message
              cr$instrumented_bindings <- cp$instrumented_bindings
              pt <- pt + 1L; phase_rows_all[[pt]] <- cr
            }
          }
        }
      }
    }
    rm(tree, trait, context); gc(FALSE)
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  phases <- if (length(phase_rows_all)) do.call(rbind, phase_rows_all) else data.frame()
  write_csv(out, "stage2c_k_permutation_timings.csv")
  write_csv(summarize_timings(out, c("method", "workload", "shape", "n", "route", "nsim", "ncores")),
            "stage2c_k_permutation_summary.csv")
  write_csv(phases, "stage2c_k_permutation_phases.csv")
  out
}

run_prepared_k_validation_audit <- function() {
  target_n <- suppressWarnings(as.integer(Sys.getenv(
    "FASTPHYLOSIG_STAGE2C_PREPARED_K_N", unset = if (smoke) "100" else "20000"
  )))
  if (length(target_n) != 1L || is.na(target_n) || target_n < 2L) {
    stop("FASTPHYLOSIG_STAGE2C_PREPARED_K_N must be an integer >= 2.", call. = FALSE)
  }
  rows <- list(); at <- 0L
  for (shape in single_shapes) {
    tree <- make_fixture(shape, target_n)
    trait <- make_continuous_trait(tree)
    context <- timed_call(function() fastphylosig::prepare_tree(tree), timeout_seconds)$value
    if (is.null(context)) next
    set.seed(860000L + target_n)
    direct <- timed_call(function() call_fast_k("prepared", tree, context, trait), timeout_seconds)
    set.seed(860000L + target_n)
    probe <- phase_probe(
      function() call_fast_k("prepared", tree, context, trait),
      c(".prepare_analysis", ".validated_context", ".validate_prepared_context",
        ".v2_context_evidence", ".tree_fingerprint", ".v2_protected_snapshot",
        ".prepared_tree_subset", ".as_trait_matrix", "fast_k_tree_batch_cpp")
    )
    z <- phase_rows(probe, "K", "prepared_v2_validation", shape, target_n,
                    "prepared", mode = "diagnostic",
                    note = "inclusive helper timings; not authoritative wall-time")
    z$direct_elapsed_s <- direct$elapsed_s
    z$probe_elapsed_s <- probe$elapsed_s
    z$probe_over_direct_ratio <- if (is.finite(direct$elapsed_s) && direct$elapsed_s > 0) {
      probe$elapsed_s / direct$elapsed_s
    } else NA_real_
    z$probe_status <- probe$status
    z$probe_warning <- probe$warning_text
    z$probe_error <- probe$error_message
    z$instrumented_bindings <- probe$instrumented_bindings
    rows[[length(rows) + 1L]] <- z
    rm(tree, trait, context); gc(FALSE)
  }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  write_csv(out, "stage2c_prepared_k_v2_phases.csv")
  if (nrow(out)) {
    summary <- aggregate(
      cbind(elapsed_s, direct_elapsed_s, probe_elapsed_s) ~ method + workload +
        shape + n + route + phase,
      out, FUN = function(x) stats::median(x, na.rm = TRUE)
    )
  } else summary <- data.frame()
  write_csv(summary, "stage2c_prepared_k_v2_summary.csv")
  out
}

run_k_correctness_guard <- function() {
  n <- if (smoke) 50L else 100L
  tree <- make_fixture("random", n)
  trait <- make_continuous_trait(tree, offset = 17L)
  before_tree <- serialize(tree, NULL, version = 3L)
  before_trait <- serialize(trait, NULL, version = 3L)
  context <- fastphylosig::prepare_tree(tree)
  rows <- list(); at <- 0L
  add <- function(check, pass, detail) {
    rows[[length(rows) + 1L]] <<- data.frame(
      check = check, pass = isTRUE(pass), detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }
  get_k <- function(value) {
    if (is.list(value) && !is.null(value$K)) return(as.numeric(value$K[[1L]]))
    if (is.data.frame(value) && "K_fast" %in% names(value)) return(as.numeric(value$K_fast[[1L]]))
    as.numeric(value[[1L]])
  }
  compare <- function(left, right, tolerance = 1e-10) {
    if (left$status != "ok" || right$status != "ok") return(FALSE)
    isTRUE(all.equal(get_k(left$value), get_k(right$value), tolerance = tolerance))
  }
  raw <- timed_call(function() call_fast_k("raw", tree, context, trait), timeout_seconds)
  prepared <- timed_call(function() call_fast_k("prepared", tree, context, trait), timeout_seconds)
  add("raw_prepared_K_parity", compare(raw, prepared),
      sprintf("raw=%s prepared=%s", get_k(raw$value), get_k(prepared$value)))
  add("raw_prepared_status_parity", identical(raw$status, prepared$status),
      paste(raw$status, prepared$status, sep = " / "))

  set.seed(861101L)
  permutations <- make_permutation_matrix(n, if (smoke) 19L else 39L)
  controlled_raw <- timed_call(function() call_fast_k(
    "raw", tree, context, trait, test = TRUE,
    nsim = nrow(permutations), permutations = permutations
  ), timeout_seconds)
  controlled_prepared <- timed_call(function() call_fast_k(
    "prepared", tree, context, trait, test = TRUE,
    nsim = nrow(permutations), permutations = permutations
  ), timeout_seconds)
  controlled_equal <- if (controlled_raw$status == "ok" && controlled_prepared$status == "ok") {
    all.equal(get_k(controlled_raw$value), get_k(controlled_prepared$value), tolerance = 1e-10) &&
      all.equal(result_p(controlled_raw$value, "K"), result_p(controlled_prepared$value, "K"), tolerance = 1e-10)
  } else FALSE
  add("controlled_permutation_raw_prepared_parity", isTRUE(controlled_equal),
      "same supplied permutation matrix")
  add("controlled_permutation_status_parity",
      identical(controlled_raw$status, controlled_prepared$status),
      paste(controlled_raw$status, controlled_prepared$status, sep = " / "))

  shuffled <- ape::reorder.phylo(tree, order = "cladewise")
  shuffled_result <- timed_call(function() call_fast_k("raw", shuffled, context, trait), timeout_seconds)
  add("edge_order_representation_parity", compare(raw, shuffled_result),
      "cladewise versus fixed random fixture")
  add("tree_input_immutability", identical(serialize(tree, NULL, version = 3L), before_tree),
      "tree serialization unchanged")
  add("trait_input_immutability", identical(serialize(trait, NULL, version = 3L), before_trait),
      "trait serialization unchanged")
  guard <- if (length(rows)) do.call(rbind, rows) else data.frame()
  write_csv(guard, "stage2c_k_correctness.csv")
  if (nrow(guard) && any(!guard$pass)) {
    stop("K correctness guard failed; benchmark interpretation stopped.", call. = FALSE)
  }
  guard
}

git_commit_for <- function() {
  explicit <- Sys.getenv("FASTPHYLOSIG_STAGE2C_COMMIT", unset = "")
  if (nzchar(explicit)) return(explicit)
  if (!dir.exists(file.path(repo, ".git")) &&
      !file.exists(file.path(repo, ".git"))) return("<not available>")
  value <- tryCatch(system2("git", c("-C", repo, "rev-parse", "HEAD"),
                            stdout = TRUE, stderr = TRUE), error = function(e) character())
  value <- trimws(value)
  if (length(value) == 1L && grepl("^[0-9a-fA-F]{40}$", value)) value else "<not available>"
}

write_provenance <- function() {
  blas <- tryCatch(extSoftVersion()[["BLAS"]], error = function(e) NA_character_)
  compiler <- if (length(R.version$compiler)) as.character(R.version$compiler[[1L]]) else NA_character_
  detected <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_)
  data.frame(
    script = "run_stage2c.R", stage = "Stage 2C K-only post-V2 reprofile",
    mode = if (smoke) "quick_smoke" else "formal",
    source_commit = git_commit_for(), package_version = as.character(utils::packageVersion("fastphylosig")),
    R_version = R.version.string, R_platform = R.version$platform, R_arch = R.version$arch,
    OS = Sys.info()[["sysname"]], OS_release = Sys.info()[["release"]], machine = Sys.info()[["machine"]],
    compiler = compiler, OpenMP_status = env_or("FASTPHYLOSIG_OPENMP_STATUS", "not independently exposed"),
    BLAS = blas, LC_ALL = Sys.getenv("LC_ALL", ""), LC_COLLATE = Sys.getlocale("LC_COLLATE"),
    hardware = paste(Sys.info()[c("sysname", "release", "machine")], collapse = "/"),
    logical_cores = detected, repo = repo, output_dir = out_dir,
    n_single = paste(single_n, collapse = ","), n_batch = paste(batch_n, collapse = ","),
    n_permutation = paste(perm_n, collapse = ","), nsim_permutation = paste(perm_nsim, collapse = ","),
    reference_dense_cap_bytes = reference_dense_cap_bytes,
    reference_batch_work_cap_bytes = reference_batch_work_cap_bytes,
    reference_single_cubic_cap = reference_single_cubic_cap,
    reference_batch_cubic_cap = reference_batch_cubic_cap,
    reference_policy = paste(
      "dense bytes plus estimated repeated cubic work are bounded;",
      "preflight failures are resource_limit and are not extrapolated"
    ),
    repeats = if (smoke) "1" else "n<10000:10; n>=10000:5",
    timeout_seconds = if (is.finite(timeout_seconds)) timeout_seconds else "Inf",
    warmup = "one warmup per route/cell", alternation = "odd pairs raw->prepared; even pairs prepared->raw",
    benchmark_parallelism = "serial; no concurrent benchmark R processes",
    source_change_scope = "harness only; no production estimator change",
    stringsAsFactors = FALSE
  ) -> out
  write_csv(out, "stage2c_provenance.csv")
  capture.output(sessionInfo(), file = file.path(out_dir, "stage2c_sessionInfo.txt"))
  out
}

median_cell <- function(dat, route, shape, n, nsim = NA_integer_, traits = NA_integer_) {
  if (is.null(dat) || !nrow(dat)) return(NA_real_)
  keep <- dat$route == route & dat$shape == shape & dat$n == n & dat$status == "ok"
  if (!is.na(nsim) && "nsim" %in% names(dat)) keep <- keep & dat$nsim == nsim
  if (!is.na(traits) && "traits" %in% names(dat)) keep <- keep & dat$traits == traits
  z <- dat$elapsed_s[keep & is.finite(dat$elapsed_s)]
  if (length(z)) stats::median(z) else NA_real_
}

build_k_unified_budget <- function(single, batch, permutation, validation) {
  rows <- list(); at <- 0L
  if (nrow(single)) for (shape in unique(single$shape)) for (n in unique(single$n)) {
    raw <- median_cell(single, "B_fast_k_raw", shape, n)
    prep <- median_cell(single, "prepare_tree", shape, n)
    prepared <- median_cell(single, "C_fast_k_prepared", shape, n)
    reference <- median_cell(single, "A_reference", shape, n)
    at <- at + 1L
    rows[[at]] <- data.frame(
      method = "K", workload = "single_trait_crossover", shape = shape, n = n,
      raw_total = raw, prepare = prep, prepared_total = prep + prepared,
      validation = NA_real_, kernel_optimizer = NA_real_, simulation = 0,
      startup = 0, packaging = NA_real_, reference_time = reference,
      speedup = if (is.finite(reference) && is.finite(raw) && raw > 0) reference / raw else NA_real_,
      memory_context_cost = NA_real_, stringsAsFactors = FALSE
    )
  }
  if (nrow(batch)) for (shape in unique(batch$shape)) for (n in unique(batch$n))
    for (p in sort(unique(batch$traits))) {
      raw <- median_cell(batch, "raw_matrix", shape, n, traits = p)
      prepared <- median_cell(batch, "prepared_matrix", shape, n, traits = p)
      reference <- median_cell(batch, "reference_column_loop", shape, n, traits = p)
      at <- at + 1L
      rows[[at]] <- data.frame(
        method = "K", workload = paste0("batch_", p, "_traits"), shape = shape, n = n,
        raw_total = raw, prepare = NA_real_, prepared_total = prepared,
        validation = NA_real_, kernel_optimizer = NA_real_, simulation = 0,
        startup = 0, packaging = NA_real_, reference_time = reference,
        speedup = if (is.finite(reference) && is.finite(raw) && raw > 0) reference / raw else NA_real_,
        memory_context_cost = NA_real_, stringsAsFactors = FALSE
      )
    }
  if (nrow(permutation)) for (shape in unique(permutation$shape))
    for (n in unique(permutation$n)) for (nsim in sort(unique(permutation$nsim))) {
      raw <- median_cell(permutation, "raw", shape, n, nsim)
      prepared <- median_cell(permutation, "prepared", shape, n, nsim)
      at <- at + 1L
      rows[[at]] <- data.frame(
        method = "K", workload = paste0("permutation_nsim_", nsim), shape = shape, n = n,
        raw_total = raw, prepare = NA_real_, prepared_total = prepared,
        validation = NA_real_, kernel_optimizer = NA_real_, simulation = NA_real_,
        startup = 0, packaging = NA_real_, reference_time = NA_real_, speedup = NA_real_,
        memory_context_cost = NA_real_, stringsAsFactors = FALSE
      )
    }
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  write_csv(out, "stage2c_unified_budget.csv")
  out
}

build_k_candidate_ranking <- function(permutation, phases, validation) {
  simulation_fraction <- NA_real_
  if (nrow(phases)) {
    total <- phases[phases$phase == "fast_k_tree_permutation_cpp" & is.finite(phases$elapsed_s), "share_total"]
    if (length(total)) simulation_fraction <- stats::median(total, na.rm = TRUE)
  }
  validation_fraction <- NA_real_
  if (nrow(validation)) {
    total <- validation[validation$phase == ".validate_prepared_context" & is.finite(validation$elapsed_s), "share_total"]
    if (length(total)) validation_fraction <- stats::median(total, na.rm = TRUE)
  }
  out <- data.frame(
    candidate = c("K_permutation_null_engine", "prepared_exact_validation_efficiency"),
    hotspot_share = c(simulation_fraction, validation_fraction),
    ideal_end_to_end_saving = 100 * c(simulation_fraction, validation_fraction),
    measured_practical_saving = NA_real_, correctness_risk = c("medium", "high"),
    engineering_risk = c("medium", "high"), authorization = c(
      if (is.finite(simulation_fraction) && simulation_fraction >= 0.5) "AUTHORIZED_FOR_CANDIDATE_DESIGN" else "NO-GO",
      "NO-GO; exact-equivalent design not established"
    ), scope = "profile/rank only; no candidate implemented",
    stringsAsFactors = FALSE
  )
  write_csv(out, "stage2c_k_candidate_ranking.csv")
  out
}

main <- function() {
  message("[stage2c] K-only ", if (smoke) "quick smoke" else "formal", " run; output=", out_dir)
  write_provenance()
  guard <- run_k_correctness_guard()
  single <- run_k_single()
  batch <- run_k_batch()
  permutation <- run_k_permutation()
  validation <- run_prepared_k_validation_audit()
  phases <- if (file.exists(file.path(out_dir, "stage2c_k_permutation_phases.csv"))) {
    utils::read.csv(file.path(out_dir, "stage2c_k_permutation_phases.csv"), stringsAsFactors = FALSE)
  } else data.frame()
  unified <- build_k_unified_budget(single, batch, permutation, validation)
  ranking <- build_k_candidate_ranking(permutation, phases, validation)
  write_csv(data.frame(
    stage = "Stage 2C K-only", status = "PASS", correctness_guard = "PASS",
    production_code_changed = "NO", candidate_implemented = "NO",
    formal_run = !smoke, next_stage = "STOP; do not implement candidate automatically",
    stringsAsFactors = FALSE
  ), "stage2c_status.csv")
  invisible(list(guard = guard, single = single, batch = batch,
                 permutation = permutation, validation = validation,
                 unified = unified, ranking = ranking))
}

main()
