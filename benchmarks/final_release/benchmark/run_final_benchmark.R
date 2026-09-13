#!/usr/bin/env Rscript

# Final 0.2.0 public benchmark and memory characterization.
# This is audit-only: it never edits package source, tests, metadata, or API.

argv <- commandArgs(trailingOnly = TRUE)
if (!length(argv) && nzchar(Sys.getenv("FASTPHYLOSIG_FINAL_ARGS", unset = ""))) {
  argv <- trimws(strsplit(Sys.getenv("FASTPHYLOSIG_FINAL_ARGS"), "[[:space:]]+")[[1L]])
}

parse_args <- function(values) {
  modes <- intersect(values, c("--smoke", "--formal", "--memory-only"))
  if (length(modes) != 1L) {
    stop("Use exactly one of --smoke, --formal, or --memory-only.",
         call. = FALSE)
  }
  positional <- values[!values %in% modes]
  if (length(positional) > 2L) {
    stop("usage: run_final_benchmark.R [--smoke|--formal|--memory-only] [repo] [out].",
         call. = FALSE)
  }
  repo <- if (length(positional)) positional[[1L]] else getwd()
  out <- if (length(positional) >= 2L) positional[[2L]] else
    file.path(repo, "benchmarks", "final_release", "benchmark", "results")
  list(mode = sub("^--", "", modes[[1L]]),
       repo = normalizePath(repo, winslash = "/", mustWork = FALSE),
       out = normalizePath(out, winslash = "/", mustWork = FALSE))
}

args <- parse_args(argv)
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
formal <- identical(args$mode, "formal")
memory_only <- identical(args$mode, "memory-only")

env_or <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

as_int <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.integer(env_or(name, as.character(default))))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(name, " must be an integer >= ", minimum, call. = FALSE)
  }
  value
}

csv_ints <- function(name, default, minimum = 1L) {
  value <- env_or(name, paste(default, collapse = ","))
  out <- suppressWarnings(as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])))
  if (!length(out) || anyNA(out) || any(out < minimum)) {
    stop(name, " must contain comma-separated integers >= ", minimum,
         call. = FALSE)
  }
  unique(out)
}

csv_chars <- function(name, default) {
  value <- env_or(name, paste(default, collapse = ","))
  out <- unique(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (!length(out) || any(!nzchar(out))) stop(name, " is empty.", call. = FALSE)
  out
}

repeats <- as_int("FASTPHYLOSIG_FINAL_REPEATS", if (formal) 10L else 1L,
                  minimum = if (formal) 10L else 1L)
sizes <- csv_ints("FASTPHYLOSIG_FINAL_SIZES",
                  if (formal) c(500L, 2000L, 5000L) else 100L, 2L)
memory_sizes <- csv_ints(
  "FASTPHYLOSIG_FINAL_MEMORY_SIZES",
  if (formal) c(500L, 1000L, 5000L, 10000L, 20000L) else 500L, 2L
)
benchmark_shapes <- csv_chars(
  "FASTPHYLOSIG_FINAL_SHAPES",
  if (formal) c("balanced", "random", "pectinate") else "balanced"
)
memory_shapes <- csv_chars("FASTPHYLOSIG_FINAL_MEMORY_SHAPES", "balanced")
matrix_traits <- csv_ints("FASTPHYLOSIG_FINAL_MATRIX_TRAITS",
                          if (formal) c(8L, 32L, 100L) else 8L, 1L)
permutation_sizes <- csv_ints(
  "FASTPHYLOSIG_FINAL_PERMUTATION_SIZES",
  if (formal) c(500L, 2000L, 5000L) else 100L, 2L
)
permutation_nsim <- as_int("FASTPHYLOSIG_FINAL_PERMUTATION_NSIM",
                           if (formal) 999L else 19L, minimum = 1L)
heavy_budget <- as.numeric(env_or("FASTPHYLOSIG_FINAL_HEAVY_BUDGET_S",
                                  if (formal) "600" else "60"))
light_budget <- as.numeric(env_or("FASTPHYLOSIG_FINAL_LIGHT_BUDGET_S",
                                  if (formal) "120" else "60"))
memory_budget <- as.numeric(env_or("FASTPHYLOSIG_FINAL_MEMORY_BUDGET_S",
                                   if (formal) "300" else "120"))
if (any(!is.finite(c(heavy_budget, light_budget, memory_budget))) ||
    any(c(heavy_budget, light_budget, memory_budget) <= 0)) {
  stop("per-cell budgets must be positive finite seconds.", call. = FALSE)
}

git_value <- function(repo, git_args) {
  value <- tryCatch(system2("git", c("-C", repo, git_args),
                            stdout = TRUE, stderr = FALSE),
                    error = function(e) character())
  if (length(value)) paste(value, collapse = "\n") else "UNAVAILABLE"
}

compiler_value <- function() {
  r_bin <- file.path(R.home("bin"), "R.exe")
  value <- tryCatch(system2(r_bin, c("CMD", "config", "CXX17"),
                            stdout = TRUE, stderr = FALSE),
                    error = function(e) character())
  if (length(value)) paste(value, collapse = " ") else "UNAVAILABLE"
}

sys <- Sys.info()
provenance <- list(
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  commit = env_or("FASTPHYLOSIG_FINAL_COMMIT", git_value(args$repo, "rev-parse HEAD")),
  branch = env_or("FASTPHYLOSIG_FINAL_BRANCH", git_value(args$repo, "branch --show-current")),
  package_version = NA_character_,
  r_version = R.version.string,
  platform = R.version$platform,
  os_type = unname(sys[["sysname"]]),
  os_release = unname(sys[["release"]]),
  machine = unname(sys[["machine"]]),
  nodename = unname(sys[["nodename"]]),
  compiler_cxx17 = compiler_value(),
  openmp = if (isTRUE(capabilities("openmp"))) "enabled" else "disabled_or_unavailable",
  blas = if ("BLAS" %in% names(extSoftVersion()))
    unname(extSoftVersion()[["BLAS"]]) else "UNAVAILABLE",
  lapack = if ("LAPACK" %in% names(extSoftVersion()))
    unname(extSoftVersion()[["LAPACK"]]) else "UNAVAILABLE",
  locale = paste(Sys.getlocale(), collapse = " | "),
  logical_cores = parallel::detectCores(logical = TRUE),
  mode = args$mode,
  repeats = repeats,
  light_budget_s = light_budget,
  heavy_budget_s = heavy_budget,
  memory_budget_s = memory_budget
)

load_package <- function(repo) {
  library_dir <- env_or("FASTPHYLOSIG_FINAL_LIBRARY")
  if (nzchar(library_dir) && dir.exists(file.path(library_dir, "fastphylosig"))) {
    library_dir <- normalizePath(library_dir, winslash = "/", mustWork = TRUE)
    .libPaths(c(library_dir, .libPaths()))
    suppressPackageStartupMessages(library(
      "fastphylosig", lib.loc = library_dir, character.only = TRUE
    ))
  } else if (requireNamespace("pkgload", quietly = TRUE)) {
    suppressPackageStartupMessages(pkgload::load_all(repo, quiet = TRUE))
  } else {
    stop("Use FASTPHYLOSIG_FINAL_LIBRARY or install pkgload.", call. = FALSE)
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("ape is required for fixed fixtures.", call. = FALSE)
  }
  provenance$package_version <<- as.character(utils::packageVersion("fastphylosig"))
  invisible(TRUE)
}

load_package(args$repo)

restore_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(code)
}

make_balanced <- function(n) {
  root <- n + 1L
  edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  q_node <- q_lo <- q_hi <- integer(n - 1L)
  q_node[[1L]] <- root; q_lo[[1L]] <- 1L; q_hi[[1L]] <- n
  head <- 1L; tail <- 1L; next_internal <- root + 1L; edge_i <- 0L
  while (head <= tail) {
    node <- q_node[[head]]; lo <- q_lo[[head]]; hi <- q_hi[[head]]; head <- head + 1L
    mid <- floor((lo + hi) / 2L)
    child_lo <- c(lo, mid + 1L); child_hi <- c(mid, hi)
    for (j in 1:2) {
      child <- if (child_lo[[j]] == child_hi[[j]]) child_lo[[j]] else {
        new_node <- next_internal; next_internal <- next_internal + 1L
        tail <- tail + 1L; q_node[[tail]] <- new_node
        q_lo[[tail]] <- child_lo[[j]]; q_hi[[tail]] <- child_hi[[j]]; new_node
      }
      edge_i <- edge_i + 1L; edge[edge_i, ] <- c(node, child)
    }
  }
  tree <- list(edge = edge[seq_len(edge_i), , drop = FALSE],
               tip.label = paste0("sp", seq_len(n)),
               edge.length = 0.5 + seq_len(edge_i) / edge_i,
               Nnode = n - 1L)
  class(tree) <- "phylo"
  ape::reorder.phylo(tree, order = "postorder")
}

make_pectinate <- function(n) {
  root <- n + 1L
  if (n == 2L) edge <- matrix(c(root, 1L, root, 2L), ncol = 2L, byrow = TRUE)
  else {
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

make_random <- function(n, seed) {
  restore_seed(seed + as.integer(n), {
    tree <- ape::reorder.phylo(ape::rtree(n), order = "postorder")
    tree$tip.label <- paste0("sp", seq_len(n))
    tree$edge.length <- 0.5 + seq_len(nrow(tree$edge)) / nrow(tree$edge)
    tree
  })
}

make_tree <- function(n, shape, seed = 206020L) {
  switch(shape, balanced = make_balanced(n), random = make_random(n, seed),
         pectinate = make_pectinate(n), stop("unknown shape: ", shape,
                                             call. = FALSE))
}

make_continuous <- function(tree, k = 1L, seed = 702020L) {
  restore_seed(seed, {
    values <- matrix(stats::rnorm(length(tree$tip.label) * k),
                     nrow = length(tree$tip.label), ncol = k)
    rownames(values) <- tree$tip.label
    colnames(values) <- paste0("trait", seq_len(k))
    if (k == 1L) { out <- as.numeric(values[, 1L]); names(out) <- tree$tip.label; out }
    else values
  })
}

make_binary <- function(tree) {
  out <- as.numeric(seq_along(tree$tip.label) %% 2L == 0L)
  names(out) <- tree$tip.label
  out
}

numeric_close <- function(x, y, tolerance = 1e-7) {
  if (length(x) != length(y)) return(FALSE)
  x <- as.numeric(x); y <- as.numeric(y)
  same_na <- is.na(x) & is.na(y); finite <- is.finite(x) & is.finite(y)
  scale <- pmax(1, abs(x), abs(y))
  all(same_na | (finite & abs(x - y) <= tolerance * scale))
}

`%||%` <- function(x, y) if (is.null(x)) y else x
extract_k <- function(value) if (is.list(value)) as.numeric(value$K %||% value$estimate) else as.numeric(value)
extract_lambda <- function(value) if (is.list(value)) as.numeric(value$lambda %||% value$estimate) else as.numeric(value)
extract_d <- function(value, name) {
  if (!is.list(value)) return(NA_real_)
  keys <- switch(name,
    observed = c("DEstimate", "estimate"), random = c("P_random", "p_random"),
    brownian = c("P_Brownian", "p_brownian"),
    mcse_random = c("MCSE_P_random", "MCSE_random", "mcse_random"),
    mcse_brownian = c("MCSE_P_Brownian", "MCSE_Brownian", "mcse_brownian"),
    character())
  for (key in keys) if (!is.null(value[[key]])) return(as.numeric(value[[key]])[[1L]])
  NA_real_
}
compare_d <- function(x, y) {
  fields <- c("observed", "random", "brownian", "mcse_random", "mcse_brownian")
  all(vapply(fields, function(f) numeric_close(extract_d(x, f), extract_d(y, f)), logical(1L)))
}

capture_call <- function(fun, seed = NULL, budget_s = Inf) {
  warnings <- character(); status <- "PASS"; error_message <- ""
  started <- proc.time()[["elapsed"]]
  body <- if (is.null(seed)) fun else function() restore_seed(seed, fun())
  value <- tryCatch({
    setTimeLimit(elapsed = budget_s, transient = TRUE)
    withCallingHandlers(body(), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    msg <- conditionMessage(e)
    status <<- if (grepl("time limit|elapsed", msg, ignore.case = TRUE))
      "RESOURCE_LIMIT" else "ERROR"
    error_message <<- msg; NULL
  }, finally = setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))
  list(value = value, elapsed_s = as.numeric(proc.time()[["elapsed"]] - started),
       status = status, warning_count = length(unique(warnings)),
       warning_text = paste(unique(warnings), collapse = " | "),
       error_message = error_message)
}

row_template <- function(workload, method, shape, n, traits, nsim, route,
                         rep_id, timed, ncores, fixture_seed, guard_id,
                         reference_name = NA_character_, reference_status = NA_character_) {
  data.frame(timestamp_utc = provenance$timestamp_utc, commit = provenance$commit,
    package_version = provenance$package_version, r_version = provenance$r_version,
    os_type = provenance$os_type, compiler_cxx17 = provenance$compiler_cxx17,
    openmp = provenance$openmp, blas = provenance$blas, lapack = provenance$lapack,
    locale = provenance$locale, workload = workload,
    method = method, shape = shape, n = as.integer(n), traits = as.integer(traits),
    nsim = as.integer(nsim), route = route, repeat_id = as.integer(rep_id),
    elapsed_s = timed$elapsed_s, status = timed$status,
    warning_count = as.integer(timed$warning_count), warning_text = timed$warning_text,
    error_message = timed$error_message, ncores = as.integer(ncores),
    fixture_seed = as.integer(fixture_seed), guard_id = guard_id,
    reference_name = reference_name, reference_status = reference_status,
    stringsAsFactors = FALSE)
}

timing_path <- file.path(args$out, "final_benchmark.csv")
memory_path <- file.path(args$out, "memory_characterization.csv")
status_path <- file.path(args$out, "benchmark_status.csv")
for (path in c(timing_path, memory_path, status_path)) if (file.exists(path)) file.remove(path)
append_rows <- function(rows, path) {
  if (!length(rows)) return(invisible(NULL))
  first <- !file.exists(path)
  utils::write.table(rows, path, sep = ",", row.names = FALSE, col.names = first,
                     append = !first, qmethod = "double", na = "NA")
}
record_status <- function(check, status, detail = "", workload = "") {
  append_rows(data.frame(timestamp_utc = provenance$timestamp_utc,
    commit = provenance$commit, package_version = provenance$package_version,
    check = check, status = status, workload = workload, detail = as.character(detail),
    stringsAsFactors = FALSE), status_path)
}
record_status("provenance", "PASS", paste(names(provenance), provenance,
                                            sep = "=", collapse = " | "))
record_status("scope", "PASS", "single-machine; cross-platform qualification not run")

run_guard <- function(tree, ctx, cont, binary, shape, n, fixture_seed) {
  bytes_before <- serialize(tree, NULL, version = 3)
  raw_k <- capture_call(function() fastphylosig::fast_k(tree, cont, test = FALSE,
    verbose = FALSE, progress = FALSE), budget_s = light_budget)
  prep_k <- capture_call(function() fastphylosig::fast_k(ctx, cont, test = FALSE,
    verbose = FALSE, progress = FALSE), budget_s = light_budget)
  raw_l <- capture_call(function() fastphylosig::fast_lambda(tree, cont,
    verbose = FALSE, progress = FALSE), budget_s = light_budget)
  prep_l <- capture_call(function() fastphylosig::fast_lambda(ctx, cont,
    verbose = FALSE, progress = FALSE), budget_s = light_budget)
  # Use fixed null matrices for the non-stochastic raw/prepared parity guard.
  # This avoids conflating a route comparison with the independent streaming
  # RNG path; stochastic replay is checked separately below.
  fixed_random <- matrix(rep(binary, 3L), nrow = length(binary), ncol = 3L)
  fixed_brownian <- matrix(rep(1 - binary, 3L), nrow = length(binary), ncol = 3L)
  raw_d <- capture_call(function() fastphylosig::fast_d(tree, binary, test = FALSE,
    nsim = 3L, random_states = fixed_random, brownian_states = fixed_brownian,
    verbose = FALSE, progress = FALSE), budget_s = light_budget)
  prep_d <- capture_call(function() fastphylosig::fast_d(ctx, binary, test = FALSE,
    nsim = 3L, random_states = fixed_random, brownian_states = fixed_brownian,
    verbose = FALSE, progress = FALSE), budget_s = light_budget)
  immutable <- identical(bytes_before, serialize(tree, NULL, version = 3))
  parity <- raw_k$status == "PASS" && prep_k$status == "PASS" &&
    raw_l$status == "PASS" && prep_l$status == "PASS" && raw_d$status == "PASS" &&
    prep_d$status == "PASS" && numeric_close(extract_k(raw_k$value), extract_k(prep_k$value)) &&
    numeric_close(extract_lambda(raw_l$value), extract_lambda(prep_l$value)) &&
    compare_d(raw_d$value, prep_d$value)
  replay_k <- capture_call(function() fastphylosig::fast_k(ctx, cont, test = TRUE,
    nsim = 19L, verbose = FALSE, progress = FALSE), seed = fixture_seed + 91L,
    budget_s = light_budget)
  replay_k2 <- capture_call(function() fastphylosig::fast_k(ctx, cont, test = TRUE,
    nsim = 19L, verbose = FALSE, progress = FALSE), seed = fixture_seed + 91L,
    budget_s = light_budget)
  replay_d <- capture_call(function() fastphylosig::fast_d(ctx, binary, test = TRUE,
    nsim = 19L, return_sim = FALSE, keep_null = FALSE, ncores = 1L,
    verbose = FALSE, progress = FALSE), seed = fixture_seed + 97L,
    budget_s = light_budget)
  replay_d2 <- capture_call(function() fastphylosig::fast_d(ctx, binary, test = TRUE,
    nsim = 19L, return_sim = FALSE, keep_null = FALSE, ncores = 1L,
    verbose = FALSE, progress = FALSE), seed = fixture_seed + 97L,
    budget_s = light_budget)
  replay <- replay_k$status == "PASS" && replay_k2$status == "PASS" &&
    numeric_close(extract_k(replay_k$value), extract_k(replay_k2$value)) &&
    identical(replay_k$value$P, replay_k2$value$P) && replay_d$status == "PASS" &&
    replay_d2$status == "PASS" && compare_d(replay_d$value, replay_d2$value)
  matrix <- make_continuous(tree, 2L, seed = fixture_seed + 11L)
  matrix_call <- capture_call(function() fastphylosig::fast_k(ctx, matrix,
    test = FALSE, verbose = FALSE, progress = FALSE), budget_s = light_budget)
  single_call <- capture_call(function() fastphylosig::fast_k(ctx, matrix[, 1L],
    test = FALSE, verbose = FALSE, progress = FALSE), budget_s = light_budget)
  matrix_ok <- matrix_call$status == "PASS" && single_call$status == "PASS" &&
    numeric_close(matrix_call$value$estimate[[1L]], extract_k(single_call$value))
  ok <- immutable && parity && replay && matrix_ok
  list(status = if (ok) "PASS" else "CORRECTNESS_FAIL",
       detail = paste0("parity=", parity, "; replay=", replay,
                       "; matrix=", matrix_ok, "; immutable=", immutable),
       id = paste(shape, n, sep = "_"))
}

reference_rows <- function(tree, cont, binary, shape, n, fixture_seed) {
  if (n > 150L) return(list())
  rows <- list()
  if (requireNamespace("phytools", quietly = TRUE)) for (method in c("K", "lambda")) {
    ref_call <- if (method == "K") function() phytools::phylosig(
      tree, cont, method = "K", test = FALSE) else function() phytools::phylosig(tree, cont, method = "lambda")
    ref <- capture_call(ref_call, budget_s = light_budget)
    fast <- if (method == "K") capture_call(function() fastphylosig::fast_k(
      tree, cont, test = FALSE, verbose = FALSE, progress = FALSE), budget_s = light_budget) else
      capture_call(function() fastphylosig::fast_lambda(tree, cont, verbose = FALSE,
                                                        progress = FALSE), budget_s = light_budget)
    ref_value <- if (ref$status == "PASS") if (method == "K") as.numeric(ref$value) else as.numeric(ref$value$lambda) else NA_real_
    fast_value <- if (fast$status == "PASS") if (method == "K") extract_k(fast$value) else extract_lambda(fast$value) else NA_real_
    guard <- ref$status == "PASS" && fast$status == "PASS" && numeric_close(ref_value, fast_value, 1e-5)
    ref_row <- row_template("reference", method, shape, n, 1L, 0L,
      "reference", 0L, ref, 1L, fixture_seed, paste(shape, n, method, sep = "_"),
      "phytools", if (guard) "PASS" else "CORRECTNESS_FAIL")
    ref_row$status <- if (guard) "PASS" else "CORRECTNESS_FAIL"
    rows[[length(rows) + 1L]] <- ref_row
  }
  if (requireNamespace("caper", quietly = TRUE)) {
    dat <- data.frame(species = names(binary), state = as.numeric(binary), stringsAsFactors = FALSE)
    # caper resolves these columns by numeric index on current releases;
    # using indices also avoids locale-dependent name matching in the guard.
    ref <- capture_call(function() caper::phylo.d(data = dat, phy = tree,
      names.col = 1L, binvar = 2L, permut = 0L), budget_s = light_budget)
    rows[[length(rows) + 1L]] <- row_template("reference", "D", shape, n, 1L, 0L,
      "reference", 0L, ref, 1L, fixture_seed, paste(shape, n, "D", sep = "_"),
      "caper", if (ref$status == "PASS") "PASS" else "UNAVAILABLE")
  }
  rows
}

run_pair <- function(tree, ctx, cont, binary, workload, method, shape, n, traits,
                     nsim, fixture_seed, fun_raw, fun_prepared, budget_s, guard,
                     routes_to_time = c("raw", "prepared")) {
  workload_for <- function(route) {
    if (length(workload) > 1L && !is.null(names(workload)))
      as.character(workload[[route]]) else as.character(workload[[1L]])
  }
  status_workload <- paste(unique(vapply(routes_to_time, workload_for,
                                         character(1L))), collapse = "+")
  if (!identical(guard$status, "PASS")) {
    record_status(paste(status_workload, method, shape, n, "guard", sep = ":"),
                  "CORRECTNESS_FAIL", guard$detail, status_workload)
    return(invisible(NULL))
  }
  warm_raw <- if ("raw" %in% routes_to_time)
    capture_call(fun_raw, seed = fixture_seed + 301L, budget_s = budget_s) else NULL
  warm_prepared <- if ("prepared" %in% routes_to_time)
    capture_call(fun_prepared, seed = fixture_seed + 302L, budget_s = budget_s) else NULL
  if ((!is.null(warm_raw) && warm_raw$status != "PASS") ||
      (!is.null(warm_prepared) && warm_prepared$status != "PASS")) {
    record_status(paste(status_workload, method, shape, n, "warmup", sep = ":"),
                  if (!is.null(warm_raw) && warm_raw$status != "PASS") warm_raw$status else warm_prepared$status,
                  "warmup failed", workload)
    return(invisible(NULL))
  }
  timing_statuses <- character()
  for (rep_id in seq_len(repeats)) {
    routes <- if (rep_id %% 2L) c("raw", "prepared") else c("prepared", "raw")
    routes <- routes[routes %in% routes_to_time]
    for (route in routes) {
      timed <- capture_call(if (route == "raw") fun_raw else fun_prepared,
                            seed = fixture_seed + 1000L + rep_id, budget_s = budget_s)
      timing_statuses <- c(timing_statuses, timed$status)
      append_rows(row_template(workload_for(route), method, shape, n, traits, nsim, route,
        rep_id, timed, 1L, fixture_seed, guard$id), timing_path)
    }
  }
  timing_status <- if (all(timing_statuses == "PASS")) "PASS" else
    if (any(timing_statuses == "ERROR")) "ERROR" else "RESOURCE_LIMIT"
  record_status(paste(status_workload, method, shape, n, "timing", sep = ":"),
                timing_status, paste0("repeats=", repeats), status_workload)
  invisible(NULL)
}

run_memory <- function() {
  helper <- file.path(args$repo, "benchmarks", "final_release", "benchmark", "memory_validate.R")
  if (!file.exists(helper)) helper <- file.path(getwd(), "memory_validate.R")
  if (!file.exists(helper)) {
    record_status("memory_helper", "ERROR", "memory_validate.R not found", "E-memory")
    return(invisible(NULL))
  }
  for (shape in memory_shapes) for (n in memory_sizes) {
    seed <- 842020L + n
    tree_call <- capture_call(function() make_tree(n, shape, seed = seed), budget_s = memory_budget)
    if (tree_call$status != "PASS") {
      row <- data.frame(timestamp_utc = provenance$timestamp_utc, commit = provenance$commit,
        package_version = provenance$package_version, r_version = provenance$r_version,
        os_type = provenance$os_type, compiler_cxx17 = provenance$compiler_cxx17,
        openmp = provenance$openmp, blas = provenance$blas, lapack = provenance$lapack,
        locale = provenance$locale, shape = shape, n = n,
        status = tree_call$status, object_size_bytes = NA_real_, rds_size_bytes = NA_real_,
        save_elapsed_s = NA_real_, fresh_read_validate_elapsed_s = NA_real_,
        rds_compression = "none", error_message = tree_call$error_message,
        stringsAsFactors = FALSE)
      append_rows(row, memory_path); next
    }
    prepared <- capture_call(function() fastphylosig::prepare_tree(tree_call$value), budget_s = memory_budget)
    if (prepared$status != "PASS") {
      row <- data.frame(timestamp_utc = provenance$timestamp_utc, commit = provenance$commit,
        package_version = provenance$package_version, r_version = provenance$r_version,
        os_type = provenance$os_type, compiler_cxx17 = provenance$compiler_cxx17,
        openmp = provenance$openmp, blas = provenance$blas, lapack = provenance$lapack,
        locale = provenance$locale, shape = shape, n = n,
        status = prepared$status, object_size_bytes = NA_real_, rds_size_bytes = NA_real_,
        save_elapsed_s = NA_real_, fresh_read_validate_elapsed_s = NA_real_,
        rds_compression = "none", error_message = prepared$error_message,
        stringsAsFactors = FALSE)
      append_rows(row, memory_path); record_status(paste("memory", shape, n, sep = ":"),
        prepared$status, prepared$error_message, "E-memory"); next
    }
    ctx <- prepared$value
    rds_path <- file.path(args$out, paste0(".memory-", shape, "-", n, ".rds"))
    fresh_path <- file.path(args$out, paste0(".memory-", shape, "-", n, ".csv"))
    save_call <- capture_call(function() saveRDS(ctx, rds_path, compress = FALSE), budget_s = memory_budget)
    rds_size <- if (file.exists(rds_path)) as.numeric(file.info(rds_path)$size) else NA_real_
    fresh_status <- "RESOURCE_LIMIT"; fresh_elapsed <- NA_real_; fresh_error <- "fresh process not started"
    if (save_call$status == "PASS" && is.finite(rds_size)) {
      rscript <- file.path(R.home("bin"), "Rscript.exe")
      fresh_call <- capture_call(function() {
        system2(rscript, c(helper, rds_path, fresh_path, args$repo), stdout = TRUE, stderr = TRUE)
        if (!file.exists(fresh_path)) stop("fresh validation did not write output", call. = FALSE)
        parsed <- utils::read.csv(fresh_path, stringsAsFactors = FALSE)
        if (!nrow(parsed) || parsed$status[[1L]] != "PASS")
          stop(if (nrow(parsed)) parsed$error_message[[1L]] else "fresh validation failed", call. = FALSE)
        as.numeric(parsed$elapsed_s[[1L]])
      }, budget_s = memory_budget)
      fresh_status <- fresh_call$status; fresh_elapsed <- if (fresh_call$status == "PASS") fresh_call$value else NA_real_
      fresh_error <- if (fresh_call$status == "PASS") "" else fresh_call$error_message
    }
    final_status <- if (save_call$status == "PASS" && fresh_status == "PASS") "PASS" else
      if (save_call$status == "RESOURCE_LIMIT" || fresh_status == "RESOURCE_LIMIT") "RESOURCE_LIMIT" else "ERROR"
    row <- data.frame(timestamp_utc = provenance$timestamp_utc, commit = provenance$commit,
      package_version = provenance$package_version, r_version = provenance$r_version,
      os_type = provenance$os_type, compiler_cxx17 = provenance$compiler_cxx17,
      openmp = provenance$openmp, blas = provenance$blas, lapack = provenance$lapack,
      locale = provenance$locale, shape = shape, n = n,
      status = final_status, object_size_bytes = as.numeric(object.size(ctx)),
      rds_size_bytes = rds_size, save_elapsed_s = save_call$elapsed_s,
      fresh_read_validate_elapsed_s = fresh_elapsed, rds_compression = "none",
      error_message = paste(c(save_call$error_message, fresh_error), collapse = " | "),
      stringsAsFactors = FALSE)
    append_rows(row, memory_path)
    record_status(paste("memory", shape, n, sep = ":"), final_status,
      paste0("object.size=", row$object_size_bytes, "; rds=", row$rds_size_bytes), "E-memory")
    if (file.exists(rds_path)) file.remove(rds_path)
    if (file.exists(fresh_path)) file.remove(fresh_path)
  }
  invisible(NULL)
}

if (!memory_only) {
  # Keep a small, fixed external-reference guard in formal runs without
  # changing the user-workload timing grid.  Reference timing is not used
  # for speed claims; it only checks the final installed source.
  if (formal) for (shape in benchmark_shapes) {
    reference_n <- 100L
    reference_seed <- 492020L + match(shape, benchmark_shapes) * 1000L
    reference_tree <- make_tree(reference_n, shape, seed = reference_seed)
    reference_cont <- make_continuous(reference_tree, 1L, reference_seed + 1L)
    reference_binary <- make_binary(reference_tree)
    refs <- reference_rows(reference_tree, reference_cont, reference_binary,
                           shape, reference_n, reference_seed)
    if (length(refs)) {
      ref_table <- do.call(rbind, refs)
      append_rows(ref_table, timing_path)
      ref_status <- if (any(ref_table$reference_status == "CORRECTNESS_FAIL"))
        "CORRECTNESS_FAIL" else "PASS"
      record_status(paste("reference_guard", shape, reference_n, sep = ":"),
                    ref_status, paste0("rows=", nrow(ref_table)), "reference")
    } else {
      record_status(paste("reference_guard", shape, reference_n, sep = ":"),
                    "UNAVAILABLE", "no optional reference package available", "reference")
    }
  }
  for (shape in benchmark_shapes) for (n in sizes) {
    seed <- 562020L + n + match(shape, benchmark_shapes) * 100000L
    tree <- make_tree(n, shape, seed); cont <- make_continuous(tree, 1L, seed + 1L)
    binary <- make_binary(tree); ctx <- prepare_tree(tree)
    guard <- run_guard(tree, ctx, cont, binary, shape, n, seed)
    record_status(paste("guard", shape, n, sep = ":"), guard$status, guard$detail, "guard")
    refs <- reference_rows(tree, cont, binary, shape, n, seed)
    if (length(refs)) append_rows(do.call(rbind, refs), timing_path)
    run_pair(tree, ctx, cont, binary, c(raw = "A-single-trait", prepared = "B-prepared-repeated"), "K", shape, n, 1L, 0L, seed,
      function() fastphylosig::fast_k(tree, cont, test = FALSE, verbose = FALSE, progress = FALSE),
      function() fastphylosig::fast_k(ctx, cont, test = FALSE, verbose = FALSE, progress = FALSE), light_budget, guard)
    run_pair(tree, ctx, cont, binary, c(raw = "A-single-trait", prepared = "B-prepared-repeated"), "lambda", shape, n, 1L, 0L, seed + 10L,
      function() fastphylosig::fast_lambda(tree, cont, verbose = FALSE, progress = FALSE),
      function() fastphylosig::fast_lambda(ctx, cont, verbose = FALSE, progress = FALSE), light_budget, guard)
    run_pair(tree, ctx, cont, binary, c(raw = "A-single-trait", prepared = "B-prepared-repeated"), "fast_signal_K", shape, n, 1L, 0L, seed + 20L,
      function() fastphylosig::fast_signal(tree, cont, method = "K", test = FALSE, verbose = FALSE, progress = FALSE),
      function() fastphylosig::fast_signal(ctx, cont, method = "K", test = FALSE, verbose = FALSE, progress = FALSE), light_budget, guard)
  }
  matrix_n <- if (formal) 500L else min(sizes)
  for (shape in benchmark_shapes) {
    tree <- make_tree(matrix_n, shape, seed = 982020L + match(shape, benchmark_shapes)); ctx <- prepare_tree(tree)
    for (traits in matrix_traits) {
      seed <- 682020L + traits + match(shape, benchmark_shapes) * 1000L; matrix <- make_continuous(tree, traits, seed)
      cont <- matrix[, 1L]; binary <- make_binary(tree); guard <- run_guard(tree, ctx, cont, binary, shape, matrix_n, seed)
      matrix_raw_guard <- capture_call(function() fastphylosig::fast_k(tree, matrix,
        test = FALSE, verbose = FALSE, progress = FALSE), budget_s = light_budget)
      matrix_prepared_guard <- capture_call(function() fastphylosig::fast_k(ctx, matrix,
        test = FALSE, verbose = FALSE, progress = FALSE), budget_s = light_budget)
      matrix_parity <- matrix_raw_guard$status == "PASS" &&
        matrix_prepared_guard$status == "PASS" &&
        identical(matrix_raw_guard$value$status, matrix_prepared_guard$value$status) &&
        numeric_close(matrix_raw_guard$value$estimate, matrix_prepared_guard$value$estimate)
      if (!matrix_parity) {
        guard$status <- "CORRECTNESS_FAIL"
        guard$detail <- paste0(guard$detail, "; matrix_raw_prepared=FALSE")
      }
      run_pair(tree, ctx, cont, binary, "C-batch-matrix", "K", shape, matrix_n, traits, 0L, seed,
        function() fastphylosig::fast_k(tree, matrix, test = FALSE, verbose = FALSE, progress = FALSE),
        function() fastphylosig::fast_k(ctx, matrix, test = FALSE, verbose = FALSE, progress = FALSE), light_budget, guard)
      run_pair(tree, ctx, cont, binary, "C-batch-matrix", "fast_signal_K", shape, matrix_n, traits, 0L, seed + 1L,
        function() fastphylosig::fast_signal(tree, matrix, method = "K", test = FALSE, verbose = FALSE, progress = FALSE),
        function() fastphylosig::fast_signal(ctx, matrix, method = "K", test = FALSE, verbose = FALSE, progress = FALSE), light_budget, guard)
    }
  }
  for (shape in benchmark_shapes) for (n in permutation_sizes) {
    seed <- 772020L + n + match(shape, benchmark_shapes) * 100000L; tree <- make_tree(n, shape, seed); ctx <- prepare_tree(tree)
    cont <- make_continuous(tree, 1L, seed + 1L); binary <- make_binary(tree); guard <- run_guard(tree, ctx, cont, binary, shape, n, seed)
    run_pair(tree, ctx, cont, binary, "D-permutation-heavy", "K", shape, n, 1L, permutation_nsim, seed,
      function() fastphylosig::fast_k(tree, cont, test = TRUE, nsim = permutation_nsim, verbose = FALSE, progress = FALSE),
      function() fastphylosig::fast_k(ctx, cont, test = TRUE, nsim = permutation_nsim, verbose = FALSE, progress = FALSE), heavy_budget, guard)
    run_pair(tree, ctx, cont, binary, "D-permutation-heavy", "D", shape, n, 1L, permutation_nsim, seed + 10L,
      function() fastphylosig::fast_d(tree, binary, test = TRUE, nsim = permutation_nsim, return_sim = FALSE, keep_null = FALSE, ncores = 1L, verbose = FALSE, progress = FALSE),
      function() fastphylosig::fast_d(ctx, binary, test = TRUE, nsim = permutation_nsim, return_sim = FALSE, keep_null = FALSE, ncores = 1L, verbose = FALSE, progress = FALSE), heavy_budget, guard)
  }
}

run_memory()

timing <- if (file.exists(timing_path)) utils::read.csv(timing_path, stringsAsFactors = FALSE) else data.frame()
memory <- if (file.exists(memory_path)) utils::read.csv(memory_path, stringsAsFactors = FALSE) else data.frame()
status_rows <- if (file.exists(status_path)) utils::read.csv(status_path, stringsAsFactors = FALSE) else data.frame()
timing_failures <- if (nrow(timing)) sum(timing$workload != "reference" & timing$status != "PASS") else 0L
memory_failures <- if (nrow(memory)) sum(memory$status != "PASS") else 0L
guard_failures <- if (nrow(status_rows))
  sum(grepl("^guard:", status_rows$check) & status_rows$status != "PASS") else 0L
reference_failures <- if (nrow(status_rows))
  sum(status_rows$workload == "reference" & status_rows$status == "CORRECTNESS_FAIL") else 0L
complete_status <- if (formal && (timing_failures || memory_failures || guard_failures || reference_failures))
  "INCOMPLETE" else "PASS"
record_status("formal_completeness", complete_status,
              paste0("timing_failures=", timing_failures,
                     "; memory_failures=", memory_failures,
                     "; guard_failures=", guard_failures,
                     "; reference_failures=", reference_failures), "summary")

memory_report <- file.path(dirname(memory_path), "MEMORY_CHARACTERIZATION.md")
memory_lines <- c(
  "# Memory Characterization",
  "",
  paste0("Status: ", if (nrow(memory) && all(memory$status == "PASS")) "PASS" else "INCOMPLETE"),
  "",
  "Single-machine measurements from the exact source commit below. `object.size` is an R object size and is not peak RSS. RDS files were written with `compress = FALSE`; temporary RDS and fresh-process CSV files were removed after validation.",
  "",
  paste0("- Commit: `", provenance$commit, "`"),
  paste0("- Package version: `", provenance$package_version, "`"),
  paste0("- R: `", provenance$r_version, "`"),
  paste0("- OS: `", provenance$os_type, " ", provenance$os_release, "`"),
  paste0("- Compiler: `", provenance$compiler_cxx17, "`"),
  paste0("- Locale: `", provenance$locale, "`"),
  "",
  "| shape | n | status | object.size (bytes) | RDS (bytes) | save (s) | fresh read + first validation (s) |",
  "|---|---:|---|---:|---:|---:|---:|"
)
if (nrow(memory)) for (i in seq_len(nrow(memory))) {
  memory_lines <- c(memory_lines, paste0(
    "| ", memory$shape[[i]], " | ", memory$n[[i]], " | ", memory$status[[i]],
    " | ", memory$object_size_bytes[[i]], " | ", memory$rds_size_bytes[[i]],
    " | ", formatC(memory$save_elapsed_s[[i]], format = "fg", digits = 6),
    " | ", formatC(memory$fresh_read_validate_elapsed_s[[i]], format = "fg", digits = 6), " |"
  ))
}
if (!nrow(memory)) memory_lines <- c(memory_lines, "| no completed rows | | | | | | |")
writeLines(memory_lines, memory_report)

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (requireNamespace("digest", quietly = TRUE)) digest::digest(file = path, algo = "sha256") else NA_character_
}
writeLines(c(paste0("FINAL_BENCHMARK_COMMIT=", provenance$commit),
  paste0("PACKAGE_VERSION=", provenance$package_version), paste0("MODE=", args$mode),
  paste0("TIMING_ROWS=", nrow(timing)), paste0("MEMORY_ROWS=", nrow(memory)),
  paste0("TIMING_SHA256=", sha256_file(timing_path)), paste0("MEMORY_SHA256=", sha256_file(memory_path)),
  paste0("STATUS_SHA256=", sha256_file(status_path)), paste0("FORMAL_COMPLETENESS=", complete_status)),
  file.path(args$out, "benchmark_manifest.txt"))
