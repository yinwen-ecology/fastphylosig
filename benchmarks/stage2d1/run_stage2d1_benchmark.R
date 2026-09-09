#!/usr/bin/env Rscript

# Stage 2D1 performance benchmark for the private K permutation prototype.
# This script is deliberately opt-in and refuses to run unless the preceding
# correctness script recorded PASS in the same output directory.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- all_args[grepl("^--file=", all_args)]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else
  file.path(getwd(), "run_stage2d1_benchmark.R")
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
source(file.path(script_dir, "stage2d1_common.R"), local = TRUE)

args <- stage2d1_args(commandArgs(trailingOnly = TRUE))
if (!isTRUE(args$smoke) && !isTRUE(args$formal)) {
  stop("Formal benchmark is explicit: add --formal or --smoke.", call. = FALSE)
}
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(args$out, "stage2d1_correctness_status.csv")
if (!file.exists(status_path)) {
  stop("Correctness status is missing; run run_stage2d1_correctness.R first.",
       call. = FALSE)
}
status <- utils::read.csv(status_path, stringsAsFactors = FALSE)
if (!nrow(status) || !identical(as.character(status$status[[1L]]), "PASS")) {
  stop("Correctness/RNG gates are not PASS; formal timing is blocked.",
       call. = FALSE)
}

Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")
assign("stage2d1_smoke", isTRUE(args$smoke), envir = .GlobalEnv)
Sys.setenv(FASTPHYLOSIG_STAGE2D1_RETURN_ORDERED = "false")
ns <- stage2d1_load_package(args$repo)
hooks <- stage2d1_load_hooks(args$repo)
if (!is.function(hooks$rng)) {
  stop("The prototype RNG hook is unavailable after correctness PASS.", call. = FALSE)
}
production_engine <- stage2d1_production_cpp(ns)
baseline_kind <- stage2d1_env(
  "FASTPHYLOSIG_STAGE2D1_BENCH_BASELINE", "production"
)
if (!baseline_kind %in% c("production", "private_oracle")) {
  stop("FASTPHYLOSIG_STAGE2D1_BENCH_BASELINE must be production or private_oracle.",
       call. = FALSE)
}
old_engine <- production_engine
if (identical(baseline_kind, "private_oracle")) {
  if (!is.function(hooks$prototype)) {
    stop("The private oracle requires the checked-in prototype entry point.",
         call. = FALSE)
  }
  oracle_target <- hooks$prototype
  old_engine <- function(...) {
    call_args <- list(...)
    call_args$mode <- "oracle"
    call_args$return_ordered <- FALSE
    call_args$collect_counters <- FALSE
    do.call(oracle_target, call_args)
  }
}

n_grid <- stage2d1_grid("FASTPHYLOSIG_STAGE2D1_BENCH_N",
                        c(500L, 2000L, 5000L, 10000L), c(500L), minimum = 2L)
nsim_grid <- stage2d1_grid("FASTPHYLOSIG_STAGE2D1_BENCH_NSIM",
                           c(199L, 999L, 9999L), c(19L), minimum = 1L)
ncores <- stage2d1_grid("FASTPHYLOSIG_STAGE2D1_BENCH_NCORES",
                        c(1L, 2L), c(1L), minimum = 1L)
if (isTRUE(args$smoke)) ncores <- 1L
detected_cores <- tryCatch(parallel::detectCores(logical = TRUE),
                           error = function(e) 1L)
if (is.finite(detected_cores) && detected_cores >= 1L) {
  ncores <- ncores[ncores <= detected_cores]
}
if (!length(ncores)) ncores <- 1L
shape_value <- stage2d1_env("FASTPHYLOSIG_STAGE2D1_BENCH_SHAPES")
shapes <- if (nzchar(shape_value)) {
  unique(trimws(strsplit(shape_value, ",", fixed = TRUE)[[1L]]))
} else "balanced"
if (any(!shapes %in% c("balanced", "random", "pectinate"))) {
  stop("FASTPHYLOSIG_STAGE2D1_BENCH_SHAPES contains an unknown shape.",
       call. = FALSE)
}
trait_chunk <- as.integer(stage2d1_env("FASTPHYLOSIG_STAGE2D1_TRAIT_CHUNK", "64"))
if (is.na(trait_chunk) || trait_chunk < 1L) stop("invalid trait chunk", call. = FALSE)
repeats_value <- suppressWarnings(as.integer(stage2d1_env(
  "FASTPHYLOSIG_STAGE2D1_BENCH_REPEATS", if (isTRUE(args$smoke)) "1" else "10"
)))
if (is.na(repeats_value) || repeats_value < 1L) stop("invalid repeats", call. = FALSE)
seed <- as.integer(stage2d1_env("FASTPHYLOSIG_STAGE2D1_BENCH_SEED", "20260902"))
if (is.na(seed)) stop("invalid benchmark seed", call. = FALSE)

timing_rows <- list()
raw_rows <- list()
batch_rows <- list()
memory_rows <- list()
add_timing <- function(x) timing_rows[[length(timing_rows) + 1L]] <<- x

run_direct <- function(fun, mode, compiled_tree, X, nsim, threads, chunk,
                       return_sim = FALSE, seed_value = seed) {
  stage2d1_capture(function() stage2d1_call_engine(
    fun, mode = mode, compiled_tree = compiled_tree, X = X, nsim = nsim,
    permutations = NULL, trait_chunk = trait_chunk, return_sim = return_sim,
    include_observed = TRUE, ncores = threads, simulation_chunk = chunk,
    seed = seed_value
  ))
}

time_pair <- function(shape, n, nsim, threads, chunk, pair) {
  tree <- stage2d1_make_tree(n, shape)
  compiled_tree <- stage2d1_compile_tree(tree, ns)
  X <- stage2d1_make_matrix(tree, 1L)
  warm_old <- run_direct(old_engine, "rng", compiled_tree, X, nsim, threads,
                         chunk, return_sim = FALSE)
  warm_new <- run_direct(hooks$rng, "rng", compiled_tree, X, nsim, threads,
                         chunk, return_sim = FALSE)
  if (warm_old$status != "ok" || warm_new$status != "ok") {
    stop("Warmup failed for ", shape, "/n=", n, "/nsim=", nsim,
         "/ncores=", threads, call. = FALSE)
  }
  for (i in seq_len(repeats_value)) {
    # Alternate first caller on every paired repeat.  Each call receives the
    # same seed, so timing order cannot change the stochastic input.
    first_old <- ((i + pair) %% 2L) == 0L
    one <- if (first_old) old_engine else hooks$rng
    two <- if (first_old) hooks$rng else old_engine
    one_name <- if (first_old) "old" else "prototype"
    two_name <- if (first_old) "prototype" else "old"
    first <- run_direct(one, "rng", compiled_tree, X, nsim, threads, chunk,
                        return_sim = FALSE)
    second <- run_direct(two, "rng", compiled_tree, X, nsim, threads, chunk,
                         return_sim = FALSE)
    add_timing(data.frame(
      route = "prepared_direct", implementation = one_name, shape = shape,
      n = n, nsim = nsim, ncores = threads, simulation_chunk = chunk,
      pair = i, order = paste(one_name, two_name, sep = ">"),
      elapsed_s = first$elapsed_s, status = first$status,
      warning_text = first$warning_text, error_message = first$error_message,
      stringsAsFactors = FALSE
    ))
    add_timing(data.frame(
      route = "prepared_direct", implementation = two_name, shape = shape,
      n = n, nsim = nsim, ncores = threads, simulation_chunk = chunk,
      pair = i, order = paste(one_name, two_name, sep = ">"),
      elapsed_s = second$elapsed_s, status = second$status,
      warning_text = second$warning_text, error_message = second$error_message,
      stringsAsFactors = FALSE
    ))
  }
  old_proxy <- stage2d1_workspace_proxy(n, nsim, 1L, trait_chunk, chunk,
                                        return_sim = FALSE, candidate = FALSE)
  new_proxy <- stage2d1_workspace_proxy(n, nsim, 1L, trait_chunk, chunk,
                                        return_sim = FALSE, candidate = TRUE)
  memory_rows[[length(memory_rows) + 1L]] <<- cbind(
    old_proxy, implementation = "old", shape = shape
  )
  memory_rows[[length(memory_rows) + 1L]] <<- cbind(
    new_proxy, implementation = "prototype", shape = shape
  )
  invisible(NULL)
}

if (isTRUE(args$smoke)) {
  for (shape in shapes) for (n in n_grid) for (nsim in nsim_grid)
    for (threads in ncores) time_pair(shape, n, nsim, threads, 128L, 1L)
} else {
  # The formal benchmark is serialized: no cells are run concurrently.
  pair <- 0L
  for (shape in shapes) for (n in n_grid) for (nsim in nsim_grid)
    for (threads in ncores) {
      pair <- pair + 1L
      message("[stage2d1] prepared prototype ", shape, "/n=", n,
              "/nsim=", nsim, "/ncores=", threads)
      time_pair(shape, n, nsim, threads, 128L, pair)
    }
}

# Public raw route is retained as an end-to-end guard.  It is not used as a
# candidate speedup because no raw production binding is changed in Stage 2D1.
raw_repeats <- if (isTRUE(args$smoke)) 1L else min(3L, repeats_value)
for (shape in shapes) for (n in n_grid) for (nsim in nsim_grid) {
  tree <- stage2d1_make_tree(n, shape)
  X <- stage2d1_make_trait(tree)
  for (threads in ncores) {
    chunk <- 128L
    warm <- stage2d1_capture(function() fastphylosig::fast_k(
      tree, X, test = TRUE, nsim = nsim, ncores = threads,
      simulation_chunk = chunk, return_sim = FALSE, verbose = FALSE,
      progress = FALSE
    ))
    if (warm$status != "ok") stop("raw warmup failed", call. = FALSE)
    for (i in seq_len(raw_repeats)) {
      measured <- stage2d1_capture(function() stage2d1_restore_seed(
        seed, function() fastphylosig::fast_k(
          tree, X, test = TRUE, nsim = nsim, ncores = threads,
          simulation_chunk = chunk, return_sim = FALSE, verbose = FALSE,
          progress = FALSE
        )
      ))
      raw_rows[[length(raw_rows) + 1L]] <- data.frame(
        route = "raw_public", implementation = "production", shape = shape,
        n = n, nsim = nsim, ncores = threads, simulation_chunk = chunk,
        pair = i, elapsed_s = measured$elapsed_s, status = measured$status,
        warning_text = measured$warning_text, error_message = measured$error_message,
        stringsAsFactors = FALSE
      )
    }
  }
}

# test=FALSE matrix protection at n=500.  The candidate is intentionally not
# substituted into this path; this records the unchanged production baseline.
batch_traits <- stage2d1_grid("FASTPHYLOSIG_STAGE2D1_TRAITS",
                              c(1L, 8L, 32L, 100L), c(1L, 8L), minimum = 1L)
batch_n <- if (isTRUE(args$smoke)) 50L else 500L
for (p in batch_traits) {
  tree <- stage2d1_make_tree(batch_n, "balanced")
  X <- stage2d1_make_matrix(tree, p)
  warm <- stage2d1_capture(function() fastphylosig::fast_k(
    tree, X, test = FALSE, verbose = FALSE, progress = FALSE
  ))
  for (i in seq_len(if (isTRUE(args$smoke)) 1L else min(5L, repeats_value))) {
    measured <- stage2d1_capture(function() fastphylosig::fast_k(
      tree, X, test = FALSE, verbose = FALSE, progress = FALSE
    ))
    batch_rows[[length(batch_rows) + 1L]] <- data.frame(
      gate = "test_false", implementation = "production", traits = p,
      n = batch_n, nsim = NA_integer_, pair = i,
      elapsed_s = measured$elapsed_s, status = measured$status,
      warning_text = measured$warning_text, comparison = "unchanged_path",
      stringsAsFactors = FALSE
    )
  }
}

timings <- if (length(timing_rows)) do.call(rbind, timing_rows) else data.frame()
raw <- if (length(raw_rows)) do.call(rbind, raw_rows) else data.frame()
batch <- if (length(batch_rows)) do.call(rbind, batch_rows) else data.frame()
memory <- if (length(memory_rows)) do.call(rbind, memory_rows) else data.frame()
stage2d1_write_csv(timings, file.path(args$out, "stage2d1_prepared_timings.csv"))
stage2d1_write_csv(raw, file.path(args$out, "stage2d1_raw_guard_timings.csv"))
stage2d1_write_csv(batch, file.path(args$out, "stage2d1_test_false_batch_guard.csv"))
stage2d1_write_csv(memory, file.path(args$out, "stage2d1_memory_proxy.csv"))

summarize <- function(dat, keys) {
  if (!nrow(dat)) return(data.frame())
  key <- do.call(paste, c(dat[keys], sep = "\r"))
  groups <- split(seq_len(nrow(dat)), key, drop = TRUE)
  out <- lapply(groups, function(ii) {
    z <- dat[ii, , drop = FALSE]
    values <- z$elapsed_s[z$status == "ok" & is.finite(z$elapsed_s)]
    if (!length(values)) values <- NA_real_
    z[1L, keys, drop = FALSE]
    data.frame(
      z[1L, keys, drop = FALSE], observations = nrow(z),
      ok_observations = sum(z$status == "ok"),
      median_s = if (all(is.na(values))) NA_real_ else median(values, na.rm = TRUE),
      IQR_s = if (all(is.na(values))) NA_real_ else IQR(values, na.rm = TRUE),
      min_s = if (all(is.na(values))) NA_real_ else min(values, na.rm = TRUE),
      max_s = if (all(is.na(values))) NA_real_ else max(values, na.rm = TRUE),
      warnings = paste(unique(z$warning_text[nzchar(z$warning_text)]), collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

summary <- summarize(timings,
                     c("route", "implementation", "shape", "n", "nsim",
                       "ncores", "simulation_chunk"))
if (nrow(summary)) {
  old <- summary[summary$implementation == "old", , drop = FALSE]
  new <- summary[summary$implementation == "prototype", , drop = FALSE]
  join_keys <- c("route", "shape", "n", "nsim", "ncores", "simulation_chunk")
  names(old)[names(old) == "median_s"] <- "old_median_s"
  names(old)[names(old) == "IQR_s"] <- "old_IQR_s"
  names(new)[names(new) == "median_s"] <- "prototype_median_s"
  names(new)[names(new) == "IQR_s"] <- "prototype_IQR_s"
  old <- old[c(join_keys, "old_median_s", "old_IQR_s")]
  new <- new[c(join_keys, "prototype_median_s", "prototype_IQR_s")]
  paired <- merge(old, new, by = join_keys, all = TRUE)
  paired$reduction_fraction <- 1 - paired$prototype_median_s / paired$old_median_s
  paired$speedup <- paired$old_median_s / paired$prototype_median_s
  stage2d1_write_csv(paired, file.path(args$out, "stage2d1_prepared_summary.csv"))
} else {
  paired <- data.frame()
  stage2d1_write_csv(paired, file.path(args$out, "stage2d1_prepared_summary.csv"))
}

# Resource and acceptance gates are intentionally conservative.  Missing
# comparable candidate rows are NOT converted to a performance PASS.
heavy <- if (nrow(paired)) paired[paired$nsim >= 999L, , drop = FALSE] else data.frame()
heavy <- heavy[is.finite(heavy$reduction_fraction) & is.finite(heavy$speedup), , drop = FALSE]
heavy_pass <- nrow(heavy) >= 2L && sum(heavy$reduction_fraction >= 0.15) >= 2L &&
  any(heavy$speedup >= 1.20)
light <- if (nrow(paired)) paired[paired$nsim < 999L, , drop = FALSE] else data.frame()
light_slowdown <- if (nrow(light)) any(light$reduction_fraction < -0.05, na.rm = TRUE) else FALSE
batch_ok <- nrow(batch) > 0L && all(batch$status == "ok") &&
  all(!nzchar(batch$warning_text))
overall <- heavy_pass && !light_slowdown && batch_ok && nrow(memory) > 0L

stage2d1_write_csv(data.frame(
  candidate = "prototype_hook",
  benchmark_baseline = baseline_kind,
  controlled_exact = "PASS (prerequisite status file)",
  same_seed_replay = "PASS (prerequisite status file)",
  heavy_performance_gate = if (heavy_pass) "PASS" else "FAIL",
  light_slowdown_gate = if (!light_slowdown) "PASS" else "FAIL",
  test_false_batch_guard = if (batch_ok) "PASS" else "FAIL",
  memory_bounded_proxy = if (nrow(memory)) "RECORDED" else "FAIL",
  K_PERMUTATION_PROTOTYPE = if (overall) "PERFORMANCE_PASS" else "FAIL",
  production_integration = "NOT_RUN",
  production_code_changed = "NO",
  stringsAsFactors = FALSE
), file.path(args$out, "stage2d1_performance_status.csv"))

stage2d1_write_csv(stage2d1_provenance(args$repo, "Stage 2D1 benchmark", args, hooks),
                   file.path(args$out, "stage2d1_benchmark_provenance.csv"))
message("[stage2d1] benchmark complete; production integration remains blocked")
