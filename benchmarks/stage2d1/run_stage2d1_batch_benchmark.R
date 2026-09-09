#!/usr/bin/env Rscript

# Stage 2D1 private batch test=TRUE timing guard.
#
# This audit compares the checked-in prototype's exact private oracle with
# modes A and B.  It is deliberately separate from the broad performance
# benchmark: the grid is fixed to a balanced n=500 tree and p=1/8/32/100.
# No package source, tests, or public API are modified by this script.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- all_args[grepl("^--file=", all_args)]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else
  file.path(getwd(), "run_stage2d1_batch_benchmark.R")
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
source(file.path(script_dir, "stage2d1_common.R"), local = TRUE)

args <- stage2d1_args(commandArgs(trailingOnly = TRUE))
if (!isTRUE(args$formal) && !isTRUE(args$smoke)) {
  stop("Batch timing is explicit: add --formal or --smoke.", call. = FALSE)
}
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)

# The formal batch guard is downstream of the exact correctness/RNG gate.
if (isTRUE(args$formal)) {
  correctness_path <- file.path(args$out, "stage2d1_correctness_status.csv")
  if (!file.exists(correctness_path)) {
    stop("Correctness status is missing; run run_stage2d1_correctness.R first.",
         call. = FALSE)
  }
  correctness <- utils::read.csv(correctness_path, stringsAsFactors = FALSE)
  if (!nrow(correctness) || !identical(as.character(correctness$status[[1L]]),
                                       "PASS")) {
    stop("Correctness/RNG gates are not PASS; formal batch timing is blocked.",
         call. = FALSE)
  }
}

Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")
ns <- stage2d1_load_package(args$repo)

# Always use the checked-in source by default.  An explicit path is allowed
# for a copied, byte-identical audit source on machines with path constraints,
# but a missing source fails closed rather than silently using another hook.
proto_cpp <- stage2d1_env(
  "FASTPHYLOSIG_STAGE2D1_PROTO_CPP",
  file.path(args$repo, "benchmarks", "stage2d1", "prototype",
            "k_permutation_prototype.cpp")
)
proto_exists <- tryCatch(isTRUE(file.exists(proto_cpp)), error = function(e) FALSE)
if (!proto_exists) {
  stop("Checked-in Stage 2D1 prototype source is unavailable: ", proto_cpp,
       call. = FALSE)
}
if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Rcpp is required to compile the checked-in Stage 2D1 prototype.",
       call. = FALSE)
}
prototype_env <- new.env(parent = globalenv())
Rcpp::sourceCpp(normalizePath(proto_cpp, winslash = "/", mustWork = TRUE),
                env = prototype_env)
prototype <- if (exists("stage2d1_k_permutation_prototype",
                         envir = prototype_env, inherits = FALSE)) {
  get("stage2d1_k_permutation_prototype", envir = prototype_env,
      inherits = FALSE)
} else NULL
if (!is.function(prototype)) {
  stop("The checked-in prototype did not expose stage2d1_k_permutation_prototype().",
       call. = FALSE)
}

hooks_for_provenance <- list(
  source_file = "",
  source_cpp = proto_cpp,
  controlled_name = "stage2d1_k_permutation_prototype(mode=controlled)",
  rng_name = "stage2d1_k_permutation_prototype(mode=oracle/A/B)",
  batch_name = "",
  prototype_mode = "oracle,A,B"
)
stage2d1_frozen_contract(args$repo)

# The formal grid is intentionally non-configurable so that these rows remain
# comparable across machines and cannot be mistaken for the broad benchmark.
if (isTRUE(args$smoke)) {
  n <- 50L
  nsim <- 19L
  trait_grid <- c(1L, 8L)
  repeats <- 2L
} else {
  n <- 500L
  nsim <- 199L
  trait_grid <- c(1L, 8L, 32L, 100L)
  repeats <- 10L
}
ncores <- 1L
trait_chunk <- 64L
simulation_chunk <- 128L
base_seed <- as.integer(stage2d1_env(
  "FASTPHYLOSIG_STAGE2D1_BATCH_SEED", "20260909"
))
if (is.na(base_seed)) stop("invalid FASTPHYLOSIG_STAGE2D1_BATCH_SEED", call. = FALSE)

call_proto <- function(compiled_tree, X, nsim, mode, seed,
                       return_sim = FALSE, return_ordered = FALSE,
                       collect_counters = FALSE) {
  invoke <- function() prototype(
    compiled_tree = compiled_tree,
    X = X,
    nsim = as.integer(nsim),
    permutations = NULL,
    trait_chunk = as.integer(trait_chunk),
    return_sim = isTRUE(return_sim),
    include_observed = TRUE,
    n_threads = as.integer(ncores),
    simulation_chunk = as.integer(simulation_chunk),
    mode = mode,
    return_ordered = isTRUE(return_ordered),
    collect_counters = isTRUE(collect_counters)
  )
  stage2d1_restore_seed(seed, invoke)
}

capture_proto <- function(compiled_tree, X, nsim, mode, seed,
                          return_sim = FALSE, return_ordered = FALSE,
                          collect_counters = FALSE) {
  stage2d1_capture(function() call_proto(
    compiled_tree = compiled_tree, X = X, nsim = nsim, mode = mode,
    seed = seed, return_sim = return_sim, return_ordered = return_ordered,
    collect_counters = collect_counters
  ))
}

timing_rows <- list()
parity_rows <- list()
memory_rows <- list()

append_timing <- function(implementation, mode, candidate_group, p, pair,
                          seed, order, captured) {
  timing_rows[[length(timing_rows) + 1L]] <<- data.frame(
    stage = "Stage 2D1 batch test=TRUE guard",
    baseline = "private_oracle",
    implementation = implementation,
    mode = mode,
    candidate_group = candidate_group,
    shape = "balanced",
    n = n,
    nsim = nsim,
    traits = p,
    ncores = ncores,
    trait_chunk = trait_chunk,
    simulation_chunk = simulation_chunk,
    pair = pair,
    seed = seed,
    order = order,
    elapsed_s = captured$elapsed_s,
    status = captured$status,
    warning_text = captured$warning_text,
    error_message = captured$error_message,
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

for (p in trait_grid) {
  message("[stage2d1] private batch guard p=", p,
          ", n=", n, ", nsim=", nsim)
  tree <- stage2d1_make_tree(n, "balanced")
  compiled_tree <- stage2d1_compile_tree(tree, ns)
  X <- stage2d1_make_matrix(tree, p)

  # Exact parity is checked with the ordered null retained.  This is separate
  # from timing so output construction cannot obscure the test=TRUE timings.
  parity_seed <- base_seed + as.integer(p)
  oracle_exact <- capture_proto(
    compiled_tree, X, nsim, "oracle", parity_seed,
    return_sim = TRUE, return_ordered = TRUE
  )
  for (candidate_mode in c("A", "B")) {
    candidate_exact <- capture_proto(
      compiled_tree, X, nsim, candidate_mode, parity_seed,
      return_sim = TRUE, return_ordered = TRUE
    )
    compared <- stage2d1_compare_captures(
      oracle_exact, candidate_exact, nsim, p, numeric_tolerance = 0
    )
    exact_pass <- isTRUE(compared$pass) &&
      identical(oracle_exact$status, "ok") &&
      identical(candidate_exact$status, "ok")
    max_abs <- compared$result$max_abs
    parity_rows[[length(parity_rows) + 1L]] <- data.frame(
      stage = "Stage 2D1 batch test=TRUE guard",
      baseline = "private_oracle",
      candidate = candidate_mode,
      shape = "balanced",
      n = n,
      nsim = nsim,
      traits = p,
      ncores = ncores,
      seed = parity_seed,
      oracle_status = oracle_exact$status,
      candidate_status = candidate_exact$status,
      oracle_warning = oracle_exact$warning_text,
      candidate_warning = candidate_exact$warning_text,
      status_pass = isTRUE(compared$status_pass),
      warning_pass = isTRUE(compared$warning_pass),
      result_pass = isTRUE(compared$result$pass),
      pass = exact_pass,
      max_abs_K = unname(max_abs[["K"]]),
      max_abs_P = unname(max_abs[["P"]]),
      max_abs_MCSE_P = unname(max_abs[["MCSE_P"]]),
      max_abs_exceedance = unname(max_abs[["exceedance"]]),
      max_abs_sim_K = unname(max_abs[["sim_K"]]),
      oracle_error = oracle_exact$error_message,
      candidate_error = candidate_exact$error_message,
      stringsAsFactors = FALSE
    )
  }

  # Warm both implementations before the paired series.  The formal run has
  # exactly ten pairs per candidate mode and alternates call order by pair.
  for (mode in c("oracle", "A", "B")) {
    warm <- capture_proto(compiled_tree, X, nsim, mode, base_seed,
                          return_sim = FALSE, return_ordered = FALSE)
    if (!identical(warm$status, "ok")) {
      stop("Warmup failed for mode=", mode, ", traits=", p, ": ",
           warm$error_message, call. = FALSE)
    }
  }

  for (candidate_mode in c("A", "B")) {
    for (pair in seq_len(repeats)) {
      pair_seed <- base_seed + 1000L * as.integer(p) + as.integer(pair)
      oracle_first <- (pair %% 2L) == 1L
      first_mode <- if (oracle_first) "oracle" else candidate_mode
      second_mode <- if (oracle_first) candidate_mode else "oracle"
      first <- capture_proto(compiled_tree, X, nsim, first_mode, pair_seed,
                             return_sim = FALSE, return_ordered = FALSE)
      second <- capture_proto(compiled_tree, X, nsim, second_mode, pair_seed,
                              return_sim = FALSE, return_ordered = FALSE)
      order <- paste(first_mode, second_mode, sep = ">")
      append_timing(first_mode, first_mode, candidate_mode, p, pair,
                    pair_seed, order, first)
      append_timing(second_mode, second_mode, candidate_mode, p, pair,
                    pair_seed, order, second)
    }
  }

  # Record the prototype-reported workspace and a conservative audit proxy;
  # this does not enter the timing acceptance decision.
  for (mode in c("oracle", "A", "B")) {
    observed <- capture_proto(compiled_tree, X, nsim, mode, base_seed,
                              return_sim = FALSE, return_ordered = FALSE,
                              collect_counters = TRUE)
    reported <- stage2d1_pick(
      observed$value,
      c("peak_workspace_proxy_bytes", "workspace_bytes"),
      NA_real_
    )
    proxy <- stage2d1_workspace_proxy(
      n = n, nsim = nsim, p = p, trait_chunk = trait_chunk,
      simulation_chunk = simulation_chunk, return_sim = FALSE,
      candidate = mode != "oracle", reported = as.numeric(reported)[[1L]]
    )
    proxy$stage <- "Stage 2D1 batch test=TRUE guard"
    proxy$baseline <- "private_oracle"
    proxy$implementation <- mode
    proxy$shape <- "balanced"
    proxy$ncores <- ncores
    memory_rows[[length(memory_rows) + 1L]] <- proxy
  }
}

timings <- stage2d1_append_rows(timing_rows)
parity <- stage2d1_append_rows(parity_rows)
memory <- stage2d1_append_rows(memory_rows)

# Summaries retain both ordinary median/IQR and paired same-seed ratios.  The
# latter is useful when the two implementations are close to timer resolution.
summary_rows <- list()
paired_rows <- list()
for (p in trait_grid) {
  for (candidate_mode in c("A", "B")) {
    one <- timings[timings$traits == p &
                     timings$implementation == "oracle" &
                     timings$candidate_group == candidate_mode, , drop = FALSE]
    two <- timings[timings$traits == p &
                     timings$implementation == candidate_mode &
                     timings$candidate_group == candidate_mode, , drop = FALSE]
    one_ok <- one[one$status == "ok" & is.finite(one$elapsed_s), , drop = FALSE]
    two_ok <- two[two$status == "ok" & is.finite(two$elapsed_s), , drop = FALSE]
    oracle_median <- if (nrow(one_ok)) median(one_ok$elapsed_s) else NA_real_
    candidate_median <- if (nrow(two_ok)) median(two_ok$elapsed_s) else NA_real_
    oracle_iqr <- if (nrow(one_ok)) IQR(one_ok$elapsed_s) else NA_real_
    candidate_iqr <- if (nrow(two_ok)) IQR(two_ok$elapsed_s) else NA_real_
    joined <- merge(
      one_ok[, c("pair", "elapsed_s"), drop = FALSE],
      two_ok[, c("pair", "elapsed_s"), drop = FALSE],
      by = "pair", suffixes = c("_oracle", "_candidate")
    )
    ratios <- if (nrow(joined)) joined$elapsed_s_oracle /
      joined$elapsed_s_candidate else numeric()
    finite_ratios <- ratios[is.finite(ratios)]
    paired_rows[[length(paired_rows) + 1L]] <- data.frame(
      shape = "balanced", n = n, nsim = nsim, traits = p,
      baseline = "private_oracle", candidate = candidate_mode,
      pairs = nrow(joined),
      paired_median_speedup = if (length(finite_ratios))
        median(finite_ratios) else NA_real_,
      paired_IQR_speedup = if (length(finite_ratios))
        IQR(finite_ratios) else NA_real_,
      paired_median_difference_s = if (nrow(joined))
        median(joined$elapsed_s_oracle - joined$elapsed_s_candidate) else NA_real_,
      stringsAsFactors = FALSE
    )
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      shape = "balanced", n = n, nsim = nsim, traits = p,
      baseline = "private_oracle", candidate = candidate_mode,
      oracle_observations = nrow(one), oracle_ok = nrow(one_ok),
      oracle_median_s = oracle_median, oracle_IQR_s = oracle_iqr,
      candidate_observations = nrow(two), candidate_ok = nrow(two_ok),
      candidate_median_s = candidate_median, candidate_IQR_s = candidate_iqr,
      speedup = oracle_median / candidate_median,
      reduction_fraction = 1 - candidate_median / oracle_median,
      stringsAsFactors = FALSE
    )
  }
}
summary <- stage2d1_append_rows(summary_rows)
paired <- stage2d1_append_rows(paired_rows)
if (nrow(summary) && nrow(paired)) {
  summary <- merge(summary, paired,
                   by = c("shape", "n", "nsim", "traits", "baseline",
                          "candidate"), all = TRUE, sort = FALSE)
}

expected_timing_rows <- length(trait_grid) * 2L * repeats * 2L
parity_pass <- nrow(parity) == length(trait_grid) * 2L &&
  nrow(parity) > 0L && all(!is.na(parity$pass) & parity$pass)
timing_complete <- nrow(timings) == expected_timing_rows &&
  all(timings$status == "ok")
warning_free <- nrow(timings) > 0L && !any(nzchar(timings$warning_text)) &&
  nrow(parity) > 0L && !any(nzchar(parity$oracle_warning)) &&
  !any(nzchar(parity$candidate_warning))
slowdown_gt_5 <- nrow(summary) > 0L && any(
  is.finite(summary$reduction_fraction) & summary$reduction_fraction < -0.05
)
run_status <- if (parity_pass && timing_complete && warning_free) "PASS" else "FAIL"

stage2d1_write_csv(timings, file.path(args$out,
                                      "stage2d1_batch_test_true_timings.csv"))
stage2d1_write_csv(parity, file.path(args$out,
                                     "stage2d1_batch_exact_parity.csv"))
stage2d1_write_csv(summary, file.path(args$out,
                                     "stage2d1_batch_test_true_summary.csv"))
stage2d1_write_csv(paired, file.path(args$out,
                                    "stage2d1_batch_test_true_paired.csv"))
stage2d1_write_csv(memory, file.path(args$out,
                                     "stage2d1_batch_memory_proxy.csv"))
stage2d1_write_csv(data.frame(
  stage = "Stage 2D1 batch test=TRUE guard",
  run_mode = if (isTRUE(args$formal)) "formal" else "smoke",
  baseline = "private_oracle",
  shape = "balanced",
  n = n,
  nsim = nsim,
  traits = paste(trait_grid, collapse = ","),
  ncores = ncores,
  repeats_per_candidate = repeats,
  warmup = "one per oracle/A/B and trait",
  alternating_order = "PASS",
  same_seed_per_pair = "PASS",
  exact_private_oracle_parity = if (parity_pass) "PASS" else "FAIL",
  timing_completeness = if (timing_complete) "PASS" else "FAIL",
  warning_free = if (warning_free) "PASS" else "FAIL",
  slowdown_gt_5pct_observed = if (slowdown_gt_5) "YES" else "NO",
  status = run_status,
  production_integration = "NOT_RUN",
  production_code_changed = "NO",
  stringsAsFactors = FALSE
), file.path(args$out, "stage2d1_batch_test_true_status.csv"))

provenance <- stage2d1_provenance(
  args$repo, "Stage 2D1 batch test=TRUE guard", args, hooks_for_provenance
)
provenance$baseline <- "private_oracle"
provenance$fixed_grid <- paste0("balanced/n=", n, "/nsim=", nsim,
                                "/traits=", paste(trait_grid, collapse = ","),
                                "/ncores=1")
provenance$repeats_per_candidate <- repeats
provenance$prototype_modes <- "oracle,A,B"
stage2d1_write_csv(provenance, file.path(args$out,
                                         "stage2d1_batch_test_true_provenance.csv"))
message("[stage2d1] private batch guard complete: ", run_status,
        "; production integration remains NOT_RUN")
