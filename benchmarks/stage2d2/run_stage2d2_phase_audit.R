#!/usr/bin/env Rscript

# Stage 2D2 private K null-engine phase decomposition.
#
# This is an audit runner only.  It compiles the same-translation-unit
# harness, keeps every fixture fixed within a cell, and writes evidence below
# benchmarks/stage2d2/results/.  It never changes package sources or bindings.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- all_args[grepl("^--file=", all_args)]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else
  file.path(getwd(), "run_stage2d2_phase_audit.R")
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))

args <- commandArgs(trailingOnly = TRUE)
formal <- "--formal" %in% args
smoke <- "--smoke" %in% args || !formal
args <- args[!args %in% c("--formal", "--smoke")]
repo <- if (length(args)) args[[1L]] else
  normalizePath(file.path(script_dir, "..", ".."), winslash = "/",
                mustWork = FALSE)
out_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path(repo, "benchmarks", "stage2d2", "results", "phase-audit")
repo <- normalizePath(repo, winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

parse_ints <- function(value, name, minimum = 1L) {
  out <- suppressWarnings(as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])))
  if (!length(out) || anyNA(out) || any(out < minimum)) {
    stop(name, " must contain comma-separated integers >= ", minimum, ".",
         call. = FALSE)
  }
  unique(out)
}

grid <- function(name, default, smoke_default, minimum = 1L) {
  value <- env_or(name, "")
  if (nzchar(value)) parse_ints(value, name, minimum) else
    if (smoke) smoke_default else default
}

lib <- env_or("FASTPHYLOSIG_STAGE2C_LIB", "")
if (!nzchar(lib) || !dir.exists(file.path(lib, "fastphylosig"))) {
  guesses <- Sys.glob("C:/Users/wenyi/AppData/Local/Temp/fastphylosig-stage2c-*/library")
  guesses <- guesses[dir.exists(file.path(guesses, "fastphylosig"))]
  if (length(guesses)) lib <- guesses[[length(guesses)]]
}
if (!nzchar(lib) || !dir.exists(file.path(lib, "fastphylosig"))) {
  stop("Set FASTPHYLOSIG_STAGE2C_LIB to the audited fastphylosig library.",
       call. = FALSE)
}

cpp <- env_or("FASTPHYLOSIG_STAGE2D2_CPP", file.path(
  repo, "benchmarks", "stage2d2", "stage2d2_phase_harness.cpp"
))
if (!file.exists(cpp)) stop("Stage 2D2 harness source is missing: ", cpp,
                             call. = FALSE)

suppressPackageStartupMessages({
  library(fastphylosig, lib.loc = lib)
  library(ape)
  library(Rcpp)
})

private_env <- new.env(parent = globalenv())
Rcpp::sourceCpp(normalizePath(cpp, winslash = "/", mustWork = TRUE),
                env = private_env, rebuild = TRUE, showOutput = FALSE,
                verbose = FALSE)
phase_profile <- private_env$stage2d2_phase_profile
if (!is.function(phase_profile)) {
  stop("The compiled Stage 2D2 harness did not expose stage2d2_phase_profile.",
       call. = FALSE)
}
compile_tree <- getFromNamespace("compile_tree_cpp", "fastphylosig")

make_fixture <- function(n, seed = 206020L) {
  set.seed(as.integer(seed + n))
  tree <- ape::reorder.phylo(ape::rtree(n), order = "postorder")
  tree$tip.label <- paste0("sp", seq_len(n))
  tree$edge.length <- 0.5 + seq_len(nrow(tree$edge)) / nrow(tree$edge)
  x <- matrix(
    sin(seq_len(n) + 1 / 11) + cos(seq_len(n) / 7 + 1) +
      seq_len(n) / (n + 17),
    nrow = n, ncol = 1L,
    dimnames = list(tree$tip.label, "trait1")
  )
  compiled <- compile_tree(tree$edge, tree$edge.length, ape::Ntip(tree))
  list(tree = tree, compiled = compiled, X = x)
}

capture <- function(fun) {
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(
    withCallingHandlers(fun(), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) structure(list(error = conditionMessage(e)),
                                  class = "stage2d2_error")
  )
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(value, "stage2d2_error")) {
    return(list(value = NULL, elapsed_s = max(0, elapsed), status = "error",
                warning_text = paste(unique(warnings), collapse = " | "),
                error_message = value$error))
  }
  list(value = value, elapsed_s = max(0, elapsed), status = "ok",
       warning_text = paste(unique(warnings), collapse = " | "),
       error_message = "")
}

call_profile <- function(fixture, phase, nsim, n_threads, simulation_chunk,
                         seed, trait_chunk = 1L, trace_rows = 3L) {
  # Reset the R stream for each observation.  Generation and full pipeline
  # therefore receive the identical seed and can be compared row-for-row.
  set.seed(as.integer(seed))
  capture(function() phase_profile(
    fixture$compiled, fixture$X, nsim = as.integer(nsim),
    trait_chunk = as.integer(trait_chunk), include_observed = TRUE,
    n_threads = as.integer(n_threads),
    simulation_chunk = as.integer(simulation_chunk), phase = phase,
    controlled_count = 8L, trace_rows = as.integer(trace_rows)
  ))
}

trace_equal <- function(left, right) {
  !is.null(left) && !is.null(right) &&
    isTRUE(all.equal(as.integer(left), as.integer(right),
                     check.attributes = TRUE))
}

source_sha256 <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256"))
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    return(paste(as.character(openssl::sha256(file(path))), collapse = ""))
  }
  NA_character_
}

git_value <- function(...) {
  value <- tryCatch(system2("git", c(...), stdout = TRUE, stderr = FALSE),
                    error = function(e) character())
  if (length(value)) paste(value, collapse = "\n") else "UNAVAILABLE"
}

n_grid <- grid("FASTPHYLOSIG_STAGE2D2_PHASE_N",
               c(500L, 2000L, 5000L, 10000L), c(50L), minimum = 2L)
nsim_grid <- grid("FASTPHYLOSIG_STAGE2D2_PHASE_NSIM",
                  c(999L, 9999L), c(19L), minimum = 1L)
threads_grid <- grid("FASTPHYLOSIG_STAGE2D2_PHASE_THREADS",
                     c(1L, 2L), c(1L), minimum = 1L)
simulation_chunk <- as.integer(env_or("FASTPHYLOSIG_STAGE2D2_PHASE_CHUNK", "128"))
if (is.na(simulation_chunk) || simulation_chunk < 1L) {
  stop("FASTPHYLOSIG_STAGE2D2_PHASE_CHUNK must be positive.", call. = FALSE)
}
repeats <- as.integer(env_or("FASTPHYLOSIG_STAGE2D2_PHASE_REPEATS",
                             if (smoke) "2" else "10"))
if (is.na(repeats) || repeats < 1L) stop("repeats must be positive.", call. = FALSE)
trace_rows <- as.integer(env_or("FASTPHYLOSIG_STAGE2D2_PHASE_TRACE_ROWS", "3"))
if (is.na(trace_rows) || trace_rows < 1L) stop("trace_rows must be positive.", call. = FALSE)
base_seed <- as.integer(env_or("FASTPHYLOSIG_STAGE2D2_PHASE_SEED", "206020"))
if (is.na(base_seed)) stop("invalid base seed", call. = FALSE)

rows <- list()
parity <- list()
phases <- c("generation", "evaluation", "full")

for (n in n_grid) {
  message("[stage2d2] fixture n=", n)
  fixture <- make_fixture(n)
  for (nsim in nsim_grid) for (threads in threads_grid) {
    # Warm each phase once.  Warmups are never included in the output.
    for (phase in phases) {
      warm <- call_profile(fixture, phase, nsim, threads, simulation_chunk,
                           base_seed + n + nsim + threads, trace_rows)
      if (warm$status != "ok") stop("warmup failed: ", phase, "/n=", n,
                                    "/nsim=", nsim, ": ", warm$error_message,
                                    call. = FALSE)
    }
    for (rep in seq_len(repeats)) {
      order <- if (rep %% 2L) phases else rev(phases)
      observations <- list()
      for (phase in order) {
        seed <- base_seed + n + nsim + rep
        measured <- call_profile(fixture, phase, nsim, threads,
                                 simulation_chunk, seed, trace_rows = trace_rows)
        value <- measured$value
        value_or <- function(name, default = NA_real_) {
          if (measured$status == "ok" && !is.null(value[[name]]) &&
              length(value[[name]]) == 1L) value[[name]] else default
        }
        if (measured$status == "ok") {
          if (phase == "generation") {
            expected_draws <- (nsim - 1L) * (n - 1L)
            draw_pass <- identical(as.numeric(value$rng_draws),
                                   as.numeric(expected_draws))
          } else if (phase == "full") {
            expected_draws <- (nsim - 1L) * (n - 1L)
            draw_pass <- identical(as.numeric(value$rng_draws),
                                   as.numeric(expected_draws))
          } else draw_pass <- TRUE
          observations[[phase]] <- list(
            value = value, draw_pass = draw_pass
          )
        }
        rows[[length(rows) + 1L]] <- data.frame(
          n = n, nsim = nsim, n_threads_requested = threads,
          n_threads_effective = value_or("n_threads_effective", NA_integer_),
          openmp_compiled = value_or("openmp_compiled", NA),
          openmp_max_threads = value_or("openmp_max_threads", NA_integer_),
          simulation_chunk = simulation_chunk, replicate = rep, phase = phase,
          elapsed_s = measured$elapsed_s, native_elapsed_s = if (
            measured$status == "ok"
          ) value$elapsed_s else NA_real_,
          status = measured$status, warning_text = measured$warning_text,
          error_message = measured$error_message,
          rng_draws = value_or("rng_draws"),
          fisher_yates_swaps = value_or("fisher_yates_swaps"),
          checksum = if (measured$status == "ok") {
            if (phase == "full") value_or("generation_checksum", NA_character_)
            else value_or("checksum", NA_character_)
          } else NA_character_,
          trace_rows = value_or("trace_rows", NA_integer_),
          peak_index_bytes = value_or("peak_index_bytes"),
          observed_s = if (phase == "full" && measured$status == "ok")
            value_or("observed_s") else NA_real_,
          generation_and_chunk_s = if (phase == "full" && measured$status == "ok")
            value_or("generation_and_chunk_s") else NA_real_,
          evaluation_s = if (phase == "full" && measured$status == "ok")
            value_or("evaluation_s") else if (phase == "evaluation" &&
                                        measured$status == "ok") value_or("evaluation_s") else NA_real_,
          accounting_s = if (phase == "full" && measured$status == "ok")
            value_or("accounting_s") else NA_real_,
          stringsAsFactors = FALSE
        )
      }
      gen <- observations$generation
      full <- observations$full
      parity[[length(parity) + 1L]] <- data.frame(
        n = n, nsim = nsim, n_threads_requested = threads, replicate = rep,
        generation_ok = !is.null(gen) && isTRUE(gen$draw_pass),
        full_ok = !is.null(full) && isTRUE(full$draw_pass),
        generation_full_trace_equal = !is.null(gen) && !is.null(full) &&
          trace_equal(gen$value$trace, full$value$trace),
        generation_full_checksum_equal = !is.null(gen) && !is.null(full) &&
          identical(as.character(gen$value$checksum),
                    as.character(full$value$generation_checksum)),
        generation_full_draws_equal = !is.null(gen) && !is.null(full) &&
          identical(as.numeric(gen$value$rng_draws),
                    as.numeric(full$value$rng_draws)),
        stringsAsFactors = FALSE
      )
    }
  }
}

timings <- do.call(rbind, rows)
parity <- do.call(rbind, parity)
utils::write.csv(timings, file.path(out_dir, "stage2d2_phase_timings.csv"),
                 row.names = FALSE)
utils::write.csv(parity, file.path(out_dir, "stage2d2_rng_parity.csv"),
                 row.names = FALSE)

key <- paste(timings$n, timings$nsim, timings$n_threads_requested,
             timings$phase, sep = "\r")
groups <- split(seq_len(nrow(timings)), key, drop = TRUE)
summary_rows <- lapply(groups, function(ii) {
  z <- timings[ii, , drop = FALSE]
  ok <- z$status == "ok" & is.finite(z$elapsed_s) &
    is.finite(z$native_elapsed_s)
  values <- z$native_elapsed_s[ok]
  wall <- z$elapsed_s[ok]
  z[1L, c("n", "nsim", "n_threads_requested", "n_threads_effective",
           "openmp_compiled", "openmp_max_threads", "simulation_chunk",
           "phase"), drop = FALSE]
  data.frame(
    z[1L, c("n", "nsim", "n_threads_requested", "n_threads_effective",
             "openmp_compiled", "openmp_max_threads", "simulation_chunk",
             "phase"), drop = FALSE],
    observations = nrow(z), ok_observations = sum(ok),
    median_s = if (length(values)) median(values) else NA_real_,
    IQR_s = if (length(values)) IQR(values) else NA_real_,
    wall_median_s = if (length(wall)) median(wall) else NA_real_,
    wall_IQR_s = if (length(wall)) IQR(wall) else NA_real_,
    native_median_s = if (length(values)) median(values) else NA_real_,
    native_IQR_s = if (length(values)) IQR(values) else NA_real_,
    warnings = paste(unique(z$warning_text[nzchar(z$warning_text)]), collapse = " | "),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
utils::write.csv(summary, file.path(out_dir, "stage2d2_phase_summary.csv"),
                 row.names = FALSE)

expected_rows <- length(n_grid) * length(nsim_grid) * length(threads_grid) *
  repeats * length(phases)
no_failures <- nrow(timings) == expected_rows && all(timings$status == "ok") &&
  all(!nzchar(timings$warning_text))
rng_pass <- nrow(parity) > 0L && all(
  parity$generation_ok & parity$full_ok & parity$generation_full_trace_equal &
    parity$generation_full_checksum_equal & parity$generation_full_draws_equal
)
requested_two <- timings$n_threads_requested >= 2L
openmp_two_thread <- if (!any(requested_two)) {
  "NOT_REQUESTED"
} else if (all(timings$openmp_compiled[requested_two] &
              timings$n_threads_effective[requested_two] >= 2L)) {
  "PASS"
} else {
  "UNSUPPORTED"
}
usable_for_requested_threads <- no_failures && openmp_two_thread %in%
  c("PASS", "NOT_REQUESTED")
utils::write.csv(data.frame(
  phase_decomposition = if (no_failures) "PASS" else "FAIL",
  rng_replay = if (rng_pass) "PASS" else "FAIL",
  openmp_two_thread = openmp_two_thread,
  requested_thread_measurement = if (usable_for_requested_threads) "COMPLETE" else
    "INCOMPLETE_OPENMP",
  formal_grid = if (formal) "YES" else "NO_SMOKE",
  production_code_changed = "NO",
  stringsAsFactors = FALSE
), file.path(out_dir, "stage2d2_phase_status.csv"), row.names = FALSE)

provenance <- c(
  paste0("repo=", repo),
  paste0("git_head=", git_value("-C", repo, "rev-parse", "HEAD")),
  paste0("git_status=", git_value("-C", repo, "status", "--short")),
  paste0("harness=", normalizePath(cpp, winslash = "/", mustWork = TRUE)),
  paste0("harness_sha256=", source_sha256(cpp)),
  paste0("production_evaluator=src/k_permutation.cpp::kperm::compute_one"),
  paste0("production_rng=src/k_permutation.cpp::R::runif + descending Fisher-Yates"),
  paste0("identity_convention=include_observed=TRUE makes replicate 1 identity; all following rows draw n-1 uniforms"),
  paste0("controlled_fixture=bounded eight-pattern deterministic cyclic/reversed permutations; prepared outside evaluator timer"),
  paste0("fixture_tree=ape::rtree with fixed seed per n, postorder, fixed positive edge lengths"),
  paste0("n_grid=", paste(n_grid, collapse = ",")),
  paste0("nsim_grid=", paste(nsim_grid, collapse = ",")),
  paste0("threads_grid=", paste(threads_grid, collapse = ",")),
  paste0("openmp_gate=", openmp_two_thread),
  paste0("repeats=", repeats),
  paste0("simulation_chunk=", simulation_chunk),
  paste0("measurement=warmup; serialized alternating phase order; median/IQR"),
  paste0("generated_at=", format(Sys.time(), tz = "UTC"))
)
writeLines(provenance, file.path(out_dir, "stage2d2_phase_provenance.txt"))
message("[stage2d2] phase audit complete: ", out_dir)
print(utils::read.csv(file.path(out_dir, "stage2d2_phase_status.csv"),
                      stringsAsFactors = FALSE))
