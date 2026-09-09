#!/usr/bin/env Rscript

# Stage 2D1 private audit harness.  This file never changes package sources or
# bindings.  It compiles the sibling C++ prototype in a private environment,
# compares it with the installed production kernel, and writes evidence only
# below benchmarks/stage2d1/results/operation-audit/.

args <- commandArgs(trailingOnly = TRUE)
repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  repo_root <- normalizePath(file.path(repo_root, "..", "..", ".."),
                             winslash = "/", mustWork = TRUE)
}

out_dir <- file.path(repo_root, "benchmarks", "stage2d1", "results",
                     "operation-audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lib <- Sys.getenv("FASTPHYLOSIG_STAGE2C_LIB", unset = "")
if (!nzchar(lib) || !dir.exists(file.path(lib, "fastphylosig"))) {
  guesses <- Sys.glob("C:/Users/wenyi/AppData/Local/Temp/fastphylosig-stage2c-*/library")
  guesses <- guesses[dir.exists(file.path(guesses, "fastphylosig"))]
  if (length(guesses)) lib <- guesses[[length(guesses)]]
}
if (!nzchar(lib) || !dir.exists(file.path(lib, "fastphylosig"))) {
  stop("Set FASTPHYLOSIG_STAGE2C_LIB to the audited fastphylosig library.")
}

suppressPackageStartupMessages({
  library(fastphylosig, lib.loc = lib)
  library(ape)
  library(Rcpp)
})

proto_env <- new.env(parent = globalenv())
proto_source <- Sys.getenv(
  "FASTPHYLOSIG_STAGE2D1_CPP",
  unset = "C:/Users/wenyi/AppData/Local/Temp/fastphylosig-stage2d1-k-prototype.cpp"
)
if (!file.exists(proto_source)) {
  stop("Stage 2D1 prototype source is missing: ", proto_source)
}
Rcpp::sourceCpp(
  proto_source,
  env = proto_env,
  rebuild = TRUE,
  showOutput = FALSE,
  verbose = FALSE
)
proto <- proto_env$stage2d1_k_permutation_prototype
prod_kernel <- getFromNamespace("fast_k_tree_permutation_cpp", "fastphylosig")
compile_tree <- getFromNamespace("compile_tree_cpp", "fastphylosig")

make_fixture <- function(n, seed = 20260909L, p = 1L) {
  set.seed(seed)
  tree <- ape::reorder.phylo(ape::rtree(n), "pruningwise")
  X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p,
              dimnames = list(tree$tip.label, paste0("trait", seq_len(p))))
  compiled <- compile_tree(tree$edge, tree$edge.length, ape::Ntip(tree))
  list(tree = tree, compiled = compiled, X = X)
}

make_permutations <- function(n, nsim, seed = 981L) {
  set.seed(seed)
  out <- matrix(seq_len(n), nrow = nsim, ncol = n, byrow = TRUE)
  if (nsim > 1L) {
    for (i in 2:nsim) out[i, ] <- sample.int(n, n, replace = FALSE)
  }
  out
}

max_diff <- function(x, y) {
  if (length(x) != length(y)) return(Inf)
  if (!length(x)) return(0)
  z <- abs(as.numeric(x) - as.numeric(y))
  if (all(is.na(z))) 0 else max(z, na.rm = TRUE)
}

same_payload <- function(a, b, tolerance = 0) {
  fields <- c("K", "P", "exceedance_count", "MCSE_P")
  ok <- vapply(fields, function(nm) {
    if (is.null(a[[nm]]) || is.null(b[[nm]])) return(FALSE)
    isTRUE(all.equal(as.numeric(a[[nm]]), as.numeric(b[[nm]]),
                     tolerance = tolerance, check.attributes = FALSE))
  }, logical(1))
  ordered_ok <- TRUE
  if (!is.null(a$sim_K) || !is.null(b$sim_K)) {
    ordered_ok <- isTRUE(all.equal(as.numeric(a$sim_K), as.numeric(b$sim_K),
                                   tolerance = tolerance,
                                   check.attributes = FALSE))
  }
  c(setNames(ok, fields), ordered_null_K = ordered_ok,
    max_abs_sim_K = max_diff(a$sim_K, b$sim_K))
}

call_prod <- function(f, compiled, X, nsim, chunk, n_threads,
                      supplied = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  f(compiled, X, nsim = nsim, permutations = supplied,
    trait_chunk = 1L, return_sim = TRUE, include_observed = is.null(supplied),
    n_threads = n_threads, simulation_chunk = chunk)
}

call_proto <- function(f, compiled, X, nsim, chunk, n_threads,
                       supplied = NULL, seed = NULL, mode = "A") {
  if (!is.null(seed)) set.seed(seed)
  f(compiled, X, nsim = nsim, permutations = supplied,
    trait_chunk = 1L, return_sim = TRUE, return_ordered = TRUE,
    include_observed = is.null(supplied), n_threads = n_threads,
    simulation_chunk = chunk, mode = mode)
}

run_controlled_audit <- function() {
  f <- make_fixture(50L, p = 2L)
  perms <- make_permutations(50L, 31L)
  rows <- list()
  for (threads in c(1L, 2L)) {
    oracle <- call_prod(prod_kernel, f$compiled, f$X, 31L, 7L, threads,
                        supplied = perms)
    for (mode in c("oracle", "A", "B")) {
      got <- call_proto(proto, f$compiled, f$X, 31L, 7L, threads,
                        supplied = perms, mode = mode)
      chk <- same_payload(oracle, got, tolerance = 0)
      rows[[length(rows) + 1L]] <- data.frame(
        audit = "controlled",
        mode = mode,
        n = 50L,
        nsim = 31L,
        n_threads = threads,
        K = chk[["K"]], P = chk[["P"]],
        exceedance_count = chk[["exceedance_count"]],
        MCSE_P = chk[["MCSE_P"]],
        ordered_null_K = chk[["ordered_null_K"]],
        max_abs_sim_K = unname(chk[["max_abs_sim_K"]]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

run_rng_audit <- function() {
  rows <- list()
  cases <- list(
    c(n = 50L, nsim = 1L), c(n = 50L, nsim = 2L),
    c(n = 50L, nsim = 199L), c(n = 500L, nsim = 199L),
    c(n = 2000L, nsim = 199L), c(n = 5000L, nsim = 199L)
  )
  chunks <- c(1L, 2L, 7L, 127L, 128L, 129L, 256L)
  for (case in cases) {
    n <- unname(case[["n"]]); nsim <- unname(case[["nsim"]])
    f <- make_fixture(n, seed = 7700L + n)
    for (chunk in chunks) {
      if (chunk > nsim + 1L) next
      for (threads in c(1L, 2L)) {
        seed <- 13000L + n + nsim + chunk + threads
        oracle <- call_prod(prod_kernel, f$compiled, f$X, nsim, chunk,
                            threads, seed = seed)
        for (mode in c("A", "B")) {
          got <- call_proto(proto, f$compiled, f$X, nsim, chunk,
                            threads, seed = seed, mode = mode)
          chk <- same_payload(oracle, got, tolerance = 1e-12)
          rows[[length(rows) + 1L]] <- data.frame(
            audit = "rng_replay",
            mode = mode,
            n = n,
            nsim = nsim,
            simulation_chunk = chunk,
            n_threads = threads,
            K = chk[["K"]], P = chk[["P"]],
            exceedance_count = chk[["exceedance_count"]],
            MCSE_P = chk[["MCSE_P"]],
            ordered_null_K = chk[["ordered_null_K"]],
            max_abs_sim_K = unname(chk[["max_abs_sim_K"]]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}

run_batch_audit <- function() {
  rows <- list()
  f <- make_fixture(500L, p = 1L)
  for (p in c(1L, 8L, 32L, 100L)) {
    f <- make_fixture(500L, seed = 707L, p = p)
    seed <- 4400L + p
    old <- call_prod(prod_kernel, f$compiled, f$X, 199L, 128L, 1L,
                     seed = seed)
    for (mode in c("A", "B")) {
      got <- call_proto(proto, f$compiled, f$X, 199L, 128L, 1L,
                        seed = seed, mode = mode)
      chk <- same_payload(old, got, tolerance = 1e-12)
      rows[[length(rows) + 1L]] <- data.frame(
        audit = "batch",
        mode = mode, n = 500L, traits = p,
        nsim = 199L, K = chk[["K"]], P = chk[["P"]],
        ordered_null_K = chk[["ordered_null_K"]],
        max_abs_sim_K = unname(chk[["max_abs_sim_K"]]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

run_counter_audit <- function() {
  f <- make_fixture(5000L, seed = 2202L, p = 1L)
  rows <- list()
  for (mode in c("oracle", "A", "B")) {
    set.seed(62001L)
    got <- proto(f$compiled, f$X, nsim = 999L, trait_chunk = 1L,
                 return_sim = FALSE, return_ordered = FALSE,
                 include_observed = TRUE, n_threads = 1L,
                 simulation_chunk = 128L, mode = mode)
    z <- got$counters
    rows[[length(rows) + 1L]] <- data.frame(
      mode = mode,
      nsim = 999L,
      null_replicates = z$null_replicates,
      index_buffer_allocations = z$index_buffer_allocations,
      index_initialization_elements = z$index_initialization_elements,
      fisher_yates_swaps = z$fisher_yates_swaps,
      generation_to_chunk_elements = z$generation_to_chunk_elements,
      chunk_to_worker_elements = z$chunk_to_worker_elements,
      controlled_to_worker_elements = z$controlled_to_worker_elements,
      index_reads = z$index_reads,
      trait_reads = z$trait_reads,
      trait_gather_allocations = z$trait_gather_allocations,
      compute_workspace_allocations = z$compute_workspace_allocations,
      compute_workspace_bytes = z$compute_workspace_bytes,
      compute_workspace_peak_bytes = z$compute_workspace_peak_bytes,
      permutation_chunk_bytes = z$permutation_chunk_bytes,
      reusable_generation_bytes = z$reusable_generation_bytes,
      worker_index_bytes = z$worker_index_bytes,
      peak_workspace_proxy_bytes = z$peak_workspace_proxy_bytes,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

run_edge_audit <- function() {
  rows <- list()
  cases <- list(
    list(label = "two_tip", n = 2L, p = 1L, values = NULL),
    list(label = "constant", n = 12L, p = 1L, values = "constant"),
    list(label = "near_constant", n = 12L, p = 1L, values = "near"),
    list(label = "large_offset", n = 12L, p = 2L, values = "offset")
  )
  for (case in cases) {
    f <- make_fixture(case$n, seed = 1001L + case$n, p = case$p)
    if (identical(case$values, "constant")) f$X[,] <- 1
    if (identical(case$values, "near")) {
      f$X[, 1L] <- 1 + seq_len(case$n) * 1e-12
    }
    if (identical(case$values, "offset")) {
      f$X[,] <- 1e12 + matrix(stats::rnorm(case$n * case$p),
                               nrow = case$n, ncol = case$p)
    }
    for (mode in c("A", "B")) {
      old <- call_prod(prod_kernel, f$compiled, f$X, 19L, 7L, 1L,
                       seed = 71000L + case$n)
      got <- call_proto(proto, f$compiled, f$X, 19L, 7L, 1L,
                        seed = 71000L + case$n, mode = mode)
      chk <- same_payload(old, got, tolerance = 1e-12)
      rows[[length(rows) + 1L]] <- data.frame(
        audit = case$label, mode = mode,
        K = chk[["K"]], P = chk[["P"]],
        exceedance_count = chk[["exceedance_count"]],
        MCSE_P = chk[["MCSE_P"]],
        ordered_null_K = chk[["ordered_null_K"]],
        max_abs_sim_K = unname(chk[["max_abs_sim_K"]]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

timed <- function(fun) {
  t <- system.time(fun)[["elapsed"]]
  as.numeric(t)
}

run_timing_audit <- function(reps = 10L) {
  rows <- list()
  for (n in c(500L, 2000L, 5000L, 10000L)) {
    f <- make_fixture(n, seed = 9000L + n, p = 1L)
    for (nsim in c(199L, 999L, 9999L)) {
      for (rep in seq_len(reps)) {
        seed <- 800000L + n + nsim + rep
        order <- if (rep %% 2L) c("oracle", "A", "B") else
          c("B", "A", "oracle")
        for (mode in order) {
          set.seed(seed)
          elapsed <- if (mode == "oracle") {
            timed(prod_kernel(f$compiled, f$X, nsim = nsim,
                              trait_chunk = 1L, return_sim = FALSE,
                              include_observed = TRUE, n_threads = 1L,
                              simulation_chunk = 128L))
          } else {
            timed(proto(f$compiled, f$X, nsim = nsim,
                        trait_chunk = 1L, return_sim = FALSE,
                        return_ordered = FALSE, include_observed = TRUE,
                        n_threads = 1L, simulation_chunk = 128L,
                        mode = mode))
          }
          rows[[length(rows) + 1L]] <- data.frame(
            audit = "timing",
            n = n, nsim = nsim, rep = rep, mode = mode,
            elapsed_seconds = elapsed,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}

quick_only <- "--quick" %in% args
timing_arg <- grep("^--timing-reps=", args, value = TRUE)
timing_reps <- if (length(timing_arg)) {
  as.integer(sub("^--timing-reps=", "", timing_arg[[1L]]))
} else 10L
controlled <- run_controlled_audit()
rng <- run_rng_audit()
batch <- run_batch_audit()
counters <- run_counter_audit()
edge <- run_edge_audit()
utils::write.csv(controlled, file.path(out_dir, "controlled_parity.csv"),
                 row.names = FALSE)
utils::write.csv(rng, file.path(out_dir, "rng_replay.csv"), row.names = FALSE)
utils::write.csv(batch, file.path(out_dir, "batch_parity.csv"), row.names = FALSE)
utils::write.csv(counters, file.path(out_dir, "operation_counters.csv"),
                 row.names = FALSE)
utils::write.csv(edge, file.path(out_dir, "edge_parity.csv"), row.names = FALSE)

if (!quick_only) {
  if (!is.finite(timing_reps) || timing_reps < 1L) {
    stop("--timing-reps must be a positive integer")
  }
  timing <- run_timing_audit(reps = timing_reps)
  utils::write.csv(timing, file.path(out_dir, "timings.csv"), row.names = FALSE)
}

status <- data.frame(
  controlled_all_pass = all(controlled$K & controlled$P &
                            controlled$exceedance_count & controlled$MCSE_P &
                            controlled$ordered_null_K),
  rng_all_pass = all(rng$K & rng$P & rng$exceedance_count & rng$MCSE_P &
                     rng$ordered_null_K),
  batch_all_pass = all(batch$K & batch$P & batch$ordered_null_K),
  edge_all_pass = all(edge$K & edge$P & edge$exceedance_count &
                      edge$MCSE_P & edge$ordered_null_K),
  production_code_changed = FALSE,
  stringsAsFactors = FALSE
)
utils::write.csv(status, file.path(out_dir, "prototype_status.csv"),
                 row.names = FALSE)
print(status)
print(counters)
