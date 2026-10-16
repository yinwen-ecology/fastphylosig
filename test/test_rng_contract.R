# Cross-install RNG contract: frozen fastphylosig 0.2.0 vs the 0.3.0 candidate.
#
# Run from the repository root:
#   Rscript --vanilla test/test_rng_contract.R <v020-library> <v030-library>
#
# The production API does not expose generated permutation matrices.  This
# script therefore checks a fixed, explicitly materialized permutation matrix
# exactly, then checks internal-RNG calls through their complete K-null output
# and final .Random.seed.  The latter is strong stream evidence, but cannot
# recover or directly compare the hidden permutation rows.

numeric_capture <- function(fit, nsim, label) {
  needed <- c("K_fast", "P_fast", "exceedance_count", "sim.K_fast")
  missing <- setdiff(needed, names(fit))
  if (length(missing)) {
    stop(label, " result is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  sim <- fit[["sim.K_fast"]]
  if (!is.list(sim) || any(lengths(sim) != nsim)) {
    stop(label, " did not return all requested K-null values.", call. = FALSE)
  }
  sim_mat <- vapply(sim, as.numeric, numeric(nsim))
  list(
    K_obs = as.numeric(fit$K_fast),
    K_null = matrix(sim_mat, nrow = nsim),
    P = as.numeric(fit$P_fast),
    exceedance = as.numeric(fit$exceedance_count)
  )
}

make_fixture <- function() {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(20260929L)
  tree <- ape::rtree(24L)
  tree$edge.length <- tree$edge.length / max(tree$edge.length)
  X <- cbind(
    continuous = stats::rnorm(24L),
    binary = as.numeric(sample.int(2L, 24L, replace = TRUE) - 1L),
    categorical = as.numeric(sample.int(4L, 24L, replace = TRUE))
  )
  rownames(X) <- tree$tip.label
  list(tree = tree, X = X)
}

worker <- function(lib, output) {
  lib <- normalizePath(lib, winslash = "/", mustWork = TRUE)
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(
    library(fastphylosig, lib.loc = lib, character.only = FALSE)
  )
  version <- as.character(utils::packageVersion("fastphylosig"))
  fixture <- make_fixture()
  nsim <- 31L

  # Fixed explicit sequence, shared by the two isolated installations.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(82731L)
  permutations <- matrix(NA_integer_, nrow = nsim,
                         ncol = ape::Ntip(fixture$tree))
  permutations[1L, ] <- seq_len(ncol(permutations))
  for (i in 2:nsim) permutations[i, ] <- sample.int(ncol(permutations))

  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(424242L)
  explicit_seed_before <- get(".Random.seed", envir = .GlobalEnv)
  explicit_fit <- suppressWarnings(suppressMessages(fast_k(
    fixture$tree, fixture$X, test = TRUE, nsim = nsim,
    permutations = permutations, return_sim = TRUE,
    verbose = FALSE, progress = FALSE
  )))
  explicit_seed_after <- get(".Random.seed", envir = .GlobalEnv)

  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(424243L)
  internal_fit <- suppressWarnings(suppressMessages(fast_k(
    fixture$tree, fixture$X, test = TRUE, nsim = nsim,
    return_sim = TRUE, verbose = FALSE, progress = FALSE
  )))
  internal_seed_after <- get(".Random.seed", envir = .GlobalEnv)

  saveRDS(list(
    version = version,
    fixture = list(
      edge = fixture$tree$edge,
      edge_length = fixture$tree$edge.length,
      tip_label = fixture$tree$tip.label,
      X = fixture$X
    ),
    permutations = permutations,
    explicit_seed_before = explicit_seed_before,
    explicit_seed_after = explicit_seed_after,
    explicit = numeric_capture(explicit_fit, nsim, "explicit-permutation"),
    internal_seed_after = internal_seed_after,
    internal = numeric_capture(internal_fit, nsim, "internal-RNG")
  ), output)
}

max_abs_difference <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (length(a) != length(b)) return(Inf)
  if (!identical(is.finite(a), is.finite(b))) return(Inf)
  finite <- is.finite(a) & is.finite(b)
  if (!any(finite)) return(0)
  max(abs(a[finite] - b[finite]))
}

max_relative_difference <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (length(a) != length(b)) return(Inf)
  if (!identical(is.finite(a), is.finite(b))) return(Inf)
  finite <- is.finite(a) & is.finite(b)
  if (!any(finite)) return(0)
  max(abs(a[finite] - b[finite]) /
        pmax(abs(a[finite]), .Machine$double.eps))
}

compare_outputs <- function(reference, candidate, label, tol = 1e-10) {
  numeric_fields <- c("K_obs", "K_null")
  relative <- max(vapply(numeric_fields, function(field) {
    max_relative_difference(reference[[field]], candidate[[field]])
  }, numeric(1)))
  absolute <- max(vapply(numeric_fields, function(field) {
    max_abs_difference(reference[[field]], candidate[[field]])
  }, numeric(1)))
  p_equal <- identical(reference$P, candidate$P)
  exceedance_equal <- identical(reference$exceedance, candidate$exceedance)
  pass <- is.finite(relative) && relative <= tol && p_equal && exceedance_equal
  list(
    label = label,
    max_abs_K_difference = absolute,
    max_relative_K_difference = relative,
    P_difference_count = sum(reference$P != candidate$P),
    exceedance_difference_count = sum(reference$exceedance !=
                                        candidate$exceedance),
    pass = pass
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) && identical(args[[1L]], "--worker")) {
  if (length(args) != 3L) stop("worker needs <library> <output.rds>")
  worker(args[[2L]], args[[3L]])
} else {
  if (length(args) != 2L) {
    stop(paste(
      "usage: Rscript --vanilla test/test_rng_contract.R",
      "<fastphylosig-0.2.0-library> <fastphylosig-0.3.0-library>"
    ), call. = FALSE)
  }
  baseline_lib <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
  candidate_lib <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) != 1L) {
    stop("Could not determine the RNG contract script path.", call. = FALSE)
  }
  script <- normalizePath(sub("^--file=", "", script_arg[[1L]]),
                          winslash = "/", mustWork = TRUE)
  baseline_rds <- tempfile(fileext = ".rds")
  candidate_rds <- tempfile(fileext = ".rds")

  run_worker <- function(lib, output_file, label) {
    rscript <- file.path(R.home("bin"), "Rscript.exe")
    if (.Platform$OS.type != "windows") {
      rscript <- file.path(R.home("bin"), "Rscript")
    }
    worker_args <- c(
      "--vanilla", shQuote(script), "--worker", shQuote(lib),
      shQuote(output_file)
    )
    log <- suppressWarnings(system2(
      rscript, args = worker_args, stdout = TRUE, stderr = TRUE
    ))
    status <- attr(log, "status")
    if (is.null(status)) status <- 0L
    if (status != 0L || !file.exists(output_file)) {
      stop(label, " capture failed (status ", status, "): ",
           paste(tail(log, 20L), collapse = "\n"), call. = FALSE)
    }
    invisible(log)
  }

  run_worker(baseline_lib, baseline_rds, "0.2.0")
  run_worker(candidate_lib, candidate_rds, "0.3.0 candidate")
  baseline <- readRDS(baseline_rds)
  candidate <- readRDS(candidate_rds)
  unlink(c(baseline_rds, candidate_rds))

  if (!identical(baseline$version, "0.2.0")) {
    stop("Baseline library must contain fastphylosig 0.2.0; found ",
         baseline$version, call. = FALSE)
  }
  if (!identical(candidate$version, "0.3.0")) {
    stop("Candidate library must contain fastphylosig 0.3.0; found ",
         candidate$version, call. = FALSE)
  }
  fixture_equal <- identical(baseline$fixture, candidate$fixture)
  sequence_equal <- identical(baseline$permutations, candidate$permutations)
  explicit_rng_unchanged <-
    identical(baseline$explicit_seed_before, baseline$explicit_seed_after) &&
    identical(candidate$explicit_seed_before, candidate$explicit_seed_after)
  explicit <- compare_outputs(baseline$explicit, candidate$explicit,
                              "fixed explicit permutation sequence")
  internal <- compare_outputs(baseline$internal, candidate$internal,
                              "internal RNG-generated sequence")
  internal_seed_equal <- identical(baseline$internal_seed_after,
                                   candidate$internal_seed_after)

  report <- data.frame(
    comparison = c(explicit$label, internal$label),
    max_abs_K_difference = c(explicit$max_abs_K_difference,
                             internal$max_abs_K_difference),
    max_relative_K_difference = c(explicit$max_relative_K_difference,
                                  internal$max_relative_K_difference),
    P_difference_count = c(explicit$P_difference_count,
                           internal$P_difference_count),
    exceedance_difference_count = c(explicit$exceedance_difference_count,
                                    internal$exceedance_difference_count),
    stringsAsFactors = FALSE
  )
  cat("STAGE4 RNG CONTRACT: fastphylosig 0.2.0 vs 0.3.0 candidate\n")
  cat("fixture identical:", fixture_equal, "\n")
  cat("explicit permutation sequence identical:", sequence_equal, "\n")
  cat("explicit-permutation calls leave .Random.seed unchanged:",
      explicit_rng_unchanged, "\n")
  cat("internal .Random.seed identical:", internal_seed_equal, "\n")
  print(report, row.names = FALSE)
  cat(paste(
    "Evidence boundary: the API does not expose internally generated",
    "permutation rows; their stream is checked indirectly by complete",
    "K-null outputs plus final .Random.seed.\n"
  ))

  passed <- fixture_equal && sequence_equal && explicit_rng_unchanged &&
    explicit$pass && internal$pass && internal_seed_equal
  cat("STAGE4 RNG CONTRACT", if (passed) "PASS" else "FAIL", "\n")
  if (!passed) stop("The 0.2.0 to candidate RNG contract did not pass.",
                    call. = FALSE)
}
