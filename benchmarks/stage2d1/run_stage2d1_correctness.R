#!/usr/bin/env Rscript

# Stage 2D1 correctness and replay gates for a private K permutation prototype.
# The production namespace is used as an old oracle.  A prototype is loaded
# only through an explicitly configured hook and is never installed here.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- all_args[grepl("^--file=", all_args)]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else
  file.path(getwd(), "run_stage2d1_correctness.R")
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
source(file.path(script_dir, "stage2d1_common.R"), local = TRUE)

args <- stage2d1_args(commandArgs(trailingOnly = TRUE))
if (!isTRUE(args$smoke) && !isTRUE(args$formal)) {
  stop("Correctness audit is explicit: add --formal or --smoke.", call. = FALSE)
}
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")
assign("stage2d1_smoke", isTRUE(args$smoke), envir = .GlobalEnv)
Sys.setenv(FASTPHYLOSIG_STAGE2D1_RETURN_ORDERED = "true")

ns <- stage2d1_load_package(args$repo)
hooks <- stage2d1_load_hooks(args$repo)
stage2d1_frozen_contract(args$repo)
prov <- stage2d1_provenance(args$repo, "Stage 2D1 correctness", args, hooks)
stage2d1_write_csv(prov, file.path(args$out, "stage2d1_correctness_provenance.csv"))

status_path <- file.path(args$out, "stage2d1_correctness_status.csv")
if (!is.function(hooks$controlled) || !is.function(hooks$rng)) {
  reason <- paste0(
    "Prototype hooks are required before correctness can run. Configure ",
    "FASTPHYLOSIG_STAGE2D1_HOOK_FILE or FASTPHYLOSIG_STAGE2D1_PROTO_CPP; ",
    "expected controlled hook `", hooks$controlled_name,
    "` and RNG hook `", hooks$rng_name, "`."
  )
  stage2d1_write_csv(data.frame(
    status = "NOT_RUN", core_gates = "NOT_RUN", reason = reason,
    production_code_changed = "NO", stringsAsFactors = FALSE
  ), status_path)
  stop(reason, call. = FALSE)
}

old_engine <- stage2d1_production_cpp(ns)

formal_n <- c(50L, 500L, 2000L, 5000L)
formal_nsim <- c(1L, 2L, 199L, 999L)
formal_seeds <- c(1L, 20260901L, 99173L)
n_grid <- stage2d1_grid("FASTPHYLOSIG_STAGE2D1_N", formal_n,
                        c(50L, 500L), minimum = 2L)
nsim_grid <- stage2d1_grid("FASTPHYLOSIG_STAGE2D1_NSIM", formal_nsim,
                           c(1L, 2L, 19L), minimum = 1L)
seed_value <- stage2d1_env("FASTPHYLOSIG_STAGE2D1_SEEDS")
seeds <- if (nzchar(seed_value)) stage2d1_ints(seed_value, "SEEDS") else
  if (isTRUE(args$smoke)) formal_seeds[[1L]] else formal_seeds
shape_value <- stage2d1_env("FASTPHYLOSIG_STAGE2D1_SHAPES")
shapes <- if (nzchar(shape_value)) {
  unique(trimws(strsplit(shape_value, ",", fixed = TRUE)[[1L]]))
} else if (isTRUE(args$smoke)) "balanced" else "balanced"
if (any(!shapes %in% c("balanced", "random", "pectinate"))) {
  stop("FASTPHYLOSIG_STAGE2D1_SHAPES contains an unknown shape.", call. = FALSE)
}
ncores <- stage2d1_grid(
  "FASTPHYLOSIG_STAGE2D1_NCORES", c(1L, 2L), c(1L), minimum = 1L
)
if (any(ncores > 1L) && isTRUE(args$smoke)) ncores <- 1L
detected_cores <- tryCatch(parallel::detectCores(logical = TRUE),
                           error = function(e) 1L)
if (is.finite(detected_cores) && detected_cores >= 1L) {
  ncores <- ncores[ncores <= detected_cores]
}
if (!length(ncores)) ncores <- 1L
trait_chunk <- as.integer(stage2d1_env("FASTPHYLOSIG_STAGE2D1_TRAIT_CHUNK", "64"))
if (is.na(trait_chunk) || trait_chunk < 1L) stop("invalid trait chunk", call. = FALSE)

rows <- list()
controlled_null_rows <- list()
rng_rows <- list()
failure_rows <- list()
na_rows <- list()
batch_rows <- list()
add_row <- function(...) rows[[length(rows) + 1L]] <<- data.frame(...,
  stringsAsFactors = FALSE)

compare_pair <- function(mode, candidate, shape, n, nsim, seed, threads,
                         chunk, X, permutations = NULL, include_observed,
                         old_fun, new_fun) {
  old <- stage2d1_capture(function() stage2d1_call_engine(
    old_fun, mode, compiled_tree, X, nsim, permutations, trait_chunk,
    return_sim = TRUE, include_observed = include_observed, ncores = threads,
    simulation_chunk = chunk, seed = seed
  ))
  new <- stage2d1_capture(function() stage2d1_call_engine(
    new_fun, mode, compiled_tree, X, nsim, permutations, trait_chunk,
    return_sim = TRUE, include_observed = include_observed, ncores = threads,
    simulation_chunk = chunk, seed = seed
  ))
  cmp <- stage2d1_compare_captures(old, new, nsim, ncol(X), numeric_tolerance = 0)
  list(old = old, new = new, comparison = cmp)
}

record_pair <- function(check, mode, candidate, shape, n, nsim, seed,
                        threads, chunk, result) {
  cmp <- result$comparison
  detail <- if (isTRUE(cmp$pass)) "exact result/status/warning parity" else
    paste0("status=", cmp$left_status, "/", cmp$right_status,
           "; warning=", cmp$warning_pass,
           "; max_sim_abs=", cmp$result$max_abs[["sim_K"]])
  add_row(
    check = check, mode = mode, candidate = candidate, shape = shape,
    n = n, nsim = nsim, seed = if (is.null(seed)) NA_integer_ else seed,
    ncores = threads, simulation_chunk = chunk, pass = isTRUE(cmp$pass),
    old_status = cmp$left_status, new_status = cmp$right_status,
    old_warning = cmp$left_warning, new_warning = cmp$right_warning,
    detail = detail
  )
  cmp
}

record_null_order <- function(check, mode, shape, n, nsim, seed, threads,
                              chunk, cmp) {
  if (!isTRUE(cmp$pass) && is.null(cmp$result$left$sim_K)) return(invisible(NULL))
  left <- cmp$result$left
  right <- cmp$result$right
  if (is.null(left$sim_K) || is.null(right$sim_K)) return(invisible(NULL))
  if (!identical(dim(left$sim_K), dim(right$sim_K))) return(invisible(NULL))
  observed <- as.numeric(left$K)
  for (i in seq_len(nsim)) {
    old_null <- left$sim_K[i, 1L]
    new_null <- right$sim_K[i, 1L]
    controlled_null_rows[[length(controlled_null_rows) + 1L]] <<- data.frame(
      check = check, mode = mode, shape = shape, n = n, nsim = nsim,
      seed = if (is.null(seed)) NA_integer_ else seed,
      ncores = threads, simulation_chunk = chunk, replicate = i,
      old_null_K = old_null, new_null_K = new_null,
      null_exact = identical(old_null, new_null),
      old_exceedance = isTRUE(old_null >= observed),
      new_exceedance = isTRUE(new_null >= as.numeric(right$K)[[1L]]),
      stringsAsFactors = FALSE
    )
  }
  invisible(NULL)
}

# Direct fixed-permutation oracle.  The same integer matrix is passed to both
# evaluators, so no RNG draw is involved and each replicate can be audited.
for (shape in shapes) for (n in n_grid) {
  tree <- stage2d1_make_tree(n, shape)
  compiled_tree <- stage2d1_compile_tree(tree, ns)
  X <- stage2d1_make_matrix(tree, 1L)
  for (nsim in nsim_grid) {
    perms <- stage2d1_make_permutations(n, nsim, seed = 610000L + n + nsim)
    result <- compare_pair(
      "controlled", "A/B/C-oracle", shape, n, nsim, NA_integer_, 1L, 1L,
      X, perms, FALSE, old_engine, hooks$controlled
    )
    cmp <- record_pair("controlled_permutation_exact", "controlled", "A/B/C",
                       shape, n, nsim, NA_integer_, 1L, 1L, result)
    record_null_order("controlled_null_order", "controlled", shape, n, nsim,
                      NA_integer_, 1L, 1L, cmp)
    if (length(ncores) > 1L) {
      result_mt <- compare_pair(
        "controlled", "A/B/C-oracle", shape, n, nsim, NA_integer_, 2L, 1L,
        X, perms, FALSE, old_engine, hooks$controlled
      )
      record_pair("controlled_thread_replay", "controlled", "A/B/C", shape,
                  n, nsim, NA_integer_, 2L, 1L, result_mt)
    }
  }
}

# Internal RNG replay.  The comparison is deliberately made at the same
# seed, ncores, and chunk.  Chunk values include the current 128 boundary and
# its neighbours plus the requested-nsim boundary cases.
chunks_for <- function(nsim) {
  unique(pmax(1L, c(1L, 2L, 127L, 128L, 129L, 255L, 256L, 257L,
                    nsim - 1L, nsim, nsim + 1L)))
}
old_rng_cache <- list()
new_rng_cache <- list()
for (shape in shapes) for (n in n_grid) {
  tree <- stage2d1_make_tree(n, shape)
  compiled_tree <- stage2d1_compile_tree(tree, ns)
  X <- stage2d1_make_matrix(tree, 1L)
  for (nsim in nsim_grid) for (seed in seeds) for (threads in ncores)
    for (chunk in chunks_for(nsim)) {
      result <- compare_pair(
        "rng", "A/B", shape, n, nsim, seed, threads, chunk, X, NULL, TRUE,
        old_engine, hooks$rng
      )
      cmp <- record_pair("same_seed_same_ncores_replay", "rng", "A/B", shape,
                         n, nsim, seed, threads, chunk, result)
      record_null_order("rng_null_order", "rng", shape, n, nsim, seed, threads,
                        chunk, cmp)
      key <- paste(shape, n, nsim, seed, threads, sep = "|")
      if (isTRUE(result$old$status == "ok")) {
        old_rng_cache[[paste(key, chunk, sep = "|")]] <-
          stage2d1_normalize_result(result$old$value, nsim, ncol(X))
      }
      if (isTRUE(result$new$status == "ok")) {
        new_rng_cache[[paste(key, chunk, sep = "|")]] <-
          stage2d1_normalize_result(result$new$value, nsim, ncol(X))
      }
    }
}

check_chunk_replay <- function(cache, label) {
  keys <- unique(sub("\\|[^|]+$", "", names(cache)))
  for (key in keys) {
    entries <- cache[startsWith(names(cache), paste0(key, "|"))]
    if (length(entries) < 2L) next
    first <- entries[[1L]]
    pass <- vapply(entries[-1L], function(z) {
      stage2d1_num_equal(first$sim_K, z$sim_K, 0) &&
        stage2d1_num_equal(first$P, z$P, 0) &&
        stage2d1_num_equal(first$MCSE_P, z$MCSE_P, 0) &&
        stage2d1_num_equal(first$exceedance, z$exceedance, 0)
    }, logical(1))
    parts <- strsplit(key, "\\|", fixed = FALSE)[[1L]]
    add_row(
      check = paste0(label, "_chunk_boundary_replay"), mode = "rng",
      candidate = label, shape = parts[[1L]], n = as.integer(parts[[2L]]),
      nsim = as.integer(parts[[3L]]), seed = as.integer(parts[[4L]]),
      ncores = as.integer(parts[[5L]]), simulation_chunk = NA_integer_,
      pass = all(pass), old_status = "ok", new_status = "ok",
      old_warning = "", new_warning = "",
      detail = paste0("compared ", length(entries), " chunk settings")
    )
  }
}
check_chunk_replay(old_rng_cache, "old_production")
check_chunk_replay(new_rng_cache, "prototype")

# Explicit malformed-input failures exercise validation/accounting without
# weakening the candidate's error contract.
for (shape in shapes) {
  n <- min(n_grid)
  tree <- stage2d1_make_tree(n, shape)
  compiled_tree <- stage2d1_compile_tree(tree, ns)
  X <- stage2d1_make_matrix(tree, 1L)
  perms <- stage2d1_make_permutations(n, 2L, 611L)
  perms[1L, 2L] <- perms[1L, 1L]
  old <- stage2d1_capture(function() stage2d1_call_engine(
    old_engine, "controlled", compiled_tree, X, 2L, perms, trait_chunk,
    TRUE, FALSE, 1L, 1L, NULL
  ))
  new <- stage2d1_capture(function() stage2d1_call_engine(
    hooks$controlled, "controlled", compiled_tree, X, 2L, perms, trait_chunk,
    TRUE, FALSE, 1L, 1L, NULL
  ))
  pass <- old$status == "error" && new$status == "error" &&
    grepl("permutation|index|once", old$error_message, ignore.case = TRUE) &&
    grepl("permutation|index|once", new$error_message, ignore.case = TRUE)
  failure_rows[[length(failure_rows) + 1L]] <- data.frame(
    shape = shape, failure = "duplicate permutation index", pass = pass,
    old_status = old$status, new_status = new$status,
    old_error = old$error_message, new_error = new$error_message,
    stringsAsFactors = FALSE
  )
}

# Scientific edge fixtures: the ordered null vector and all accounting fields
# are compared exactly for ordinary, offset, constant, near-constant, and
# two-tip traits.  Constant traits are expected to retain the production
# non-finite/NA semantics, not to be coerced into a finite P value.
edge_specs <- list(
  ordinary = function(tree) stage2d1_make_trait(tree),
  large_offset = function(tree) stage2d1_make_trait(tree, 700000L) + 1e14,
  constant = function(tree) stats::setNames(rep(3, length(tree$tip.label)), tree$tip.label),
  near_constant = function(tree) {
    y <- stats::setNames(rep(3, length(tree$tip.label)), tree$tip.label)
    y[[1L]] <- y[[1L]] + 1e-8
    y
  }
)
for (shape in shapes) for (spec_name in names(edge_specs)) {
  n <- if (isTRUE(args$smoke)) 8L else max(50L, min(n_grid))
  tree <- stage2d1_make_tree(n, shape)
  compiled_tree <- stage2d1_compile_tree(tree, ns)
  X <- matrix(edge_specs[[spec_name]](tree), ncol = 1L,
              dimnames = list(tree$tip.label, spec_name))
  perms <- stage2d1_make_permutations(n, 2L, 612L)
  result <- compare_pair("controlled", spec_name, shape, n, 2L, NA_integer_,
                         1L, 1L, X, perms, FALSE, old_engine,
                         hooks$controlled)
  record_pair(paste0("edge_", spec_name), "controlled", spec_name, shape, n,
              2L, NA_integer_, 1L, 1L, result)
}
tree_two <- stage2d1_make_tree(2L, "balanced")
compiled_tree <- stage2d1_compile_tree(tree_two, ns)
X_two <- matrix(stage2d1_make_trait(tree_two), ncol = 1L,
                dimnames = list(tree_two$tip.label, "two_tip"))
perms_two <- stage2d1_make_permutations(2L, 2L, 613L)
two_result <- compare_pair("controlled", "two_tip", "balanced", 2L, 2L,
                           NA_integer_, 1L, 1L, X_two, perms_two, FALSE,
                           old_engine, hooks$controlled)
record_pair("edge_two_tip", "controlled", "two_tip", "balanced", 2L, 2L,
            NA_integer_, 1L, 1L, two_result)

# Public NA and input-immutability guard.  Direct old/new equality on each
# retained subset is paired with a public call so removed tips cannot return
# through the candidate path.
na_kinds <- if (isTRUE(args$smoke)) "shared" else
  c("none", "shared", "repeated", "near_unique")
for (kind in na_kinds) {
  n <- if (isTRUE(args$smoke)) 12L else max(50L, min(n_grid))
  tree <- stage2d1_make_tree(n, "balanced")
  X <- stage2d1_make_na_matrix(tree, kind, traits = 4L)
  before_tree <- serialize(tree, NULL, version = 3L)
  before_X <- serialize(X, NULL, version = 3L)
  public <- stage2d1_capture(function() fastphylosig::fast_k(
    tree, X, test = TRUE, nsim = if (isTRUE(args$smoke)) 5L else 19L,
    verbose = FALSE, progress = FALSE, ncores = 1L
  ))
  immutable <- identical(serialize(tree, NULL, version = 3L), before_tree) &&
    identical(serialize(X, NULL, version = 3L), before_X)
  na_rows[[length(na_rows) + 1L]] <- data.frame(
    mask = kind, n = n, public_status = public$status,
    warning_text = public$warning_text, input_immutable = immutable,
    retained_values_only = public$status == "ok", pass = immutable &&
      public$status == "ok", detail = "public K handled mask without input mutation",
    stringsAsFactors = FALSE
  )
  # One direct retained-subtree audit per mask type proves the prototype sees
  # exactly the same legal analysis set as production.
  j <- 1L
  keep <- is.finite(X[, j])
  retained_tree <- ape::drop.tip(tree, tree$tip.label[!keep])
  compiled_tree <- stage2d1_compile_tree(retained_tree, ns)
  retained_X <- matrix(X[keep, j], ncol = 1L,
                       dimnames = list(retained_tree$tip.label, "trait1"))
  perms <- stage2d1_make_permutations(sum(keep), 2L, 614L)
  result <- compare_pair("controlled", paste0("NA_", kind), shape = "balanced",
                         n = sum(keep), nsim = 2L, seed = NA_integer_, threads = 1L,
                         chunk = 1L, X = retained_X, permutations = perms,
                         include_observed = FALSE, old_fun = old_engine,
                         new_fun = hooks$controlled)
  record_pair(paste0("NA_retained_", kind), "controlled", paste0("NA_", kind),
              "balanced", sum(keep), 2L, NA_integer_, 1L, 1L, result)
}

# test=FALSE and the required batch trait sizes are a protection baseline.  A
# candidate batch hook is optional because Stage 2D1 changes only the null
# engine; when absent this remains a non-comparative production guard.
batch_traits <- stage2d1_grid("FASTPHYLOSIG_STAGE2D1_TRAITS",
                              c(1L, 8L, 32L, 100L), c(1L, 8L), minimum = 1L)
for (p in batch_traits) {
  n <- if (isTRUE(args$smoke)) 50L else 500L
  tree <- stage2d1_make_tree(n, "balanced")
  X <- stage2d1_make_matrix(tree, p)
  before_tree <- serialize(tree, NULL, version = 3L)
  before_X <- serialize(X, NULL, version = 3L)
  public <- stage2d1_capture(function() fastphylosig::fast_k(
    tree, X, test = FALSE, verbose = FALSE, progress = FALSE
  ))
  pass <- public$status == "ok" && !nzchar(public$warning_text) &&
    identical(serialize(tree, NULL, version = 3L), before_tree) &&
    identical(serialize(X, NULL, version = 3L), before_X)
  if (is.function(hooks$batch) && public$status == "ok") {
    compiled_tree <- stage2d1_compile_tree(tree, ns)
    candidate <- stage2d1_capture(function() do.call(hooks$batch, list(
      compiled_tree = compiled_tree, X = X, trait_chunk = trait_chunk,
      n_threads = 1L
    )))
    pass <- pass && candidate$status == "ok"
    candidate_status <- candidate$status
  } else candidate_status <- "NOT_RUN"
  batch_rows[[length(batch_rows) + 1L]] <- data.frame(
    traits = p, n = n, public_status = public$status,
    candidate_status = candidate_status, warning_text = public$warning_text,
    pass = pass, comparison = if (is.function(hooks$batch)) "hook" else
      "production_guard_only", stringsAsFactors = FALSE
  )
}

failure <- if (length(failure_rows)) do.call(rbind, failure_rows) else data.frame()
na_guard <- if (length(na_rows)) do.call(rbind, na_rows) else data.frame()
batch_guard <- if (length(batch_rows)) do.call(rbind, batch_rows) else data.frame()
null_order <- if (length(controlled_null_rows)) {
  do.call(rbind, controlled_null_rows)
} else data.frame()
all_rows <- if (length(rows)) do.call(rbind, rows) else data.frame()

stage2d1_write_csv(all_rows, file.path(args$out, "stage2d1_correctness_rows.csv"))
stage2d1_write_csv(null_order, file.path(args$out, "stage2d1_ordered_nulls.csv"))
stage2d1_write_csv(failure, file.path(args$out, "stage2d1_failure_accounting.csv"))
stage2d1_write_csv(na_guard, file.path(args$out, "stage2d1_na_mask_guard.csv"))
stage2d1_write_csv(batch_guard, file.path(args$out, "stage2d1_batch_guard.csv"))

core_pass <- nrow(all_rows) > 0L && all(all_rows$pass) &&
  (nrow(null_order) == 0L || all(null_order$null_exact)) &&
  (nrow(failure) == 0L || all(failure$pass)) &&
  (nrow(na_guard) == 0L || all(na_guard$pass)) &&
  (nrow(batch_guard) == 0L || all(batch_guard$pass))
final_status <- if (core_pass) "PASS" else "FAIL"
stage2d1_write_csv(data.frame(
  status = final_status,
  core_gates = if (core_pass) "PASS" else "FAIL",
  controlled_permutation_exact = if (nrow(all_rows)) all(all_rows$pass) else FALSE,
  ordered_nulls_exact = if (nrow(null_order)) all(null_order$null_exact) else TRUE,
  failure_accounting = if (nrow(failure)) all(failure$pass) else TRUE,
  NA_masks = if (nrow(na_guard)) all(na_guard$pass) else TRUE,
  batch_trait_guard = if (nrow(batch_guard)) all(batch_guard$pass) else TRUE,
  production_code_changed = "NO",
  formal_benchmark_authorized = core_pass,
  detail = if (core_pass) "Correctness/RNG gates passed; benchmark may be run explicitly." else
    "Do not run formal benchmark; inspect failed gate rows.",
  stringsAsFactors = FALSE
), status_path)

if (!core_pass) stop("Stage 2D1 correctness gate failed.", call. = FALSE)
message("[stage2d1] correctness PASS; formal benchmark is separately opt-in")
