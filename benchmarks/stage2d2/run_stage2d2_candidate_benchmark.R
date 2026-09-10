#!/usr/bin/env Rscript

# Stage 2D2 private candidate benchmark.
#
# This is audit infrastructure only.  The oracle and candidate are loaded
# from FASTPHYLOSIG_STAGE2D2_PROTO_CPP in one sourceCpp environment.  The
# runner never changes package code, package bindings, prototypes, or tests.
# Formal mode is deliberately explicit because the largest cells are costly.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- all_args[grepl("^--file=", all_args)]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else
  file.path(getwd(), "run_stage2d2_candidate_benchmark.R")
script_dir <- dirname(normalizePath(script_file, winslash = "/",
                                    mustWork = FALSE))
source(file.path(script_dir, "stage2d2_common.R"), local = TRUE)

parse_cli <- function(values) {
  formal <- "--formal" %in% values
  smoke <- "--smoke" %in% values
  if (formal == smoke) {
    stop("Use exactly one of --formal or --smoke.", call. = FALSE)
  }
  positional <- values[!values %in% c("--formal", "--smoke")]
  if (length(positional) > 2L) {
    stop("usage: run_stage2d2_candidate_benchmark.R [--formal|--smoke] " ,
         "[repo] [out].", call. = FALSE)
  }
  repo <- if (length(positional)) positional[[1L]] else getwd()
  out <- if (length(positional) >= 2L) positional[[2L]] else
    file.path(repo, "benchmarks", "stage2d2", "results",
              "candidate-benchmark")
  list(
    formal = isTRUE(formal), smoke = isTRUE(smoke),
    repo = normalizePath(repo, winslash = "/", mustWork = FALSE),
    out = normalizePath(out, winslash = "/", mustWork = FALSE)
  )
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(args$out, "stage2d2_candidate_benchmark_status.csv")

write_early_status <- function(status, reason) {
  stage2d2_write_csv(data.frame(
    status = as.character(status), candidate_mode = NA_character_,
    same_translation_unit = FALSE, source_cpp = NA_character_,
    reason = as.character(reason), production_code_changed = "NO",
    stringsAsFactors = FALSE
  ), status_path)
  stop(reason, call. = FALSE)
}

Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")

env_or <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

parse_grid <- function(name, formal_default, smoke_default, minimum = 1L) {
  value <- env_or(name, "")
  if (nzchar(value)) {
    return(stage2d2_ints(value, name, minimum = minimum))
  }
  if (isTRUE(args$smoke)) smoke_default else formal_default
}

parse_shapes <- function() {
  value <- env_or("FASTPHYLOSIG_STAGE2D2_BENCH_SHAPES", "")
  out <- if (nzchar(value))
    unique(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])) else
    if (isTRUE(args$smoke)) "balanced" else
      c("balanced", "random", "pectinate")
  if (!length(out) || any(!nzchar(out)) ||
      any(!out %in% c("balanced", "random", "pectinate"))) {
    stop("FASTPHYLOSIG_STAGE2D2_BENCH_SHAPES must contain only balanced, ",
         "random, or pectinate.", call. = FALSE)
  }
  out
}

parse_scalar_int <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.integer(env_or(name, as.character(default))))
  if (is.na(value) || value < minimum) {
    stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
  }
  value
}

shapes <- parse_shapes()
n_grid <- parse_grid(
  "FASTPHYLOSIG_STAGE2D2_BENCH_N",
  c(5000L, 10000L), c(50L), minimum = 2L
)
nsim_grid <- parse_grid(
  "FASTPHYLOSIG_STAGE2D2_BENCH_NSIM",
  c(9999L), c(19L), minimum = 1L
)
threads_grid <- parse_grid(
  "FASTPHYLOSIG_STAGE2D2_BENCH_THREADS",
  c(1L, 2L), c(1L), minimum = 1L
)
block_grid <- parse_grid(
  "FASTPHYLOSIG_STAGE2D2_BENCH_BLOCKS",
  c(2L, 4L, 8L, 16L, 32L), c(2L, 4L), minimum = 1L
)
if (any(!block_grid %in% c(2L, 4L, 8L, 16L, 32L))) {
  stop("FASTPHYLOSIG_STAGE2D2_BENCH_BLOCKS must use 2, 4, 8, 16, or 32.",
       call. = FALSE)
}
trait_chunk <- parse_scalar_int(
  "FASTPHYLOSIG_STAGE2D2_BENCH_TRAIT_CHUNK", 64L, minimum = 1L
)
simulation_chunk <- parse_scalar_int(
  "FASTPHYLOSIG_STAGE2D2_BENCH_SIMULATION_CHUNK", 128L, minimum = 1L
)
guard_n <- parse_scalar_int(
  "FASTPHYLOSIG_STAGE2D2_BENCH_GUARD_N",
  if (isTRUE(args$smoke)) 50L else 500L, minimum = 2L
)
light_nsim <- parse_scalar_int(
  "FASTPHYLOSIG_STAGE2D2_BENCH_LIGHT_NSIM", 199L, minimum = 1L
)
batch_traits <- parse_grid(
  "FASTPHYLOSIG_STAGE2D2_BENCH_TRAITS",
  c(1L, 8L, 32L, 100L), c(1L, 8L, 32L, 100L), minimum = 1L
)
guard_threads <- parse_grid(
  "FASTPHYLOSIG_STAGE2D2_BENCH_GUARD_THREADS", c(1L), c(1L), minimum = 1L
)
repeats <- parse_scalar_int(
  "FASTPHYLOSIG_STAGE2D2_BENCH_REPEATS",
  if (isTRUE(args$smoke)) 2L else 10L, minimum = 1L
)
if (isTRUE(args$formal) && repeats < 10L) {
  stop("Formal heavy timing requires at least 10 paired repeats.",
       call. = FALSE)
}
light_repeats <- parse_scalar_int(
  "FASTPHYLOSIG_STAGE2D2_BENCH_LIGHT_REPEATS",
  if (isTRUE(args$smoke)) 1L else min(3L, repeats), minimum = 1L
)
guard_repeats <- parse_scalar_int(
  "FASTPHYLOSIG_STAGE2D2_BENCH_TEST_FALSE_REPEATS",
  if (isTRUE(args$smoke)) 1L else min(3L, repeats), minimum = 1L
)
base_seed <- parse_scalar_int(
  "FASTPHYLOSIG_STAGE2D2_BENCH_SEED", 20260910L, minimum = 1L
)

# The sourceCpp file is the only permitted private implementation input.  A
# single environment is used for every export so a claimed oracle/candidate
# pair cannot silently come from two separately compiled translation units.
find_export <- function(env, candidates) {
  for (name in candidates) {
    if (exists(name, envir = env, inherits = FALSE)) {
      value <- get(name, envir = env, inherits = FALSE)
      if (is.function(value)) return(list(name = name, fun = value))
    }
  }
  NULL
}

make_mode_wrapper <- function(target, mode) {
  force(target)
  force(mode)
  function(...) {
    call_args <- list(...)
    call_args$mode <- mode
    do.call(target, stage2d2_filter_call_args(target, call_args))
  }
}

load_candidate_hooks <- function() {
  cpp_file <- env_or("FASTPHYLOSIG_STAGE2D2_PROTO_CPP", "")
  if (!nzchar(cpp_file)) {
    stop("Set FASTPHYLOSIG_STAGE2D2_PROTO_CPP to one audit translation unit.",
         call. = FALSE)
  }
  cpp_file <- normalizePath(cpp_file, winslash = "/", mustWork = TRUE)
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Rcpp is required for FASTPHYLOSIG_STAGE2D2_PROTO_CPP.",
         call. = FALSE)
  }
  hook_env <- new.env(parent = globalenv())
  before <- ls(envir = hook_env, all.names = TRUE)
  Rcpp::sourceCpp(cpp_file, env = hook_env, rebuild = TRUE,
                  showOutput = FALSE, verbose = FALSE)
  exports <- setdiff(ls(envir = hook_env, all.names = TRUE), before)
  cpp_bound <- function(entry) {
    !is.null(entry) && entry$name %in% exports &&
      exists(entry$name, envir = hook_env, inherits = FALSE) &&
      identical(get(entry$name, envir = hook_env, inherits = FALSE), entry$fun)
  }

  generic <- find_export(hook_env, c(
    "stage2d2_k_permutation_prototype", "stage2d2_k_null_prototype"
  ))
  oracle_raw <- find_export(hook_env, c(
    "stage2d2_k_permutation_oracle", "stage2d2_k_null_oracle", "oracle"
  ))
  candidate_a_raw <- find_export(hook_env, c(
    "stage2d2_k_permutation_candidate_a", "stage2d2_k_null_candidate_a",
    "candidate_a"
  ))
  candidate_b_raw <- find_export(hook_env, c(
    "stage2d2_k_permutation_candidate_b", "stage2d2_k_null_candidate_b",
    "candidate_b"
  ))

  requested <- toupper(trimws(env_or(
    "FASTPHYLOSIG_STAGE2D2_CANDIDATE_MODE", "AUTO"
  )))
  if (!requested %in% c("AUTO", "A", "B", "C")) {
    stop("FASTPHYLOSIG_STAGE2D2_CANDIDATE_MODE must be AUTO, A, B, or C.",
         call. = FALSE)
  }
  mode <- requested
  if (mode == "AUTO") {
    mode <- if (!is.null(generic)) "A" else if (
      !is.null(oracle_raw) && !is.null(candidate_b_raw)
    ) "B" else if (!is.null(oracle_raw) && !is.null(candidate_a_raw)) "A" else
      "UNAVAILABLE"
  }
  if (mode == "UNAVAILABLE") {
    stop("The sourceCpp TU exposed neither generic mode A nor separate B ",
         "oracle/candidate exports.", call. = FALSE)
  }
  b_gate <- toupper(trimws(env_or(
    "FASTPHYLOSIG_STAGE2D2_CANDIDATE_B_GATE", "NOT_AUTHORIZED"
  )))
  if (mode == "B" && !b_gate %in% c("GO", "AUTHORIZED", "RUN", "PASS")) {
    stop("Candidate B timing is disabled until its exactness gate is ",
         "authorized (set FASTPHYLOSIG_STAGE2D2_CANDIDATE_B_GATE).",
         call. = FALSE)
  }

  if (mode %in% c("A", "C")) {
    if (!is.null(generic)) {
      if (!cpp_bound(generic)) stop("Generic export is not sourceCpp-bound.",
                                    call. = FALSE)
      oracle <- make_mode_wrapper(generic$fun, "oracle")
      candidate <- make_mode_wrapper(generic$fun, mode)
      oracle_name <- paste0(generic$name, "(mode=oracle)")
      candidate_name <- paste0(generic$name, "(mode=", mode, ")")
      same_tu <- cpp_bound(generic)
      layout <- "generic_mode"
    } else {
      if (mode != "A" || is.null(oracle_raw) || is.null(candidate_a_raw) ||
          !cpp_bound(oracle_raw) || !cpp_bound(candidate_a_raw)) {
        stop("Candidate A requires a sourceCpp-bound generic export or ",
             "same-TU oracle/candidate_a exports.", call. = FALSE)
      }
      oracle <- oracle_raw$fun
      candidate <- candidate_a_raw$fun
      oracle_name <- oracle_raw$name
      candidate_name <- candidate_a_raw$name
      same_tu <- cpp_bound(oracle_raw) && cpp_bound(candidate_a_raw)
      layout <- "separate_exports"
    }
  } else {
    if (is.null(oracle_raw) || is.null(candidate_b_raw) ||
        !cpp_bound(oracle_raw) || !cpp_bound(candidate_b_raw)) {
      stop("Candidate B requires sourceCpp-bound oracle and candidate_b ",
           "exports from the same TU.", call. = FALSE)
    }
    oracle <- oracle_raw$fun
    candidate <- candidate_b_raw$fun
    oracle_name <- oracle_raw$name
    candidate_name <- candidate_b_raw$name
    same_tu <- cpp_bound(oracle_raw) && cpp_bound(candidate_b_raw)
    layout <- "separate_exports"
  }
  if (!isTRUE(same_tu)) stop("The oracle/candidate pair is not same-TU bound.",
                             call. = FALSE)
  list(
    cpp_file = cpp_file, env = hook_env, exports = exports,
    mode = mode, layout = layout, oracle = oracle, candidate = candidate,
    oracle_name = oracle_name, candidate_name = candidate_name,
    generic_name = if (is.null(generic)) "" else generic$name,
    candidate_b_gate = b_gate,
    same_translation_unit = same_tu
  )
}

hooks <- tryCatch(load_candidate_hooks(), error = function(e) {
  write_early_status("NOT_RUN", conditionMessage(e))
})

ns <- tryCatch(stage2d2_load_package(args$repo), error = function(e) {
  write_early_status("NOT_RUN", paste0("Cannot load fastphylosig: ",
                                        conditionMessage(e)))
})
if (!is.environment(ns)) write_early_status("NOT_RUN",
                                             "fastphylosig namespace unavailable.")

compile_tree <- function(tree) stage2d2_compile_tree(tree, ns)

source_sha256 <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256"))
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    return(paste(as.character(openssl::sha256(file(path))), collapse = ""))
  }
  NA_character_
}

git_value <- function(repo) {
  value <- tryCatch(system2("git", c("-C", repo, "rev-parse", "HEAD"),
                            stdout = TRUE, stderr = FALSE),
                    error = function(e) character())
  if (length(value)) paste(value, collapse = "\n") else "UNAVAILABLE"
}

fixture <- function(n, shape, traits, seed) {
  tree <- stage2d2_make_tree(n, shape)
  X <- stage2d2_make_matrix(tree, traits = traits, kind = "normal", seed = seed)
  list(tree = tree, compiled_tree = compile_tree(tree), X = X,
       n = as.integer(n), shape = shape, traits = as.integer(traits))
}

trait_chunk_for <- function(p) as.integer(min(trait_chunk, max(1L, p)))

# Use a reproducible, bounded integer seed without overflowing R's signed
# integer representation when a user supplies a large grid.
cell_seed <- function(stream, n, nsim, traits, threads, block, pair = 0L) {
  block_value <- if (is.na(block)) 0 else as.double(block)
  raw <- as.double(base_seed) + as.double(stream) * 1000003 +
    as.double(n) * 1009 + as.double(nsim) * 9176 +
    as.double(traits) * 101 + as.double(threads) * 17 +
    block_value * 13 + as.double(pair)
  as.integer((raw %% 2147483646) + 1)
}

call_native <- function(fun, compiled_tree, X, nsim, n_threads, block_size,
                        seed, return_sim = FALSE, include_observed = TRUE,
                        permutations = NULL, return_ordered = FALSE,
                        collect_counters = FALSE) {
  call_args <- list(
    compiled_tree = compiled_tree, X = X, nsim = as.integer(nsim),
    permutations = permutations, trait_chunk = trait_chunk_for(ncol(X)),
    return_sim = isTRUE(return_sim), include_observed = isTRUE(include_observed),
    n_threads = as.integer(n_threads), simulation_chunk = simulation_chunk,
    return_ordered = isTRUE(return_ordered)
  )
  if (!is.na(block_size)) call_args$block_size <- as.integer(block_size)
  if (isTRUE(collect_counters)) call_args$collect_counters <- TRUE
  stage2d2_capture(function() stage2d2_call_with_seed(fun, call_args, seed))
}

one_value <- function(value, default = NA) {
  if (is.null(value) || !length(value)) return(default)
  value[[1L]]
}

pick_value <- function(value, names, default = NA) {
  if (!is.list(value)) return(default)
  for (name in names) if (!is.null(value[[name]])) return(value[[name]])
  default
}

as_number <- function(value, default = NA_real_) {
  one <- one_value(value, default)
  out <- suppressWarnings(as.numeric(one))
  if (!length(out) || is.na(out)) default else out[[1L]]
}

safe_object_size <- function(value) {
  if (is.null(value)) return(NA_real_)
  tryCatch(as.numeric(utils::object.size(value)), error = function(e) NA_real_)
}

memory_row <- function(captured, implementation, cfg, phase, timed) {
  value <- captured$value
  memory <- if (is.list(value)) value$memory else NULL
  candidate_metadata <- if (is.list(value)) value$candidate_metadata else NULL
  reported <- function(names, default = NA_real_) {
    top <- pick_value(value, names, NULL)
    nested <- pick_value(memory, names, NULL)
    as_number(if (!is.null(top)) top else nested, default)
  }
  text_reported <- function(names, default = NA_character_) {
    top <- pick_value(value, names, NULL)
    nested <- pick_value(memory, names, NULL)
    out <- if (!is.null(top)) top else nested
    if (is.null(out) || !length(out)) default else as.character(out[[1L]])
  }
  candidate_reported <- function(names, default = NA_real_) {
    as_number(pick_value(candidate_metadata, names, NULL), default)
  }
  metadata_present <- is.list(value) && (is.list(memory) ||
    any(c("stage2d2_workspace_bytes", "stage2d2_workspace_peak_bytes",
          "peak_workspace_bytes", "peak_index_bytes") %in% names(value)))
  data.frame(
    stage = "Stage 2D2 candidate benchmark",
    route = cfg$route, phase = phase, timed = isTRUE(timed),
    implementation = implementation, candidate = hooks$mode,
    shape = cfg$shape, n = as.integer(cfg$n), nsim = as.integer(cfg$nsim),
    traits = as.integer(cfg$traits), n_threads = as.integer(cfg$n_threads),
    block_size = as.integer(cfg$block_size), pair = as.integer(cfg$pair),
    seed = as.integer(cfg$seed), elapsed_s = captured$elapsed_s,
    status = captured$status, metadata_status = if (metadata_present)
      "reported" else if (captured$status == "ok") "not_reported" else "error",
    workspace_bytes = reported(c("stage2d2_workspace_bytes", "workspace_bytes")),
    workspace_peak_bytes = reported(c(
      "stage2d2_workspace_peak_bytes", "workspace_peak_bytes"
    )),
    workspace_allocations = reported(c("stage2d2_workspace_allocations",
                                       "workspace_allocations")),
    workspace_reuses = reported(c("stage2d2_workspace_reuses",
                                  "workspace_reuses")),
    peak_index_bytes = reported(c("peak_index_bytes")),
    peak_workspace_bytes = reported(c("peak_workspace_bytes")),
    requested_sim_output_bytes = reported(c("requested_sim_output_bytes")),
    blocks_processed = candidate_reported(c("blocks_processed")),
    generated_permutations = candidate_reported(c("generated_permutations")),
    fisher_yates_swaps = candidate_reported(c("fisher_yates_swaps")),
    compute_one_equivalents = candidate_reported(c("compute_one_equivalents")),
    memory_bounded = text_reported(c("bounded"), NA_character_),
    working_memory_model = text_reported(c("working_memory_model"),
                                         NA_character_),
    memory_note = text_reported(c("note"), NA_character_),
    result_object_bytes = safe_object_size(value),
    warning_text = captured$warning_text, error_message = captured$error_message,
    stringsAsFactors = FALSE
  )
}

timing_rows <- list()
memory_rows <- list()
parity_rows <- list()

record_capture <- function(captured, implementation, cfg, phase, timed,
                           order) {
  timing_rows[[length(timing_rows) + 1L]] <<- data.frame(
    stage = "Stage 2D2 candidate benchmark", route = cfg$route,
    phase = phase, timed = isTRUE(timed), implementation = implementation,
    candidate = hooks$mode, shape = cfg$shape, n = as.integer(cfg$n),
    nsim = as.integer(cfg$nsim), traits = as.integer(cfg$traits),
    n_threads = as.integer(cfg$n_threads), block_size = as.integer(cfg$block_size),
    trait_chunk = trait_chunk_for(cfg$traits),
    simulation_chunk = as.integer(simulation_chunk), pair = as.integer(cfg$pair),
    seed = as.integer(cfg$seed), order = as.character(order),
    elapsed_s = captured$elapsed_s, status = captured$status,
    warning_text = captured$warning_text,
    warning_classes = captured$warning_classes,
    error_class = captured$error_class,
    error_message = captured$error_message,
    result_status = if (is.list(captured$value)) as.character(
      one_value(captured$value$status, NA_character_)
    ) else NA_character_,
    stringsAsFactors = FALSE
  )
  memory_rows[[length(memory_rows) + 1L]] <<- memory_row(
    captured, implementation, cfg, phase, timed
  )
  invisible(captured)
}

blocks_for_mode <- function() if (identical(hooks$mode, "B")) block_grid else
  NA_integer_

run_timed_cell <- function(route, f, nsim, n_threads, block_size, repeats,
                           stream) {
  p <- ncol(f$X)
  cfg <- list(route = route, shape = f$shape, n = f$n, nsim = nsim,
             traits = p, n_threads = n_threads, block_size = block_size,
             pair = 0L, seed = cell_seed(stream, f$n, nsim, p, n_threads,
                                          block_size, 0L))
  warm_seed <- cfg$seed
  warm_oracle <- call_native(
    hooks$oracle, f$compiled_tree, f$X, nsim, n_threads, block_size,
    warm_seed, return_sim = FALSE, include_observed = TRUE
  )
  warm_candidate <- call_native(
    hooks$candidate, f$compiled_tree, f$X, nsim, n_threads, block_size,
    warm_seed, return_sim = FALSE, include_observed = TRUE
  )
  record_capture(warm_oracle, "oracle", cfg, "warmup", FALSE, "warmup")
  record_capture(warm_candidate, "candidate", cfg, "warmup", FALSE, "warmup")
  if (warm_oracle$status != "ok" || warm_candidate$status != "ok") {
    stop("Warmup failed for ", route, "/", f$shape, "/n=", f$n,
         "/nsim=", nsim, "/threads=", n_threads, "/block=", block_size,
         ": ", warm_oracle$error_message, " / ",
         warm_candidate$error_message, call. = FALSE)
  }

  for (pair in seq_len(repeats)) {
    pair_seed <- cell_seed(stream, f$n, nsim, p, n_threads, block_size, pair)
    oracle_first <- pair %% 2L == 1L
    first_fun <- if (oracle_first) hooks$oracle else hooks$candidate
    second_fun <- if (oracle_first) hooks$candidate else hooks$oracle
    first_name <- if (oracle_first) "oracle" else "candidate"
    second_name <- if (oracle_first) "candidate" else "oracle"
    order <- paste(first_name, second_name, sep = ">")
    cfg$pair <- pair
    cfg$seed <- pair_seed
    first <- call_native(first_fun, f$compiled_tree, f$X, nsim, n_threads,
                         block_size, pair_seed, return_sim = FALSE,
                         include_observed = TRUE)
    second <- call_native(second_fun, f$compiled_tree, f$X, nsim, n_threads,
                          block_size, pair_seed, return_sim = FALSE,
                          include_observed = TRUE)
    record_capture(first, first_name, cfg, "timed", TRUE, order)
    record_capture(second, second_name, cfg, "timed", TRUE, order)
  }
  invisible(NULL)
}

# A small controlled-permutation parity gate catches a mismatched export or a
# candidate that changes result shape before expensive timing begins.
run_parity <- function() {
  n <- max(2L, min(17L, guard_n))
  tree <- stage2d2_make_tree(n, "balanced")
  X <- stage2d2_make_matrix(tree, traits = 1L, kind = "normal", seed = 20261017L)
  compiled <- compile_tree(tree)
  perms <- stage2d2_make_permutations(n, 7L, seed = 20261018L)
  parity_blocks <- blocks_for_mode()
  for (block_size in parity_blocks) {
    seed <- cell_seed(700L, n, 7L, 1L, 1L, block_size, 1L)
    oracle <- call_native(
      hooks$oracle, compiled, X, 7L, 1L, block_size, seed,
      return_sim = TRUE, include_observed = FALSE, permutations = perms,
      return_ordered = TRUE
    )
    candidate <- call_native(
      hooks$candidate, compiled, X, 7L, 1L, block_size, seed,
      return_sim = TRUE, include_observed = FALSE, permutations = perms,
      return_ordered = TRUE
    )
    compared <- stage2d2_compare_captures(oracle, candidate, 7L, 1L)
    max_field <- function(name) {
      value <- compared$result$max_abs
      if (is.null(value) || is.null(value[[name]])) NA_real_ else
        as.numeric(value[[name]])
    }
    parity_rows[[length(parity_rows) + 1L]] <<- data.frame(
      stage = "Stage 2D2 candidate benchmark parity",
      candidate = hooks$mode, shape = "balanced", n = n, nsim = 7L,
      traits = 1L, n_threads = 1L, block_size = as.integer(block_size),
      oracle_status = oracle$status, candidate_status = candidate$status,
      status_pass = compared$status_pass, warning_pass = compared$warning_pass,
      result_pass = compared$result$pass, pass = isTRUE(compared$pass),
      max_abs_K = max_field("K"), max_abs_P = max_field("P"),
      max_abs_MCSE_P = max_field("MCSE_P"),
      max_abs_exceedance = max_field("exceedance"),
      max_abs_sim_K = max_field("sim_K"),
      oracle_warning = oracle$warning_text,
      candidate_warning = candidate$warning_text,
      oracle_error = oracle$error_message,
      candidate_error = candidate$error_message,
      stringsAsFactors = FALSE
    )
  }
  parity <- stage2d2_rbind(parity_rows)
  if (!nrow(parity) || !all(parity$pass)) {
    stage2d2_write_csv(parity, file.path(args$out,
                                        "stage2d2_candidate_benchmark_parity.csv"))
    write_early_status("FAIL", "same-TU oracle/candidate parity gate failed.")
  }
  invisible(NULL)
}

run_parity()

# Formal cells are serial in this process.  No future, fork, or parallel R
# worker is created by the runner; only the requested native n_threads value
# is passed to the private evaluator.
for (shape in shapes) for (n in n_grid) {
  message("[stage2d2] heavy fixture shape=", shape, ", n=", n)
  f <- fixture(n, shape, 1L, seed = cell_seed(101L, n, 1L, 1L, 1L,
                                              NA_integer_))
  for (nsim in nsim_grid) for (n_threads in threads_grid) {
    for (block_size in blocks_for_mode()) {
      run_timed_cell(
        "formal_heavy", f, nsim, n_threads, block_size, repeats, stream = 111L
      )
    }
  }
}

# The light path keeps the requested nsim=199 and exercises all batch trait
# widths.  A smaller guard tree/repeat count can be selected for smoke runs.
for (p in batch_traits) {
  message("[stage2d2] light batch traits=", p, ", n=", guard_n,
          ", nsim=", light_nsim)
  f <- fixture(guard_n, "balanced", p,
               seed = cell_seed(201L, guard_n, light_nsim, p, 1L,
                                NA_integer_))
  for (n_threads in guard_threads) for (block_size in blocks_for_mode()) {
    run_timed_cell(
      "light_nsim199", f, light_nsim, n_threads, block_size,
      light_repeats, stream = 211L
    )
  }
}

# Public test=FALSE is a production guard, recorded but never used as a
# candidate speedup.  It confirms the unchanged batch path for p=1/8/32/100.
run_public_test_false <- function() {
  for (p in batch_traits) {
    f <- fixture(guard_n, "balanced", p,
                 seed = cell_seed(301L, guard_n, 0L, p, 1L, NA_integer_))
    cfg <- list(route = "test_false", shape = "balanced", n = guard_n,
                nsim = NA_integer_, traits = p, n_threads = 1L,
                block_size = NA_integer_, pair = 0L,
                seed = NA_integer_)
    public_call <- function() stage2d2_capture(function() fastphylosig::fast_k(
      f$tree, x = f$X, test = FALSE, ncores = 1L,
      trait_chunk = trait_chunk_for(p), simulation_chunk = simulation_chunk,
      return_sim = FALSE, keep_null = FALSE, verbose = FALSE, progress = FALSE
    ))
    warm <- public_call()
    record_capture(warm, "production", cfg, "warmup", FALSE, "warmup")
    if (warm$status != "ok") {
      stop("test=FALSE warmup failed for traits=", p, ": ",
           warm$error_message, call. = FALSE)
    }
    for (pair in seq_len(guard_repeats)) {
      cfg$pair <- pair
      measured <- public_call()
      record_capture(measured, "production", cfg, "timed", TRUE, "public")
    }
  }
  invisible(NULL)
}
run_public_test_false()

raw <- stage2d2_rbind(timing_rows)
memory <- stage2d2_rbind(memory_rows)
parity <- stage2d2_rbind(parity_rows)

summary_rows <- list()
summary_keys <- c("route", "candidate", "shape", "n", "nsim", "traits",
                  "n_threads", "block_size", "trait_chunk",
                  "simulation_chunk")
timed <- raw[raw$timed & raw$route %in% c("formal_heavy", "light_nsim199"),
             , drop = FALSE]
if (nrow(timed)) {
  group_key <- do.call(paste, c(timed[summary_keys], sep = "\r"))
  groups <- split(seq_len(nrow(timed)), group_key, drop = TRUE)
  for (indices in groups) {
    z <- timed[indices, , drop = FALSE]
    first <- z[z$implementation == "oracle", , drop = FALSE]
    second <- z[z$implementation == "candidate", , drop = FALSE]
    pair_ids <- intersect(first$pair, second$pair)
    pair_rows <- lapply(pair_ids, function(pair) {
      a <- first[first$pair == pair, , drop = FALSE]
      b <- second[second$pair == pair, , drop = FALSE]
      data.frame(pair = pair, oracle_s = a$elapsed_s[[1L]],
                 candidate_s = b$elapsed_s[[1L]],
                 oracle_status = a$status[[1L]],
                 candidate_status = b$status[[1L]], stringsAsFactors = FALSE)
    })
    joined <- if (length(pair_rows)) do.call(rbind, pair_rows) else data.frame()
    ok_pairs <- if (nrow(joined)) joined$oracle_status == "ok" &
      joined$candidate_status == "ok" & is.finite(joined$oracle_s) &
      is.finite(joined$candidate_s) else logical()
    ratios <- if (any(ok_pairs)) joined$oracle_s[ok_pairs] /
      joined$candidate_s[ok_pairs] else numeric()
    oracle_values <- first$elapsed_s[first$status == "ok" &
                                      is.finite(first$elapsed_s)]
    candidate_values <- second$elapsed_s[second$status == "ok" &
                                          is.finite(second$elapsed_s)]
    oracle_median <- if (length(oracle_values)) median(oracle_values) else NA_real_
    candidate_median <- if (length(candidate_values)) median(candidate_values) else
      NA_real_
    speedup <- if (is.finite(oracle_median) && is.finite(candidate_median) &&
                   candidate_median > 0) oracle_median / candidate_median else
      NA_real_
    reduction <- if (is.finite(speedup)) 1 - candidate_median / oracle_median else
      NA_real_
    expected_orders <- ifelse(z$pair %% 2L == 1L, "oracle>candidate",
                              "candidate>oracle")
    order_ok <- all(z$order == expected_orders)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      z[1L, summary_keys, drop = FALSE], pairs = length(pair_ids),
      ok_pairs = sum(ok_pairs), oracle_observations = nrow(first),
      candidate_observations = nrow(second),
      oracle_median_s = oracle_median,
      oracle_IQR_s = if (length(oracle_values)) IQR(oracle_values) else NA_real_,
      candidate_median_s = candidate_median,
      candidate_IQR_s = if (length(candidate_values)) IQR(candidate_values) else
        NA_real_,
      speedup = speedup, reduction_fraction = reduction,
      reduction_percent = if (is.finite(reduction)) 100 * reduction else NA_real_,
      paired_median_speedup = if (length(ratios)) median(ratios) else NA_real_,
      paired_mean_speedup = if (length(ratios)) mean(ratios) else NA_real_,
      paired_min_speedup = if (length(ratios)) min(ratios) else NA_real_,
      paired_max_speedup = if (length(ratios)) max(ratios) else NA_real_,
      alternating_order = order_ok,
      warnings_clear = all(!nzchar(z$warning_text)),
      timing_status = if (length(pair_ids) && all(ok_pairs) && order_ok &&
                           all(!nzchar(z$warning_text))) "PASS" else "FAIL",
      stringsAsFactors = FALSE
    )
  }
}
summary <- stage2d2_rbind(summary_rows)

guards <- raw[raw$route == "test_false" & raw$timed, , drop = FALSE]
if (!nrow(guards)) guards <- data.frame()

stage2d2_write_csv(raw, file.path(args$out,
                                  "stage2d2_candidate_benchmark_raw.csv"))
stage2d2_write_csv(summary, file.path(args$out,
                                     "stage2d2_candidate_benchmark_summary.csv"))
stage2d2_write_csv(memory, file.path(args$out,
                                    "stage2d2_candidate_benchmark_memory.csv"))
stage2d2_write_csv(parity, file.path(args$out,
                                    "stage2d2_candidate_benchmark_parity.csv"))
stage2d2_write_csv(guards, file.path(args$out,
                                    "stage2d2_candidate_benchmark_guards.csv"))

source_cpp_hash <- source_sha256(hooks$cpp_file)
provenance <- data.frame(
  stage = "Stage 2D2 candidate benchmark",
  benchmark_mode = if (isTRUE(args$formal)) "formal" else "smoke",
  candidate_mode = hooks$mode, export_layout = hooks$layout,
  source_cpp = hooks$cpp_file, source_cpp_sha256 = source_cpp_hash,
  source_cpp_exports = paste(hooks$exports, collapse = ","),
  oracle_export = hooks$oracle_name, candidate_export = hooks$candidate_name,
  candidate_b_gate = hooks$candidate_b_gate,
  same_translation_unit = isTRUE(hooks$same_translation_unit),
  source_commit = git_value(args$repo), package_version = tryCatch(
    as.character(utils::packageVersion("fastphylosig")),
    error = function(e) NA_character_
  ), R_version = R.version.string, platform = R.version$platform,
  os = Sys.info()[["sysname"]], machine = Sys.info()[["machine"]],
  n_grid = paste(n_grid, collapse = ","), nsim_grid = paste(nsim_grid,
                                                             collapse = ","),
  shapes = paste(shapes, collapse = ","),
  threads = paste(threads_grid, collapse = ","),
  block_grid = paste(block_grid, collapse = ","), repeats = repeats,
  light_nsim = light_nsim, batch_traits = paste(batch_traits, collapse = ","),
  timing_protocol = "one serialized warmup; paired same-seed calls; alternating order",
  return_sim = FALSE, include_observed = TRUE,
  test_false_route = "public fast_k(test=FALSE), unchanged production path",
  production_code_changed = "NO", stringsAsFactors = FALSE
)
stage2d2_write_csv(provenance, file.path(args$out,
                                         "stage2d2_candidate_benchmark_provenance.csv"))

heavy_summary <- summary[summary$route == "formal_heavy", , drop = FALSE]
light_summary <- summary[summary$route == "light_nsim199", , drop = FALSE]
required_pairs <- if (isTRUE(args$formal)) 10L else repeats
heavy_pass <- nrow(heavy_summary) > 0L &&
  all(heavy_summary$pairs >= required_pairs & heavy_summary$ok_pairs >=
      required_pairs & heavy_summary$timing_status == "PASS")
light_pass <- nrow(light_summary) > 0L &&
  all(light_summary$timing_status == "PASS")
test_false_pass <- nrow(guards) > 0L && all(guards$status == "ok") &&
  all(!nzchar(guards$warning_text)) && all(guards$traits %in% batch_traits)
memory_pass <- nrow(memory) > 0L && all(memory$stage ==
                                        "Stage 2D2 candidate benchmark")
parity_pass <- nrow(parity) > 0L && all(parity$pass)
overall <- isTRUE(hooks$same_translation_unit) && parity_pass && heavy_pass &&
  light_pass && test_false_pass && memory_pass
reason <- paste0(
  "parity=", parity_pass, "; heavy=", heavy_pass,
  "; light_nsim199=", light_pass, "; test_false=", test_false_pass,
  "; memory=", memory_pass
)
status <- data.frame(
  status = if (overall) "PASS" else "FAIL", candidate_mode = hooks$mode,
  same_translation_unit = isTRUE(hooks$same_translation_unit),
  source_cpp = hooks$cpp_file, export_layout = hooks$layout,
  candidate_b_gate = hooks$candidate_b_gate,
  formal_heavy_cells = nrow(heavy_summary),
  formal_heavy_min_pairs = if (nrow(heavy_summary)) min(heavy_summary$pairs) else
    NA_integer_, formal_heavy_pair_requirement = required_pairs,
  heavy_timing_gate = if (heavy_pass) "PASS" else "FAIL",
  light_nsim199_guard = if (light_pass) "PASS" else "FAIL",
  test_false_batch_guard = if (test_false_pass) "PASS" else "FAIL",
  parity_gate = if (parity_pass) "PASS" else "FAIL",
  memory_metadata = if (memory_pass) "RECORDED" else "FAIL",
  speedup_reduction_summary = if (nrow(summary)) "WRITTEN" else "MISSING",
  warnings_or_errors = if (all(raw$status == "ok") &&
                            all(!nzchar(raw$warning_text))) "NONE" else "PRESENT",
  reason = reason, production_code_changed = "NO", stringsAsFactors = FALSE
)
stage2d2_write_csv(status, status_path)
message("[stage2d2] candidate benchmark complete; production code unchanged")
if (!overall) stop("Stage 2D2 candidate benchmark gate failed: ", reason,
                   call. = FALSE)
