#!/usr/bin/env Rscript

# Stage 2D2B true blocked-evaluator benchmark.
#
# This file is audit infrastructure only.  It never edits package source,
# package bindings, tests, or the public API.  It fails closed until the
# independent exactness harness has written status=PASS.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- all_args[grepl("^--file=", all_args)]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else
  file.path(getwd(), "run_stage2d2b_benchmark.R")
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
common_file <- file.path(dirname(script_dir), "stage2d2", "stage2d2_common.R")
if (!file.exists(common_file)) {
  stop("Stage 2D2 shared helper file is unavailable: ", common_file,
       call. = FALSE)
}
source(common_file, local = TRUE)

parse_args <- function(values) {
  modes <- intersect(values, c("--smoke", "--formal", "--full", "--guards"))
  if (length(modes) != 1L) {
    stop("Use exactly one of --smoke, --formal, --full, or --guards.",
         call. = FALSE)
  }
  positional <- values[!values %in% modes]
  if (length(positional) > 2L) {
    stop("usage: run_stage2d2b_benchmark.R [--smoke|--formal|--full|--guards] ",
         "[repo] [out].", call. = FALSE)
  }
  repo <- if (length(positional)) positional[[1L]] else getwd()
  repo <- normalizePath(repo, winslash = "/", mustWork = FALSE)
  mode <- sub("^--", "", modes[[1L]])
  default_out <- file.path(repo, "benchmarks", "stage2d2b", "results", mode)
  out <- if (length(positional) >= 2L) positional[[2L]] else default_out
  list(
    mode = mode,
    repo = repo,
    out = normalizePath(out, winslash = "/", mustWork = FALSE)
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(args$out, "stage2d2b_benchmark_status.csv")

env_or <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

write_status <- function(status, reason, extra = list()) {
  base <- data.frame(
    status = as.character(status), stage = "Stage 2D2B benchmark",
    mode = args$mode, candidate = "B", production_code_changed = "NO",
    reason = as.character(reason), stringsAsFactors = FALSE
  )
  if (length(extra)) {
    for (nm in names(extra)) base[[nm]] <- extra[[nm]]
  }
  stage2d2_write_csv(base, status_path)
  invisible(base)
}

not_run <- function(reason, extra = list()) {
  write_status("NOT_RUN", reason, extra)
  stop(reason, call. = FALSE)
}

read_csv_gate <- function(path, label) {
  if (!file.exists(path)) {
    not_run(paste0(label, " status file is missing: ", path))
  }
  value <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE,
                                    check.names = FALSE),
                    error = function(e) not_run(
                      paste0("Cannot read ", label, " status file: ",
                             conditionMessage(e))))
  if (!nrow(value) || !("status" %in% names(value)) ||
      !identical(toupper(trimws(as.character(value$status[[1L]]))), "PASS")) {
    seen <- if (nrow(value) && "status" %in% names(value))
      as.character(value$status[[1L]]) else "MISSING"
    not_run(paste0(label, " gate is not PASS (status=", seen, ")."))
  }
  value
}

default_exactness_path <- file.path(
  args$repo, "benchmarks", "stage2d2b", "results", "exactness",
  "stage2d2b_exactness_status.csv"
)
compat_exactness_path <- file.path(
  args$repo, "benchmarks", "stage2d2b", "results", "correctness",
  "stage2d2b_correctness_status.csv"
)
configured_exactness <- env_or("FASTPHYLOSIG_STAGE2D2B_EXACTNESS_STATUS")
exactness_path <- if (nzchar(configured_exactness)) configured_exactness else
  if (file.exists(default_exactness_path)) default_exactness_path else
    compat_exactness_path
exactness_path <- normalizePath(exactness_path, winslash = "/", mustWork = FALSE)
# This is the first gate.  No candidate or oracle is compiled before it.
exactness_status <- read_csv_gate(exactness_path, "Stage 2D2B exactness")

if (args$mode %in% c("full", "guards")) {
  prereq_name <- if (args$mode == "full") "evaluation" else "full-pipeline"
  default_prereq <- file.path(args$repo, "benchmarks", "stage2d2b", "results",
                               if (args$mode == "full") "formal" else "full",
                               "stage2d2b_benchmark_status.csv")
  prereq_path <- normalizePath(env_or(
    if (args$mode == "full") "FASTPHYLOSIG_STAGE2D2B_EVALUATION_STATUS" else
      "FASTPHYLOSIG_STAGE2D2B_FULL_STATUS", default_prereq
  ), winslash = "/", mustWork = FALSE)
  prereq_status <- read_csv_gate(prereq_path, prereq_name)
} else {
  prereq_path <- ""
  prereq_status <- data.frame()
}

Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")

parse_ints <- function(name, default, minimum = 1L, allowed = NULL) {
  value <- env_or(name, paste(default, collapse = ","))
  out <- suppressWarnings(as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])))
  if (!length(out) || anyNA(out) || any(out < minimum) ||
      (!is.null(allowed) && any(!out %in% allowed))) {
    suffix <- if (is.null(allowed)) paste0("integers >= ", minimum) else
      paste0("values in {", paste(allowed, collapse = ","), "}")
    stop(name, " must contain comma-separated ", suffix, ".", call. = FALSE)
  }
  unique(out)
}

parse_shapes <- function(default) {
  value <- env_or("FASTPHYLOSIG_STAGE2D2B_SHAPES", paste(default, collapse = ","))
  out <- unique(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (!length(out) || any(!nzchar(out)) ||
      any(!out %in% c("balanced", "random", "pectinate"))) {
    stop("FASTPHYLOSIG_STAGE2D2B_SHAPES must contain only balanced, random, ",
         "or pectinate.", call. = FALSE)
  }
  out
}

parse_scalar <- function(name, default, minimum = 1L) {
  out <- suppressWarnings(as.integer(env_or(name, as.character(default))))
  if (length(out) != 1L || is.na(out) || out < minimum) {
    stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
  }
  out
}

is_smoke <- identical(args$mode, "smoke")
is_formal <- identical(args$mode, "formal")
is_full <- identical(args$mode, "full")
is_guards <- identical(args$mode, "guards")

default_shapes <- if (is_smoke) "balanced" else c("balanced", "random", "pectinate")
shapes <- parse_shapes(default_shapes)
if (is_formal) {
  n_grid <- parse_ints("FASTPHYLOSIG_STAGE2D2B_N",
                       c(2000L, 5000L, 10000L), minimum = 2L)
  nsim_grid <- parse_ints("FASTPHYLOSIG_STAGE2D2B_NSIM",
                          c(999L, 9999L), minimum = 1L)
  block_grid <- parse_ints("FASTPHYLOSIG_STAGE2D2B_BLOCKS",
                           c(2L, 4L, 8L, 16L), minimum = 1L,
                           allowed = c(2L, 4L, 8L, 16L))
  thread_grid <- 1L
} else if (is_smoke) {
  n_grid <- parse_ints("FASTPHYLOSIG_STAGE2D2B_N", c(50L, 500L), minimum = 2L)
  nsim_grid <- parse_ints("FASTPHYLOSIG_STAGE2D2B_NSIM", c(19L, 199L), minimum = 1L)
  block_grid <- parse_ints("FASTPHYLOSIG_STAGE2D2B_BLOCKS",
                           c(2L, 4L, 8L, 16L), minimum = 1L,
                           allowed = c(2L, 4L, 8L, 16L))
  thread_grid <- 1L
} else if (is_full) {
  n_grid <- parse_ints("FASTPHYLOSIG_STAGE2D2B_N",
                       c(2000L, 5000L, 10000L), minimum = 2L)
  nsim_grid <- parse_ints("FASTPHYLOSIG_STAGE2D2B_NSIM", 9999L, minimum = 1L)
  block_grid <- NA_integer_
  thread_grid <- 1L
} else {
  n_grid <- parse_ints("FASTPHYLOSIG_STAGE2D2B_GUARD_N",
                       parse_scalar("FASTPHYLOSIG_STAGE2D2B_GUARD_N", 500L,
                                    minimum = 2L), minimum = 2L)
  nsim_grid <- parse_ints("FASTPHYLOSIG_STAGE2D2B_GUARD_NSIM", 199L, minimum = 1L)
  block_grid <- NA_integer_
  thread_grid <- parse_ints("FASTPHYLOSIG_STAGE2D2B_GUARD_THREADS",
                            c(1L, 2L, 4L, 8L), minimum = 1L)
}

repeats <- parse_scalar(
  "FASTPHYLOSIG_STAGE2D2B_REPEATS",
  if (is_formal) 10L else if (is_smoke) 2L else 3L,
  minimum = 1L
)
if (is_formal && repeats < 10L) {
  stop("Formal evaluation-only timing requires at least 10 paired repeats.",
       call. = FALSE)
}
trait_chunk <- parse_scalar("FASTPHYLOSIG_STAGE2D2B_TRAIT_CHUNK", 64L)
simulation_chunk <- parse_scalar("FASTPHYLOSIG_STAGE2D2B_SIMULATION_CHUNK", 128L)
base_seed <- parse_scalar("FASTPHYLOSIG_STAGE2D2B_SEED", 20260910L)
guard_traits <- parse_ints("FASTPHYLOSIG_STAGE2D2B_GUARD_TRAITS",
                           c(1L, 8L, 32L, 100L), minimum = 1L)

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256"))
  }
  if (requireNamespace("openssl", quietly = TRUE)) {
    return(paste(as.character(openssl::sha256(file(path))), collapse = ""))
  }
  NA_character_
}

git_head <- function(repo) {
  value <- tryCatch(system2("git", c("-C", repo, "rev-parse", "HEAD"),
                            stdout = TRUE, stderr = FALSE),
                    error = function(e) character())
  if (length(value)) paste(value, collapse = "\n") else "UNAVAILABLE"
}

find_export <- function(env, candidates) {
  for (name in candidates) {
    if (exists(name, envir = env, inherits = FALSE)) {
      value <- get(name, envir = env, inherits = FALSE)
      if (is.function(value)) return(list(name = name, fun = value))
    }
  }
  NULL
}

mode_wrapper <- function(target, mode) {
  force(target)
  force(mode)
  function(...) {
    call_args <- list(...)
    call_args$mode <- mode
    do.call(target, stage2d2_filter_call_args(target, call_args))
  }
}

load_hooks <- function() {
  cpp_file <- env_or("FASTPHYLOSIG_STAGE2D2B_PROTO_CPP")
  if (!nzchar(cpp_file)) {
    not_run("Set FASTPHYLOSIG_STAGE2D2B_PROTO_CPP to one candidate B sourceCpp translation unit.")
  }
  cpp_file <- normalizePath(cpp_file, winslash = "/", mustWork = TRUE)
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    not_run("Rcpp is required to load the Stage 2D2B translation unit.")
  }
  hook_env <- new.env(parent = globalenv())
  before <- ls(envir = hook_env, all.names = TRUE)
  Rcpp::sourceCpp(cpp_file, env = hook_env, rebuild = TRUE,
                  showOutput = FALSE, verbose = FALSE)
  exports <- setdiff(ls(envir = hook_env, all.names = TRUE), before)
  bound <- function(entry) {
    !is.null(entry) && entry$name %in% exports &&
      exists(entry$name, envir = hook_env, inherits = FALSE) &&
      identical(get(entry$name, envir = hook_env, inherits = FALSE), entry$fun)
  }
  generic <- find_export(hook_env, c(
    "stage2d2b_k_null_prototype", "stage2d2b_k_permutation_prototype",
    "stage2d2b_true_blocked_k_null"
  ))
  oracle_raw <- find_export(hook_env, c(
    "stage2d2b_k_null_oracle", "stage2d2b_k_permutation_oracle",
    "stage2d2_k_null_oracle", "stage2d2_k_permutation_oracle", "oracle"
  ))
  candidate_raw <- find_export(hook_env, c(
    "stage2d2b_k_null_candidate_b", "stage2d2b_k_permutation_candidate_b",
    "stage2d2b_candidate_b", "stage2d2_k_null_candidate_b",
    "stage2d2_k_permutation_candidate_b", "run_candidate_b", "candidate_b"
  ))
  if (!is.null(generic)) {
    if (!bound(generic)) not_run("Generic Stage 2D2B export is not sourceCpp-bound.")
    oracle <- mode_wrapper(generic$fun, "oracle")
    candidate <- mode_wrapper(generic$fun, "B")
    same_tu <- TRUE
    layout <- "generic_mode"
    oracle_name <- paste0(generic$name, "(mode=oracle)")
    candidate_name <- paste0(generic$name, "(mode=B)")
  } else {
    if (is.null(oracle_raw) || is.null(candidate_raw) ||
        !bound(oracle_raw) || !bound(candidate_raw)) {
      not_run("Stage 2D2B requires sourceCpp-bound oracle and Candidate B exports from one TU.")
    }
    oracle <- oracle_raw$fun
    candidate <- candidate_raw$fun
    same_tu <- TRUE
    layout <- "separate_exports"
    oracle_name <- oracle_raw$name
    candidate_name <- candidate_raw$name
  }
  list(
    cpp_file = cpp_file, cpp_sha256 = sha256_file(cpp_file), env = hook_env,
    exports = exports, oracle = oracle, candidate = candidate,
    oracle_name = oracle_name, candidate_name = candidate_name,
    same_translation_unit = same_tu, layout = layout
  )
}

hooks <- tryCatch(load_hooks(), error = function(e) {
  not_run(paste0("Cannot load Stage 2D2B translation unit: ", conditionMessage(e)))
})

ns <- tryCatch(stage2d2_load_package(args$repo), error = function(e) {
  not_run(paste0("Cannot load audited fastphylosig package: ", conditionMessage(e)))
})
if (!is.environment(ns)) not_run("fastphylosig namespace is unavailable.")

# If exactness provenance carries the prototype hash, do not benchmark a
# different candidate by accident.
status_hash_names <- c("source_cpp_sha256", "candidate_cpp_sha256",
                       "prototype_cpp_sha256")
hash_columns <- intersect(status_hash_names, names(exactness_status))
gate_hash <- if (length(hash_columns))
  as.character(exactness_status[[hash_columns[[1L]]]][[1L]]) else NA_character_
if (length(hash_columns) && !is.na(gate_hash) && nzchar(gate_hash) &&
    !is.na(hooks$cpp_sha256) &&
    !identical(tolower(gate_hash), tolower(hooks$cpp_sha256))) {
  not_run("Exactness status source hash does not match the candidate TU.")
}

cell_seed <- function(stream, n, nsim, traits, threads, block, pair = 0L) {
  block_value <- if (is.na(block)) 0 else as.double(block)
  raw <- as.double(base_seed) + as.double(stream) * 1000003 +
    as.double(n) * 1009 + as.double(nsim) * 9176 +
    as.double(traits) * 101 + as.double(threads) * 17 +
    block_value * 13 + as.double(pair)
  as.integer((raw %% 2147483646) + 1)
}

pick <- function(value, names, default = NULL) {
  if (is.null(value) || !is.list(value)) return(default)
  for (nm in names) if (!is.null(value[[nm]])) return(value[[nm]])
  default
}

number <- function(value, names, default = NA_real_) {
  out <- pick(value, names, NULL)
  if (is.null(out) || !length(out)) return(default)
  out <- suppressWarnings(as.numeric(out[[1L]]))
  if (is.na(out)) default else out
}

text_value <- function(value, names, default = NA_character_) {
  out <- pick(value, names, NULL)
  if (is.null(out) || !length(out)) return(default)
  as.character(out[[1L]])
}

result_meta <- function(value) {
  if (!is.list(value)) return(list())
  meta <- value$candidate_metadata
  memory <- value$memory
  list(
    old_compute_one_calls = number(value, c("candidate_old_compute_one_calls"),
                                   number(meta, c("candidate_old_compute_one_calls",
                                                 "old_compute_one_calls"))),
    oracle_call_count = number(value, c("oracle_call_count"),
                               number(meta, c("oracle_call_count"))),
    candidate_compute_count = number(value, c("candidate_compute_count"),
                                     number(meta, c("candidate_compute_count",
                                                   "candidate_compute_calls"))),
    blocks_processed = number(value, c("blocks_processed"),
                              number(meta, c("blocks_processed"))),
    generated_permutations = number(value, c("generated_permutations"),
                                     number(meta, c("generated_permutations"))),
    workspace_bytes = number(value, c("workspace_bytes"),
                             number(memory, c("workspace_bytes",
                                              "stage2d2b_workspace_bytes"))),
    workspace_peak_bytes = number(value, c("workspace_peak_bytes",
                                           "stage2d2b_workspace_peak_bytes"),
                                  number(memory, c("workspace_peak_bytes",
                                                   "stage2d2b_workspace_peak_bytes",
                                                   "peak_workspace_bytes"))),
    permutation_workspace_bytes = number(value, c("permutation_workspace_bytes"),
                                          number(memory, c("permutation_workspace_bytes",
                                                           "peak_index_bytes"))),
    message_workspace_bytes = number(value, c("message_workspace_bytes"),
                                      number(memory, c("message_workspace_bytes"))),
    state_workspace_bytes = number(value, c("state_workspace_bytes"),
                                    number(memory, c("state_workspace_bytes"))),
    bounded = text_value(value, c("bounded", "memory_bounded"),
                         text_value(memory, c("bounded", "memory_bounded"))),
    working_memory_model = text_value(
      value, c("working_memory_model"),
      text_value(memory, c("working_memory_model"))
    ),
    memory_note = text_value(value, c("memory_note", "note"),
                             text_value(memory, c("memory_note", "note")))
  )
}

counter_gate <- function(capture, require_positive = TRUE) {
  if (!identical(capture$status, "ok")) {
    return(list(pass = FALSE, old = NA_real_, oracle = NA_real_, compute = NA_real_,
                blocks = NA_real_, reason = "candidate call did not return"))
  }
  meta <- result_meta(capture$value)
  positive <- !require_positive || (is.finite(meta$candidate_compute_count) &&
                                    meta$candidate_compute_count > 0)
  pass <- is.finite(meta$old_compute_one_calls) &&
    meta$old_compute_one_calls == 0 &&
    is.finite(meta$oracle_call_count) && meta$oracle_call_count == 0 && positive
  list(
    pass = isTRUE(pass), old = meta$old_compute_one_calls,
    oracle = meta$oracle_call_count, compute = meta$candidate_compute_count,
    blocks = meta$blocks_processed,
    reason = if (pass) "hard counters pass" else
      "missing/non-zero old_compute_one/oracle counter or non-positive candidate count"
  )
}

memory_record <- function(capture, cfg, implementation, phase) {
  meta <- result_meta(capture$value)
  n_total <- max(1, 2 * as.double(cfg$n) - 1)
  b <- if (is.na(cfg$block_size)) 0 else as.double(cfg$block_size)
  chunk <- max(1, min(as.double(cfg$trait_chunk), as.double(cfg$traits)))
  permutation_bound <- as.double(cfg$n) * b * 4
  four_array_bound <- n_total * b * chunk * 4 * 8 + b * chunk * 2 * 8
  model <- tolower(ifelse(is.na(meta$working_memory_model), "", meta$working_memory_model))
  bounded <- tolower(ifelse(is.na(meta$bounded), "", meta$bounded))
  forbidden <- grepl("n.?times.?nsim|n.?x.?nsim|quadratic", model) ||
    grepl("n.?times.?nsim|n.?x.?nsim|quadratic", bounded)
  has_memory <- any(is.finite(c(meta$workspace_bytes, meta$workspace_peak_bytes,
                               meta$permutation_workspace_bytes,
                               meta$message_workspace_bytes,
                               meta$state_workspace_bytes))) ||
    nzchar(model) || nzchar(bounded)
  data.frame(
    stage = "Stage 2D2B memory audit", phase = phase,
    implementation = implementation, shape = cfg$shape, n = as.integer(cfg$n),
    nsim = as.integer(cfg$nsim), traits = as.integer(cfg$traits),
    n_threads = as.integer(cfg$n_threads), block_size = as.integer(cfg$block_size),
    pair = as.integer(cfg$pair), workspace_bytes = meta$workspace_bytes,
    workspace_peak_bytes = meta$workspace_peak_bytes,
    permutation_workspace_bytes = meta$permutation_workspace_bytes,
    message_workspace_bytes = meta$message_workspace_bytes,
    state_workspace_bytes = meta$state_workspace_bytes,
    independent_permutation_bound_bytes = permutation_bound,
    independent_four_array_bound_bytes = four_array_bound,
    bounded = meta$bounded, working_memory_model = meta$working_memory_model,
    forbidden_n_x_nsim_workspace = forbidden,
    memory_metadata_present = has_memory,
    memory_status = if (implementation != "candidate") "NOT_APPLICABLE" else
      if (has_memory && !forbidden) "PASS" else "FAIL",
    memory_note = meta$memory_note, stringsAsFactors = FALSE
  )
}

raw_rows <- list()
memory_rows <- list()

record_capture <- function(capture, implementation, cfg, phase, order,
                           generation_s = NA_real_) {
  meta <- result_meta(capture$value)
  gate <- if (implementation == "candidate") counter_gate(capture) else
    list(pass = NA, old = NA_real_, oracle = NA_real_, compute = NA_real_,
         blocks = NA_real_, reason = "not applicable")
  raw_rows[[length(raw_rows) + 1L]] <<- data.frame(
    stage = "Stage 2D2B benchmark", phase = phase, route = cfg$route,
    implementation = implementation, candidate = "B", shape = cfg$shape,
    n = as.integer(cfg$n), nsim = as.integer(cfg$nsim), traits = as.integer(cfg$traits),
    n_threads = as.integer(cfg$n_threads), block_size = as.integer(cfg$block_size),
    trait_chunk = as.integer(cfg$trait_chunk), simulation_chunk = as.integer(cfg$simulation_chunk),
    trait_kind = cfg$trait_kind, mask_kind = cfg$mask_kind,
    pair = as.integer(cfg$pair), seed = as.integer(cfg$seed), order = order,
    permutation_generation_s = generation_s, elapsed_s = capture$elapsed_s,
    status = capture$status, warning_text = capture$warning_text,
    warning_classes = capture$warning_classes, error_class = capture$error_class,
    error_message = capture$error_message, old_compute_one_calls = gate$old,
    oracle_call_count = gate$oracle, candidate_compute_count = gate$compute,
    blocks_processed = gate$blocks, counter_status = if (implementation == "candidate")
      if (gate$pass) "PASS" else "FAIL" else "NOT_APPLICABLE",
    workspace_peak_bytes = meta$workspace_peak_bytes,
    working_memory_model = meta$working_memory_model,
    stringsAsFactors = FALSE
  )
  memory_rows[[length(memory_rows) + 1L]] <<-
    memory_record(capture, cfg, implementation, phase)
  invisible(capture)
}

call_engine <- function(fun, compiled_tree, X, nsim, permutations, include_observed,
                        n_threads, block_size, seed) {
  fml <- tryCatch(names(formals(fun)), error = function(e) character())
  call_args <- list(
    compiled_tree = compiled_tree, X = X, nsim = as.integer(nsim),
    permutations = permutations, trait_chunk = as.integer(trait_chunk),
    return_sim = FALSE, include_observed = isTRUE(include_observed),
    n_threads = as.integer(n_threads), simulation_chunk = as.integer(simulation_chunk),
    # Ordered null payloads are required by the independent exactness gate.
    # Timing isolates the evaluator and therefore does not request that large
    # result matrix again; candidate counters and P/MCSE are still returned.
    return_ordered = FALSE, block_size = as.integer(block_size),
    collect_counters = TRUE
  )
  if (!is.na(block_size) && length(fml) && !"block_size" %in% fml &&
      !"..." %in% fml) call_args$block_size <- NULL
  stage2d2_capture(function() stage2d2_call_with_seed(fun, call_args, seed))
}

make_fixture <- function(n, shape, traits = 1L, kind = "normal", seed = base_seed) {
  tree <- stage2d2_make_tree(n, shape)
  X <- stage2d2_make_matrix(tree, traits = traits, kind = kind, seed = seed)
  list(tree = tree, X = X, compiled_tree = stage2d2_compile_tree(tree, ns),
       n = as.integer(n), shape = shape, traits = as.integer(traits), kind = kind)
}

run_pair <- function(fixture, route, nsim, block_size, n_threads = 1L,
                     permutations = NULL, include_observed = TRUE,
                     stream = 1L, pair_repeats = repeats, trait_kind = "normal",
                     mask_kind = "none", generation_s = NA_real_) {
  cfg <- list(route = route, shape = fixture$shape, n = fixture$n, nsim = nsim,
              traits = ncol(fixture$X), n_threads = n_threads,
              block_size = block_size, trait_chunk = trait_chunk,
              simulation_chunk = simulation_chunk, trait_kind = trait_kind,
              mask_kind = mask_kind, pair = 0L,
              seed = cell_seed(stream, fixture$n, nsim, ncol(fixture$X),
                               n_threads, block_size, 0L))
  warm_oracle <- call_engine(hooks$oracle, fixture$compiled_tree, fixture$X,
                             nsim, permutations, include_observed, n_threads,
                             block_size, cfg$seed)
  warm_candidate <- call_engine(hooks$candidate, fixture$compiled_tree, fixture$X,
                                nsim, permutations, include_observed, n_threads,
                                block_size, cfg$seed)
  record_capture(warm_oracle, "oracle", cfg, "warmup", "warmup", generation_s)
  record_capture(warm_candidate, "candidate", cfg, "warmup", "warmup", generation_s)
  for (pair in seq_len(pair_repeats)) {
    cfg$pair <- pair
    cfg$seed <- cell_seed(stream, fixture$n, nsim, ncol(fixture$X),
                          n_threads, block_size, pair)
    oracle_first <- pair %% 2L == 1L
    first <- if (oracle_first) "oracle" else "candidate"
    second <- if (oracle_first) "candidate" else "oracle"
    first_fun <- if (oracle_first) hooks$oracle else hooks$candidate
    second_fun <- if (oracle_first) hooks$candidate else hooks$oracle
    first_capture <- call_engine(first_fun, fixture$compiled_tree, fixture$X,
                                 nsim, permutations, include_observed, n_threads,
                                 block_size, cfg$seed)
    second_capture <- call_engine(second_fun, fixture$compiled_tree, fixture$X,
                                  nsim, permutations, include_observed, n_threads,
                                  block_size, cfg$seed)
    record_capture(first_capture, first, cfg, "timed", first, generation_s)
    record_capture(second_capture, second, cfg, "timed", second, generation_s)
  }
  invisible(NULL)
}

generate_controlled <- function(n, nsim, seed) {
  started <- stage2d2_clock()
  perms <- stage2d2_make_permutations(n, nsim, seed = seed)
  list(permutations = perms, elapsed_s = stage2d2_clock() - started)
}

if (is_formal || is_smoke) {
  for (shape in shapes) for (n in n_grid) {
    message("[stage2d2b] evaluation shape=", shape, ", n=", n)
    for (nsim in nsim_grid) {
      generated <- generate_controlled(n, nsim,
                                       cell_seed(101L, n, nsim, 1L, 1L,
                                                 block_grid[[1L]]))
      fixture <- make_fixture(n, shape, traits = 1L, seed = cell_seed(
        102L, n, nsim, 1L, 1L, block_grid[[1L]]
      ))
      for (block_size in block_grid) {
        run_pair(fixture, "evaluation", nsim, block_size,
                 n_threads = 1L, permutations = generated$permutations,
                 include_observed = FALSE, stream = 111L,
                 pair_repeats = repeats, generation_s = generated$elapsed_s)
      }
    }
  }
}

raw <- stage2d2_rbind(raw_rows)
memory <- stage2d2_rbind(memory_rows)

summary_rows <- list()
make_summary <- function(data, route_filter = NULL) {
  if (!nrow(data)) return(data.frame(stringsAsFactors = FALSE))
  data <- data[data$phase == "timed" &
                 data$implementation %in% c("oracle", "candidate"), , drop = FALSE]
  if (!is.null(route_filter)) data <- data[data$route %in% route_filter, , drop = FALSE]
  if (!nrow(data)) return(data.frame(stringsAsFactors = FALSE))
  keys <- c("route", "shape", "n", "nsim", "traits", "n_threads",
            "block_size", "trait_kind", "mask_kind")
  key <- do.call(paste, c(data[keys], sep = "\r"))
  groups <- split(seq_len(nrow(data)), key, drop = TRUE)
  rows <- lapply(groups, function(idx) {
    z <- data[idx, , drop = FALSE]
    o <- z[z$implementation == "oracle", , drop = FALSE]
    c <- z[z$implementation == "candidate", , drop = FALSE]
    pairs <- intersect(o$pair, c$pair)
    if (length(pairs)) {
      oa <- o[match(pairs, o$pair), , drop = FALSE]
      ca <- c[match(pairs, c$pair), , drop = FALSE]
      ok <- oa$status == "ok" & ca$status == "ok" &
        !nzchar(oa$warning_text) & !nzchar(ca$warning_text) &
        ca$counter_status == "PASS" &
        is.finite(oa$elapsed_s) & is.finite(ca$elapsed_s)
      os <- oa$elapsed_s[ok]
      cs <- ca$elapsed_s[ok]
    } else {
      oa <- ca <- data.frame()
      ok <- logical()
      os <- cs <- numeric()
    }
    median_o <- if (length(os)) stats::median(os) else NA_real_
    median_c <- if (length(cs)) stats::median(cs) else NA_real_
    data.frame(
      route = z$route[[1L]], shape = z$shape[[1L]], n = z$n[[1L]],
      nsim = z$nsim[[1L]], traits = z$traits[[1L]],
      n_threads = z$n_threads[[1L]], block_size = z$block_size[[1L]],
      trait_kind = z$trait_kind[[1L]], mask_kind = z$mask_kind[[1L]],
      pairs = length(pairs), ok_pairs = sum(ok),
      oracle_median_s = median_o, oracle_iqr_s = if (length(os)) stats::IQR(os) else NA_real_,
      candidate_median_s = median_c, candidate_iqr_s = if (length(cs)) stats::IQR(cs) else NA_real_,
      speedup = if (is.finite(median_o) && is.finite(median_c) && median_c > 0)
        median_o / median_c else NA_real_,
      reduction = if (is.finite(median_o) && median_o > 0 && is.finite(median_c))
        1 - median_c / median_o else NA_real_,
      timing_status = if (length(pairs) >= if (is_formal) 10L else pair_repeats_for_summary() &&
                          sum(ok) == length(pairs)) "PASS" else "FAIL",
      warning_or_error = if (nrow(z) && (any(z$status != "ok") || any(nzchar(z$warning_text))))
        "PRESENT" else "NONE",
      counter_status = if (nrow(c) && all(c$counter_status == "PASS")) "PASS" else "FAIL",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

pair_repeats_for_summary <- function() {
  if (is_formal) 10L else repeats
}

# R resolves function bodies at call time, so the helper above is available
# when make_summary is invoked below.
summary <- if (is_formal || is_smoke) make_summary(raw, "evaluation") else
  data.frame(stringsAsFactors = FALSE)

select_block <- function(summary) {
  if (!nrow(summary)) return(list(block = NA_integer_, reason = "no valid evaluation timings"))
  heavy <- summary[summary$route == "evaluation" & summary$n >= 5000L &
                     summary$nsim == 9999L & summary$n_threads == 1L &
                     summary$timing_status == "PASS" & is.finite(summary$candidate_median_s),
                   , drop = FALSE]
  if (!nrow(heavy)) return(list(block = NA_integer_, reason = "no valid heavy evaluation cells"))
  by_block <- stats::aggregate(candidate_median_s ~ block_size, heavy, stats::median)
  names(by_block)[2L] <- "median_candidate_s"
  best <- by_block$block_size[[which.min(by_block$median_candidate_s)]]
  reason <- "minimum median candidate evaluator time"
  if (all(c(8L, 16L) %in% by_block$block_size)) {
    t8 <- by_block$median_candidate_s[by_block$block_size == 8L]
    t16 <- by_block$median_candidate_s[by_block$block_size == 16L]
    if (abs(t16 - t8) / min(t8, t16) < 0.03) {
      best <- 8L
      reason <- "blocks 8 and 16 differ by <3%; selected smaller block 8"
    }
  }
  list(block = as.integer(best), reason = reason, by_block = by_block)
}

selection <- select_block(summary)

evaluation_gate <- function(summary, selection) {
  if (!nrow(summary) || is.na(selection$block)) return(FALSE)
  heavy <- summary[summary$route == "evaluation" & summary$n >= 5000L &
                     summary$nsim == 9999L & summary$n_threads == 1L &
                     summary$block_size == selection$block, , drop = FALSE]
  if (!nrow(heavy)) return(FALSE)
  sum(heavy$timing_status == "PASS" & heavy$counter_status == "PASS" &
        is.finite(heavy$speedup) & heavy$speedup >= 1.35) >= 2L
}

eval_pass <- if (is_formal) evaluation_gate(summary, selection) else NA

if (is_full) {
  selected <- suppressWarnings(as.integer(env_or("FASTPHYLOSIG_STAGE2D2B_BLOCK_SIZE", "NA")))
  if (is.na(selected)) {
    selected_col <- intersect(c("selected_block_size", "block_size"), names(prereq_status))
    if (!length(selected_col)) not_run("Evaluation status does not record selected block size.")
    selected <- suppressWarnings(as.integer(prereq_status[[selected_col[[1L]]]][[1L]]))
  }
  if (is.na(selected) || !selected %in% c(2L, 4L, 8L, 16L)) {
    not_run("Full-pipeline benchmark requires a selected block size in {2,4,8,16}.")
  }
  block_grid <- selected
  for (shape in shapes) for (n in n_grid) {
    message("[stage2d2b] full pipeline shape=", shape, ", n=", n)
    for (nsim in nsim_grid) {
      fixture <- make_fixture(n, shape, traits = 1L,
                              seed = cell_seed(301L, n, nsim, 1L, 1L, selected))
      run_pair(fixture, "full_pipeline", nsim, selected, n_threads = 1L,
               permutations = NULL, include_observed = TRUE, stream = 311L,
               pair_repeats = repeats)
    }
  }
  raw <- stage2d2_rbind(raw_rows)
  memory <- stage2d2_rbind(memory_rows)
  summary <- make_summary(raw, "full_pipeline")
  selection <- list(block = selected, reason = "inherited from evaluation gate")
}

run_public_test_false <- function() {
  if (!is_guards) return(invisible(NULL))
  for (p in guard_traits) {
    fixture <- make_fixture(n_grid[[1L]], "balanced", traits = p,
                            seed = cell_seed(401L, n_grid[[1L]], 0L, p, 1L,
                                              selection$block))
    cfg <- list(route = "test_false", shape = "balanced", n = fixture$n,
                nsim = NA_integer_, traits = p, n_threads = 1L,
                block_size = NA_integer_, trait_chunk = trait_chunk,
                simulation_chunk = simulation_chunk, trait_kind = "normal",
                mask_kind = "none", pair = 0L,
                seed = cell_seed(402L, fixture$n, 0L, p, 1L, NA_integer_))
    call <- function() stage2d2_capture(function() fastphylosig::fast_k(
      fixture$tree, x = fixture$X, test = FALSE, ncores = 1L,
      trait_chunk = trait_chunk, simulation_chunk = simulation_chunk,
      return_sim = FALSE, keep_null = FALSE, verbose = FALSE, progress = FALSE
    ))
    warm <- call()
    record_capture(warm, "production", cfg, "warmup", "public")
    if (warm$status != "ok") next
    for (pair in seq_len(max(1L, min(repeats, 3L)))) {
      cfg$pair <- pair
      measured <- call()
      record_capture(measured, "production", cfg, "timed", "public")
    }
  }
  invisible(NULL)
}

run_guards <- function() {
  if (!is_guards) return(invisible(NULL))
  selected <- selection$block
  if (is.na(selected)) not_run("Guards require a selected block size from the prior gate.")
  guard_n <- n_grid[[1L]]
  guard_nsim <- nsim_grid[[1L]]
  # Light evaluator and batch-width guard.
  for (p in guard_traits) {
    fixture <- make_fixture(guard_n, "balanced", traits = p,
                            seed = cell_seed(501L, guard_n, guard_nsim, p, 1L, selected))
    generated <- generate_controlled(guard_n, guard_nsim,
                                     cell_seed(502L, guard_n, guard_nsim, p, 1L, selected))
    run_pair(fixture, "light_nsim199", guard_nsim, selected,
             n_threads = 1L, permutations = generated$permutations,
             include_observed = FALSE, stream = 511L, pair_repeats = repeats,
             generation_s = generated$elapsed_s)
  }
  # Numerical edge traits, including constant-result failure semantics.
  for (kind in c("large_offset", "near_constant", "constant")) {
    fixture <- make_fixture(guard_n, "balanced", traits = 1L, kind = kind,
                            seed = cell_seed(521L, guard_n, guard_nsim, 1L, 1L, selected))
    generated <- generate_controlled(guard_n, guard_nsim,
                                     cell_seed(522L, guard_n, guard_nsim, 1L, 1L, selected))
    run_pair(fixture, paste0("guard_", kind), guard_nsim, selected,
             n_threads = 1L, permutations = generated$permutations,
             include_observed = FALSE, stream = 531L, pair_repeats = repeats,
             trait_kind = kind, generation_s = generated$elapsed_s)
  }
  # Two-tip and retained-subset/NA-mask cases.
  two_tip <- make_fixture(2L, "balanced", traits = 1L,
                          seed = cell_seed(541L, 2L, guard_nsim, 1L, 1L, selected))
  two_perm <- generate_controlled(2L, guard_nsim,
                                  cell_seed(542L, 2L, guard_nsim, 1L, 1L, selected))
  run_pair(two_tip, "guard_two_tip", guard_nsim, selected, 1L,
           permutations = two_perm$permutations, include_observed = FALSE,
           stream = 551L, pair_repeats = repeats, generation_s = two_perm$elapsed_s)
  masked <- stage2d2_masked_cases(
    stage2d2_make_tree(guard_n, "balanced"), traits = 8L, kind = "shared"
  )
  if (length(masked)) {
    for (i in seq_len(min(length(masked), 2L))) {
      f <- masked[[i]]
      fixture <- list(tree = f$tree, X = f$X,
                      compiled_tree = stage2d2_compile_tree(f$tree, ns),
                      n = nrow(f$X), shape = "balanced", traits = ncol(f$X), kind = "normal")
      generated <- generate_controlled(fixture$n, guard_nsim,
                                       cell_seed(561L + i, fixture$n, guard_nsim,
                                                 fixture$traits, 1L, selected))
      run_pair(fixture, "guard_retained_subset", guard_nsim, selected, 1L,
               permutations = generated$permutations, include_observed = FALSE,
               stream = 571L + i, pair_repeats = repeats, mask_kind = "shared",
               generation_s = generated$elapsed_s)
    }
  }
  # Native thread guard uses the same controlled permutations and chosen block.
  thread_fixture <- make_fixture(guard_n, "balanced", traits = 1L,
                                 seed = cell_seed(581L, guard_n, guard_nsim, 1L,
                                                   thread_grid[[1L]], selected))
  for (threads in thread_grid) {
    generated <- generate_controlled(guard_n, guard_nsim,
                                     cell_seed(582L, guard_n, guard_nsim, 1L,
                                               threads, selected))
    run_pair(thread_fixture, "guard_threads", guard_nsim, selected, threads,
             permutations = generated$permutations, include_observed = FALSE,
             stream = 591L, pair_repeats = repeats, generation_s = generated$elapsed_s)
  }
  run_public_test_false()
  invisible(NULL)
}

if (is_full) {
  # Full pipeline uses the inherited selected block; the call above already ran it.
  invisible(NULL)
} else if (is_guards) {
  # Guards are intentionally run only after the full-pipeline prerequisite gate.
  selected_from_full <- suppressWarnings(as.integer(
    env_or("FASTPHYLOSIG_STAGE2D2B_BLOCK_SIZE", "NA")
  ))
  if (is.na(selected_from_full)) {
    selected_col <- intersect(c("selected_block_size", "block_size"),
                              names(prereq_status))
    if (!length(selected_col)) {
      not_run("Full-pipeline status does not record selected block size.")
    }
    selected_from_full <- suppressWarnings(as.integer(
      prereq_status[[selected_col[[1L]]]][[1L]]
    ))
  }
  if (is.na(selected_from_full) ||
      !selected_from_full %in% c(2L, 4L, 8L, 16L)) {
    not_run("Guards require a selected block size in {2,4,8,16}.")
  }
  selection <- list(block = selected_from_full,
                    reason = "inherited from full-pipeline gate")
  run_guards()
  raw <- stage2d2_rbind(raw_rows)
  memory <- stage2d2_rbind(memory_rows)
  summary <- make_summary(raw, c("light_nsim199", "guard_large_offset",
                                  "guard_near_constant", "guard_constant",
                                  "guard_two_tip", "guard_retained_subset",
                                  "guard_threads"))
}

if (!exists("raw", inherits = FALSE)) raw <- stage2d2_rbind(raw_rows)
if (!exists("memory", inherits = FALSE)) memory <- stage2d2_rbind(memory_rows)
if (!exists("summary", inherits = FALSE)) summary <- data.frame(stringsAsFactors = FALSE)

stage2d2_write_csv(raw, file.path(args$out, "stage2d2b_benchmark_raw.csv"))
stage2d2_write_csv(summary, file.path(args$out, "stage2d2b_benchmark_summary.csv"))
stage2d2_write_csv(memory, file.path(args$out, "stage2d2b_benchmark_memory.csv"))

memory_candidate <- if (nrow(memory)) memory[memory$implementation == "candidate", , drop = FALSE] else
  data.frame()
memory_pass <- if (!nrow(memory_candidate)) FALSE else
  all(memory_candidate$memory_status == "PASS") &&
  all(!memory_candidate$forbidden_n_x_nsim_workspace)

heavy_full_pass <- NA
full_speedup_gate <- NA
if (is_full && nrow(summary)) {
  heavy <- summary[summary$n >= 5000L & summary$nsim == 9999L &
                     summary$n_threads == 1L, , drop = FALSE]
  heavy_full_pass <- sum(heavy$timing_status == "PASS" & heavy$counter_status == "PASS" &
                           is.finite(heavy$reduction) & heavy$reduction >= 0.20) >= 2L
  full_speedup_gate <- any(heavy$timing_status == "PASS" & is.finite(heavy$speedup) &
                             heavy$speedup >= 1.25)
}

guards_timed <- if (is_guards && nrow(raw)) raw[raw$phase == "timed", , drop = FALSE] else
  data.frame()
guard_pair_pass <- if (!nrow(guards_timed)) FALSE else {
  cand <- guards_timed[guards_timed$implementation == "candidate", , drop = FALSE]
  all(cand$counter_status == "PASS") && all(cand$status == "ok") &&
    all(!nzchar(cand$warning_text))
}
test_false <- if (is_guards && nrow(raw)) raw[raw$route == "test_false" &
                                                  raw$phase == "timed", , drop = FALSE] else
  data.frame()
test_false_pass <- if (!is_guards) NA else nrow(test_false) > 0L &&
  all(test_false$implementation == "production") &&
  all(test_false$status == "ok") && all(!nzchar(test_false$warning_text))

if (is_formal || is_smoke) {
  overall <- isTRUE(hooks$same_translation_unit) &&
    if (is_smoke) {
      nrow(summary) > 0L && all(summary$timing_status == "PASS") && memory_pass
    } else {
      isTRUE(eval_pass) && memory_pass
    }
  timing_gate <- if (is_smoke) "SMOKE" else if (isTRUE(eval_pass)) "PASS" else "FAIL"
  reason <- paste0("exactness=PASS; same_tu=", hooks$same_translation_unit,
                   "; evaluation_gate=", timing_gate, "; memory=", memory_pass)
} else if (is_full) {
  overall <- isTRUE(heavy_full_pass) && isTRUE(full_speedup_gate) && memory_pass
  timing_gate <- if (overall) "PASS" else "FAIL"
  reason <- paste0("evaluation_prerequisite=PASS; full_reduction_gate=", heavy_full_pass,
                   "; full_speedup_gate=", full_speedup_gate,
                   "; memory=", memory_pass)
} else {
  overall <- isTRUE(guard_pair_pass) && isTRUE(test_false_pass)
  timing_gate <- if (overall) "PASS" else "FAIL"
  reason <- paste0("full_prerequisite=PASS; evaluator_guards=", guard_pair_pass,
                   "; test_false=", test_false_pass)
}

selected_block_value <- if (length(selection$block)) selection$block else NA_integer_
status <- data.frame(
  status = if (overall) "PASS" else "FAIL",
  stage = "Stage 2D2B benchmark", mode = args$mode, candidate = "B",
  exactness_status = "PASS", exactness_status_file = exactness_path,
  same_translation_unit = isTRUE(hooks$same_translation_unit),
  candidate_source_cpp = hooks$cpp_file,
  candidate_source_cpp_sha256 = hooks$cpp_sha256,
  oracle_export = hooks$oracle_name, candidate_export = hooks$candidate_name,
  selected_block_size = selected_block_value,
  block_selection_reason = if (length(selection$reason)) selection$reason else NA_character_,
  evaluation_heavy_pass_cells = if (is_formal && nrow(summary)) {
    heavy <- summary[summary$n >= 5000L & summary$nsim == 9999L &
                       summary$n_threads == 1L & summary$block_size == selected_block_value,
                     , drop = FALSE]
    sum(heavy$timing_status == "PASS" & is.finite(heavy$speedup) & heavy$speedup >= 1.35)
  } else NA_integer_,
  evaluation_gate = if (is_formal) if (isTRUE(eval_pass)) "PASS" else "FAIL" else
    if (is_smoke) "NOT_APPLICABLE" else "INHERITED",
  full_reduction_gate = if (is_full) if (isTRUE(heavy_full_pass)) "PASS" else "FAIL" else
    "NOT_RUN",
  full_speedup_gate = if (is_full) if (isTRUE(full_speedup_gate)) "PASS" else "FAIL" else
    "NOT_RUN",
  light_batch_parallel_guard = if (is_guards) if (isTRUE(guard_pair_pass)) "PASS" else "FAIL" else
    "NOT_RUN",
  test_false_guard = if (is_guards) if (isTRUE(test_false_pass)) "PASS" else "FAIL" else
    "NOT_RUN",
  memory_gate = if (memory_pass) "PASS" else "FAIL",
  old_compute_one_counter = if (nrow(raw) && any(raw$implementation == "candidate"))
    if (all(raw$counter_status[raw$implementation == "candidate"] == "PASS")) "0" else "FAIL" else
    "MISSING",
  oracle_call_counter = if (nrow(raw) && any(raw$implementation == "candidate"))
    if (all(raw$oracle_call_count[raw$implementation == "candidate"] == 0)) "0" else "FAIL" else
    "MISSING",
  warnings_or_errors = if (nrow(raw) && all(raw$status == "ok") &&
                            all(!nzchar(raw$warning_text))) "NONE" else "PRESENT",
  timing_protocol = "same-TU; pre-generated controlled permutations for evaluation; warmup; alternating; serialized",
  repeats = repeats, n_grid = paste(n_grid, collapse = ","),
  nsim_grid = paste(nsim_grid, collapse = ","),
  shapes = paste(shapes, collapse = ","),
  threads = paste(thread_grid, collapse = ","),
  reason = reason, production_code_changed = "NO", stringsAsFactors = FALSE
)
stage2d2_write_csv(status, status_path)

provenance <- data.frame(
  stage = "Stage 2D2B benchmark", mode = args$mode, candidate = "B",
  source_commit = git_head(args$repo), package_version = tryCatch(
    as.character(utils::packageVersion("fastphylosig")), error = function(e) NA_character_
  ), R_version = R.version.string, platform = R.version$platform,
  os = Sys.info()[["sysname"]], machine = Sys.info()[["machine"]],
  locale_collate = Sys.getlocale("LC_COLLATE"), omp_threads = env_or("OMP_NUM_THREADS"),
  exactness_status_file = exactness_path,
  exactness_status_sha256 = sha256_file(exactness_path),
  candidate_source_cpp = hooks$cpp_file,
  candidate_source_cpp_sha256 = hooks$cpp_sha256,
  oracle_export = hooks$oracle_name, candidate_export = hooks$candidate_name,
  same_translation_unit = isTRUE(hooks$same_translation_unit),
  n_grid = paste(n_grid, collapse = ","), nsim_grid = paste(nsim_grid, collapse = ","),
  block_grid = paste(block_grid, collapse = ","), shapes = paste(shapes, collapse = ","),
  thread_grid = paste(thread_grid, collapse = ","), repeats = repeats,
  trait_chunk = trait_chunk, simulation_chunk = simulation_chunk,
  fixed_fixture = TRUE, warmup = TRUE, alternating_order = TRUE,
  serialized = TRUE, controlled_permutations_pre_generated = is_formal || is_smoke || is_guards,
  return_sim = FALSE, production_code_changed = "NO", stringsAsFactors = FALSE
)
stage2d2_write_csv(provenance, file.path(args$out, "stage2d2b_benchmark_provenance.csv"))

message("[stage2d2b] benchmark complete; status=", status$status[[1L]],
        "; production code unchanged")
