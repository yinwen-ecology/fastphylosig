# Stage 2D2B true blocked K null evaluator correctness gate.
#
# This audit runner is opt-in and writes evidence only under its requested
# output directory.  It does not modify package source, tests, public API, or
# the frozen Stage 2D2 evidence.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- all_args[grepl("^--file=", all_args)]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else
  file.path(getwd(), "run_stage2d2b_correctness.R")
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
source(file.path(script_dir, "stage2d2b_common.R"), local = TRUE)

args <- stage2d2b_parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(args$out, "stage2d2b_correctness_status.csv")

write_status <- function(status, reason, extra = list()) {
  row <- c(list(
    status = as.character(status), mode = args$mode,
    reason = as.character(reason), production_code_changed = "NO",
    public_api_changed = "NO", rng_policy_changed = "NO"
  ), extra)
  stage2d2b_write_csv(as.data.frame(row, stringsAsFactors = FALSE), status_path)
}

fail_before_run <- function(reason) {
  write_status("NOT_RUN", reason)
  stop(reason, call. = FALSE)
}

Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")

hooks <- tryCatch(stage2d2b_load_hooks(args$repo), error = function(e) e)
if (inherits(hooks, "error")) fail_before_run(conditionMessage(hooks))
if (!isTRUE(hooks$same_translation_unit)) {
  fail_before_run("oracle and Candidate B are not sourceCpp-bound in one TU.")
}

ns <- tryCatch(stage2d2b_load_package(args$repo), error = function(e) e)
if (inherits(ns, "error")) {
  fail_before_run(paste0("Cannot load audited fastphylosig: ",
                         conditionMessage(ns)))
}
if (!is.environment(ns)) fail_before_run("fastphylosig namespace unavailable.")

parse_grid <- function(name, formal_default, smoke_default, minimum = 1L) {
  value <- stage2d2b_env(name)
  if (nzchar(value)) return(stage2d2b_ints(value, name, minimum))
  if (isTRUE(args$smoke)) smoke_default else formal_default
}

parse_scalar <- function(name, formal_default, smoke_default, minimum = 1L) {
  value <- stage2d2b_env(name)
  if (!nzchar(value)) value <- as.character(
    if (isTRUE(args$smoke)) smoke_default else formal_default
  )
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < minimum) {
    stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
  }
  out
}

parse_shapes <- function() {
  value <- stage2d2b_env("FASTPHYLOSIG_STAGE2D2B_SHAPES")
  if (nzchar(value)) return(stage2d2b_strings(
    value, "FASTPHYLOSIG_STAGE2D2B_SHAPES",
    c("balanced", "random", "pectinate")
  ))
  c("balanced", "random", "pectinate")
}

shapes <- parse_shapes()
n_grid <- parse_grid("FASTPHYLOSIG_STAGE2D2B_N",
                     c(50L, 500L, 2000L, 5000L), c(50L), minimum = 2L)
blocks <- parse_grid("FASTPHYLOSIG_STAGE2D2B_BLOCKS",
                     c(2L, 4L, 8L, 16L), c(2L, 4L), minimum = 1L)
if (any(!blocks %in% c(2L, 4L, 8L, 16L))) {
  stop("FASTPHYLOSIG_STAGE2D2B_BLOCKS must use 2, 4, 8, or 16.",
       call. = FALSE)
}
nsim_controlled <- parse_scalar("FASTPHYLOSIG_STAGE2D2B_NSIM",
                                19L, 7L, minimum = 1L)
traits_grid <- parse_grid("FASTPHYLOSIG_STAGE2D2B_TRAITS",
                          if (isTRUE(args$smoke)) c(1L, 8L) else c(1L),
                          c(1L, 8L), minimum = 1L)
replay_threads <- parse_grid("FASTPHYLOSIG_STAGE2D2B_REPLAY_THREADS",
                             c(1L, 2L), c(1L), minimum = 1L)
replay_nsim <- parse_scalar("FASTPHYLOSIG_STAGE2D2B_REPLAY_NSIM",
                            19L, 7L, minimum = 1L)
base_seed <- parse_scalar("FASTPHYLOSIG_STAGE2D2B_SEED",
                          20260910L, 20260910L, minimum = 1L)

patterns <- if (isTRUE(args$smoke)) {
  c("identity", "single_swap", "reverse")
} else stage2d2b_permutation_patterns
pattern_value <- stage2d2b_env("FASTPHYLOSIG_STAGE2D2B_PATTERNS")
if (nzchar(pattern_value)) patterns <- stage2d2b_strings(
  pattern_value, "FASTPHYLOSIG_STAGE2D2B_PATTERNS",
  stage2d2b_permutation_patterns
)

compile_tree <- function(tree) stage2d2_compile_tree(tree, ns)

fixture <- function(n, shape, traits, seed) {
  tree <- stage2d2_make_tree(n, shape)
  X <- stage2d2_make_matrix(tree, traits = traits, kind = "normal", seed = seed)
  list(tree = tree, compiled = compile_tree(tree), X = X,
       n = as.integer(n), shape = shape, traits = as.integer(traits))
}

cell_seed <- function(stream, n, nsim, traits, block, thread = 1L, pair = 0L) {
  raw <- as.double(base_seed) + as.double(stream) * 1000003 +
    as.double(n) * 1009 + as.double(nsim) * 9176 +
    as.double(traits) * 101 + as.double(block) * 13 +
    as.double(thread) * 17 + as.double(pair)
  as.integer((raw %% 2147483646) + 1)
}

trait_chunk_for <- function(p) as.integer(min(64L, max(1L, p)))

max_field <- function(compared, name) {
  value <- compared$result$max_abs
  if (is.null(value) || is.null(value[[name]])) NA_real_ else
    as.numeric(value[[name]])
}

ordered_payload_ok <- function(captured, nsim, p) {
  if (captured$status != "ok") return(FALSE)
  normalized <- stage2d2b_normalize(captured$value, nsim, p)
  !is.null(normalized) && !is.null(normalized$sim_K) &&
    identical(dim(normalized$sim_K), c(as.integer(nsim), as.integer(p)))
}

counter_cfg <- function(check, implementation, f, nsim, pattern, block,
                        threads, seed) {
  list(check = check, implementation = implementation, shape = f$shape,
       n = f$n, traits = f$traits, nsim = nsim, pattern = pattern,
       block_size = block, n_threads = threads, seed = seed)
}

summary_row <- function(check, f, nsim, pattern, block, threads, seed,
                        compared, old_accounting, new_accounting,
                        old_capture, new_capture, counter, ordered, pass,
                        detail) {
  data.frame(
    check = check, shape = f$shape, n = f$n, traits = f$traits,
    nsim = nsim, pattern = pattern, block_size = block,
    n_threads = threads, seed = seed, pass = isTRUE(pass),
    bitwise_payload_pass = isTRUE(compared$pass),
    contract_fields_present = isTRUE(compared$contract_fields_present),
    ordered_return_pass = isTRUE(ordered),
    old_accounting_pass = isTRUE(old_accounting$pass),
    candidate_accounting_pass = isTRUE(new_accounting$pass),
    requested_pass = isTRUE(new_accounting$requested_ok),
    successful_failed_pass = isTRUE(new_accounting$counts_ok),
    identity_first_pass = isTRUE(new_accounting$identity_first_ok),
    P_formula_pass = isTRUE(new_accounting$p_formula_ok),
    old_status = old_capture$status, candidate_status = new_capture$status,
    old_warning = old_capture$warning_text,
    candidate_warning = new_capture$warning_text,
    max_abs_K = max_field(compared, "K"),
    max_abs_P = max_field(compared, "P"),
    max_abs_MCSE_P = max_field(compared, "MCSE_P"),
    max_abs_exceedance = max_field(compared, "exceedance"),
    max_abs_ordered_null_K = max_field(compared, "sim_K"),
    counter_gate_pass = isTRUE(counter$counter_gate_pass),
    old_compute_one_calls = counter$old_compute_one_calls,
    oracle_call_count = counter$oracle_call_count,
    candidate_compute_count = counter$candidate_compute_count,
    detail = as.character(detail), stringsAsFactors = FALSE
  )
}

counter_rows <- list()
summary_rows <- list()
replicate_rows <- list()
failure_rows <- list()
immutability_rows <- list()

run_pair <- function(check, f, X, nsim, perms, pattern, block, threads,
                     seed, include_observed = FALSE, require_identity = TRUE) {
  old <- stage2d2b_capture(function() stage2d2b_call_engine(
    hooks$oracle, f$compiled, X, nsim, permutations = perms,
    block_size = block, n_threads = threads,
    trait_chunk = trait_chunk_for(ncol(X)), include_observed = include_observed,
    seed = seed, return_sim = TRUE
  ))
  candidate <- stage2d2b_capture(function() stage2d2b_call_engine(
    hooks$candidate, f$compiled, X, nsim, permutations = perms,
    block_size = block, n_threads = threads,
    trait_chunk = trait_chunk_for(ncol(X)), include_observed = include_observed,
    seed = seed, return_sim = TRUE
  ))
  compared <- stage2d2b_compare(old, candidate, nsim, ncol(X))
  old_accounting <- stage2d2b_accounting(
    old, nsim, ncol(X), include_observed, controlled = !is.null(perms),
    require_identity = require_identity
  )
  new_accounting <- stage2d2b_accounting(
    candidate, nsim, ncol(X), include_observed, controlled = !is.null(perms),
    require_identity = require_identity
  )
  cfg <- counter_cfg(check, "candidate", f, nsim, pattern, block, threads, seed)
  counter <- stage2d2b_counter_row(candidate, cfg)
  counter_rows[[length(counter_rows) + 1L]] <<- counter
  ordered <- ordered_payload_ok(candidate, nsim, ncol(X))
  pass <- isTRUE(compared$pass) && isTRUE(old_accounting$pass) &&
    isTRUE(new_accounting$pass) && isTRUE(ordered) &&
    isTRUE(counter$counter_gate_pass)
  detail <- if (pass) "oracle/candidate exact with valid blocked provenance" else
    paste0("comparison=", compared$pass,
           "; accounting=", old_accounting$pass, "/", new_accounting$pass,
           "; ordered=", ordered, "; counter=", counter$counter_gate_pass,
           "; candidate_error=", candidate$error_message)
  summary_rows[[length(summary_rows) + 1L]] <<- summary_row(
    check, f, nsim, pattern, block, threads, seed, compared,
    old_accounting, new_accounting, old, candidate, counter, ordered, pass,
    detail
  )
  rows <- stage2d2b_replicate_rows(compared, list(
    check = check, shape = f$shape, n = f$n, pattern = pattern,
    block_size = block, n_threads = threads
  ))
  if (nrow(rows)) replicate_rows[[length(replicate_rows) + 1L]] <<- rows
  list(old = old, candidate = candidate, compared = compared,
       old_accounting = old_accounting, candidate_accounting = new_accounting,
       counter = counter, ordered = ordered, pass = pass)
}

# Gate 1: controlled permutations.  This loop is intentionally first; no
# internal RNG replay or failure fixture is attempted until every controlled
# ordered-null replicate has passed bitwise and accounting comparison.
message("[stage2d2b] controlled permutation exactness")
for (shape in shapes) for (n in n_grid) for (p in traits_grid) {
  f <- fixture(n, shape, p,
               cell_seed(101L, n, nsim_controlled, p, 1L))
  for (pattern in patterns) {
    perms <- stage2d2b_make_permutations(
      n, nsim_controlled, pattern,
      seed = cell_seed(102L, n, nsim_controlled, p, 1L)
    )
    for (block in blocks) {
      run_pair("controlled_permutation", f, f$X, nsim_controlled, perms,
               pattern, block, 1L,
               cell_seed(103L, n, nsim_controlled, p, block),
               include_observed = FALSE, require_identity = TRUE)
    }
  }
}

controlled <- stage2d2b_rbind(summary_rows)
controlled_pass <- nrow(controlled) > 0L && all(controlled$pass)
stage2d2b_write_csv(controlled, file.path(args$out,
                                          "stage2d2b_controlled_summary.csv"))
stage2d2b_write_csv(stage2d2b_rbind(replicate_rows), file.path(
  args$out, "stage2d2b_controlled_replicates.csv"
))
stage2d2b_write_csv(stage2d2b_rbind(counter_rows), file.path(
  args$out, "stage2d2b_controlled_counters.csv"
))
if (!controlled_pass) {
  write_status("FAIL", "controlled permutation exactness/accounting gate failed",
               list(controlled_rows = nrow(controlled),
                    controlled_failures = sum(!controlled$pass)))
  stop("Stage 2D2B controlled permutation gate failed.", call. = FALSE)
}

# Gate 2: identity-first and same-seed RNG replay.  The first returned row is
# required to be the observed identity result; the frozen accounting fields
# carry the no-draw identity contract.
message("[stage2d2b] identity-first and RNG replay")
rng_rows <- list()
for (shape in shapes) for (n in n_grid) {
  f <- fixture(n, shape, 1L,
               cell_seed(201L, n, replay_nsim, 1L, 1L))
  for (threads in replay_threads) for (block in blocks) {
    seed <- cell_seed(202L, n, replay_nsim, 1L, block, threads)
    first <- run_pair("identity_first_rng", f, f$X, replay_nsim, NULL,
                      "internal_rng", block, threads, seed,
                      include_observed = TRUE, require_identity = TRUE)
    repeated <- stage2d2b_capture(function() stage2d2b_call_engine(
      hooks$candidate, f$compiled, f$X, replay_nsim,
      permutations = NULL, block_size = block, n_threads = threads,
      trait_chunk = 1L, include_observed = TRUE, seed = seed,
      return_sim = TRUE
    ))
    replay_cmp <- if (first$candidate$status == "ok" &&
                      repeated$status == "ok") {
      # Do not treat the original oracle as the replay reference: both
      # candidate calls are restored to exactly the same R seed.
      stage2d2b_compare(first$candidate, repeated, replay_nsim, 1L)
    } else list(pass = FALSE, result = list(max_abs = numeric()),
                contract_fields_present = FALSE)
    replay_max <- if (is.null(replay_cmp$result$max_abs) ||
                      !length(replay_cmp$result$max_abs)) NA_real_ else {
      z <- as.numeric(replay_cmp$result$max_abs)
      if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE)
    }
    rng_rows[[length(rng_rows) + 1L]] <- data.frame(
      check = "same_seed_rng_replay", shape = shape, n = n,
      nsim = replay_nsim, block_size = block, n_threads = threads,
      identity_first = isTRUE(first$candidate_accounting$identity_first_ok),
      requested_ok = isTRUE(first$candidate_accounting$requested_ok),
      successful_failed_ok = isTRUE(first$candidate_accounting$counts_ok),
      candidate_vs_oracle = isTRUE(first$pass),
      candidate_replay = isTRUE(replay_cmp$pass),
      max_abs_replay = replay_max, pass = isTRUE(first$pass) &&
        isTRUE(replay_cmp$pass), stringsAsFactors = FALSE
    )
  }
}
rng_summary <- stage2d2b_rbind(rng_rows)
stage2d2b_write_csv(rng_summary, file.path(args$out,
                                           "stage2d2b_rng_replay.csv"))

# Gate 3: failure semantics.  The candidate must reject the same malformed
# supplied permutations as the oracle.  Invalid block size is a candidate
# self-check because the oracle intentionally ignores that private knob.
message("[stage2d2b] failure semantics and input immutability")
bad_tree <- stage2d2_make_tree(50L, "balanced")
bad_f <- list(
  tree = bad_tree, compiled = compile_tree(bad_tree),
  X = stage2d2_make_matrix(bad_tree, 1L, "normal", 301050L),
  n = 50L, shape = "balanced", traits = 1L
)
good <- stage2d2b_make_permutations(50L, 2L, "identity", seed = 301051L)
bad_perms <- list(
  duplicate = good,
  out_of_range = good,
  wrong_dimension = good[, -1L, drop = FALSE]
)
bad_perms$duplicate[1L, 2L] <- bad_perms$duplicate[1L, 1L]
bad_perms$out_of_range[1L, 1L] <- 0L
for (bad_name in names(bad_perms)) {
  old <- stage2d2b_capture(function() stage2d2b_call_engine(
    hooks$oracle, bad_f$compiled, bad_f$X, 2L,
    permutations = bad_perms[[bad_name]], block_size = 2L,
    n_threads = 1L, trait_chunk = 1L, include_observed = FALSE,
    return_sim = TRUE
  ))
  candidate <- stage2d2b_capture(function() stage2d2b_call_engine(
    hooks$candidate, bad_f$compiled, bad_f$X, 2L,
    permutations = bad_perms[[bad_name]], block_size = 2L,
    n_threads = 1L, trait_chunk = 1L, include_observed = FALSE,
    return_sim = TRUE
  ))
  pass <- old$status == "error" && candidate$status == "error" &&
    identical(old$error_class, candidate$error_class) &&
    identical(old$warning_classes, candidate$warning_classes) &&
    identical(old$warning_text, candidate$warning_text) &&
    identical(old$error_message, candidate$error_message)
  failure_rows[[length(failure_rows) + 1L]] <- data.frame(
    check = "malformed_permutation", case = bad_name,
    oracle_status = old$status, candidate_status = candidate$status,
    oracle_error_class = old$error_class,
    candidate_error_class = candidate$error_class,
    oracle_error = old$error_message, candidate_error = candidate$error_message,
    pass = isTRUE(pass), stringsAsFactors = FALSE
  )
}

nonfinite <- bad_f$X
nonfinite[1L, 1L] <- NA_real_
old <- stage2d2b_capture(function() stage2d2b_call_engine(
  hooks$oracle, bad_f$compiled, nonfinite, 2L, permutations = good,
  block_size = 2L, n_threads = 1L, trait_chunk = 1L,
  include_observed = FALSE, return_sim = TRUE
))
candidate <- stage2d2b_capture(function() stage2d2b_call_engine(
  hooks$candidate, bad_f$compiled, nonfinite, 2L, permutations = good,
  block_size = 2L, n_threads = 1L, trait_chunk = 1L,
  include_observed = FALSE, return_sim = TRUE
))
nonfinite_pass <- old$status == "error" && candidate$status == "error" &&
  identical(old$error_class, candidate$error_class) &&
  identical(old$warning_classes, candidate$warning_classes) &&
  identical(old$error_message, candidate$error_message)
failure_rows[[length(failure_rows) + 1L]] <- data.frame(
  check = "nonfinite_trait", case = "NA", oracle_status = old$status,
  candidate_status = candidate$status, oracle_error_class = old$error_class,
  candidate_error_class = candidate$error_class,
  oracle_error = old$error_message, candidate_error = candidate$error_message,
  pass = isTRUE(nonfinite_pass), stringsAsFactors = FALSE
)

bad_block <- stage2d2b_capture(function() stage2d2b_call_engine(
  hooks$candidate, bad_f$compiled, bad_f$X, 2L, permutations = good,
  block_size = 3L, n_threads = 1L, trait_chunk = 1L,
  include_observed = FALSE, return_sim = TRUE
))
failure_rows[[length(failure_rows) + 1L]] <- data.frame(
  check = "invalid_block_size", case = "3", oracle_status = NA_character_,
  candidate_status = bad_block$status, oracle_error_class = NA_character_,
  candidate_error_class = bad_block$error_class, oracle_error = NA_character_,
  candidate_error = bad_block$error_message,
  pass = bad_block$status == "error", stringsAsFactors = FALSE
)

# Input immutability is tested independently of the oracle comparison.  The
# object snapshots include edge, edge.length, labels, Nnode, trait values, and
# the supplied permutation matrix.
immutable_tree <- stage2d2_make_tree(50L, "random")
immutable_X <- stage2d2_make_matrix(immutable_tree, 8L, "normal", 302050L)
immutable_perms <- stage2d2b_make_permutations(50L, 7L, "adversarial", 302051L)
before_tree <- immutable_tree
before_X <- immutable_X
before_perms <- immutable_perms
invisible(stage2d2b_capture(function() stage2d2b_call_engine(
  hooks$candidate, compile_tree(immutable_tree), immutable_X, 7L,
  permutations = immutable_perms, block_size = 4L, n_threads = 1L,
  trait_chunk = 8L, include_observed = FALSE, return_sim = TRUE
)))
immutability_rows[[1L]] <- data.frame(
  check = "input_immutability",
  tree_unchanged = identical(before_tree, immutable_tree),
  X_unchanged = identical(before_X, immutable_X),
  permutations_unchanged = identical(before_perms, immutable_perms),
  pass = identical(before_tree, immutable_tree) &&
    identical(before_X, immutable_X) && identical(before_perms, immutable_perms),
  stringsAsFactors = FALSE
)

# A bounded batch guard exercises multiple trait widths without introducing
# stochastic-generation cost.  It also ensures ordered nulls are retained for
# every trait column, not just the one-trait heavy fixture.
batch_rows <- list()
for (p in c(1L, 8L, 32L, 100L)) {
  batch_f <- fixture(50L, "balanced", p, 303050L + p)
  perms <- stage2d2b_make_permutations(50L, 7L, "fixed_random", 303051L + p)
  one <- run_pair("batch_controlled", batch_f, batch_f$X, 7L, perms,
                  "fixed_random", 8L, 1L, 303052L + p,
                  include_observed = FALSE, require_identity = TRUE)
  batch_rows[[length(batch_rows) + 1L]] <- data.frame(
    traits = p, pass = isTRUE(one$pass), stringsAsFactors = FALSE
  )
}

failures <- stage2d2b_rbind(failure_rows)
immutability <- stage2d2b_rbind(immutability_rows)
batch <- stage2d2b_rbind(batch_rows)
summary <- stage2d2b_rbind(summary_rows)
replicates <- stage2d2b_rbind(replicate_rows)
counters <- stage2d2b_rbind(counter_rows)

stage2d2b_write_csv(summary, file.path(args$out,
                                       "stage2d2b_correctness_summary.csv"))
stage2d2b_write_csv(replicates, file.path(args$out,
                                          "stage2d2b_correctness_replicates.csv"))
stage2d2b_write_csv(counters, file.path(args$out,
                                        "stage2d2b_correctness_counters.csv"))
stage2d2b_write_csv(failures, file.path(args$out,
                                        "stage2d2b_correctness_failures.csv"))
stage2d2b_write_csv(immutability, file.path(args$out,
                                            "stage2d2b_correctness_immutability.csv"))
stage2d2b_write_csv(batch, file.path(args$out,
                                     "stage2d2b_batch_guard.csv"))
stage2d2b_write_csv(rng_summary, file.path(args$out,
                                           "stage2d2b_rng_replay.csv"))

replicate_pass <- nrow(replicates) > 0L &&
  all(replicates$bitwise_exact & replicates$success_exact)
rng_pass <- nrow(rng_summary) > 0L && all(rng_summary$pass)
failure_pass <- nrow(failures) > 0L && all(failures$pass)
immutability_pass <- nrow(immutability) > 0L && all(immutability$pass)
batch_pass <- nrow(batch) > 0L && all(batch$pass)
counter_pass <- nrow(counters) > 0L && all(counters$counter_gate_pass)
overall <- nrow(summary) > 0L && all(summary$pass) && replicate_pass &&
  rng_pass && failure_pass && immutability_pass && batch_pass && counter_pass

provenance <- data.frame(
  stage = "fastphylosig 0.2.0 Stage 2D2B true blocked evaluator correctness",
  mode = args$mode, source_commit = stage2d2b_git_value(args$repo),
  package_version = tryCatch(as.character(utils::packageVersion(
    "fastphylosig")), error = function(e) NA_character_),
  R_version = R.version.string, platform = R.version$platform,
  os = Sys.info()[["sysname"]], machine = Sys.info()[["machine"]],
  prototype_cpp = hooks$cpp_file,
  prototype_cpp_sha256 = stage2d2b_source_sha256(hooks$cpp_file),
  oracle_export = hooks$oracle_name, candidate_export = hooks$candidate_name,
  same_translation_unit = hooks$same_translation_unit, candidate_gate = hooks$gate,
  shapes = paste(shapes, collapse = ","), n_grid = paste(n_grid, collapse = ","),
  patterns = paste(patterns, collapse = ","), blocks = paste(blocks, collapse = ","),
  nsim_controlled = nsim_controlled, replay_nsim = replay_nsim,
  replay_threads = paste(replay_threads, collapse = ","),
  controlled_rows = nrow(controlled), controlled_failures = sum(!controlled$pass),
  replicate_rows = nrow(replicates), rng_rows = nrow(rng_summary),
  failure_rows = nrow(failures), immutability_rows = nrow(immutability),
  batch_rows = nrow(batch),
  old_compute_one_call_contract = "candidate counter == 0",
  oracle_call_contract = "candidate counter == 0",
  candidate_compute_contract = "candidate counter > 0",
  production_code_changed = "NO", public_api_changed = "NO",
  estimator_changed = "NO", RNG_policy_changed = "NO",
  stringsAsFactors = FALSE
)
stage2d2b_write_csv(provenance, file.path(args$out,
                                          "stage2d2b_correctness_provenance.csv"))

status <- if (overall) "PASS" else "FAIL"
reason <- if (overall) "all Stage 2D2B correctness/provenance gates passed" else
  paste0("controlled=", controlled_pass, "; replicate=", replicate_pass,
         "; rng=", rng_pass, "; failures=", failure_pass,
         "; immutability=", immutability_pass, "; batch=", batch_pass,
         "; counters=", counter_pass)
write_status(status, reason, list(
  controlled_rows = nrow(controlled), controlled_failures = sum(!controlled$pass),
  replicate_rows = nrow(replicates), rng_rows = nrow(rng_summary),
  failure_rows = nrow(failures), immutability_rows = nrow(immutability),
  batch_rows = nrow(batch), counter_rows = nrow(counters),
  candidate_old_compute_one_calls_zero = counter_pass &&
    all(counters$old_compute_one_calls == 0),
  candidate_oracle_call_count_zero = counter_pass &&
    all(counters$oracle_call_count == 0),
  candidate_compute_count_positive = counter_pass &&
    all(counters$candidate_compute_count > 0),
  controlled_bitwise_exact = controlled_pass,
  rng_replay_exact = rng_pass
))

if (!overall) stop("Stage 2D2B correctness/provenance gate failed.", call. = FALSE)
