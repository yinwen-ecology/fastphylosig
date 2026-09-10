#!/usr/bin/env Rscript

# Opt-in Stage 2D2 exactness runner.
#
# This runner performs bounded correctness checks only.  It does not run the
# Stage 2D2 formal timing grid, bind a private prototype to production, or
# modify package files.  Set one of:
#
#   FASTPHYLOSIG_STAGE2D2_PROTO_CPP=/ascii/path/to/one_translation_unit.cpp
#   FASTPHYLOSIG_STAGE2D2_HOOK_FILE=/ascii/path/to/hook.R
#
# The hook must provide an exact same-translation-unit oracle plus Candidate A
# and Candidate B.  A sourceCpp function named
# `stage2d2_k_permutation_prototype()` may expose all three through
# `mode = "oracle"`, `"A"`, and `"B"`; separate exported functions are also
# accepted.  Missing or incomplete hooks fail closed with NOT_RUN.
# The public `test = FALSE` and batch paths are independent performance guards;
# they are deliberately not substituted by this internal exactness runner and
# are recorded as NOT_RUN_HERE.

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- all_args[grepl("^--file=", all_args)]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else
  file.path(getwd(), "run_stage2d2_correctness.R")
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
source(file.path(script_dir, "stage2d2_common.R"), local = TRUE)

args <- stage2d2_parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(args$out, "stage2d2_correctness_status.csv")

write_not_run <- function(reason) {
  stage2d2_write_csv(data.frame(
    status = "NOT_RUN", checks = 0L, failures = 0L,
    reason = as.character(reason), production_code_changed = "NO",
    stringsAsFactors = FALSE
  ), status_path)
  stop(reason, call. = FALSE)
}

Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")

hooks <- tryCatch(
  stage2d2_load_hooks(args$repo),
  error = function(e) list(ready = FALSE, reason = conditionMessage(e))
)
if (!isTRUE(hooks$ready)) write_not_run(hooks$reason)

ns <- tryCatch(
  stage2d2_load_package(args$repo),
  error = function(e) {
    write_not_run(paste0("Cannot load audited fastphylosig package: ",
                         conditionMessage(e)))
  }
)
if (!is.environment(ns)) write_not_run("fastphylosig namespace is unavailable.")

prov <- tryCatch(stage2d2_provenance(args$repo, args, hooks),
                 error = function(e) data.frame(
                   stage = "Stage 2D2 exactness correctness",
                   level = args$level, source_commit = "<unavailable>",
                   provenance_error = conditionMessage(e),
                   stringsAsFactors = FALSE
                 ))
stage2d2_write_csv(prov, file.path(args$out,
                                   "stage2d2_correctness_provenance.csv"))
stage2d2_write_csv(stage2d2_memory_contract(), file.path(
  args$out, "stage2d2_memory_contract.csv"
))

oracle <- hooks$oracle
if (identical(hooks$candidate_mode, "B")) {
  candidates <- list(B = hooks$candidate_b)
} else if (identical(hooks$candidate_mode, "C")) {
  candidates <- list(C = hooks$candidate_c)
} else {
  candidates <- list(A = hooks$candidate_a)
}
if (!identical(hooks$candidate_mode, "B") &&
    isTRUE(hooks$candidate_b_authorized)) {
  if (!is.function(hooks$candidate_b)) {
    write_not_run("Candidate B was authorized but no Candidate B hook exists.")
  }
  candidates$B <- hooks$candidate_b
}
compiled_tree <- function(tree) stage2d2_compile_tree(tree, ns)

summary_rows <- list()
replicate_rows <- list()
input_rows <- list()
failure_rows <- list()

config_row <- function(check, candidate, shape, n, nsim, seed, n_threads,
                       chunk, trait_kind, mask_kind, p, pass, detail,
                       cmp = NULL, accounting = NA) {
  max_abs <- if (!is.null(cmp) && !is.null(cmp$result)) cmp$result$max_abs else
    c(K = NA_real_, P = NA_real_, MCSE_P = NA_real_,
      exceedance = NA_real_, sim_K = NA_real_)
  data.frame(
    check = check, candidate = candidate, shape = shape, n = as.integer(n),
    traits = as.integer(p), nsim = as.integer(nsim),
    seed = if (is.null(seed)) NA_integer_ else as.integer(seed),
    n_threads = as.integer(n_threads), simulation_chunk = as.integer(chunk),
    trait_kind = trait_kind, mask_kind = mask_kind,
    pass = isTRUE(pass), accounting_pass = accounting,
    max_abs_K = unname(max_abs[["K"]]),
    max_abs_P = unname(max_abs[["P"]]),
    max_abs_MCSE_P = unname(max_abs[["MCSE_P"]]),
    max_abs_exceedance = unname(max_abs[["exceedance"]]),
    max_abs_sim_K = unname(max_abs[["sim_K"]]),
    detail = as.character(detail), stringsAsFactors = FALSE
  )
}

accounting_pass <- function(capture, nsim, p) {
  if (capture$status != "ok") return(NA)
  value <- stage2d2_normalize_result(capture$value, nsim, p)
  if (is.null(value) || is.null(value$P) || is.null(value$exceedance) ||
      is.null(value$requested) || is.null(value$successful) ||
      is.null(value$failed) || is.null(value$n_randomizations) ||
      is.null(value$include_observed) || length(value$P) != p ||
      length(value$exceedance) != p || length(value$requested) != 1L ||
      length(value$n_randomizations) != 1L ||
      length(value$successful) %in% c(0L) ||
      length(value$failed) %in% c(0L)) {
    return(FALSE)
  }
  expected <- as.numeric(value$exceedance) / as.numeric(value$requested)
  count_len <- max(length(value$successful), length(value$failed))
  successful <- rep(as.numeric(value$successful), length.out = count_len)
  failed <- rep(as.numeric(value$failed), length.out = count_len)
  count_pass <- stage2d2_exact_equal(successful + failed,
                                     rep(as.numeric(value$requested),
                                         length.out = count_len))
  random_expected <- if (isTRUE(value$include_observed))
    as.numeric(value$requested) - 1 else as.numeric(value$requested)
  random_pass <- stage2d2_exact_equal(as.numeric(value$n_randomizations),
                                      random_expected)
  p_pass <- if (all(is.na(value$P))) TRUE else
    stage2d2_exact_equal(value$P, expected)
  isTRUE(count_pass) && isTRUE(random_pass) && isTRUE(p_pass)
}

run_pair <- function(check, candidate_name, candidate_fun, tree, X, nsim,
                     permutations = NULL, seed = NULL, n_threads = 1L,
                     chunk = 128L, trait_kind = "normal", mask_kind = "none",
                     require_payload = TRUE, shape = "unknown") {
  compiled <- compiled_tree(tree)
  include_observed <- is.null(permutations)
  old <- stage2d2_capture(function() stage2d2_call_engine(
    oracle, compiled, X, nsim, permutations = permutations,
    trait_chunk = if (ncol(X) == 1L) 1L else 8L,
    return_sim = TRUE, include_observed = include_observed,
    n_threads = n_threads, simulation_chunk = chunk, seed = seed
  ))
  new <- stage2d2_capture(function() stage2d2_call_engine(
    candidate_fun, compiled, X, nsim, permutations = permutations,
    trait_chunk = if (ncol(X) == 1L) 1L else 8L,
    return_sim = TRUE, include_observed = include_observed,
    n_threads = n_threads, simulation_chunk = chunk, seed = seed
  ))
  cmp <- stage2d2_compare_captures(old, new, nsim, ncol(X))
  old_valid <- stage2d2_validate_payload(old, nsim, ncol(X))
  new_valid <- stage2d2_validate_payload(new, nsim, ncol(X))
  payload_pass <- if (require_payload) old_valid && new_valid else TRUE
  acc_old <- accounting_pass(old, nsim, ncol(X))
  acc_new <- accounting_pass(new, nsim, ncol(X))
  acc_pass <- if (require_payload) isTRUE(acc_old) && isTRUE(acc_new) else NA
  pass <- isTRUE(cmp$pass) && isTRUE(payload_pass) &&
    (is.na(acc_pass) || isTRUE(acc_pass))
  detail <- if (pass) "exact oracle/candidate payload and accounting" else
    paste0("comparison=", cmp$pass, "; payload=", payload_pass,
           "; accounting=", acc_old, "/", acc_new,
           "; warning=", cmp$warning_pass,
           "; errors=", old$error_class, "/", new$error_class)
  config <- list(
    check = check, candidate = candidate_name, shape = shape, n = nrow(X),
    nsim = nsim, seed = seed, n_threads = n_threads, chunk = chunk,
    trait_kind = trait_kind, mask_kind = mask_kind
  )
  summary_rows[[length(summary_rows) + 1L]] <<- config_row(
    check, candidate_name, shape, config$n, nsim, seed, n_threads, chunk,
    trait_kind, mask_kind, ncol(X), pass, detail, cmp, acc_pass
  )
  if (old$status == "ok" && new$status == "ok") {
    replicate_rows[[length(replicate_rows) + 1L]] <<-
      stage2d2_replicate_rows(cmp, config)
  }
  list(old = old, new = new, comparison = cmp, pass = pass)
}

run_exact_case <- function(check, tree, shape, X, nsim, permutations = NULL,
                           seed = NULL, n_threads = 1L, chunk = 2L,
                           trait_kind = "normal", mask_kind = "none",
                           require_payload = TRUE) {
  for (candidate_name in names(candidates)) {
    run_pair(check, candidate_name, candidates[[candidate_name]], tree, X,
             nsim, permutations, seed, n_threads, chunk, trait_kind,
             mask_kind, require_payload, shape)
  }
}

# The smoke level is deliberately small but covers every requested fixture
# family.  --correctness adds one larger bounded tree and more replay edges;
# neither mode is the Stage 2D2 formal timing grid.
if (isTRUE(args$smoke)) {
  n_grid <- c(2L, 7L, 17L)
  controlled_nsim <- c(1L, 2L, 7L)
  rng_nsim <- c(1L, 7L)
  seed_grid <- 20260909L
} else {
  n_grid <- c(2L, 7L, 17L, 50L)
  controlled_nsim <- c(1L, 2L, 7L, 19L)
  rng_nsim <- c(2L, 7L, 19L)
  seed_grid <- c(20260909L, 99173L)
}
shapes <- c("balanced", "random", "pectinate")
threads <- c(1L, 2L)

# Controlled exactness and batch trait coverage.
for (shape in shapes) for (n in n_grid) {
  tree <- stage2d2_make_tree(n, shape)
  for (p in c(1L, 8L, 32L, 100L)) {
    X <- stage2d2_make_matrix(tree, traits = p, kind = "normal",
                              seed = 203000L + n + p)
    for (nsim in controlled_nsim) for (n_threads in threads) {
      perms <- stage2d2_make_permutations(n, nsim,
                                          seed = 204000L + n + nsim)
      run_exact_case("controlled_permutation", tree, shape, X, nsim,
                     permutations = perms, n_threads = n_threads,
                     chunk = min(7L, nsim))
    }
  }
}

# Numerical edge traits, including constant-result failure/NA semantics.
edge_kinds <- c("large_offset", "near_constant", "constant")
for (shape in shapes) {
  tree <- stage2d2_make_tree(17L, shape)
  for (kind in edge_kinds) {
    X <- stage2d2_make_matrix(tree, traits = 1L, kind = kind,
                              seed = 205000L)
    perms <- stage2d2_make_permutations(17L, 7L, seed = 205017L)
    run_exact_case("edge_trait", tree, shape, X, 7L, permutations = perms,
                   n_threads = 1L, chunk = 2L, trait_kind = kind)
  }
}

# Two-tip boundary.
two_tip <- stage2d2_make_tree(2L, "balanced")
two_tip_X <- stage2d2_make_matrix(two_tip, traits = 2L, kind = "normal",
                                  seed = 206002L)
two_tip_perms <- stage2d2_make_permutations(2L, 2L, seed = 206003L)
run_exact_case("two_tip", two_tip, "two_tip", two_tip_X, 2L,
               permutations = two_tip_perms, n_threads = 1L, chunk = 1L)

# NA masks are represented by their retained, complete-case subsets.  Shared
# masks exercise a batch route; multiple masks exercise separate retained
# subsets and therefore more than one compiled tree.
na_tree <- stage2d2_make_tree(17L, "balanced")
for (mask_kind in c("shared", "multiple")) {
  mask_cases <- stage2d2_masked_cases(na_tree, traits = 8L, kind = mask_kind)
  if (!length(mask_cases)) stop("NA mask fixture retained no valid subset.",
                               call. = FALSE)
  for (case in mask_cases) {
    for (n_threads in threads) {
      nsim <- 7L
      perms <- stage2d2_make_permutations(nrow(case$X), nsim,
                                          seed = 207000L + length(case$X))
      run_exact_case("na_retained_subset", case$tree, "balanced", case$X,
                     nsim, permutations = perms, n_threads = n_threads,
                     chunk = 2L, mask_kind = mask_kind)
    }
  }
}

# Internal RNG replay: same seed/configuration, chunk boundaries, both thread
# settings, and a second invocation of each candidate.  Permutations remain
# serially generated by the frozen RNG contract; only evaluator layout is
# under test.
for (shape in shapes) for (n in n_grid) {
  tree <- stage2d2_make_tree(n, shape)
  X <- stage2d2_make_matrix(tree, traits = 1L, kind = "normal",
                            seed = 208000L + n)
  for (nsim in rng_nsim) for (seed in if (length(seed_grid) > 1L)
    seed_grid else list(seed_grid)) for (n_threads in threads) {
    chunks <- unique(pmax(1L, c(1L, 2L, min(7L, nsim), nsim)))
    replay_reference <- list()
    replay_pass <- list()
    for (chunk in chunks) {
      for (candidate_name in names(candidates)) {
        candidate_fun <- candidates[[candidate_name]]
        pair <- run_pair(
          "rng_oracle_replay", candidate_name, candidate_fun, tree, X, nsim,
          permutations = NULL, seed = seed, n_threads = n_threads,
          chunk = chunk, shape = shape
        )
        key <- candidate_name
        current <- if (pair$new$status == "ok") pair$new$value else NULL
        if (is.null(replay_reference[[key]])) {
          replay_reference[[key]] <- current
          replay_pass[[key]] <- pair$new$status == "ok"
        } else {
          replay_pass[[key]] <- isTRUE(replay_pass[[key]]) &&
            !is.null(current) && stage2d2_compare_results(
              replay_reference[[key]], current, nsim, ncol(X)
            )$pass
        }
        # Repeat the same call after restoring the seed; this is an explicit
        # same-seed replay check independent of the oracle comparison.
        repeated <- stage2d2_capture(function() stage2d2_call_engine(
          candidate_fun, compiled_tree(tree), X, nsim,
          permutations = NULL, trait_chunk = 1L, return_sim = TRUE,
          include_observed = TRUE, n_threads = n_threads,
          simulation_chunk = chunk, seed = seed
        ))
        same <- pair$new$status == "ok" && repeated$status == "ok" &&
          stage2d2_compare_results(pair$new$value, repeated$value,
                                   nsim, ncol(X))$pass
        summary_rows[[length(summary_rows) + 1L]] <- config_row(
          "rng_same_seed_replay", candidate_name, shape, n, nsim, seed,
          n_threads, chunk, "normal", "none", ncol(X), same,
          if (same) "candidate repeated under restored seed" else
            "candidate changed under restored seed", NULL, NA
        )
      }
    }
    # A candidate's ordered null and all accounting fields must not depend on
    # the chunk boundary under the supported fixed configuration.
    for (candidate_name in names(candidates)) {
      pass <- isTRUE(replay_pass[[candidate_name]]) &&
        length(chunks) > 1L
      summary_rows[[length(summary_rows) + 1L]] <- config_row(
        "rng_chunk_boundary_replay", candidate_name, shape, n, nsim, seed,
        n_threads, NA_integer_, "normal", "none", ncol(X), pass,
        if (pass) "candidate chunk settings replay exactly" else
          "candidate chunk settings differ", NULL, NA
      )
    }
  }
}

# Invalid controlled permutations must retain the same failure class/status.
bad_tree <- stage2d2_make_tree(7L, "balanced")
bad_X <- stage2d2_make_matrix(bad_tree, traits = 1L, kind = "normal",
                              seed = 209007L)
good_perms <- stage2d2_make_permutations(7L, 2L, seed = 209008L)
bad_perms <- list(
  duplicate = good_perms,
  out_of_range = good_perms,
  wrong_dimension = good_perms[, -1L, drop = FALSE]
)
bad_perms$duplicate[1L, 2L] <- bad_perms$duplicate[1L, 1L]
bad_perms$out_of_range[1L, 1L] <- 0L
for (bad_name in names(bad_perms)) {
  bad <- bad_perms[[bad_name]]
  for (candidate_name in names(candidates)) {
    candidate_fun <- candidates[[candidate_name]]
    compiled <- compiled_tree(bad_tree)
    old <- stage2d2_capture(function() stage2d2_call_engine(
      oracle, compiled, bad_X, 2L, permutations = bad,
      trait_chunk = 1L, return_sim = TRUE, include_observed = FALSE,
      n_threads = 1L, simulation_chunk = 2L
    ))
    new <- stage2d2_capture(function() stage2d2_call_engine(
      candidate_fun, compiled, bad_X, 2L, permutations = bad,
      trait_chunk = 1L, return_sim = TRUE, include_observed = FALSE,
      n_threads = 1L, simulation_chunk = 2L
    ))
    cmp <- stage2d2_compare_captures(old, new, 2L, 1L)
    pass <- old$status == "error" && new$status == "error" &&
      cmp$status_pass && cmp$error_pass && cmp$warning_pass
    failure_rows[[length(failure_rows) + 1L]] <- data.frame(
      check = "malformed_permutation", failure = bad_name,
      candidate = candidate_name, old_status = old$status,
      new_status = new$status, old_error_class = old$error_class,
      new_error_class = new$error_class, old_error = old$error_message,
      new_error = new$error_message, pass = pass,
      stringsAsFactors = FALSE
    )
  }
}

# Input immutability is a separate gate for ordinary controlled data.
immutable_tree <- stage2d2_make_tree(17L, "random")
immutable_X <- stage2d2_make_matrix(immutable_tree, traits = 8L,
                                    kind = "normal", seed = 210017L)
immutable_perms <- stage2d2_make_permutations(17L, 7L, seed = 210018L)
for (candidate_name in names(candidates)) {
  before_tree <- immutable_tree
  before_X <- immutable_X
  before_perms <- immutable_perms
  compiled <- compiled_tree(immutable_tree)
  stage2d2_capture(function() stage2d2_call_engine(
    candidates[[candidate_name]], compiled, immutable_X, 7L,
    permutations = immutable_perms, trait_chunk = 8L, return_sim = TRUE,
    include_observed = FALSE, n_threads = 1L, simulation_chunk = 2L
  ))
  pass <- identical(before_tree, immutable_tree) &&
    identical(before_X, immutable_X) && identical(before_perms, immutable_perms)
  input_rows[[length(input_rows) + 1L]] <- data.frame(
    check = "input_immutability", candidate = candidate_name,
    tree_unchanged = identical(before_tree, immutable_tree),
    X_unchanged = identical(before_X, immutable_X),
    permutations_unchanged = identical(before_perms, immutable_perms),
    pass = pass, stringsAsFactors = FALSE
  )
}

summary <- stage2d2_rbind(summary_rows)
replicates <- stage2d2_rbind(replicate_rows)
inputs <- stage2d2_rbind(input_rows)
failures <- stage2d2_rbind(failure_rows)
stage2d2_write_csv(summary, file.path(args$out,
                                      "stage2d2_correctness_summary.csv"))
stage2d2_write_csv(replicates, file.path(args$out,
                                         "stage2d2_correctness_replicates.csv"))
stage2d2_write_csv(inputs, file.path(args$out,
                                     "stage2d2_correctness_immutability.csv"))
stage2d2_write_csv(failures, file.path(args$out,
                                      "stage2d2_correctness_failures.csv"))

summary_pass <- nrow(summary) > 0L && all(summary$pass)
replicate_pass <- nrow(replicates) == 0L || all(replicates$null_exact &
                                                replicates$success_exact)
input_pass <- nrow(inputs) > 0L && all(inputs$pass)
failure_pass <- nrow(failures) > 0L && all(failures$pass)
overall <- summary_pass && replicate_pass && input_pass && failure_pass
status <- if (overall) "PASS" else "FAIL"
reason <- if (overall) "all bounded exactness gates passed" else paste0(
  "summary=", summary_pass, "; replicate=", replicate_pass,
  "; input_immutability=", input_pass, "; failure_semantics=", failure_pass
)
stage2d2_write_csv(data.frame(
  status = status, checks = nrow(summary), failures = sum(!summary$pass),
  replicate_rows = nrow(replicates), immutability_rows = nrow(inputs),
  malformed_rows = nrow(failures), same_translation_unit_oracle =
    isTRUE(hooks$same_translation_unit), reason = reason,
  candidate_a_status = "RUN",
  candidate_b_gate = hooks$candidate_b_gate %||% "NOT_AUTHORIZED",
  candidate_b_status = hooks$candidate_b_status %||% "NOT_AUTHORIZED",
  public_test_false_batch_guard = "INDEPENDENT_NOT_RUN_HERE",
  formal_performance_grid = "NOT_RUN",
  production_code_changed = "NO", stringsAsFactors = FALSE
), status_path)

if (!overall) stop("Stage 2D2 bounded correctness audit failed.", call. = FALSE)
