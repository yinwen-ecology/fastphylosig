#!/usr/bin/env Rscript

# Opt-in Stage 2D2 workspace-safety audit for Candidate A.
#
# This runner is audit infrastructure only.  It loads one sourceCpp
# translation unit supplied explicitly by the caller, and never changes
# package sources, package bindings, tests, or the public API.

# Set the locale before reading the Rscript --file path.  This keeps the
# ASCII-only runner usable when its parent workspace path is not ASCII.
Sys.setenv(LC_ALL = "C", LANG = "C", LANGUAGE = "C")

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- all_args[grepl("^--file=", all_args)]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else
  file.path(getwd(), "run_stage2d2_workspace_safety.R")
script_dir <- tryCatch(
  dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE)),
  error = function(e) {
    # Windows may pass a non-ASCII --file path in a non-UTF-8 code page.
    # When launched from the repository, the working directory is reliable.
    cwd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    candidates <- unique(c(
      file.path(cwd, "benchmarks", "stage2d2"),
      file.path(cwd, "..", "benchmarks", "stage2d2"),
      file.path(cwd, "..", "..", "benchmarks", "stage2d2")
    ))
    candidates <- candidates[file.exists(file.path(
      candidates, "stage2d2_common.R"
    ))]
    if (length(candidates)) candidates[[1L]] else stop(
      "Cannot locate stage2d2_common.R after --file path conversion failed.",
      call. = FALSE
    )
  }
)
common_file <- file.path(script_dir, "stage2d2_common.R")
if (!file.exists(common_file)) {
  stop("stage2d2_common.R is missing beside the runner.", call. = FALSE)
}
source(common_file, local = TRUE)

parse_cli <- function(values) {
  smoke <- "--smoke" %in% values
  correctness <- "--correctness" %in% values
  if (isTRUE(smoke) == isTRUE(correctness)) {
    stop("Use exactly one of --smoke or --correctness.", call. = FALSE)
  }
  positional <- values[!values %in% c("--smoke", "--correctness")]
  if (length(positional) > 2L) {
    stop("usage: run_stage2d2_workspace_safety.R [--smoke|--correctness] ",
         "[repo] [out].", call. = FALSE)
  }
  default_repo <- normalizePath(file.path(script_dir, "..", ".."),
                                winslash = "/", mustWork = FALSE)
  repo <- if (length(positional)) positional[[1L]] else default_repo
  env_out <- stage2d2_env("FASTPHYLOSIG_STAGE2D2_WORKSPACE_SAFETY_OUT")
  out <- if (length(positional) >= 2L) positional[[2L]] else if (nzchar(env_out))
    env_out else file.path(repo, "benchmarks", "stage2d2", "results",
                           "workspace-safety")
  list(
    level = if (isTRUE(correctness)) "targeted" else "smoke",
    smoke = isTRUE(smoke), correctness = isTRUE(correctness),
    repo = normalizePath(repo, winslash = "/", mustWork = FALSE),
    out = normalizePath(out, winslash = "/", mustWork = FALSE)
  )
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(args$out, "stage2d2_workspace_safety_status.csv")
evidence_path <- file.path(args$out, "stage2d2_workspace_safety_evidence.csv")

write_status <- function(status, reason, checks = 0L, failures = 0L,
                         source_cpp = NA_character_,
                         same_translation_unit = FALSE, extra = list()) {
  base <- list(
    status = as.character(status), checks = as.integer(checks),
    failures = as.integer(failures), reason = as.character(reason),
    source_cpp_opt_in = FALSE, source_cpp = source_cpp,
    same_translation_unit = same_translation_unit, candidate_mode = "A",
    candidate_debug_workspace = TRUE, shapes = NA_character_,
    trait_kinds = NA_character_, batch8_cells = 0L,
    ordered_sim_exact = NA, K_exact = NA, P_exact = NA,
    MCSE_P_exact = NA, accounting_exact = NA,
    workspace_reuse_positive = NA, debug_checks_positive = NA,
    canary_checks_positive = NA, stale_state_failures = NA_real_,
    canary_failures = NA_real_, production_code_changed = "NO"
  )
  for (nm in names(extra)) base[[nm]] <- extra[[nm]]
  stage2d2_write_csv(as.data.frame(base, stringsAsFactors = FALSE),
                     status_path)
  invisible(status_path)
}

fail_early <- function(reason, source_cpp = NA_character_,
                       same_translation_unit = FALSE) {
  write_status("NOT_RUN", reason, source_cpp = source_cpp,
               same_translation_unit = same_translation_unit)
  stop(reason, call. = FALSE)
}

# A prototype is never selected implicitly.  The target-specific variable is
# preferred; the common Stage 2D2 variable remains an explicit opt-in alias.
proto_cpp <- stage2d2_env("FASTPHYLOSIG_STAGE2D2_WORKSPACE_SAFETY_PROTO_CPP")
if (!nzchar(proto_cpp)) {
  proto_cpp <- stage2d2_env("FASTPHYLOSIG_STAGE2D2_PROTO_CPP")
}
if (!nzchar(proto_cpp)) {
  fail_early(paste0(
    "Workspace-safety audit is opt-in: set ",
    "FASTPHYLOSIG_STAGE2D2_WORKSPACE_SAFETY_PROTO_CPP or ",
    "FASTPHYLOSIG_STAGE2D2_PROTO_CPP to one sourceCpp translation unit."
  ))
}
proto_cpp <- tryCatch(normalizePath(proto_cpp, winslash = "/", mustWork = TRUE),
                      error = function(e) fail_early(
                        paste0("Prototype source is unavailable: ",
                               conditionMessage(e))))
if (!requireNamespace("Rcpp", quietly = TRUE)) {
  fail_early("Rcpp is required for the explicit sourceCpp prototype.",
             source_cpp = proto_cpp)
}

hook_env <- new.env(parent = globalenv())
exports_before <- ls(envir = hook_env, all.names = TRUE)
source_error <- NULL
tryCatch(
  Rcpp::sourceCpp(proto_cpp, env = hook_env, rebuild = TRUE,
                  showOutput = FALSE, verbose = FALSE),
  error = function(e) source_error <<- e
)
if (!is.null(source_error)) {
  fail_early(paste0("sourceCpp prototype failed: ",
                    conditionMessage(source_error)), source_cpp = proto_cpp)
}
exports <- setdiff(ls(envir = hook_env, all.names = TRUE), exports_before)
generic_name <- NULL
for (nm in c("stage2d2_k_permutation_prototype",
             "stage2d2_k_null_prototype")) {
  if (nm %in% exports && exists(nm, envir = hook_env, inherits = FALSE) &&
      is.function(get(nm, envir = hook_env, inherits = FALSE))) {
    generic_name <- nm
    break
  }
}
if (is.null(generic_name)) {
  fail_early(paste0(
    "The sourceCpp translation unit must export the generic Stage 2D2 ",
    "prototype with oracle and A modes."), source_cpp = proto_cpp)
}
prototype <- get(generic_name, envir = hook_env, inherits = FALSE)
prototype_bound <- identical(
  get(generic_name, envir = hook_env, inherits = FALSE), prototype
)
prototype_formals <- tryCatch(names(formals(prototype)),
                              error = function(e) character())
if (!prototype_bound || !all(c("mode", "debug_workspace", "return_ordered") %in%
                             prototype_formals)) {
  fail_early(paste0(
    "The sourceCpp prototype must expose mode, debug_workspace, and ",
    "return_ordered formals."), source_cpp = proto_cpp,
    same_translation_unit = prototype_bound)
}

ns <- tryCatch(stage2d2_load_package(args$repo), error = function(e) {
  fail_early(paste0("Cannot load audited fastphylosig package: ",
                    conditionMessage(e)), source_cpp = proto_cpp,
             same_translation_unit = prototype_bound)
})
if (!is.environment(ns)) {
  fail_early("fastphylosig namespace is unavailable.", source_cpp = proto_cpp,
             same_translation_unit = prototype_bound)
}

compile_tree <- function(tree) stage2d2_compile_tree(tree, ns)

if (isTRUE(args$smoke)) {
  n_grid <- c(2L, 7L)
  nsim_grid <- c(7L)
  chunk_grid <- c(2L)
} else {
  n_grid <- c(2L, 7L, 17L)
  nsim_grid <- c(7L)
  chunk_grid <- c(2L)
}
shapes <- c("balanced", "random", "pectinate")
trait_kinds <- c("normal", "large_offset", "near_constant")
trait_grid <- c(1L, 8L)

cell_seed <- function(kind_index, n, nsim, traits, chunk, pair = 0L) {
  value <- 20260910 + as.double(kind_index) * 100003 +
    as.double(n) * 1009 + as.double(nsim) * 9176 +
    as.double(traits) * 101 + as.double(chunk) * 17 + as.double(pair)
  as.integer((value %% 2147483646) + 1)
}

call_mode <- function(mode, compiled, X, nsim, permutations, trait_chunk,
                      simulation_chunk, seed, debug_workspace) {
  call_args <- list(
    compiled_tree = compiled, X = X, nsim = as.integer(nsim),
    permutations = permutations, trait_chunk = as.integer(trait_chunk),
    return_sim = TRUE, include_observed = FALSE, n_threads = 1L,
    simulation_chunk = as.integer(simulation_chunk), return_ordered = TRUE,
    mode = mode, debug_workspace = isTRUE(debug_workspace)
  )
  stage2d2_capture(function() stage2d2_call_with_seed(
    prototype, call_args, seed = seed
  ))
}

counter <- function(value, names, default = NA_real_) {
  if (is.null(value) || !is.list(value)) return(default)
  for (nm in names) {
    if (!is.null(value[[nm]]) && length(value[[nm]])) {
      out <- suppressWarnings(as.numeric(value[[nm]][[1L]]))
      if (length(out) && !is.na(out)) return(out[[1L]])
    }
  }
  default
}

field_equal <- function(left, right, nsim, p, field) {
  if (left$status != "ok" || right$status != "ok") return(FALSE)
  a <- stage2d2_normalize_result(left$value, nsim, p)
  b <- stage2d2_normalize_result(right$value, nsim, p)
  if (is.null(a) || is.null(b)) return(FALSE)
  stage2d2_exact_equal(a[[field]], b[[field]])
}

accounting_exact <- function(capture, nsim, p) {
  if (capture$status != "ok") return(FALSE)
  value <- stage2d2_normalize_result(capture$value, nsim, p)
  if (is.null(value) || is.null(value$P) || is.null(value$exceedance) ||
      is.null(value$requested) || is.null(value$successful) ||
      is.null(value$failed) || is.null(value$n_randomizations) ||
      is.null(value$include_observed)) return(FALSE)
  if (length(value$P) != p || length(value$exceedance) != p ||
      length(value$requested) != 1L || length(value$n_randomizations) != 1L ||
      !identical(as.logical(value$include_observed), FALSE)) return(FALSE)
  requested <- as.numeric(value$requested)
  if (!stage2d2_exact_equal(requested, as.numeric(nsim)) ||
      !stage2d2_exact_equal(as.numeric(value$n_randomizations),
                            as.numeric(nsim))) return(FALSE)
  count_len <- max(length(value$successful), length(value$failed), p)
  successful <- rep(as.numeric(value$successful), length.out = count_len)
  failed <- rep(as.numeric(value$failed), length.out = count_len)
  if (!stage2d2_exact_equal(successful + failed,
                            rep(as.numeric(nsim), length.out = count_len))) {
    return(FALSE)
  }
  expected_p <- as.numeric(value$exceedance) / as.numeric(nsim)
  stage2d2_exact_equal(as.numeric(value$P), expected_p)
}

evidence_rows <- list()
cell_index <- 0L

run_cell <- function(shape, n, trait_kind, traits, nsim, chunk,
                     kind_index) {
  tree <- stage2d2_make_tree(n, shape)
  X <- stage2d2_make_matrix(
    tree, traits = traits, kind = trait_kind,
    seed = 203000L + as.integer(n) + as.integer(traits) + kind_index
  )
  compiled <- compile_tree(tree)
  permutations <- stage2d2_make_permutations(
    n, nsim, seed = 204000L + as.integer(n) + as.integer(nsim) + kind_index
  )
  seed <- cell_seed(kind_index, n, nsim, traits, chunk)
  oracle_first <- ((kind_index + n + traits + chunk) %% 2L) == 0L
  if (isTRUE(oracle_first)) {
    oracle <- call_mode("oracle", compiled, X, nsim, permutations,
                       min(8L, traits), chunk, seed, FALSE)
    candidate <- call_mode("A", compiled, X, nsim, permutations,
                          min(8L, traits), chunk, seed, TRUE)
    call_order <- "oracle>A"
  } else {
    candidate <- call_mode("A", compiled, X, nsim, permutations,
                          min(8L, traits), chunk, seed, TRUE)
    oracle <- call_mode("oracle", compiled, X, nsim, permutations,
                       min(8L, traits), chunk, seed, FALSE)
    call_order <- "A>oracle"
  }

  p <- as.integer(traits)
  compared <- stage2d2_compare_captures(oracle, candidate, nsim, p)
  oracle_valid <- stage2d2_validate_payload(oracle, nsim, p)
  candidate_valid <- stage2d2_validate_payload(candidate, nsim, p)
  exact_fields <- c(
    ordered_sim = field_equal(oracle, candidate, nsim, p, "sim_K"),
    K = field_equal(oracle, candidate, nsim, p, "K"),
    P = field_equal(oracle, candidate, nsim, p, "P"),
    MCSE_P = field_equal(oracle, candidate, nsim, p, "MCSE_P"),
    accounting_fields = isTRUE(compared$pass) &&
      accounting_exact(oracle, nsim, p) && accounting_exact(candidate, nsim, p)
  )
  workspace_reuses <- counter(candidate$value,
                              c("stage2d2_workspace_reuses", "workspace_reuses"))
  debug_checks <- counter(candidate$value,
                          c("stage2d2_debug_checks", "debug_checks"))
  canary_checks <- counter(candidate$value,
                           c("stage2d2_canary_checks", "canary_checks"))
  stale_failures <- counter(candidate$value, c(
    "stage2d2_stale_state_failures", "stale_state_failures"
  ))
  canary_failures <- counter(candidate$value, c(
    "stage2d2_canary_failures", "canary_failures"
  ))
  counters_present <- all(is.finite(c(workspace_reuses, debug_checks,
                                      canary_checks, stale_failures,
                                      canary_failures)))
  counter_checks <- c(
    workspace_reuse_positive = counters_present && workspace_reuses > 0,
    debug_checks_positive = counters_present && debug_checks > 0,
    canary_checks_positive = counters_present && canary_checks > 0,
    stale_state_failures_zero = counters_present && stale_failures == 0,
    canary_failures_zero = counters_present && canary_failures == 0
  )
  pass <- oracle_valid && candidate_valid && all(exact_fields) &&
    all(counter_checks)
  detail <- if (pass) "exact ordered payload, accounting, and clean workspace diagnostics" else
    paste0("payload=", oracle_valid, "/", candidate_valid,
           "; exact=", paste(names(exact_fields)[!exact_fields], collapse = ","),
           "; counters=", paste(names(counter_checks)[!counter_checks],
                                 collapse = ","),
           "; oracle_status=", oracle$status,
           "; candidate_status=", candidate$status)
  data.frame(
    check = "workspace_safety", level = args$level,
    shape = shape, n = as.integer(n), traits = p, nsim = as.integer(nsim),
    trait_kind = trait_kind, permutation_mode = "controlled",
    simulation_chunk = as.integer(chunk), trait_chunk = min(8L, p),
    seed = as.integer(seed), call_order = call_order,
    oracle_mode = "oracle", candidate_mode = "A",
    candidate_debug_workspace = TRUE,
    oracle_status = oracle$status, candidate_status = candidate$status,
    oracle_payload_valid = oracle_valid,
    candidate_payload_valid = candidate_valid,
    ordered_sim_exact = exact_fields[["ordered_sim"]],
    K_exact = exact_fields[["K"]], P_exact = exact_fields[["P"]],
    MCSE_P_exact = exact_fields[["MCSE_P"]],
    accounting_exact = exact_fields[["accounting_fields"]],
    oracle_warning = oracle$warning_text,
    candidate_warning = candidate$warning_text,
    oracle_error = oracle$error_message,
    candidate_error = candidate$error_message,
    workspace_reuses = workspace_reuses, debug_checks = debug_checks,
    canary_checks = canary_checks,
    stale_state_failures = stale_failures,
    canary_failures = canary_failures,
    workspace_reuse_positive = counter_checks[["workspace_reuse_positive"]],
    debug_checks_positive = counter_checks[["debug_checks_positive"]],
    canary_checks_positive = counter_checks[["canary_checks_positive"]],
    stale_state_failures_zero = counter_checks[["stale_state_failures_zero"]],
    canary_failures_zero = counter_checks[["canary_failures_zero"]],
    pass = pass, detail = detail, stringsAsFactors = FALSE
  )
}

for (shape in shapes) for (n in n_grid) for (traits in trait_grid) {
  for (nsim in nsim_grid) for (chunk in chunk_grid) {
    # Adjacent cells cycle through the three numerical regimes.
    for (kind_index in seq_along(trait_kinds)) {
      trait_kind <- trait_kinds[[kind_index]]
      cell_index <- cell_index + 1L
      message("[stage2d2] workspace cell ", cell_index,
              ": shape=", shape, ", n=", n, ", traits=", traits,
              ", kind=", trait_kind, ", nsim=", nsim,
              ", chunk=", chunk)
      evidence_rows[[length(evidence_rows) + 1L]] <- run_cell(
        shape, n, trait_kind, traits, nsim, chunk, kind_index
      )
    }
  }
}

evidence <- stage2d2_rbind(evidence_rows)
stage2d2_write_csv(evidence, evidence_path)

batch8_cells <- if (nrow(evidence)) sum(evidence$traits == 8L) else 0L
all_or <- function(name) nrow(evidence) > 0L && all(evidence[[name]])
counter_total <- function(name) {
  if (!nrow(evidence) || !(name %in% names(evidence))) return(NA_real_)
  sum(evidence[[name]], na.rm = TRUE)
}
overall <- nrow(evidence) > 0L && all(evidence$pass) &&
  batch8_cells > 0L && all(evidence$candidate_debug_workspace)
status_reason <- if (overall) {
  "all controlled workspace-safety cells passed"
} else {
  paste0("cells=", nrow(evidence), "; failures=", sum(!evidence$pass),
         "; batch8_cells=", batch8_cells)
}
write_status(
  if (overall) "PASS" else "FAIL", status_reason,
  checks = nrow(evidence), failures = sum(!evidence$pass),
  extra = list(
    source_cpp_opt_in = TRUE, source_cpp = proto_cpp,
    same_translation_unit = prototype_bound, candidate_mode = "A",
    candidate_debug_workspace = TRUE,
    shapes = paste(shapes, collapse = ","),
    trait_kinds = paste(trait_kinds, collapse = ","),
    batch8_cells = as.integer(batch8_cells),
    ordered_sim_exact = all_or("ordered_sim_exact"),
    K_exact = all_or("K_exact"), P_exact = all_or("P_exact"),
    MCSE_P_exact = all_or("MCSE_P_exact"),
    accounting_exact = all_or("accounting_exact"),
    workspace_reuse_positive = all_or("workspace_reuse_positive"),
    debug_checks_positive = all_or("debug_checks_positive"),
    canary_checks_positive = all_or("canary_checks_positive"),
    stale_state_failures = counter_total("stale_state_failures"),
    canary_failures = counter_total("canary_failures")
  )
)

message("[stage2d2] workspace-safety audit complete; production code unchanged")
if (!overall) {
  stop("Stage 2D2 workspace-safety audit failed: ", status_reason,
       call. = FALSE)
}
