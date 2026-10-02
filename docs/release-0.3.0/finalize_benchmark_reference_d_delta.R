options(stringsAsFactors = FALSE, warn = 1)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "docs", "release-0.3.0")
raw_path <- file.path(out, "benchmark_reference_d_delta_runs.csv")
pairs_path <- file.path(out, "benchmark_reference_d_delta_pairs.csv")
summary_path <- file.path(out, "benchmark_reference_d_delta_summary.csv")
metadata_path <- file.path(out, "benchmark_reference_d_delta_metadata.csv")
qa_path <- file.path(out, "benchmark_reference_d_delta_qa.md")
runner_path <- file.path(out, "benchmark_reference_d_delta.R")
snapshot_path <- file.path(root, "benchmarks", "stage5",
  "reference_d_delta_postprocess", "pairs_before_order_audit_repair.csv")
smoke_meta_path <- file.path(out,
  "benchmark_reference_d_delta_smoke_20261001T181645Z_metadata.csv")
delta_source_path <- file.path(root, "paper", "reference_sources",
  "delta_borges2019_code.R")

hash_file <- function(path) {
  digest::digest(path, algo = "sha256", serialize = FALSE, file = TRUE)
}
if (!all(file.exists(c(raw_path, pairs_path, runner_path, snapshot_path,
                       smoke_meta_path, delta_source_path)))) {
  stop("Missing raw runs, pair snapshot, runner, smoke metadata, or reference source.")
}
expected <- list(
  delta_sha = "0dc1e987c13f0e972b0e7c3b4829b4d0fdc41ff0de3856bae4803f9ce0115175",
  tar_sha = "909b7ad27485dcc2d96a7bb1fc70b65b54e93245d75761aca406b9c52f06cbff",
  snapshot_sha = "96164c7e9e26ed0b720f7ca402c9b1997ba7a153a14c2e92daa803fb7c94bda9",
  executed_runner_sha = "440cb42d4c18ebc35863cc65e0ba771afb1cc9e918e4f9888333c4e1ab31349f"
)
delta_sha <- hash_file(delta_source_path)
snapshot_sha <- hash_file(snapshot_path)
old_pairs_sha <- hash_file(pairs_path)
current_raw_sha <- hash_file(raw_path)
existing_postprocess <- file.exists(c(summary_path, metadata_path, qa_path))
if (any(existing_postprocess)) {
  if (!all(existing_postprocess)) {
    stop("Partial prior post-processing outputs found; preserve and inspect them.")
  }
  prior_metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)
  prior_raw_sha <- prior_metadata$value[
    match("raw_runs_sha256", prior_metadata$field)
  ]
  if (length(prior_raw_sha) != 1L ||
      !identical(prior_raw_sha, current_raw_sha)) {
    stop("Existing derived outputs do not match the immutable raw-run hash.")
  }
}
current_pairs <- read.csv(pairs_path, stringsAsFactors = FALSE)
repaired_pair_table <- nrow(current_pairs) == 240L &&
  "first_implementation" %in% names(current_pairs)
if (!identical(delta_sha, expected$delta_sha) ||
    !identical(snapshot_sha, expected$snapshot_sha) ||
    (!identical(old_pairs_sha, expected$snapshot_sha) &&
       !repaired_pair_table)) {
  stop("Pinned source or pre-repair pair snapshot SHA mismatch.")
}

raw_sha_before <- current_raw_sha
raw <- read.csv(raw_path, stringsAsFactors = FALSE, check.names = FALSE)
if (nrow(raw) != 480L || length(unique(raw$run_id)) != 1L ||
    !identical(unique(raw$profile), "full")) {
  stop("Raw file is not the expected single full run.")
}
if (any(!raw$success) || any(!is.finite(raw$estimate)) ||
    any(!is.finite(raw$elapsed_seconds) | raw$elapsed_seconds <= 0) ||
    any(!is.na(raw$error) & nzchar(raw$error))) {
  stop("A timed call failed or returned a non-finite result; stop.")
}
if (any(raw$n_species_used != raw$n_tip) ||
    any(raw$matched_species != raw$n_tip) || any(raw$n_removed_na != 0L)) {
  stop("Actual species-count or NA audit failed.")
}
if (any(raw$source_tarball_sha256 != expected$tar_sha) ||
    any(raw$package_version != "0.3.0") ||
    length(unique(raw$dll_file_sha256)) != 1L) {
  stop("Package, tarball, or DLL provenance mismatch.")
}
if (any(raw$method == "D" & (!raw$test | raw$d_nsim != 199L)) ||
    any(raw$method == "Delta" & (raw$test | raw$delta_mcmc_sim != 10000L |
      raw$delta_thin != 10L | raw$delta_burn != 100L |
      raw$delta_lambda0 != 0.1 | raw$delta_proposal_sd != 0.5 |
      raw$delta_chains != 2L | raw$delta_entropy != "LSE" |
      raw$delta_model != "ARD"))) {
  stop("Recorded parameter contract differs from the approved design.")
}

fast_d <- raw[raw$method == "D" &
  raw$implementation == "fastphylosig", , drop = FALSE]
fast_delta <- raw[raw$method == "Delta" &
  raw$implementation == "fastphylosig", , drop = FALSE]
if (nrow(fast_d) != 120L || any(fast_d$reported_status != "ok") ||
    any(fast_d$nsim_requested != 199L) ||
    any(fast_d$nsim_successful_random != 199L) ||
    any(fast_d$nsim_successful_brownian != 199L)) {
  stop("Fast D status or actual permutation count audit failed.")
}
delta_controls <- c(
  reported_delta_control_mcmc_sim = "10000",
  reported_delta_control_thin = "10",
  reported_delta_control_burn = "100",
  reported_delta_control_lambda0 = "0.1",
  reported_delta_control_proposal_sd = "0.5",
  reported_delta_control_entropy = "LSE",
  reported_delta_control_model = "ARD",
  reported_delta_control_ace_engine = "fast",
  reported_delta_control_ncores = "1"
)
if (nrow(fast_delta) != 120L || any(fast_delta$reported_status != "ok") ||
    any(fast_delta$delta_expected_total_iterations != 20000) ||
    any(fast_delta$delta_reported_requested_iterations != 20000) ||
    any(fast_delta$delta_reported_successful_iterations != 20000) ||
    any(!fast_delta$delta_requested_iterations_match) ||
    any(!fast_delta$delta_successful_iterations_match)) {
  stop("Fast Delta status or actual iteration counts failed.")
}
for (field in names(delta_controls)) {
  if (any(as.character(fast_delta[[field]]) != delta_controls[[field]])) {
    stop("Fast Delta control mismatch: ", field)
  }
}

ids <- unique(raw$comparison_id)
if (length(ids) != 240L) stop("Expected 240 matched comparison IDs.")
pair_rows <- lapply(ids, function(id) {
  z <- raw[raw$comparison_id == id, , drop = FALSE]
  f <- z[z$implementation == "fastphylosig", , drop = FALSE]
  r <- z[z$implementation == "reference", , drop = FALSE]
  if (nrow(z) != 2L || nrow(f) != 1L || nrow(r) != 1L ||
      f$first_implementation != r$first_implementation) {
    stop("Invalid matched raw rows: ", id)
  }
  data.frame(
    run_id = f$run_id, profile = f$profile, comparison_id = id,
    first_implementation = f$first_implementation,
    method = f$method, reference = f$reference, scenario = f$scenario,
    signal_mixture_h2 = f$signal_mixture_h2, n_tip = f$n_tip,
    replicate = f$replicate,
    tree_seed = f$tree_seed, component_seed = f$component_seed,
    bm_seed = f$bm_seed, noise_seed = f$noise_seed,
    tree_hash = f$tree_hash, bm_hash = f$bm_hash,
    noise_hash = f$noise_hash, trait_hash = f$trait_hash,
    reference_input_hash = f$reference_input_hash,
    d_nsim = f$d_nsim, delta_mcmc_sim = f$delta_mcmc_sim,
    delta_thin = f$delta_thin, delta_burn = f$delta_burn,
    delta_lambda0 = f$delta_lambda0, delta_proposal_sd = f$delta_proposal_sd,
    delta_chains = f$delta_chains, delta_entropy = f$delta_entropy,
    delta_model = f$delta_model,
    delta_expected_total_iterations = f$delta_expected_total_iterations,
    delta_reported_requested_iterations = f$delta_reported_requested_iterations,
    delta_reported_successful_iterations = f$delta_reported_successful_iterations,
    delta_requested_iterations_match = f$delta_requested_iterations_match,
    delta_successful_iterations_match = f$delta_successful_iterations_match,
    fast_reported_delta_control_mcmc_sim = f$reported_delta_control_mcmc_sim,
    fast_reported_delta_control_thin = f$reported_delta_control_thin,
    fast_reported_delta_control_burn = f$reported_delta_control_burn,
    fast_reported_delta_control_lambda0 = f$reported_delta_control_lambda0,
    fast_reported_delta_control_proposal_sd = f$reported_delta_control_proposal_sd,
    fast_reported_delta_control_entropy = f$reported_delta_control_entropy,
    fast_reported_delta_control_model = f$reported_delta_control_model,
    fast_reported_delta_control_ace_engine = f$reported_delta_control_ace_engine,
    fast_reported_delta_control_ncores = f$reported_delta_control_ncores,
    package_version = f$package_version,
    source_tarball_sha256 = f$source_tarball_sha256,
    dll_file_sha256 = f$dll_file_sha256,
    reference_version = f$reference_version,
    test = f$test, timing_clock = f$timing_clock,
    timing_iterations = f$timing_iterations,
    n_fast_species = f$n_species_used, n_reference_species = r$n_species_used,
    fast_matched_species = f$matched_species,
    reference_matched_species = r$matched_species,
    fast_n_removed_na = f$n_removed_na, reference_n_removed_na = r$n_removed_na,
    fast_elapsed_seconds = f$elapsed_seconds,
    reference_elapsed_seconds = r$elapsed_seconds,
    fast_estimate = f$estimate, reference_estimate = r$estimate,
    fast_status = f$reported_status, reference_status = r$reported_status,
    fast_success = f$success, reference_success = r$success,
    complete_pair = f$success && r$success,
    same_tree_hash = identical(f$tree_hash, r$tree_hash),
    same_trait_hash = identical(f$trait_hash, r$trait_hash),
    same_reference_input_hash = identical(f$reference_input_hash,
                                            r$reference_input_hash),
    species_match = f$n_species_used == r$n_species_used &&
      f$n_species_used == f$n_tip && r$n_species_used == r$n_tip &&
      f$matched_species == f$n_tip && r$matched_species == r$n_tip &&
      f$n_removed_na == 0L && r$n_removed_na == 0L,
    fast_warning = f$warning, reference_warning = r$warning,
    fast_error = f$error, reference_error = r$error,
    reference_over_fast_ratio = r$elapsed_seconds / f$elapsed_seconds,
    stringsAsFactors = FALSE
  )
})
pairs <- do.call(rbind, pair_rows)
pairs <- pairs[order(pairs$method, pairs$scenario, pairs$n_tip,
                     pairs$replicate), , drop = FALSE]
rownames(pairs) <- NULL
if (nrow(pairs) != 240L || any(!pairs$complete_pair) ||
    any(!pairs$species_match) || any(!pairs$same_tree_hash) ||
    any(!pairs$same_trait_hash) || any(!pairs$same_reference_input_hash)) {
  stop("Reconstructed pair completeness/species/hash audit failed.")
}
cell_counts <- table(interaction(pairs$method, pairs$scenario, pairs$n_tip,
                                 drop = TRUE))
if (length(cell_counts) != 24L || any(cell_counts != 10L)) {
  stop("Expected 10 pairs in each of 24 cells.")
}
order_counts <- aggregate(
  comparison_id ~ method + scenario + n_tip + first_implementation,
  data = pairs, FUN = length
)
if (any(order_counts$comparison_id != 5L)) {
  stop("Call order is not balanced 5/5 within every cell.")
}

summary_rows <- lapply(split(
  pairs, interaction(pairs$method, pairs$scenario, pairs$n_tip, drop = TRUE)
), function(z) {
  data.frame(
    run_id = z$run_id[[1L]], method = z$method[[1L]],
    scenario = z$scenario[[1L]], signal_mixture_h2 = z$signal_mixture_h2[[1L]],
    n_tip = z$n_tip[[1L]], n_pairs = nrow(z),
    fast_median_seconds = median(z$fast_elapsed_seconds),
    fast_q1_seconds = unname(quantile(z$fast_elapsed_seconds, 0.25)),
    fast_q3_seconds = unname(quantile(z$fast_elapsed_seconds, 0.75)),
    reference_median_seconds = median(z$reference_elapsed_seconds),
    reference_q1_seconds = unname(quantile(z$reference_elapsed_seconds, 0.25)),
    reference_q3_seconds = unname(quantile(z$reference_elapsed_seconds, 0.75)),
    reference_over_fast_ratio_of_medians =
      median(z$reference_elapsed_seconds) / median(z$fast_elapsed_seconds),
    paired_ratio_median = median(z$reference_over_fast_ratio),
    paired_ratio_q1 = unname(quantile(z$reference_over_fast_ratio, 0.25)),
    paired_ratio_q3 = unname(quantile(z$reference_over_fast_ratio, 0.75)),
    warning_calls_fast = sum(!is.na(z$fast_warning) & nzchar(z$fast_warning)),
    warning_calls_reference =
      sum(!is.na(z$reference_warning) & nzchar(z$reference_warning)),
    test = z$test[[1L]], timing_clock = "Sys.time",
    timing_iterations = 1L, preparation_included = FALSE,
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
summary <- summary[order(summary$method, summary$scenario, summary$n_tip), ]
rownames(summary) <- NULL
if (nrow(summary) != 24L || any(summary$n_pairs != 10L)) {
  stop("Unexpected cell summary.")
}

warnings <- raw[!is.na(raw$warning) & nzchar(raw$warning), , drop = FALSE]
warning_components <- do.call(rbind, lapply(seq_len(nrow(warnings)), function(i) {
  parts <- strsplit(warnings$warning[[i]], " || ", fixed = TRUE)[[1L]]
  data.frame(
    method = warnings$method[[i]], implementation = warnings$implementation[[i]],
    scenario = warnings$scenario[[i]], n_tip = warnings$n_tip[[i]],
    warning_component = parts, stringsAsFactors = FALSE
  )
}))
warning_calls <- aggregate(
  rep(1L, nrow(warnings)),
  list(method = warnings$method, implementation = warnings$implementation), sum
)
names(warning_calls)[[3L]] <- "warning_calls"
warning_counts <- sort(table(warning_components$warning_component),
                       decreasing = TRUE)
reference_warning_components <- warning_components[
  warning_components$method == "Delta" &
    warning_components$implementation == "reference", , drop = FALSE
]
reference_warning_counts <- sort(
  table(reference_warning_components$warning_component), decreasing = TRUE
)
reference_warning_by_cell <- xtabs(
  ~ warning_component + n_tip + scenario, reference_warning_components
)
reference_delta <- raw[raw$method == "Delta" &
  raw$implementation == "reference", , drop = FALSE]
diag_ok <- tolower(as.character(fast_delta$diagnostics_available)) == "true"
fast_delta_warning <- fast_delta[!is.na(fast_delta$warning) &
  nzchar(fast_delta$warning), , drop = FALSE]
if (nrow(reference_delta) != 120L ||
    sum(is.finite(reference_delta$estimate)) != 120L ||
    any(!is.na(reference_delta$error) & nzchar(reference_delta$error)) ||
    sum(diag_ok, na.rm = TRUE) != 120L || nrow(fast_delta_warning) != 1L) {
  stop("Delta return-value or diagnostic audit failed.")
}

smoke_meta <- read.csv(smoke_meta_path, stringsAsFactors = FALSE)
executed_runner_sha <- smoke_meta$value[
  match("benchmark_script_sha256", smoke_meta$field)
]
if (!identical(executed_runner_sha, expected$executed_runner_sha)) {
  stop("Cannot verify executed runner SHA from smoke metadata.")
}
corrected_runner_sha <- hash_file(runner_path)
postprocessor_sha <- hash_file(file.path(out, "finalize_benchmark_reference_d_delta.R"))
raw_sha_after <- hash_file(raw_path)
if (!identical(raw_sha_before, raw_sha_after)) stop("Raw CSV changed during audit.")

metadata_list <- c(
  run_id = unique(raw$run_id),
  profile = "full",
  original_benchmark_process_exit_code = "1",
  original_process_outcome = paste(
    "All timed calls and raw rows completed; exit 1 in post-run pair-order",
    "aggregate because first_implementation was omitted from pair projection.",
    "No timed calls were rerun; final pairs and summaries were reconstructed."
  ),
  executed_runner_sha256 = executed_runner_sha,
  corrected_runner_sha256 = corrected_runner_sha,
  postprocessor_sha256 = postprocessor_sha,
  raw_runs_sha256 = raw_sha_before,
  pre_repair_pairs_snapshot_path =
    "benchmarks/stage5/reference_d_delta_postprocess/pairs_before_order_audit_repair.csv",
  pre_repair_pairs_snapshot_sha256 = snapshot_sha,
  source_tarball_sha256 = unique(raw$source_tarball_sha256),
  package_version = unique(raw$package_version),
  dll_file_sha256 = unique(raw$dll_file_sha256),
  D_call = "fast_d(test=TRUE, nsim=199) vs caper::phylo.d(permut=199)",
  D_fast_observed_counts = "199/199 random and 199/199 Brownian in all 120 calls",
  Delta_fast_call =
    "test=FALSE,mcmc_sim=10000,thin=10,burn=100,lambda0=0.1,proposal_sd=0.5,entropy=LSE,model=ARD",
  Delta_reference_source_sha256 = delta_sha,
  Delta_reference_contract =
    "ape::ace(model=ARD), nentropy() LSE-equivalent transform, two emcmc chains",
  Delta_fast_observed_iterations = "requested=20000; successful=20000 in all 120 calls",
  species_grid = "50,100,500,2000",
  scenarios = "none=h2 0; moderate=h2 0.5; strong=h2 0.95",
  paired_replicates_per_cell = "10",
  pair_count = as.character(nrow(pairs)),
  timing = paste(
    "one public function call, Sys.time; tree/trait generation, prepare_tree,",
    "caper comparative.data, matching, seeding and gc excluded; Delta ACE/MCMC included"
  ),
  call_order = "5 fast-first and 5 reference-first per cell",
  warning_calls = as.character(nrow(warnings)),
  delta_reference_warning_calls = as.character(sum(
    warnings$method == "Delta" & warnings$implementation == "reference"
  )),
  delta_fast_warning_calls = as.character(sum(
    warnings$method == "Delta" & warnings$implementation == "fastphylosig"
  )),
  inference_scope =
    "runtime only; not equal-accuracy or final-fit convergence validation"
)
metadata <- data.frame(field = names(metadata_list),
                       value = as.character(metadata_list),
                       stringsAsFactors = FALSE)

show <- function(x) paste(capture.output(print(x, row.names = FALSE)),
                          collapse = "\n")
ratio_text <- show(summary[, c(
  "method", "scenario", "n_tip", "n_pairs", "fast_median_seconds",
  "reference_median_seconds", "reference_over_fast_ratio_of_medians",
  "paired_ratio_median", "warning_calls_fast", "warning_calls_reference"
)])
status_text <- paste(capture.output(print(with(raw, table(
  method, implementation, reported_status, success
)))), collapse = "\n")
order_text <- show(order_counts[order(order_counts$method,
  order_counts$scenario, order_counts$n_tip,
  order_counts$first_implementation), ])
warn_call_text <- show(warning_calls)
warn_component_text <- paste(capture.output(print(warning_counts)),
                             collapse = "\n")
warn_reference_component_text <- paste(
  capture.output(print(reference_warning_counts)), collapse = "\n"
)
warn_cell_text <- paste(capture.output(print(reference_warning_by_cell)),
                        collapse = "\n")
diag_text <- show(data.frame(
  calls = nrow(fast_delta), diagnostics_available = sum(diag_ok, na.rm = TRUE),
  ESS_alpha_min = min(fast_delta$ESS_alpha, na.rm = TRUE),
  ESS_beta_min = min(fast_delta$ESS_beta, na.rm = TRUE),
  Rhat_alpha_max = max(fast_delta$Rhat_alpha, na.rm = TRUE),
  Rhat_beta_max = max(fast_delta$Rhat_beta, na.rm = TRUE),
  MCSE_Delta_max = max(fast_delta$MCSE_Delta, na.rm = TRUE)
))
fast_warn_text <- show(fast_delta_warning[, c(
  "scenario", "n_tip", "replicate", "warning", "estimate",
  "ESS_alpha", "ESS_beta", "Rhat_alpha", "Rhat_beta", "MCSE_Delta"
)])
d_ratios <- summary$reference_over_fast_ratio_of_medians[summary$method == "D"]
delta_ratios <- summary$reference_over_fast_ratio_of_medians[
  summary$method == "Delta"
]
qa <- c(
  "# D/Delta full-run raw audit and post-processing recovery", "",
  "## Run integrity",
  "The original benchmark runner exited 1 after timing completed: its final pair-order aggregate referenced first_implementation, which was absent from the pair projection. Raw timed rows were saved and left unchanged. The previous pairs CSV was snapshotted before repair. No benchmark calls were repeated.",
  paste("Executed runner SHA-256, from the matching immediately preceding smoke run:", executed_runner_sha),
  paste("Corrected runner SHA-256, adding the omitted pair field:", corrected_runner_sha),
  paste("Postprocessor SHA-256:", postprocessor_sha),
  paste("Raw runs SHA-256, verified unchanged:", raw_sha_before),
  paste("Pre-repair pair snapshot SHA-256:", snapshot_sha),
  "The metadata records the original process exit code; it does not claim the benchmark runner exited 0.",
  "", "## Design and provenance",
  paste("Package fastphylosig", unique(raw$package_version),
        "; tarball SHA-256", unique(raw$source_tarball_sha256),
        "; DLL SHA-256", unique(raw$dll_file_sha256)),
  "D: fast_d(test=TRUE, nsim=199) versus caper::phylo.d(permut=199); all fast calls report 199/199 random and 199/199 Brownian simulations.",
  paste("Delta reference source SHA-256:", delta_sha,
        "; source uses ape::ace(model=ARD), nentropy() and two emcmc chains."),
  "Fast Delta actual controls and requested/successful total iterations were audited from returned fields: 10000, thin 10, burn 100, lambda0 0.1, proposal_sd 0.5, LSE, ARD, two-chain contract, total 20000.",
  "Each timed interval measured one public call with Sys.time. Tree/trait generation, prepare_tree, comparative.data, matching, seeding and gc were excluded; Delta ACE/reconstruction and MCMC remained inside the call.",
  "The four panels compare runtime for different configured analyses: D uses a permutation test, Delta uses two-chain MCMC, and K/lambda use test=FALSE point estimates. They do not imply equal statistical work.",
  "", "## Raw audit",
  paste("Raw calls:", nrow(raw), "(240 fast, 240 reference); pairs:",
        nrow(pairs), "; cells: 24; paired replicates/cell: 10."),
  paste("Finite successful calls:", sum(raw$success), "/", nrow(raw),
        "; errors:", sum(!is.na(raw$error) & nzchar(raw$error))),
  paste("Observed actual species match:", sum(pairs$species_match), "/",
        nrow(pairs), "; all counts equal n_tip; NA removals: 0."),
  paste("Tree, trait and reference-input hash matches:",
        sum(pairs$same_tree_hash), "/", nrow(pairs), ";",
        sum(pairs$same_trait_hash), "/", nrow(pairs), ";",
        sum(pairs$same_reference_input_hash), "/", nrow(pairs)),
  "First-call order balance by cell:",
  "~~~text", order_text, "~~~",
  "Statuses (reference APIs do not expose a status field):",
  "~~~text", status_text, "~~~",
  "", "## Warning and Delta diagnostic audit",
  paste("Warning-bearing calls:", nrow(warnings), "/480; D=0; Delta reference=112/120; fast Delta=1/120. All warnings remain in raw rows; no retries or row removals."),
  "Warning-call counts by method and implementation:",
  "~~~text", warn_call_text, "~~~",
  "All warning-component frequencies (components may co-occur):",
  "~~~text", warn_component_text, "~~~",
  "Reference Delta warning-component frequencies:",
  "~~~text", warn_reference_component_text, "~~~",
  "Reference Delta warning components by species count and scenario:",
  "~~~text", warn_cell_text, "~~~",
  "The pinned reference delta() returns only a numeric scalar and discards the ACE convergence object; result_fields are <numeric>. Call-level warning capture cannot distinguish an optimizer trial warning from an accepted/final-fit warning. Finite estimates do not certify convergence.",
  "Fast Delta chain diagnostics:",
  "~~~text", diag_text, "~~~",
  "One fast Delta warning row:",
  "~~~text", fast_warn_text, "~~~",
  "Delta elapsed times describe calls returning finite values only. They do not validate equal accuracy, effective chains, or final-fit convergence.",
  "", "## Observed runtime",
  "Ratio is median(reference seconds) / median(fast seconds), descriptive only:",
  "~~~text", ratio_text, "~~~",
  paste("D ratio range:", format(min(d_ratios), digits=4), "to",
        format(max(d_ratios), digits=4), "x."),
  paste("Delta ratio range:", format(min(delta_ratios), digits=4), "to",
        format(max(delta_ratios), digits=4), "x."),
  "The Delta figure panel must disclose reference warnings and the absence of a convergence certification."
)

write.csv(pairs, pairs_path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(summary, summary_path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(metadata, metadata_path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
writeLines(qa, qa_path, useBytes = TRUE)
if (!identical(hash_file(raw_path), raw_sha_before)) {
  stop("Raw timed-run CSV changed during post-processing.")
}
message("Reconstructed 240 pairs and 24 summaries from immutable raw rows; no benchmark calls were run.")
