#!/usr/bin/env Rscript

# Rebuild summaries from the completed Stage 2C formal timing files. This
# script is deliberately read-only with respect to package production code and
# never invokes a benchmark workload.

run_args <- commandArgs(trailingOnly = TRUE)
result_arg <- grep("^--result-root=", run_args, value = TRUE)
if (length(result_arg)) {
  result_root <- sub("^--result-root=", "", result_arg[[1L]])
} else {
  # The repository-root invocation avoids translating a non-ASCII Windows
  # script path during startup under the benchmark's C locale.
  result_root <- file.path("benchmarks", "stage2c", "results",
                           "2026-09-09-post-v2")
}

if (!dir.exists(result_root)) stop("Formal result directory is missing: ", result_root)

make_summary <- function(timings) {
  key_names <- c("method", "shape", "n", "workload", "path", "nsim", "ncores")
  key_fields <- timings[key_names]
  key_fields[] <- lapply(key_fields, function(value) {
    value <- as.character(value)
    value[is.na(value)] <- "<NA>"
    value
  })
  key <- do.call(paste, c(key_fields, sep = "\r"))
  rows <- lapply(split(seq_len(nrow(timings)), key), function(ii) {
    z <- timings[ii, , drop = FALSE]
    finite <- z$elapsed_ms[is.finite(z$elapsed_ms)]
    med <- if (length(finite)) stats::median(finite) else NA_real_
    data.frame(
      method = z$method[[1L]], shape = z$shape[[1L]], n = z$n[[1L]],
      workload = z$workload[[1L]], path = z$path[[1L]], nsim = z$nsim[[1L]],
      ncores = z$ncores[[1L]], repeats_requested = z$repeats[[1L]],
      observations = nrow(z), ok_observations = sum(z$status == "ok"),
      resource_limit_observations = sum(z$status == "RESOURCE_LIMIT"),
      error_observations = sum(z$status == "error"),
      median_ms = med, IQR_ms = if (length(finite)) stats::IQR(finite) else NA_real_,
      traits_per_sec = if (is.finite(med) && med > 0) 1000 / med else NA_real_,
      per_trait_ms = med, warning_calls = sum(z$warning_count),
      statuses = paste(sort(unique(z$status)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$method, out$shape, out$n, out$workload, out$ncores, out$path), ]
}

method_summaries <- list()
for (name in c("lambda", "d", "delta")) {
  result_dir <- file.path(result_root, paste0(name, "-formal"))
  timings <- utils::read.csv(file.path(result_dir, "stage2c_methods_timings.csv"),
                             stringsAsFactors = FALSE)
  rebuilt <- make_summary(timings)
  utils::write.csv(rebuilt,
                   file.path(result_dir, "stage2c_methods_summary_rebuilt.csv"),
                   row.names = FALSE, na = "")
  method_summaries[[name]] <- rebuilt
}

first_value <- function(x, default = NA_real_) {
  if (length(x)) x[[1L]] else default
}

route_budget <- function(summary) {
  analysis <- summary[summary$path %in% c("raw", "prepared"), , drop = FALSE]
  keys <- unique(analysis[c("method", "shape", "n", "workload", "nsim", "ncores")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, , drop = FALSE]
    same <- analysis$method == k$method & analysis$shape == k$shape &
      analysis$n == k$n & analysis$workload == k$workload &
      analysis$ncores == k$ncores &
      (is.na(analysis$nsim) & is.na(k$nsim) | analysis$nsim == k$nsim)
    z <- analysis[same, , drop = FALSE]
    prep <- summary[summary$method == k$method & summary$shape == k$shape &
                      summary$n == k$n & summary$path == "prepare_tree", , drop = FALSE]
    ref <- summary[summary$method == k$method & summary$shape == k$shape &
                     summary$n == k$n & summary$path == "reference", , drop = FALSE]
    raw_ms <- first_value(z$median_ms[z$path == "raw"])
    prepared_ms <- first_value(z$median_ms[z$path == "prepared"])
    prepare_ms <- first_value(prep$median_ms)
    reference_ms <- first_value(ref$median_ms)
    reference_status <- first_value(ref$statuses, "not_run")
    data.frame(
      method = k$method, shape = k$shape, n = k$n, workload = k$workload,
      nsim = k$nsim, ncores = k$ncores, raw_total_ms = raw_ms,
      prepare_ms = prepare_ms, prepared_compute_ms = prepared_ms,
      prepared_end_to_end_ms = prepare_ms + prepared_ms,
      reference_ms = reference_ms,
      reference_over_raw_speedup = if (is.finite(reference_ms) && is.finite(raw_ms) && raw_ms > 0)
        reference_ms / raw_ms else NA_real_,
      reference_status = reference_status,
      raw_status = first_value(z$statuses[z$path == "raw"], "missing"),
      prepared_status = first_value(z$statuses[z$path == "prepared"], "missing"),
      warning_calls = sum(z$warning_calls),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

routes <- do.call(rbind, lapply(method_summaries, route_budget))
routes <- routes[order(routes$method, routes$shape, routes$n, routes$workload,
                       routes$ncores), ]
utils::write.csv(routes, file.path(result_root, "stage2c_method_route_budget.csv"),
                 row.names = FALSE, na = "")

lookup_median <- function(summary, method, n, path, nsim, ncores = 1L) {
  hit <- summary$method == method & summary$n == n & summary$path == path &
    summary$ncores == ncores & summary$nsim == nsim
  hit[is.na(hit)] <- FALSE
  first_value(summary$median_ms[hit])
}

k_phase <- utils::read.csv(file.path(result_root, "k-formal",
                                     "stage2c_k_permutation_phases.csv"),
                           stringsAsFactors = FALSE)
k_hit <- k_phase[k_phase$n == 10000 & k_phase$route == "prepared" &
                   k_phase$nsim == 9999 & k_phase$mode == "internal_rng" &
                   k_phase$phase == "fast_k_tree_permutation_cpp", , drop = FALSE]

d_sum <- method_summaries$d
d_199 <- lookup_median(d_sum, "D", 2000, "prepared", 199, 1L)
d_9999 <- lookup_median(d_sum, "D", 2000, "prepared", 9999, 1L)
delta_sum <- method_summaries$delta
delta_19 <- lookup_median(delta_sum, "Delta", 500, "prepared", 19, 1L)
delta_199 <- lookup_median(delta_sum, "Delta", 500, "prepared", 199, 1L)

hotspots <- data.frame(
  method = c("K", "D", "Delta"),
  representative_workload = c(
    "prepared n=10000, nsim=9999, internal RNG",
    "prepared n=2000, nsim=9999, ncores=1",
    "prepared n=500, nsim=199, ncores=1"
  ),
  hotspot = c("permutation C++ engine", "random/Brownian null simulation",
              "MCMC plus permutation ACE"),
  hotspot_share_pct = c(100 * first_value(k_hit$share_total),
                        100 * (d_9999 - d_199) / d_9999,
                        100 * (delta_199 - delta_19) / delta_199),
  share_interpretation = c(
    "single low-overhead phase probe; inclusive share",
    "lower bound from added 9800 simulations",
    "lower bound from added 180 test simulations"
  ),
  ideal_end_to_end_ceiling_pct = c(100 * first_value(k_hit$share_total),
                                    100 * (d_9999 - d_199) / d_9999,
                                    100 * (delta_199 - delta_19) / delta_199),
  stringsAsFactors = FALSE
)
utils::write.csv(hotspots, file.path(result_root, "stage2c_hotspot_gates.csv"),
                 row.names = FALSE, na = "")

candidate_decisions <- data.frame(
  rank = c(1L, 2L),
  candidate = c("K_PERMUTATION_NULL_ENGINE_DESIGN",
                "D_NULL_ENGINE_ORDER_STATISTIC_DESIGN"),
  gate_scope = c("heavy prepared K: n=10000, nsim=9999",
                 "heavy prepared D: n=2000, nsim=9999, ncores=1"),
  measured_hotspot_share_pct = c(hotspots$hotspot_share_pct[hotspots$method == "K"],
                                 hotspots$hotspot_share_pct[hotspots$method == "D"]),
  authorization = c("GO_HIGH_PRIORITY_DESIGN_ONLY", "GO_DESIGN_ONLY"),
  interpretation = c(
    paste("This heavy-workload gate supersedes the broad all-grid median as",
          "the Stage 2C candidate decision; the old all-grid table remains",
          "diagnostic evidence for light workloads."),
    "The share is a conservative lower bound from paired simulation-count differences."
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(candidate_decisions,
                 file.path(result_root, "stage2c_final_candidate_decisions.csv"),
                 row.names = FALSE, na = "")

cat("STAGE2C_SUMMARY_REBUILD=PASS\n")
cat("METHOD_ROUTE_ROWS=", nrow(routes), "\n", sep = "")
cat("HOTSPOT_ROWS=", nrow(hotspots), "\n", sep = "")
cat("FINAL_CANDIDATES=", nrow(candidate_decisions), "\n", sep = "")
