#!/usr/bin/env Rscript

# Combine independently serialized tree-shape runs without re-running timings.

args <- commandArgs(trailingOnly = TRUE)
results_dir <- normalizePath(
  if (length(args)) args[[1L]] else file.path("benchmarks", "stage2b2a", "results"),
  winslash = "/", mustWork = TRUE
)
shapes <- c("balanced", "random", "pectinate")

read_shape <- function(file) {
  do.call(rbind, lapply(shapes, function(shape) {
    utils::read.csv(file.path(results_dir, shape, file), check.names = FALSE)
  }))
}

prepare <- read_shape("prepare_tree_summary.csv")
helpers <- read_shape("helper_phase_summary.csv")
operations <- read_shape("operation_counts.csv")
calls <- read_shape("helper_call_counts.csv")
calibration <- read_shape("instrumentation_calibration_summary.csv")

# One early shape run used integer multiplication for this diagnostic upper
# bound. Recompute it in double precision from the preserved exact counts.
operations$estimated_name_comparisons_linear_upper <-
  as.double(operations$named_adjacency_lookup_total) *
  as.double(operations$named_adjacency_parent_groups)

empirical_exponents <- function(data, groups) {
  split_key <- interaction(data[groups], drop = TRUE, lex.order = TRUE)
  rows <- list()
  id <- 0L
  for (key in unique(split_key)) {
    z <- data[split_key == key, , drop = FALSE]
    z <- z[order(z$n), , drop = FALSE]
    if (nrow(z) < 2L) next
    for (i in seq_len(nrow(z) - 1L)) {
      id <- id + 1L
      rows[[id]] <- cbind(
        z[i, groups, drop = FALSE],
        data.frame(
          n_from = z$n[[i]], n_to = z$n[[i + 1L]],
          elapsed_from_ms = z$median_ms[[i]],
          elapsed_to_ms = z$median_ms[[i + 1L]],
          empirical_p = log(z$median_ms[[i + 1L]] / z$median_ms[[i]]) /
            log(z$n[[i + 1L]] / z$n[[i]]),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  do.call(rbind, rows)
}

key_phases <- c(
  "inspect_tree_core", "safe_canonicalize_core", "canonical_tree_signature",
  "preorder_connectivity", "root_distance", "signature_descendant_keys"
)
key_helpers <- helpers[helpers$phase %in% key_phases, , drop = FALSE]

phase_wide <- reshape(
  key_helpers[, c("shape", "n", "phase", "median_ms")],
  idvar = c("shape", "n"), timevar = "phase", direction = "wide"
)
names(phase_wide) <- sub("^median_ms\\.", "", names(phase_wide))
budget <- merge(
  prepare[, c("shape", "n", "median_ms")], phase_wide,
  by = c("shape", "n"), all.x = TRUE
)
names(budget)[names(budget) == "median_ms"] <- "prepare_tree_ms"
budget$inspection_share_pct <- 100 * budget$inspect_tree_core / budget$prepare_tree_ms
budget$canonicalization_share_pct <-
  100 * budget$safe_canonicalize_core / budget$prepare_tree_ms
budget$signature_share_pct <-
  100 * budget$canonical_tree_signature / budget$prepare_tree_ms
budget$inspection_traversal_ceiling_pct <-
  100 * (budget$preorder_connectivity + budget$root_distance) /
  budget$prepare_tree_ms
budget$descendant_key_ceiling_pct <-
  100 * budget$signature_descendant_keys / budget$prepare_tree_ms

utils::write.csv(prepare, file.path(results_dir, "combined_prepare_tree_summary.csv"),
                 row.names = FALSE)
utils::write.csv(helpers, file.path(results_dir, "combined_helper_phase_summary.csv"),
                 row.names = FALSE)
utils::write.csv(empirical_exponents(prepare, "shape"),
                 file.path(results_dir, "combined_prepare_tree_empirical_exponent.csv"),
                 row.names = FALSE)
utils::write.csv(empirical_exponents(key_helpers, c("shape", "phase")),
                 file.path(results_dir, "combined_helper_empirical_exponent.csv"),
                 row.names = FALSE)
utils::write.csv(operations, file.path(results_dir, "combined_operation_counts.csv"),
                 row.names = FALSE)
utils::write.csv(calls, file.path(results_dir, "combined_helper_call_counts.csv"),
                 row.names = FALSE)
utils::write.csv(calibration,
                 file.path(results_dir, "combined_instrumentation_calibration.csv"),
                 row.names = FALSE)
utils::write.csv(budget[order(budget$shape, budget$n), ],
                 file.path(results_dir, "combined_hotspot_budget.csv"),
                 row.names = FALSE)
