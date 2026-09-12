#!/usr/bin/env Rscript

# Derive the compact Stage 2E1A qualification tables from the completed
# formal timing evidence. This script is audit-only and never loads or edits
# package production code.

args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[[1L]] else
  file.path("benchmarks", "stage2e1a", "results", "formal-authoritative")
output <- if (length(args) >= 2L) args[[2L]] else input
dir.create(output, recursive = TRUE, showWarnings = FALSE)

raw <- utils::read.csv(file.path(input, "stage2e1a_phase_timings.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
summary <- utils::read.csv(file.path(input, "stage2e1a_phase_summary.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(raw) != 540L || any(raw$status != "PASS")) {
  stop("Expected 540 passing formal timing rows.", call. = FALSE)
}

heavy <- summary[summary$n >= 2000L & summary$nsim == 9999L &
                   summary$timing_status == "PASS", , drop = FALSE]
if (nrow(heavy) != 18L) {
  stop("Expected 18 passing representative heavy cells.", call. = FALSE)
}

ratio_median <- function(x, numerator, denominator) {
  stats::median(x[[numerator]] / x[[denominator]], na.rm = TRUE)
}

groups <- split(raw, raw$cell_id, drop = TRUE)
paired <- do.call(rbind, lapply(groups, function(x) {
  data.frame(
    cell_id = x$cell_id[[1L]],
    paired_random_fraction = ratio_median(
      x, "random_null_s", "complete_D_ncores1_s"
    ),
    paired_brownian_generation_fraction = ratio_median(
      x, "brownian_generation_s", "complete_D_ncores1_s"
    ),
    paired_sort_threshold_fraction = ratio_median(
      x, "brownian_sort_threshold_s", "complete_D_ncores1_s"
    ),
    paired_binary_D_fraction = ratio_median(
      x, "binary_D_s", "complete_D_ncores1_s"
    ),
    paired_sort_fraction_of_complete_brownian = ratio_median(
      x, "brownian_sort_s", "complete_brownian_s"
    ),
    paired_ncores2_ratio = ratio_median(
      x, "complete_D_ncores2_s", "complete_D_ncores1_s"
    ),
    stringsAsFactors = FALSE
  )
}))

keep <- c(
  "cell_id", "shape", "n", "prevalence", "nsim",
  "complete_D_ncores1_s_median_s", "complete_D_ncores1_s_iqr_s",
  "random_fraction_of_total_D", "Brownian_generation_fraction",
  "sort_fraction_of_complete_Brownian",
  "sort_threshold_fraction_of_total_D", "binary_D_fraction",
  "theoretical_free_sort_speedup", "speedup_old_sort_20pct",
  "speedup_old_sort_30pct", "speedup_old_sort_50pct", "ncores2_ratio"
)
heavy <- merge(heavy[, keep, drop = FALSE], paired, by = "cell_id",
               all.x = TRUE, sort = FALSE)
heavy <- heavy[order(heavy$prevalence, heavy$n,
                     match(heavy$shape, c("balanced", "random", "pectinate"))), ]
rownames(heavy) <- NULL
utils::write.csv(heavy, file.path(output, "stage2e1a_heavy_qualification.csv"),
                 row.names = FALSE, na = "NA")

range_row <- function(metric, values) {
  data.frame(metric = metric, minimum = min(values, na.rm = TRUE),
             median = stats::median(values, na.rm = TRUE),
             maximum = max(values, na.rm = TRUE), stringsAsFactors = FALSE)
}

aggregate <- do.call(rbind, list(
  range_row("heavy_total_s", heavy$complete_D_ncores1_s_median_s),
  range_row("random_fraction", heavy$random_fraction_of_total_D),
  range_row("brownian_generation_fraction", heavy$Brownian_generation_fraction),
  range_row("full_sort_fraction_of_complete_brownian",
            heavy$sort_fraction_of_complete_Brownian),
  range_row("sort_threshold_fraction_of_total_D",
            heavy$sort_threshold_fraction_of_total_D),
  range_row("binary_D_fraction", heavy$binary_D_fraction),
  range_row("paired_sort_threshold_fraction",
            heavy$paired_sort_threshold_fraction),
  range_row("theoretical_free_sort_speedup",
            heavy$theoretical_free_sort_speedup),
  range_row("speedup_with_20pct_old_sort_cost",
            heavy$speedup_old_sort_20pct),
  range_row("speedup_with_30pct_old_sort_cost",
            heavy$speedup_old_sort_30pct),
  range_row("speedup_with_50pct_old_sort_cost",
            heavy$speedup_old_sort_50pct),
  range_row("ncores2_over_ncores1", heavy$ncores2_ratio)
))
utils::write.csv(aggregate,
                 file.path(output, "stage2e1a_heavy_aggregate.csv"),
                 row.names = FALSE, na = "NA")

decision <- data.frame(
  stage = "Stage 2E1A Final D Candidate Qualification",
  source_commit = "03d72813f9da08b27f924aa37437ab2715cf8dfe",
  correctness_checks = 1299L,
  correctness_failures = 0L,
  formal_cells = nrow(summary),
  formal_repeats_per_cell = 10L,
  formal_timing_rows = nrow(raw),
  formal_timing_failures = sum(raw$status != "PASS"),
  representative_heavy_cells = nrow(heavy),
  heavy_sort_threshold_ge_30pct = sum(
    heavy$sort_threshold_fraction_of_total_D >= 0.30
  ),
  heavy_sort_threshold_ge_20pct = sum(
    heavy$sort_threshold_fraction_of_total_D >= 0.20
  ),
  D_EXACT_ORDER_STATISTIC = "NO_GO",
  D_OPTIMIZATION_0_2_0 = "CLOSED",
  NEXT_STAGE = "FINAL_0_2_0_VALIDATION",
  production_code_changed = "NO",
  stringsAsFactors = FALSE
)
utils::write.csv(decision, file.path(output, "FINAL_DECISION.csv"),
                 row.names = FALSE, na = "NA")

print(decision)
