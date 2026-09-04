args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: run_candidate_paired.R BEFORE_LIB AFTER_LIB LABEL OUTPUT_DIR")
}
before_lib <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
after_lib <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
label <- args[[3L]]
output_dir <- normalizePath(args[[4L]], winslash = "/", mustWork = TRUE)

time_case <- function(lib, workload, tree, trait, warm = FALSE,
                      iterations = 1L) {
  callr::r(function(workload, tree, trait, warm, iterations) {
    library(fastphylosig)
    prepared <- grepl("prepared", workload, fixed = TRUE)
    context <- if (prepared) fastphylosig::prepare_tree(tree) else tree
    run <- function() {
      if (grepl("lambda", workload, fixed = TRUE)) {
        fastphylosig::fast_lambda(
          context, trait, test = FALSE, lambda_profile = FALSE,
          verbose = FALSE, progress = FALSE
        )
      } else {
        fastphylosig::fast_k(
          context, trait, test = FALSE, verbose = FALSE, progress = FALSE
        )
      }
    }
    if (isTRUE(warm)) invisible(run())
    gc(FALSE)
    elapsed <- system.time({
      for (iteration in seq_len(iterations)) result <- run()
    })[["elapsed"]] / iterations
    statistic <- if (is.data.frame(result) && "K_fast" %in% names(result)) {
      unname(result$K_fast[[1L]])
    } else if (is.list(result) && !is.null(result$lambda)) {
      unname(result$lambda[[1L]])
    } else if (is.list(result) && !is.null(result$K)) {
      unname(result$K[[1L]])
    } else {
      unname(as.numeric(result)[[1L]])
    }
    list(elapsed = unname(elapsed), statistic = statistic)
  }, list(workload = workload, tree = tree, trait = trait, warm = warm,
          iterations = iterations),
  libpath = c(lib, .libPaths()), spinner = FALSE)
}

set.seed(20260904)
tree_500 <- ape::rtree(500L)
trait_500 <- stats::setNames(stats::rnorm(500L), tree_500$tip.label)
matrix_500 <- matrix(stats::rnorm(500L * 100L), nrow = 500L, ncol = 100L,
                     dimnames = list(tree_500$tip.label,
                                     paste0("trait", seq_len(100L))))
tree_5000 <- ape::rtree(5000L)
trait_5000 <- stats::setNames(stats::rnorm(5000L), tree_5000$tip.label)

cases <- list(
  list(workload = "K prepared n=500", tree = tree_500, trait = trait_500,
       iterations = 20L),
  list(workload = "K prepared n=5000", tree = tree_5000,
       trait = trait_5000, iterations = 1L),
  list(workload = "K raw n=5000", tree = tree_5000, trait = trait_5000,
       iterations = 1L),
  list(workload = "K prepared matrix n=500 traits=100", tree = tree_500,
       trait = matrix_500, iterations = 10L),
  list(workload = "lambda prepared n=500", tree = tree_500,
       trait = trait_500, iterations = 20L)
)

all_rows <- list()
row_index <- 0L
for (case in cases) {
  # Equal workload-specific warmup outside formal timing.
  invisible(time_case(before_lib, case$workload, case$tree, case$trait,
                      warm = TRUE, iterations = case$iterations))
  invisible(time_case(after_lib, case$workload, case$tree, case$trait,
                      warm = TRUE, iterations = case$iterations))
  for (pair in seq_len(10L)) {
    versions <- if (pair %% 2L) c("before", "after") else c("after", "before")
    for (position in seq_along(versions)) {
      version <- versions[[position]]
      result <- time_case(
        if (version == "before") before_lib else after_lib,
        case$workload,
        case$tree,
        case$trait,
        warm = grepl("prepared", case$workload, fixed = TRUE),
        iterations = case$iterations
      )
      row_index <- row_index + 1L
      all_rows[[row_index]] <- data.frame(
        workload = case$workload,
        pair = pair,
        position = position,
        version = version,
        iterations = case$iterations,
        elapsed_ms = 1000 * result$elapsed,
        statistic = result$statistic,
        stringsAsFactors = FALSE
      )
    }
  }
}
timings <- do.call(rbind, all_rows)
utils::write.csv(
  timings,
  file.path(output_dir, paste0(label, "_paired_timings.csv")),
  row.names = FALSE
)

summaries <- lapply(unique(timings$workload), function(workload) {
  rows <- timings[timings$workload == workload, , drop = FALSE]
  before <- rows[rows$version == "before", , drop = FALSE]
  after <- rows[rows$version == "after", , drop = FALSE]
  before_median <- stats::median(before$elapsed_ms)
  after_median <- stats::median(after$elapsed_ms)
  data.frame(
    workload = workload,
    repeats = nrow(before),
    before_median_ms = before_median,
    before_iqr_ms = stats::IQR(before$elapsed_ms),
    after_median_ms = after_median,
    after_iqr_ms = stats::IQR(after$elapsed_ms),
    speedup = before_median / after_median,
    max_statistic_abs_diff = max(abs(before$statistic - after$statistic)),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summaries)
utils::write.csv(
  summary,
  file.path(output_dir, paste0(label, "_paired_summary.csv")),
  row.names = FALSE
)
print(summary)
