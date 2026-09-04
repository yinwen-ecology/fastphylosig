args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: run_candidate1_paired.R BEFORE_LIB AFTER_LIB OUTPUT_DIR")
}
before_lib <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
after_lib <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[3L]], winslash = "/", mustWork = TRUE)

time_raw_k <- function(lib, tree, trait) {
  callr::r(function(tree, trait) {
    library(fastphylosig)
    gc(FALSE)
    elapsed <- system.time({
      result <- fastphylosig::fast_k(
        tree,
        trait,
        test = FALSE,
        verbose = FALSE,
        progress = FALSE
      )
    })[["elapsed"]]
    list(elapsed = unname(elapsed), statistic = unname(as.numeric(result)[[1L]]))
  }, list(tree = tree, trait = trait),
  libpath = c(lib, .libPaths()), spinner = FALSE)
}

time_raw_lambda <- function(lib, tree, trait) {
  callr::r(function(tree, trait) {
    library(fastphylosig)
    gc(FALSE)
    tryCatch({
      elapsed <- system.time({
        result <- fastphylosig::fast_lambda(
          tree,
          trait,
          test = FALSE,
          lambda_profile = FALSE,
          verbose = FALSE,
          progress = FALSE
        )
      })[["elapsed"]]
      list(ok = TRUE, elapsed = unname(elapsed), statistic = result$lambda,
           message = NA_character_)
    }, error = function(error) {
      list(ok = FALSE, elapsed = NA_real_, statistic = NA_real_,
           message = conditionMessage(error))
    })
  }, list(tree = tree, trait = trait),
  libpath = c(lib, .libPaths()), spinner = FALSE)
}

make_pectinate <- function(n) {
  edge <- matrix(0L, nrow = 2L * (n - 1L), ncol = 2L)
  if (n > 2L) {
    for (i in seq_len(n - 2L)) {
      parent <- n + i
      edge[2L * i - 1L, ] <- c(parent, i)
      edge[2L * i, ] <- c(parent, parent + 1L)
    }
  }
  last <- 2L * n - 1L
  edge[2L * n - 3L, ] <- c(last, n - 1L)
  edge[2L * n - 2L, ] <- c(last, n)
  structure(
    list(edge = edge, edge.length = rep(1, nrow(edge)),
         tip.label = paste0("tip", seq_len(n)), Nnode = n - 1L),
    class = "phylo"
  )
}

set.seed(20260904)
k_tree <- ape::rtree(5000L)
k_trait <- stats::setNames(stats::rnorm(5000L), k_tree$tip.label)
lambda_tree <- make_pectinate(1000L)
lambda_trait <- stats::setNames(seq_len(1000L) / 1000,
                                lambda_tree$tip.label)

# Equal warmup outside formal timing.
invisible(time_raw_k(before_lib, k_tree, k_trait))
invisible(time_raw_k(after_lib, k_tree, k_trait))

rows <- vector("list", 20L)
index <- 0L
for (pair in seq_len(10L)) {
  versions <- if (pair %% 2L) c("before", "after") else c("after", "before")
  for (position in seq_along(versions)) {
    version <- versions[[position]]
    result <- time_raw_k(if (version == "before") before_lib else after_lib,
                         k_tree, k_trait)
    index <- index + 1L
    rows[[index]] <- data.frame(
      workload = "K raw n=5000",
      pair = pair,
      position = position,
      version = version,
      elapsed_ms = 1000 * result$elapsed,
      statistic = result$statistic,
      stringsAsFactors = FALSE
    )
  }
}
timings <- do.call(rbind, rows)
utils::write.csv(
  timings,
  file.path(output_dir, "candidate1_paired_timings.csv"),
  row.names = FALSE
)

summarize <- function(values) {
  c(median_ms = stats::median(values),
    iqr_ms = stats::IQR(values))
}
before_summary <- summarize(timings$elapsed_ms[timings$version == "before"])
after_summary <- summarize(timings$elapsed_ms[timings$version == "after"])
summary <- data.frame(
  workload = "K raw n=5000",
  repeats = 10L,
  before_median_ms = before_summary[["median_ms"]],
  before_iqr_ms = before_summary[["iqr_ms"]],
  after_median_ms = after_summary[["median_ms"]],
  after_iqr_ms = after_summary[["iqr_ms"]],
  speedup = before_summary[["median_ms"]] / after_summary[["median_ms"]],
  max_statistic_abs_diff = max(abs(
    timings$statistic[timings$version == "before"] -
      timings$statistic[timings$version == "after"]
  )),
  stringsAsFactors = FALSE
)
utils::write.csv(
  summary,
  file.path(output_dir, "candidate1_paired_summary.csv"),
  row.names = FALSE
)

lambda_evidence <- rbind(
  data.frame(version = "before", time_raw_lambda(before_lib, lambda_tree,
                                                  lambda_trait)),
  data.frame(version = "after", time_raw_lambda(after_lib, lambda_tree,
                                                 lambda_trait))
)
utils::write.csv(
  lambda_evidence,
  file.path(output_dir, "candidate1_pectinate_lambda.csv"),
  row.names = FALSE
)
print(summary)
print(lambda_evidence)
