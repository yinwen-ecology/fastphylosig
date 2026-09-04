args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: collect_final_call_counts.R PACKAGE_LIB OUTPUT_CSV")
}
.libPaths(c(normalizePath(args[[1L]], winslash = "/", mustWork = TRUE),
            .libPaths()))
library(fastphylosig)

set.seed(20260904)
tree <- ape::rtree(500L)
trait <- stats::setNames(stats::rnorm(500L), tree$tip.label)
context <- prepare_tree(tree)
counts <- new.env(parent = emptyenv())
counts$fingerprint <- 0L
counts$inspection <- 0L
counts$subset <- 0L
original_fingerprint <- fastphylosig:::.tree_fingerprint
original_inspection <- fastphylosig:::.inspect_tree_core
original_subset <- fastphylosig:::.prepared_tree_subset

testthat::local_mocked_bindings(
  .tree_fingerprint = function(...) {
    counts$fingerprint <- counts$fingerprint + 1L
    original_fingerprint(...)
  },
  .inspect_tree_core = function(...) {
    counts$inspection <- counts$inspection + 1L
    original_inspection(...)
  },
  .prepared_tree_subset = function(...) {
    counts$subset <- counts$subset + 1L
    original_subset(...)
  },
  .package = "fastphylosig"
)

invisible(fast_k(context, trait, test = FALSE,
                 verbose = FALSE, progress = FALSE))
result <- data.frame(
  path = "prepared_vector_n500",
  function_name = c(".tree_fingerprint", ".inspect_tree_core",
                    ".prepared_tree_subset"),
  calls = c(counts$fingerprint, counts$inspection, counts$subset),
  stringsAsFactors = FALSE
)
utils::write.csv(result, args[[2L]], row.names = FALSE)
print(result)
