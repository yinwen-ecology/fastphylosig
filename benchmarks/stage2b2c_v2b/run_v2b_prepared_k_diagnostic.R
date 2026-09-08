#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("usage: run_v2b_prepared_k_diagnostic.R LIBRARY LABEL OUTPUT_CSV")
}
lib <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
label <- args[[2L]]
output <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(fastphylosig, lib.loc = lib))

elapsed_ms <- function(expr) {
  started <- proc.time()[["elapsed"]]
  force(expr)
  1000 * (proc.time()[["elapsed"]] - started)
}
make_tree <- function(shape, n) {
  if (shape == "pectinate") {
    tree <- ape::stree(n, type = "left")
  } else {
    set.seed(20260908L + n)
    tree <- ape::rtree(n)
  }
  tree$tip.label <- paste0("sp", seq_len(n))
  tree$edge.length <- 0.5 + seq_len(nrow(tree$edge)) / nrow(tree$edge)
  tree
}

ns <- asNamespace("fastphylosig")
validate <- get(".validate_prepared_context", ns)
validated <- get(".validated_context", ns)
subset_tree <- get(".prepared_tree_subset", ns)
kernel <- get("fast_k_tree_batch_cpp", ns)

rows <- list()
at <- 0L
for (shape in c("random", "pectinate")) {
  n <- 20000L
  tree <- make_tree(shape, n)
  trait <- stats::setNames(sin(seq_len(n)) + cos(seq_len(n) / 7), tree$tip.label)
  context <- prepare_tree(tree)
  token <- validated(context, verify = TRUE)
  group <- subset_tree(token, need_matrix = FALSE)
  X <- matrix(trait[group$tree$tip.label], ncol = 1L)
  invisible(kernel(group$compiled_tree, X, trait_chunk = 64L))
  for (repeat_id in seq_len(10L)) {
    phases <- c(
      context_validation = elapsed_ms(validate(context)),
      validated_subset_lookup = elapsed_ms(subset_tree(token, need_matrix = FALSE)),
      k_kernel = elapsed_ms(kernel(group$compiled_tree, X, trait_chunk = 64L))
    )
    for (phase in names(phases)) {
      at <- at + 1L
      rows[[at]] <- data.frame(
        version = label, shape = shape, n = n, repeat_id = repeat_id,
        phase = phase, elapsed_ms = unname(phases[[phase]]),
        stringsAsFactors = FALSE
      )
    }
  }
}
raw <- do.call(rbind, rows)
summary <- do.call(rbind, lapply(
  split(raw, interaction(raw$version, raw$shape, raw$phase, drop = TRUE)),
  function(z) data.frame(
    version = z$version[[1L]], shape = z$shape[[1L]], n = z$n[[1L]],
    phase = z$phase[[1L]], observations = nrow(z),
    median_ms = stats::median(z$elapsed_ms), IQR_ms = stats::IQR(z$elapsed_ms),
    stringsAsFactors = FALSE
  )
))
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(summary, output, row.names = FALSE)
print(summary, row.names = FALSE)
