#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("usage: run_v2b_prepared_k_confirmation.R BEFORE_LIB AFTER_LIB OUTPUT")
}
libs <- c(candidate1 = normalizePath(args[[1L]], mustWork = TRUE),
          v2 = normalizePath(args[[2L]], mustWork = TRUE))
output <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
sessions <- lapply(libs, function(lib) {
  session <- callr::r_session$new(options = callr::r_session_options(
    user_profile = FALSE, system_profile = FALSE
  ))
  session$run(function(lib) {
    .libPaths(c(lib, .libPaths()))
    suppressPackageStartupMessages(library(fastphylosig, lib.loc = lib))
  }, args = list(lib = lib))
  session
})
on.exit(for (session in sessions) try(session$close(), silent = TRUE), add = TRUE)

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

rows <- list()
at <- 0L
for (shape in c("random", "pectinate")) {
  tree <- make_tree(shape, 20000L)
  trait <- stats::setNames(
    sin(seq_len(20000L)) + cos(seq_len(20000L) / 7), tree$tip.label
  )
  for (version in names(sessions)) {
    sessions[[version]]$run(function(tree, trait) {
      assign(".v2b_ctx", fastphylosig::prepare_tree(tree), .GlobalEnv)
      assign(".v2b_trait", trait, .GlobalEnv)
      invisible(fastphylosig::fast_k(
        .v2b_ctx, .v2b_trait, test = FALSE, verbose = FALSE, progress = FALSE
      ))
    }, args = list(tree = tree, trait = trait))
  }
  for (pair in seq_len(20L)) {
    version_order <- if (pair %% 2L) names(sessions) else rev(names(sessions))
    for (position in seq_along(version_order)) {
      version <- version_order[[position]]
      elapsed <- sessions[[version]]$run(function() {
        started <- proc.time()[["elapsed"]]
        for (i in seq_len(5L)) {
          invisible(fastphylosig::fast_k(
            .v2b_ctx, .v2b_trait, test = FALSE,
            verbose = FALSE, progress = FALSE
          ))
        }
        1000 * (proc.time()[["elapsed"]] - started) / 5
      })
      at <- at + 1L
      rows[[at]] <- data.frame(
        version = version, shape = shape, n = 20000L, pair = pair,
        position = position, inner_iterations = 5L,
        elapsed_ms = as.numeric(elapsed), stringsAsFactors = FALSE
      )
    }
  }
}
raw <- do.call(rbind, rows)
summary <- do.call(rbind, lapply(
  split(raw, interaction(raw$version, raw$shape, drop = TRUE)),
  function(z) data.frame(
    version = z$version[[1L]], shape = z$shape[[1L]], n = z$n[[1L]],
    observations = nrow(z), inner_iterations = z$inner_iterations[[1L]],
    median_ms = stats::median(z$elapsed_ms), IQR_ms = stats::IQR(z$elapsed_ms),
    stringsAsFactors = FALSE
  )
))
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(raw, sub("\\.csv$", "_raw.csv", output), row.names = FALSE)
utils::write.csv(summary, output, row.names = FALSE)
print(summary, row.names = FALSE)
