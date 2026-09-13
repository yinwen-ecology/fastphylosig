# Final release clean-install smoke runner.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: run_clean_smoke.R <library> <output>", call. = FALSE)
}
smoke_library <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- args[[2L]]
.libPaths(c(smoke_library, .libPaths()))

connection <- file(output, open = "wt", encoding = "UTF-8")
sink(connection)
sink(connection, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(connection)
}, add = TRUE)

library(fastphylosig)

set.seed(20260913)
tree <- ape::rtree(36)
tips <- tree$tip.label
trait <- setNames(seq_len(36) / 7 + stats::rnorm(36, sd = 0.05), tips)
binary <- setNames(rep(c(0L, 1L), 18), tips)
categorical <- ape::rTraitDisc(
  tree, model = "ER", k = 3, rate = 1.2, states = letters[1:3]
)
categorical <- setNames(as.character(categorical[tips]), tips)

checked <- check_tree(tree)
context <- prepare_tree(tree)
results <- list(
  raw_k = fast_k(tree, trait, test = FALSE, progress = FALSE, verbose = FALSE),
  prepared_k = fast_k(context, trait, test = FALSE, progress = FALSE, verbose = FALSE),
  raw_lambda = fast_lambda(tree, trait, test = FALSE, progress = FALSE, verbose = FALSE),
  prepared_lambda = fast_lambda(context, trait, test = FALSE, progress = FALSE, verbose = FALSE),
  raw_d = fast_d(tree, binary, test = FALSE, progress = FALSE, verbose = FALSE),
  prepared_d = fast_d(context, binary, test = FALSE, progress = FALSE, verbose = FALSE),
  raw_delta = fast_delta(tree, categorical, test = FALSE, model = "ER", progress = FALSE, verbose = FALSE),
  prepared_delta = fast_delta(context, categorical, test = FALSE, model = "ER", progress = FALSE, verbose = FALSE),
  raw_ace = fast_ace(categorical, phy = tree, CI = FALSE, model = "ER", progress = FALSE),
  prepared_ace = fast_ace(categorical, prepared = context, CI = FALSE, model = "ER", progress = FALSE),
  raw_signal = fast_signal(tree, data = trait, method = "K", test = FALSE, progress = FALSE, verbose = FALSE),
  prepared_signal = fast_signal(context, data = trait, method = "K", test = FALSE, progress = FALSE, verbose = FALSE)
)

stopifnot(
  isTRUE(checked$ready),
  inherits(context, "fastphylosig_tree"),
  all(vapply(results, function(result) !is.null(result), logical(1))),
  isTRUE(all.equal(as.numeric(results$raw_k),
                   as.numeric(results$prepared_k), tolerance = 0)),
  isTRUE(all.equal(unname(results$raw_lambda$lambda),
                   unname(results$prepared_lambda$lambda), tolerance = 0))
)

cat("CLEAN_INSTALL_SMOKE=PASS\n")
cat("PACKAGE_VERSION=", as.character(packageVersion("fastphylosig")), "\n", sep = "")
cat("PATHS=raw,prepared\n")
cat("ENTRY_POINTS=check_tree,prepare_tree,fast_k,fast_lambda,fast_d,fast_delta,fast_ace,fast_signal\n")
