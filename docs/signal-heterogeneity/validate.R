# Reproduce numerical tests and bounded large-tree smoke workloads.
Sys.setlocale("LC_CTYPE", "English_United States.utf8")
task_library <- readLines("docs/signal-heterogeneity/library-path.txt", warn = FALSE)
.libPaths(c(task_library, .libPaths()))
library(fastphylosig)
out <- "docs/signal-heterogeneity"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
testthat::set_max_fails(Inf)
tests <- testthat::test_dir("tests/testthat", reporter = "silent",
  env = new.env(parent = asNamespace("fastphylosig")),
  stop_on_failure = FALSE, load_package = "none")
results <- as.data.frame(tests)
totals <- data.frame(passed = sum(results$passed), failed = sum(results$failed),
                     warnings = sum(results$warning), skipped = sum(results$skipped))
write.csv(totals, file.path(out, "test-totals.csv"), row.names = FALSE)
write.csv(results[, c("file", "test", "nb", "failed", "skipped", "error", "warning", "passed")],
          file.path(out, "test-results.csv"), row.names = FALSE)
print(totals)
if (any(totals[c("failed", "warnings")] > 0)) {
  print(results[results$failed > 0 | results$warning > 0 | results$error,
                c("file", "test", "failed", "error", "warning")])
  stop("Regression tests did not pass.")
}

timings <- list()
for (n in c(2000L, 50000L)) {
  set.seed(100 + n)
  tree <- ape::rtree(n)
  x <- setNames(rnorm(n), tree$tip.label)
  gc()
  timing <- system.time(fit <- fast_hetero(
    tree, x, min_tips = 20, max_tips = 200, progress = FALSE))
  stopifnot(!fit$selection$overlap,
            all(fit$clades$n_tip_total >= 20 & fit$clades$n_tip_total <= 200),
            all(fit$clades$status == "ok"),
            fit$selection$total_clade_tips <= n)
  timings[[length(timings) + 1L]] <- data.frame(n_tip = n, test = FALSE, nsim = 0L,
    n_clades = length(fit$selection$nodes), coverage = fit$selection$coverage,
    elapsed_seconds = unname(timing[["elapsed"]]))
  if (n == 2000L) {
    set.seed(88)
    timing <- system.time(fit <- fast_hetero(
      tree, x, min_tips = 20, max_tips = 200, test = TRUE, nsim = 199,
      simulation_chunk = 32, progress = FALSE))
    stopifnot(fit$heterogeneity$status == "ok", fit$heterogeneity$nsim_successful == 199)
    timings[[length(timings) + 1L]] <- data.frame(n_tip = n, test = TRUE, nsim = 199L,
      n_clades = length(fit$selection$nodes), coverage = fit$selection$coverage,
      elapsed_seconds = unname(timing[["elapsed"]]))
  }
}
timings <- do.call(rbind, timings)
write.csv(timings, file.path(out, "smoke-timings.csv"), row.names = FALSE)
print(timings)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))

# One fully reproducible illustration; these are simulated, not empirical data.
tree <- ape::stree(128, type = "balanced")
tree$edge.length <- rep(1, nrow(tree$edge))
base <- fast_hetero(tree, setNames(sin(1:128), tree$tip.label),
  min_tips = 20, max_tips = 32, progress = FALSE)
set.seed(18)
bm_tree <- ape::extract.clade(tree, base$selection$nodes[[2L]])
bm_core <- fastphylosig:::compile_tree_cpp(bm_tree$edge, bm_tree$edge.length, 32L)
reference_env <- new.env(parent = asNamespace("fastphylosig"))
sys.source("tests/testthat/helper-heterogeneity-reference.R", reference_env)
x <- setNames(c(rep(c(-1, 1), each = 16), as.numeric(reference_env$.hetero_bm_reference(bm_core, 1)),
                rnorm(32), rep(c(-1, 1), 16)), tree$tip.label)
set.seed(19)
demo <- fast_hetero(tree, x, min_tips = 20, max_tips = 32,
  test = TRUE, nsim = 199, keep_null = TRUE, progress = FALSE)
saveRDS(list(tree = tree, trait = x, result = demo), file.path(out, "demo.rds"))
write.csv(demo$clades, file.path(out, "demo-clades.csv"), row.names = FALSE)
write.csv(demo$heterogeneity, file.path(out, "demo-heterogeneity.csv"), row.names = FALSE)
grDevices::png(file.path(out, "demo.png"), width = 1200, height = 700, res = 140)
graphics::par(mar = c(4.5, 5, 5.5, 1))
plot(demo, xlab = "Blomberg's K", main = "")
graphics::mtext("Simulated clade signal heterogeneity", side = 3, line = 3.7,
                font = 2, cex = 1.2)
graphics::mtext(sprintf("H = %.3f; BM reference P = %.3f; 4 clades, 32 tips each",
  demo$heterogeneity$dispersion, demo$heterogeneity$P_BM), side = 3, line = 2.3,
  cex = 0.9)
grDevices::dev.off()
print(demo$heterogeneity)
