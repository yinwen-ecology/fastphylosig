# Benchmark single-trait multicore paths.
#
# This script compares ncores = 1, 2, and 4 for the parts of the package where
# one trait still creates many independent numerical jobs:
#   - fast_signal(..., method = "K", test = TRUE)
#   - fast_d(..., nsim = ...)
#
# It intentionally uses precomputed permutation/state matrices so the timing
# focuses on the C++ kernels affected by ncores.

library(ape)
library(fastphylosig)

time_expr <- function(expr) {
  gc()
  unname(system.time(force(expr))[["elapsed"]])
}

make_perms <- function(n, nsim) {
  out <- matrix(NA_integer_, nrow = nsim, ncol = n)
  out[1L, ] <- seq_len(n)
  if (nsim > 1L) {
    for (i in 2L:nsim) out[i, ] <- sample.int(n)
  }
  out
}

bench_one <- function(tips = 500, nsim = 1000, cores = c(1, 2, 4)) {
  set.seed(1000 + tips)
  tree <- ape::rtree(tips)
  x <- stats::rnorm(tips)
  names(x) <- tree$tip.label
  xb <- as.integer(x > stats::median(x))
  names(xb) <- tree$tip.label

  perms <- make_perms(tips, nsim)
  random_states <- matrix(xb[t(perms)], nrow = tips, ncol = nsim)
  C <- unclass(caper::VCV.array(tree))
  samples <- t(chol(C)) %*% matrix(stats::rnorm(tips * nsim), nrow = tips)
  prop_state1 <- mean(xb == sort(unique(xb))[[1]])
  brownian_states <- fastphylosig:::brownian_threshold_cpp(samples, prop_state1)

  rows <- list()
  for (nc in cores) {
    rows[[length(rows) + 1L]] <- data.frame(
      statistic = "K_test",
      tips = tips,
      nsim = nsim,
      ncores = nc,
      elapsed_sec = time_expr(
        fast_signal(
          tree, x, method = "K", test = TRUE, nsim = nsim,
          permutations = perms, return_sim = FALSE, verbose = FALSE,
          ncores = nc
        )
      )
    )
    rows[[length(rows) + 1L]] <- data.frame(
      statistic = "D",
      tips = tips,
      nsim = nsim,
      ncores = nc,
      elapsed_sec = time_expr(
        fast_d(
          tree, xb, nsim = nsim, random_states = random_states,
          brownian_states = brownian_states, return_sim = FALSE,
          verbose = FALSE, ncores = nc
        )
      )
    )
  }
  out <- do.call(rbind, rows)
  out$speedup_vs_1core <- ave(
    out$elapsed_sec, interaction(out$statistic, out$tips, drop = TRUE),
    FUN = function(z) z[1L] / z
  )
  out
}

tips <- as.integer(strsplit(Sys.getenv("FASTPHYLOSIG_BENCH_TIPS", "100,500"),
                            ",", fixed = TRUE)[[1]])
nsim <- as.integer(Sys.getenv("FASTPHYLOSIG_BENCH_NSIM", "1000"))
cores <- as.integer(strsplit(Sys.getenv("FASTPHYLOSIG_BENCH_CORES", "1,2,4"),
                             ",", fixed = TRUE)[[1]])

result <- do.call(rbind, lapply(tips, bench_one, nsim = nsim, cores = cores))
print(result)

out_dir <- file.path("benchmarks", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  result,
  file.path(out_dir, "multicore_single_trait_benchmark.csv"),
  row.names = FALSE
)

png(
  file.path(out_dir, "multicore_single_trait_speedup.png"),
  width = 1400, height = 800, res = 150
)
op <- par(mfrow = c(1, length(unique(result$statistic))), mar = c(5, 4, 3, 1))
for (stat in unique(result$statistic)) {
  z <- result[result$statistic == stat, ]
  labels <- paste0(z$tips, " tips\n", z$ncores, " cores")
  barplot(
    z$speedup_vs_1core,
    names.arg = labels,
    ylim = c(0, max(result$speedup_vs_1core, na.rm = TRUE) * 1.15),
    ylab = "speedup vs 1 core",
    main = stat,
    col = "#4C78A8",
    las = 1
  )
  abline(h = 1, lty = 2, col = "gray40")
}
par(op)
dev.off()
