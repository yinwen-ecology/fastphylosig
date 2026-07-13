#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(fastphylosig)
})

out_dir <- file.path(getwd(), "benchmarks", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tips <- as.integer(strsplit(Sys.getenv("FAST_ACE_TIPS", "5000,10000,20000"), ",")[[1]])
mcmc_sim <- as.integer(Sys.getenv("FAST_ACE_MCMC_SIM", "1000"))
thin <- as.integer(Sys.getenv("FAST_ACE_THIN", "10"))
burn <- as.integer(Sys.getenv("FAST_ACE_BURN", "100"))

scale_tree <- function(tree) {
  depth <- ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)]
  tree$edge.length <- tree$edge.length / max(depth)
  tree
}

qbase <- matrix(
  c(0, 1.2, 0.25,
    0.45, 0, 1.8,
    0.15, 0.7, 0),
  nrow = 3, byrow = TRUE
)

bench_one <- function(n, model) {
  set.seed(100 + n)
  tree <- scale_tree(ape::rtree(n))
  if (model == "ER") {
    x <- ape::rTraitDisc(
      tree, model = "ER", k = 3, rate = 2, states = letters[1:3]
    )
  } else {
    x <- ape::rTraitDisc(tree, model = qbase, states = letters[1:3])
  }
  names(x) <- tree$tip.label

  gc()
  t_fast_delta <- system.time({
    set.seed(1)
    delta_fast <- suppressWarnings(fast_delta(
      tree, x, model = model, ace_engine = "fast", test = FALSE,
      mcmc_sim = mcmc_sim, thin = thin, burn = burn, verbose = FALSE
    ))
  })[["elapsed"]]

  gc()
  t_ape_delta <- system.time({
    set.seed(1)
    delta_ape <- suppressWarnings(fast_delta(
      tree, x, model = model, ace_engine = "ape", test = FALSE,
      mcmc_sim = mcmc_sim, thin = thin, burn = burn, verbose = FALSE
    ))
  })[["elapsed"]]

  y <- factor(x[tree$tip.label])
  gc()
  t_fast_ace <- system.time({
    ace_fast <- fast_ace(y, tree, model = model)
  })[["elapsed"]]

  gc()
  t_ape_ace <- system.time({
    ace_ape <- suppressWarnings(ape::ace(
      y, tree, type = "discrete", method = "ML", model = model
    ))
  })[["elapsed"]]

  data.frame(
    tips = n,
    model = model,
    fast_delta_fastace_s = t_fast_delta,
    fast_delta_apeake_s = t_ape_delta,
    fast_delta_speedup = t_ape_delta / t_fast_delta,
    fast_ace_s = t_fast_ace,
    ape_ace_s = t_ape_ace,
    ace_speedup = t_ape_ace / t_fast_ace,
    delta_fast = as.numeric(delta_fast$delta),
    delta_ape = as.numeric(delta_ape$delta),
    delta_abs_diff = abs(as.numeric(delta_fast$delta) -
                           as.numeric(delta_ape$delta)),
    loglik_abs_diff = abs(ace_fast$loglik - ace_ape$loglik),
    lik_anc_max_abs_diff = max(abs(ace_fast$lik.anc - ace_ape$lik.anc)),
    stringsAsFactors = FALSE
  )
}

result <- do.call(rbind, lapply(tips, function(n) {
  do.call(rbind, lapply(c("ER", "ARD"), function(model) bench_one(n, model)))
}))

result_file <- file.path(out_dir, "fast_ace_full_workflow_benchmark.csv")
utils::write.csv(result, result_file, row.names = FALSE)

png(file.path(out_dir, "fast_ace_full_workflow_speed.png"),
    width = 1700, height = 950, res = 150)
op <- par(mfrow = c(1, 2), mar = c(5, 5, 3, 1), las = 1)
for (model in c("ER", "ARD")) {
  z <- result[result$model == model, ]
  ylim <- range(c(z$fast_delta_apeake_s, z$fast_delta_fastace_s))
  plot(z$tips, z$fast_delta_apeake_s, type = "b", pch = 16, lwd = 2,
       log = "xy", col = "#E45756", xlab = "tips", ylab = "seconds",
       ylim = c(min(ylim) * 0.7, max(ylim) * 1.3),
       main = paste0(model, ": full fast_delta workflow"))
  lines(z$tips, z$fast_delta_fastace_s, type = "b", pch = 16, lwd = 2,
        col = "#4C78A8")
  legend("topleft", legend = c("ape::ace engine", "fast_ace engine"),
         col = c("#E45756", "#4C78A8"), pch = 16, lwd = 2, bty = "n")
}
par(op)
dev.off()

print(result)
message("Wrote: ", result_file)
