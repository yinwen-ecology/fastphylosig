#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(fastphylosig)
})

root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "benchmarks", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

python_bin <- Sys.getenv(
  "PYTHON_BIN",
  "C:/Users/wenyi/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
)
if (!file.exists(python_bin)) python_bin <- "python"
python_ref <- file.path(root, "benchmarks", "delta_python_reference.py")

tips <- as.integer(strsplit(Sys.getenv("DELTA_TIPS", "5000,10000,20000"), ",")[[1]])
mcmc_sim <- as.integer(Sys.getenv("DELTA_MCMC_SIM", "10000"))
thin <- as.integer(Sys.getenv("DELTA_THIN", "10"))
burn <- as.integer(Sys.getenv("DELTA_BURN", "100"))
repeats <- as.integer(Sys.getenv("DELTA_REPEATS", "3"))
lambda0 <- as.numeric(Sys.getenv("DELTA_LAMBDA0", "0.1"))
proposal_sd <- as.numeric(Sys.getenv("DELTA_PROPOSAL_SD", "0.5"))
entropy <- Sys.getenv("DELTA_ENTROPY", "LSE")
model_families <- strsplit(Sys.getenv("DELTA_MODEL_FAMILIES", "ER,ARD"), ",")[[1]]
model_families <- toupper(trimws(model_families))
seed0 <- as.integer(Sys.getenv("DELTA_SEED", "20260705"))

scenarios <- data.frame(
  scenario = c("slow_rate_strong_signal",
               "medium_rate_moderate_signal",
               "fast_rate_weak_signal"),
  er_rate = c(0.8, 2.0, 8.0),
  stringsAsFactors = FALSE
)

ard_base <- matrix(
  c(0.0, 1.2, 0.25,
    0.45, 0.0, 1.8,
    0.15, 0.7, 0.0),
  nrow = 3, byrow = TRUE
)
ard_base <- ard_base / mean(ard_base[row(ard_base) != col(ard_base)])

scale_tree_height <- function(tree) {
  depth <- ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)]
  tree$edge.length <- tree$edge.length / max(depth)
  tree
}

simulate_trait <- function(tree, rate, model_family, seed) {
  set.seed(seed)
  if (model_family == "ER") {
    x <- ape::rTraitDisc(
      tree, model = "ER", k = 3, rate = rate, states = letters[1:3]
    )
  } else if (model_family == "ARD") {
    x <- ape::rTraitDisc(
      tree, model = rate * ard_base, states = letters[1:3]
    )
  } else {
    stop("Unsupported model family: ", model_family, call. = FALSE)
  }
  names(x) <- tree$tip.label
  x
}

ace_prob <- function(tree, x, model) {
  y <- factor(x[tree$tip.label])
  ar <- suppressWarnings(
    ape::ace(y, tree, type = "discrete", method = "ML", model = model)$lik.anc
  )
  if (is.complex(ar)) ar <- Re(ar)
  ar <- as.matrix(ar)
  storage.mode(ar) <- "double"
  ar
}

time_expr <- function(expr) {
  gc()
  t0 <- Sys.time()
  value <- force(expr)
  list(value = value, time_s = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

run_python_reference <- function(prob_file, out_prefix, seed) {
  py_csv <- paste0(out_prefix, "_python_delta.csv")
  py_entropy_file <- paste0(out_prefix, "_python_entropy.csv")
  py_args <- c(
    shQuote(normalizePath(python_ref, winslash = "/", mustWork = TRUE)),
    "--prob", shQuote(normalizePath(prob_file, winslash = "/", mustWork = TRUE)),
    "--out", shQuote(normalizePath(py_csv, winslash = "/", mustWork = FALSE)),
    "--entropy-out", shQuote(normalizePath(py_entropy_file, winslash = "/", mustWork = FALSE)),
    "--lambda0", lambda0,
    "--se", proposal_sd,
    "--sim", mcmc_sim,
    "--thin", thin,
    "--burn", burn,
    "--entropy", entropy,
    "--repeats", repeats,
    "--seed", seed
  )
  t0 <- Sys.time()
  status <- system2(python_bin, py_args)
  process_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!identical(status, 0L)) {
    stop("Python benchmark failed for ", basename(out_prefix), call. = FALSE)
  }
  list(
    rows = utils::read.csv(py_csv, stringsAsFactors = FALSE),
    entropy = scan(py_entropy_file, quiet = TRUE, sep = ","),
    process_time_s = process_time
  )
}

plot_scenarios <- function(result) {
  png(file.path(out_dir, "delta_scenarios_speed.png"),
      width = 1900, height = 1050, res = 150)
  op <- par(mfrow = c(length(unique(result$model_family)), 1),
            mar = c(7, 5, 3, 1), las = 2)
  on.exit(par(op), add = TRUE)
  for (fam in unique(result$model_family)) {
    tmp <- result[result$model_family == fam, ]
    labels <- paste(tmp$tips, sub("_signal", "", tmp$scenario), sep = "\n")
    mat <- t(as.matrix(tmp[, c("fast_delta_full_time_s",
                               "fast_delta_core_median_s",
                               "python_core_median_s")]))
    barplot(pmax(mat, 1e-4), beside = TRUE, log = "y", names.arg = labels,
            ylab = "seconds (log scale)",
            col = c("#4C78A8", "#54A24B", "#E45756"),
            main = paste0("Delta runtime: ", fam,
                          " simulation and ACE model"))
    legend("topleft",
           legend = c("fast_delta full", "fast_delta Rcpp core",
                      "original Python core"),
           fill = c("#4C78A8", "#54A24B", "#E45756"), bty = "n")
  }
  dev.off()

  png(file.path(out_dir, "delta_scenarios_accuracy.png"),
      width = 1700, height = 1000, res = 150)
  op <- par(mfrow = c(length(unique(result$model_family)), 1),
            mar = c(5, 5, 3, 1), las = 1)
  on.exit(par(op), add = TRUE)
  cols <- c("#4C78A8", "#F58518", "#54A24B")
  for (fam in unique(result$model_family)) {
    fam_result <- result[result$model_family == fam, ]
    plot(NA, xlim = range(fam_result$tips),
         ylim = range(fam_result$delta_mean_abs_diff),
         log = "x", xlab = "tips",
         ylab = "absolute difference in mean Delta",
         main = paste0("Rcpp vs original Python Delta estimates: ", fam))
    for (i in seq_along(unique(fam_result$scenario))) {
      sc <- unique(fam_result$scenario)[i]
      tmp <- fam_result[fam_result$scenario == sc, ]
      lines(tmp$tips, tmp$delta_mean_abs_diff, type = "b",
            pch = 16, lwd = 2, col = cols[i])
    }
    legend("topright", legend = unique(fam_result$scenario), col = cols,
           pch = 16, lwd = 2, bty = "n")
    abline(h = 0, col = "grey70")
  }
  dev.off()
}

tree_cache <- new.env(parent = emptyenv())
rows <- list()
row_id <- 0L

for (n in tips) {
  key <- as.character(n)
  if (!exists(key, tree_cache, inherits = FALSE)) {
    set.seed(seed0 + n)
    assign(key, scale_tree_height(ape::rtree(n)), tree_cache)
  }
  tree <- get(key, tree_cache, inherits = FALSE)

  for (model_family in model_families) {
  for (s in seq_len(nrow(scenarios))) {
    scenario <- scenarios$scenario[s]
    rate <- scenarios$er_rate[s]
    seed <- seed0 + n * 100L + s * 10L + match(model_family, c("ER", "ARD"))
    message("Running Delta scenario benchmark: ", n, " tips, ",
            model_family, ", ", scenario)

    x <- simulate_trait(tree, rate, model_family, seed)
    state_counts <- table(x)
    ace_model <- model_family

    full <- time_expr(
      suppressWarnings(
        fast_delta(
          tree, x, test = FALSE, mcmc_sim = mcmc_sim,
          thin = thin, burn = burn, lambda0 = lambda0,
          proposal_sd = proposal_sd, entropy = entropy, model = ace_model,
          verbose = FALSE
        )
      )
    )

    prob_time <- time_expr(ace_prob(tree, x, ace_model))
    prob <- prob_time$value
    prefix <- file.path(
      out_dir,
      paste0("delta_scenario_", n, "_", model_family, "_", scenario)
    )
    prob_file <- paste0(prefix, "_prob.csv")
    utils::write.table(prob, prob_file, sep = ",", row.names = FALSE,
                       col.names = FALSE)

    entropy_code <- match(entropy, c("LSE", "SE", "GINI"))
    r_entropy <- fastphylosig:::delta_entropy_cpp(prob, entropy_code)
    r_values <- numeric(repeats)
    r_times <- numeric(repeats)
    for (i in seq_len(repeats)) {
      set.seed(seed + 1000L + i)
      core <- time_expr(
        fastphylosig:::delta_mcmc_cpp(
          probabilities = prob, lambda0 = lambda0,
          proposal_sd = proposal_sd, sim = mcmc_sim, thin = thin,
          burn = burn, entropy_type = entropy_code
        )
      )
      r_values[i] <- core$value$delta
      r_times[i] <- core$time_s
    }

    py <- run_python_reference(prob_file, prefix, seed + 2000L)

    r_mean <- mean(r_values)
    py_mean <- mean(py$rows$delta)
    r_sd <- stats::sd(r_values)
    py_sd <- stats::sd(py$rows$delta)
    mc_se_diff <- sqrt(r_sd^2 / repeats + py_sd^2 / repeats)
    if (!is.finite(mc_se_diff)) mc_se_diff <- NA_real_

    row_id <- row_id + 1L
    rows[[row_id]] <- data.frame(
      tips = n,
      model_family = model_family,
      scenario = scenario,
      simulation_model = model_family,
      simulation_rate = rate,
      ace_model = ace_model,
      nodes = nrow(prob),
      states = ncol(prob),
      state_counts = paste(names(state_counts), as.integer(state_counts),
                           sep = ":", collapse = ";"),
      mcmc_sim = mcmc_sim,
      thin = thin,
      burn = burn,
      repeats = repeats,
      fast_delta_full = as.numeric(full$value$delta),
      fast_delta_full_time_s = full$time_s,
      ace_time_s = prob_time$time_s,
      fast_delta_core_mean = r_mean,
      fast_delta_core_sd = r_sd,
      fast_delta_core_median_s = stats::median(r_times),
      python_core_mean = py_mean,
      python_core_sd = py_sd,
      python_core_median_s = stats::median(py$rows$time_s),
      python_process_time_s = py$process_time_s,
      python_first_call_s = py$rows$first_time_s[[1]],
      delta_mean_abs_diff = abs(r_mean - py_mean),
      relative_delta_diff = abs(r_mean - py_mean) / abs(py_mean),
      delta_mean_diff_over_mc_se = if (is.finite(mc_se_diff) && mc_se_diff > 0) {
        abs(r_mean - py_mean) / mc_se_diff
      } else {
        NA_real_
      },
      entropy_max_abs_diff = max(abs(as.numeric(r_entropy) - py$entropy)),
      core_speedup_python_over_rcpp =
        stats::median(py$rows$time_s) / stats::median(r_times),
      full_speedup_python_core_over_full =
        stats::median(py$rows$time_s) / full$time_s,
      stringsAsFactors = FALSE
    )

    utils::write.csv(do.call(rbind, rows),
                     file.path(out_dir, "delta_scenarios_benchmark_partial.csv"),
                     row.names = FALSE)
  }
  }
}

result <- do.call(rbind, rows)
result_file <- file.path(out_dir, "delta_scenarios_benchmark.csv")
utils::write.csv(result, result_file, row.names = FALSE)
plot_scenarios(result)

print(result)
message("Wrote: ", result_file)
message("Wrote: ", file.path(out_dir, "delta_scenarios_speed.png"))
message("Wrote: ", file.path(out_dir, "delta_scenarios_accuracy.png"))
