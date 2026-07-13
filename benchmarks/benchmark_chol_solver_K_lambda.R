# Benchmark Cholesky solver changes for K and lambda.
#
# K comparison:
#   old path = ape::vcv.phylo + solve(C) + original invC K kernel
#   new path = ape::vcv.phylo + chol(C) + Cholesky K kernel
#
# Lambda comparison:
#   fixed likelihood compares old explicit solve/determinant in R against the
#   current C++ Cholesky likelihood. Full optimization compares the current
#   fast_signal(lambda) against phytools on moderate trees.

library(ape)
library(fastphylosig)

out_dir <- file.path("benchmarks", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

time_call <- function(expr) {
  gc()
  elapsed <- system.time(value <- force(expr))[["elapsed"]]
  list(value = value, elapsed = unname(elapsed))
}

old_k_from_inv <- function(C, x) {
  invC <- solve(C)
  n <- nrow(C)
  sum_invC <- sum(invC)
  a <- sum(invC %*% x) / sum_invC
  centered <- x - a
  norm_const <- (sum(diag(C)) - n / sum_invC) / (n - 1)
  as.numeric(
    (sum(centered * centered) /
       drop(t(centered) %*% invC %*% centered)) / norm_const
  )
}

new_k_from_chol <- function(C, x) {
  cholC <- chol(C)
  fastphylosig:::fast_k_chol_batch_cpp(
    matrix(x, ncol = 1), cholC, sum(diag(C))
  )[[1]]
}

lambda_loglik_solve_R <- function(lambda, C, y) {
  n <- nrow(C)
  Cl <- lambda * C
  diag(Cl) <- diag(C)
  invCl <- solve(Cl)
  sum_invCl <- sum(invCl)
  a <- sum(invCl %*% y) / sum_invCl
  centered <- y - a
  quad <- drop(t(centered) %*% invCl %*% centered)
  sig2 <- quad / n
  logdetCl <- as.numeric(determinant(Cl, logarithm = TRUE)$modulus)
  -0.5 * quad / sig2 -
    0.5 * n * log(2 * pi) -
    0.5 * (n * log(sig2) + logdetCl)
}

bench_k <- function(tips) {
  set.seed(6000 + tips)
  tree <- ape::rtree(tips)
  x <- ape::rTraitCont(tree)

  t_vcv <- time_call(ape::vcv.phylo(tree))
  C <- t_vcv$value
  t_old <- time_call(old_k_from_inv(C, x))
  t_new <- time_call(new_k_from_chol(C, x))

  data.frame(
    statistic = "K",
    tips = tips,
    vcv_sec = t_vcv$elapsed,
    old_solver_sec = t_old$elapsed,
    new_chol_sec = t_new$elapsed,
    old_total_sec = t_vcv$elapsed + t_old$elapsed,
    new_total_sec = t_vcv$elapsed + t_new$elapsed,
    speedup_total = (t_vcv$elapsed + t_old$elapsed) /
      (t_vcv$elapsed + t_new$elapsed),
    old_value = t_old$value,
    new_value = t_new$value,
    abs_error = abs(t_old$value - t_new$value),
    stringsAsFactors = FALSE
  )
}

bench_lambda_fixed <- function(tips, lambda = 0.5) {
  set.seed(7000 + tips)
  tree <- ape::rtree(tips)
  y <- ape::rTraitCont(tree)
  C <- ape::vcv.phylo(tree)

  t_old <- time_call(lambda_loglik_solve_R(lambda, C, y))
  t_new <- time_call(fastphylosig:::lambda_loglik_cpp(lambda, C, y))

  data.frame(
    statistic = "lambda_fixed_loglik",
    tips = tips,
    lambda = lambda,
    old_solve_sec = t_old$elapsed,
    new_chol_sec = t_new$elapsed,
    speedup = t_old$elapsed / t_new$elapsed,
    old_loglik = t_old$value,
    new_loglik = t_new$value,
    abs_error = abs(t_old$value - t_new$value),
    stringsAsFactors = FALSE
  )
}

bench_lambda_full <- function(tips) {
  if (!requireNamespace("phytools", quietly = TRUE)) {
    return(NULL)
  }
  set.seed(8000 + tips)
  tree <- ape::rtree(tips)
  y <- ape::rTraitCont(tree)

  t_fast <- time_call(fast_signal(
    tree, y, method = "lambda", test = TRUE, verbose = FALSE
  ))
  t_ref <- time_call(phytools::phylosig(
    tree, y, method = "lambda", test = TRUE, se = NULL
  ))

  data.frame(
    statistic = "lambda_full_optimize",
    tips = tips,
    fast_sec = t_fast$elapsed,
    phytools_sec = t_ref$elapsed,
    speedup = t_ref$elapsed / t_fast$elapsed,
    lambda_fast = t_fast$value$lambda,
    lambda_phytools = t_ref$value$lambda,
    logL_fast = t_fast$value$logL,
    logL_phytools = t_ref$value$logL,
    P_fast = t_fast$value$P,
    P_phytools = t_ref$value$P,
    abs_lambda_error = abs(t_fast$value$lambda - t_ref$value$lambda),
    abs_logL_error = abs(t_fast$value$logL - t_ref$value$logL),
    abs_P_error = abs(t_fast$value$P - t_ref$value$P),
    stringsAsFactors = FALSE
  )
}

k_tips <- as.integer(strsplit(Sys.getenv("FASTPHYLOSIG_K_CHOL_TIPS",
                                         "1000,2000,5000"),
                              ",", fixed = TRUE)[[1]])
lambda_fixed_tips <- as.integer(strsplit(
  Sys.getenv("FASTPHYLOSIG_LAMBDA_FIXED_TIPS", "500,1000,2000"),
  ",", fixed = TRUE
)[[1]])
lambda_full_tips <- as.integer(strsplit(
  Sys.getenv("FASTPHYLOSIG_LAMBDA_FULL_TIPS", "100,300,500"),
  ",", fixed = TRUE
)[[1]])

k_res <- do.call(rbind, lapply(k_tips, bench_k))
lambda_fixed_res <- do.call(rbind, lapply(lambda_fixed_tips,
                                          bench_lambda_fixed))
lambda_full_res <- do.call(rbind, lapply(lambda_full_tips,
                                         bench_lambda_full))

write.csv(k_res, file.path(out_dir, "chol_solver_K_benchmark.csv"),
          row.names = FALSE)
write.csv(lambda_fixed_res,
          file.path(out_dir, "chol_solver_lambda_fixed_benchmark.csv"),
          row.names = FALSE)
if (!is.null(lambda_full_res)) {
  write.csv(lambda_full_res,
            file.path(out_dir, "chol_solver_lambda_full_benchmark.csv"),
            row.names = FALSE)
}

png(file.path(out_dir, "chol_solver_K_lambda_speedup.png"),
    width = 1500, height = 800, res = 150)
op <- par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))
barplot(k_res$speedup_total, names.arg = paste0(k_res$tips, " tips"),
        main = "K total", ylab = "speedup", col = "#4C78A8", las = 1,
        ylim = c(0, max(k_res$speedup_total, lambda_fixed_res$speedup,
                        if (!is.null(lambda_full_res)) lambda_full_res$speedup,
                        na.rm = TRUE) * 1.15))
abline(h = 1, lty = 2, col = "gray40")
barplot(lambda_fixed_res$speedup,
        names.arg = paste0(lambda_fixed_res$tips, " tips"),
        main = "lambda fixed logLik", ylab = "speedup", col = "#4C78A8",
        las = 1)
abline(h = 1, lty = 2, col = "gray40")
if (!is.null(lambda_full_res)) {
  barplot(lambda_full_res$speedup,
          names.arg = paste0(lambda_full_res$tips, " tips"),
          main = "lambda optimize", ylab = "speedup", col = "#4C78A8",
          las = 1)
  abline(h = 1, lty = 2, col = "gray40")
} else {
  plot.new()
  title("lambda optimize")
  text(0.5, 0.5, "phytools unavailable")
}
par(op)
dev.off()

print(k_res)
print(lambda_fixed_res)
if (!is.null(lambda_full_res)) print(lambda_full_res)
