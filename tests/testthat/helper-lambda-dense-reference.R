# Independent covariance oracle: one eigendecomposition of D^(-1/2) C D^(-1/2)
# evaluates Pagel's likelihood in O(n) per lambda, with no production kernels.
.lambda_dense_reference <- function(tree, y) {
  C <- ape::vcv.phylo(tree)
  if (!is.null(names(y))) y <- y[rownames(C)]
  d <- diag(C)
  eig <- eigen(C / sqrt(outer(d, d)), symmetric = TRUE)
  origin <- mean(y)
  z <- drop(crossprod(eig$vectors, (y - origin) / sqrt(d)))
  one <- drop(crossprod(eig$vectors, 1 / sqrt(d)))
  n <- length(y)
  evaluate <- function(lambda) {
    v <- 1 + lambda * (eig$values - 1)
    if (any(!is.finite(v) | v <= 0)) {
      return(c(logL = -Inf, gls_mean = NA_real_, sigma2 = NA_real_))
    }
    a <- sum(one * z / v) / sum(one^2 / v)
    sigma2 <- sum((z - a * one)^2 / v) / n
    ll <- -.5 * (n * (1 + log(2*pi) + log(sigma2)) +
                   sum(log(d)) + sum(log(v)))
    c(logL = ll, gls_mean = origin + a, sigma2 = sigma2)
  }
  depth <- ape::node.depth.edgelength(tree)
  terminal <- tree$edge[, 2] <= n
  parent_depth <- depth[tree$edge[terminal, 1]]
  upper <- if (ape::is.ultrametric(tree) && any(parent_depth > 0)) {
    min(depth[tree$edge[terminal, 2]][parent_depth > 0] /
          parent_depth[parent_depth > 0])
  } else 1
  grid <- sort(unique(c(seq(0, upper, length.out = 2049),
                        upper * 2^-(1:35), upper * (1 - 2^-(1:35)))))
  objective <- function(lambda) evaluate(lambda)[['logL']]
  values <- vapply(grid, objective, numeric(1))
  m <- length(grid)
  peaks <- which(values[2:(m-1)] >= values[1:(m-2)] &
                   values[2:(m-1)] >= values[3:m] &
                   (values[2:(m-1)] > values[1:(m-2)] |
                      values[2:(m-1)] > values[3:m])) + 1L
  refined <- vapply(peaks, function(k) stats::optimize(
    objective, grid[c(k-1L, k+1L)], maximum = TRUE, tol = 1e-10
  )$maximum, numeric(1))
  candidates <- c(grid, refined)
  lambda <- candidates[which.max(vapply(candidates, objective, numeric(1)))]
  best <- evaluate(lambda)
  logL0 <- objective(0)
  LR <- max(0, 2 * (best[['logL']] - logL0))
  list(lambda = lambda, logL = best[['logL']],
       gls_mean = best[['gls_mean']], sigma2 = best[['sigma2']],
       logL0 = logL0, LR = LR, P = stats::pchisq(LR, 1, lower.tail = FALSE),
       evaluate = evaluate)
}
