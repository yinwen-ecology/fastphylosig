test_that("lambda searches distinct likelihood modes and legal endpoints", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((t4:4.602493115,(t1:0.3598982585,t3:1.399696597):1.633107075):8.712189423,t2:0.8149253698);")
  x <- c(1.34191599368243, -1.14167599611392, -1.47283661624142, .745508116974182)
  C <- ape::vcv(tree)
  dense_ll <- function(lambda, y) {
    V <- lambda * C
    diag(V) <- diag(C)
    W <- solve(V)
    a <- sum(W %*% y) / sum(W)
    q <- sum((y - a) * (W %*% (y - a)))
    -.5 * (length(y) * (1 + log(2*pi) + log(q/length(y))) +
             as.numeric(determinant(V, logarithm = TRUE)$modulus))
  }
  compiled <- fastphylosig:::compile_tree_cpp(
    matrix(as.integer(tree$edge), ncol=2), tree$edge.length, ape::Ntip(tree))
  X <- cbind(x, x * 3 + 2, c(1, -1, 1, -1))
  got <- fastphylosig:::fast_lambda_tree_optimize_cpp(
    compiled, X, max_lambda=1, profile=TRUE, tol=1e-9)
  grid <- seq(0, 1, length.out=10001)
  for (j in seq_len(ncol(X))) {
    vals <- vapply(grid, dense_ll, numeric(1), y=X[,j])
    expect_gte(got$logLik[j] + 1e-8, max(vals))
    expect_equal(got$logLik[j], dense_ll(got$lambda[j], X[,j]), tolerance=1e-8)
    expect_equal(got$LR[j], 2*(got$logLik[j]-dense_ll(0, X[,j])), tolerance=1e-8)
    expect_equal(got$P[j], pchisq(got$LR[j], 1, lower.tail=FALSE), tolerance=1e-8)
    expect_gte(got$logLik[j] + 1e-8, max(got$profile$logLik[,j]))
    solo <- fastphylosig:::fast_lambda_tree_optimize_cpp(
      compiled, X[,j,drop=FALSE], max_lambda=1, tol=1e-9)
    expect_equal(got$lambda[j], solo$lambda[1], tolerance=1e-10)
  }
  expect_equal(got$lambda[1], 1)
})

test_that("lambda agrees with dense multistart searches across 20 random trees", {
  skip_if_not_installed("ape")
  for (seed in seq_len(20)) {
    set.seed(73000 + seed)
    tree <- ape::rtree(8)
    tree$edge.length <- 10^runif(nrow(tree$edge), -2, 2)
    C <- ape::vcv(tree)
    X <- matrix(rnorm(24), nrow=8)
    X[,2] <- X[,2] * 10^runif(8, -1, 1)
    compiled <- fastphylosig:::compile_tree_cpp(
      matrix(as.integer(tree$edge), ncol=2), tree$edge.length, 8L)
    got <- fastphylosig:::fast_lambda_tree_optimize_cpp(
      compiled, X, max_lambda=1, tol=1e-8)
    dense <- function(z, y) {
      V <- z*C
      diag(V) <- diag(C)
      W <- solve(V)
      a <- sum(W %*% y)/sum(W)
      q <- sum((y-a)*(W %*% (y-a)))
      -.5*(8*(1+log(2*pi)+log(q/8)) +
        as.numeric(determinant(V, logarithm=TRUE)$modulus))
    }
    grid <- sort(unique(c(seq(0,1,length.out=1001), 2^-(1:35), 1-2^-(1:35))))
    for (j in 1:3) {
      fun <- function(z) dense(z, X[,j])
      vals <- vapply(grid, fun, numeric(1))
      n <- length(grid)
      peaks <- which(vals[2:(n-1)] >= vals[1:(n-2)] &
                     vals[2:(n-1)] >= vals[3:n]) + 1L
      refined <- vapply(peaks, function(k)
        optimize(fun, grid[c(k-1L,k+1L)], maximum=TRUE,
                 tol=1e-10)$objective, numeric(1))
      expect_equal(got$logLik[j], max(c(vals, refined)), tolerance=1e-7,
                   info=paste("seed", seed, "trait", j))
      expect_equal(got$logLik[j], fun(got$lambda[j]), tolerance=1e-8)
    }
  }
})

test_that("pooled ML and REML reference resolves narrow boundary modes", {
  set.seed(89378)
  tree <- ape::rtree(8)
  tree$edge.length <- 10^runif(nrow(tree$edge), -4, 4)
  X <- matrix(rnorm(24), 8)
  X[,2] <- X[,2]*10^runif(8, -2, 2)
  y <- X[,1]
  C <- ape::vcv(tree)
  core <- fastphylosig:::compile_tree_cpp(tree$edge, tree$edge.length, 8L)
  for (reml in c(FALSE, TRUE)) {
    dense <- function(z) {
      V <- z*C; diag(V) <- diag(C)
      W <- solve(V)
      a <- sum(W %*% y)/sum(W)
      q <- sum((y-a)*(W %*% (y-a)))
      df <- 8-as.integer(reml)
      # Two identical crowns: twice the separately profiled likelihood.
      -(df*(1+log(2*pi)+log(q/df)) +
        as.numeric(determinant(V, logarithm=TRUE)$modulus) +
        if (reml) log(sum(W)) else 0)
    }
    grid <- sort(unique(c(seq(0,1,length.out=4001),1-2^-(1:35))))
    vals <- vapply(grid, dense, numeric(1)); n <- length(grid)
    peaks <- which(vals[2:(n-1)] >= vals[1:(n-2)] &
                   vals[2:(n-1)] >= vals[3:n]) + 1L
    optima <- vapply(peaks, function(k)
      optimize(dense, grid[c(k-1L,k+1L)], maximum=TRUE,
               tol=1e-10)$objective, numeric(1))
    got <- fastphylosig:::heterogeneity_lambda_calibrate_cpp(
      list(core,core), list(y,y), .2, 2L, FALSE, reml)
    expect_equal(got$null_logLik, max(c(vals,optima)), tolerance=1e-8)
    expect_equal(got$null_logLik, dense(got$lambda_common), tolerance=1e-8)
    expect_lt(got$lambda_common, 1)
    if (!reml) {
      single <- fastphylosig:::fast_lambda_tree_optimize_cpp(
        core, matrix(y,ncol=1), tol=1e-7)
      expect_equal(got$null_logLik, 2*single$logLik[1], tolerance=1e-8)
      expect_equal(got$lambda_common, single$lambda[1], tolerance=1e-7)
    }
  }
})
