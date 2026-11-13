test_that("shared lambda grid weights preserve per-trait optimizer results", {
  set.seed(20261112)
  for (n in c(8L, 64L)) for (scale in c(1e-100, 1, 1e100)) {
    tree <- ape::rtree(n)
    tree$edge.length <- tree$edge.length * scale
    compiled <- compile_tree_cpp(matrix(as.integer(tree$edge), ncol = 2L),
                                 tree$edge.length, n)
    X <- cbind(stats::rnorm(n), stats::rnorm(n) * 1e50,
               stats::rnorm(n) * 1e-50, rep(1, n))
    before <- .Random.seed
    batch <- fast_lambda_tree_optimize_cpp(compiled, X, test = TRUE,
                                          profile = TRUE, profile_points = 15L)
    expect_identical(.Random.seed, before)
    for (j in seq_len(ncol(X))) {
      single <- fast_lambda_tree_optimize_cpp(compiled, X[, j, drop=FALSE],
                                             test=TRUE, profile=FALSE)
      for (field in c("lambda", "logLik", "gls_mean", "sigma2", "logLik0",
                      "LR", "P", "valid", "status")) {
        expect_identical(batch[[field]][j], single[[field]][1L])
      }
      fixed <- fast_lambda_tree_fixed_cpp(compiled, X[,j,drop=FALSE],
                                         batch$profile$lambda)
      expect_equal(batch$profile$logLik[,j], fixed$logLik[,1], tolerance=0)
    }
  }
})