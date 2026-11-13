test_that("sharp boundary LR endpoints use likelihood rather than linear interpolation", {
  tr <- ape::read.tree(text = "((a:1e-6,b:1e-6):.999999,(c:1e-6,d:.100001):.999999);")
  x <- c(a = -1, b = -1, c = 1, d = 1)
  V <- ape::vcv(tr); y <- x[rownames(V)]
  ell <- function(l) {
    C <- V * l; diag(C) <- diag(V); iv <- solve(C)
    mu <- sum(iv %*% y) / sum(iv); r <- y - mu
    s2 <- drop(t(r) %*% iv %*% r) / length(y)
    -.5 * (length(y) * (log(2*pi) + 1 + log(s2)) +
      as.numeric(determinant(C, logarithm = TRUE)$modulus))
  }
  cut <- ell(1) - .5 * qchisq(.95, 1)
  root <- uniroot(function(z) ell(z) - cut, c(.99, 1), tol = 1e-13)$root
  fits <- lapply(c(101L, 1001L), function(n) fast_lambda(
    tr, x, test = TRUE, lambda_profile_points = n, progress = FALSE))
  for (fit in fits) {
    expect_equal(unname(fit$lambda_CI), c(root, 1), tolerance = 1e-10)
    expect_lt(abs(ell(fit$lambda_CI[[1]]) - cut), 1e-6)
  }
  expect_equal(fits[[1]]$lambda_CI, fits[[2]]$lambda_CI, tolerance = 1e-10)
  mat <- cbind(a = x, b = x * 2, missing = c(-1, -1, 1, NA), bad = rep(1, 4))
  tab <- fast_lambda(tr, mat, test = TRUE, lambda_profile = TRUE, progress = FALSE)
  expect_equal(tab$lambda_CI_lower_fast[1:2], rep(root, 2), tolerance = 1e-10)
  expect_true(is.na(tab$lambda_CI_lower_fast[4]))
  prepared <- prepare_tree(tr)
  fit <- fast_lambda(prepared, x[c(4, 1, 3, 2)], test = TRUE, progress = FALSE)
  expect_equal(fit$lambda_CI, fits[[1]]$lambda_CI, tolerance = 1e-10)
})

test_that("CI refinement respects boundary, flat and non-finite brackets", {
  ci_fun <- fastphylosig:::.lambda_profile_ci
  cut <- -.5*qchisq(.95, 1)
  p <- function(z) data.frame(lambda = c(0, .5, 1), logL = cut + z)
  calls <- 0L
  ell <- function(z) { calls <<- calls + 1L; cut + 1 }
  expect_equal(ci_fun(p(c(1,1,1)), .5, 0, loglik = ell)$intervals,
               data.frame(lower = 0, upper = 1))
  expect_equal(calls, 0L)
  expect_equal(ci_fun(p(c(-1,0,-1)), .5, 0, loglik = ell)$intervals,
               data.frame(lower = .5, upper = .5))
  expect_equal(calls, 0L)
  expect_equal(ci_fun(p(c(1,NA,1)), 0, 0, loglik = ell)$intervals,
               data.frame(lower = c(0,1), upper = c(0,1)))
  expect_equal(calls, 0L)
  expect_warning(ci <- ci_fun(p(c(-1,1,1)), .5, 0,
                              loglik = function(z) NA_real_), "refinement failed")
  expect_equal(ci$lower, .25)
})

test_that("all detected disconnected components receive nonlinear root refinement", {
  cut <- -.5*qchisq(.95, 1)
  ell <- function(z) cut + (z-.1)*(z-.4)*(z-.6)*(z-.9)
  grid <- seq(0, 1, length.out=5)
  ci <- fastphylosig:::.lambda_profile_ci(
    data.frame(lambda=grid, logL=ell(grid)), 0, 0, loglik=ell)
  expect_equal(ci$intervals, data.frame(lower=c(0,.4,.9), upper=c(.1,.6,1)),
               tolerance=1e-10)
  expect_true(ci$disconnected)
})
