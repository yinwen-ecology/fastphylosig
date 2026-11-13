test_that("profile LR set preserves all connected components", {
  cutoff <- -0.5 * qchisq(.95, 1)
  prof <- data.frame(lambda = seq(0, 1, length.out = 5),
                     logL = cutoff + c(1, -1, -1, 1, 1))
  ci <- fastphylosig:::.lambda_profile_ci(prof, 1, 0)
  expect_equal(ci$intervals, data.frame(lower = c(0, .625), upper = c(.125, 1)))
  expect_equal(c(ci$lower, ci$upper), c(0, 1))
  expect_true(ci$disconnected)
  expect_equal(ci$cutoff, cutoff)
})

test_that("profile LR set handles connected, flat, boundary and missing profiles", {
  ci_fun <- fastphylosig:::.lambda_profile_ci
  cut <- -.5 * qchisq(.95, 1)
  p <- function(z) data.frame(lambda = seq(0, 1, length.out = length(z)), logL = cut + z)
  ci <- ci_fun(p(c(-1, 1, -1)), .5, 0)
  expect_equal(ci$intervals, data.frame(lower = .25, upper = .75))
  expect_false(ci$disconnected)
  flat <- ci_fun(p(rep(1, 5)), .5, 0)
  expect_equal(flat$intervals, data.frame(lower = 0, upper = 1))
  expect_false(flat$disconnected)
  point <- ci_fun(p(c(-1, 0, -1)), .5, 0)
  expect_equal(point$intervals, data.frame(lower = .5, upper = .5))
  gap <- ci_fun(p(c(1, NA, 1)), 0, 0)
  expect_equal(gap$intervals, data.frame(lower = c(0, 1), upper = c(0, 1)))
  expect_true(gap$disconnected)
  bad <- ci_fun(data.frame(lambda = c(0, 1), logL = c(-Inf, -Inf)), NA_real_, -Inf)
  expect_true(is.na(bad$lower))
  expect_equal(nrow(bad$intervals), 0L)
})

test_that("multimodal lambda CI retains zero and plots separated intervals", {
  set.seed(9123)
  for (i in 1:4) {
    tr <- ape::rtree(8)
    tr$edge.length <- 10^runif(nrow(tr$edge), -3, 3)
    for (j in 1:30) x <- setNames(rnorm(8) * 10^runif(8, -1, 1), tr$tip.label)
  }
  fit <- fast_lambda(tr, x, test = TRUE, lambda_profile_points = 1001, progress = FALSE)
  expect_gt(fit$P, .05)
  expect_equal(unname(fit$lambda_CI), c(0, 1))
  expect_true(fit$lambda_CI_disconnected)
  ints <- fit$lambda_CI_intervals
  expect_equal(nrow(ints), 2L)
  expect_equal(ints$lower[1], 0)
  expect_equal(ints$upper[2], 1)
  # Independent dense Gaussian likelihood, rather than another package search.
  V <- ape::vcv(tr); y <- x[rownames(V)]
  ell <- function(l) {
    C <- V * l; diag(C) <- diag(V); iv <- solve(C)
    mu <- sum(iv %*% y) / sum(iv); r <- y - mu
    s2 <- drop(t(r) %*% iv %*% r) / length(y)
    -.5 * (length(y) * (log(2 * pi) + 1 + log(s2)) +
      as.numeric(determinant(C, logarithm = TRUE)$modulus))
  }
  cut <- fit$lambda_CI_cutoff
  expect_gt(ell(0), cut)
  expect_lt(ell(.1), cut)
  roots <- c(uniroot(function(z) ell(z) - cut, c(0, .1), tol = 1e-12)$root,
             uniroot(function(z) ell(z) - cut, c(.1, 1), tol = 1e-12)$root)
  expect_equal(c(ints$upper[1], ints$lower[2]), roots, tolerance = 2e-5)
  mat <- cbind(a = x, b = x)
  tab <- fast_lambda(tr, mat, test = TRUE, lambda_profile = TRUE, lambda_profile_points = 1001, progress = FALSE)
  expect_true(all(tab$lambda_CI_disconnected_fast))
  expect_equal(tab$lambda_CI_intervals_fast[[1]], ints)
  expect_equal(tab$lambda_CI_intervals_fast[[2]], ints)
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  drawn <- plot_signal(fit)
  grDevices::dev.off()
  unlink(file)
  expect_equal(drawn$ci_intervals[[1]], ints)
})