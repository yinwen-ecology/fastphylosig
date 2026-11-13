test_that("valid ACE information retains its standard-error arithmetic", {
  for (h in list(matrix(3, 1L), matrix(c(4, 1, 1, 3), 2L))) {
    for (unit in c(1, .01, 100)) {
      expected <- unit * sqrt(2 * diag(solve(h)))
      expect_identical(.fast_ace_rate_standard_errors(h, nrow(h), unit), expected)
    }
  }
})

test_that("invalid ACE curvature cannot leave selectively finite standard errors", {
  # This indefinite Hessian has ALL positive inverse diagonal entries.
  # Testing only diag(solve(H)) therefore misses the invalid covariance.
  q <- matrix(c(1, 1, 1, 1, 1, -1, 1, -1,
                1, 1, -1, -1, 1, -1, -1, 1), 4L) / 2
  indefinite <- q %*% diag(c(.1, 10, 10, -1)) %*% t(q)
  expect_true(all(diag(indefinite) > 0))
  expect_true(all(diag(solve(indefinite)) > 0))
  invalid <- list(indefinite, diag(c(1, -1)), matrix(1, 2L, 2L),
                  diag(c(1, 0)), matrix(NA_real_, 2L, 2L),
                  diag(c(1, Inf)), matrix(c(1, 0, 2, 1), 2L), NULL)
  for (h in invalid) {
    np <- if (is.matrix(h)) nrow(h) else 2L
    messages <- list()
    se <- withCallingHandlers(.fast_ace_rate_standard_errors(h, np),
      fastphylosig_ace_se_warning = function(w) {
        messages[[length(messages) + 1L]] <<- w
        invokeRestart("muffleWarning")
      })
    expect_length(messages, 1L)
    expect_match(conditionMessage(messages[[1L]]), "standard errors are unavailable")
    expect_true(all(is.nan(se)))
    expect_length(se, np)
  }
})

test_that("ACE curvature guard handles singular inversion and invalid scaling", {
  h <- diag(c(1, 1e-20)) # positive definite, but solve() rejects its conditioning
  expect_warning(se <- .fast_ace_rate_standard_errors(h, 2L),
                 class = "fastphylosig_ace_se_warning")
  expect_true(all(is.nan(se)))
  expect_warning(se <- .fast_ace_rate_standard_errors(diag(2), 2L, Inf),
                 class = "fastphylosig_ace_se_warning")
  expect_true(all(is.nan(se)))
})

test_that("unavailable ACE standard errors do not change estimates or RNG", {
  set.seed(1941)
  tree <- ape::rtree(32L)
  x <- setNames(rep(c("a", "b"), 16L), tree$tip.label)
  before <- .Random.seed
  reference <- fast_ace(x, tree, model = "ER", progress = FALSE)
  calls <- 0L
  original_nlm <- stats::nlm
  local({
  testthat::local_mocked_bindings(nlm = function(...) {
    calls <<- calls + 1L
    list(hessian = matrix(-1, 1L))
  }, .package = "stats")
  expect_warning(actual <- fast_ace(x, tree, model = "ER", progress = FALSE),
                 class = "fastphylosig_ace_se_warning")
  expect_true(all(is.nan(actual$se)))
  expect_equal(actual$rates, reference$rates, tolerance = 1e-12)
  expect_equal(actual$loglik, reference$loglik, tolerance = 1e-12)
  expect_identical(actual$convergence, reference$convergence)
  expect_equal(actual$lik.anc, reference$lik.anc, tolerance = 1e-12)
  expect_identical(.Random.seed, before)
  expect_identical(calls, 1L)
  expect_no_warning(no_se <- fast_ace(x, tree, model = "ER", CI = FALSE,
                                     progress = FALSE))
  expect_identical(calls, 1L)
  expect_true(all(is.nan(no_se$se)))
  })
  expect_identical(stats::nlm, original_nlm)
})
