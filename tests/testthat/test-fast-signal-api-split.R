# Unified dispatcher and fixed-method specialist API -------------------------

test_that("continuous specialists and dispatcher preserve result contracts", {
  set.seed(901)
  tree <- ape::rtree(12)
  x <- stats::rnorm(12)
  names(x) <- tree$tip.label

  k <- suppressMessages(fast_k(tree, data = x, verbose = FALSE))
  high_k <- suppressMessages(
    fast_signal(tree, data = x, method = "k", verbose = FALSE)
  )
  expect_s3_class(k, "phylosig")
  expect_s3_class(high_k, "phylosig")
  expect_equal(as.numeric(high_k), as.numeric(k), tolerance = 1e-12)
  expect_equal(attr(high_k, "method"), "K")
  expect_equal(attr(high_k, "workflow")$production_function, "fast_k")

  lambda <- suppressMessages(fast_lambda(tree, x = x, verbose = FALSE))
  high_lambda <- suppressMessages(
    fast_signal(tree, x = x, method = "Lambda", verbose = FALSE)
  )
  expect_equal(high_lambda$lambda, lambda$lambda, tolerance = 1e-10)
  expect_equal(attr(high_lambda, "workflow")$production_function,
               "fast_lambda")
})

test_that("dispatcher routes D and Delta without nesting target results", {
  set.seed(902)
  tree <- ape::rtree(10)
  z <- as.integer(stats::rnorm(10) > 0)
  names(z) <- tree$tip.label

  d <- suppressMessages(
    fast_signal(tree, data = z, method = "D", nsim = 2,
                return_sim = FALSE, verbose = FALSE)
  )
  expect_s3_class(d, "phylo.d")
  expect_false(is.list(d$DEstimate))
  expect_equal(attr(d, "workflow")$production_function, "fast_d")

  cat <- rep(c("a", "b", "c"), length.out = 10)
  names(cat) <- tree$tip.label
  delta <- suppressWarnings(suppressMessages(
    fast_signal(tree, data = cat, method = "delta", test = FALSE,
                mcmc_sim = 40, thin = 10, burn = 10, verbose = FALSE)
  ))
  expect_s3_class(delta, "phylo_delta")
  expect_true(is.numeric(delta$delta))
  expect_equal(attr(delta, "workflow")$production_function, "fast_delta")
})

test_that("method and argument validation is explicit", {
  tree <- ape::rtree(4)
  x <- setNames(stats::rnorm(4), tree$tip.label)

  expect_error(fast_signal(tree, data = x), "method")
  expect_error(fast_signal(tree, data = x, method = c("K", "lambda")),
               "exactly one")
  expect_error(fast_signal(tree, data = x, method = "unknown"),
               "Unknown method")
  expect_error(fast_signal(tree, data = x, method = "D",
                           lambda_profile = TRUE), "not applicable")
  expect_error(fast_signal(tree, data = x, x = x, method = "K"),
               "only one")
})
