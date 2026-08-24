test_that("D and Delta reject simultaneously supplied x and X", {
  expect_error(
    fast_d(NULL, x = NULL, X = NULL, progress = FALSE),
    "only one of x or X"
  )
  expect_error(
    fast_delta(NULL, x = NULL, X = NULL, progress = FALSE),
    "only one of x or X"
  )
})

test_that("fast_ace validates method as a scalar supported value", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- c(a = "A", b = "A", c = "B", d = "B")

  expect_error(
    fast_ace(x, tree, method = NA_character_, progress = FALSE),
    "single, non-missing character string"
  )
  expect_error(
    fast_ace(x, tree, method = c("ML", "ML"), progress = FALSE),
    "single, non-missing character string"
  )
  expect_error(
    fast_ace(x, tree, method = "REML", progress = FALSE),
    "supports only method = \"ML\""
  )
})

test_that("Delta ignores nsim when permutation testing is disabled", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- c(a = "A", b = "A", c = "B", d = "B")
  unused_nsim <- new.env(parent = emptyenv())

  fit <- suppressWarnings(fast_delta(
    tree, x, test = FALSE, nsim = unused_nsim,
    mcmc_sim = 2L, thin = 1L, burn = 1L,
    verbose = FALSE, progress = FALSE
  ))

  expect_s3_class(fit, "fastphylosig_result")
  expect_true(all(is.na(fit$requested_simulations)))
  expect_error(
    fast_delta(
      tree, x, test = TRUE, nsim = unused_nsim,
      mcmc_sim = 2L, thin = 1L, burn = 1L,
      verbose = FALSE, progress = FALSE
    ),
    "positive integer when test = TRUE"
  )
})

test_that("trait tables share one species-name contract", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")

  unnamed <- data.frame(trait = c(0, 0, 1, 1))
  fit <- fast_d(
    tree, unnamed, test = FALSE, verbose = FALSE, progress = FALSE
  )
  expect_s3_class(fit, "fastphylosig_result")

  named <- stats::setNames(c(0, 0, 1, 1), tree$tip.label)
  matched <- match_tree_data(tree, data = named, verbose = FALSE)
  expect_named(matched$data, "x")

  wrong_length <- data.frame(trait = c(0, 1, 0))
  expect_error(
    fast_d(
      tree, wrong_length, test = FALSE, verbose = FALSE, progress = FALSE
    ),
    "non-empty species names"
  )

  duplicated <- matrix(
    c("A", "A", "B", "B"), ncol = 1L,
    dimnames = list(c("a", "a", "c", "d"), "trait")
  )
  expect_error(
    fastphylosig:::.as_named_trait_table(
      duplicated, tree, verbose = FALSE, input_name = "data"
    ),
    "species names must be unique"
  )
})
