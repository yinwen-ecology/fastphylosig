test_that("progress output is method-level and reports elapsed runtime", {
  set.seed(811)
  tree <- ape::rtree(8)
  continuous <- stats::setNames(stats::rnorm(8), tree$tip.label)
  binary <- stats::setNames(rep(0:1, length.out = 8), tree$tip.label)
  categorical <- stats::setNames(rep(letters[1:3], length.out = 8),
                                  tree$tip.label)

  expect_message(
    fast_signal(tree, continuous, method = "K", progress = TRUE,
                verbose = TRUE),
    "\\[fastphylosig: K\\] Preparing data\\.\\.\\."
  )
  expect_message(
    fast_d(tree, binary, nsim = 2, return_sim = FALSE,
           progress = TRUE, verbose = TRUE),
    "\\[fastphylosig: D\\] Done \\| [0-9.]+ s"
  )
  expect_message(
    suppressWarnings(fast_delta(
      tree, categorical, mcmc_sim = 20, thin = 5, burn = 5,
      progress = TRUE, verbose = TRUE
    )),
    "\\[fastphylosig: Delta\\] Done \\| [0-9.]+ s"
  )
  expect_no_warning(
    expect_message(
      fast_ace(categorical, tree, CI = FALSE, progress = TRUE),
      "\\[fastphylosig: ACE\\] Done \\| [0-9.]+ s"
    )
  )
})

test_that("progress FALSE is fully silent and timing is stable", {
  set.seed(812)
  tree <- ape::rtree(7)
  x <- stats::setNames(stats::rnorm(7), tree$tip.label)

  scalar <- expect_silent(fast_signal(
    tree, x, method = "K", progress = FALSE, verbose = TRUE
  ))
  timing_scalar <- attr(scalar, "timing")
  expect_true(is.list(timing_scalar))
  expect_true(is.finite(timing_scalar$total_elapsed))
  expect_gte(timing_scalar$total_elapsed, 0)

  lambda <- expect_silent(fast_signal(
    tree, x, method = "lambda", progress = FALSE, verbose = TRUE
  ))
  expect_true(is.list(lambda$timing))
  expect_identical(attr(lambda, "timing"), lambda$timing)

  X <- cbind(a = x, b = x + 1)
  batch <- expect_silent(fast_signal(
    tree, X, method = "K", progress = FALSE, verbose = TRUE
  ))
  expect_true(is.list(attr(batch, "timing")))
  expect_true(is.finite(attr(batch, "timing")$total_elapsed))
})

test_that("batch progress has no per-trait stage spam", {
  set.seed(813)
  tree <- ape::rtree(8)
  X <- matrix(stats::rnorm(8 * 4), nrow = 8,
              dimnames = list(tree$tip.label, paste0("trait_", 1:4)))
  messages <- character()
  withCallingHandlers(
    fast_signal(tree, X, method = "K", progress = TRUE, verbose = TRUE),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_equal(sum(grepl("Calculating K...", messages, fixed = TRUE)), 1L)
  expect_equal(sum(grepl("Done \\|", messages)), 1L)
})

test_that("elapsed time formatting is readable for long jobs", {
  expect_identical(fastphylosig:::.runtime_format_elapsed(0.423), "0.42 s")
  expect_identical(fastphylosig:::.runtime_format_elapsed(134.6),
                   "2 min 14.6 s")
  expect_identical(fastphylosig:::.runtime_format_elapsed(3723.4),
                   "1 h 2 min 3.4 s")
})
