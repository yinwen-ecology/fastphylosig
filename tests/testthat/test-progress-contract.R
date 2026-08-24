# User-facing progress contract ---------------------------------------------

.progress_contract_capture <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    eval(substitute(expr), envir = parent.frame()),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, messages = messages)
}

.progress_contract_expect <- function(captured, method, stages) {
  prefix <- paste0("^\\[fastphylosig: ", method, "\\] ")
  for (stage in stages) {
    expect_true(
      any(grepl(paste0(prefix, stage), captured$messages)),
      info = paste(method, "missing stage", stage)
    )
  }
  expect_true(
    any(grepl(paste0(prefix, "Done \\| "), captured$messages)),
    info = paste(method, "missing Done message")
  )
  # No C++ symbol, pointer, or per-trait implementation detail is part of the
  # user-facing progress contract.
  forbidden <- paste(
    c("\\.Call", "compiled_tree", "fast_[A-Za-z0-9_]*_cpp",
      "\\btrait(s)?\\b", "\\bchunk(s)?\\b", "\\bnode(s)?\\b"),
    collapse = "|"
  )
  expect_false(any(grepl(forbidden, captured$messages, ignore.case = TRUE)))
  invisible(captured)
}

test_that("successful K/lambda/D/Delta calls expose stable method stages", {
  set.seed(20260812L)
  tree <- ape::rtree(8)
  continuous <- stats::setNames(stats::rnorm(8), tree$tip.label)
  binary <- stats::setNames(rep(0:1, length.out = 8), tree$tip.label)
  categorical <- stats::setNames(rep(c("a", "b", "c"), length.out = 8),
                                  tree$tip.label)

  k <- .progress_contract_capture(fast_k(
    tree, continuous, test = TRUE, nsim = 2, verbose = FALSE,
    progress = TRUE
  ))
  .progress_contract_expect(k, "K", c("Checking tree\\.\\.\\.",
                                      "Preparing data\\.\\.\\.",
                                      "Calculating K\\.\\.\\.",
                                      "Running 2 randomizations\\.\\.\\."))

  lambda <- .progress_contract_capture(fast_lambda(
    tree, continuous, test = TRUE, verbose = FALSE,
    progress = TRUE, lambda_profile = FALSE
  ))
  .progress_contract_expect(lambda, "lambda", c("Checking tree\\.\\.\\.",
                                                "Preparing data\\.\\.\\.",
                                                "Optimizing lambda\\.\\.\\.",
                                                "Testing lambda = 0\\.\\.\\."))

  d <- .progress_contract_capture(fast_d(
    tree, binary, test = TRUE, nsim = 2, return_sim = FALSE,
    verbose = FALSE, progress = TRUE
  ))
  .progress_contract_expect(d, "D", c("Checking tree\\.\\.\\.",
                                      "Preparing data\\.\\.\\.",
                                      "Calculating observed D\\.\\.\\.",
                                      "Simulating random and Brownian nulls\\.\\.\\."))

  delta <- .progress_contract_capture(suppressWarnings(fast_delta(
    tree, categorical, test = TRUE, nsim = 2,
    mcmc_sim = 20, thin = 5, burn = 5,
    verbose = FALSE, progress = TRUE
  )))
  .progress_contract_expect(delta, "Delta", c("Checking tree\\.\\.\\.",
                                              "Preparing data\\.\\.\\.",
                                              "Estimating ancestral states\\.\\.\\.",
                                              "Running MCMC\\.\\.\\.",
                                              "Running permutations\\.\\.\\.",
                                              "Finalizing diagnostics\\.\\.\\."))
})

test_that("progress FALSE is silent and batch calls do not spam per trait", {
  set.seed(20260813L)
  tree <- ape::rtree(8)
  continuous <- stats::setNames(stats::rnorm(8), tree$tip.label)
  binary <- stats::setNames(rep(0:1, length.out = 8), tree$tip.label)
  categorical <- stats::setNames(rep(c("a", "b", "c"), length.out = 8),
                                  tree$tip.label)

  silent_calls <- list(
    .progress_contract_capture(fast_k(
      tree, continuous, verbose = TRUE, progress = FALSE
    )),
    .progress_contract_capture(fast_lambda(
      tree, continuous, verbose = TRUE, progress = FALSE,
      lambda_profile = FALSE
    )),
    .progress_contract_capture(fast_d(
      tree, binary, test = FALSE, nsim = 2, return_sim = FALSE,
      verbose = TRUE, progress = FALSE
    )),
    .progress_contract_capture(suppressWarnings(fast_delta(
      tree, categorical, test = FALSE, mcmc_sim = 20, thin = 5, burn = 5,
      verbose = TRUE, progress = FALSE
    )))
  )
  expect_true(all(vapply(silent_calls, function(z) !length(z$messages),
                         logical(1))))

  X <- matrix(stats::rnorm(8 * 4), nrow = 8,
              dimnames = list(tree$tip.label, paste0("trait_", 1:4)))
  batch <- .progress_contract_capture(fast_k(
    tree, X, verbose = FALSE, progress = TRUE
  ))
  expect_equal(sum(grepl("Checking tree\\.\\.\\.", batch$messages)), 1L)
  expect_equal(sum(grepl("Preparing data\\.\\.\\.", batch$messages)), 1L)
  expect_equal(sum(grepl("Calculating K\\.\\.\\.", batch$messages)), 1L)
  expect_equal(sum(grepl("Done \\|", batch$messages)), 1L)
})

test_that("simulation and likelihood-test stages are conditional", {
  set.seed(20260814L)
  tree <- ape::rtree(8)
  continuous <- stats::setNames(stats::rnorm(8), tree$tip.label)
  binary <- stats::setNames(rep(0:1, length.out = 8), tree$tip.label)
  categorical <- stats::setNames(rep(c("a", "b", "c"), length.out = 8),
                                  tree$tip.label)

  k <- .progress_contract_capture(fast_k(
    tree, continuous, test = FALSE, verbose = FALSE, progress = TRUE
  ))
  expect_false(any(grepl("Running [0-9]+ randomizations", k$messages)))

  lambda <- .progress_contract_capture(fast_lambda(
    tree, continuous, test = FALSE, verbose = FALSE, progress = TRUE,
    lambda_profile = FALSE
  ))
  expect_false(any(grepl("Testing lambda = 0", lambda$messages,
                         fixed = TRUE)))

  d <- .progress_contract_capture(fast_d(
    tree, binary, test = FALSE, nsim = 2, return_sim = FALSE,
    verbose = FALSE, progress = TRUE
  ))
  expect_false(any(grepl("Simulating random and Brownian nulls", d$messages,
                         fixed = TRUE)))

  delta <- .progress_contract_capture(suppressWarnings(fast_delta(
    tree, categorical, test = FALSE, mcmc_sim = 20, thin = 5, burn = 5,
    verbose = FALSE, progress = TRUE
  )))
  expect_false(any(grepl("Running permutations", delta$messages,
                         fixed = TRUE)))
  expect_true(any(grepl("Finalizing diagnostics", delta$messages,
                        fixed = TRUE)))
})
