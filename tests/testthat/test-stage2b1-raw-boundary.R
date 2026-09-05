# Stage 2B1 raw/prepared boundary contracts

.stage2b1_tree <- function() {
  ape::read.tree(text = paste0(
    "(((a:0.5,b:0.7):0.8,(c:0.4,d:0.6):0.9):0.3,",
    "((e:0.2,f:0.3):0.4,(g:0.5,h:0.6):0.7):0.8);"
  ))
}

.stage2b1_traits <- function(tree) {
  tips <- tree$tip.label
  list(
    continuous = stats::setNames(c(-0.5, -0.2, 0.1, 0.4, 0.7, 1, 1.3, 1.6), tips),
    binary = stats::setNames(c(0L, 0L, 0L, 1L, 1L, 1L, 0L, 1L), tips),
    categorical = stats::setNames(c("a", "b", "c", "a", "b", "c", "a", "b"), tips)
  )
}

.stage2b1_expect_tree_error <- function(tree, label, regexp) {
  traits <- .stage2b1_traits(tree)
  testthat::expect_error(
    fast_k(tree, traits$continuous, test = FALSE, verbose = FALSE,
           progress = FALSE), regexp, info = paste(label, "K")
  )
  testthat::expect_error(
    fast_lambda(tree, traits$continuous, lambda_profile = FALSE,
                verbose = FALSE, progress = FALSE), regexp,
    info = paste(label, "lambda")
  )
  testthat::expect_error(
    fast_d(tree, traits$binary, test = FALSE, return_sim = FALSE,
           verbose = FALSE, progress = FALSE), regexp,
    info = paste(label, "D")
  )
  testthat::expect_error(
    fast_delta(tree, traits$categorical, test = FALSE, mcmc_sim = 8L,
               thin = 2L, burn = 2L, verbose = FALSE, progress = FALSE),
    regexp, info = paste(label, "Delta")
  )
}

test_that("raw and prepared specialist routes remain numerically equivalent", {
  testthat::skip_if_not_installed("ape")
  tree <- .stage2b1_tree()
  traits <- .stage2b1_traits(tree)
  ctx <- prepare_tree(tree)

  raw_k <- suppressMessages(fast_k(
    tree, traits$continuous, test = FALSE, verbose = FALSE, progress = FALSE
  ))
  prepared_k <- suppressMessages(fast_k(
    ctx, traits$continuous, test = FALSE, verbose = FALSE, progress = FALSE
  ))
  testthat::expect_equal(as.numeric(raw_k), as.numeric(prepared_k),
                         tolerance = 1e-12)

  raw_lambda <- suppressMessages(fast_lambda(
    tree, traits$continuous, lambda_profile = FALSE,
    verbose = FALSE, progress = FALSE
  ))
  prepared_lambda <- suppressMessages(fast_lambda(
    ctx, traits$continuous, lambda_profile = FALSE,
    verbose = FALSE, progress = FALSE
  ))
  for (field in c("lambda", "logL", "gls_mean", "sig2")) {
    testthat::expect_equal(raw_lambda[[field]], prepared_lambda[[field]],
                           tolerance = 1e-10, info = field)
  }

  d_nsim <- 4L
  random_states <- matrix(rep(traits$binary, d_nsim), nrow = length(traits$binary))
  brownian_states <- matrix(
    rep(1L - as.integer(traits$binary), d_nsim),
    nrow = length(traits$binary)
  )
  raw_d <- suppressMessages(fast_d(
    tree, traits$binary, test = FALSE, nsim = d_nsim,
    random_states = random_states, brownian_states = brownian_states,
    return_sim = FALSE,
    verbose = FALSE, progress = FALSE
  ))
  prepared_d <- suppressMessages(fast_d(
    ctx, traits$binary, test = FALSE, nsim = d_nsim,
    random_states = random_states, brownian_states = brownian_states,
    return_sim = FALSE,
    verbose = FALSE, progress = FALSE
  ))
  testthat::expect_equal(raw_d$DEstimate, prepared_d$DEstimate,
                         tolerance = 1e-12)
  testthat::expect_equal(raw_d$observed, prepared_d$observed,
                         tolerance = 1e-12)

  set.seed(20260905L)
  raw_delta <- suppressWarnings(fast_delta(
    tree, traits$categorical, test = FALSE, mcmc_sim = 20L,
    thin = 5L, burn = 5L, verbose = FALSE, progress = FALSE
  ))
  set.seed(20260905L)
  prepared_delta <- suppressWarnings(fast_delta(
    ctx, traits$categorical, test = FALSE, mcmc_sim = 20L,
    thin = 5L, burn = 5L, verbose = FALSE, progress = FALSE
  ))
  for (field in c("Delta", "alpha_mean", "beta_mean", "status")) {
    testthat::expect_equal(raw_delta[[field]], prepared_delta[[field]],
                           tolerance = 1e-10, info = field)
  }
})

test_that("raw specialist calls leave tree and trait inputs unchanged", {
  testthat::skip_if_not_installed("ape")
  tree <- .stage2b1_tree()
  traits <- .stage2b1_traits(tree)
  tree_before <- serialize(tree, NULL)
  continuous_before <- serialize(traits$continuous, NULL)
  binary_before <- serialize(traits$binary, NULL)
  categorical_before <- serialize(traits$categorical, NULL)

  suppressMessages(fast_k(tree, traits$continuous, verbose = FALSE,
                          progress = FALSE))
  suppressMessages(fast_lambda(tree, traits$continuous,
                               lambda_profile = FALSE, verbose = FALSE,
                               progress = FALSE))
  suppressMessages(fast_d(tree, traits$binary, test = FALSE,
                          return_sim = FALSE, verbose = FALSE,
                          progress = FALSE))
  suppressWarnings(fast_delta(tree, traits$categorical, test = FALSE,
                              mcmc_sim = 8L, thin = 2L, burn = 2L,
                              verbose = FALSE, progress = FALSE))

  testthat::expect_identical(serialize(tree, NULL), tree_before)
  testthat::expect_identical(serialize(traits$continuous, NULL), continuous_before)
  testthat::expect_identical(serialize(traits$binary, NULL), binary_before)
  testthat::expect_identical(serialize(traits$categorical, NULL), categorical_before)
})

test_that("raw invalid trees retain method-specific failure contracts", {
  testthat::skip_if_not_installed("ape")
  valid <- .stage2b1_tree()

  negative <- valid
  negative$edge.length[[1L]] <- -1
  .stage2b1_expect_tree_error(
    negative, "negative branch", "branch|negative|finite|non-negative|positive|check_tree"
  )

  zero_terminal <- valid
  terminal_row <- which(valid$edge[, 2L] <= ape::Ntip(valid))[[1L]]
  zero_terminal$edge.length[[terminal_row]] <- 0
  .stage2b1_expect_tree_error(
    zero_terminal, "zero terminal branch",
    "branch|zero|positive|USER_ACTION_REQUIRED|check_tree"
  )

  invalid_nnode <- valid
  invalid_nnode$Nnode <- NA_integer_
  .stage2b1_expect_tree_error(
    invalid_nnode, "invalid Nnode", "Nnode|topology|invalid|check_tree|compile"
  )

  malformed_edge <- valid
  malformed_edge$edge <- malformed_edge$edge[, 1L, drop = FALSE]
  .stage2b1_expect_tree_error(
    malformed_edge, "malformed edge", "edge|topology|invalid|check_tree|compile"
  )
})

test_that("raw method-specific root and topology gates remain distinct", {
  testthat::skip_if_not_installed("ape")
  rooted <- .stage2b1_tree()
  polytomy <- ape::read.tree(text = "((a:1,b:1,c:1):1,(d:1,e:1):1);")
  categorical <- stats::setNames(c("a", "b", "c", "a", "b"),
                                 polytomy$tip.label)

  checked <- check_tree(polytomy)
  testthat::expect_true(checked$ready_by_signal[["K"]])
  testthat::expect_true(checked$ready_by_signal[["lambda"]])
  testthat::expect_true(checked$ready_by_signal[["D"]])
  testthat::expect_false(checked$ready_by_signal[["Delta"]])
  testthat::expect_error(
    fast_delta(polytomy, categorical, test = FALSE, mcmc_sim = 8L,
               thin = 2L, burn = 2L, verbose = FALSE, progress = FALSE),
    "binary|polytom|root"
  )

  unrooted <- ape::unroot(rooted)
  testthat::expect_error(
    fast_delta(unrooted,
               stats::setNames(c("a", "b", "c", "a", "b", "c", "a", "b"),
                                unrooted$tip.label),
               test = FALSE, mcmc_sim = 8L, thin = 2L, burn = 2L,
               verbose = FALSE, progress = FALSE),
    "root|Delta|binary"
  )

  unary <- list(
    edge = matrix(c(3L, 4L, 4L, 1L, 3L, 2L), ncol = 2L, byrow = TRUE),
    edge.length = c(1, 1, 1), Nnode = 2L, tip.label = c("a", "b")
  )
  class(unary) <- "phylo"
  unary_binary <- stats::setNames(c(0L, 1L), unary$tip.label)
  unary_categorical <- stats::setNames(c("a", "b"), unary$tip.label)
  testthat::expect_error(
    fast_d(unary, unary_binary, test = FALSE, return_sim = FALSE,
           verbose = FALSE, progress = FALSE), "single-child|unary|D"
  )
  testthat::expect_error(
    fast_delta(unary, unary_categorical, test = FALSE, mcmc_sim = 8L,
               thin = 2L, burn = 2L, verbose = FALSE, progress = FALSE),
    "single-child|unary|Delta|binary"
  )
})
