test_that("mixed exponent caches and batches retain ACE values and RNG", {
  set.seed(916)
  # Below/at the uncached-work threshold, multiple batches and cache saturation.
  for (case in list(c(2048L, 4L, 2047L), c(2048L, 4L, 2048L),
                    c(8192L, 20L, 12341L), c(8192L, 4L, -1L))) {
    n <- case[1L]; k <- case[2L]; uncached <- case[3L]
    tr <- ape::rtree(n)
    if (uncached < 0L) {
      # More than 4096 repeated groups: uncached repeated values also batch.
      tr$edge.length <- rep(seq_len(n - 1L) / n, each = 2L)
    } else {
      tr$edge.length[seq_len(nrow(tr$edge) - uncached)] <- .2
    }
    w <- .fast_ace_workspace(.fast_ace_tree_workspace(tr), "ARD", k)
    state <- sample.int(k, n, TRUE); state[seq(3, n, 17)] <- NA_integer_
    par <- runif(w$rate$np, .01, .3)
    for (threads in c(1L, 2L, 4L)) {
      ptr <- fast_ace_fit_prepare_cpp(w$edge, w$edge_length, state, w$rate$rate_index, threads)
      for (p in list(par, par * 100, rep(0, length(par)), rep(-1, length(par)), par)) {
        set.seed(381); seed <- .Random.seed
        expected <- fast_ace_discrete_deviance_cpp(w$edge, w$edge_length, state, w$rate$rate_index, p, TRUE)
        expect_identical(fast_ace_fit_deviance_cpp(ptr, p), expected)
        expect_identical(.Random.seed, seed)
      }
      fast_ace_fit_release_cpp(ptr)
    }
  }
})

test_that("mixed-length thread choices retain ARD objective trajectories", {
  set.seed(619); tr <- ape::rtree(4096)
  rows <- sample.int(nrow(tr$edge), nrow(tr$edge) %/% 2L)
  tr$edge.length[rows] <- pmax(.01, round(tr$edge.length[rows], 2))
  state <- sample.int(4L, 4096, TRUE)
  w <- .fast_ace_workspace(.fast_ace_tree_workspace(tr), "ARD", 4L)
  run <- function(threads) {
    set.seed(13)
    ptr <- fast_ace_fit_prepare_cpp(w$edge, w$edge_length, state, w$rate$rate_index, threads)
    on.exit(fast_ace_fit_release_cpp(ptr))
    calls <- list()
    dev <- function(p) {
      value <- fast_ace_fit_deviance_cpp(ptr, p)
      calls[[length(calls) + 1L]] <<- list(par = p[], value = value); value
    }
    fit <- nlminb(rep(.1, 12), dev, lower = rep(0, 12), upper = rep(1e50, 12),
      control = list(iter.max = 15L))
    list(fit = fit, calls = calls, seed = .Random.seed, next_normal = rnorm(5), next_uniform = runif(5))
  }
  expected <- run(1L)
  expect_identical(run(2L), expected)
  expect_identical(run(4L), expected)
})
