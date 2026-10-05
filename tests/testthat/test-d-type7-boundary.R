# Independent reference: R's type-7 quantile and the documented strict cutoff.
d_type7_states <- function(raw, p) {
  matrix(vapply(seq_len(ncol(raw)), function(j) {
    as.integer(raw[, j] < stats::quantile(raw[, j], p, type = 7,
                                         names = FALSE))
  }, integer(nrow(raw))), nrow = nrow(raw), ncol = ncol(raw))
}

# Reproduce branch draws only; thresholding is entirely in R.
d_type7_raw <- function(tree, nsim) {
  n <- length(tree$tip.label)
  children <- split(seq_len(nrow(tree$edge)), tree$edge[, 1])
  vapply(seq_len(nsim), function(i) {
    values <- numeric(max(tree$edge))
    stack <- n + 1L
    while (length(stack)) {
      parent <- tail(stack, 1L)
      stack <- head(stack, -1L)
      for (e in children[[as.character(parent)]]) {
        child <- tree$edge[e, 2]
        values[child] <- values[parent] + sqrt(tree$edge.length[e]) * rnorm(1)
        if (child > n) stack <- c(stack, child)
      }
    }
    values[seq_len(n)]
  }, numeric(n))
}

test_that('D matrix thresholds agree independently with R type-7 quantiles', {
  for (v in list(5, c(-2, 3), 0:3, c(-7, -1, 0, 2, 9, 12),
                 c(0, 0, 1, 1), c(-3, -3, -3))) {
    raw <- cbind(v, rev(v))
    for (p in unique(c(0, 1, 1 / length(v), .01, .1, .25, .5, .73))) {
      expect_equal(fastphylosig:::brownian_threshold_cpp(raw, p),
                   d_type7_states(raw, p), info = paste(length(v), p))
    }
  }
  expect_equal(as.vector(fastphylosig:::brownian_threshold_cpp(
    matrix(0:3, 4L), .25)), c(1L, 0L, 0L, 0L))
})

test_that('D tree thresholding preserves prevalence, branch draws and strict ties', {
  tree <- ape::reorder.phylo(ape::read.tree(text =
    '((a:0.2,b:1.1):0.4,((c:0.7,d:1.8):0.3,(e:0.6,f:1.3):0.9):0.5);'),
    'pruningwise')
  n <- length(tree$tip.label)
  for (p in c(0, 1/n, .5, 1 - 1/n, 1)) {
    set.seed(811)
    raw <- d_type7_raw(tree, 64L)
    reference_seed <- .Random.seed
    expected <- d_type7_states(raw, p)
    set.seed(811)
    actual <- fastphylosig:::brownian_tree_threshold_cpp(
      tree$edge, tree$edge.length, n, 64L, p)
    expect_equal(actual, expected)
    expect_identical(.Random.seed, reference_seed)
    if (p == 1/n) expect_equal(colSums(actual), rep(1, 64))
    if (p == 1 - 1/n) expect_equal(colSums(actual), rep(n - 1, 64))
  }
})

test_that('rare-state fast_d streaming and compatibility use independently calibrated nulls', {
  tree <- ape::reorder.phylo(ape::read.tree(text =
    '((a:0.2,b:1.1):0.4,((c:0.7,d:1.8):0.3,(e:0.6,f:1.3):0.9):0.5);'),
    'pruningwise')
  n <- length(tree$tip.label)
  nsim <- 64L
  for (rare_first in c(TRUE, FALSE)) {
    x <- setNames(factor(c('rare', rep('common', n - 1)),
      levels = if (rare_first) c('rare', 'common') else c('common', 'rare')),
      tree$tip.label)
    ds <- fastphylosig:::.binary_state(x)$values
    p <- if (rare_first) 1/n else 1 - 1/n
    set.seed(912)
    # Streaming draws all random permutations in each chunk before Brownian.
    raw <- do.call(cbind, lapply(1:2, function(chunk) {
      for (sim in 1:32) for (tip in (n - 1):1) runif(1, 0, tip + 1)
      d_type7_raw(tree, 32L)
    }))
    reference_seed <- .Random.seed
    brownian <- d_type7_states(raw, p)
    expect_equal(colSums(brownian), rep(if (rare_first) 1 else n - 1, nsim))
    expected <- fastphylosig:::phylo_d_sums_cpp(
      brownian, tree$edge, tree$edge.length, n)
    set.seed(912)
    fit <- fastphylosig::fast_d(tree, x, nsim = nsim, chunk_size = 32L,
                               verbose = FALSE, progress = FALSE)
    expect_identical(.Random.seed, reference_seed)
    expect_equal(fit$Permutations$brownian, as.numeric(expected), tolerance = 1e-12)
    expect_gt(fit$Parameters$MeanBrownian, 0)
    expect_true(is.finite(fit$DEstimate))
    expect_equal(fit$DEstimate,
      (fit$Parameters$Observed - mean(expected)) /
        (fit$Parameters$MeanRandom - mean(expected)), tolerance = 1e-12)
    expect_equal(fit$P_Brownian, mean(expected > fit$Parameters$Observed))
    # Force compatibility generation with supplied random states.
    random <- matrix(rep(ds, nsim), nrow = n)
    set.seed(913)
    raw <- d_type7_raw(tree, nsim)
    reference_seed <- .Random.seed
    brownian <- d_type7_states(raw, p)
    set.seed(913)
    compat <- fastphylosig::fast_d(tree, x, nsim = nsim, random_states = random,
                                  verbose = FALSE, progress = FALSE)
    expect_identical(.Random.seed, reference_seed)
    controlled <- fastphylosig::fast_d(tree, x, nsim = nsim,
      random_states = random, brownian_states = brownian,
      verbose = FALSE, progress = FALSE)
    expect_equal(compat$Parameters, controlled$Parameters, tolerance = 1e-12)
    expect_equal(compat$Permutations, controlled$Permutations, tolerance = 1e-12)
    expect_equal(compat$DEstimate, controlled$DEstimate, tolerance = 1e-12)
  }
})

test_that('two-tip D retains natural calibration degeneracy', {
  tree <- ape::read.tree(text = '(a:0.4,b:1.3);')
  x <- setNames(c(0, 1), tree$tip.label)
  set.seed(914)
  fit <- fastphylosig::fast_d(tree, x, nsim = 32L,
                             verbose = FALSE, progress = FALSE)
  expect_equal(fit$Parameters$MeanRandom, fit$Parameters$MeanBrownian)
  expect_equal(fit$Parameters$MeanBrownian, 1)
  expect_true(is.na(fit$DEstimate))
  expect_match(fit$note, 'identical')
})
