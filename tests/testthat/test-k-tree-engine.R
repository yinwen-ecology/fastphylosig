test_that("exact tree K engine matches dense GLS on representative trees", {
  skip_if_not(exists("fast_k_tree_batch_cpp", asNamespace("fastphylosig")))
  skip_if_not_installed("ape")

  dense <- function(tree, X) {
    C <- ape::vcv.phylo(tree)
    Q <- solve(C)
    one <- rep(1, nrow(C))
    sum_inv <- as.numeric(crossprod(one, Q %*% one))
    a <- as.numeric(crossprod(one, Q %*% X)) / sum_inv
    centered <- sweep(X, 2L, a, "-")
    numerator <- colSums(centered^2)
    denominator <- colSums(centered * (Q %*% centered))
    normalization <- (sum(diag(C)) - nrow(C) / sum_inv) /
      (nrow(C) - 1)
    list(
      K = (numerator / denominator) / normalization,
      gls_mean = a,
      numerator = numerator,
      denominator = denominator,
      sum_inv = sum_inv,
      normalization = normalization
    )
  }

  compile <- function(tree) {
    d_tree <- ape::reorder.phylo(tree, "pruningwise")
    fastphylosig:::compile_tree_cpp(
      matrix(as.integer(d_tree$edge), ncol = 2L),
      as.numeric(d_tree$edge.length), ape::Ntip(d_tree)
    )
  }

  make_tree <- function(text) {
    tree <- ape::read.tree(text = text)
    tree$edge.length <- as.numeric(tree$edge.length)
    tree
  }

  trees <- list(
    balanced = make_tree(
      "((t1:0.3,t2:0.8):0.6,(t3:0.4,t4:0.9):1.1);"
    ),
    pectinate = make_tree(
      "((((t1:0.2,t2:0.7):0.5,t3:1.2):0.4,t4:0.3):0.8,t5:1.4);"
    ),
    polytomy = make_tree(
      "((t1:0.3,t2:0.4,t3:0.5):0.6,(t4:0.2,t5:0.7):0.8);"
    ),
    two_tip = make_tree("(t1:0.5,t2:1.7);")
  )
  set.seed(20260809)
  random_tree <- ape::rtree(9)
  random_tree$edge.length <- runif(nrow(random_tree$edge), 0.05, 2.5)
  trees$random_nonultrametric <- random_tree

  for (tree in trees) {
    X <- matrix(rnorm(ape::Ntip(tree) * 7L), ncol = 7L)
    rownames(X) <- tree$tip.label
    X <- X[tree$tip.label, , drop = FALSE]
    ref <- dense(tree, X)
    got <- fastphylosig:::fast_k_tree_batch_cpp(
      compile(tree), X, trait_chunk = 2L
    )
    expect_equal(got$K, ref$K, tolerance = 1e-8)
    expect_equal(got$gls_mean, ref$gls_mean, tolerance = 1e-8)
    expect_equal(got$numerator, ref$numerator, tolerance = 1e-8)
    expect_equal(got$denominator, ref$denominator, tolerance = 1e-8)
    expect_equal(got$sum_inv, ref$sum_inv, tolerance = 1e-8)
    expect_equal(got$normalization, ref$normalization, tolerance = 1e-8)

    # Chunking is an implementation detail and must not alter the result.
    got_one <- fastphylosig:::fast_k_tree_batch_cpp(
      compile(tree), X, trait_chunk = 1L
    )
    expect_equal(got_one, got, tolerance = 1e-12)
  }
})

test_that("tree K engine rejects unsupported zero branches and non-finite traits", {
  skip_if_not(exists("fast_k_tree_batch_cpp", asNamespace("fastphylosig")))
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "(a:0.5,b:1.0);")
  compiled <- fastphylosig:::compile_tree_cpp(
    matrix(as.integer(tree$edge), ncol = 2L),
    as.numeric(tree$edge.length), ape::Ntip(tree)
  )
  X <- matrix(c(1, 2), ncol = 1L)
  expect_error(
    fastphylosig:::fast_k_tree_batch_cpp(compiled, X, trait_chunk = 0L),
    "trait_chunk"
  )

  zero <- tree
  zero$edge.length[[1L]] <- 0
  zero_compiled <- fastphylosig:::compile_tree_cpp(
    matrix(as.integer(zero$edge), ncol = 2L),
    as.numeric(zero$edge.length), ape::Ntip(zero)
  )
  expect_error(
    fastphylosig:::fast_k_tree_batch_cpp(zero_compiled, X),
    "strictly positive"
  )
  X[1, 1] <- NA_real_
  expect_error(
    fastphylosig:::fast_k_tree_batch_cpp(compiled, X),
    "finite"
  )
})

test_that("public tree and validate engines preserve dense K", {
  set.seed(20260810)
  tree <- ape::rtree(30)
  X <- matrix(stats::rnorm(30 * 5), 30, 5,
              dimnames = list(tree$tip.label, paste0("x", 1:5)))
  ctx <- prepare_tree(tree)

  tree_fit <- fastphylosig:::.fast_signal_with_engine(
    ctx, X, method = "K", engine = "tree", test = FALSE,
    trait_chunk = 2L, verbose = FALSE
  )
  expect_equal(cache_info(ctx)$n_numerical_entries, 0L)
  dense_fit <- fastphylosig:::.fast_signal_with_engine(
    tree, X, method = "K", engine = "dense", test = FALSE,
    verbose = FALSE
  )
  expect_equal(tree_fit$K_fast, dense_fit$K_fast, tolerance = 1e-8)

  checked <- fastphylosig:::.fast_signal_with_engine(
    tree, X, method = "K", engine = "validate", test = FALSE,
    validate_tolerance = 1e-8, verbose = FALSE
  )
  expect_true(all(checked$validation_pass))
  expect_lt(max(checked$absolute_error), 1e-8)
  expect_equal(checked$K_fast, checked$K_tree, tolerance = 0)
  tested <- fastphylosig:::.fast_signal_with_engine(
    tree, X[, 1L], method = "K", engine = "tree",
    test = TRUE, nsim = 5L, verbose = FALSE
  )
  expect_true(is.finite(tested$K))
  expect_true(is.finite(tested$P))
  expect_equal(tested$nsim_successful, 5L)
})

test_that("tree K residual pass is stable under large common offsets", {
  skip_if_not(exists("fast_k_tree_batch_cpp", asNamespace("fastphylosig")))
  skip_if_not_installed("ape")

  tree <- ape::read.tree(
    text = "((a:0.3,b:0.8):0.6,(c:0.4,d:0.9):1.1);"
  )
  d_tree <- ape::reorder.phylo(tree, "pruningwise")
  compiled <- fastphylosig:::compile_tree_cpp(
    matrix(as.integer(d_tree$edge), ncol = 2L),
    as.numeric(d_tree$edge.length), ape::Ntip(d_tree)
  )

  dense_centered <- function(X) {
    C <- ape::vcv.phylo(d_tree)
    Q <- solve(C)
    one <- rep(1, nrow(C))
    sum_inv <- as.numeric(crossprod(one, Q %*% one))
    baseline <- X[1L, ]
    centered <- sweep(X, 2L, baseline, "-")
    delta <- as.numeric(crossprod(one, Q %*% centered)) / sum_inv
    residual <- sweep(centered, 2L, delta, "-")
    numerator <- colSums(residual^2)
    denominator <- colSums(residual * (Q %*% residual))
    normalization <- (sum(diag(C)) - nrow(C) / sum_inv) /
      (nrow(C) - 1L)
    list(
      K = (numerator / denominator) / normalization,
      gls_mean = baseline + delta,
      numerator = numerator,
      denominator = denominator,
      sum_inv = sum_inv,
      normalization = normalization
    )
  }

  # Differences are deliberately much larger than the ulp at 1e15, so the
  # translated matrix still contains the same observable residuals.
  centered <- matrix(c(
    0, 1200, -800,
    500, -1500, 700,
    -700, 900, 2200,
    1800, -400, -1300
  ), nrow = 4L, byrow = TRUE)
  reference_k <- NULL
  for (offset in c(0, 1e12, 1e15)) {
    X <- sweep(centered, 2L, offset, "+")
    ref <- dense_centered(X)
    got <- fastphylosig:::fast_k_tree_batch_cpp(
      compiled, X, trait_chunk = 1L
    )
    expect_equal(got$K, ref$K, tolerance = 1e-8)
    expect_equal(got$gls_mean, ref$gls_mean, tolerance = 1e-8)
    expect_equal(got$numerator, ref$numerator, tolerance = 1e-8)
    expect_equal(got$denominator, ref$denominator, tolerance = 1e-8)
    expect_equal(got$sum_inv, ref$sum_inv, tolerance = 1e-8)
    expect_equal(got$normalization, ref$normalization, tolerance = 1e-8)
    if (is.null(reference_k)) reference_k <- got$K
    expect_equal(got$K, reference_k, tolerance = 1e-8)
  }

  constant <- matrix(c(rep(1e15, 4L), rep(-1e15, 4L)), nrow = 4L)
  constant_fit <- fastphylosig:::fast_k_tree_batch_cpp(compiled, constant)
  expect_equal(constant_fit$gls_mean, c(1e15, -1e15), tolerance = 0)
  expect_equal(constant_fit$numerator, c(0, 0), tolerance = 0)
  expect_equal(constant_fit$denominator, c(0, 0), tolerance = 0)
  expect_true(all(is.nan(constant_fit$K)))

  near_constant <- sweep(matrix(c(
    0, 0.5,
    1, -0.5,
    -1, 1.5,
    2, -1.5
  ), nrow = 4L, byrow = TRUE), 2L, 1e15, "+")
  near_ref <- dense_centered(near_constant)
  near_fit <- fastphylosig:::fast_k_tree_batch_cpp(
    compiled, near_constant, trait_chunk = 2L
  )
  expect_equal(near_fit$K, near_ref$K, tolerance = 1e-8)
  expect_equal(near_fit$numerator, near_ref$numerator, tolerance = 1e-8)
  expect_equal(near_fit$denominator, near_ref$denominator, tolerance = 1e-8)
})
