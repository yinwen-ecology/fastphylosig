test_that("fixed lambda tree kernel matches dense C_lambda GLS oracle", {
  skip_if_not_installed("ape")
  skip_if_not(exists("fast_lambda_tree_fixed_cpp", asNamespace("fastphylosig")))

  dense <- function(tree, X, lambda) {
    C <- ape::vcv.phylo(tree)
    n <- nrow(C)
    p <- ncol(X)
    out <- list(
      gls_mean = matrix(NA_real_, length(lambda), p),
      sigma2 = matrix(NA_real_, length(lambda), p),
      logLik = matrix(-Inf, length(lambda), p),
      valid = matrix(FALSE, length(lambda), p)
    )
    for (i in seq_along(lambda)) {
      z <- lambda[[i]]
      if (!is.finite(z) || z < 0) next
      Cl <- z * (C - diag(diag(C))) + diag(diag(C))
      Q <- tryCatch(solve(Cl), error = function(e) NULL)
      if (is.null(Q)) next
      one <- rep(1, n)
      q1 <- drop(crossprod(one, Q %*% one))
      if (!is.finite(q1) || q1 <= 0) next
      for (j in seq_len(p)) {
        y <- X[, j]
        mu <- drop(crossprod(one, Q %*% y)) / q1
        r <- y - mu
        q <- drop(crossprod(r, Q %*% r))
        sig <- q / n
        if (!is.finite(sig) || sig <= 0) next
        ld <- determinant(Cl, logarithm = TRUE)$modulus[[1L]]
        ll <- -0.5 * q / sig - 0.5 * n * log(2 * pi) -
          0.5 * (n * log(sig) + ld)
        out$gls_mean[i, j] <- mu
        out$sigma2[i, j] <- sig
        out$logLik[i, j] <- ll
        out$valid[i, j] <- is.finite(ll)
      }
    }
    out
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
    balanced = make_tree("((a:0.3,b:0.8):0.6,(c:0.4,d:0.9):1.1);"),
    pectinate = make_tree(
      "((((a:0.2,b:0.7):0.5,c:1.2):0.4,d:0.3):0.8,e:1.4);"
    ),
    polytomy = make_tree("((a:0.3,b:0.4,c:0.5):0.6,(d:0.2,e:0.7):0.8);"),
    two_tip = make_tree("(a:0.5,b:1.7);")
  )
  set.seed(20260809)
  random_tree <- ape::rtree(9)
  random_tree$edge.length <- stats::runif(nrow(random_tree$edge), 0.05, 2.5)
  trees$random_nonultrametric <- random_tree

  for (tree in trees) {
    n <- ape::Ntip(tree)
    X <- matrix(stats::rnorm(n * 4L), ncol = 4L,
                dimnames = list(tree$tip.label, paste0("x", seq_len(4L))))
    X <- X[tree$tip.label, , drop = FALSE]
    lambda <- c(0, 1e-12, 0.2, 0.5, 1)
    ref <- dense(tree, X, lambda)
    got <- fastphylosig:::fast_lambda_tree_fixed_cpp(
      compile(tree), X, lambda, trait_chunk = 2L
    )
    expect_equal(got$gls_mean[got$valid], ref$gls_mean[ref$valid],
                 tolerance = 1e-8)
    expect_equal(got$sigma2[got$valid], ref$sigma2[ref$valid],
                 tolerance = 1e-8)
    expect_equal(got$logLik[got$valid], ref$logLik[ref$valid],
                 tolerance = 1e-8)
    expect_true(all(got$valid == ref$valid))

    # Trait chunking must not alter any fixed-surface value.
    got_one <- fastphylosig:::fast_lambda_tree_fixed_cpp(
      compile(tree), X, lambda, trait_chunk = 1L
    )
    expect_equal(got_one$gls_mean, got$gls_mean, tolerance = 1e-12)
    expect_equal(got_one$sigma2, got$sigma2, tolerance = 1e-12)
    expect_equal(got_one$logLik, got$logLik, tolerance = 1e-12)
  }
})

test_that("lambda zero, max boundary, and invalid rows are diagnosable", {
  skip_if_not_installed("ape")
  skip_if_not(exists("fast_lambda_tree_fixed_cpp", asNamespace("fastphylosig")))

  # Ultrametric trees can have a phytools upper bound above one.  The exact
  # boundary has a zero transformed terminal variance and is singular.
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  compiled <- fastphylosig:::compile_tree_cpp(
    matrix(as.integer(tree$edge), ncol = 2L),
    as.numeric(tree$edge.length), ape::Ntip(tree)
  )
  X <- matrix(c(0, 1, 2, 4), ncol = 1L)
  max_lam <- fastphylosig:::fast_lambda_tree_fixed_cpp(
    compiled, X, 0.5
  )$max_lambda
  expect_gt(max_lam, 1)
  got <- fastphylosig:::fast_lambda_tree_fixed_cpp(
    compiled, X, c(-1, 0, max_lam * (1 - 1e-8), max_lam, max_lam + 1),
    trait_chunk = 1L
  )
  expect_true(got$valid[2, 1])
  expect_true(got$valid[3, 1])
  expect_false(got$valid[1, 1])
  expect_false(got$valid[4, 1])
  expect_false(got$valid[5, 1])
  expect_equal(got$logLik[!got$valid], rep(-Inf, sum(!got$valid)))
  expect_true(all(nzchar(got$status[!got$valid])))

  # Constant traits have a valid GLS mean but no positive ML residual scale.
  constant <- matrix(rep(1e15, ape::Ntip(tree)), ncol = 1L)
  const_fit <- fastphylosig:::fast_lambda_tree_fixed_cpp(
    compiled, constant, c(0, 0.5, 1), trait_chunk = 1L
  )
  expect_true(all(!const_fit$valid[, 1]))
  expect_true(all(is.na(const_fit$sigma2[, 1])))
  expect_true(all(const_fit$logLik[, 1] == -Inf))
  expect_equal(const_fit$gls_mean[, 1], rep(1e15, 3L), tolerance = 0)
})

test_that("fixed lambda kernel rejects malformed dimensions and non-finite traits", {
  skip_if_not_installed("ape")
  skip_if_not(exists("fast_lambda_tree_fixed_cpp", asNamespace("fastphylosig")))
  tree <- ape::read.tree(text = "(a:0.5,b:1.0);")
  compiled <- fastphylosig:::compile_tree_cpp(
    matrix(as.integer(tree$edge), ncol = 2L),
    as.numeric(tree$edge.length), ape::Ntip(tree)
  )
  expect_error(
    fastphylosig:::fast_lambda_tree_fixed_cpp(compiled, matrix(1:3, ncol = 1), 0.5),
    "one row per"
  )
  bad <- matrix(c(1, NA_real_), ncol = 1L)
  expect_error(
    fastphylosig:::fast_lambda_tree_fixed_cpp(compiled, bad, 0.5),
    "finite"
  )
  expect_error(
    fastphylosig:::fast_lambda_tree_fixed_cpp(compiled, matrix(c(1, 2), ncol = 1),
                                              0.5, trait_chunk = 0L),
    "trait_chunk"
  )
})

