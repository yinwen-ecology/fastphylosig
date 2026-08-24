test_that("tree lambda optimizer matches R bounded optimize and fixed kernel", {
  skip_if_not_installed("ape")
  skip_if_not(exists("fast_lambda_tree_fixed_cpp", asNamespace("fastphylosig")))
  skip_if_not(exists("fast_lambda_tree_optimize_cpp", asNamespace("fastphylosig")))

  tree <- ape::read.tree(
    text = "((a:0.3,b:0.8):0.6,(c:0.4,d:0.9):1.1);"
  )
  compiled <- fastphylosig:::compile_tree_cpp(
    matrix(as.integer(tree$edge), ncol = 2L),
    as.numeric(tree$edge.length), ape::Ntip(tree)
  )
  set.seed(20260811)
  X <- matrix(stats::rnorm(ape::Ntip(tree) * 3L), ncol = 3L,
              dimnames = list(tree$tip.label, paste0("x", 1:3)))
  X <- X[tree$tip.label, , drop = FALSE]

  got <- fastphylosig:::fast_lambda_tree_optimize_cpp(
    compiled, X, test = TRUE, profile = TRUE, profile_points = 31L,
    trait_chunk = 2L, tol = 1e-8
  )
  expect_true(all(got$valid))
  expect_equal(length(got$lambda), ncol(X))
  expect_true(all(is.finite(got$lambda)))

  for (j in seq_len(ncol(X))) {
    objective <- function(z) {
      fixed <- fastphylosig:::fast_lambda_tree_fixed_cpp(
        compiled, X[, j, drop = FALSE], z, trait_chunk = 1L
      )
      as.numeric(fixed$logLik[1, 1])
    }
    ref <- stats::optimize(objective, c(0, got$max_lambda), maximum = TRUE,
                           tol = 1e-8)
    expect_equal(got$lambda[j], ref$maximum, tolerance = 1e-8)
    expect_equal(got$logLik[j], ref$objective, tolerance = 1e-8)

    fixed_hat <- fastphylosig:::fast_lambda_tree_fixed_cpp(
      compiled, X[, j, drop = FALSE], got$lambda[j], trait_chunk = 1L
    )
    expect_equal(got$gls_mean[j], fixed_hat$gls_mean[1, 1], tolerance = 1e-8)
    expect_equal(got$sigma2[j], fixed_hat$sigma2[1, 1], tolerance = 1e-8)
    expect_equal(got$logLik0[j], objective(0), tolerance = 1e-8)
    expect_equal(got$LR[j], max(0, 2 * (got$logLik[j] - got$logLik0[j])),
                 tolerance = 1e-8)
    expect_equal(got$P[j], stats::pchisq(got$LR[j], df = 1,
                                         lower.tail = FALSE),
                 tolerance = 1e-8)
  }

  prof <- got$profile
  expect_true(is.list(prof))
  expect_true(any(abs(prof$lambda) < 1e-15))
  expect_true(any(abs(prof$lambda - got$max_lambda) < 1e-15))
  if (got$max_lambda >= 1) expect_true(any(abs(prof$lambda - 1) < 1e-15))
  for (j in seq_len(ncol(X))) {
    expect_true(any(abs(prof$lambda - got$lambda[j]) < 1e-15))
  }
  expect_equal(dim(prof$logLik), c(length(prof$lambda), ncol(X)))

  if (requireNamespace("phytools", quietly = TRUE)) {
    got_phytools <- fastphylosig:::fast_lambda_tree_optimize_cpp(
      compiled, X[, 1L, drop = FALSE], test = FALSE, profile = FALSE
    )
    ref_phytools <- phytools::phylosig(
      tree, X[, 1L], method = "lambda", test = FALSE, se = NULL
    )
    expect_equal(got_phytools$lambda[[1L]], ref_phytools$lambda,
                 tolerance = 5e-4)
    expect_equal(got_phytools$logLik[[1L]], ref_phytools$logL,
                 tolerance = 1e-5)
  }
})

test_that("optimizer handles ultrametric boundaries, two tips, and batches", {
  skip_if_not_installed("ape")
  skip_if_not(exists("fast_lambda_tree_fixed_cpp", asNamespace("fastphylosig")))
  skip_if_not(exists("fast_lambda_tree_optimize_cpp", asNamespace("fastphylosig")))

  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  compiled <- fastphylosig:::compile_tree_cpp(
    matrix(as.integer(tree$edge), ncol = 2L),
    as.numeric(tree$edge.length), ape::Ntip(tree)
  )
  X <- cbind(c(0, 1, 2, 4), c(5, 5.1, 4.8, 5.3))
  got <- fastphylosig:::fast_lambda_tree_optimize_cpp(
    compiled, X, max_lambda = 1.5, test = TRUE, profile = TRUE,
    profile_points = 17L, trait_chunk = 1L
  )
  expect_true(all(got$valid))
  expect_lte(max(got$lambda), 1.5 + 1e-12)
  expect_true(any(abs(got$profile$lambda - 1) < 1e-15))
  expect_true(any(abs(got$profile$lambda - got$max_lambda) < 1e-15))
  expect_true(all(is.finite(got$logLik0)))
  got_threads <- fastphylosig:::fast_lambda_tree_optimize_cpp(
    compiled, X, max_lambda = 1.5, test = TRUE, profile = FALSE,
    trait_chunk = 2L, n_threads = 3L
  )
  expect_equal(got_threads$lambda, got$lambda, tolerance = 1e-12)
  expect_equal(got_threads$logLik, got$logLik, tolerance = 1e-12)

  # At the full phytools boundary the transformed terminal edge is singular;
  # it must remain a profile -Inf point and must not become lambda_hat.
  full <- fastphylosig:::fast_lambda_tree_optimize_cpp(
    compiled, X, test = TRUE, profile = TRUE, profile_points = 17L,
    trait_chunk = 2L, tol = 1e-8
  )
  expect_true(all(full$valid))
  boundary_row <- which.min(abs(full$profile$lambda - full$max_lambda))
  expect_true(all(full$profile$logLik[boundary_row, ] == -Inf))
  expect_true(all(full$lambda < full$max_lambda))

  two <- ape::read.tree(text = "(a:0.5,b:1.7);")
  two_compiled <- fastphylosig:::compile_tree_cpp(
    matrix(as.integer(two$edge), ncol = 2L),
    as.numeric(two$edge.length), ape::Ntip(two)
  )
  y <- matrix(c(1, 4), ncol = 1L)
  two_fit <- fastphylosig:::fast_lambda_tree_optimize_cpp(
    two_compiled, y, test = TRUE, profile = FALSE
  )
  expect_true(two_fit$valid[[1]])
  expect_true(is.finite(two_fit$lambda[[1]]))
  expect_true(is.finite(two_fit$logLik[[1]]))
})

test_that("optimizer diagnoses constant traits and invalid optimization bounds", {
  skip_if_not_installed("ape")
  skip_if_not(exists("fast_lambda_tree_optimize_cpp", asNamespace("fastphylosig")))
  tree <- ape::read.tree(text = "((a:0.5,b:0.9):0.7,(c:0.4,d:1.2):0.8);")
  compiled <- fastphylosig:::compile_tree_cpp(
    matrix(as.integer(tree$edge), ncol = 2L),
    as.numeric(tree$edge.length), ape::Ntip(tree)
  )
  constant <- matrix(rep(1e15, ape::Ntip(tree)), ncol = 1L)
  bad <- fastphylosig:::fast_lambda_tree_optimize_cpp(
    compiled, constant, test = TRUE, profile = TRUE, profile_points = 7L
  )
  expect_false(bad$valid[[1]])
  expect_true(is.na(bad$lambda[[1]]))
  expect_true(is.na(bad$gls_mean[[1]]) || is.finite(bad$gls_mean[[1]]))
  expect_true(is.na(bad$logLik0[[1]]) || bad$logLik0[[1]] == -Inf)
  expect_true(nzchar(bad$status[[1]]))

  expect_error(
    fastphylosig:::fast_lambda_tree_optimize_cpp(
      compiled, matrix(c(1, 2, 3, 4), ncol = 1L), max_lambda = 0,
      test = TRUE
    ),
    "max_lambda"
  )
  expect_error(
    fastphylosig:::fast_lambda_tree_optimize_cpp(
      compiled, matrix(c(1, 2, 3, 4), ncol = 1L), profile = TRUE,
      profile_points = 1L
    ),
    "profile_points"
  )
})
