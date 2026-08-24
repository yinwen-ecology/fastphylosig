# Acceptance tests for the internal Pagel-lambda optimizer.
#
# The optimizer is intentionally an internal milestone while the tree kernel
# is being validated.  Releases that do not expose
# fast_lambda_tree_optimize_cpp() skip these tests cleanly.  The oracle is the
# dense C_lambda GLS likelihood used by phytools::phylosig(method = "lambda")
# followed by stats::optimize() over [0, max_lambda].

.lov_has_engine <- function() {
  exists("fast_lambda_tree_optimize_cpp", envir = asNamespace("fastphylosig"),
         inherits = FALSE)
}

.lov_get <- function(x, candidates) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) {
    hit <- candidates[candidates %in% names(x)]
    if (length(hit)) return(x[[hit[[1L]]]])
  }
  nm <- names(x)
  if (!is.null(nm)) {
    hit <- candidates[candidates %in% nm]
    if (length(hit)) return(unname(x[[hit[[1L]]]]))
  }
  NULL
}

.lov_scalar <- function(x, p, label, default = NULL) {
  if (is.null(x)) return(default)
  x <- as.vector(x)
  if (length(x) == p) return(x)
  if (length(x) == 1L && p == 1L) return(x)
  # A few development wrappers return a p x 1 or 1 x p matrix.
  d <- dim(x)
  if (!is.null(d) && length(x) == p) return(as.vector(x))
  stop(label, " has length ", length(x), "; expected ", p, call. = FALSE)
}

.lov_call <- function(compiled_tree, X, max_lambda = NULL, test = TRUE,
                      profile = FALSE, profile_points = 101L,
                      trait_chunk = 64L, n_threads = 1L) {
  fn <- get("fast_lambda_tree_optimize_cpp", envir = asNamespace("fastphylosig"),
            inherits = FALSE)
  fml <- names(formals(fn))
  args <- list()
  put <- function(candidates, value) {
    hit <- candidates[candidates %in% fml]
    if (length(hit)) args[[hit[[1L]]]] <<- value
  }
  put(c("compiled_tree", "compiled", "tree_core"), compiled_tree)
  put(c("X", "x", "Y", "traits"), X)
  if (!is.null(max_lambda)) {
    put(c("max_lambda", "lambda_max", "upper", "upper_lambda"),
        as.numeric(max_lambda))
  }
  put(c("test", "return_test", "compute_test"), isTRUE(test))
  put(c("profile", "lambda_profile", "return_profile"), isTRUE(profile))
  put(c("profile_points", "lambda_profile_points", "n_profile"),
      as.integer(profile_points))
  put(c("trait_chunk", "chunk_size", "trait_block"), as.integer(trait_chunk))
  put(c("n_threads", "ncores", "threads"), as.integer(n_threads))
  if (!any(names(args) %in% c("compiled_tree", "compiled", "tree_core"))) {
    stop("optimizer has no compiled-tree argument.", call. = FALSE)
  }
  if (!any(names(args) %in% c("X", "x", "Y", "traits"))) {
    stop("optimizer has no trait-matrix argument.", call. = FALSE)
  }
  do.call(fn, args)
}

.lov_unpack <- function(result, p) {
  out <- list(
    lambda = .lov_scalar(.lov_get(result,
                                 c("lambda", "lambda_hat", "lambdas")),
                         p, "lambda"),
    logLik = .lov_scalar(.lov_get(result,
                                  c("logLik", "loglik", "logL", "log_lik")),
                         p, "logLik"),
    gls_mean = .lov_scalar(.lov_get(result,
                                    c("gls_mean", "gls.mean", "mean", "mu")),
                           p, "gls_mean"),
    sigma2 = .lov_scalar(.lov_get(result,
                                  c("sigma2", "sigma_sq", "sig2", "variance")),
                         p, "sigma2"),
    logLik0 = .lov_scalar(.lov_get(result,
                                   c("logLik0", "logL0", "log_lik0",
                                     "loglik0")),
                          p, "logLik0"),
    LR = .lov_scalar(.lov_get(result, c("LR", "lr", "likelihood_ratio")),
                     p, "LR"),
    P = .lov_scalar(.lov_get(result, c("P", "p", "p_value", "pvalue")),
                    p, "P"),
    valid = .lov_scalar(.lov_get(result, c("valid", "is_valid", "ok")),
                        p, "valid"),
    status = .lov_scalar(.lov_get(result,
                                  c("status", "statuses", "message")),
                         p, "status"),
    profile = .lov_get(result, c("profile", "lambda_profile",
                                 "lambda_profile_fast"))
  )
  out
}

.lov_tree <- function(kind, n = 8L) {
  make_stree <- function(type) {
    out <- tryCatch(
      ape::stree(n = n, type = type,
                 tip.label = paste0("sp", seq_len(n))),
      error = function(e) ape::rtree(n)
    )
    out$tip.label <- paste0("sp", seq_len(ape::Ntip(out)))
    out
  }
  tree <- switch(
    kind,
    balanced = make_stree("balanced"),
    pectinate = make_stree("left"),
    unbalanced = {
      out <- make_stree("left")
      out$edge.length <- 10^seq(-4, 2, length.out = nrow(out$edge))
      out
    },
    ultrametric = ape::compute.brlen(make_stree("balanced"), method = "Grafen"),
    nonultrametric = {
      out <- ape::rtree(n)
      out$tip.label <- paste0("sp", seq_len(n))
      out
    },
    very_short = {
      out <- make_stree("balanced")
      out$edge.length <- seq(1e-8, 4e-7, length.out = nrow(out$edge))
      out
    },
    heterogeneous = {
      out <- make_stree("left")
      out$edge.length <- 10^seq(-5, 2, length.out = nrow(out$edge))
      out
    },
    two_tip = ape::read.tree(text = "(sp1:0.4,sp2:1.7);"),
    three_tip = ape::read.tree(text = "((sp1:0.2,sp2:1.1):0.7,sp3:2.3);"),
    polytomy = ape::read.tree(text = "(sp1:0.4,sp2:0.7,sp3:1.2,sp4:1.8);"),
    stop("unknown optimizer tree fixture", call. = FALSE)
  )
  if (is.null(tree$tip.label)) {
    tree$tip.label <- paste0("sp", seq_len(ape::Ntip(tree)))
  }
  tree$tip.label <- paste0("sp", seq_len(ape::Ntip(tree)))
  if (is.null(tree$edge.length) ||
      length(tree$edge.length) != nrow(tree$edge)) {
    tree$edge.length <- rep(1, nrow(tree$edge))
  }
  tree$edge.length <- pmax(as.numeric(tree$edge.length), 1e-12)
  ape::reorder.phylo(tree, "cladewise")
}

.lov_max_lambda <- function(tree) {
  if (exists(".max_lambda", envir = asNamespace("fastphylosig"),
             inherits = FALSE)) {
    z <- tryCatch(fastphylosig:::.max_lambda(tree), error = function(e) NA_real_)
    if (length(z) == 1L && is.finite(z) && z > 0) return(as.numeric(z))
  }
  1
}

.lov_traits <- function(tree, p = 6L, seed = 1L) {
  set.seed(seed)
  n <- ape::Ntip(tree)
  X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p,
              dimnames = list(tree$tip.label,
                              paste0("trait_", seq_len(p))))
  if (p >= 2L) X[, 2L] <- 0.5 + stats::rnorm(n, sd = 1e-7)
  if (p >= 3L) X[, 3L] <- 3
  if (p >= 4L) X[, 4L] <- 1e12 + stats::rnorm(n)
  if (p >= 5L) X[, 5L] <- seq_len(n) / max(1, n)
  if (p >= 6L) X[, 6L] <- -2 + stats::rnorm(n, sd = 1e-4)
  X
}

.lov_compiled <- function(tree) {
  ctx <- fastphylosig::prepare_tree(tree)
  fastphylosig:::.prepared_tree_subset(
    ctx, seq_len(ape::Ntip(tree)), need_lambda = FALSE, need_matrix = TRUE
  )
}

.lov_oracle <- function(tree, y, lambda) {
  C <- ape::vcv.phylo(tree)
  n <- nrow(C)
  Cl <- lambda * C
  diag(Cl) <- diag(C)
  inv <- tryCatch(solve(Cl), error = function(e) NULL)
  if (is.null(inv)) {
    return(list(mean = NA_real_, sigma2 = NA_real_, logLik = -Inf,
                valid = FALSE))
  }
  one <- rep(1, n)
  q1 <- sum(inv)
  if (!is.finite(q1) || q1 <= 0) {
    return(list(mean = NA_real_, sigma2 = NA_real_, logLik = -Inf,
                valid = FALSE))
  }
  baseline <- as.numeric(y[[1L]])
  shifted <- as.numeric(y) - baseline
  delta <- as.numeric(crossprod(one, inv %*% shifted) / q1)
  mu <- baseline + delta
  centered <- shifted - delta
  q <- as.numeric(crossprod(centered, inv %*% centered))
  sig <- q / n
  ld <- tryCatch(as.numeric(determinant(Cl, logarithm = TRUE)$modulus),
                 error = function(e) NA_real_)
  ll <- if (is.finite(sig) && sig > 0 && is.finite(ld)) {
    -0.5 * q / sig - 0.5 * n * log(2 * pi) -
      0.5 * (n * log(sig) + ld)
  } else -Inf
  list(mean = mu, sigma2 = sig, logLik = as.numeric(ll),
       valid = is.finite(ll))
}

.lov_opt_oracle <- function(tree, y, max_lambda = .lov_max_lambda(tree)) {
  lik <- function(z) .lov_oracle(tree, y, z)$logLik
  fit <- stats::optimize(lik, interval = c(0, max_lambda),
                         maximum = TRUE)
  at0 <- lik(0)
  LR <- if (is.finite(fit$objective) && is.finite(at0)) {
    max(0, 2 * (fit$objective - at0))
  } else NA_real_
  P <- if (is.finite(LR)) stats::pchisq(LR, df = 1, lower.tail = FALSE) else NA_real_
  at_hat <- .lov_oracle(tree, y, fit$maximum)
  list(lambda = fit$maximum, logLik = fit$objective, logLik0 = at0,
       LR = LR, P = P, mean = at_hat$mean, sigma2 = at_hat$sigma2)
}

.lov_close <- function(got, ref, abs_tol = 1e-5, rel_tol = 1e-5,
                       label = "value") {
  if (is.finite(ref)) {
    testthat::expect_true(is.finite(got), label = label)
    testthat::expect_lte(abs(got - ref), abs_tol + rel_tol * max(1, abs(ref)),
                         label = label)
  } else if (is.infinite(ref) && ref < 0) {
    testthat::expect_true(is.infinite(got) && got < 0, label = label)
  } else {
    testthat::expect_true(is.na(got), label = label)
  }
}


test_that("lambda optimizer matches dense optimize across tree fixtures", {
  testthat::skip_if_not_installed("ape")
  if (!.lov_has_engine()) {
    testthat::skip("fast_lambda_tree_optimize_cpp is not integrated yet")
  }

  families <- list(
    balanced = .lov_tree("balanced", 8L),
    pectinate = .lov_tree("pectinate", 8L),
    unbalanced = .lov_tree("unbalanced", 8L),
    ultrametric = .lov_tree("ultrametric", 8L),
    nonultrametric = .lov_tree("nonultrametric", 8L),
    very_short = .lov_tree("very_short", 8L),
    heterogeneous = .lov_tree("heterogeneous", 8L),
    two_tip = .lov_tree("two_tip", 2L),
    three_tip = .lov_tree("three_tip", 3L),
    polytomy = .lov_tree("polytomy", 4L)
  )

  for (nm in names(families)) {
    tree <- families[[nm]]
    X <- .lov_traits(tree, p = 2L, seed = 20260840L + nchar(nm))
    # The optimizer is called once for a trait batch.  The near-constant column
    # is retained to exercise stable baseline/residual arithmetic; the dense
    # oracle still determines whether its surface is finite.
    sub <- .lov_compiled(tree)
    Xg <- X[sub$tree$tip.label, , drop = FALSE]
    upper <- .lov_max_lambda(sub$tree)
    got <- .lov_unpack(
      .lov_call(sub$compiled_tree, Xg, max_lambda = upper, test = TRUE,
                profile = FALSE, trait_chunk = 1L, n_threads = 1L),
      ncol(Xg)
    )
    ref <- lapply(seq_len(ncol(Xg)), function(j) {
      .lov_opt_oracle(sub$tree, Xg[, j], upper)
    })
    for (j in seq_len(ncol(Xg))) {
      r <- ref[[j]]
      .lov_close(got$logLik[[j]], r$logLik, label = paste(nm, "logLik", j))
      .lov_close(got$logLik0[[j]], r$logLik0, label = paste(nm, "logLik0", j))
      .lov_close(got$LR[[j]], r$LR, label = paste(nm, "LR", j))
      .lov_close(got$P[[j]], r$P, label = paste(nm, "P", j))
      .lov_close(got$gls_mean[[j]], r$mean, abs_tol = 2e-5,
                 label = paste(nm, "gls_mean", j))
      .lov_close(got$sigma2[[j]], r$sigma2, abs_tol = 2e-5,
                 label = paste(nm, "sigma2", j))
      if (is.finite(r$logLik) && is.finite(got$lambda[[j]])) {
        # A flat likelihood can place the numerical maximum at a different
        # nearby point.  Require lambda agreement on curved surfaces and use
        # the log-likelihood/LR checks above as the acceptance criterion for
        # flat surfaces.
        left <- .lov_oracle(sub$tree, Xg[, j], max(0, r$lambda - 1e-3))$logLik
        right <- .lov_oracle(sub$tree, Xg[, j],
                             min(upper, r$lambda + 1e-3))$logLik
        curvature <- abs(r$logLik - max(left, right))
        if (is.finite(curvature) && curvature > 1e-5) {
          .lov_close(got$lambda[[j]], r$lambda, abs_tol = 1e-5,
                     rel_tol = 1e-5, label = paste(nm, "lambda", j))
        }
      }
    }
    if (!is.null(got$valid)) {
      testthat::expect_true(all(is.logical(got$valid) | is.numeric(got$valid)))
    }
    if (!is.null(got$status)) {
      testthat::expect_true(all(nzchar(as.character(got$status))))
    }
  }
})


test_that("optimizer agrees with phytools and preserves batch/chunk/thread results", {
  testthat::skip_if_not_installed("ape")
  if (!.lov_has_engine()) {
    testthat::skip("fast_lambda_tree_optimize_cpp is not integrated yet")
  }
  tree <- .lov_tree("nonultrametric", 10L)
  X <- .lov_traits(tree, p = 4L, seed = 20260841L)
  sub <- .lov_compiled(tree)
  Xg <- X[sub$tree$tip.label, , drop = FALSE]
  upper <- .lov_max_lambda(sub$tree)
  one <- .lov_unpack(
    .lov_call(sub$compiled_tree, Xg, max_lambda = upper, test = TRUE,
              profile = FALSE, trait_chunk = 1L, n_threads = 1L), ncol(Xg)
  )
  many <- .lov_unpack(
    .lov_call(sub$compiled_tree, Xg, max_lambda = upper, test = TRUE,
              profile = FALSE, trait_chunk = 3L, n_threads = 2L), ncol(Xg)
  )
  for (nm in c("lambda", "logLik", "gls_mean", "sigma2", "logLik0", "LR", "P")) {
    if (!is.null(one[[nm]]) && !is.null(many[[nm]])) {
      testthat::expect_equal(one[[nm]], many[[nm]], tolerance = 2e-8,
                             label = paste("chunk/thread", nm))
    }
  }

  if (requireNamespace("phytools", quietly = TRUE)) {
    for (j in seq_len(ncol(Xg))) {
      y <- Xg[, j]
      phy <- tryCatch(
        suppressWarnings(phytools::phylosig(sub$tree, y, method = "lambda",
                                            test = TRUE, se = FALSE)),
        error = function(e) NULL
      )
      if (is.null(phy) || !is.finite(phy$lambda) || !is.finite(phy$logL)) next
      .lov_close(one$lambda[[j]], phy$lambda, abs_tol = 1e-5,
                 rel_tol = 1e-5, label = paste("phytools lambda", j))
      .lov_close(one$logLik[[j]], phy$logL, abs_tol = 1e-5,
                 rel_tol = 1e-5, label = paste("phytools logL", j))
      if (!is.null(one$logLik0) && is.finite(phy$logL0)) {
        .lov_close(one$logLik0[[j]], phy$logL0, abs_tol = 1e-5,
                   rel_tol = 1e-5, label = paste("phytools logL0", j))
      }
      if (!is.null(one$LR) && is.finite(phy$logL0)) {
        expected_lr <- max(0, 2 * (phy$logL - phy$logL0))
        .lov_close(one$LR[[j]], expected_lr, abs_tol = 1e-5,
                   rel_tol = 1e-5, label = paste("phytools LR", j))
      }
      if (!is.null(one$P) && is.finite(phy$P)) {
        .lov_close(one$P[[j]], phy$P, abs_tol = 1e-5,
                   rel_tol = 1e-5, label = paste("phytools P", j))
      }
    }
  }
})


test_that("optimizer reports boundary and degenerate traits explicitly", {
  testthat::skip_if_not_installed("ape")
  if (!.lov_has_engine()) {
    testthat::skip("fast_lambda_tree_optimize_cpp is not integrated yet")
  }
  tree <- .lov_tree("balanced", 6L)
  n <- ape::Ntip(tree)
  X <- cbind(
    constant = rep(7, n),
    near_constant = 7 + stats::rnorm(n, sd = 1e-9),
    random = stats::rnorm(n),
    offset = 1e12 + stats::rnorm(n)
  )
  rownames(X) <- tree$tip.label
  sub <- .lov_compiled(tree)
  upper <- .lov_max_lambda(sub$tree)
  out <- .lov_unpack(
    .lov_call(sub$compiled_tree, X[sub$tree$tip.label, , drop = FALSE],
              max_lambda = upper, test = TRUE, profile = FALSE,
              trait_chunk = 2L, n_threads = 1L), ncol(X)
  )
  if (!is.null(out$valid)) {
    testthat::expect_false(isTRUE(out$valid[[1L]]))
    testthat::expect_true(isTRUE(out$valid[[2L]]) || isTRUE(out$valid[[3L]]))
  }
  if (!is.null(out$status)) {
    testthat::expect_true(all(nzchar(as.character(out$status))))
  }
  for (j in 2:4) {
    if (!is.null(out$lambda) && is.finite(out$lambda[[j]])) {
      testthat::expect_gte(out$lambda[[j]], -1e-10)
      testthat::expect_lte(out$lambda[[j]], upper + 1e-8)
    }
  }
})


test_that("optional optimizer profile contains endpoints, optimum, and CI surface", {
  testthat::skip_if_not_installed("ape")
  if (!.lov_has_engine()) {
    testthat::skip("fast_lambda_tree_optimize_cpp is not integrated yet")
  }
  tree <- .lov_tree("balanced", 8L)
  X <- .lov_traits(tree, p = 1L, seed = 20260842L)
  sub <- .lov_compiled(tree)
  upper <- .lov_max_lambda(sub$tree)
  raw <- .lov_call(sub$compiled_tree, X[sub$tree$tip.label, , drop = FALSE],
                   max_lambda = upper, test = TRUE, profile = TRUE,
                   profile_points = 31L, trait_chunk = 1L, n_threads = 1L)
  out <- .lov_unpack(raw, 1L)
  prof <- out$profile
  if (is.null(prof)) {
    testthat::skip("optimizer profile is optional and not returned by this build")
  }
  if (is.list(prof) && !is.data.frame(prof) && length(prof) == 1L) {
    prof <- prof[[1L]]
  }
  if (is.list(prof) && !is.data.frame(prof) &&
      all(c("lambda", "logLik") %in% names(prof))) {
    prof <- data.frame(lambda = prof$lambda, logLik = prof$logLik)
  }
  testthat::expect_true(is.data.frame(prof) || is.matrix(prof))
  if (is.matrix(prof)) prof <- as.data.frame(prof)
  lam <- .lov_get(prof, c("lambda", "lambda_grid", "x"))
  ll <- .lov_get(prof, c("logLik", "logL", "loglik", "y"))
  testthat::expect_true(length(lam) >= 3L && length(ll) == length(lam))
  testthat::expect_true(any(abs(as.numeric(lam)) <= 1e-10))
  testthat::expect_true(any(abs(as.numeric(lam) - upper) <=
                             1e-8 * max(1, upper)))
  testthat::expect_true(any(is.finite(as.numeric(ll))))
  if (is.finite(out$lambda[[1L]])) {
    testthat::expect_true(any(abs(as.numeric(lam) - out$lambda[[1L]]) <=
                               1e-8 * max(1, upper)))
  }
  # CI fields, when present, must be ordered and lie on the profiled interval.
  ci <- .lov_get(raw, c("CI", "lambda_CI", "ci"))
  if (!is.null(ci) && length(ci) >= 2L && all(is.finite(as.numeric(ci)))) {
    testthat::expect_lte(as.numeric(ci[[1L]]), as.numeric(ci[[2L]]))
    testthat::expect_gte(as.numeric(ci[[1L]]), min(as.numeric(lam)) - 1e-8)
    testthat::expect_lte(as.numeric(ci[[2L]]), max(as.numeric(lam)) + 1e-8)
  }
})
