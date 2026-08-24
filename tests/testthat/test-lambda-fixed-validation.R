# Acceptance tests for the tree-linear fixed-lambda likelihood milestone.
#
# The production entry point is deliberately internal until the likelihood
# surface has been validated.  These tests therefore skip cleanly on releases
# that do not yet expose fast_lambda_tree_fixed_cpp().  Once present, every
# finite cell is compared with the dense Pagel covariance operation order used
# by phytools::phylosig(): C_lambda = lambda * offdiag(C) + diag(C).

.lfv_has_engine <- function() {
  exists("fast_lambda_tree_fixed_cpp", envir = asNamespace("fastphylosig"),
         inherits = FALSE)
}

.lfv_get <- function(x, candidates) {
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

.lfv_call <- function(compiled_tree, X, lambda, trait_chunk = 4L,
                      n_threads = 1L) {
  fn <- get("fast_lambda_tree_fixed_cpp", envir = asNamespace("fastphylosig"),
            inherits = FALSE)
  fml <- names(formals(fn))
  args <- list()
  put <- function(candidates, value) {
    hit <- candidates[candidates %in% fml]
    if (length(hit)) args[[hit[[1L]]]] <<- value
  }
  # The milestone contract uses these names.  Aliases keep the acceptance
  # test useful while an Rcpp wrapper is being regenerated during development.
  put(c("compiled_tree", "compiled", "tree_core"), compiled_tree)
  put(c("X", "x", "Y", "traits"), X)
  put(c("lambda", "lambdas", "lambda_grid"), as.numeric(lambda))
  put(c("trait_chunk", "chunk_size"), as.integer(trait_chunk))
  put(c("n_threads", "ncores", "threads"), as.integer(n_threads))
  if (!length(args) || !any(names(args) %in% c(
    "compiled_tree", "compiled", "tree_core"))) {
    stop("fast_lambda_tree_fixed_cpp has no compiled-tree argument.",
         call. = FALSE)
  }
  if (!any(names(args) %in% c("X", "x", "Y", "traits"))) {
    stop("fast_lambda_tree_fixed_cpp has no trait-matrix argument.",
         call. = FALSE)
  }
  if (!any(names(args) %in% c("lambda", "lambdas", "lambda_grid"))) {
    stop("fast_lambda_tree_fixed_cpp has no lambda-grid argument.",
         call. = FALSE)
  }
  do.call(fn, args)
}

.lfv_matrix <- function(x, n_lambda, n_trait, label) {
  if (is.null(x)) stop("fixed-lambda result is missing ", label, call. = FALSE)
  d <- dim(x)
  if (is.null(d)) {
    if (length(x) != n_lambda * n_trait) {
      stop(label, " has an incompatible length.", call. = FALSE)
    }
    return(matrix(x, nrow = n_lambda, ncol = n_trait))
  }
  if (all(as.integer(d) == c(n_lambda, n_trait))) {
    return(matrix(x, nrow = n_lambda, ncol = n_trait))
  }
  if (all(as.integer(d) == c(n_trait, n_lambda))) {
    return(t(matrix(x, nrow = n_trait, ncol = n_lambda)))
  }
  stop(label, " has dimensions ", paste(d, collapse = " x "),
       "; expected lambda x trait.", call. = FALSE)
}

.lfv_unpack <- function(result, n_lambda, n_trait) {
  lam <- .lfv_get(result, c("lambda", "lambdas", "lambda_grid"))
  list(
    lambda = if (is.null(lam)) NULL else as.numeric(lam),
    gls_mean = .lfv_matrix(
      .lfv_get(result, c("gls_mean", "gls.mean", "mean", "mu", "a")),
      n_lambda, n_trait, "gls_mean"
    ),
    sigma2 = .lfv_matrix(
      .lfv_get(result, c("sigma2", "sigma_sq", "sig2", "variance")),
      n_lambda, n_trait, "sigma2"
    ),
    logLik = .lfv_matrix(
      .lfv_get(result, c("logLik", "loglik", "logL", "log_lik")),
      n_lambda, n_trait, "logLik"
    ),
    valid = {
      z <- .lfv_get(result, c("valid", "is_valid", "ok"))
      if (is.null(z)) NULL else .lfv_matrix(z, n_lambda, n_trait, "valid")
    },
    status = {
      z <- .lfv_get(result, c("status", "statuses", "message"))
      if (is.null(z)) NULL else .lfv_matrix(z, n_lambda, n_trait, "status")
    }
  )
}

.lfv_tree <- function(kind, n = 8L) {
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
    strongly_unbalanced = make_stree("left"),
    ultrametric = ape::compute.brlen(make_stree("balanced"), method = "Grafen"),
    nonultrametric = ape::rtree(n),
    very_short_positive = make_stree("balanced"),
    heterogeneous_positive = make_stree("left"),
    two_tip = ape::read.tree(text = "(sp1:0.4,sp2:1.7);"),
    three_tip = ape::read.tree(text = "((sp1:0.2,sp2:1.1):0.7,sp3:2.3);"),
    polytomy = ape::read.tree(
      text = "(sp1:0.4,sp2:0.7,sp3:1.2,sp4:1.8);"
    ),
    stop("unknown fixed-lambda tree fixture", call. = FALSE)
  )
  if (is.null(tree$tip.label)) {
    tree$tip.label <- paste0("sp", seq_len(ape::Ntip(tree)))
  }
  tree$tip.label <- paste0("sp", seq_len(ape::Ntip(tree)))
  if (is.null(tree$edge.length)) {
    tree$edge.length <- rep(1, nrow(tree$edge))
  }
  tree$edge.length <- as.numeric(tree$edge.length)
  if (kind == "very_short_positive") {
    tree$edge.length <- seq(1e-8, 4e-7, length.out = nrow(tree$edge))
  }
  if (kind == "heterogeneous_positive") {
    tree$edge.length <- 10^seq(-5, 2, length.out = nrow(tree$edge))
  }
  tree$edge.length <- pmax(tree$edge.length, 1e-12)
  ape::reorder.phylo(tree, "cladewise")
}

.lfv_traits <- function(tree, p = 6L, seed = 1L) {
  set.seed(seed)
  n <- ape::Ntip(tree)
  p <- max(1L, as.integer(p))
  X <- matrix(
    stats::rnorm(n * p), nrow = n, ncol = p,
    dimnames = list(tree$tip.label, paste0("trait_", seq_len(p)))
  )
  if (p >= 2L) X[, 2L] <- 0.75 + stats::rnorm(n, sd = 1e-7)
  if (p >= 3L) X[, 3L] <- 4
  if (p >= 4L) X[, 4L] <- 1e12 + stats::rnorm(n)
  if (p >= 5L) X[, 5L] <- seq_len(n) / max(1, n)
  if (p >= 6L) X[, 6L] <- -2 + stats::rnorm(n, sd = 1e-4)
  X
}

.lfv_grid <- function(tree) {
  upper <- if (exists(".max_lambda", envir = asNamespace("fastphylosig"),
                     inherits = FALSE)) {
    fastphylosig:::.max_lambda(tree)
  } else 1
  upper <- as.numeric(upper[[1L]])
  if (!is.finite(upper) || upper <= 0) upper <- 1
  out <- c(0, 1e-12, 1e-8, 0.25, 0.5, 0.9,
           max(0, upper - 1e-10), upper)
  if (upper >= 1) out <- c(out, 1 - 1e-10, 1)
  sort(unique(out[is.finite(out) & out >= 0 & out <= upper]))
}

.lfv_oracle <- function(tree, y, lambda) {
  y <- as.numeric(y)
  n <- length(y)
  C <- ape::vcv.phylo(tree)
  Cl <- lambda * C
  diag(Cl) <- diag(C)
  inv <- tryCatch(solve(Cl), error = function(e) NULL)
  if (is.null(inv)) {
    return(list(mean = NA_real_, sigma2 = NA_real_, logLik = -Inf,
                valid = FALSE))
  }
  one <- rep(1, n)
  sum_inv <- sum(inv)
  if (!is.finite(sum_inv) || sum_inv <= 0) {
    return(list(mean = NA_real_, sigma2 = NA_real_, logLik = -Inf,
                valid = FALSE))
  }
  # Keep the operation order stable for common offsets of 1e12 and above.
  baseline <- y[[1L]]
  shifted <- y - baseline
  delta <- as.numeric(crossprod(one, inv %*% shifted) / sum_inv)
  gls_mean <- baseline + delta
  centered <- shifted - delta
  quad <- as.numeric(crossprod(centered, inv %*% centered))
  sigma2 <- quad / n
  if (!is.finite(sigma2) || sigma2 <= 0) {
    return(list(mean = gls_mean, sigma2 = NA_real_, logLik = -Inf,
                valid = FALSE))
  }
  logdet <- tryCatch(
    as.numeric(determinant(Cl, logarithm = TRUE)$modulus),
    error = function(e) NA_real_
  )
  if (!is.finite(logdet)) {
    return(list(mean = NA_real_, sigma2 = NA_real_, logLik = -Inf,
                valid = FALSE))
  }
  logLik <- -0.5 * quad / sigma2 -
    0.5 * n * log(2 * pi) -
    0.5 * (n * log(sigma2) + logdet)
  if (!is.finite(logLik)) {
    return(list(mean = gls_mean, sigma2 = sigma2, logLik = -Inf,
                valid = FALSE))
  }
  list(mean = gls_mean, sigma2 = sigma2, logLik = as.numeric(logLik),
       valid = TRUE)
}

.lfv_expect_close <- function(got, ref, abs_tol = 1e-7, rel_tol = 1e-7,
                              label = "value") {
  if (is.finite(ref)) {
    testthat::expect_true(is.finite(got), label = label)
    limit <- abs_tol + rel_tol * max(1, abs(ref))
    testthat::expect_lte(abs(got - ref), limit, label = label)
  } else if (is.infinite(ref) && ref < 0) {
    testthat::expect_true(is.infinite(got) && got < 0, label = label)
  } else {
    testthat::expect_true(is.na(got), label = label)
  }
}

.lfv_compiled_subset <- function(tree) {
  ctx <- fastphylosig::prepare_tree(tree)
  fastphylosig:::.prepared_tree_subset(
    ctx, seq_len(ape::Ntip(tree)), need_lambda = FALSE, need_matrix = FALSE
  )
}


test_that("fixed lambda tree likelihood matches dense GLS over the allowed grid", {
  testthat::skip_if_not_installed("ape")
  if (!.lfv_has_engine()) {
    testthat::skip("fast_lambda_tree_fixed_cpp is not integrated yet")
  }

  families <- list(
    balanced = .lfv_tree("balanced", 8L),
    pectinate = .lfv_tree("pectinate", 8L),
    strongly_unbalanced = .lfv_tree("strongly_unbalanced", 9L),
    ultrametric = .lfv_tree("ultrametric", 8L),
    nonultrametric = .lfv_tree("nonultrametric", 8L),
    very_short_positive = .lfv_tree("very_short_positive", 8L),
    heterogeneous_positive = .lfv_tree("heterogeneous_positive", 8L),
    two_tip = .lfv_tree("two_tip", 2L),
    three_tip = .lfv_tree("three_tip", 3L),
    polytomy = .lfv_tree("polytomy", 4L)
  )

  for (nm in names(families)) {
    tree <- families[[nm]]
    X <- .lfv_traits(tree, p = 6L, seed = 20260820L + nchar(nm))
    lambda <- .lfv_grid(tree)
    subset <- .lfv_compiled_subset(tree)
    Xg <- X[subset$tree$tip.label, , drop = FALSE]
    result <- .lfv_unpack(
      .lfv_call(subset$compiled_tree, Xg, lambda,
                trait_chunk = 2L, n_threads = 1L),
      length(lambda), ncol(Xg)
    )
    if (!is.null(result$lambda)) {
      testthat::expect_equal(result$lambda, lambda, tolerance = 1e-12,
                             label = paste(nm, "lambda grid"))
    }
    for (i in seq_along(lambda)) {
      for (j in seq_len(ncol(Xg))) {
        ref <- .lfv_oracle(subset$tree, Xg[, j], lambda[[i]])
        # A constant trait has a finite GLS mean but an undefined ML variance;
        # the fixed kernel must report NA sigma2 and -Inf logLik, explicitly.
        .lfv_expect_close(
          result$gls_mean[i, j], ref$mean,
          abs_tol = if (j == 4L) 2e-4 else 2e-6,
          rel_tol = if (j == 4L) 2e-12 else 2e-8,
          label = paste(nm, "gls_mean", i, j)
        )
        .lfv_expect_close(
          result$sigma2[i, j], ref$sigma2,
          abs_tol = 2e-8, rel_tol = 2e-6,
          label = paste(nm, "sigma2", i, j)
        )
        .lfv_expect_close(
          result$logLik[i, j], ref$logLik,
          abs_tol = 2e-7, rel_tol = 2e-7,
          label = paste(nm, "logLik", i, j)
        )
        if (!is.null(result$valid)) {
          testthat::expect_identical(
            isTRUE(result$valid[i, j]), isTRUE(ref$valid),
            label = paste(nm, "valid", i, j)
          )
        }
      }
    }
    # The retained dense C++ likelihood is an independent implementation of
    # the same operation order.  Check it for the ordinary (non-degenerate)
    # trait as a second oracle on every topology.
    # The retained Cholesky oracle is also checked on well-conditioned trees.
    # Very-short and highly heterogeneous fixtures intentionally stress the
    # tree kernel beyond the dense solver's conditioning envelope, so the
    # stable R dense oracle above is the acceptance reference there.
    if (nm %in% c("balanced", "pectinate", "strongly_unbalanced",
                  "ultrametric", "nonultrametric", "two_tip",
                  "three_tip", "polytomy") &&
        exists("lambda_loglik_cpp", envir = asNamespace("fastphylosig"),
               inherits = FALSE)) {
      dense_cpp <- vapply(lambda, function(z) {
        fastphylosig:::lambda_loglik_cpp(z, ape::vcv.phylo(subset$tree),
                                         Xg[, 1L])
      }, numeric(1))
      for (i in seq_along(lambda)) {
        .lfv_expect_close(result$logLik[i, 1L], dense_cpp[[i]],
                          abs_tol = 2e-6, rel_tol = 2e-7,
                          label = paste(nm, "dense C++ logLik", i))
      }
    }
  }
})


test_that("fixed lambda handles many traits, repeated NA masks, and chunking", {
  testthat::skip_if_not_installed("ape")
  if (!.lfv_has_engine()) {
    testthat::skip("fast_lambda_tree_fixed_cpp is not integrated yet")
  }
  tree <- .lfv_tree("nonultrametric", 24L)
  X <- .lfv_traits(tree, p = 24L, seed = 20260821L)
  masks <- list(
    integer(), 1:2, 3:5, c(2L, 8L, 17L), 6:10,
    c(1L, 12L, 24L), 14:16, 19:22
  )
  for (j in seq_len(ncol(X))) {
    drop <- masks[[(j - 1L) %% length(masks) + 1L]]
    if (length(drop)) X[drop, j] <- NA_real_
  }
  present <- !is.na(X)
  # Grouping is intentionally done in R here; the production caller can use
  # packed masks, but this test checks that every resulting subset is correct.
  keys <- apply(present, 2L, paste, collapse = ",")
  for (key in unique(keys)) {
    cols <- which(keys == key)
    keep <- which(present[, cols[[1L]]])
    if (length(keep) < 2L) next
    drop_labels <- setdiff(tree$tip.label, tree$tip.label[keep])
    tree_j <- if (length(drop_labels)) {
      ape::drop.tip(tree, drop_labels)
    } else tree
    subset <- .lfv_compiled_subset(tree_j)
    Xg <- X[subset$tree$tip.label, cols, drop = FALSE]
    lambda <- .lfv_grid(tree_j)
    one <- .lfv_unpack(
      .lfv_call(subset$compiled_tree, Xg, lambda,
                trait_chunk = 1L, n_threads = 1L),
      length(lambda), ncol(Xg)
    )
    many <- .lfv_unpack(
      .lfv_call(subset$compiled_tree, Xg, lambda,
                trait_chunk = 7L, n_threads = 2L),
      length(lambda), ncol(Xg)
    )
    testthat::expect_equal(one$gls_mean, many$gls_mean, tolerance = 2e-9,
                           label = paste("gls chunk", key))
    testthat::expect_equal(one$sigma2, many$sigma2, tolerance = 2e-9,
                           label = paste("sigma2 chunk", key))
    testthat::expect_equal(one$logLik, many$logLik, tolerance = 2e-9,
                           label = paste("logLik chunk", key))
    for (i in seq_along(lambda)) {
      for (j in seq_len(ncol(Xg))) {
        ref <- .lfv_oracle(tree_j, Xg[, j], lambda[[i]])
        .lfv_expect_close(one$gls_mean[i, j], ref$mean,
                          abs_tol = if (cols[[j]] == 4L) 2e-4 else 2e-6,
                          rel_tol = if (cols[[j]] == 4L) 2e-12 else 2e-8,
                          label = paste("NA gls_mean", key, i, j))
        .lfv_expect_close(one$sigma2[i, j], ref$sigma2,
                          abs_tol = 2e-8, rel_tol = 2e-6,
                          label = paste("NA sigma2", key, i, j))
        .lfv_expect_close(one$logLik[i, j], ref$logLik,
                          abs_tol = 2e-7, rel_tol = 2e-7,
                          label = paste("NA logLik", key, i, j))
      }
    }
  }
})


test_that("fixed lambda marks degenerate and boundary cells explicitly", {
  testthat::skip_if_not_installed("ape")
  if (!.lfv_has_engine()) {
    testthat::skip("fast_lambda_tree_fixed_cpp is not integrated yet")
  }
  tree <- .lfv_tree("balanced", 6L)
  X <- cbind(
    constant = rep(7, ape::Ntip(tree)),
    near_constant = 7 + stats::rnorm(ape::Ntip(tree), sd = 1e-9),
    random = stats::rnorm(ape::Ntip(tree))
  )
  rownames(X) <- tree$tip.label
  lambda <- c(0, 1e-10, 0.5, 1 - 1e-10, 1)
  subset <- .lfv_compiled_subset(tree)
  out <- .lfv_unpack(
    .lfv_call(subset$compiled_tree, X, lambda,
              trait_chunk = 3L, n_threads = 1L),
    length(lambda), ncol(X)
  )
  testthat::expect_true(all(is.na(out$sigma2[, 1L])))
  testthat::expect_true(all(is.infinite(out$logLik[, 1L]) &
                              out$logLik[, 1L] < 0))
  testthat::expect_true(all(is.finite(out$gls_mean[, 1L])))
  testthat::expect_true(all(is.finite(out$sigma2[, 2:3, drop = FALSE])))
  testthat::expect_true(all(is.finite(out$logLik[, 2:3, drop = FALSE])))
  if (!is.null(out$valid)) {
    testthat::expect_false(any(out$valid[, 1L]))
    testthat::expect_true(all(out$valid[, 2:3, drop = FALSE]))
  }
  if (!is.null(out$status)) {
    testthat::expect_true(all(nzchar(as.character(out$status))))
  }
})
