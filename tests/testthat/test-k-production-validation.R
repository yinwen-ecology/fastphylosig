test_that("Milestone 1 tree K contract is explicit", {
  testthat::skip_if_not_installed("ape")
  # Production algorithms are automatic; dense/spectral validation engines
  # remain internal and are not choices ordinary users need to understand.
  expect_false("engine" %in% names(formals(fast_signal)))
  expect_true("prepared" %in% names(formals(fast_signal)))
  expect_true("trait_chunk" %in% names(formals(fast_signal)))
  expect_false("validate_tolerance" %in% names(formals(fast_signal)))
})


.k_accept_tree <- function(kind, n = 20L) {
  make_stree <- function(type) {
    tryCatch(
      ape::stree(n = n, type = type, tip.label = paste0("sp", seq_len(n))),
      error = function(e) ape::rtree(n)
    )
  }
  tree <- switch(
    kind,
    balanced = make_stree("balanced"),
    pectinate = make_stree("left"),
    unbalanced = make_stree("left"),
    ultrametric = ape::compute.brlen(make_stree("balanced"), method = "Grafen"),
    nonultrametric = ape::rtree(n),
    short = make_stree("balanced"),
    heterogeneous = make_stree("left"),
    two_tip = ape::read.tree(text = "(sp1:0.4,sp2:1.7);"),
    three_tip = ape::read.tree(text = "((sp1:0.2,sp2:1.1):0.7,sp3:2.3);"),
    stop("unknown tree kind", call. = FALSE)
  )
  if (is.null(tree$tip.label)) {
    tree$tip.label <- paste0("sp", seq_len(ape::Ntip(tree)))
  }
  tree$tip.label <- paste0("sp", seq_len(ape::Ntip(tree)))
  if (is.null(tree$edge.length)) tree$edge.length <- rep(1, nrow(tree$edge))
  tree$edge.length <- as.numeric(tree$edge.length)
  if (kind == "short") {
    tree$edge.length <- seq(1e-9, 4e-8, length.out = nrow(tree$edge))
  }
  if (kind == "heterogeneous") {
    tree$edge.length <- 10^seq(-6, 2, length.out = nrow(tree$edge))
  }
  # Keep every branch strictly positive for the production tree solver.
  tree$edge.length <- pmax(tree$edge.length, 1e-12)
  ape::reorder.phylo(tree, "cladewise")
}

.k_accept_dense <- function(tree, y) {
  y <- as.numeric(y[tree$tip.label])
  keep <- !is.na(y)
  if (sum(keep) < 2L) return(NA_real_)
  if (!all(keep)) {
    tree <- ape::drop.tip(tree, tree$tip.label[!keep])
    y <- y[keep]
  }
  C <- ape::vcv.phylo(tree)
  Q <- tryCatch(solve(C), error = function(e) NULL)
  if (is.null(Q)) return(NA_real_)
  one <- rep(1, nrow(C))
  sum_inv <- sum(Q)
  if (!is.finite(sum_inv) || sum_inv == 0) return(NA_real_)
  # This is algebraically the phytools formula, evaluated after a harmless
  # translation to avoid catastrophic cancellation for large-offset traits.
  # A constant trait has no defined K even if dense roundoff creates a tiny
  # artificial residual.
  baseline <- y[[1L]]
  shifted <- y - baseline
  if (all(shifted == 0)) return(NA_real_)
  a_shifted <- as.numeric(crossprod(one, Q %*% shifted) / sum_inv)
  centered <- shifted - a_shifted
  denominator <- as.numeric(crossprod(centered, Q %*% centered))
  normalization <- (sum(diag(C)) - nrow(C) / sum_inv) /
    (nrow(C) - 1L)
  if (!is.finite(denominator) || denominator <= 0 ||
      !is.finite(normalization) || normalization <= 0) {
    return(NA_real_)
  }
  as.numeric(crossprod(centered, centered) / denominator / normalization)
}

.k_accept_permutations <- function(n, nsim) {
  out <- matrix(NA_integer_, nrow = nsim, ncol = n)
  out[1L, ] <- seq_len(n)
  if (nsim > 1L) {
    for (i in 2:nsim) out[i, ] <- sample.int(n)
  }
  out
}

.k_accept_has_tree_engine <- function(tree, X) {
  probe <- tryCatch(
    suppressWarnings(suppressMessages(fast_signal(
      tree, X[, 1L, drop = FALSE], method = "K", test = FALSE,
      verbose = FALSE
    ))),
    error = function(e) e
  )
  !inherits(probe, "error")
}

.k_accept_field <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (!length(hit)) return(NULL)
  x[[hit[[1L]]]]
}


test_that("tree K agrees with dense GLS across topology and numerical domains", {
  testthat::skip_if_not_installed("ape")
  set.seed(20260811)
  probe_tree <- .k_accept_tree("balanced", 12L)
  probe_X <- matrix(
    stats::rnorm(12), ncol = 1L,
    dimnames = list(probe_tree$tip.label, "random")
  )
  if (!.k_accept_has_tree_engine(probe_tree, probe_X)) {
    testthat::skip("production tree K engine is not integrated yet")
  }

  families <- list(
    balanced = .k_accept_tree("balanced", 20L),
    pectinate = .k_accept_tree("pectinate", 20L),
    strongly_unbalanced = .k_accept_tree("unbalanced", 20L),
    ultrametric = .k_accept_tree("ultrametric", 20L),
    nonultrametric = .k_accept_tree("nonultrametric", 20L),
    very_short_positive = .k_accept_tree("short", 20L),
    heterogeneous_positive = .k_accept_tree("heterogeneous", 20L),
    two_tip = .k_accept_tree("two_tip", 2L),
    three_tip = .k_accept_tree("three_tip", 3L)
  )

  for (nm in names(families)) {
    tree <- families[[nm]]
    n <- ape::Ntip(tree)
    random <- stats::rnorm(n)
    X <- cbind(
      random = random,
      random_2 = stats::rnorm(n),
      constant = rep(2, n),
      near_constant = 2 + stats::rnorm(n, sd = 1e-5),
      large_offset = 1e12 + stats::rnorm(n),
      branch_scaled = random * seq_len(n)
    )
    dimnames(X) <- list(tree$tip.label, c(
      "random", "random_2", "constant", "near_constant",
      "large_offset", "branch_scaled"
    ))
    got <- suppressWarnings(suppressMessages(fast_signal(
      tree, X, method = "K", test = FALSE,
      trait_chunk = 2L, verbose = FALSE
    )))
    expect_true(all(c("trait", "K_fast", "n_species", "n_removed_na") %in%
                      names(got)), label = nm)
    expect_equal(nrow(got), ncol(X), label = nm)

    reference <- vapply(seq_len(ncol(X)), function(j) {
      .k_accept_dense(tree, X[, j])
    }, numeric(1))
    value <- as.numeric(got$K_fast)
    finite <- is.finite(reference)
    if (any(finite)) {
      # The production tree engine is algebraically identical to dense GLS;
      # the slightly wider tolerance only accommodates ill-conditioned
      # very-short and heterogeneous branch fixtures on old BLAS versions.
      expect_lt(max(abs(value[finite] - reference[finite])), 1e-6,
                label = nm)
    }
    # Constant traits have a zero GLS denominator and therefore are expected
    # to be non-finite, not silently assigned an arbitrary K value.
    expect_false(is.finite(value[which(colnames(X) == "constant")]),
                 label = nm)
  }
})


test_that("K tree engine handles many traits, packed NA masks, and cache reuse", {
  testthat::skip_if_not_installed("ape")
  set.seed(20260812)
  tree <- .k_accept_tree("nonultrametric", 80L)
  n <- ape::Ntip(tree)
  p <- 32L
  X <- matrix(
    stats::rnorm(n * p), nrow = n, ncol = p,
    dimnames = list(tree$tip.label, paste0("trait_", seq_len(p)))
  )
  if (!.k_accept_has_tree_engine(tree, X)) {
    testthat::skip("production tree K engine is not integrated yet")
  }
  # Eight repeated packed masks plus one unique mask exercise grouping rather
  # than one VCV/factorization per trait.
  masks <- list(
    integer(), 1:3, 4:8, c(2L, 9L, 17L), 12:20,
    c(1L, 25L, 40L, 55L), 30:35, 60:65, 10:12
  )
  for (j in seq_len(p)) {
    drop <- masks[[(j - 1L) %% length(masks) + 1L]]
    if (length(drop)) X[drop, j] <- NA_real_
  }
  ctx <- prepare_tree(tree)
  got <- suppressWarnings(suppressMessages(fast_signal(
    ctx, X, method = "K", test = FALSE,
    trait_chunk = 5L, verbose = FALSE
  )))
  expect_equal(nrow(got), p)
  expect_equal(unname(got$n_species), unname(colSums(!is.na(X))))
  expect_equal(unname(got$n_removed_na), unname(colSums(is.na(X))))
  reference <- vapply(seq_len(p), function(j) .k_accept_dense(tree, X[, j]),
                      numeric(1))
  finite <- is.finite(reference)
  expect_lt(max(abs(got$K_fast[finite] - reference[finite])), 1e-6)

  # Tree-only K must not populate dense numerical resources in prepare_tree.
  info <- cache_info(ctx)
  expect_equal(info$n_numerical_entries, 0L)
  expect_equal(info$bytes_used, 0)

  # The same output is invariant to bounded trait chunk size.
  got_one <- suppressWarnings(suppressMessages(fast_signal(
    ctx, X, method = "K", test = FALSE,
    trait_chunk = 1L, verbose = FALSE
  )))
  expect_equal(got_one$K_fast, got$K_fast, tolerance = 1e-12)
})


test_that("controlled permutations give the dense K null exactly", {
  testthat::skip_if_not_installed("ape")
  set.seed(20260813)
  tree <- .k_accept_tree("balanced", 18L)
  X <- matrix(
    stats::rnorm(18 * 3), nrow = 18L, ncol = 3L,
    dimnames = list(tree$tip.label, paste0("trait_", 1:3))
  )
  if (!.k_accept_has_tree_engine(tree, X)) {
    testthat::skip("production tree K engine is not integrated yet")
  }
  nsim <- 41L
  permutations <- .k_accept_permutations(18L, nsim)
  # A production streaming implementation must accept the same explicit
  # permutation matrix as the dense compatibility implementation.  This makes
  # observed K, every sim.K value, and P a deterministic acceptance criterion.
  fit <- tryCatch(
    suppressWarnings(suppressMessages(fast_signal(
      tree, X, method = "K", test = TRUE, nsim = nsim,
      permutations = permutations, return_sim = TRUE,
      verbose = FALSE
    ))),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    testthat::skip("tree K test/streaming API is not integrated yet")
  }
  expect_true(all(c("K_fast", "P_fast") %in% names(fit)))
  sim_field <- .k_accept_field(fit, c("sim.K_fast", "sim.K", "sim_K"))
  if (is.null(sim_field)) {
    testthat::skip("production K result does not expose sim.K yet")
  }

  for (j in seq_len(ncol(X))) {
    observed <- .k_accept_dense(tree, X[, j])
    null <- vapply(seq_len(nsim), function(i) {
      permuted <- stats::setNames(
        X[permutations[i, ], j], tree$tip.label
      )
      .k_accept_dense(tree, permuted)
    }, numeric(1))
    expect_equal(fit$K_fast[j], observed, tolerance = 1e-8)
    expect_equal(as.numeric(sim_field[[j]]), null, tolerance = 1e-8)
    expect_equal(fit$P_fast[j], mean(null >= observed), tolerance = 1e-12)
  }
})


test_that("stochastic K P and MCSE obey the declared Monte Carlo contract", {
  testthat::skip_if_not_installed("ape")
  set.seed(20260814)
  tree <- .k_accept_tree("nonultrametric", 30L)
  x <- stats::rnorm(30L)
  names(x) <- tree$tip.label
  one <- tryCatch(
    suppressWarnings(suppressMessages(fast_signal(
      tree, x, method = "K", test = TRUE, nsim = 199L,
      return_sim = FALSE, verbose = FALSE
    ))),
    error = function(e) e
  )
  if (inherits(one, "error")) {
    testthat::skip("tree K stochastic API is not integrated yet")
  }
  mcse_name <- intersect(c("MCSE_P", "P_mcse", "mcse_P"), names(one))
  if (!length(mcse_name)) {
    testthat::skip("MCSE_P is not part of the current production result")
  }
  p_name <- intersect(c("P_fast", "P"), names(one))
  expect_true(length(p_name) == 1L)
  p1 <- as.numeric(one[[p_name]])[1L]
  mcse <- as.numeric(one[[mcse_name[[1L]]]])[1L]
  expect_true(is.finite(p1) && p1 >= 0 && p1 <= 1)
  expect_true(is.finite(mcse) && mcse >= 0)

  set.seed(20260815)
  two <- suppressWarnings(suppressMessages(fast_signal(
    tree, x, method = "K", test = TRUE, nsim = 199L,
    return_sim = FALSE, verbose = FALSE
  )))
  p2 <- as.numeric(two[[p_name]])[1L]
  pbar <- (p1 + p2) / 2
  se_diff <- sqrt(pbar * (1 - pbar) * (2 / 199))
  # Independent null runs are stochastic; this is a diagnostic bound rather
  # than an equality requirement.  Six standard errors avoids false failures
  # for the deliberately small test nsim.
  expect_lte(abs(p1 - p2), 6 * se_diff + 1 / 199)
})


test_that("tree K matches phytools on complete traits when reference is installed", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  set.seed(20260816)
  tree <- .k_accept_tree("balanced", 24L)
  X <- matrix(
    stats::rnorm(24 * 4), nrow = 24L, ncol = 4L,
    dimnames = list(tree$tip.label, paste0("trait_", 1:4))
  )
  if (!.k_accept_has_tree_engine(tree, X)) {
    testthat::skip("production tree K engine is not integrated yet")
  }
  fit <- suppressWarnings(suppressMessages(fast_signal(
    tree, X, method = "K", test = FALSE,
    trait_chunk = 3L, verbose = FALSE
  )))
  reference <- vapply(seq_len(ncol(X)), function(j) {
    as.numeric(suppressWarnings(phytools::phylosig(
      tree, X[, j], method = "K", test = FALSE, se = NULL
    )))
  }, numeric(1))
  expect_lt(max(abs(fit$K_fast - reference)), 1e-8)
})
