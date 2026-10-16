.hetero_dense_k <- function(tree, x) {
  C <- ape::vcv.phylo(tree)
  y <- x[rownames(C)]
  inv <- solve(C)
  a <- sum(inv %*% y) / sum(inv)
  z <- y - a
  sum(z^2) / as.numeric(crossprod(z, inv %*% z)) /
    ((sum(diag(C)) - length(y) / sum(inv)) / (length(y) - 1))
}

test_that("automatic selection is disjoint, bounded and has an independent K oracle", {
  set.seed(912)
  tree <- ape::rtree(150)
  x <- setNames(rnorm(150), tree$tip.label)
  fit <- fast_hetero(tree, x, min_tips = 10, max_tips = 40,
                                   progress = FALSE)
  expect_false(fit$selection$overlap)
  expect_true(all(fit$clades$n_tip_total >= 10 & fit$clades$n_tip_total <= 40))
  tips <- unlist(fit$selection$tips, use.names = FALSE)
  expect_equal(anyDuplicated(tips), 0L)
  expect_setequal(c(tips, fit$selection$uncovered_tips), tree$tip.label)
  expect_equal(fit$selection$coverage, length(tips) / 150)
  ref <- vapply(fit$selection$nodes, function(node) {
    sub <- ape::extract.clade(tree, node)
    .hetero_dense_k(sub, x)
  }, numeric(1))
  expect_equal(fit$clades$estimate, ref, tolerance = 1e-9)
  expect_equal(fit$heterogeneity$dispersion, sd(log(ref)))
  expect_equal(fit$heterogeneity$status, "descriptive")
  expect_true(all(is.na(fit$heterogeneity$P_BM)))
})

test_that("explicit nodes use source IDs and nested scans remain descriptive", {
  z <- .hetero_fixture()
  tree <- z$tree
  n <- length(tree$tip.label)
  # A representation change that canonicalization repairs by renumbering.
  ids <- c(seq_len(n), n + rev(seq_len(tree$Nnode)))
  tree$edge[] <- ids[tree$edge]
  tree$edge <- tree$edge[nrow(tree$edge):1L, , drop = FALSE]
  nodes <- ids[z$nodes]
  fit <- fast_hetero(tree, z$X, nodes = nodes,
                                   min_tips = 10, progress = FALSE)
  expect_identical(fit$selection$nodes, as.integer(nodes))
  for (i in seq_along(nodes)) {
    sub <- ape::extract.clade(z$tree, z$nodes[[i]])
    expect_setequal(fit$selection$tips[[i]], sub$tip.label)
    expect_equal(fit$clades$estimate[[i]], .hetero_dense_k(sub, z$X[, 1]), tolerance = 1e-9)
  }
  nested <- c(65L, z$nodes[[1L]])
  fit <- fast_hetero(z$tree, z$X, nodes = nested, min_tips = 10,
                                   progress = FALSE)
  expect_true(fit$selection$overlap)
  expect_error(fast_hetero(z$tree, z$X, nodes = nested,
    min_tips = 10, test = TRUE, progress = FALSE), "non-overlapping")
  expect_error(fast_hetero(z$tree, z$X, nodes = nested,
    min_tips = 10, max_total_tips = 64, progress = FALSE), "max_total_tips")
})

test_that("lambda scan agrees with the specialist and accepts prepared contexts", {
  z <- .hetero_fixture()
  fit <- fast_hetero(z$tree, z$X, nodes = z$nodes,
    method = "lambda", min_tips = 10, progress = FALSE)
  ref <- vapply(z$nodes, function(node) {
    sub <- ape::extract.clade(z$tree, node)
    fast_lambda(sub, z$X[, 1], lambda_profile = FALSE, progress = FALSE)$lambda
  }, numeric(1))
  expect_equal(fit$clades$estimate, ref)
  ctx <- prepare_tree(z$tree)
  prepared <- fast_hetero(ctx, z$X, min_tips = 10, max_tips = 32,
                                        progress = FALSE)
  raw <- fast_hetero(ctx$tree, z$X, min_tips = 10, max_tips = 32,
                                   progress = FALSE)
  expect_equal(prepared$clades, raw$clades)
  set.seed(30)
  calibrated <- fast_hetero(z$tree, z$X, method = "lambda", nodes = z$nodes,
    min_tips = 10, test = TRUE, nsim = 9, progress = FALSE)
  expect_equal(calibrated$heterogeneity$status, "ok")
  expect_true(is.finite(calibrated$heterogeneity$P_heterogeneity))
})

test_that("matching, missingness and failures are trait-wise", {
  z <- .hetero_fixture()
  X <- cbind(good = z$X[, 1], missing = z$X[, 1], constant = 1, invalid = z$X[, 1])
  X[1:25, "missing"] <- NA_real_
  X[1, "invalid"] <- Inf
  X <- rbind(X[-64, , drop = FALSE], extra = 5)
  fit <- fast_hetero(z$tree, X[nrow(X):1, ], nodes = z$nodes,
    min_tips = 10, test = TRUE, nsim = 9, progress = FALSE)
  expect_identical(fit$selection$tree_tips_without_data, z$tree$tip.label[64])
  expect_identical(fit$selection$data_rows_removed, "extra")
  expect_equal(fit$clades$n_species[fit$clades$trait == "missing"], c(7L, 31L))
  expect_equal(fit$clades$status[fit$clades$trait == "missing"], c("insufficient_data", "ok"))
  expect_true(all(!is.finite(fit$clades$estimate[fit$clades$trait == "constant"])))
  expect_equal(fit$clades$status[fit$clades$trait == "invalid"], c("invalid_trait", "ok"))
  expect_equal(fit$heterogeneity$status, c("ok", rep("insufficient_clades", 3)))
  expect_true(is.finite(fit$heterogeneity$P_BM[[1L]]))
})

test_that("BM simulator agrees with branch increments and dense covariance", {
  z <- .hetero_fixture()
  sub <- ape::extract.clade(z$tree, z$nodes[[1L]])
  core <- compile_tree_cpp(sub$edge, sub$edge.length, length(sub$tip.label))
  set.seed(45)
  draws <- .hetero_bm_reference(core, 3L)
  set.seed(45)
  children <- core$preorder[-1L]
  increments <- matrix(rnorm(length(children) * 3), nrow = length(children)) *
    sqrt(core$branch_length_by_node[children])
  # Independent path-summing oracle, with no upward-message implementation.
  ref <- matrix(0, nrow = core$n_tip, ncol = 3L)
  for (tip in seq_len(core$n_tip)) {
    node <- tip
    while (node != core$root) {
      ref[tip, ] <- ref[tip, ] + increments[match(node, children), ]
      node <- core$parent[[node]]
    }
  }
  expect_equal(draws, ref)
  set.seed(46)
  many <- .hetero_bm_reference(core, 20000L)
  C <- ape::vcv.phylo(sub)
  C <- C[sub$tip.label, sub$tip.label]
  expect_lt(max(abs(stats::cov(t(many)) - C)), 0.22)
})

test_that("BM dispersion and probability agree with a dense Monte Carlo oracle", {
  z <- .hetero_fixture()
  X <- z$X
  X[c(2, 4), 1] <- NA_real_
  set.seed(64)
  fit <- fast_hetero(z$tree, X, nodes = z$nodes, min_tips = 10,
    test = TRUE, nsim = 29, keep_null = TRUE, simulation_chunk = 4, progress = FALSE)
  set.seed(64)
  dense <- matrix(NA_real_, nrow = 29, ncol = 2L)
  for (i in seq_along(z$nodes)) {
    sub <- ape::extract.clade(z$tree, z$nodes[[i]])
    keep <- sub$tip.label[!is.na(X[sub$tip.label, 1])]
    if (length(keep) != length(sub$tip.label)) sub <- ape::keep.tip(sub, keep)
    core <- compile_tree_cpp(sub$edge, sub$edge.length, length(sub$tip.label))
    draws <- .hetero_bm_reference(core, 29L)
    dense[, i] <- apply(draws, 2L, function(x) .hetero_dense_k(sub, setNames(x, sub$tip.label)))
  }
  H <- apply(log(dense), 1L, sd)
  expect_equal(fit$null[[1]], H, tolerance = 1e-9)
  count <- sum(H >= fit$heterogeneity$dispersion)
  expect_equal(fit$heterogeneity$P_BM, (count + 1) / 30)
  expect_equal(fit$heterogeneity$exceedance_count, count)
  expect_equal(fit$heterogeneity$nsim_successful, 29L)
  for (i in 1:2) {
    expect_equal(fit$clades$BM_mean[[i]], mean(dense[, i]), tolerance = 1e-9)
    expect_equal(c(fit$clades$BM_lower[[i]], fit$clades$BM_upper[[i]]),
                 unname(quantile(dense[, i], c(.025, .975))), tolerance = 1e-9)
  }
  seed_after <- .Random.seed
  set.seed(64)
  other <- fast_hetero(z$tree, X, nodes = z$nodes, min_tips = 10,
    test = TRUE, nsim = 29, keep_null = TRUE, simulation_chunk = 11, progress = FALSE)
  expect_identical(.Random.seed, seed_after)
  expect_equal(other$null, fit$null, tolerance = 1e-9)
})

test_that("contrasting clade structure gives excess BM dispersion and K is scale invariant", {
  z <- .hetero_fixture()
  # First crown clade varies at its deepest split, second alternates at tips.
  x <- setNames(c(rep(c(-1, 1), each = 16), rep(c(-1, 1), 16)), z$tree$tip.label)
  set.seed(77)
  fit <- fast_hetero(z$tree, x, nodes = z$nodes, min_tips = 10,
    test = TRUE, nsim = 99, progress = FALSE)
  expect_gt(fit$clades$estimate[[1]], fit$clades$estimate[[2]] * 10)
  expect_lte(fit$heterogeneity$P_BM, .05)
  set.seed(77)
  scaled <- fast_hetero(z$tree, 30 + 4 * x, nodes = z$nodes,
    min_tips = 10, test = TRUE, nsim = 99, progress = FALSE)
  expect_equal(scaled$clades, fit$clades, tolerance = 1e-9)
  expect_equal(scaled$heterogeneity, fit$heterogeneity, tolerance = 1e-9)
})

test_that("contracts, print and plot are usable", {
  z <- .hetero_fixture()
  expect_error(fast_hetero(z$tree, z$X, min_tips = 2), "min_tips")
  expect_error(fast_hetero(z$tree, z$X, nodes = 1), "internal node")
  expect_error(fast_hetero(z$tree, z$X, nodes = rep(65, 2)), "unique")
  expect_error(fast_hetero(z$tree, z$X, test = NA), "TRUE or FALSE")
  expect_error(fast_hetero(z$tree, z$X, test = TRUE, nsim = 1), "nsim")
  expect_error(fast_hetero(z$tree, z$X, min_tips = 30, max_tips = 31), "No clades")
  fit <- fast_hetero(z$tree, z$X, nodes = z$nodes, min_tips = 10,
    progress = FALSE)
  expect_output(print(fit), "Clade-scale")
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  expect_equal(plot(fit)$estimate, fit$clades$estimate)
  expect_error(plot(fit, trait = "unknown"), "trait must name")
  X <- cbind(z$X, 2 * z$X)
  colnames(X) <- c("trait", "other")
  multi <- fast_hetero(z$tree, X,
    nodes = z$nodes, min_tips = 10, progress = FALSE)
  expect_error(plot(multi), "Specify trait")
})

test_that("direct clade kernels preserve specialist results for NA groups and failures", {
  z <- .hetero_fixture()
  X <- cbind(good = z$X[, 1], shifted = 1e6 + z$X[, 1],
             missing = z$X[, 1], same_mask = 2 * z$X[, 1],
             other_mask = z$X[, 1], constant = 1, infinite = z$X[, 1],
             nan = z$X[, 1])
  X[c(2, 5, 36), c("missing", "same_mask")] <- NA_real_
  X[c(6, 10, 48), "other_mask"] <- NA_real_
  X[1, "infinite"] <- Inf
  X[33, "nan"] <- NaN
  for (method in c("K", "lambda")) {
    scan <- fast_hetero(z$tree, X, nodes = z$nodes,
      method = method, min_tips = 10, progress = FALSE)
    for (i in seq_along(z$nodes)) {
      sub <- ape::extract.clade(z$tree, z$nodes[[i]])
      ref <- if (method == "K") fast_k(sub, X, progress = FALSE) else
        fast_lambda(sub, X, lambda_profile = FALSE, progress = FALSE)
      rows <- scan$clades$node == z$nodes[[i]]
      expect_equal(scan$clades$estimate[rows],
        if (method == "K") ref$K_fast else ref$lambda_fast, tolerance = 1e-8)
      expect_identical(scan$clades$status[rows], ref$status)
      expect_identical(scan$clades$message[rows], ref$message)
    }
  }
})

test_that("lightweight scanning preserves input and prepared snapshot checks", {
  z <- .hetero_fixture()
  ctx <- prepare_tree(z$tree)
  ctx$tree$edge.length[[1L]] <- 2
  expect_error(fast_hetero(ctx, z$X, min_tips = 10, progress = FALSE))
  for (change in list(
      function(t) { t$tip.label[1] <- t$tip.label[2]; t },
      function(t) { t$edge.length[1] <- -1; t },
      function(t) { t$edge[1, 2] <- t$edge[2, 2]; t },
      function(t) { t$Nnode <- t$Nnode + 1L; t })) {
    expect_error(fast_hetero(change(z$tree), z$X,
      nodes = z$nodes, min_tips = 10, progress = FALSE), "cannot compile")
  }
  expect_error(fast_hetero(list(), z$X, progress = FALSE), "cannot compile")
})

test_that("local zero branch diagnostics survive the direct kernel path", {
  z <- .hetero_fixture()
  for (terminal in c(TRUE, FALSE)) {
    tree <- z$tree
    branch <- which((tree$edge[, 2] <= 64L) == terminal &
                    tree$edge[, 1] != 65L)[[1L]]
    tree$edge.length[branch] <- 0
    for (method in c("K", "lambda")) {
      fit <- fast_hetero(tree, z$X, nodes = z$nodes,
        method = method, min_tips = 10, progress = FALSE)
      expect_equal(fit$clades$status[[2L]], "ok")
      if (method == "K" || terminal) {
        expect_equal(fit$clades$status[[1L]], "invalid_subtree")
        expect_match(fit$clades$message[[1L]], "zero")
      } else {
        sub <- ape::extract.clade(tree, z$nodes[[1L]])
        ref <- fast_lambda(sub, z$X[, 1], lambda_profile = FALSE, progress = FALSE)
        expect_equal(fit$clades$estimate[[1L]], ref$lambda)
      }
    }
  }
})
