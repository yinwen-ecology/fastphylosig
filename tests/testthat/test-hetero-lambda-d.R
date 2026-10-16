.hetero_dense_ll <- function(tree, x, lambda) {
  C <- ape::vcv.phylo(tree)[tree$tip.label, tree$tip.label]
  C <- lambda * (C - diag(diag(C))) + diag(diag(C))
  ans <- tryCatch({
    inv <- solve(C)
    mean <- sum(inv %*% x) / sum(inv)
    s2 <- as.numeric(crossprod(x-mean, inv %*% (x-mean)))/length(x)
    -.5*(length(x)*(log(2*pi)+log(s2)+1) + as.numeric(determinant(C, TRUE)$modulus))
  }, error = function(e) -Inf)
  if (is.finite(ans)) ans else -Inf
}

test_that("shared lambda and its bootstrap match independent dense/R references", {
  z <- .hetero_fixture()
  z$X[c(3, 36), 1] <- NA_real_
  trees <- lapply(z$nodes, function(node) {
    s <- ape::extract.clade(z$tree, node)
    ape::keep.tip(s, s$tip.label[!is.na(z$X[s$tip.label, 1])])
  })
  cores <- lapply(trees, function(s) compile_tree_cpp(s$edge, s$edge.length, length(s$tip.label)))
  pooled <- function(lambda) sum(vapply(trees, function(s)
    .hetero_dense_ll(s, z$X[s$tip.label, 1], lambda), numeric(1)))
  bound <- min(vapply(cores, function(core) {
    if (diff(range(core$root_distance[seq_len(core$n_tip)])) > 1e-7) return(1)
    min(core$root_distance[seq_len(core$n_tip)] / core$root_distance[core$parent[seq_len(core$n_tip)]])
  }, numeric(1)))
  ref <- optimize(pooled, c(0, bound), maximum = TRUE, tol = 1e-7)
  candidates <- c(0, ref$maximum)
  common <- candidates[which.max(vapply(candidates, pooled, numeric(1)))]
  set.seed(531)
  h <- fast_hetero(z$tree, z$X, method = "lambda", nodes = z$nodes,
    min_tips = 10, test = TRUE, nsim = 19, keep_null = TRUE, progress = FALSE,
    lambda_reference = "ML")
  expect_equal(h$heterogeneity$lambda_common, common, tolerance = 2e-6)
  seed <- .Random.seed
  set.seed(531)
  reference <- numeric(19)
  for (b in seq_len(19)) {
    estimates <- vapply(seq_along(cores), function(i) {
      core <- cores[[i]]; nt <- core$n_tip
      variance <- h$heterogeneity$lambda_common*core$branch_length_by_node
      variance[seq_len(nt)] <- variance[seq_len(nt)] +
        (1-h$heterogeneity$lambda_common)*core$root_distance[seq_len(nt)]
      values <- numeric(core$n_total)
      for (node in core$preorder[-1]) values[node] <- values[core$parent[node]] + rnorm(1)*sqrt(variance[node])
      fast_lambda(trees[[i]], setNames(values[seq_len(nt)], trees[[i]]$tip.label),
        lambda_profile = FALSE, progress = FALSE, verbose = FALSE)$estimate
    }, numeric(1))
    reference[b] <- sd(estimates)
  }
  expect_equal(h$null[[1]], reference, tolerance = 1e-10)
  expect_identical(.Random.seed, seed)
  expect_equal(h$heterogeneity$P_heterogeneity, (1+sum(reference >= h$heterogeneity$dispersion))/20)
  expect_true(is.na(h$heterogeneity$P_BM))
  set.seed(531)
  scaled <- fast_hetero(z$tree, 50+3*z$X, method = "lambda", nodes = z$nodes,
    min_tips = 10, test = TRUE, nsim = 19, keep_null = TRUE, progress = FALSE,
    lambda_reference = "ML")
  expect_equal(scaled$null, h$null, tolerance = 1e-8)
})

# Direct R contrast recursion and full count-by-count recalculation oracle.

test_that("D normalization matches independent repeated contrasts", {
  set.seed(923)
  tree <- ape::rtree(8)
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 8L)
  x <- c(0, 0, 0, 1, 1, 1, 0, 1)
  set.seed(45)
  fit <- heterogeneity_d_calibrate_cpp(list(core), list(x), 41L, 0L, FALSE, FALSE)
  seed <- .Random.seed
  rs <- bs <- 0
  set.seed(45)
  for (b in 1:41) {
    order <- 1:8
    for (k in 8:2) { j <- floor(runif(1, 0, k))+1L; order[c(k, j)] <- order[c(j, k)] }
    bm <- numeric(core$n_total)
    for (node in core$preorder[-1]) bm[node] <- bm[core$parent[node]]+rnorm(1)*sqrt(core$branch_length_by_node[node])
    bm_order <- order(bm[1:8])
    r <- v <- numeric(8); r[order[1:4]] <- 1; v[bm_order[1:4]] <- 1
    rs <- rs+.hetero_d_contrast(core, r)
    bs <- bs+.hetero_d_contrast(core, v)
  }
  rs <- rs/41; bs <- bs/41
  expect_equal(fit$D_mean_random, rs, tolerance = 1e-12)
  expect_equal(fit$D_mean_Brownian, bs, tolerance = 1e-12)
  expect_identical(.Random.seed, seed)
  expect_equal(fit$estimate, (.hetero_d_contrast(core, x)-bs)/(rs-bs), tolerance = 1e-12)
})

test_that("D supports binary labels, missing tips, tests, and distinct overall tests", {
  z <- .hetero_fixture()
  x <- setNames(rep(c(0, 1, 0, 0, 1, 1, 0, 1), 8), z$tree$tip.label)
  x[c(2, 38)] <- NA
  set.seed(932)
  h <- fast_hetero(z$tree, x, method = "D", nodes = z$nodes,
    min_tips = 10, test = TRUE, nsim = 29, d_nsim = 99, keep_null = TRUE,
    overall_test = TRUE, progress = FALSE)
  expect_equal(h$clades$n_species, c(31L, 31L))
  expect_identical(h$heterogeneity$status, "ok")
  expect_equal(h$heterogeneity$P_heterogeneity,
    (1+sum(h$null[[1]] >= h$heterogeneity$dispersion))/30)
  expect_equal(h$clades$estimate, (h$clades$D_contrast-h$clades$D_mean_Brownian)/
    (h$clades$D_mean_random-h$clades$D_mean_Brownian))
  seed <- .Random.seed
  set.seed(932)
  labelled <- fast_hetero(z$tree, setNames(ifelse(is.na(x), NA, ifelse(x == 0, "absent", "present")), names(x)),
    method = "D", nodes = z$nodes, min_tips = 10, test = TRUE, nsim = 29,
    d_nsim = 99, keep_null = TRUE, overall_test = TRUE, simulation_chunk = 11,
    progress = FALSE)
  # Display labels differ deliberately; numerical coding/estimators do not.
  comparable <- setdiff(names(h$clades), c("state_counts", "state_frequencies"))
  expect_equal(labelled$clades[, comparable], h$clades[, comparable])
  expect_equal(labelled$heterogeneity, h$heterogeneity)
  expect_equal(labelled$null, h$null)
  expect_identical(.Random.seed, seed)
  file <- tempfile(fileext = ".pdf"); grDevices::pdf(file)
  p <- plot(h, tree = z$tree)
  expect_match(p$summary_text[2], "random.*Brownian")
  expect_match(p$summary_text[4], "Conditional ER heterogeneity P")
  expect_identical(.Random.seed, seed)
  grDevices::dev.off(); unlink(file)
  expect_error(fast_hetero(z$tree, 1:64, method = "D", progress = FALSE), "more than two states")
  expect_error(fast_hetero(z$tree, x, method = "D", d_nsim = 1, progress = FALSE), "d_nsim")
  flat <- fast_hetero(z$tree, rep(1, 64), method = "D", nodes = z$nodes,
    min_tips = 10, test = TRUE, nsim = 9, progress = FALSE)
  expect_identical(flat$clades$status, rep("undefined", 2))
  expect_identical(flat$heterogeneity$status, "insufficient_clades")
  expect_true(is.na(flat$heterogeneity$P_heterogeneity))
})

# Enumerate ancestral states directly: this does not use a pruning likelihood.
.hetero_er_pattern_prob <- function(core, tips, q) {
  internal <- seq.int(core$n_tip+1, core$n_total)
  combinations <- as.matrix(expand.grid(rep(list(0:1), length(internal))))
  probability <- 0
  for (row in seq_len(nrow(combinations))) {
    states <- c(tips, combinations[row, ])
    children <- core$preorder[-1]
    change <- -.5*expm1(-2*q*core$branch_length_by_node[children])
    transition <- ifelse(states[children] == states[core$parent[children]], 1-change, change)
    probability <- probability+.5*prod(transition)
  }
  probability
}

test_that("conditional ER fit and simulated D dispersion agree with exact enumeration", {
  tree <- ape::read.tree(text = "((a:0.2,b:0.4):0.3,(c:0.7,d:0.5):0.6);")
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 4L)
  # A pair of monophyletic and non-monophyletic patterns gives an interior rate.
  data <- list(c(0, 0, 1, 1), c(0, 1, 0, 1))
  set.seed(62)
  fit <- heterogeneity_d_calibrate_cpp(list(core, core), data, 499L, 9999L, TRUE, TRUE)
  expect_true(is.finite(fit$q_common))
  expect_gt(fit$q_common, 0)
  patterns <- as.matrix(expand.grid(rep(list(0:1), 4)))
  patterns <- patterns[rowSums(patterns) == 2, ]
  loglik <- function(r) {
    q <- r/(1-r)
    count_mass <- sum(apply(patterns, 1, function(x) .hetero_er_pattern_prob(core, x, q)))
    sum(vapply(data, function(x) log(.hetero_er_pattern_prob(core, x, q)/count_mass), numeric(1)))
  }
  reference <- optimize(loglik, c(.000001, .999999), maximum = TRUE, tol = 1e-8)
  expect_equal(fit$q_common, reference$maximum/(1-reference$maximum), tolerance = 1e-5)
  probabilities <- apply(patterns, 1, function(x) .hetero_er_pattern_prob(core, x, fit$q_common))
  probabilities <- probabilities/sum(probabilities)
  D <- vapply(1:2, function(i) apply(patterns, 1, function(x) {
    (.hetero_d_contrast(core, x)-fit$D_mean_Brownian[i]) /
      (fit$D_mean_random[i]-fit$D_mean_Brownian[i])
  }), numeric(nrow(patterns)))
  exact_h <- abs(outer(D[, 1], D[, 2], "-"))/sqrt(2)
  weights <- outer(probabilities, probabilities)
  expected <- sum(exact_h*weights)
  variance <- sum((exact_h-expected)^2*weights)
  expect_lt(abs(mean(fit$null)-expected), 5*sqrt(variance/9999))
  # All draws contain both states: no constant-trait draws were dropped.
  expect_identical(fit$nsim_successful, 9999L)
})

test_that("D rejects nested testing and polytomies, and preserves state-label symmetry", {
  z <- .hetero_fixture()
  x <- setNames(rep(c(0, 1), 32), z$tree$tip.label)
  expect_error(fast_hetero(z$tree, x, method = "D", nodes = c(65L, z$nodes[1]),
    min_tips = 10, test = TRUE, progress = FALSE), "non-overlapping")
  star <- ape::stree(8); star$edge.length <- rep(1, 8)
  bad <- fast_hetero(star, rep(c(0, 1), 4), method = "D", nodes = 9L,
    min_tips = 3, d_nsim = 9, progress = FALSE)
  expect_identical(bad$clades$status, "invalid_subtree")
  expect_match(bad$clades$message, "bifurcating")
  set.seed(173)
  a <- fast_hetero(z$tree, x, method = "D", nodes = z$nodes, min_tips = 10,
    test = TRUE, nsim = 39, d_nsim = 79, keep_null = TRUE, progress = FALSE)
  seed <- .Random.seed
  set.seed(173)
  b <- fast_hetero(z$tree, 1-x, method = "D", nodes = z$nodes, min_tips = 10,
    test = TRUE, nsim = 39, d_nsim = 79, keep_null = TRUE, progress = FALSE)
  expect_equal(a$clades$estimate, b$clades$estimate, tolerance = 1e-12)
  expect_equal(a$heterogeneity, b$heterogeneity, tolerance = 1e-12)
  expect_equal(a$null, b$null)
  expect_identical(.Random.seed, seed)
})
