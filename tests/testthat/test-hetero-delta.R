.hetero_delta_control <- list(model = "ER", mcmc_sim = 1200L, thin = 5L, burn = 100L)

test_that("Delta crowns reuse the specialist estimator and serial RNG", {
  z <- .hetero_fixture()
  X <- data.frame(three = rep(c("A", "B", "C"), length.out = 64),
    four = rep(c("w", "x", "y", "z"), 16), row.names = z$tree$tip.label)
  X[c(1, 40), 1] <- NA
  set.seed(121)
  h <- suppressWarnings(fast_hetero(z$tree, X, method = "Delta", nodes = z$nodes,
    min_tips = 20, delta_control = .hetero_delta_control, progress = FALSE))
  seed <- .Random.seed
  set.seed(121)
  expected <- unlist(lapply(z$nodes, function(node) {
    sub <- ape::extract.clade(z$tree, node)
    unlist(lapply(X, function(y) {
      names(y) <- rownames(X)
      # Public Delta drops unobserved states and prunes each NA mask.
      suppressWarnings(do.call(fast_delta, c(list(tree = sub, x = y,
        verbose = FALSE, progress = FALSE), .hetero_delta_control)))$estimate
    }), use.names = FALSE)
  }), use.names = FALSE)
  expect_equal(h$clades$estimate, expected, tolerance = 1e-12)
  expect_identical(.Random.seed, seed)
  expect_identical(h$clades$n_species, c(31L, 32L, 31L, 32L))
  expect_identical(h$clades$n_states, c(3L, 4L, 3L, 4L))
  expect_identical(names(h$clades$state_counts[[1L]]), c("A", "B", "C"))
  expect_equal(sum(h$clades$state_counts[[1L]]), 31)
  expect_equal(sum(h$clades$state_frequencies[[1L]]), 1)
  expect_equal(h$heterogeneity$dispersion,
    vapply(c("three", "four"), function(nm) sd(log(h$clades$estimate[h$clades$trait == nm])), numeric(1)),
    ignore_attr = TRUE)
  expect_true(all(is.na(h$heterogeneity$P_heterogeneity)))
})

test_that("random-tree canonicalization and interleaved NA groups preserve Delta RNG", {
  set.seed(619); tree <- ape::rtree(150)
  X <- data.frame(a = sample(letters[1:4], 150, TRUE),
    b = sample(letters[1:3], 150, TRUE), c = sample(letters[1:4], 150, TRUE),
    row.names = tree$tip.label)
  X[c(2, 12, 48, 83, 99), 2] <- NA
  set.seed(459)
  h <- suppressWarnings(fast_hetero(tree, X, method = "Delta", min_tips = 15,
    max_tips = 45, delta_control = .hetero_delta_control, progress = FALSE))
  seed <- .Random.seed
  set.seed(459)
  expected <- unlist(lapply(h$selection$nodes, function(node) {
    sub <- ape::extract.clade(tree, node)
    # a and c share a mask, but b must still run before c in the RNG stream.
    ref <- suppressWarnings(do.call(fast_delta, c(list(tree = sub, x = X,
      verbose = FALSE, progress = FALSE), .hetero_delta_control)))
    ref$estimate
  }), use.names = FALSE)
  expect_equal(h$clades$estimate, expected, tolerance = 1e-10)
  expect_identical(.Random.seed, seed)
})

test_that("overall Delta, contrasts and plots describe the same statistic", {
  z <- .hetero_fixture()
  x <- setNames(rep(c(1, 2, 3, 4), 16), z$tree$tip.label)
  args <- list(tree = z$tree, data = x, method = "Delta", nodes = z$nodes,
    min_tips = 20, delta_control = .hetero_delta_control, progress = FALSE)
  set.seed(81)
  local <- suppressWarnings(do.call(fast_hetero, args))
  ref <- suppressWarnings(do.call(fast_delta, c(list(tree = z$tree, x = x,
    verbose = FALSE, progress = FALSE), .hetero_delta_control)))
  set.seed(81)
  h <- suppressWarnings(do.call(fast_hetero, c(args, list(overall = TRUE))))
  expect_equal(h$overall$estimate, ref$estimate, tolerance = 1e-12)
  expect_equal(h$clades$estimate, local$clades$estimate)
  ratios <- log(h$clades$estimate) - log(h$overall$estimate)
  expect_equal(h$clades$log_ratio_overall, ratios)
  expect_equal(h$heterogeneity$local_global_deviation, sqrt(mean(ratios^2)))
  expect_equal(h$clades$difference_overall, h$clades$estimate - h$overall$estimate)
  file <- tempfile(fileext = ".pdf"); grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  seed <- .Random.seed
  p <- plot(h, tree = z$tree, node_labels = TRUE)
  expect_identical(.Random.seed, seed)
  expect_identical(p$color_scale, "log")
  expect_match(p$summary_text[1], "Overall Delta")
  expect_true(any(grepl("RMS", p$summary_text)))
  expect_identical(plot(h, view = "clades"), h$clades)
})

test_that("native multistate simulation matches analytic ER and chunk RNG", {
  set.seed(20); tree <- ape::rtree(35)
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 35L)
  for (k in c(3L, 4L, 6L)) {
    index <- matrix(1L, k, k); diag(index) <- 0L
    workspace <- heterogeneity_mk_prepare_cpp(core, index, .4)
    set.seed(48); native <- heterogeneity_mk_simulate_cpp(workspace, 11L); seed <- .Random.seed
    set.seed(48); expected <- .hetero_er_states(core, k, .4, 11L)
    expect_identical(native, expected)
    expect_identical(.Random.seed, seed)
    set.seed(48)
    split <- cbind(heterogeneity_mk_simulate_cpp(workspace, 3L),
      heterogeneity_mk_simulate_cpp(workspace, 8L))
    expect_identical(split, native)
    expect_identical(.Random.seed, seed)
  }
  expect_error(heterogeneity_mk_simulate_cpp(NULL, 1L), "invalid")
})

test_that("ARD simulation matches an independent uniformization oracle", {
  set.seed(34); tree <- ape::rtree(18)
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 18L)
  for (k in c(3L, 4L, 6L)) {
    rate <- .fast_ace_rate_index("ARD", k)
    rates <- seq(.05, .8, length.out = k*(k-1L))
    Q <- matrix(0, k, k)
    off <- row(Q) != col(Q); Q[off] <- rates[rate$rate_index[off]]
    diag(Q) <- -rowSums(Q)
    mu <- max(-diag(Q)); jump <- diag(k) + Q/mu
    transitions <- lapply(core$branch_length_by_node, function(t) {
      term <- diag(k); P <- term
      for (n in 1:180) { term <- term %*% jump * (mu*t/n); P <- P + term }
      P * exp(-mu*t)
    })
    workspace <- heterogeneity_mk_prepare_cpp(core, rate$rate_index, rates)
    set.seed(94); native <- heterogeneity_mk_simulate_cpp(workspace, 13L); seed <- .Random.seed
    set.seed(94); expected <- matrix(0L, core$n_tip, 13L)
    for (b in 1:13) {
      states <- integer(core$n_total)
      states[core$root] <- as.integer(floor(runif(1, 0, k))) + 1L
      for (node in core$preorder[-1L]) {
        probabilities <- transitions[[node]][states[core$parent[node]], ]
        states[node] <- which(runif(1) < cumsum(probabilities))[1L]
      }
      expected[,b] <- states[seq_len(core$n_tip)]
    }
    expect_identical(native, expected)
    expect_identical(.Random.seed, seed)
  }
})

test_that("Delta bootstrap matches whole-tree R simulation and specialist refits", {
  z <- .hetero_fixture(); k <- 4L
  x <- setNames(rep(seq_len(k), 16L), z$tree$tip.label)
  set.seed(171)
  h <- suppressWarnings(fast_hetero(z$tree, x, method = "Delta", nodes = z$nodes,
    min_tips = 20, test = TRUE, nsim = 3L, keep_null = TRUE,
    delta_control = .hetero_delta_control, progress = FALSE))
  seed <- .Random.seed
  specialist <- function(tree, x) suppressWarnings(do.call(fast_delta,
    c(list(tree = tree, x = x, verbose = FALSE, progress = FALSE), .hetero_delta_control)))$estimate
  set.seed(171)
  subtrees <- lapply(z$nodes, function(node) ape::extract.clade(z$tree, node))
  local <- vapply(subtrees, specialist, numeric(1), x = x)
  global <- specialist(z$tree, x)
  expect_equal(h$heterogeneity$statistic, sqrt(mean((log(local) - log(global))^2)))
  core <- compile_tree_cpp(z$tree$edge, z$tree$edge.length, 64L)
  common <- h$heterogeneity$reference_model[[1L]]$rates
  expected <- numeric(3L)
  for (b in 1:3) {
    y <- setNames(.hetero_er_states(core, k, common, 1L)[,1L], z$tree$tip.label)
    global <- specialist(z$tree, y)
    local <- vapply(subtrees, specialist, numeric(1), x = y)
    expected[b] <- sqrt(mean((log(local) - log(global))^2))
  }
  expect_equal(h$null[[1L]], expected, tolerance = 1e-9)
  expect_identical(.Random.seed, seed)
})

test_that("Delta reference tests retain all draws and test signal deviations", {
  z <- .hetero_fixture()
  for (k in c(3L, 4L, 6L)) {
    x <- setNames(rep(seq_len(k), length.out = 64), z$tree$tip.label)
    if (k == 6L) {
      core <- compile_tree_cpp(z$tree$edge, z$tree$edge.length, 64L)
      set.seed(381)
      x[] <- .hetero_er_states(core, k, .15, 1L)[,1L]
    }
    x[] <- LETTERS[as.integer(x)]
    args <- list(tree = z$tree, data = x, method = "Delta", nodes = z$nodes,
      min_tips = 20, test = TRUE, nsim = 3L, keep_null = TRUE,
      delta_control = .hetero_delta_control, progress = FALSE)
    set.seed(213)
    h <- suppressWarnings(do.call(fast_hetero, args)); seed <- .Random.seed
    expect_identical(h$heterogeneity$status, "ok")
    expect_identical(h$heterogeneity$nsim_successful, 3L)
    expect_equal(h$heterogeneity$statistic, h$heterogeneity$local_global_deviation)
    expect_equal(h$heterogeneity$P_heterogeneity,
      (1+sum(h$null[[1]] >= h$heterogeneity$statistic))/4)
    expect_match(h$heterogeneity$test_statistic, "Delta_local / Delta_overall", fixed = TRUE)
    expect_equal(length(h$heterogeneity$reference_model[[1]]$states), k)
    expect_identical(h$heterogeneity$reference_model[[1]]$states, LETTERS[seq_len(k)])
    set.seed(213)
    again <- suppressWarnings(do.call(fast_hetero, c(args, list(simulation_chunk = 1L))))
    expect_identical(again$null, h$null)
    expect_identical(.Random.seed, seed)
  }
  tail <- heterogeneity_bootstrap_summary_cpp(c(1, NA, 3), 2, FALSE)
  expect_equal(tail$nsim_successful, 2L)
  expect_true(is.na(tail$P_heterogeneity))
  expect_true(is.na(tail$null_mean))
})

test_that("categorical state support, missingness and invalid controls are explicit", {
  z <- .hetero_fixture()
  x <- setNames(c(rep(1, 32), rep(c(2, 3), 16)), z$tree$tip.label)
  h <- suppressWarnings(fast_hetero(z$tree, x, method = "Δ", nodes = z$nodes,
    min_tips = 20, test = TRUE, nsim = 2L, delta_control = .hetero_delta_control, progress = FALSE))
  expect_identical(h$clades$n_states, c(1L, 2L))
  expect_identical(h$clades$status, c("undefined", "ok"))
  expect_true(is.na(h$clades$estimate[1]))
  expect_identical(h$heterogeneity$status, "unavailable_clades")
  expect_true(is.na(h$heterogeneity$P_heterogeneity))
  sparse <- x; sparse[1:25] <- NA
  h <- suppressWarnings(fast_hetero(z$tree, sparse, method = "Delta", nodes = z$nodes,
    min_tips = 20, delta_control = .hetero_delta_control, progress = FALSE))
  expect_identical(h$clades$status[1], "insufficient_data")
  expect_equal(h$clades$n_states[1], 1L)
  expect_error(fast_hetero(z$tree, x, method = "Delta", delta_control = list(thin = 0), progress = FALSE), "thin")
  expect_error(fast_hetero(z$tree, x, method = "Delta", delta_control = list(unknown = 1), progress = FALSE), "Unknown")
  expect_error(fast_hetero(z$tree, x, method = "Delta", delta_control = list(NULL), progress = FALSE), "named")
})

test_that("ARD reference handles missing tips on the original root geometry", {
  z <- .hetero_fixture()
  x <- setNames(rep(c(1, 2, 3), length.out = 64), z$tree$tip.label)
  x[c(1, 8, 35)] <- NA
  control <- .hetero_delta_control; control$model <- "ARD"
  set.seed(791)
  h <- suppressWarnings(fast_hetero(z$tree, x, method = "Delta", nodes = z$nodes,
    min_tips = 20, test = TRUE, nsim = 2L, keep_null = TRUE,
    delta_control = control, progress = FALSE))
  model <- h$heterogeneity$reference_model[[1L]]
  expect_identical(model$model, "ARD")
  expect_identical(model$root, "uniform")
  expect_length(model$rates, 6L)
  workspace <- .fast_ace_workspace(.fast_ace_tree_workspace(z$tree), "ARD", 3L)
  expected <- suppressWarnings(.fast_ace_fit_workspace(as.integer(factor(x)),
    levels(factor(x)), workspace, ip = .1, CI = FALSE, estimate_se = FALSE))
  expect_equal(model$rates, expected$rates)
  expect_equal(model$loglik, expected$loglik)
  expect_identical(h$clades$n_species, c(30L, 31L))
  # Finite and failed replicate statuses both obey complete-draw accounting.
  if (all(is.finite(h$null[[1L]]))) expect_true(is.finite(h$heterogeneity$P_heterogeneity)) else
    expect_true(is.na(h$heterogeneity$P_heterogeneity))
})

test_that("optional local signal tests are separate from heterogeneity", {
  z <- .hetero_fixture()
  x <- setNames(rep(c(1, 2, 3, 4), 16), z$tree$tip.label)
  h <- suppressWarnings(fast_hetero(z$tree, x, method = "Delta", nodes = z$nodes,
    min_tips = 20, clade_test = TRUE, clade_nsim = 2L,
    overall_test = TRUE, overall_nsim = 2L, delta_control = .hetero_delta_control,
    progress = FALSE))
  expect_true(all(is.finite(h$clades$P_signal)))
  expect_true(is.finite(h$overall$P))
  expect_true(is.na(h$heterogeneity$P_heterogeneity))
  expect_match(h$overall$test_method, "randomization")
  h$overall$P <- 0; h$overall$nsim_successful <- 2L
  expect_match(.heterogeneity_summary_text(h, "x")[2L], "0/2", fixed = TRUE)
  for (method in c("K", "lambda", "D")) {
    y <- if (method == "D") setNames(rep(c(0, 1), 32), z$tree$tip.label) else z$X
    h <- fast_hetero(z$tree, y, method = method, nodes = z$nodes, min_tips = 20,
      clade_test = TRUE, clade_nsim = 3, d_nsim = 41, overall = TRUE, progress = FALSE)
    expect_true(all(is.finite(h$clades$P_signal)))
    expect_true(is.finite(h$heterogeneity$local_global_deviation))
    if (method == "D") expect_true(all(is.na(h$clades$log_ratio_overall)))
  }
})
