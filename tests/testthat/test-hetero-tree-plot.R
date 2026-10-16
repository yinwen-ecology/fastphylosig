test_that("fast_hetero is the sole heterogeneity API and lambda can be calibrated", {
  expect_false("fast_signal_heterogeneity" %in% getNamespaceExports("fastphylosig"))
  expect_false(exists("fast_signal_heterogeneity", envir = asNamespace("fastphylosig"), inherits = FALSE))
  z <- .hetero_fixture()
  lambda <- fast_hetero(z$tree, z$X, nodes = z$nodes, min_tips = 10L,
    method = "lambda", progress = FALSE)
  expect_equal(lambda$heterogeneity$dispersion, sd(lambda$clades$estimate))
  expect_identical(lambda$heterogeneity$scale, "SD(lambda)")
  expect_true(all(is.na(lambda$heterogeneity$P_BM)))
  calibrated <- fast_hetero(z$tree, z$X, nodes = z$nodes, min_tips = 10L,
    method = "lambda", test = TRUE, nsim = 9, progress = FALSE)
  expect_identical(calibrated$heterogeneity$status, "ok")
})

test_that("native annotation matches an independent ancestor walk and source IDs", {
  tree <- ape::stree(64L, "balanced")
  tree$edge.length <- seq(0.1, 2, length.out = nrow(tree$edge))
  nt <- length(tree$tip.label)
  ids <- c(seq_len(nt), nt + rev(seq_len(tree$Nnode)))
  tree$edge[] <- ids[tree$edge]
  tree$edge <- tree$edge[nrow(tree$edge):1L, , drop = FALSE]
  tree$edge.length <- rev(tree$edge.length)
  core <- compile_tree_cpp(tree$edge, tree$edge.length, nt)
  index <- .hetero_index_reference(core)
  # Deliberately reverse nested selection order.
  nodes <- rev(head(core$preorder[core$preorder > nt], 4L))
  tips <- lapply(nodes, function(node) {
    ref <- .hetero_subtree_reference(tree, core, index, node)
    attr(ref, "source_tip_index")
  })
  mapped <- heterogeneity_tree_map_cpp(core, nodes, unlist(tips, use.names = FALSE),
    as.integer(c(0, cumsum(lengths(tips)))))
  reference <- integer(core$n_total)
  for (node in seq_len(core$n_total)) {
    ancestor <- node
    while (ancestor > 0L) {
      selected <- match(ancestor, nodes, nomatch = 0L)
      if (selected) { reference[[node]] <- selected; break }
      ancestor <- core$parent[[ancestor]]
    }
  }
  expect_identical(mapped$node_owner[mapped$source_to_plot], reference)
  inverse <- integer(core$n_total)
  inverse[mapped$source_to_plot] <- seq_len(core$n_total)
  child <- inverse[mapped$edge[, 2L]]
  parent <- inverse[mapped$edge[, 1L]]
  expect_identical(parent, core$parent[child])
  expect_identical(mapped$edge_length, core$branch_length_by_node[child])
  expected_edge <- reference[child]
  stems <- child %in% nodes
  expected_edge[stems] <- reference[parent[stems]]
  expect_identical(mapped$edge_owner, expected_edge)
  expect_identical(mapped$source_to_plot[[core$root]], nt + 1L)
  expect_error(heterogeneity_tree_map_cpp(core, nodes, rep(NA_integer_, sum(lengths(tips))),
    as.integer(c(0, cumsum(lengths(tips))))), "does not match")
})

test_that("tree plots support K/lambda, unavailable crowns and prepared trees", {
  z <- .hetero_fixture()
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  for (method in c("K", "lambda")) {
    h <- fast_hetero(z$tree, z$X, nodes = z$nodes, min_tips = 10L,
      method = method, progress = FALSE)
    set.seed(129)
    seed <- .Random.seed
    p <- plot(h, tree = z$tree, node_labels = TRUE)
    expect_identical(.Random.seed, seed)
    expect_identical(p$color_scale, if (method == "K") "log" else "raw")
    expect_length(p$edge_colors, nrow(z$tree$edge))
    expect_identical(p$clades, h$clades)
    expect_equal(sum(p$edge_owner == 0L), 2L) # The two crown stems.
    fan <- plot(h, tree = prepare_tree(z$tree), type = "fan", show_legend = FALSE)
    expect_identical(fan$edge_colors, p$edge_colors)
    expect_error(plot(h, view = "tree"), "Supply the original tree")
    wrong <- z$tree
    wrong$tip.label[c(1, 64)] <- wrong$tip.label[c(64, 1)]
    expect_error(plot(h, tree = wrong), "does not match")
  }
  h <- fast_hetero(z$tree, z$X * NA_real_, nodes = z$nodes, min_tips = 10L,
    progress = FALSE)
  p <- plot(h, tree = z$tree)
  expect_true(all(p$edge_colors == "grey80"))
  expect_length(p$legend_values, 0L)
  expect_error(plot(h, tree = z$tree, palette = "red"), "at least two")
  expect_error(plot(h, tree = z$tree, na_col = NA_character_), "one colour")
  expect_error(plot(h, tree = z$tree, show_legend = NA), "TRUE or FALSE")
})

test_that("unavailable nested crowns override enclosing estimates on tree plots", {
  z <- .hetero_fixture()
  core <- compile_tree_cpp(z$tree$edge, z$tree$edge.length, 64L)
  nodes <- c(z$nodes[[1L]], core$root)
  h <- fast_hetero(z$tree, z$X, nodes = nodes, min_tips = 10L, progress = FALSE)
  # Simulate an unestimable inner crown, with a valid enclosing estimate.
  h$clades$status[[1L]] <- "undefined"
  h$clades$estimate[[1L]] <- NA_real_
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  p <- plot(h, tree = z$tree, node_labels = TRUE)
  expect_true(all(p$edge_colors[p$edge_owner == 1L] == "grey80"))
  expect_true(all(p$edge_colors[p$edge_owner == 2L] != "grey80"))
})

test_that("overall signal tests reuse public engines on all observed tree tips", {
  z <- .hetero_fixture()
  X <- cbind(x = z$X[, 1L], missing = z$X[, 1L])
  rownames(X) <- rownames(z$X)
  X[c(1, 7, 40), 2L] <- NA_real_
  for (method in c("K", "lambda")) {
    set.seed(318)
    h <- fast_hetero(z$tree, X, nodes = z$nodes[[1L]], min_tips = 10L,
      method = method, overall_test = TRUE, overall_nsim = 19L, progress = FALSE)
    seed <- .Random.seed
    set.seed(318)
    ref <- if (method == "K") fast_k(z$tree, X, test = TRUE, nsim = 19L,
      simulation_chunk = 32L, verbose = FALSE, progress = FALSE) else
      fast_lambda(z$tree, X, test = TRUE, lambda_profile = FALSE,
        verbose = FALSE, progress = FALSE)
    expect_identical(h$overall$estimate, ref$estimate)
    expect_identical(h$overall$P, ref$P_fast)
    expect_identical(h$overall$n_species, c(64L, 61L))
    expect_identical(h$overall$status, ref$status)
    expect_identical(.Random.seed, seed)
    if (method == "lambda") expect_identical(h$overall$LR, ref$LR_fast)
    file <- tempfile(fileext = ".pdf")
    grDevices::pdf(file)
    p <- plot(h, tree = z$tree, trait = "missing", type = "fan")
    expect_identical(.Random.seed, seed)
    expect_match(p$summary_text[[1L]], "n = 61", fixed = TRUE)
    expect_match(p$summary_text[[2L]], if (method == "K") "Randomization P" else "LRT P")
    if (method == "lambda")
      expect_match(p$summary_text[[4L]], "Shared-lambda heterogeneity test: not performed", fixed = TRUE)
    expect_length(plot(h, tree = z$tree, trait = "x", show_summary = FALSE)$summary_text, 0L)
    expect_error(plot(h, tree = z$tree, trait = "x", show_summary = NA), "TRUE or FALSE")
    grDevices::dev.off()
    unlink(file)
  }
})

test_that("overall analysis is optional and cannot change the BM reference", {
  z <- .hetero_fixture()
  call <- list(tree = z$tree, data = z$X, nodes = z$nodes, min_tips = 10L,
    test = TRUE, nsim = 19L, keep_null = TRUE, progress = FALSE)
  set.seed(910)
  local <- do.call(fast_hetero, call)
  expect_null(local$overall)
  set.seed(910)
  full <- do.call(fast_hetero, c(call, list(overall_test = TRUE, overall_nsim = 7L)))
  unchanged <- setdiff(names(local$clades), c("overall_estimate", "difference_overall", "log_ratio_overall"))
  expect_identical(full$clades[, unchanged], local$clades[, unchanged])
  unchanged <- setdiff(names(local$heterogeneity), "local_global_deviation")
  expect_identical(full$heterogeneity[, unchanged], local$heterogeneity[, unchanged])
  expect_identical(full$null, local$null)
  expect_identical(full$overall$nsim_requested, 7L)
  h <- fast_hetero(z$tree, z$X, nodes = z$nodes, min_tips = 10L,
    overall = TRUE, progress = FALSE)
  expect_true(is.finite(h$overall$estimate))
  expect_true(is.na(h$overall$P))
  expect_identical(h$overall$test_method, "none")
  expect_match(.heterogeneity_summary_text(h, colnames(z$X)[[1L]])[[2L]], "not performed")
  expect_match(.heterogeneity_summary_text(h, colnames(z$X)[[1L]])[[4L]], "not performed")
  unavailable <- fast_hetero(z$tree, z$X * NA_real_, nodes = z$nodes,
    min_tips = 10L, overall_test = TRUE, progress = FALSE)
  expect_true(is.na(unavailable$overall$estimate))
  expect_true(is.na(unavailable$overall$P))
  expect_match(.heterogeneity_summary_text(unavailable, colnames(z$X)[[1L]])[[1L]], "unavailable")
  expect_error(fast_hetero(z$tree, z$X, overall = NA), "TRUE or FALSE")
  expect_error(fast_hetero(z$tree, z$X, overall_test = NA), "TRUE or FALSE")
  expect_error(fast_hetero(z$tree, z$X, overall_test = TRUE, overall_nsim = 0), "overall_nsim")
})
