test_that("check_tree reports method-specific readiness without mutation", {
  testthat::skip_if_not_installed("ape")

  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  before <- serialize(tree, NULL)
  got <- check_tree(tree)
  after <- serialize(tree, NULL)

  expect_s3_class(got, "fastphylosig_tree_check")
  expect_true(got$ready)
  expect_equal(got$ready_by_signal,
               c(K = TRUE, lambda = TRUE, D = TRUE, Delta = TRUE))
  expect_equal(got$tree_summary$n_tip, 4L)
  expect_true(got$tree_summary$rooted)
  expect_true(got$tree_summary$ultrametric)
  expect_equal(got$tree_summary$n_polytomies, 0L)
  expect_equal(got$tree_summary$n_single_child, 0L)
  expect_equal(before, after)
  expect_true(all(c("code", "severity", "signal", "message") %in%
                    names(got$issues)))
})

test_that("polytomies are not blanket-rejected", {
  testthat::skip_if_not_installed("ape")

  tree <- ape::read.tree(text = "((a:1,b:1,c:1):1,d:2);")
  got <- check_tree(tree)

  expect_true(got$ready_by_signal[["K"]])
  expect_true(got$ready_by_signal[["lambda"]])
  expect_true(got$ready_by_signal[["D"]])
  expect_false(got$ready_by_signal[["Delta"]])
  expect_true(got$tree_summary$has_polytomies)
  expect_true(got$tree_summary$rooted)
  expect_true(any(got$issues$code == "delta_requires_binary_tree"))
  expect_true(any(got$issues$code == "polytomy_supported"))
})

test_that("zero, negative, near-zero, and invalid labels are diagnosed", {
  testthat::skip_if_not_installed("ape")

  zero <- ape::read.tree(text = "((a:0,b:1):1,(c:1,d:1):1);")
  z <- check_tree(zero)
  expect_false(z$ready_by_signal[["K"]])
  expect_false(z$ready_by_signal[["D"]])
  expect_false(z$ready_by_signal[["Delta"]])
  expect_false(z$ready_by_signal[["lambda"]])
  expect_equal(z$tree_summary$n_zero_terminal, 1L)
  expect_true(any(z$issues$code == "k_requires_positive_branches"))
  expect_true(any(z$issues$code == "lambda_zero_terminal_branch"))

  negative <- zero
  negative$edge.length[[2L]] <- -1
  n <- check_tree(negative, signal = c("K", "D"))
  expect_false(n$ready)
  expect_true(any(n$issues$code == "negative_branch_lengths"))
  expect_equal(names(n$ready_by_signal), c("K", "D"))

  near <- zero
  near$edge.length <- rep(1, nrow(near$edge))
  near$edge.length[[1L]] <- 1e-12
  near$tip.label[[2L]] <- ""
  near$tip.label[[3L]] <- near$tip.label[[1L]]
  q <- check_tree(near)
  expect_true(any(q$issues$code == "near_zero_branch"))
  expect_true(any(q$issues$code == "empty_tip_label"))
  expect_true(any(q$issues$code == "duplicate_tip_labels"))
})

test_that("single-child nodes block D and Delta but do not impose a global binary rule", {
  tree <- list(
    edge = matrix(c(3L, 4L, 4L, 1L, 3L, 2L), ncol = 2L, byrow = TRUE),
    edge.length = c(1, 1, 1),
    Nnode = 2L,
    tip.label = c("a", "b")
  )
  class(tree) <- "phylo"
  got <- check_tree(tree)
  expect_true(got$ready_by_signal[["K"]])
  expect_true(got$ready_by_signal[["lambda"]])
  expect_false(got$ready_by_signal[["D"]])
  expect_false(got$ready_by_signal[["Delta"]])
  expect_equal(got$tree_summary$n_single_child, 1L)
  expect_true(any(got$issues$code == "d_requires_no_single_child"))
})

test_that("mixed unary and polytomous topology cannot pass Delta by cancellation", {
  tree <- list(
    edge = matrix(c(
      5L, 6L, 5L, 1L, 5L, 2L, 6L, 7L, 7L, 3L, 7L, 4L
    ), ncol = 2L, byrow = TRUE),
    edge.length = rep(1, 6),
    Nnode = 3L,
    tip.label = c("a", "b", "c", "d")
  )
  class(tree) <- "phylo"
  got <- check_tree(tree)
  expect_true(got$tree_summary$has_polytomies)
  expect_true(got$tree_summary$has_single_child)
  expect_false(got$ready_by_signal[["Delta"]])
  expect_true(any(got$issues$code == "delta_requires_binary_tree"))
})

test_that("invalid inputs return a check object rather than mutating or stopping", {
  got <- check_tree(NULL, signal = "K")
  expect_s3_class(got, "fastphylosig_tree_check")
  expect_false(got$ready)
  expect_false(got$ready_by_signal[["K"]])
  expect_true(any(got$issues$code == "invalid_input"))
  expect_error(check_tree(NULL, signal = "not-a-signal"), "signal")

  tree <- ape::rtree(5)
  tree$Nnode <- NA_real_
  invalid_nnode <- check_tree(tree, signal = "K")
  expect_false(invalid_nnode$ready)
  expect_true(any(invalid_nnode$issues$code == "invalid_nnode"))
})

test_that("print groups repeated multi-signal diagnostics", {
  tree <- ape::unroot(ape::rtree(8))
  checked <- check_tree(tree)
  output <- capture.output(print(checked))
  rooted_lines <- output[grepl("conventional rooted flag", output, fixed = TRUE)]
  expect_length(rooted_lines, 1L)
  expect_match(rooted_lines, "K/lambda/D", fixed = TRUE)
})

test_that("unrooted encodings are distinguished from structural rootedness", {
  testthat::skip_if_not_installed("ape")

  rooted <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  unrooted <- ape::unroot(rooted)
  before <- serialize(unrooted, NULL)
  got <- check_tree(unrooted)
  after <- serialize(unrooted, NULL)

  expect_false(got$tree_summary$rooted)
  expect_true(got$tree_summary$structural_rooted)
  expect_true(any(got$issues$code == "conventional_unrooted"))
  expect_true(got$ready_by_signal[["K"]])
  expect_true(got$ready_by_signal[["lambda"]])
  # D can use the structural compatibility path when ape retains the
  # canonical n_tip + 1 root; otherwise the explicit canonical-root issue is
  # the method-specific blocking diagnostic.
  if (identical(got$tree_summary$root, got$tree_summary$n_tip + 1L)) {
    expect_true(got$ready_by_signal[["D"]])
  } else {
    expect_false(got$ready_by_signal[["D"]])
    expect_true(any(got$issues$code == "d_requires_canonical_root"))
  }
  expect_false(got$ready_by_signal[["Delta"]])
  expect_true(any(got$issues$code == "delta_requires_conventional_root"))
  expect_equal(before, after)
})

test_that("non-ultrametric trees remain ready for every supported method", {
  testthat::skip_if_not_installed("ape")

  tree <- ape::read.tree(text = "((a:1,b:1):1,c:3);")
  before <- serialize(tree, NULL)
  got <- check_tree(tree)
  after <- serialize(tree, NULL)

  expect_false(got$tree_summary$ultrametric)
  expect_true(all(got$ready_by_signal))
  expect_true(any(got$issues$code == "nonultrametric_tree"))
  expect_true(all(got$issues$severity[got$issues$code ==
                                        "nonultrametric_tree"] == "INFO"))
  expect_equal(before, after)
})

test_that("missing branch lengths block each method explicitly", {
  testthat::skip_if_not_installed("ape")

  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  tree$edge.length <- NULL
  before <- serialize(tree, NULL)
  got <- check_tree(tree)
  after <- serialize(tree, NULL)

  expect_false(got$ready)
  expect_false(got$ready_by_signal[["K"]])
  expect_false(got$ready_by_signal[["lambda"]])
  expect_false(got$ready_by_signal[["D"]])
  expect_false(got$ready_by_signal[["Delta"]])
  expect_true(any(got$issues$code == "missing_branch_lengths"))
  expect_true(all(c("K", "lambda", "D", "Delta") %in%
                    got$issues$signal))
  expect_equal(before, after)
})

test_that("zero internal branches are conditional for lambda but block other kernels", {
  testthat::skip_if_not_installed("ape")

  tree <- ape::read.tree(text = "((a:1,b:1):0,(c:1,d:1):1);")
  before <- serialize(tree, NULL)
  got <- check_tree(tree)
  after <- serialize(tree, NULL)

  expect_equal(got$tree_summary$n_zero_internal, 1L)
  expect_false(got$ready_by_signal[["K"]])
  expect_true(got$ready_by_signal[["lambda"]])
  expect_false(got$ready_by_signal[["D"]])
  expect_false(got$ready_by_signal[["Delta"]])
  expect_true(any(got$issues$code == "k_requires_positive_branches" &
                    got$issues$signal == "K"))
  expect_true(any(got$issues$code == "lambda_zero_internal_branch" &
                    got$issues$signal == "lambda"))
  expect_true(any(got$issues$code == "d_requires_positive_branches" &
                    got$issues$signal == "D"))
  expect_true(any(got$issues$code == "delta_requires_positive_branches" &
                    got$issues$signal == "Delta"))
  expect_equal(before, after)
})

test_that("two-tip trees are supported and retain method-specific readiness", {
  testthat::skip_if_not_installed("ape")

  tree <- ape::read.tree(text = "(a:1,b:2);")
  before <- serialize(tree, NULL)
  got <- check_tree(tree)
  after <- serialize(tree, NULL)

  expect_equal(got$tree_summary$n_tip, 2L)
  expect_equal(got$tree_summary$n_internal, 1L)
  expect_true(got$tree_summary$fully_dichotomous)
  expect_true(all(got$ready_by_signal))
  expect_equal(before, after)
})

test_that("pectinate and strongly unbalanced trees are not rejected as unresolved", {
  testthat::skip_if_not_installed("ape")

  pectinate <- ape::read.tree(
    text = "((((a:1,b:1):1,c:1):1,d:1):1,e:1);"
  )
  strongly_unbalanced <- ape::read.tree(
    text = "(((((a:1,b:1):1,c:1):1,d:1):1,e:1):1,f:1);"
  )

  for (tree in list(pectinate = pectinate,
                    strongly_unbalanced = strongly_unbalanced)) {
    before <- serialize(tree, NULL)
    got <- check_tree(tree)
    after <- serialize(tree, NULL)
    expect_true(got$tree_summary$rooted)
    expect_equal(got$tree_summary$n_polytomies, 0L)
    expect_equal(got$tree_summary$n_single_child, 0L)
    expect_true(got$tree_summary$fully_dichotomous)
    expect_true(all(got$ready_by_signal))
    expect_equal(before, after)
  }
})
