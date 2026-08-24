test_that("the internal tree requirement registry covers all production methods", {
  registry <- fastphylosig:::.tree_requirement_registry()
  expect_equal(names(registry), c("K", "lambda", "D", "Delta", "ACE"))
  expect_true(all(vapply(registry, function(x) isTRUE(x$min_tips >= 2L), logical(1))))
  expect_true(all(vapply(registry, function(x) isTRUE(x$finite_branch_lengths), logical(1))))
  expect_true(all(vapply(registry, function(x) isTRUE(x$structural_root), logical(1))))
  expect_false(fastphylosig:::.tree_requirement("K")$conventional_root)
  expect_true(fastphylosig:::.tree_requirement("Delta")$conventional_root)
  expect_equal(fastphylosig:::.tree_requirement("lambda")$zero_internal,
               "warning")
  expect_true(fastphylosig:::.tree_requirement("ACE")$binary)
})

test_that("check_tree still uses shared issue metadata", {
  testthat::skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:0,b:1):1,(c:1,d:1):1);")
  checked <- check_tree(tree, signal = "K")
  expect_true(all(c("check", "problem", "action", "auto_fixable") %in%
                    names(checked$issues)))
  row <- checked$issues[checked$issues$code == "k_requires_positive_branches", , drop = FALSE]
  expect_equal(row$check[[1L]], "branch lengths")
  expect_false(row$auto_fixable[[1L]])
})
