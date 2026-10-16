test_that("explicit large crowns override automatic bounds but retain observed minima", {
  set.seed(1512)
  tree <- ape::rtree(600)
  x <- setNames(rnorm(600), tree$tip.label)
  root <- 601L
  fit <- fast_hetero(tree, x, nodes = root, min_tips = 500,
                     progress = FALSE)
  expect_identical(fit$selection$nodes, root)
  expect_equal(fit$clades$n_tip_total, 600)
  expect_equal(fit$clades$status, "ok")
  expect_true(is.finite(fit$clades$estimate))

  x[seq_len(101)] <- NA_real_
  sparse <- fast_hetero(tree, x, nodes = root, min_tips = 500,
                        progress = FALSE)
  expect_equal(sparse$clades$status, "insufficient_data")
  expect_true(is.na(sparse$clades$estimate))

  expect_error(fast_hetero(tree, x, min_tips = 500, progress = FALSE), "max_tips")
  expect_error(fast_hetero(tree, x, min_tips = 20, max_tips = 10,
                          progress = FALSE), "max_tips")
  expect_error(fast_hetero(tree, x, nodes = root, min_tips = 2,
                          progress = FALSE), "min_tips")
  expect_error(fast_hetero(tree, x, nodes = root, max_tips = 2,
                          progress = FALSE), "max_tips")
  expect_error(fast_hetero(tree, x, nodes = root, max_tips = NA_real_,
                          progress = FALSE), "max_tips")
  expect_error(fast_hetero(tree, x, nodes = root, max_tips = 3.5,
                          progress = FALSE), "max_tips")
  expect_error(fast_hetero(tree, x, nodes = root, min_tips = 500,
                          max_total_tips = 599, progress = FALSE), "max_total_tips")
})

test_that("native explicit plans ignore only the automatic bounds relationship", {
  set.seed(1513)
  tree <- ape::rtree(600)
  core <- .heterogeneity_compile(tree, "K")
  plan <- heterogeneity_plan_cpp(core, tree$tip.label, 601L, 500L, 200L, 1e6)
  expect_identical(plan$nodes, 601L)
  expect_equal(plan$n_tip, 600L)
  expect_error(heterogeneity_plan_cpp(core, tree$tip.label, NULL, 500L, 200L, 1e6),
               "invalid clade plan dimensions")
  expect_error(heterogeneity_plan_cpp(core, tree$tip.label, 601L, 2L, 200L, 1e6),
               "invalid clade plan dimensions")
  expect_error(heterogeneity_plan_cpp(core, tree$tip.label, 601L, 500L, 2L, 1e6),
               "invalid clade plan dimensions")
  expect_error(heterogeneity_plan_cpp(core, tree$tip.label, 601L, 500L, 200L, 599),
               "max_total_tips")
})
