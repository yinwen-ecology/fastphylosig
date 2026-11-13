.review_capture_messages <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    eval(substitute(expr), parent.frame()),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, messages = messages)
}

test_that("unpruned matching reports only removals that actually happened", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  x <- c(a = 0, b = 1, c = 0)
  captured <- .review_capture_messages(match_tree_data(tree, x, prune = FALSE))
  fit <- captured$value
  expect_setequal(fit$tree$tip.label, tree$tip.label)
  expect_equal(rownames(fit$data), names(x))
  expect_identical(fit$tree_tips_removed, character())
  expect_equal(fit$report$removed_tree_tips, 0L)
  expect_length(captured$messages, 0L)

  extra <- .review_capture_messages(match_tree_data(
    tree, c(x, outside = 1), prune = FALSE
  ))
  expect_identical(extra$value$tree_tips_removed, character())
  expect_identical(extra$value$data_rows_removed, "outside")
  expect_equal(extra$value$report$removed_tree_tips, 0L)
  expect_equal(extra$value$report$removed_data_rows, 1L)
  expect_length(extra$messages, 1L)
  expect_match(extra$messages, "removed 0 tree tips and 1 data rows", fixed = TRUE)
})

test_that("pruned matching still reports the removed tree tips", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  captured <- .review_capture_messages(match_tree_data(
    tree, c(a = 0, b = 1, c = 0, outside = 1), prune = TRUE
  ))
  fit <- captured$value
  expect_setequal(fit$tree$tip.label, c("a", "b", "c"))
  expect_identical(fit$tree_tips_removed, "d")
  expect_identical(fit$data_rows_removed, "outside")
  expect_equal(fit$report$removed_tree_tips, 1L)
  expect_equal(fit$report$removed_data_rows, 1L)
  expect_length(captured$messages, 1L)
})

test_that("D emits just one matching notice for partial matches", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  captured <- .review_capture_messages(fast_d(
    tree, c(a = 0, b = 1, c = 0, outside = 1),
    test = FALSE, nsim = 2L, verbose = TRUE, progress = TRUE
  ))
  notices <- captured$messages[grepl("^Matched [0-9]+ species;", captured$messages)]
  expect_length(notices, 1L)
  expect_match(notices[[1L]], "removed 1 tree tips and 1 data rows", fixed = TRUE)
  expect_equal(attr(captured$value, "match_report")$removed_tree_tips, 1L)
  quiet <- .review_capture_messages(fast_d(
    tree, c(a = 0, b = 1, c = 0, outside = 1),
    test = FALSE, nsim = 2L, verbose = TRUE, progress = FALSE
  ))
  expect_length(quiet$messages, 0L)
})
