test_that("integer-adjacency inspection preserves the baseline object", {
  testthat::skip_if_not_installed("ape")

  fixtures <- .stage2b2b_fixture_set()
  signals <- c("K", "lambda", "D", "Delta")
  for (fixture_name in names(fixtures)) {
    tree <- fixtures[[fixture_name]]
    old_tree <- .stage2b2b_clone(tree)
    current_tree <- .stage2b2b_clone(tree)
    old_before <- serialize(old_tree, NULL, version = 2L)
    current_before <- serialize(current_tree, NULL, version = 2L)
    old <- .stage2b2b_capture(
      .stage2b2b_old_inspect_tree_core, old_tree, signals
    )
    current <- .stage2b2b_capture(
      .stage2b2b_current_inspect_tree_core, current_tree, signals
    )

    .stage2b2b_compare_capture(
      old, current, paste("fixture", fixture_name)
    )
    expect_identical(
      serialize(old_tree, NULL, version = 2L), old_before,
      info = paste(fixture_name, "old oracle input mutation")
    )
    expect_identical(
      serialize(current_tree, NULL, version = 2L), current_before,
      info = paste(fixture_name, "candidate input mutation")
    )
  }
})

test_that("inspection preserves signal selection and condition classes", {
  testthat::skip_if_not_installed("ape")

  fixtures <- .stage2b2b_fixture_set()
  selected <- list(
    k_only = "K",
    lambda_d = c("lambda", "D"),
    delta_k = c("Delta", "K")
  )
  for (fixture_name in names(fixtures)) {
    for (selection_name in names(selected)) {
      tree <- fixtures[[fixture_name]]
      old <- .stage2b2b_capture(
        .stage2b2b_old_inspect_tree_core, .stage2b2b_clone(tree),
        selected[[selection_name]]
      )
      current <- .stage2b2b_capture(
        .stage2b2b_current_inspect_tree_core, .stage2b2b_clone(tree),
        selected[[selection_name]]
      )
      .stage2b2b_compare_capture(
        old, current,
        paste("fixture", fixture_name, "selection", selection_name)
      )
    }
  }
})

