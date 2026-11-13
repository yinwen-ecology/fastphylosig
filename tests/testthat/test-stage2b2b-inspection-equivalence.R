test_that("inspection reports current readiness for valid and malformed trees", {
  skip_if_not_installed("ape")
  set.seed(20260904)
  fixtures <- .inspection_fixture_set()
  signals <- c("K", "lambda", "D", "Delta")
  all_ready <- c("balanced", "random", "pectinate", "shuffled_edge_order",
                 "safe_internal_renumbering", "different_root_numbering",
                 "two_tip", "large_pectinate")
  all_blocked <- c("cycle", "disconnected", "malformed_edge", "invalid_nnode",
                   "negative_branch", "zero_terminal_branch")
  expected_codes <- list(
    polytomy = c("polytomy_supported", "d_polytomy_compatibility",
                 "delta_requires_binary_tree"),
    cycle = "invalid_topology", disconnected = "invalid_topology",
    unary = c("single_child_internal", "d_requires_no_single_child",
              "delta_requires_binary_tree"),
    malformed_edge = "missing_edge", invalid_nnode = "invalid_nnode",
    negative_branch = "negative_branch_lengths",
    zero_terminal_branch = c("k_requires_positive_branches",
                            "lambda_zero_terminal_branch",
                            "d_requires_positive_branches",
                            "delta_requires_positive_branches"),
    invalid_root_representation = c("conventional_unrooted",
                                    "delta_requires_conventional_root")
  )
  for (name in names(fixtures)) {
    tree <- fixtures[[name]]
    before <- serialize(tree, NULL, version = 2L)
    expect_silent(check <- fastphylosig:::.inspect_tree_core(tree, signals))
    expect_s3_class(check, "fastphylosig_tree_check")
    expect_identical(names(check$ready_by_signal), signals)
    expect_identical(check$ready, all(check$ready_by_signal))
    expect_identical(serialize(tree, NULL, version = 2L), before, info = name)
    if (name %in% all_ready) expect_true(all(check$ready_by_signal), info = name)
    if (name %in% all_blocked) expect_false(any(check$ready_by_signal), info = name)
    if (name %in% c("polytomy", "invalid_root_representation")) {
      expect_identical(unname(check$ready_by_signal), c(TRUE, TRUE, TRUE, FALSE))
    }
    if (name == "unary") {
      expect_identical(unname(check$ready_by_signal), c(TRUE, TRUE, FALSE, FALSE))
    }
    if (name %in% names(expected_codes)) {
      expect_true(all(expected_codes[[name]] %in% check$issues$code), info = name)
    }
    expect_true(all(check$issues$severity %in%
                      c("INFO", "WARNING", "ERROR", "USER_ACTION_REQUIRED")))
    expect_true(all(check$issues$signal %in% c("all", signals)))
    expect_true(all(nzchar(check$issues$message)))
  }
})

test_that("inspection signal selection preserves diagnostics and condition contracts", {
  skip_if_not_installed("ape")
  set.seed(20260904)
  fixtures <- .inspection_fixture_set()
  for (name in names(fixtures)) {
    tree <- fixtures[[name]]
    full <- fastphylosig:::.inspect_tree_core(tree)
    for (selection in list("K", c("lambda", "D"), c("Delta", "K"))) {
      expect_silent(current <- fastphylosig:::.inspect_tree_core(tree, selection))
      expect_identical(current$ready_by_signal, full$ready_by_signal[selection])
      expect_identical(current$tree_summary, full$tree_summary)
      expected <- full$issues[full$issues$signal %in% c("all", selection), , drop = FALSE]
      # Selection controls the per-method emission order, not the issue content.
      key <- function(x) order(x$code, x$signal, x$severity, x$message)
      expected <- expected[key(expected), , drop = FALSE]
      actual <- current$issues[key(current$issues), , drop = FALSE]
      rownames(expected) <- rownames(actual) <- NULL
      expect_identical(actual, expected, info = name)
    }
  }
  tree <- fixtures$balanced
  for (selection in list(NULL, character(), NA_character_, 1)) {
    expect_error(fastphylosig:::.inspect_tree_core(tree, selection),
                 "signal must contain", class = "error")
  }
  expect_error(fastphylosig:::.inspect_tree_core(tree, "unknown"),
               "signal must contain only", class = "error")
  aliases <- fastphylosig:::.inspect_tree_core(tree, c("pagel", "phylo.d", "k", "k"))
  expect_identical(names(aliases$ready_by_signal), c("lambda", "D", "K"))
})
