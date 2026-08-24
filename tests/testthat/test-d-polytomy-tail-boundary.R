test_that("D polytomy one-ULP strict-tail boundary remains explicit", {
  skip_if_not_installed("ape")

  fixture <- testthat::test_path(
    "fixtures", "d_polytomy_tail_boundary_audit.csv"
  )
  skip_if_not(
    file.exists(fixture),
    "fixed D polytomy tail-boundary audit fixture is unavailable"
  )

  audit <- utils::read.csv(
    fixture, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = c(shared_state_bits = "character")
  )
  nsim <- 199L
  classification <- "D_IEEE_754_STRICT_BOUNDARY_ROUNDING"

  expect_equal(nrow(audit), nsim)
  expect_identical(audit$draw, seq_len(nsim))
  expect_true(all(audit$audit_status == "RUN"))
  expect_true(all(audit$source_classification == classification))
  expect_true(all(!audit$source_A_shared_state_mismatch))
  expect_true(all(!audit$source_B_tree_or_order_mismatch))
  expect_true(all(!audit$source_C_nonboundary_statistic_mismatch))
  expect_true(all(audit$source_D_ieee_boundary_rounding))
  expect_true(all(!audit$source_E_tail_semantics_mismatch))
  expect_true(all(audit$shared_state_identical))
  expect_true(all(audit$topology_identical))
  expect_true(all(grepl(
    "literal > tail; no epsilon or tail change.", audit$note, fixed = TRUE
  )))

  # These columns record the audit's literal `null > observed` operations.
  # Raw hex and ULP assertions below prevent decimal parsing or rounding from
  # being substituted for that comparison.
  caper_strict_gt <- audit$caper_strict_gt
  fast_strict_gt <- audit$fast_strict_gt
  expect_identical(
    audit$strict_tail_disagreement,
    caper_strict_gt != fast_strict_gt
  )
  expect_identical(
    which(caper_strict_gt != fast_strict_gt),
    132L
  )

  boundary <- audit[audit$draw == 132L, , drop = FALSE]
  expect_false(boundary$caper_strict_gt)
  expect_true(boundary$fast_strict_gt)
  expect_identical(boundary$observed_ulp_distance, 1L)
  expect_identical(boundary$null_ulp_distance, 0L)
  expect_identical(boundary$caper_observed_hex, boundary$caper_null_hex)
  expect_identical(boundary$caper_null_hex, boundary$fast_null_hex)
  expect_identical(boundary$caper_replay_null_hex, boundary$caper_null_hex)
  expect_identical(boundary$fast_observed_hex, "40133f1a53b535dc")
  expect_identical(boundary$caper_observed_hex, "40133f1a53b535dd")

  branch <- function(text) paste0(text, ":1")
  binary_subtree <- function(labels) {
    if (length(labels) == 1L) return(branch(labels))
    branch(paste0(
      "(", branch(labels[[1L]]), ",", binary_subtree(labels[-1L]), ")"
    ))
  }
  labels <- paste0("t", seq_len(20L))
  polytomy <- branch(paste0(
    "(",
    paste(
      c(
        branch(labels[[1L]]),
        branch(labels[[2L]]),
        binary_subtree(labels[3L:16L])
      ),
      collapse = ","
    ),
    ")"
  ))
  tree <- ape::read.tree(text = paste0(
    "(", polytomy, ",", binary_subtree(labels[17L:20L]), ");"
  ))
  trait <- stats::setNames(integer(20L), labels)
  trait[c(15L, 8L, 13L, 16L)] <- 1L

  state_columns <- strsplit(audit$shared_state_bits, "", fixed = TRUE)
  expect_true(all(lengths(state_columns) == 20L))
  brownian_states <- matrix(
    as.numeric(unlist(state_columns, use.names = FALSE)),
    nrow = 20L, ncol = nsim
  )
  rownames(brownian_states) <- tree$tip.label

  production <- suppressWarnings(suppressMessages(fast_d(
    tree, trait, test = TRUE, nsim = nsim,
    random_states = brownian_states,
    brownian_states = brownian_states,
    return_sim = TRUE, keep_null = TRUE,
    verbose = FALSE, progress = FALSE
  )))

  caper_count <- sum(caper_strict_gt)
  production_count <- sum(
    production$Permutations$brownian > production$Parameters$Observed
  )
  expect_identical(caper_count, 102L)
  expect_identical(sum(fast_strict_gt), 103L)
  expect_identical(production_count, 103L)
  expect_identical(production_count, caper_count + 1L)
  expect_identical(production$nsim_successful_brownian, nsim)
  expect_equal(production$P_Brownian, production_count / nsim, tolerance = 0)
})
