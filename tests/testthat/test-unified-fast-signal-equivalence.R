# End-to-end equivalence checks for the fixed-method specialists and the
# unified fast_signal() dispatcher.  The fixture is deliberately fixed (and
# all stochastic inputs are supplied explicitly) so this file is stable on
# small CI machines as well as with the compiled .rlib_api package.

.fast_signal_e2e_fixture <- function() {
  tree <- ape::read.tree(text = paste0(
    "(((a:0.5,b:0.7):0.8,(c:0.4,d:0.6):0.9):0.3,",
    "((e:0.2,f:0.3):0.4,(g:0.5,h:0.6):0.7):0.8);"
  ))
  tree <- ape::reorder.phylo(tree, "postorder")
  tips <- tree$tip.label
  continuous <- stats::setNames(
    c(-0.5, -0.2, 0.1, 0.4, 0.7, 1.0, 1.3, 1.6), tips
  )
  binary <- stats::setNames(c(0, 0, 0, 1, 1, 1, 0, 1), tips)
  categorical <- stats::setNames(
    c("a", "b", "c", "a", "b", "c", "a", "b"), tips
  )
  permutations <- rbind(
    seq_along(tips),
    rev(seq_along(tips)),
    c(2L, 1L, 4L, 3L, 6L, 5L, 8L, 7L)
  )
  random_states <- rbind(
    c(0, 1), c(1, 0), c(0, 1), c(1, 0),
    c(0, 1), c(1, 0), c(0, 1), c(1, 0)
  )
  brownian_states <- rbind(
    c(0, 0), c(0, 0), c(1, 1), c(1, 1),
    c(0, 0), c(0, 0), c(1, 1), c(1, 1)
  )
  list(
    tree = tree, tips = tips, continuous = continuous, binary = binary,
    categorical = categorical, permutations = permutations,
    random_states = random_states, brownian_states = brownian_states
  )
}

.expect_fast_signal_workflow <- function(x, method, production,
                                         input = "raw") {
  workflow <- attr(x, "workflow", exact = TRUE)
  required <- c(
    "method", "production_function", "input_tree_type", "raw_or_prepared",
    "tree_auto_normalized", "matching_performed", "tree_tips_removed",
    "data_rows_removed", "NA_patterns", "total_elapsed"
  )
  testthat::expect_true(is.list(workflow))
  testthat::expect_true(all(required %in% names(workflow)))
  testthat::expect_identical(workflow$method, method)
  testthat::expect_identical(workflow$production_function, production)
  testthat::expect_identical(workflow$input_tree_type, input)
  testthat::expect_identical(workflow$raw_or_prepared, input)
  testthat::expect_true(is.logical(workflow$tree_auto_normalized))
  testthat::expect_true(isTRUE(workflow$matching_performed))
  testthat::expect_true(is.numeric(workflow$total_elapsed))
  testthat::expect_true(is.finite(workflow$total_elapsed))
}

.expect_named_numeric_fields <- function(a, b, fields, tolerance = 1e-10) {
  for (field in fields) {
    testthat::expect_equal(a[[field]], b[[field]], tolerance = tolerance,
                           info = field)
  }
}

test_that("K dispatcher equals fast_k for controlled permutations and NA masks", {
  testthat::skip_if_not_installed("ape")
  f <- .fast_signal_e2e_fixture()
  specialist <- suppressMessages(fast_k(
    f$tree, f$continuous, test = TRUE, nsim = nrow(f$permutations),
    permutations = f$permutations, return_sim = TRUE,
    verbose = FALSE, progress = FALSE
  ))
  unified <- suppressMessages(fast_signal(
    f$tree, data = f$continuous, method = "K", test = TRUE,
    nsim = nrow(f$permutations), permutations = f$permutations,
    return_sim = TRUE, verbose = FALSE, progress = FALSE
  ))

  .expect_named_numeric_fields(
    unified, specialist,
    c("K_fast", "n_species", "n_removed_na", "P_fast",
      "nsim_requested", "nsim_successful", "exceedance_count"),
    tolerance = 1e-12
  )
  testthat::expect_equal(unified$sim.K_fast, specialist$sim.K_fast,
                         tolerance = 1e-12)
  metadata <- attr(unified, "analysis_metadata")
  if (length(metadata$matched_species) == 1L &&
      is.numeric(metadata$matched_species)) {
    testthat::expect_equal(metadata$matched_species, 8L)
  } else {
    testthat::expect_equal(metadata$matched_species, f$tips)
  }
  testthat::expect_length(metadata$tree_tips_removed, 0L)
  testthat::expect_length(metadata$data_rows_removed, 0L)
  X <- cbind(trait_1 = f$continuous, trait_2 = 2 * f$continuous)
  X["a", "trait_2"] <- NA_real_
  na_specialist <- suppressMessages(fast_k(
    f$tree, X, test = FALSE, verbose = FALSE, progress = FALSE
  ))
  na_fit <- suppressMessages(fast_signal(
    f$tree, data = X, method = "K", test = FALSE,
    verbose = FALSE, progress = FALSE
  ))
  testthat::expect_equal(na_fit$K_fast, na_specialist$K_fast,
                         tolerance = 1e-12)
  testthat::expect_equal(na_fit$n_species, na_specialist$n_species)
  testthat::expect_equal(na_fit$n_removed_na, na_specialist$n_removed_na)
  testthat::expect_equal(
    attr(na_fit, "analysis_metadata")$na_pattern_count,
    attr(na_specialist, "analysis_metadata")$na_pattern_count
  )
  testthat::expect_equal(na_fit$n_removed_na, c(0L, 1L))
  testthat::expect_equal(attr(na_fit, "workflow")$NA_patterns, 2L)

  mismatch <- c(f$continuous[f$tips[1:7]], outside = 9)
  mismatch_specialist <- suppressMessages(fast_k(
    f$tree, mismatch, verbose = FALSE, progress = FALSE
  ))
  mismatch_unified <- suppressMessages(fast_signal(
    f$tree, data = mismatch, method = "K",
    verbose = FALSE, progress = FALSE
  ))
  testthat::expect_equal(as.numeric(mismatch_unified),
                         as.numeric(mismatch_specialist), tolerance = 1e-12)
  mismatch_meta <- attr(mismatch_unified, "analysis_metadata")
  testthat::expect_equal(mismatch_meta$matching$matched_species, 7L)
  testthat::expect_equal(mismatch_meta$matching$removed_tree_tips, 1L)
  testthat::expect_equal(mismatch_meta$matching$removed_data_rows, 1L)
  testthat::expect_equal(attr(mismatch_specialist, "matching"),
                         attr(mismatch_unified, "matching"))
  testthat::expect_equal(attr(mismatch_unified, "workflow")$tree_tips_removed,
                         1L)
  testthat::expect_equal(attr(mismatch_unified, "workflow")$data_rows_removed,
                         1L)
  .expect_fast_signal_workflow(unified, "K", "fast_k")
  .expect_fast_signal_workflow(specialist, "K", "fast_k")
})

test_that("lambda profile is identical through fast_lambda and fast_signal", {
  testthat::skip_if_not_installed("ape")
  f <- .fast_signal_e2e_fixture()
  specialist <- suppressMessages(fast_lambda(
    f$tree, x = f$continuous, lambda_profile = TRUE,
    lambda_profile_points = 11L, verbose = FALSE, progress = FALSE
  ))
  unified <- suppressMessages(fast_signal(
    f$tree, data = f$continuous, method = "Lambda",
    lambda_profile = TRUE, lambda_profile_points = 11L,
    verbose = FALSE, progress = FALSE
  ))
  .expect_named_numeric_fields(
    unified, specialist,
    c("lambda", "logL", "gls_mean", "sig2", "lambda_CI",
      "lambda_CI_level", "lambda_CI_cutoff"),
    tolerance = 1e-10
  )
  testthat::expect_equal(unified$lambda_profile, specialist$lambda_profile,
                         tolerance = 1e-10)
  .expect_fast_signal_workflow(unified, "lambda", "fast_lambda")
  .expect_fast_signal_workflow(specialist, "lambda", "fast_lambda")
  testthat::expect_equal(attr(unified, "production_function"), "fast_lambda")
})

test_that("D dispatcher equals fast_d for fixed random and Brownian states", {
  testthat::skip_if_not_installed("ape")
  f <- .fast_signal_e2e_fixture()
  specialist <- suppressMessages(fast_d(
    f$tree, f$binary, test = TRUE, nsim = 2L,
    random_states = f$random_states, brownian_states = f$brownian_states,
    return_sim = TRUE, verbose = FALSE, progress = FALSE
  ))
  unified <- suppressMessages(fast_signal(
    f$tree, data = f$binary, method = "d", test = TRUE, nsim = 2L,
    random_states = f$random_states, brownian_states = f$brownian_states,
    return_sim = TRUE, verbose = FALSE, progress = FALSE
  ))
  .expect_named_numeric_fields(
    unified, specialist,
    c("DEstimate", "Pval1", "Pval0", "P_random", "P_Brownian",
      "nsim_requested", "nsim_successful_random",
      "nsim_successful_brownian"),
    tolerance = 1e-12
  )
  testthat::expect_equal(unified$Permutations, specialist$Permutations,
                         tolerance = 1e-12)
  .expect_fast_signal_workflow(unified, "D", "fast_d")
  .expect_fast_signal_workflow(specialist, "D", "fast_d")
})

test_that("Delta dispatcher is reproducible with fixed seed and permutations", {
  testthat::skip_if_not_installed("ape")
  f <- .fast_signal_e2e_fixture()
  controls <- list(
    test = TRUE, nsim = 2L, permutations = f$permutations[1:2, , drop = FALSE],
    return_sim = TRUE, mcmc_sim = 8L, thin = 2L, burn = 2L,
    verbose = FALSE, progress = FALSE
  )
  set.seed(8128)
  specialist <- suppressWarnings(do.call(
    fast_delta, c(list(tree = f$tree, x = f$categorical), controls)
  ))
  set.seed(8128)
  unified <- suppressWarnings(do.call(
    fast_signal, c(list(tree = f$tree, data = f$categorical,
                        method = "Delta"), controls)
  ))
  .expect_named_numeric_fields(
    unified, specialist,
    c("delta", "alpha_mean", "beta_mean", "n_saved",
      "n_saved_successful", "MCSE_Delta", "P", "P_MCSE",
      "n_failed_sim"),
    tolerance = 1e-10
  )
  testthat::expect_equal(unified$sim.delta, specialist$sim.delta,
                         tolerance = 1e-10)
  testthat::expect_equal(unified$stochastic$requested_permutations,
                         2L)
  .expect_fast_signal_workflow(unified, "Delta", "fast_delta")
  .expect_fast_signal_workflow(specialist, "Delta", "fast_delta")

  ctx <- prepare_tree(f$tree)
  prepared_controls <- controls
  prepared_controls$mcmc_sim <- 6L
  set.seed(8129)
  prepared_specialist <- suppressWarnings(do.call(
    fast_delta, c(list(tree = ctx, x = f$categorical), prepared_controls)
  ))
  set.seed(8129)
  prepared_unified <- suppressWarnings(do.call(
    fast_signal, c(list(tree = f$tree, data = f$categorical,
                        method = "Delta", prepared = ctx),
                   prepared_controls)
  ))
  testthat::expect_equal(prepared_unified$delta,
                         prepared_specialist$delta, tolerance = 1e-10)
  testthat::expect_equal(prepared_unified$P, prepared_specialist$P,
                         tolerance = 1e-10)
  testthat::expect_equal(prepared_unified$sim.delta,
                         prepared_specialist$sim.delta, tolerance = 1e-10)
  .expect_fast_signal_workflow(prepared_specialist, "Delta", "fast_delta",
                               "prepared")
  .expect_fast_signal_workflow(prepared_unified, "Delta", "fast_delta",
                               "prepared")
})

test_that("prepared contexts preserve specialist results and report provenance", {
  testthat::skip_if_not_installed("ape")
  f <- .fast_signal_e2e_fixture()
  ctx <- prepare_tree(f$tree)
  raw_k <- suppressMessages(fast_k(
    f$tree, f$continuous, verbose = FALSE, progress = FALSE
  ))
  prepared_k <- suppressMessages(fast_k(
    ctx, f$continuous, verbose = FALSE, progress = FALSE
  ))
  prepared_high_k <- suppressMessages(fast_signal(
    f$tree, data = f$continuous, method = "K", prepared = ctx,
    verbose = FALSE, progress = FALSE
  ))
  raw_lambda <- suppressMessages(fast_lambda(
    f$tree, f$continuous, lambda_profile = TRUE,
    lambda_profile_points = 11L, verbose = FALSE, progress = FALSE
  ))
  prepared_lambda <- suppressMessages(fast_lambda(
    ctx, f$continuous, lambda_profile = TRUE,
    lambda_profile_points = 11L, verbose = FALSE, progress = FALSE
  ))
  prepared_high_lambda <- suppressMessages(fast_signal(
    f$tree, data = f$continuous, method = "lambda", prepared = ctx,
    lambda_profile = TRUE, lambda_profile_points = 11L,
    verbose = FALSE, progress = FALSE
  ))
  raw_d <- suppressMessages(fast_d(
    f$tree, f$binary, test = TRUE, nsim = 2L,
    random_states = f$random_states, brownian_states = f$brownian_states,
    return_sim = TRUE, verbose = FALSE, progress = FALSE
  ))
  prepared_d <- suppressMessages(fast_d(
    ctx, f$binary, test = TRUE, nsim = 2L,
    random_states = f$random_states, brownian_states = f$brownian_states,
    return_sim = TRUE, verbose = FALSE, progress = FALSE
  ))
  raw_high_d <- suppressMessages(fast_signal(
    f$tree, data = f$binary, method = "D", test = TRUE, nsim = 2L,
    random_states = f$random_states, brownian_states = f$brownian_states,
    return_sim = TRUE, verbose = FALSE, progress = FALSE
  ))
  prepared_high_d <- suppressMessages(fast_signal(
    f$tree, data = f$binary, method = "D", prepared = ctx,
    test = TRUE, nsim = 2L, random_states = f$random_states,
    brownian_states = f$brownian_states, return_sim = TRUE,
    verbose = FALSE, progress = FALSE
  ))
  testthat::expect_equal(as.numeric(prepared_k), as.numeric(raw_k),
                         tolerance = 1e-12)
  testthat::expect_equal(as.numeric(prepared_high_k), as.numeric(raw_k),
                         tolerance = 1e-12)
  testthat::expect_equal(prepared_lambda$lambda, raw_lambda$lambda,
                         tolerance = 1e-10)
  testthat::expect_equal(prepared_lambda$logL, raw_lambda$logL,
                         tolerance = 1e-10)
  testthat::expect_equal(prepared_high_lambda$lambda, raw_lambda$lambda,
                         tolerance = 1e-10)
  testthat::expect_equal(prepared_high_lambda$lambda_profile,
                         raw_lambda$lambda_profile, tolerance = 1e-10)
  testthat::expect_equal(prepared_d$DEstimate, raw_d$DEstimate,
                         tolerance = 1e-12)
  testthat::expect_equal(raw_high_d$DEstimate, raw_d$DEstimate,
                         tolerance = 1e-12)
  .expect_fast_signal_workflow(prepared_k, "K", "fast_k", "prepared")
  .expect_fast_signal_workflow(prepared_high_k, "K", "fast_k", "prepared")
  .expect_fast_signal_workflow(prepared_lambda, "lambda", "fast_lambda",
                               "prepared")
  .expect_fast_signal_workflow(prepared_high_lambda, "lambda", "fast_lambda",
                               "prepared")
  .expect_fast_signal_workflow(raw_high_d, "D", "fast_d")
  .expect_fast_signal_workflow(prepared_high_d, "D", "fast_d", "prepared")
  testthat::expect_false(isTRUE(attr(prepared_k, "workflow")$tree_auto_normalized))
})

test_that("safe canonicalization is representation-only and input objects stay immutable", {
  testthat::skip_if_not_installed("ape")
  f <- .fast_signal_e2e_fixture()
  ord <- rev(seq_len(nrow(f$tree$edge)))
  scrambled <- f$tree
  scrambled$edge <- f$tree$edge[ord, , drop = FALSE]
  scrambled$edge.length <- f$tree$edge.length[ord]
  tree_before <- serialize(scrambled, NULL)
  data_before <- serialize(f$continuous, NULL)
  fit_scrambled <- suppressMessages(fast_signal(
    scrambled, data = f$continuous, method = "K",
    verbose = FALSE, progress = FALSE
  ))
  fit_reference <- suppressMessages(fast_signal(
    f$tree, data = f$continuous, method = "K",
    verbose = FALSE, progress = FALSE
  ))
  testthat::expect_equal(as.numeric(fit_scrambled), as.numeric(fit_reference),
                         tolerance = 1e-12)
  testthat::expect_identical(serialize(scrambled, NULL), tree_before)
  testthat::expect_identical(serialize(f$continuous, NULL), data_before)
  testthat::expect_identical(scrambled$tip.label, f$tree$tip.label)
  testthat::expect_equal(sort(scrambled$edge.length),
                         sort(f$tree$edge.length))

  resolved <- resolve_tree(scrambled, signal = "K")
  again <- resolve_tree(resolved, signal = "K")
  testthat::expect_equal(ape::cophenetic.phylo(resolved),
                         ape::cophenetic.phylo(f$tree))
  testthat::expect_length(
    attr(again, "fastphylosig_resolution")$changes, 0L
  )
})

test_that("dispatcher aliases, required method, and inapplicable arguments fail clearly", {
  testthat::skip_if_not_installed("ape")
  f <- .fast_signal_e2e_fixture()
  k <- suppressMessages(fast_signal(
    f$tree, f$continuous, method = "k", verbose = FALSE, progress = FALSE
  ))
  lambda <- suppressMessages(fast_signal(
    f$tree, f$continuous, method = intToUtf8(0x03bb), lambda_profile = FALSE,
    verbose = FALSE, progress = FALSE
  ))
  d <- suppressMessages(fast_signal(
    f$tree, f$binary, method = "d", test = FALSE,
    verbose = FALSE, progress = FALSE
  ))
  delta <- suppressWarnings(suppressMessages(fast_signal(
    f$tree, f$categorical, method = "delta", test = FALSE,
    mcmc_sim = 6L, thin = 2L, burn = 2L,
    verbose = FALSE, progress = FALSE
  )))
  testthat::expect_identical(attr(k, "method"), "K")
  testthat::expect_identical(attr(lambda, "method"), "lambda")
  testthat::expect_identical(attr(d, "method"), "D")
  testthat::expect_identical(attr(delta, "method"), "Delta")

  testthat::expect_error(fast_signal(f$tree, f$continuous), "method")
  testthat::expect_error(
    fast_signal(f$tree, f$continuous, method = c("K", "lambda")),
    "exactly one"
  )
  testthat::expect_error(
    fast_signal(f$tree, f$continuous, method = "unknown"),
    "Unknown method"
  )
  testthat::expect_error(
    fast_signal(f$tree, f$binary, method = "D", lambda_profile = TRUE),
    "not applicable"
  )
  testthat::expect_error(
    fast_signal(f$tree, f$continuous, method = "K", unknown_argument = 1),
    "not applicable"
  )
  testthat::expect_error(
    fast_signal(f$tree, data = f$continuous, x = f$continuous, method = "K"),
    "only one"
  )

  star <- ape::stree(8L, type = "star")
  star$edge.length <- rep(1, nrow(star$edge))
  testthat::expect_error(
    fast_signal(star, f$categorical, method = "Delta", test = FALSE,
                mcmc_sim = 6L, thin = 2L, burn = 2L,
                verbose = FALSE, progress = FALSE),
    "USER_ACTION_REQUIRED|polytom|check_tree"
  )
})
