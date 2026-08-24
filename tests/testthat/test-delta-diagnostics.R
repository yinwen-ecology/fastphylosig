test_that("Delta public MCMC defaults use the documented production controls", {
  defaults <- formals(fast_delta)

  expect_equal(defaults$mcmc_sim, 10000)
  expect_equal(defaults$thin, 10)
  expect_equal(defaults$burn, 100)
})

test_that("Delta exposes saved-chain diagnostics without changing the estimator fields", {
  set.seed(701)
  tree <- ape::rtree(18)
  x <- rep(c("a", "b", "c"), length.out = ape::Ntip(tree))
  names(x) <- tree$tip.label

  fit <- suppressWarnings(fast_delta(
    tree, x, mcmc_sim = 80, thin = 5, burn = 20, verbose = FALSE
  ))

  expect_s3_class(fit, "phylo_delta")
  expect_true(is.finite(fit$delta))
  expect_equal(fit$n_saved, 26) # two chains: 1 + floor((80 - 20) / 5)
  expect_equal(fit$n_saved_requested, 26)
  expect_equal(fit$n_saved_successful, 26)
  expect_equal(fit$requested_iterations, 160)
  expect_equal(fit$successful_iterations, 160)
  expect_true(is.finite(fit$alpha_sd))
  expect_true(is.finite(fit$beta_sd))
  expect_true(is.finite(fit$ESS_alpha))
  expect_true(is.finite(fit$ESS_beta))
  expect_true(is.finite(fit$MCSE_Delta))
  expect_true(isTRUE(fit$diagnostics_available))
  expect_true(is.list(fit$stochastic))
  expect_true(length(fit$stochastic$rng_kind) >= 1L)
  expect_true(is.numeric(fit$stochastic$seed$checksum) ||
              is.na(fit$stochastic$seed$checksum))
  expect_false("raw_seed" %in% names(fit$stochastic))
  expect_false(".Random.seed" %in% names(fit$stochastic))
})

test_that("Permutation P MCSE uses finite successful permutations", {
  set.seed(702)
  tree <- ape::rtree(14)
  x <- rep(c("a", "b", "c"), length.out = ape::Ntip(tree))
  names(x) <- tree$tip.label
  fit <- suppressWarnings(fast_delta(
    tree, x, test = TRUE, nsim = 5, mcmc_sim = 40, thin = 5,
    burn = 10, return_sim = TRUE, verbose = FALSE
  ))

  finite <- is.finite(fit$sim.delta)
  expect_equal(fit$successful_simulations, sum(finite))
  if (is.finite(fit$P) && sum(finite) > 0) {
    expect_equal(
      fit$P_MCSE,
      sqrt(fit$P * (1 - fit$P) / sum(finite)),
      tolerance = 1e-12
    )
  } else {
    expect_true(is.na(fit$P_MCSE))
  }
  expect_equal(fit$stochastic$requested_permutations, 5L)
})

test_that("R diagnostics fallback accounts for alpha-beta covariance", {
  alpha <- cbind(c(1, 2, 3, 4, 5, 6), c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5))
  beta <- 2 * alpha
  out <- fastphylosig:::.delta_diagnostics_fallback(alpha, beta)

  expect_true(isTRUE(out$available))
  expect_equal(out$alpha_mean, mean(alpha))
  expect_equal(out$beta_mean, mean(beta))
  expect_equal(out$alpha_beta_cov, out$mean_covariance[1, 2])
  expect_true(is.finite(out$MCSE_Delta))
  expect_true(out$ESS_alpha <= length(alpha))
  expect_true(out$ESS_beta <= length(beta))
})

test_that("Delta reuses compact shared preparation metadata", {
  set.seed(703)
  tree <- ape::rtree(12)
  x <- rep(c("a", "b", "c"), length.out = ape::Ntip(tree))
  names(x) <- tree$tip.label
  X <- data.frame(
    trait_a = x,
    trait_b = x,
    row.names = names(x),
    stringsAsFactors = FALSE
  )
  X$trait_b[c(2, 5)] <- NA_character_

  fit <- suppressWarnings(fast_delta(
    tree, X, mcmc_sim = 30, thin = 5, burn = 10, verbose = FALSE
  ))
  metadata <- attr(fit, "analysis_metadata")

  expect_true(is.list(metadata))
  expect_identical(metadata$signal, "Delta")
  expect_identical(metadata$raw_or_prepared, "raw")
  expect_true(is.numeric(metadata$na_pattern_count))
  expect_length(metadata$na_pattern_keys, metadata$na_pattern_count)
  expect_true(is.data.frame(metadata$matching))
  expect_false(any(c("ctx", "groups", "matched_data") %in%
                   names(metadata)))
  expect_identical(
    attr(fit, "match_report")$analysis_metadata,
    metadata
  )
})

test_that("Delta reports actionable errors for a blocking raw tree", {
  tree <- ape::stree(6, type = "star")
  tree$edge.length <- rep(1, nrow(tree$edge))
  x <- setNames(rep(c("a", "b", "c"), length.out = 6), tree$tip.label)

  expect_error(
    fast_delta(tree, x, mcmc_sim = 20, thin = 5, burn = 5,
               verbose = FALSE),
    "USER_ACTION_REQUIRED|check_tree"
  )
})

test_that("Delta keeps explicit notes for single-state and short traits", {
  set.seed(704)
  tree <- ape::rtree(8)
  tips <- tree$tip.label
  X <- data.frame(
    valid = rep(c("a", "b", "c"), length.out = length(tips)),
    single_state = rep("a", length(tips)),
    short = c("a", rep(NA_character_, length(tips) - 1L)),
    row.names = tips,
    stringsAsFactors = FALSE
  )

  fit <- suppressWarnings(fast_delta(
    tree, X, mcmc_sim = 20, thin = 5, burn = 5, verbose = FALSE
  ))

  expect_true(is.finite(fit$Delta_fast[[1L]]))
  expect_true(is.na(fit$Delta_fast[[2L]]))
  expect_match(fit$note[[2L]], "single state")
  expect_true(is.na(fit$Delta_fast[[3L]]))
  expect_match(fit$note[[3L]], "fewer than 2")
  metadata <- attr(fit, "analysis_metadata")$delta_retained_species
  expect_equal(metadata$by_trait$sample_size_status,
               c("warning", "warning", "error"))
})

test_that("Delta warns once for a batch with fewer than 20 retained species", {
  set.seed(705)
  tree <- ape::rtree(25)
  tips <- tree$tip.label
  categorical <- rep(c("a", "b", "c"), length.out = length(tips))
  X <- data.frame(
    full = categorical,
    retained_19 = categorical,
    retained_15 = categorical,
    row.names = tips,
    stringsAsFactors = FALSE
  )
  X$retained_19[1:6] <- NA_character_
  X$retained_15[1:10] <- NA_character_

  warnings <- character()
  fit <- withCallingHandlers(
    fast_delta(
      tree, X, test = FALSE, mcmc_sim = 20, thin = 5, burn = 5,
      verbose = FALSE, progress = FALSE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  delta_warnings <- warnings[grepl("^Delta retained", warnings)]
  expect_length(delta_warnings, 1L)
  expect_match(delta_warnings, "Delta retained 2 to 19 species")
  expect_match(delta_warnings, "retained_19 \\(n = 19\\)")
  expect_match(delta_warnings, "retained_15 \\(n = 15\\)")
  expect_true(all(is.finite(fit$Delta_fast)))
  expect_equal(fit$n_species, c(25L, 19L, 15L))

  metadata <- attr(fit, "analysis_metadata")$delta_retained_species
  expect_true(is.list(metadata))
  expect_identical(metadata$warning_threshold, 20L)
  expect_true(metadata$warning_emitted)
  expect_equal(metadata$by_trait$sample_size_status,
               c("ok", "warning", "warning"))
  expect_equal(metadata$affected_traits$trait, c("retained_19", "retained_15"))
  expect_equal(metadata$affected_traits$n_species, c(19L, 15L))
  expect_equal(metadata$affected_traits$n_removed_na, c(6L, 10L))
})

test_that("Delta does not warn when every retained trait has at least 20 species", {
  set.seed(706)
  tree <- ape::rtree(20)
  x <- setNames(rep(c("a", "b", "c"), length.out = 20), tree$tip.label)

  warnings <- character()
  fit <- withCallingHandlers(
    fast_delta(
      tree, x, test = FALSE, mcmc_sim = 20, thin = 5, burn = 5,
      verbose = FALSE, progress = FALSE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_length(warnings[grepl("^Delta retained", warnings)], 0L)
  metadata <- attr(fit, "analysis_metadata")$delta_retained_species
  expect_false(metadata$warning_emitted)
  expect_equal(nrow(metadata$affected_traits), 0L)
})

test_that("Delta MCMC diagnostic warnings are aggregated and retain diagnostics", {
  set.seed(707)
  tree <- ape::rtree(20)
  values <- rep(c("a", "b", "c"), length.out = ape::Ntip(tree))
  X <- data.frame(
    trait_a = values,
    trait_b = values,
    row.names = tree$tip.label,
    stringsAsFactors = FALSE
  )

  warnings <- character()
  fit <- withCallingHandlers(
    fast_delta(
      tree, X, test = FALSE, mcmc_sim = 20, thin = 5, burn = 5,
      verbose = FALSE, progress = FALSE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  mcmc_warnings <- warnings[grepl("^Delta MCMC diagnostics", warnings)]
  expect_length(mcmc_warnings, 1L)
  expect_match(mcmc_warnings, "Estimates are returned unchanged")
  expect_true(all(is.finite(fit$Delta_fast)))
  expect_equal(fit$requested_iterations, rep(40, 2L))
  expect_equal(fit$successful_iterations, rep(40, 2L))
  expect_true(all(is.finite(fit$ESS_alpha)))
  expect_true(all(is.finite(fit$ESS_beta)))

  metadata <- attr(fit, "analysis_metadata")$delta_mcmc_diagnostics
  expect_true(is.list(metadata))
  expect_identical(metadata$ess_minimum, 20)
  expect_identical(metadata$split_rhat_maximum, 1.1)
  expect_true(metadata$warning_emitted)
  expect_equal(metadata$affected_traits$trait, c("trait_a", "trait_b"))
  expect_true(all(metadata$affected_traits$warning))
})

test_that("Delta MCMC warning summary only considers finite estimates", {
  summary <- fastphylosig:::.delta_mcmc_warning_summary(
    trait = c("ess", "rhat", "mcse", "failed"),
    delta = c(0.5, 0.5, 0.5, NA_real_),
    ess_alpha = c(19, 20, 20, NA_real_),
    ess_beta = c(20, 20, 20, NA_real_),
    rhat_alpha = c(1, 1.11, 1, NA_real_),
    rhat_beta = c(1, 1, 1, NA_real_),
    mcse_delta = c(0.01, 0.01, NA_real_, NA_real_)
  )

  expect_true(summary$warning_emitted)
  expect_equal(summary$affected_traits$trait, c("ess", "rhat", "mcse"))
  expect_match(summary$by_trait$reason[[1L]], "ESS_alpha < 20")
  expect_match(summary$by_trait$reason[[2L]], "split_Rhat_alpha > 1.1")
  expect_match(summary$by_trait$reason[[3L]], "MCSE_Delta is non-finite")
  expect_false(summary$by_trait$warning[[4L]])
  expect_true(is.na(summary$by_trait$reason[[4L]]))
})
