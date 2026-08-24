candidate2_production_fixture <- function() {
  tree <- ape::read.tree(text = paste0(
    "(((a:0.5,b:0.7):0.8,(c:0.4,d:0.6):0.9):0.3,",
    "((e:0.2,f:0.3):0.4,(g:0.5,h:0.6):0.7):0.8);"
  ))
  tree <- ape::reorder.phylo(tree, "postorder")
  y <- stats::setNames(
    factor(c("A", "B", "C", "A", "B", "C", "A", "B")),
    tree$tip.label
  )
  permutations <- rbind(
    seq_along(y), rev(seq_along(y)),
    c(2L, 1L, 4L, 3L, 6L, 5L, 8L, 7L)
  )
  tree_workspace <- fastphylosig:::.fast_ace_tree_workspace(tree)
  workspace <- fastphylosig:::.fast_ace_workspace(
    tree_workspace, model = "ARD", nl = nlevels(y)
  )
  list(tree = tree, y = y, permutations = permutations,
       workspace = workspace)
}

candidate2_fit <- function(y, workspace, estimate_se) {
  y <- y[workspace$tip_label]
  suppressWarnings(fastphylosig:::.fast_ace_fit_workspace(
    tip_state = as.integer(y), lvls = levels(y), workspace = workspace,
    ip = 0.1, CI = TRUE, marginal = FALSE, estimate_se = estimate_se
  ))
}

candidate2_legacy_delta_one <- function(fixture, y, keep_chains) {
  fit <- candidate2_fit(y, fixture$workspace, estimate_se = TRUE)
  fastphylosig:::delta_mcmc_cpp(
    probabilities = fit$lik.anc, lambda0 = 0.1, proposal_sd = 0.5,
    sim = 12L, thin = 2L, burn = 2L, entropy_type = 1L,
    return_chains = keep_chains
  )
}

test_that("Candidate 2 preserves public ACE SE and fitted results", {
  skip_if_not_installed("ape")
  f <- candidate2_production_fixture()

  expect_false("engine" %in% names(formals(fast_ace)))
  public <- suppressWarnings(fast_ace(
    f$y, f$tree, model = "ARD", CI = TRUE, marginal = FALSE,
    progress = FALSE
  ))
  expect_named(public, c(
    "loglik", "rates", "se", "index.matrix", "convergence", "message",
    "lik.anc", "call", "engine", "timing"
  ))
  expect_length(public$se, length(public$rates))

  with_se <- candidate2_fit(f$y, f$workspace, estimate_se = TRUE)
  without_se <- candidate2_fit(f$y, f$workspace, estimate_se = FALSE)
  fields <- c(
    "loglik", "rates", "index.matrix", "convergence", "message", "lik.anc"
  )
  expect_equal(without_se[fields], with_se[fields], tolerance = 0)
  expect_length(without_se$se, length(without_se$rates))
  expect_true(all(is.nan(without_se$se)))
})

test_that("Candidate 2 preserves fixed-seed observed and permutation Delta", {
  skip_if_not_installed("ape")
  f <- candidate2_production_fixture()
  new_delta_one <- function(y, keep_chains) {
    fastphylosig:::.delta_one(
      tree = f$tree, y = y, lambda0 = 0.1, proposal_sd = 0.5,
      mcmc_sim = 12L, thin = 2L, burn = 2L, entropy_code = 1L,
      model = "ARD", ace_engine = "fast", keep_chains = keep_chains,
      ace_workspace = f$workspace
    )
  }

  set.seed(2050)
  legacy_observed <- candidate2_legacy_delta_one(f, f$y, TRUE)
  set.seed(2050)
  production_observed <- suppressWarnings(new_delta_one(f$y, TRUE))
  expect_equal(production_observed, legacy_observed, tolerance = 0)

  run_permutations <- function(use_legacy) {
    vapply(seq_len(nrow(f$permutations)), function(i) {
      yp <- f$y[f$permutations[i, ]]
      names(yp) <- names(f$y)
      result <- if (use_legacy) {
        candidate2_legacy_delta_one(f, yp, FALSE)
      } else {
        suppressWarnings(new_delta_one(yp, FALSE))
      }
      result$delta
    }, numeric(1))
  }
  set.seed(2051)
  legacy_permutations <- run_permutations(TRUE)
  set.seed(2051)
  production_permutations <- run_permutations(FALSE)
  expect_equal(production_permutations, legacy_permutations, tolerance = 0)
})
