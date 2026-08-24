workspace_delta_fixture <- function() {
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
  list(tree = tree, y = y, permutations = permutations)
}

test_that("ACE workspace is numerically identical to public fast_ace", {
  testthat::skip_if_not_installed("ape")
  f <- workspace_delta_fixture()
  transformed_tree <- f$tree
  transformed_tree$edge.length <- transformed_tree$edge.length^2
  for (model in c("ER", "ARD")) {
    direct <- suppressWarnings(fast_ace(
      f$y, f$tree, model = model, kappa = 2, CI = TRUE, marginal = FALSE,
      progress = FALSE
    ))
    tree_workspace <- fastphylosig:::.fast_ace_tree_workspace(
      transformed_tree, include_fingerprint = TRUE
    )
    workspace <- fastphylosig:::.fast_ace_workspace(
      tree_workspace, model = model, nl = nlevels(f$y)
    )
    cached <- suppressWarnings(fastphylosig:::.fast_ace_lik_anc_workspace(
      f$y, workspace
    ))
    testthat::expect_equal(cached, direct$lik.anc, tolerance = 0,
                           info = model)
  }
})

test_that("Delta ACE workspace preserves fixed-seed observed and permutation paths", {
  testthat::skip_if_not_installed("ape")
  f <- workspace_delta_fixture()
  tree_workspace <- fastphylosig:::.fast_ace_tree_workspace(
    f$tree, include_fingerprint = TRUE
  )
  workspace <- fastphylosig:::.fast_ace_workspace(
    tree_workspace, model = "ARD", nl = nlevels(f$y)
  )
  args <- list(
    tree = f$tree, y = f$y, lambda0 = 0.1, proposal_sd = 0.5,
    mcmc_sim = 12L, thin = 2L, burn = 2L, entropy_code = 1L,
    model = "ARD", ace_engine = "fast", keep_chains = TRUE
  )

  set.seed(2048)
  legacy_observed <- suppressWarnings(do.call(
    fastphylosig:::.delta_one, c(args, list(ace_workspace = NULL))
  ))
  set.seed(2048)
  cached_observed <- suppressWarnings(do.call(
    fastphylosig:::.delta_one, c(args, list(ace_workspace = workspace))
  ))
  testthat::expect_equal(cached_observed, legacy_observed, tolerance = 0)

  run_permutations <- function(ace_workspace) {
    vapply(seq_len(nrow(f$permutations)), function(i) {
      fastphylosig:::.delta_permutation_worker(
        i = i, y = f$y, perms = f$permutations, tree = f$tree,
        lambda0 = 0.1, proposal_sd = 0.5, mcmc_sim = 12L, thin = 2L,
        burn = 2L, entropy_code = 1L, model = "ARD", ace_engine = "fast",
        ace_workspace = ace_workspace
      )
    }, numeric(1))
  }
  set.seed(2049)
  legacy_permutations <- suppressWarnings(run_permutations(NULL))
  set.seed(2049)
  cached_permutations <- suppressWarnings(run_permutations(workspace))
  testthat::expect_equal(cached_permutations, legacy_permutations,
                         tolerance = 0)
})
