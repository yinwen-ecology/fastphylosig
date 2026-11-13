test_that("ACE prepared deviance preserves stateless fits including SE", {
  set.seed(101)
  tree <- ape::rtree(48)
  for (k in c(2L, 3L)) {
    x <- ape::rTraitDisc(tree, model = "ER", k = k, rate = .3,
                        states = LETTERS[seq_len(k)])
    x <- droplevels(factor(x[tree$tip.label]))
    if (nlevels(x) < 2L) next
    for (model in c("ER", "ARD")) for (scale in c(1e-9, 1, 1e3)) {
      scaled <- tree
      scaled$edge.length <- tree$edge.length * scale
      tw <- .fast_ace_tree_workspace(scaled)
      ws <- .fast_ace_workspace(tw, model, nlevels(x))
      # The native workspace is an implementation optimization, not a new
      # optimizer. Both paths must retain the full result and diagnostics.
      fit <- function(reuse) suppressWarnings(.fast_ace_fit_workspace(
        as.integer(x), levels(x), ws, ip = .1,
        CI = TRUE, estimate_se = TRUE, reuse_buffers = reuse))
      before <- .Random.seed
      prepared <- fit(TRUE)
      expect_identical(.Random.seed, before)
      reference <- fit(FALSE)
      expect_equal(prepared, reference, tolerance = 1e-10)
    }
  }
})