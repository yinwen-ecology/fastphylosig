test_that("odd middle observations and folded ties retain independent diagnostics", {
  for (n in c(5L, 8L, 31L, 32L, 991L)) for (chains in c(1L, 2L, 4L)) {
    a <- matrix(rep(c(1, 2, 2, 4, 8, 8, 16), length.out = n*chains), n, chains)
    b <- matrix(rep(c(16, 8, 4, 4, 2, 1, 1), length.out = n*chains), n, chains)
    if (n %% 2L) {
      a[n%/%2L + 1L, ] <- rep(c(.01, 100), length.out = chains)
      b[n%/%2L + 1L, ] <- rep(c(100, .01), length.out = chains)
    }
    native <- fastphylosig:::delta_diagnostics_cpp(a, b)
    reference <- fastphylosig:::.delta_diagnostics_fallback(a, b)
    for (nm in fastphylosig:::.delta_rank_fields())
      expect_equal(native[[nm]], reference[[nm]], tolerance = 1e-9,
        info = paste(n, chains, nm))
  }
})

test_that("burn and thinning endpoints match the independent posterior sampler", {
  p <- matrix(c(.2, .8, .35, .65, .6, .4, .8, .2), 4L, 2L, byrow = TRUE)
  controls <- rbind(c(1L, 1L, 1L), c(2L, 1L, 1L), c(17L, 3L, 4L),
    c(17L, 9L, 17L), c(17L, 1000000L, 1L), c(1000L, 7L, 101L))
  x <- 1 - rowSums(p*p)
  x <- x*2
  for (i in seq_len(nrow(controls))) {
    z <- controls[i, ]
    set.seed(28640000L)
    native <- fastphylosig:::delta_mcmc_cpp(p, .1, .5, z[1L], z[2L], z[3L], 3L, TRUE)
    seed <- .Random.seed
    set.seed(28640000L)
    reference <- .delta_posterior_reference(x, .1, .5, z[1L], z[2L], z[3L])
    expect_equal(unname(native$alpha_chain), reference$alpha, tolerance = 1e-10)
    expect_equal(unname(native$beta_chain), reference$beta, tolerance = 1e-10)
    expect_identical(.Random.seed, seed)
    expect_equal(nrow(native$alpha_chain), 1L + (z[1L] - z[3L]) %/% z[2L])
  }
})
