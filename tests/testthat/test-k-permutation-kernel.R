.k_inclusive_upper_tail <- function(value, observed) {
  direct <- value >= observed
  direct[is.na(direct)] <- FALSE
  finite <- is.finite(value) & is.finite(observed)
  scale <- pmax(1, abs(value), abs(observed))
  tie <- finite & (observed - value <=
    8 * .Machine$double.eps * scale)
  direct | tie
}

test_that("bounded K permutation kernel matches dense controlled permutations", {
  skip_if_not(exists("fast_k_tree_permutation_cpp",
                     asNamespace("fastphylosig")))

  tree <- ape::read.tree(
    text = "((a:0.5,b:0.7):0.8,(c:0.4,d:0.6):0.9);"
  )
  d_tree <- ape::reorder.phylo(tree, "pruningwise")
  compiled <- fastphylosig:::compile_tree_cpp(
    d_tree$edge, d_tree$edge.length, ape::Ntip(d_tree)
  )
  x <- matrix(c(1.2, -0.7, 2.3, 0.4), ncol = 1,
              dimnames = list(tree$tip.label, "trait"))
  perms <- rbind(
    c(1L, 2L, 3L, 4L), c(2L, 1L, 4L, 3L),
    c(3L, 4L, 1L, 2L), c(4L, 3L, 2L, 1L)
  )

  got <- fastphylosig:::fast_k_tree_permutation_cpp(
    compiled, x, nsim = nrow(perms), permutations = perms,
    trait_chunk = 1L, return_sim = TRUE, n_threads = 1L
  )
  dense <- fastphylosig:::fast_k_chol_permutation_cpp(
    x, chol(ape::vcv.phylo(d_tree)), sum(diag(ape::vcv.phylo(d_tree))),
    perms, n_threads = 1L
  )

  expect_equal(as.numeric(got$K), as.numeric(dense$K), tolerance = 1e-8)
  expect_equal(as.numeric(got$sim_K), as.numeric(dense$sim_K),
               tolerance = 1e-8)
  expected_exceedance <- sum(.k_inclusive_upper_tail(
    as.numeric(got$sim_K), as.numeric(got$K)
  ))
  expect_equal(as.numeric(got$P), expected_exceedance / nrow(perms),
               tolerance = 0)
  expect_equal(as.numeric(dense$P), expected_exceedance / nrow(perms),
               tolerance = 0)
  expect_equal(got$nsim_requested, nrow(perms))
  expect_equal(got$nsim_successful, nrow(perms))
  expect_identical(got$permutation_mode, "controlled")
  expect_false(got$include_observed)
})

test_that("finite controlled K grids preserve the inclusive tail and P lower bound", {
  skip_if_not(exists("fast_k_tree_permutation_cpp",
                     asNamespace("fastphylosig")))

  tree <- ape::read.tree(text = "(((a:0.5,b:0.7):0.8,c:0.6):0.9,(d:0.4,e:0.6):1.0);")
  d_tree <- ape::reorder.phylo(tree, "pruningwise")
  compiled <- fastphylosig:::compile_tree_cpp(
    d_tree$edge, d_tree$edge.length, ape::Ntip(d_tree)
  )
  n <- ape::Ntip(d_tree)
  x <- matrix(seq_len(n) / 3, nrow = n, ncol = 1,
              dimnames = list(d_tree$tip.label, "trait"))
  chol_c <- chol(ape::vcv.phylo(d_tree))
  trace_c <- sum(diag(ape::vcv.phylo(d_tree)))

  make_permutations <- function(nsim) {
    do.call(rbind, lapply(seq_len(nsim), function(i) {
      if (i == 1L) return(seq_len(n))
      if (i %% 3L == 0L) return(rev(seq_len(n)))
      shift <- (i - 1L) %% n
      c((shift + 1L):n, seq_len(shift))
    }))
  }

  for (nsim in c(5L, 19L, 99L)) {
    perms <- make_permutations(nsim)
    fast <- fastphylosig:::fast_k_tree_permutation_cpp(
      compiled, x, nsim = nsim, permutations = perms,
      trait_chunk = 1L, return_sim = TRUE, n_threads = 1L
    )
    dense <- fastphylosig:::fast_k_chol_permutation_cpp(
      x, chol_c, trace_c, perms, n_threads = 1L
    )
    expect_equal(as.numeric(fast$K), as.numeric(dense$K), tolerance = 1e-10,
                 info = paste("observed", nsim))
    expect_equal(as.numeric(fast$sim_K), as.numeric(dense$sim_K), tolerance = 1e-10,
                 info = paste("simulated", nsim))
    expected_exceedance <- sum(.k_inclusive_upper_tail(
      as.numeric(fast$sim_K), as.numeric(fast$K)
    ))
    expect_equal(fast$exceedance_count, expected_exceedance)
    expect_equal(fast$P, expected_exceedance / nsim, tolerance = 0)
    expect_equal(as.numeric(dense$P), expected_exceedance / nsim,
                 tolerance = 0)
    expect_true(fast$P >= 1 / nsim,
                info = paste("identity permutation lower bound", nsim))
    expect_equal(fast$nsim_requested, nsim)
    expect_equal(fast$nsim_successful, nsim)
  }
})

test_that("controlled K ties use the inclusive >= tail", {
  skip_if_not(exists("fast_k_tree_permutation_cpp",
                     asNamespace("fastphylosig")))
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  d_tree <- ape::reorder.phylo(tree, "pruningwise")
  compiled <- fastphylosig:::compile_tree_cpp(
    d_tree$edge, d_tree$edge.length, ape::Ntip(d_tree)
  )
  x <- matrix(c(1, 2, 4, 8), nrow = 4, ncol = 1,
              dimnames = list(d_tree$tip.label, "trait"))
  perms <- matrix(rep(seq_len(4L), times = 5L), nrow = 5L, byrow = TRUE)
  got <- fastphylosig:::fast_k_tree_permutation_cpp(
    compiled, x, nsim = 5L, permutations = perms,
    return_sim = TRUE, n_threads = 1L
  )
  expect_true(all(abs(got$sim_K - got$K) < 1e-12))
  expect_equal(got$exceedance_count, 5)
  expect_equal(got$P, 1, tolerance = 0)
})

test_that("internal K permutation mode is streamed and reports MCSE", {
  skip_if_not(exists("fast_k_tree_permutation_cpp",
                     asNamespace("fastphylosig")))

  set.seed(17)
  tree <- ape::rtree(12)
  d_tree <- ape::reorder.phylo(tree, "pruningwise")
  compiled <- fastphylosig:::compile_tree_cpp(
    d_tree$edge, d_tree$edge.length, ape::Ntip(d_tree)
  )
  x <- matrix(stats::rnorm(24), nrow = 12, ncol = 2,
              dimnames = list(tree$tip.label, c("a", "b")))
  got <- fastphylosig:::fast_k_tree_permutation_cpp(
    compiled, x, nsim = 32L, trait_chunk = 1L,
    return_sim = FALSE, include_observed = TRUE
  )

  expect_false("sim_K" %in% names(got))
  expect_equal(got$nsim_requested, 32L)
  expect_equal(got$nsim_successful, c(32L, 32L))
  expect_equal(got$n_randomizations, 31L)
  expect_identical(got$permutation_mode, "internal_rng")
  expect_true(all(is.finite(got$MCSE_P)))
  expect_true(all(got$P >= 0 & got$P <= 1))
})

test_that("internal RNG chunks are reproducible across chunk sizes and threads", {
  skip_if_not(exists("fast_k_tree_permutation_cpp",
                     asNamespace("fastphylosig")))

  set.seed(17)
  tree <- ape::rtree(14)
  d_tree <- ape::reorder.phylo(tree, "pruningwise")
  compiled <- fastphylosig:::compile_tree_cpp(
    d_tree$edge, d_tree$edge.length, ape::Ntip(d_tree)
  )
  x <- matrix(stats::rnorm(28), nrow = 14, ncol = 2,
              dimnames = list(tree$tip.label, c("a", "b")))

  set.seed(123)
  one <- fastphylosig:::fast_k_tree_permutation_cpp(
    compiled, x, nsim = 64L, trait_chunk = 1L,
    return_sim = TRUE, include_observed = TRUE,
    n_threads = 1L, simulation_chunk = 7L
  )
  set.seed(123)
  two <- fastphylosig:::fast_k_tree_permutation_cpp(
    compiled, x, nsim = 64L, trait_chunk = 1L,
    return_sim = TRUE, include_observed = TRUE,
    n_threads = 1L, simulation_chunk = 64L
  )
  expect_equal(one$sim_K, two$sim_K, tolerance = 1e-12)
  expect_equal(one$P, two$P, tolerance = 0)

  set.seed(123)
  threaded <- fastphylosig:::fast_k_tree_permutation_cpp(
    compiled, x, nsim = 64L, trait_chunk = 1L,
    return_sim = TRUE, include_observed = TRUE,
    n_threads = 2L, simulation_chunk = 7L
  )
  # Worker floating-point environments can differ from the main R thread by
  # a few ulps, but the generated permutations, exceedance counts, and P are
  # deterministic.
  expect_equal(one$sim_K, threaded$sim_K, tolerance = 1e-12)
  expect_equal(one$P, threaded$P, tolerance = 0)
  expect_equal(one$simulation_chunk, 7L)

  expect_error(
    fastphylosig:::fast_k_tree_permutation_cpp(
      compiled, x, nsim = 4L, simulation_chunk = 0L
    ),
    "simulation_chunk must be a positive integer"
  )
})

test_that("constant observed traits do not produce artificial K P values", {
  skip_if_not(exists("fast_k_tree_permutation_cpp",
                     asNamespace("fastphylosig")))

  tree <- ape::rtree(6)
  d_tree <- ape::reorder.phylo(tree, "pruningwise")
  compiled <- fastphylosig:::compile_tree_cpp(
    d_tree$edge, d_tree$edge.length, ape::Ntip(d_tree)
  )
  x <- matrix(1, nrow = 6, ncol = 1,
              dimnames = list(tree$tip.label, "constant"))
  got <- fastphylosig:::fast_k_tree_permutation_cpp(
    compiled, x, nsim = 8L, return_sim = FALSE
  )

  expect_true(is.na(got$K[[1L]]))
  expect_true(is.na(got$P[[1L]]))
  expect_equal(got$nsim_successful[[1L]], 0L)
})

test_that("K permutation contrasts remain stable under large common offsets", {
  skip_if_not(exists("fast_k_tree_permutation_cpp",
                     asNamespace("fastphylosig")))

  set.seed(19)
  tree <- ape::rtree(25)
  d_tree <- ape::reorder.phylo(tree, "pruningwise")
  compiled <- fastphylosig:::compile_tree_cpp(
    d_tree$edge, d_tree$edge.length, ape::Ntip(d_tree)
  )
  # Keep the trait contrast scale large enough that adding 1e15 does not
  # erase it at double precision; the test then isolates the kernel's
  # baseline-centred arithmetic rather than input quantization.
  x <- matrix(round(stats::rnorm(25) * 1e8, 3), nrow = 25, ncol = 1,
              dimnames = list(tree$tip.label, "trait"))
  perms <- t(replicate(40L, sample.int(25L)))
  unshifted <- fastphylosig:::fast_k_tree_permutation_cpp(
    compiled, x, nsim = nrow(perms), permutations = perms,
    trait_chunk = 1L, return_sim = TRUE
  )

  for (offset in c(1e12, 1e15)) {
    shifted <- fastphylosig:::fast_k_tree_permutation_cpp(
      compiled, x + offset, nsim = nrow(perms), permutations = perms,
      trait_chunk = 1L, return_sim = TRUE
    )
    expect_equal(shifted$K, unshifted$K, tolerance = 1e-8)
    expect_equal(shifted$sim_K, unshifted$sim_K, tolerance = 1e-7)
    expect_equal(shifted$P, unshifted$P, tolerance = 0)
  }
})
