# Root-sensitivity contracts.  These tests deliberately distinguish a
# representation-equivalent reroot from a biological topology repair.  ape's
# resolve.root=TRUE can create a zero-length root branch for a non-current
# root; the helper below instead splits an existing positive edge and keeps
# its total length unchanged.

.root_contract_tree <- function(kind, n) {
  n <- as.integer(n)
  if (kind == "balanced") {
    size <- 2L^ceiling(log2(n))
    tree <- ape::stree(
      size, type = "balanced",
      tip.label = paste0("sp", seq_len(size))
    )
    if (size > n) {
      tree <- ape::drop.tip(tree, paste0("sp", seq.int(n + 1L, size)))
    }
  } else if (kind %in% c("pectinate", "strongly_unbalanced")) {
    tree <- ape::stree(
      n, type = "left", tip.label = paste0("sp", seq_len(n))
    )
  } else if (kind == "nonultrametric") {
    size <- 2L^ceiling(log2(n))
    tree <- ape::stree(
      size, type = "balanced",
      tip.label = paste0("sp", seq_len(size))
    )
    if (size > n) {
      tree <- ape::drop.tip(tree, paste0("sp", seq.int(n + 1L, size)))
    }
  } else {
    stop("unknown root-contract tree kind", call. = FALSE)
  }

  tree$tip.label <- paste0("sp", seq_len(ape::Ntip(tree)))
  tree$edge.length <- rep(1, nrow(tree$edge))
  if (kind == "strongly_unbalanced") {
    # Keep the topology pectinate while making root-to-tip lengths strongly
    # heterogeneous; all branches remain strictly positive.
    tree$edge.length <- seq(0.5, 2.5, length.out = nrow(tree$edge))
  } else if (kind == "nonultrametric") {
    tree$edge.length <- seq(0.35, 2.25, length.out = nrow(tree$edge))
  }
  ape::reorder.phylo(tree, "cladewise")
}

.root_contract_insert_root <- function(tree, node, root_edge = 1L,
                                        fraction = 0.5) {
  # First remove the conventional root without changing unrooted distances.
  unrooted <- ape::root(tree, node = node, resolve.root = FALSE)
  n <- ape::Ntip(unrooted)
  old_root <- n + 1L
  root_rows <- which(unrooted$edge[, 1L] == old_root)
  if (length(root_rows) != 3L || root_edge < 1L || root_edge > 3L) {
    stop("root construction did not expose a trichotomous root", call. = FALSE)
  }
  idx <- root_rows[[root_edge]]
  child <- unrooted$edge[idx, 2L]
  branch <- unrooted$edge.length[idx]
  if (!is.finite(branch) || branch <= 0 || fraction <= 0 || fraction >= 1) {
    stop("root construction requires a positive edge", call. = FALSE)
  }

  # Move the old root to a fresh internal id.  The new conventional root is
  # inserted halfway along the selected old-root edge, so the old root keeps
  # its other two children and both new edges are positive.
  old_internal <- 2L * n - 1L
  edge <- unrooted$edge
  edge[edge == old_root] <- old_internal
  edge[idx, ] <- c(n + 1L, child)
  edge <- rbind(c(n + 1L, old_internal), edge)
  edge_length <- c(fraction * branch, unrooted$edge.length)
  edge_length[idx + 1L] <- (1 - fraction) * branch

  out <- unrooted
  out$edge <- edge
  out$edge.length <- edge_length
  out$Nnode <- n - 1L
  out$root.edge <- NULL
  ape::reorder.phylo(out, "cladewise")
}

.root_contract_rootings <- function(tree) {
  n <- ape::Ntip(tree)
  # n + 2 is an internal node for every binary fixture above.  The two
  # selected root edges represent distinct legal placements on the same
  # unrooted branch network.
  list(
    root_a = .root_contract_insert_root(tree, n + 2L, root_edge = 1L),
    root_b = .root_contract_insert_root(tree, n + 2L, root_edge = 2L)
  )
}

.root_contract_continuous <- function(tree) {
  n <- ape::Ntip(tree)
  tips <- tree$tip.label
  stats::setNames(
    sin(seq_len(n) / 3) + seq_len(n) / n,
    tips
  )
}

.root_contract_binary <- function(tree) {
  n <- ape::Ntip(tree)
  set.seed(9137 + n)
  stats::setNames(as.integer(stats::runif(n) > 0.45), tree$tip.label)
}

.root_contract_caper <- function(tree, binary, nsim = 2L, seed = 9139L) {
  data <- data.frame(
    sp = tree$tip.label,
    trait = as.integer(binary[tree$tip.label]),
    stringsAsFactors = FALSE
  )
  comparative <- caper::comparative.data(
    tree, data, names.col = sp, warn.dropped = FALSE
  )
  set.seed(seed)
  reference <- suppressWarnings(caper::phylo.d(
    comparative, binvar = trait, permut = nsim
  ))

  # Reproduce caper's random/Brownian null states so the complete D result,
  # not only the observed contrast, is compared under identical draws.
  set.seed(seed)
  ds <- as.integer(binary[tree$tip.label])
  prop_state1 <- mean(ds == sort(unique(ds))[[1L]])
  random <- replicate(nsim, sample(ds))
  vcv <- unclass(caper::VCV.array(tree))
  brownian <- mvtnorm::rmvnorm(nsim, sigma = vcv)
  brownian <- as.data.frame(t(brownian))
  threshold <- apply(brownian, 2L, quantile, prop_state1)
  brownian <- sweep(brownian, 2L, threshold, "<")
  brownian <- matrix(
    as.numeric(brownian), nrow = length(ds), ncol = nsim
  )
  list(reference = reference, random = random, brownian = brownian)
}

test_that("legal root placements preserve topology and match K/lambda references", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")

  kinds <- c("balanced", "pectinate", "strongly_unbalanced",
             "nonultrametric")
  k_differences <- numeric()
  for (n in c(20L, 100L)) {
    for (kind in kinds) {
      tree <- .root_contract_tree(kind, n)
      roots <- .root_contract_rootings(tree)
      reference_dist <- ape::cophenetic.phylo(roots[[1L]])
      testthat::expect_equal(
        ape::cophenetic.phylo(roots[[2L]]), reference_dist,
        tolerance = 1e-12, info = paste(kind, n, "root distance")
      )
      for (root_name in names(roots)) {
        rooted <- roots[[root_name]]
        checked <- check_tree(rooted)
        testthat::expect_true(
          all(checked$ready_by_signal[c("K", "lambda", "D")]),
          info = paste(kind, n, root_name)
        )
        testthat::expect_true(all(rooted$edge.length > 0))

        x <- .root_contract_continuous(rooted)
        got_k <- suppressMessages(fast_k(
          rooted, x, verbose = FALSE, progress = FALSE
        ))
        ref_k <- suppressWarnings(phytools::phylosig(
          rooted, x, method = "K", test = FALSE, se = NULL
        ))
        testthat::expect_equal(
          as.numeric(got_k), as.numeric(ref_k), tolerance = 1e-6,
          info = paste(kind, n, root_name, "K")
        )

        got_lambda <- suppressMessages(fast_lambda(
          rooted, x, lambda_profile = FALSE,
          verbose = FALSE, progress = FALSE
        ))
        ref_lambda <- suppressWarnings(phytools::phylosig(
          rooted, x, method = "lambda", test = FALSE, se = NULL
        ))
        testthat::expect_equal(
          got_lambda$lambda, unname(ref_lambda$lambda), tolerance = 1e-5,
          info = paste(kind, n, root_name, "lambda")
        )
        testthat::expect_equal(
          got_lambda$logL, unname(ref_lambda$logL), tolerance = 5e-4,
          info = paste(kind, n, root_name, "lambda logL")
        )
      }
      x <- .root_contract_continuous(roots[[1L]])
      k_a <- as.numeric(suppressMessages(fast_k(
        roots[[1L]], x, verbose = FALSE, progress = FALSE
      )))
      k_b <- as.numeric(suppressMessages(fast_k(
        roots[[2L]], x, verbose = FALSE, progress = FALSE
      )))
      k_differences <- c(k_differences, abs(k_a - k_b))
    }
  }
  # Root placement is part of the K contract: the two legal roots need not
  # yield the same K, so the test must compare each root to its reference,
  # rather than silently asserting root invariance.
  testthat::expect_true(any(k_differences > 1e-8))
})

test_that("D root placements match caper observed contrasts under controlled nulls", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("caper")

  for (kind in c("balanced", "pectinate", "strongly_unbalanced",
                 "nonultrametric")) {
    tree <- .root_contract_tree(kind, 20L)
    roots <- .root_contract_rootings(tree)
    binary <- .root_contract_binary(roots[[1L]])
    for (root_name in names(roots)) {
      rooted <- roots[[root_name]]
      caper_bundle <- .root_contract_caper(rooted, binary, nsim = 2L)
      got <- suppressWarnings(suppressMessages(fast_d(
        rooted, binary[rooted$tip.label], test = TRUE, nsim = 2L,
        random_states = caper_bundle$random,
        brownian_states = caper_bundle$brownian,
        return_sim = TRUE, verbose = FALSE, progress = FALSE
      )))
      reference <- caper_bundle$reference
      testthat::expect_equal(
        got$DEstimate, unname(reference$DEstimate),
        tolerance = 1e-10, info = paste(kind, root_name, "observed")
      )
      testthat::expect_equal(got$Pval1, reference$Pval1, tolerance = 1e-12)
      testthat::expect_equal(got$Pval0, reference$Pval0, tolerance = 1e-12)
      testthat::expect_equal(
        got$Parameters$Observed, unname(reference$Parameters$Observed),
        tolerance = 1e-10
      )
      testthat::expect_equal(
        got$Parameters$MeanRandom, reference$Parameters$MeanRandom,
        tolerance = 1e-10
      )
      testthat::expect_equal(
        got$Parameters$MeanBrownian, reference$Parameters$MeanBrownian,
        tolerance = 1e-10
      )
      testthat::expect_equal(
        got$Permutations$random, unname(reference$Permutations$random),
        tolerance = 1e-10
      )
      testthat::expect_equal(
        got$Permutations$brownian, unname(reference$Permutations$brownian),
        tolerance = 1e-10
      )
      testthat::expect_equal(got$nsim_requested, 2L)
      testthat::expect_equal(got$nsim_successful_random, 2L)
      testthat::expect_equal(got$nsim_successful_brownian, 2L)
      testthat::expect_true(is.finite(got$DEstimate))
    }
  }
})

test_that("zero terminal and internal branches agree between check_tree and fast_* boundaries", {
  testthat::skip_if_not_installed("ape")
  zero_terminal <- ape::read.tree(
    text = "((a:0,b:1):1,(c:1,d:1):1);"
  )
  zero_internal <- ape::read.tree(
    text = "((a:1,b:1):0,(c:1,d:1):1);"
  )
  continuous <- c(a = 0.1, b = 0.2, c = 0.3, d = 0.4)
  binary <- c(a = 0L, b = 0L, c = 1L, d = 1L)
  categorical <- c(a = "a", b = "b", c = "a", d = "b")

  terminal_check <- check_tree(zero_terminal)
  testthat::expect_false(any(terminal_check$ready_by_signal))
  for (method in c("K", "lambda", "D", "Delta")) {
    testthat::expect_error(
      switch(
        method,
        K = fast_signal(zero_terminal, continuous, method = "K",
                        verbose = FALSE, progress = FALSE),
        lambda = fast_signal(zero_terminal, continuous, method = "lambda",
                             verbose = FALSE, progress = FALSE),
        D = fast_signal(zero_terminal, binary, method = "D", test = FALSE,
                        verbose = FALSE, progress = FALSE),
        Delta = fast_signal(
          zero_terminal, categorical, method = "Delta", test = FALSE,
          mcmc_sim = 4L, thin = 2L, burn = 2L,
          verbose = FALSE, progress = FALSE
        )
      ),
      "USER_ACTION_REQUIRED|zero|branch|check_tree",
      info = paste("zero terminal", method)
    )
  }

  internal_check <- check_tree(zero_internal)
  testthat::expect_false(internal_check$ready_by_signal[["K"]])
  testthat::expect_true(internal_check$ready_by_signal[["lambda"]])
  testthat::expect_false(internal_check$ready_by_signal[["D"]])
  testthat::expect_false(internal_check$ready_by_signal[["Delta"]])
  testthat::expect_error(
    fast_signal(zero_internal, continuous, method = "K",
                verbose = FALSE, progress = FALSE),
    "USER_ACTION_REQUIRED|zero|branch|check_tree"
  )
  testthat::expect_s3_class(
    suppressMessages(fast_lambda(
      zero_internal, continuous, lambda_profile = FALSE,
      verbose = FALSE, progress = FALSE
    )),
    "phylosig"
  )
  testthat::expect_error(
    fast_signal(zero_internal, binary, method = "D", test = FALSE,
                verbose = FALSE, progress = FALSE),
    "USER_ACTION_REQUIRED|zero|branch|check_tree"
  )
  testthat::expect_error(
    fast_signal(zero_internal, categorical, method = "Delta", test = FALSE,
                mcmc_sim = 4L, thin = 2L, burn = 2L,
                verbose = FALSE, progress = FALSE),
    "USER_ACTION_REQUIRED|zero|branch|check_tree"
  )
})

test_that("D controlled ties use caper strict tails and expose a failed denominator", {
  testthat::skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  observed <- c(a = 0, b = 0, c = 1, d = 1)
  states <- matrix(rep(observed[tree$tip.label], times = 3L),
                   nrow = ape::Ntip(tree), ncol = 3L,
                   dimnames = list(tree$tip.label, NULL))

  got <- suppressWarnings(suppressMessages(fast_d(
    tree, observed, test = TRUE, nsim = 3L,
    random_states = states, brownian_states = states,
    return_sim = TRUE, verbose = FALSE, progress = FALSE
  )))

  # caper::phylo.d() defines P_random as sum(random < observed)/nsim
  # and P_Brownian as sum(Brownian > observed)/nsim.  Equality is not
  # an exceedance, so both controlled tails are exactly zero here.
  testthat::expect_equal(got$P_random, 0, tolerance = 0)
  testthat::expect_equal(got$P_Brownian, 0, tolerance = 0)
  testthat::expect_equal(got$nsim_successful_random, 3L)
  testthat::expect_equal(got$nsim_successful_brownian, 3L)
  testthat::expect_true(is.na(got$DEstimate))
  testthat::expect_match(
    got$note,
    "random and Brownian null means are identical",
    fixed = TRUE
  )
  testthat::expect_equal(got$Permutations$random,
                         rep(got$Parameters$Observed, 3L),
                         tolerance = 1e-12)
  testthat::expect_equal(got$Permutations$brownian,
                         rep(got$Parameters$Observed, 3L),
                         tolerance = 1e-12)
})
