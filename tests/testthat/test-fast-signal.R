.test_permutation_matrix <- function(n, nsim) {
  out <- matrix(NA_integer_, nrow = nsim, ncol = n)
  out[1L, ] <- seq_len(n)
  if (nsim > 1L) {
    for (i in 2L:nsim) {
      out[i, ] <- sample.int(n)
    }
  }
  out
}

.test_k_with_perms <- function(tree, x, permutations) {
  C <- ape::vcv.phylo(tree)
  x <- x[rownames(C)]
  invC <- solve(C)
  n <- nrow(C)
  norm_const <- (sum(diag(C)) - n / sum(invC)) / (n - 1)
  k_fun <- function(y) {
    a <- sum(invC %*% y) / sum(invC)
    as.numeric(
      (t(y - a) %*% (y - a) / (t(y - a) %*% invC %*% (y - a))) /
        norm_const
    )
  }
  K <- k_fun(x)
  simK <- vapply(seq_len(nrow(permutations)), function(i) {
    k_fun(x[permutations[i, ]])
  }, numeric(1))
  list(K = K, P = mean(simK >= K), sim.K = simK)
}

test_that("K batch matches phytools for complete traits", {
  set.seed(10)
  tree <- ape::rtree(30)
  X <- matrix(
    stats::rnorm(30 * 6),
    nrow = 30,
    dimnames = list(tree$tip.label, paste0("trait_", 1:6))
  )

  res <- fast_signal(
    tree, X, method = "K", test = FALSE, verbose = FALSE
  )
  ref <- vapply(seq_len(ncol(X)), function(j) {
    as.numeric(phytools::phylosig(
      tree, X[, j], method = "K", test = FALSE, se = NULL
    ))
  }, numeric(1))

  expect_lt(max(abs(res$K_fast - ref), na.rm = TRUE), 1e-8)
})

test_that("fast_signal vector input returns a phylosig-like object", {
  set.seed(101)
  tree <- ape::rtree(30)
  x <- stats::rnorm(30)
  names(x) <- tree$tip.label

  fast <- fast_signal(tree, x, method = "K", test = FALSE,
                      verbose = FALSE)
  ref <- phytools::phylosig(tree, x, method = "K", test = FALSE,
                            se = NULL)

  expect_s3_class(fast, "phylosig")
  expect_equal(attr(fast, "method"), "K")
  expect_lt(abs(as.numeric(fast) - as.numeric(ref)), 1e-8)
})

test_that("K handles extra rows, missing tree tips, and trait-wise NAs", {
  set.seed(11)
  tree <- ape::rtree(35)
  kept <- tree$tip.label[-c(1, 2)]
  X <- matrix(
    stats::rnorm((length(kept) + 1) * 3),
    nrow = length(kept) + 1,
    dimnames = list(c(kept, "not_in_tree"), paste0("trait_", 1:3))
  )
  X[kept[3:5], 2] <- NA_real_

  res <- fast_signal(
    tree, X, method = "K", test = FALSE, verbose = FALSE
  )
  ref <- vapply(seq_len(ncol(X)), function(j) {
    x <- X[, j]
    names(x) <- rownames(X)
    as.numeric(phytools::phylosig(
      tree, x, method = "K", test = FALSE, se = NULL
    ))
  }, numeric(1))

  expect_equal(unique(res$removed_tree_tips), 2)
  expect_equal(unique(res$removed_data_rows), 1)
  expect_equal(res$n_removed_na[2], 3)
  expect_lt(max(abs(res$K_fast - ref), na.rm = TRUE), 1e-8)
})

test_that("match_tree_data reports matched and removed species", {
  set.seed(111)
  tree <- ape::rtree(12)
  kept <- tree$tip.label[-1]
  X <- data.frame(
    trait = stats::rnorm(length(kept) + 1),
    row.names = c(kept, "not_in_tree")
  )

  m <- match_tree_data(tree, X, verbose = FALSE)

  expect_s3_class(m$tree, "phylo")
  expect_equal(ape::Ntip(m$tree), 11)
  expect_equal(nrow(m$X), 11)
  expect_identical(m$data, m$X)
  expect_equal(m$report$matched_species, 11)
  expect_equal(m$report$removed_tree_tips, 1)
  expect_equal(m$report$removed_data_rows, 1)
  expect_equal(rownames(m$X), m$tree$tip.label)
  expect_equal(rownames(m$data), m$tree$tip.label)
})

test_that("match_phylo_data is a deprecated compatibility wrapper", {
  set.seed(112)
  tree <- ape::rtree(8)
  X <- data.frame(
    trait = stats::rnorm(8),
    row.names = tree$tip.label
  )

  canonical <- match_tree_data(tree, X, verbose = FALSE)
  warning_seen <- FALSE
  legacy <- withCallingHandlers(
    match_phylo_data(tree, X, verbose = FALSE),
    warning = function(w) {
      if (grepl("deprecated", conditionMessage(w), fixed = TRUE)) {
        warning_seen <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )
  expect_true(warning_seen)

  expect_equal(legacy, canonical)
  expect_identical(legacy$data, legacy$X)
})

test_that("K test path matches a shared permutation reference", {
  set.seed(12)
  tree <- ape::rtree(25)
  X <- matrix(
    stats::rnorm(25 * 4),
    nrow = 25,
    dimnames = list(tree$tip.label, paste0("trait_", 1:4))
  )
  perms <- .test_permutation_matrix(n = 25, nsim = 40)

  res <- fast_signal(
    tree, X, method = "K", test = TRUE, nsim = 40,
    permutations = perms, return_sim = TRUE, verbose = FALSE
  )
  refs <- lapply(seq_len(ncol(X)), function(j) {
    .test_k_with_perms(tree, X[, j], perms)
  })

  expect_lt(max(abs(res$K_fast - vapply(refs, `[[`, numeric(1), "K"))), 1e-8)
  expect_lt(max(abs(res$P_fast - vapply(refs, `[[`, numeric(1), "P"))), 1e-12)
  sim_err <- max(vapply(seq_len(ncol(X)), function(j) {
    max(abs(res$sim.K_fast[[j]] - refs[[j]]$sim.K))
  }, numeric(1)))
  expect_lt(sim_err, 1e-8)

  no_sim <- fast_signal(
    tree, X, method = "K", test = TRUE, nsim = 40,
    permutations = perms, return_sim = FALSE, verbose = FALSE
  )
  expect_false("sim.K_fast" %in% names(no_sim))
  expect_lt(max(abs(no_sim$K_fast - vapply(refs, `[[`, numeric(1), "K"))), 1e-8)
  expect_lt(max(abs(no_sim$P_fast - vapply(refs, `[[`, numeric(1), "P"))), 1e-12)
})

test_that("K randomization is identical with multiple C++ threads", {
  set.seed(120)
  tree <- ape::rtree(35)
  x <- stats::rnorm(35)
  names(x) <- tree$tip.label
  perms <- .test_permutation_matrix(n = 35, nsim = 80)

  one <- fast_signal(
    tree, x, method = "K", test = TRUE, nsim = 80,
    permutations = perms, return_sim = TRUE, verbose = FALSE, ncores = 1
  )
  two <- fast_signal(
    tree, x, method = "K", test = TRUE, nsim = 80,
    permutations = perms, return_sim = TRUE, verbose = FALSE, ncores = 2
  )

  expect_equal(one$K, two$K, tolerance = 1e-12)
  expect_equal(one$P, two$P, tolerance = 1e-12)
  expect_equal(one$sim.K, two$sim.K, tolerance = 1e-12)
})

test_that("custom permutation matrices are validated", {
  set.seed(121)
  tree <- ape::rtree(10)
  x <- stats::rnorm(10)
  names(x) <- tree$tip.label

  wrong_rows <- matrix(rep(seq_len(10), 2), nrow = 2, byrow = TRUE)
  expect_error(
    fast_signal(tree, x, method = "K", test = TRUE, nsim = 3,
                permutations = wrong_rows, verbose = FALSE),
    "exactly nsim rows"
  )

  duplicated_index <- matrix(rep(seq_len(10), 3), nrow = 3, byrow = TRUE)
  duplicated_index[1, 10] <- 9L
  expect_error(
    fast_signal(tree, x, method = "K", test = TRUE, nsim = 3,
                permutations = duplicated_index, verbose = FALSE),
    "every index"
  )

  fractional <- matrix(rep(seq_len(10), 3), nrow = 3, byrow = TRUE)
  fractional[1, 1] <- 1.5
  expect_error(
    fast_signal(tree, x, method = "K", test = TRUE, nsim = 3,
                permutations = fractional, verbose = FALSE),
    "finite integer"
  )
})

test_that("lambda batch matches phytools", {
  set.seed(13)
  tree <- ape::rtree(25)
  X <- matrix(
    stats::rnorm(25 * 3),
    nrow = 25,
    dimnames = list(tree$tip.label, paste0("trait_", 1:3))
  )

  res <- fast_signal(
    tree, X, method = "lambda", test = TRUE, verbose = FALSE
  )
  refs <- lapply(seq_len(ncol(X)), function(j) {
    phytools::phylosig(
      tree, X[, j], method = "lambda", test = TRUE, se = NULL
    )
  })
  ref_lambda <- vapply(refs, `[[`, numeric(1), "lambda")
  ref_logL <- vapply(refs, `[[`, numeric(1), "logL")
  ref_P <- vapply(refs, `[[`, numeric(1), "P")

  expect_lt(max(abs(res$lambda_fast - ref_lambda), na.rm = TRUE), 1e-5)
  # Upstream phytools changed its optimizer path in current R releases;
  # retain a reference tolerance that covers the documented floating-point
  # difference without weakening the production result contract.
  expect_lt(max(abs(res$logL_fast - ref_logL), na.rm = TRUE), 5e-4)
  expect_lt(max(abs(res$P_fast - ref_P), na.rm = TRUE), 1e-5)
})

test_that("phylo.d C++ contrast sums match caper contrCalc", {
  set.seed(21)
  tree <- ape::rtree(30)
  tree$node.label <- as.character(seq_len(tree$Nnode) + length(tree$tip.label))
  phy <- ape::reorder.phylo(tree, "pruningwise")
  states <- matrix(
    sample(0:1, 30 * 6, replace = TRUE),
    nrow = 30,
    ncol = 6,
    dimnames = list(tree$tip.label, c("Obs", paste0("V", 1:5)))
  )
  edge <- matrix(as.integer(phy$edge), ncol = 2)

  fast <- fastphylosig:::phylo_d_sums_cpp(
    states = states,
    edge = edge,
    edge_length = phy$edge.length,
    n_tip = length(phy$tip.label)
  )
  ref <- caper::contrCalc(
    vals = states,
    phy = phy,
    ref.var = "V1",
    picMethod = "phylo.d",
    crunch.brlen = 0
  )

  expect_lt(max(abs(as.numeric(fast) - colSums(ref$contrMat))), 1e-10)
})

test_that("fast_d supports vector and matrix inputs", {
  set.seed(22)
  tree <- ape::rtree(40)
  z <- ape::rTraitCont(tree)
  x <- as.integer(z > stats::median(z))
  names(x) <- tree$tip.label

  one <- fast_d(
    tree, x, nsim = 50, return_sim = FALSE, verbose = FALSE
  )
  expect_s3_class(one, "phylo.d")
  expect_true(is.finite(one$DEstimate))
  expect_true(is.finite(one$Pval1))
  expect_true(is.finite(one$Pval0))

  X <- data.frame(trait_a = x, trait_b = 1L - x, row.names = names(x))
  many <- fast_d(
    tree, X, nsim = 50, return_sim = FALSE, verbose = FALSE
  )
  expect_equal(nrow(many), 2)
  expect_true(all(is.finite(many$D_fast)))
})

test_that("fast_d gives identical controlled results with multiple C++ threads", {
  set.seed(220)
  tree <- ape::rtree(35)
  z <- ape::rTraitCont(tree)
  x <- as.integer(z > stats::median(z))
  names(x) <- tree$tip.label
  ds <- x[tree$tip.label]
  nsim <- 60
  random_states <- replicate(nsim, sample(ds))
  C <- unclass(caper::VCV.array(tree))
  samples <- t(chol(C)) %*% matrix(stats::rnorm(35 * nsim), nrow = 35)
  brownian_states <- fastphylosig:::brownian_threshold_cpp(
    samples, mean(ds == sort(unique(ds))[[1]])
  )

  one <- fast_d(
    tree, x, nsim = nsim, random_states = random_states,
    brownian_states = brownian_states, return_sim = TRUE,
    verbose = FALSE, ncores = 1
  )
  two <- fast_d(
    tree, x, nsim = nsim, random_states = random_states,
    brownian_states = brownian_states, return_sim = TRUE,
    verbose = FALSE, ncores = 2
  )

  expect_equal(one$DEstimate, two$DEstimate, tolerance = 1e-12)
  expect_equal(one$Pval1, two$Pval1, tolerance = 1e-12)
  expect_equal(one$Pval0, two$Pval0, tolerance = 1e-12)
  expect_equal(one$Permutations, two$Permutations, tolerance = 1e-12)
})

test_that("fast_d is invariant to binary state labels", {
  set.seed(222)
  tree <- ape::rtree(30)
  x01 <- rep(0:1, length.out = 30)
  names(x01) <- tree$tip.label
  x_labels <- ifelse(x01 == 0, 10, 20)
  names(x_labels) <- names(x01)
  nsim <- 40
  random_states <- replicate(nsim, sample(x01))
  brownian_states <- matrix(
    sample(0:1, 30 * nsim, replace = TRUE), nrow = 30, ncol = nsim
  )

  fit01 <- fast_d(
    tree, x01, nsim = nsim, random_states = random_states,
    brownian_states = brownian_states, return_sim = TRUE, verbose = FALSE
  )
  fit_labels <- fast_d(
    tree, x_labels, nsim = nsim, random_states = random_states,
    brownian_states = brownian_states, return_sim = TRUE, verbose = FALSE
  )

  expect_equal(fit01$DEstimate, fit_labels$DEstimate, tolerance = 1e-12)
  expect_equal(fit01$Pval1, fit_labels$Pval1, tolerance = 1e-12)
  expect_equal(fit01$Pval0, fit_labels$Pval0, tolerance = 1e-12)
})

test_that("direct tree Brownian threshold simulation has expected shape", {
  set.seed(221)
  tree <- ape::rtree(40)
  phy <- ape::reorder.phylo(tree, "cladewise")
  states <- fastphylosig:::brownian_tree_threshold_cpp(
    edge = matrix(as.integer(phy$edge), ncol = 2),
    edge_length = phy$edge.length,
    n_tip = 40,
    nsim = 20,
    prop_state1 = 0.5
  )

  expect_equal(dim(states), c(40L, 20L))
  expect_true(all(states %in% 0:1))
  expect_true(all(colSums(states) == 20))
})

test_that("fast_ace matches ape::ace for discrete ER and ARD likelihoods", {
  set.seed(31)
  tree <- ape::rtree(40)

  x_er <- ape::rTraitDisc(
    tree, model = "ER", k = 3, rate = 1.5, states = letters[1:3]
  )
  x_er <- factor(x_er[tree$tip.label])
  ref_er <- ape::ace(
    x_er, tree, type = "discrete", method = "ML", model = "ER"
  )
  fast_er <- fast_ace(x_er, tree, model = "ER")
  expect_lt(abs(ref_er$loglik - fast_er$loglik), 1e-8)
  expect_lt(max(abs(ref_er$rates - fast_er$rates)), 1e-6)
  se_error_er <- max(abs(ref_er$se - fast_er$se))
  se_scale_er <- max(1, abs(ref_er$se))
  # The SE is obtained from a numerical Hessian; BLAS/LAPACK and compiler
  # choices can perturb its last digits without changing the fitted model.
  expect_lt(se_error_er, 1e-7 + 1e-5 * se_scale_er)
  expect_lt(max(abs(ref_er$lik.anc - fast_er$lik.anc)), 1e-7)

  q <- matrix(
    c(0, 1.2, 0.25,
      0.45, 0, 1.8,
      0.15, 0.7, 0),
    nrow = 3, byrow = TRUE
  )
  x_ard <- ape::rTraitDisc(tree, model = q, states = letters[1:3])
  x_ard <- factor(x_ard[tree$tip.label])
  ref_ard <- suppressWarnings(ape::ace(
    x_ard, tree, type = "discrete", method = "ML", model = "ARD"
  ))
  fast_ard <- fast_ace(x_ard, tree, model = "ARD")
  expect_lt(abs(ref_ard$loglik - fast_ard$loglik), 1e-8)
  expect_lt(max(abs(ref_ard$rates - fast_ard$rates)), 1e-7)
  expect_lt(
    max(abs(ref_ard$se - fast_ard$se) / pmax(abs(ref_ard$se), 1)),
    5e-5
  )
  expect_lt(max(abs(ref_ard$lik.anc - fast_ard$lik.anc)), 1e-6)
})

test_that("fast_ace requires exact unique tip names", {
  set.seed(311)
  tree <- ape::rtree(12)
  x <- rep(c("a", "b"), length.out = 12)
  names(x) <- tree$tip.label

  bad <- x
  names(bad)[1] <- "not_in_tree"
  expect_error(fast_ace(bad, tree), "match phy\\$tip.label exactly")

  duplicated <- x
  names(duplicated)[1] <- names(duplicated)[2]
  expect_error(fast_ace(duplicated, tree), "must be unique")
})

test_that("Delta entropy and wrapper work for categorical traits", {
  prob <- matrix(
    c(0.8, 0.1, 0.1,
      0.34, 0.33, 0.33,
      0.2, 0.5, 0.3),
    ncol = 3, byrow = TRUE
  )
  ref_lse <- rowSums(ifelse(prob <= 1 / 3, prob,
                            prob / (1 - 3) - 1 / (1 - 3)))
  expect_equal(
    as.numeric(fastphylosig:::delta_entropy_cpp(prob, 1L)),
    ref_lse,
    tolerance = 1e-12
  )

  set.seed(24)
  tree <- ape::rtree(30)
  x <- rep(c("a", "b", "c"), length.out = 30)
  names(x) <- tree$tip.label

  one <- suppressWarnings(fast_delta(
    tree, x, mcmc_sim = 120, thin = 10, burn = 20, verbose = FALSE
  ))
  expect_s3_class(one, "phylo_delta")
  expect_true(is.finite(one$delta))

  X <- data.frame(
    trait_a = x,
    trait_b = rev(x),
    row.names = names(x),
    stringsAsFactors = FALSE
  )
  many <- suppressWarnings(fast_delta(
    tree, X, test = TRUE, nsim = 3, mcmc_sim = 80, thin = 10,
    burn = 20, return_sim = FALSE, verbose = FALSE
  ))
  expect_equal(nrow(many), 2)
  expect_true(all(is.finite(many$Delta_fast)))
  expect_true(all(is.finite(many$P_fast)))
})

test_that("plot_signal returns plotted data for supported results", {
  set.seed(241)
  tree <- ape::rtree(20)
  X <- matrix(
    stats::rnorm(20 * 2),
    nrow = 20,
    dimnames = list(tree$tip.label, c("a", "b"))
  )
  res <- fast_signal(tree, X, method = "K", test = TRUE, nsim = 10,
                     return_sim = TRUE, verbose = FALSE)
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  plotted <- plot_signal(res)
  grDevices::dev.off()

  expect_true(file.exists(f))
  expect_equal(nrow(plotted), 2)
  expect_equal(plotted$method, rep("K", 2))
  expect_true(all(is.finite(plotted$estimate)))
  expect_true(all(plotted$n_sim == 10))
  expect_true(all(lengths(plotted$sim) == 10))
  expect_equal(plotted$tail_direction, rep("right", 2))
})

test_that("plot_signal requires P values", {
  set.seed(2411)
  tree <- ape::rtree(20)
  x <- stats::rnorm(20)
  names(x) <- tree$tip.label
  fit <- fast_signal(tree, x, method = "K", test = FALSE,
                     verbose = FALSE)

  expect_error(plot_signal(fit), "return_sim")
})

test_that("plot_signal draws method-specific D and lambda plots", {
  set.seed(242)
  tree <- ape::rtree(30)
  z <- ape::rTraitCont(tree)
  x <- as.integer(z > stats::median(z))
  names(x) <- tree$tip.label

  fit <- fast_d(tree, x, nsim = 20, return_sim = TRUE, verbose = FALSE)
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  plotted <- plot_signal(fit)
  grDevices::dev.off()

  expect_true(file.exists(f))
  expect_equal(plotted$method, "D")
  expect_equal(plotted$plot_type, "D_calibration")
  expect_true(is.finite(plotted$estimate))
  expect_true(is.finite(plotted$n_sim))
  expect_true(length(plotted$sim_random[[1]]) == 20)
  expect_true(length(plotted$sim_brownian[[1]]) == 20)

  denom <- fit$Parameters$MeanRandom - fit$Parameters$MeanBrownian
  expected_random_d <- (fit$Permutations$random - fit$Parameters$MeanBrownian) /
    denom
  expected_brownian_d <- (fit$Permutations$brownian -
      fit$Parameters$MeanBrownian) / denom
  expect_equal(plotted$sim_random[[1]], expected_random_d, tolerance = 1e-12)
  expect_equal(plotted$sim_brownian[[1]], expected_brownian_d,
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(plotted$sim_random[[1]],
                                fit$Permutations$random)))
  expect_equal(plotted$extreme_random,
               sum(expected_random_d < fit$DEstimate))
  expect_equal(plotted$extreme_brownian,
               sum(expected_brownian_d > fit$DEstimate))
  expect_false(identical(plotted$P_random_display, "0"))
  expect_false(identical(plotted$P_Brownian_display, "0"))

  fit_lambda <- fast_signal(tree, z, method = "lambda", test = TRUE,
                            verbose = FALSE)
  expect_true(is.data.frame(fit_lambda$lambda_profile))
  expect_true(is.finite(fit_lambda$lambda_CI[["lower"]]))
  expect_true(is.finite(fit_lambda$lambda_CI[["upper"]]))
  lambda_hit <- which.min(abs(fit_lambda$lambda_profile$lambda -
                                fit_lambda$lambda))
  expect_equal(fit_lambda$lambda_profile$logL[[lambda_hit]], fit_lambda$logL,
               tolerance = 1e-10)
  f_lambda <- tempfile(fileext = ".png")
  grDevices::png(f_lambda)
  plotted_lambda <- plot_signal(fit_lambda)
  grDevices::dev.off()
  expect_true(file.exists(f_lambda))
  expect_equal(plotted_lambda$method, "lambda")
  expect_equal(plotted_lambda$plot_type, "lambda_profile")
  expect_true(is.data.frame(plotted_lambda$profile[[1]]))
  expect_true(is.finite(plotted_lambda$LR))
})

test_that("fast_delta permutation test can use multiple workers", {
  set.seed(240)
  tree <- ape::rtree(20)
  x <- rep(c("a", "b", "c"), length.out = 20)
  names(x) <- tree$tip.label

  res <- suppressWarnings(fast_delta(
    tree, x, test = TRUE, nsim = 2, mcmc_sim = 60, thin = 10,
    burn = 20, return_sim = TRUE, verbose = FALSE, ncores = 2
  ))

  expect_s3_class(res, "phylo_delta")
  expect_true(is.finite(res$delta))
  expect_true(is.finite(res$P) || is.na(res$P))
  expect_equal(length(res$sim.delta), 2)
})

test_that("fast_signal dispatches to the specialist entry points", {
  set.seed(23)
  tree <- ape::rtree(30)
  y <- stats::rnorm(30)
  names(y) <- tree$tip.label
  z <- as.integer(y > stats::median(y))
  names(z) <- tree$tip.label
  cat_trait <- rep(c("a", "b", "c"), length.out = 30)
  names(cat_trait) <- tree$tip.label

  k <- fast_signal(tree, y, method = "K", test = FALSE, verbose = FALSE)
  lambda <- fast_signal(tree, y, method = "lambda", test = FALSE,
                        verbose = FALSE)
  d <- fast_d(tree, z, nsim = 20, return_sim = FALSE, verbose = FALSE)
  delta <- suppressWarnings(fast_delta(
    tree, cat_trait, mcmc_sim = 100, thin = 10, burn = 20, verbose = FALSE
  ))

  expect_s3_class(k, "phylosig")
  expect_s3_class(lambda, "phylosig")
  expect_s3_class(d, "phylo.d")
  expect_s3_class(delta, "phylo_delta")
  high_d <- fast_signal(
    tree, z, method = "D", nsim = 20, return_sim = FALSE,
    verbose = FALSE
  )
  expect_s3_class(high_d, "phylo.d")
  expect_identical(attr(high_d, "workflow")$production_function, "fast_d")
})
