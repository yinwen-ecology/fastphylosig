# Mathematical-oracle checks for the in-place ACE transition-vector kernel.
#
# The production implementation and this oracle intentionally use the same
# continuous-time Markov-chain definition, but the oracle performs each
# transition with base-R eigen algebra.  Keeping the reference independent of
# the C++ expression-template path catches both allocation-oriented rewrites
# and accidental changes to pruning/reconstruction semantics.

.ace_opt_rate_index <- function(model, k) {
  out <- matrix(0L, k, k)
  if (identical(model, "ER")) {
    out[row(out) != col(out)] <- 1L
  } else {
    out[row(out) != col(out)] <- seq_len(k * (k - 1L))
  }
  out
}

.ace_opt_oracle <- function(edge, edge_length, tip_state, rate_index, par,
                            marginal = FALSE) {
  n_tip <- length(tip_state)
  n_node <- n_tip - 1L
  n_total <- n_tip + n_node
  k <- nrow(rate_index)

  Q <- matrix(0, k, k)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (i != j && rate_index[i, j] > 0L) {
        Q[i, j] <- par[rate_index[i, j]]
      }
    }
    Q[i, i] <- -sum(Q[i, ])
  }
  eg <- eigen(Q)
  inv_v <- solve(eg$vectors)
  transition <- lapply(edge_length, function(t) {
    Re(eg$vectors %*% diag(exp(eg$values * t), k, k) %*% inv_v)
  })

  liks <- matrix(0, n_total, k)
  for (i in seq_len(n_tip)) {
    liks[i, tip_state[i]] <- 1
  }
  log_comp <- 0
  for (e in seq.int(1L, nrow(edge), by = 2L)) {
    anc <- edge[e, 1L]
    des_left <- edge[e, 2L]
    des_right <- edge[e + 1L, 2L]
    vl <- transition[[e]] %*% liks[des_left, ]
    vr <- transition[[e + 1L]] %*% liks[des_right, ]
    v <- as.numeric(vl) * as.numeric(vr)
    comp <- sum(v)
    liks[anc, ] <- v / comp
    log_comp <- log_comp + log(comp)
  }

  lik_anc <- liks[(n_tip + 1L):n_total, , drop = FALSE]
  if (!marginal) {
    for (e in rev(seq.int(1L, nrow(edge), by = 2L))) {
      anc <- edge[e, 1L] - n_tip
      des_left <- edge[e, 2L] - n_tip
      if (des_left > 0L) {
        denom <- as.numeric(lik_anc[des_left, , drop = FALSE] %*%
                              transition[[e]])
        tmp <- lik_anc[anc, ] / denom
        lik_anc[des_left, ] <-
          as.numeric(tmp %*% transition[[e]]) * lik_anc[des_left, ]
      }

      des_right <- edge[e + 1L, 2L] - n_tip
      if (des_right > 0L) {
        denom <- as.numeric(lik_anc[des_right, , drop = FALSE] %*%
                              transition[[e + 1L]])
        tmp <- lik_anc[anc, ] / denom
        lik_anc[des_right, ] <-
          as.numeric(tmp %*% transition[[e + 1L]]) * lik_anc[des_right, ]
      }

      rowsums <- rowSums(lik_anc)
      for (r in seq_len(nrow(lik_anc))) {
        if (rowsums[r] != 0) lik_anc[r, ] <- lik_anc[r, ] / rowsums[r]
      }
    }
  }

  list(deviance = -2 * log_comp, lik.anc = lik_anc)
}

.ace_opt_payload <- function(k, model, scale, seed) {
  set.seed(seed)
  tree <- ape::rtree(8L)
  tree <- ape::reorder.phylo(tree, "postorder")
  tree$edge.length <- tree$edge.length * scale
  edge <- matrix(as.integer(tree$edge), ncol = 2L)
  tip_state <- rep(seq_len(k), length.out = ape::Ntip(tree))
  rate_index <- .ace_opt_rate_index(model, k)
  npar <- max(rate_index)
  par <- if (identical(model, "ER")) {
    0.7
  } else {
    0.15 + 0.025 * seq_len(npar)
  }
  list(edge = edge, edge_length = tree$edge.length,
       tip_state = as.integer(tip_state), rate_index = rate_index, par = par)
}

test_that("in-place ACE transitions preserve ER/ARD likelihoods", {
  skip_if_not_installed("ape")

  # Include two independently generated postorder topologies and both short
  # and long branches.  The direct fixed-rate call keeps this test focused on
  # the transition/pruning kernel rather than optimizer tolerances.
  cases <- expand.grid(
    k = c(2L, 3L, 5L), model = c("ER", "ARD"),
    scale = c(0.01, 1, 10), seed = c(41L, 73L),
    stringsAsFactors = FALSE
  )

  for (row in seq_len(nrow(cases))) {
    p <- .ace_opt_payload(
      k = cases$k[row], model = cases$model[row],
      scale = cases$scale[row], seed = cases$seed[row]
    )
    ref <- .ace_opt_oracle(
      p$edge, p$edge_length, p$tip_state, p$rate_index, p$par
    )
    got_dev <- fastphylosig:::fast_ace_discrete_deviance_cpp(
      p$edge, p$edge_length, p$tip_state, p$rate_index, p$par
    )
    got <- fastphylosig:::fast_ace_discrete_liks_cpp(
      p$edge, p$edge_length, p$tip_state, p$rate_index, p$par,
      marginal = FALSE
    )

    expect_equal(got_dev, ref$deviance, tolerance = 1e-8,
                 info = paste(cases$model[row], cases$k[row],
                              cases$scale[row], cases$seed[row]))
    expect_equal(got$loglik, -0.5 * ref$deviance, tolerance = 1e-8,
                 info = paste(cases$model[row], cases$k[row],
                              cases$scale[row], cases$seed[row]))
    expect_equal(got$lik.anc, ref$lik.anc, tolerance = 1e-8,
                 info = paste(cases$model[row], cases$k[row],
                              cases$scale[row], cases$seed[row]))
  }
})

test_that("ACE transition/pruning calls are repeatable", {
  skip_if_not_installed("ape")
  p <- .ace_opt_payload(k = 5L, model = "ARD", scale = 1, seed = 101L)
  one <- fastphylosig:::fast_ace_discrete_liks_cpp(
    p$edge, p$edge_length, p$tip_state, p$rate_index, p$par,
    marginal = FALSE
  )
  for (i in seq_len(8L)) {
    again <- fastphylosig:::fast_ace_discrete_liks_cpp(
      p$edge, p$edge_length, p$tip_state, p$rate_index, p$par,
      marginal = FALSE
    )
    expect_equal(again$deviance, one$deviance, tolerance = 1e-12)
    expect_equal(again$loglik, one$loglik, tolerance = 1e-12)
    expect_equal(again$lik.anc, one$lik.anc, tolerance = 1e-12)
  }
})

