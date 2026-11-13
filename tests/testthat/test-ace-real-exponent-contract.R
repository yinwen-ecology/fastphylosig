test_that("real-exponent specialization retains numerical ACE equivalence", {
  set.seed(9008)
  for(k in c(2L,3L,5L,8L))for(model in c("ER","ARD"))
  for(repeated in c(FALSE,TRUE)) {
    tree <- ape::rtree(64L)
    if(repeated)tree$edge.length <- rep(c(.05,.1,.5,2),length.out=nrow(tree$edge))
    w <- .fast_ace_workspace(.fast_ace_tree_workspace(tree),model,k)
    states <- rep(seq_len(k),length.out=64L);states[c(3,17)] <- NA_integer_
    ptr <- fast_ace_fit_prepare_cpp(w$edge,w$edge_length,states,w$rate$rate_index,1L)
    params <- list(rep(.1,w$rate$np),runif(w$rate$np,.01,.4),
      rep(1000,w$rate$np),rep(0,w$rate$np),rep(.1,w$rate$np))
    for(y in list(states,rev(states))) {
      fast_ace_fit_update_cpp(ptr,y)
      for(par in params) {
        seed <- .Random.seed
        old <- fast_ace_discrete_deviance_cpp(w$edge,w$edge_length,y,
          w$rate$rate_index,par,FALSE)
        new <- fast_ace_fit_deviance_cpp(ptr,par)
        .expect_ace_numerical_equivalent(as.numeric(new), as.numeric(old))
        expect_identical(.Random.seed,seed)
        if(is.finite(new)&&new<1e49)for(marginal in c(FALSE,TRUE)) {
          old_liks <- fast_ace_discrete_liks_cpp(w$edge,w$edge_length,y,
            w$rate$rate_index,par,marginal,FALSE)
          new_liks <- fast_ace_fit_liks_cpp(ptr,par,marginal)
          .expect_ace_numerical_equivalent(new_liks, old_liks)
          expect_identical(.Random.seed,seed)
        }
      }
    }
    fast_ace_fit_release_cpp(ptr)
  }
})


test_that("native ACE paths agree with independent matrix-exponential likelihoods", {
  skip_if_not_installed("expm")
  # Fixed rates test mathematical values, independently of optimizer trajectories.
  # Include real/complex eigensystems, large rates, missing states and the
  # small-cache / serial-or-batched boundary. No optimized native result is the
  # oracle: .ace_math_loglik performs R pruning with independent expm transitions.
  for (k in c(2L, 3L, 5L, 8L)) {
    tree <- ape::read.tree(text =
      "(((a:0.2,b:0.7):0.5,(c:0.3,d:0.9):0.4):0.2,((e:0.5,f:0.4):0.3,(g:0.7,h:0.2):0.6):0.8);")
    w <- .fast_ace_workspace(.fast_ace_tree_workspace(tree), "ARD", k)
    states <- rep(seq_len(k), length.out = 8L)
    states[8L] <- NA_integer_
    local({
      ptr <- fast_ace_fit_prepare_cpp(w$edge, w$edge_length, states,
        w$rate$rate_index, 1L)
      on.exit(fast_ace_fit_release_cpp(ptr))
      for (rates in list(rep(.1, w$rate$np),
                         seq_len(w$rate$np) * .01, rep(1000, w$rate$np))) {
        truth <- .ace_math_loglik(tree, states,
          .ace_math_q(w$rate$rate_index, rates)) + log(k)
        expect_true(is.finite(truth))
        for (reuse in c(FALSE, TRUE)) {
          got <- fast_ace_discrete_liks_cpp(w$edge, w$edge_length, states,
            w$rate$rate_index, rates, FALSE, reuse)
          expect_equal(got$loglik, truth, tolerance = 1e-10)
        }
        got <- fast_ace_fit_liks_cpp(ptr, rates, FALSE)
        expect_equal(got$loglik, truth, tolerance = 1e-10)
      }
    })
  }
  Q <- matrix(c(0, .4, .01, .02, 0, .5, .6, .03, 0), 3L, byrow = TRUE)
  diag(Q) <- -rowSums(Q)
  rates <- Q[row(Q) != col(Q)]
  for (n in c(32L, 33L, 1366L, 1367L)) {
    set.seed(906 + n)
    tree <- ape::rtree(n)
    tree$edge.length <- seq(.05, 2, length.out = nrow(tree$edge))
    states <- rep(1:3, length.out = n); states[7L] <- NA_integer_
    w <- .fast_ace_workspace(.fast_ace_tree_workspace(tree), "ARD", 3L)
    truth <- .ace_math_loglik(tree, states, Q) + log(3)
    expect_true(is.finite(truth))
    for (reuse in c(FALSE, TRUE)) {
      got <- fast_ace_discrete_liks_cpp(w$edge, w$edge_length, states,
        w$rate$rate_index, rates, FALSE, reuse)
      expect_equal(got$loglik, truth, tolerance = 1e-10)
    }
    for (threads in c(1L, 2L)) local({
      ptr <- fast_ace_fit_prepare_cpp(w$edge, w$edge_length, states,
        w$rate$rate_index, threads)
      on.exit(fast_ace_fit_release_cpp(ptr))
      got <- fast_ace_fit_liks_cpp(ptr, rates, FALSE)
      expect_equal(got$loglik, truth, tolerance = 1e-10)
    })
  }
})
