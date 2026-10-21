test_that("exact eigenvalue reuse preserves changing generators and state data", {
  set.seed(20261020)
  tr <- ape::rtree(96)
  tr$edge.length <- exp(runif(nrow(tr$edge),-2,1))
  w <- .fast_ace_workspace(.fast_ace_tree_workspace(tr),"ARD",3L)
  states <- rep(1:3,length.out=96)
  ptr <- fast_ace_fit_prepare_cpp(w$edge,w$edge_length,states,w$rate$rate_index,1L)
  on.exit(fast_ace_fit_release_cpp(ptr))
  rates <- function(a,b,d) {
    Q <- matrix(c(0,a,b,0,0,d,0,0,0),3,byrow=TRUE)
    Q[row(Q)!=col(Q)]
  }
  # Changing a+b's split keeps all eigenvalues fixed while changing Q and its
  # eigenvectors. More than 24 subsequent generators exercise bounded eviction.
  base <- rates(.2,0,.6);other <- rates(.15,.05,.6)
  parameters <- c(list(base,other),lapply(1:60,function(i)rates(.05+i/1000,0,.6)),list(base,other))
  bytes <- fast_ace_fit_bytes_cpp(ptr)
  for(z in list(states,rev(states),replace(states,c(3,17,40),NA_integer_))) {
    fast_ace_fit_update_cpp(ptr,z);seed <- .Random.seed
    probabilities <- list()
    for(i in seq_along(parameters)) {
      par <- parameters[[i]]
      # Raw reuse_buffers calls have fresh, uncached scratch. Their scalar
      # evaluator and matrix products supply the exact reference on this BLAS.
      expected <- fast_ace_discrete_deviance_cpp(w$edge,w$edge_length,z,w$rate$rate_index,par,TRUE)
      expect_true(identical(fast_ace_fit_deviance_cpp(ptr,par),as.numeric(expected),
                            num.eq=FALSE,single.NA=FALSE))
      ref <- fast_ace_discrete_liks_cpp(w$edge,w$edge_length,z,w$rate$rate_index,par,FALSE,TRUE)
      got <- fast_ace_fit_liks_cpp(ptr,par,FALSE)
      expect_true(identical(got,ref[c("deviance","loglik","lik.anc")],
                            num.eq=FALSE,single.NA=FALSE))
      expect_identical(.Random.seed,seed)
      if(i<=2L)probabilities[[i]] <- got$lik.anc
    }
    expect_false(identical(probabilities[[1L]],probabilities[[2L]]))
    # A zero generator has no support for these observations and can stop
    # pruning early. A later valid trial must still see complete exponent rows.
    expect_identical(fast_ace_fit_deviance_cpp(ptr,rep(0,6)),1e50)
    expect_identical(fast_ace_fit_deviance_cpp(ptr,base),as.numeric(
      fast_ace_discrete_deviance_cpp(w$edge,w$edge_length,z,w$rate$rate_index,base,TRUE)))
  }
  expect_lt(fast_ace_fit_bytes_cpp(ptr)-bytes,1024^2)
})

test_that("large serial fits retain bounded exponent storage", {
  set.seed(20261021)
  tr <- ape::rtree(4096)
  w <- .fast_ace_workspace(.fast_ace_tree_workspace(tr),"ARD",3L)
  states <- rep(1:3,length.out=4096)
  ptr <- fast_ace_fit_prepare_cpp(w$edge,w$edge_length,states,w$rate$rate_index,1L)
  on.exit(fast_ace_fit_release_cpp(ptr))
  seed <- .Random.seed
  for(i in 1:30) {
    Q <- matrix(c(0,.05+i/1000,0,0,0,.6,0,0,0),3,byrow=TRUE)
    par <- Q[row(Q)!=col(Q)]
    expect_identical(fast_ace_fit_deviance_cpp(ptr,par),as.numeric(
      fast_ace_discrete_deviance_cpp(w$edge,w$edge_length,states,w$rate$rate_index,par,TRUE)))
  }
  expect_lt(fast_ace_fit_bytes_cpp(ptr),1024^2)
  expect_identical(.Random.seed,seed)
})

test_that("cache eligibility boundaries preserve complex ARD transitions", {
  set.seed(20261022)
  Q <- matrix(c(0,.4,.01,.02,0,.5,.6,.03,0),3,byrow=TRUE)
  par <- Q[row(Q)!=col(Q)]
  # 64 edges enables the cache; >8191/3 edges retains the batch/serial paths.
  for(n in c(32L,33L,1366L,1367L)) {
    tr <- ape::rtree(n)
    tr$edge.length <- exp(runif(nrow(tr$edge),-2,1))
    w <- .fast_ace_workspace(.fast_ace_tree_workspace(tr),"ARD",3L)
    states <- rep(1:3,length.out=n);states[7L] <- NA_integer_
    expected <- fast_ace_discrete_liks_cpp(w$edge,w$edge_length,states,
      w$rate$rate_index,par,FALSE,TRUE)[c("deviance","loglik","lik.anc")]
    for(threads in c(1L,2L))local({
      ptr <- fast_ace_fit_prepare_cpp(w$edge,w$edge_length,states,w$rate$rate_index,threads)
      on.exit(fast_ace_fit_release_cpp(ptr))
      seed <- .Random.seed
      for(i in 1:2) {
        expect_true(identical(fast_ace_fit_liks_cpp(ptr,par,FALSE),expected,
                              num.eq=FALSE,single.NA=FALSE))
        expect_identical(.Random.seed,seed)
      }
    })
  }
})
