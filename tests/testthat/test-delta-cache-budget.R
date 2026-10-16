test_that("prepared ACE caches beyond the historical group cap exactly", {
  set.seed(912)
  tr <- ape::rtree(6000L)
  tr$edge.length <- rep(seq(.01,1,length.out=5000L),length.out=nrow(tr$edge))
  for(k in c(2L,4L,8L)) {
    state <- rep(seq_len(k),length.out=6000L);state[seq(3,6000L,19)]<-NA_integer_
    w <- .fast_ace_workspace(.fast_ace_tree_workspace(tr),'ARD',k)
    rates <- seq(.01,.3,length.out=w$rate$np)
    expected <- fast_ace_discrete_deviance_cpp(w$edge,w$edge_length,state,
      w$rate$rate_index,rates,FALSE)
    for(threads in c(1L,2L)) {
      ptr <- fast_ace_fit_prepare_cpp(w$edge,w$edge_length,state,w$rate$rate_index,threads)
      actual <- fast_ace_fit_deviance_cpp(ptr,rates)
      fast_ace_fit_release_cpp(ptr)
      expect_identical(actual,expected)
    }
  }
})
