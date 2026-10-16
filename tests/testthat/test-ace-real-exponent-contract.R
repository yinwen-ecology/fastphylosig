test_that("real-exponent specialization retains the original ACE path", {
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
        expect_true(identical(as.numeric(new),as.numeric(old),num.eq=FALSE))
        expect_identical(.Random.seed,seed)
        if(is.finite(new)&&new<1e49)for(marginal in c(FALSE,TRUE)) {
          old_liks <- fast_ace_discrete_liks_cpp(w$edge,w$edge_length,y,
            w$rate$rate_index,par,marginal,FALSE)
          new_liks <- fast_ace_fit_liks_cpp(ptr,par,marginal)
          expect_true(identical(new_liks,old_liks,num.eq=FALSE,single.NA=FALSE))
          expect_identical(.Random.seed,seed)
        }
      }
    }
    fast_ace_fit_release_cpp(ptr)
  }
})

