test_that("updating prepared states preserves pruning, likelihoods and RNG", {
  set.seed(202610)
  tr <- ape::rtree(64)
  for(k in c(2L,4L,8L)) {
    w <- .fast_ace_workspace(.fast_ace_tree_workspace(tr),"ARD",k)
    states <- rep(seq_len(k),length.out=64)
    ptr <- fast_ace_fit_prepare_cpp(w$edge,w$edge_length,states,w$rate$rate_index,1L)
    for(z in list(states,rev(states),replace(states,c(3,19),NA_integer_))) {
      seed <- .Random.seed
      fast_ace_fit_update_cpp(ptr,z)
      expect_identical(.Random.seed,seed)
      for(par in list(rep(.1,w$rate$np),seq(.01,.4,length.out=w$rate$np))) {
        expect_identical(fast_ace_fit_deviance_cpp(ptr,par),
          fast_ace_discrete_deviance_cpp(w$edge,w$edge_length,z,w$rate$rate_index,par,TRUE))
        for(marginal in c(FALSE,TRUE)) expect_identical(
          fast_ace_fit_liks_cpp(ptr,par,marginal),fast_ace_discrete_liks_cpp(
            w$edge,w$edge_length,z,w$rate$rate_index,par,marginal,TRUE))
      }
    }
    before <- fast_ace_fit_deviance_cpp(ptr,rep(.1,w$rate$np))
    expect_error(fast_ace_fit_update_cpp(ptr,rep(k+1L,64)),"state indices")
    expect_error(fast_ace_fit_update_cpp(ptr,states[-1]),"state count")
    expect_identical(fast_ace_fit_deviance_cpp(ptr,rep(.1,w$rate$np)),before)
    expect_gt(fast_ace_fit_bytes_cpp(ptr),0)
    fast_ace_fit_release_cpp(ptr)
    fast_ace_fit_release_cpp(ptr)
    expect_identical(fast_ace_fit_bytes_cpp(ptr),0)
    expect_error(fast_ace_fit_deviance_cpp(ptr,rep(.1,w$rate$np)),"released")
  }
})

test_that("pooled fits retain fresh optimizer results and cleanup on eviction", {
  old <- options(fastphylosig.ace_cache_bytes=64*1024^2)
  on.exit(options(old))
  pool <- .fast_ace_native_pool();on.exit(.fast_ace_pool_close(pool),add=TRUE)
  set.seed(778);tr <- ape::rtree(32)
  tw <- .fast_ace_tree_workspace(tr)
  pooled <- .fast_ace_workspace(tw,"ARD",3L,pool)
  fresh <- .fast_ace_workspace(tw,"ARD",3L)
  run <- function(w,z) {
    set.seed(709)
    fit <- suppressWarnings(.fast_ace_fit_workspace(z,letters[1:3],w,.1,
      CI=TRUE,estimate_se=FALSE))
    list(fit=fit,seed=.Random.seed,next_normal=rnorm(4),next_uniform=runif(4))
  }
  states <- rep(1:3,length.out=32)
  for(z in list(states,rev(states),replace(states,7,NA_integer_)))
    expect_identical(run(pooled,z),run(fresh,z))
  expect_length(ls(pool$entries),1)
  entry <- pool$entries[[ls(pool$entries)[1]]]
  expect_identical(pool$bytes,entry$bytes)
  pool$budget <- 1
  .fast_ace_pool_trim(pool)
  expect_length(ls(pool$entries),0)
  expect_identical(pool$bytes,0)
  expect_identical(fast_ace_fit_bytes_cpp(entry$pointer),0)
  expect_identical(run(pooled,states),run(fresh,states))
  expect_length(ls(pool$entries),0)
  pool$budget <- 0
  expect_identical(run(pooled,states),run(fresh,states))
  expect_length(ls(pool$entries),0)
  pool$budget <- 64*1024^2;run(pooled,states)
  # A serialized pool from another process has unusable pointers. Never use it.
  foreign <- unserialize(serialize(pooled,NULL));foreign$native_pool$pid <- -1L
  expect_identical(run(foreign,states),run(fresh,states))
  .fast_ace_pool_close(pool);.fast_ace_pool_close(pool)
  expect_length(ls(pool$entries),0)
  expect_identical(pool$bytes,0)
})

test_that("public results omit analysis environments and native pointers", {
  tr <- ape::stree(64,"balanced");tr$edge.length <- rep(1,nrow(tr$edge))
  x <- setNames(rep(letters[1:3],length.out=64),tr$tip.label)
  set.seed(412)
  h <- suppressWarnings(fast_hetero(tr,x,method="Delta",min_tips=20,
    max_tips=32,test=TRUE,nsim=2,delta_control=list(mcmc_sim=30L,thin=2L,burn=4L),
    progress=FALSE))
  has_native <- function(z) {
    if(is.environment(z)||typeof(z)=="externalptr")return(TRUE)
    if(is.list(z)&&any(vapply(z,has_native,logical(1))))return(TRUE)
    a <- attributes(z)
    !is.null(a)&&any(vapply(a,has_native,logical(1)))
  }
  expect_false(has_native(h))
  expect_identical(unserialize(serialize(h,NULL)),h)
  expect_null(attr(h$delta_control,"ace_pool"))
})

