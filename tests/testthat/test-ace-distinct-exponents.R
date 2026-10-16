test_that("distinct-length exponent batches preserve ACE values and RNG", {
  set.seed(621)
  # Include the threshold boundary and a tree spanning several 2 MiB batches.
  for(case in list(c(512L,8L),c(513L,8L),c(2048L,4L),c(8192L,20L))) {
    n<-case[1L];k<-case[2L];tr<-ape::rtree(n)
    w<-.fast_ace_workspace(.fast_ace_tree_workspace(tr),"ARD",k)
    state<-sample.int(k,n,TRUE);state[seq(3,n,17)]<-NA_integer_
    par<-runif(w$rate$np,.02,.3)
    for(threads in c(1L,2L,4L)) {
      ptr<-fast_ace_fit_prepare_cpp(w$edge,w$edge_length,state,w$rate$rate_index,threads)
      for(p in list(par,par*10,rep(0,length(par)),rep(-1,length(par)),par)) {
        set.seed(139);seed<-.Random.seed
        expected<-fast_ace_discrete_deviance_cpp(w$edge,w$edge_length,state,w$rate$rate_index,p,TRUE)
        expect_identical(fast_ace_fit_deviance_cpp(ptr,p),expected)
        expect_identical(.Random.seed,seed)
      }
      fast_ace_fit_release_cpp(ptr)
    }
  }
  expect_error(fast_ace_fit_prepare_cpp(w$edge,w$edge_length,state,w$rate$rate_index,-1L),"non-negative")
})

test_that("thread choices retain the complete ARD optimizer trajectory", {
  set.seed(714);tr<-ape::rtree(2048);state<-sample.int(4L,2048,TRUE)
  w<-.fast_ace_workspace(.fast_ace_tree_workspace(tr),"ARD",4L)
  run<-function(threads) {
    set.seed(18)
    ptr<-fast_ace_fit_prepare_cpp(w$edge,w$edge_length,state,w$rate$rate_index,threads)
    on.exit(fast_ace_fit_release_cpp(ptr))
    calls<-list()
    dev<-function(p){v<-fast_ace_fit_deviance_cpp(ptr,p);calls[[length(calls)+1L]]<<-list(par=p[],value=v);v}
    fit<-nlminb(rep(.1,12),dev,lower=rep(0,12),upper=rep(1e50,12),control=list(iter.max=15L))
    list(fit=fit,calls=calls,seed=.Random.seed,next_normal=rnorm(5),next_uniform=runif(5))
  }
  expected<-run(1L)
  expect_identical(run(2L),expected)
  expect_identical(run(4L),expected)
})

test_that("Delta honors an explicit serial ACE option without changing results", {
  old<-getOption("fastphylosig.ace_threads")
  on.exit(options(fastphylosig.ace_threads=old))
  set.seed(13);tr<-ape::rtree(1100);y<-setNames(factor(sample(letters[1:4],1100,TRUE)),tr$tip.label)
  w<-.fast_ace_workspace(.fast_ace_tree_workspace(tr),"ARD",4L)
  run<-function(threads) {
    options(fastphylosig.ace_threads=threads);set.seed(78)
    result<-suppressWarnings(.delta_one(tr,y,.1,.5,300L,5L,20L,3L,"ARD","fast",TRUE,w))
    list(result=result,seed=.Random.seed,next_normal=rnorm(5))
  }
  expect_identical(run(0L),run(1L))
})
