# Independent matrix pruning; includes the uniform-root normalization.
.er_reference_oracle <- function(w, states, q) {
  k <- w$nl; n <- length(states); lik <- matrix(0,2*n-1,k); total <- 0
  for(i in seq_len(n)) if(is.na(states[i]))lik[i,] <- 1 else lik[i,states[i]] <- 1
  for(e in seq(1,nrow(w$edge),by=2)) {
    v <- lapply(e:(e+1L),function(r) {
      retain <- exp(-k*q*w$edge_length[r])
      P <- diag(retain,k)+matrix((1-retain)/k,k,k)
      as.numeric(P %*% lik[w$edge[r,2L],])
    })
    value <- v[[1]]*v[[2]]; scale <- sum(value)
    if(scale<=0)return(Inf)
    lik[w$edge[e,1L],] <- value/scale; total <- total+log(scale)
  }
  -2*(total-log(k))
}

test_that("ER reference pruning agrees with independent analytic matrix pruning", {
  set.seed(851); tr <- ape::rtree(32)
  for(k in c(2L,4L,8L)) for(scale in c(1e-6,1,1e6)) {
    tree <- tr; tree$edge.length <- tree$edge.length*scale
    w <- .fast_ace_workspace(.fast_ace_tree_workspace(tree),"ER",k)
    x <- rep(seq_len(k),length.out=32); x[c(4,21)] <- NA_integer_
    ptr <- fast_ace_fit_prepare_cpp(w$edge,w$edge_length,x,w$rate$rate_index,1L)
    seed <- .Random.seed
    for(q in c(0,.005,.1,1,10,Inf)/scale) {
      expect_equal(fast_ace_er_deviance_cpp(ptr,q)+2*log(k),
        .er_reference_oracle(w,x,q),tolerance=1e-9)
      if(is.finite(q)&&q>0) expect_equal(fast_ace_er_deviance_cpp(ptr,q),
        fast_ace_fit_deviance_cpp(ptr,q),tolerance=1e-8)
    }
    fit <- fast_ace_er_reference_cpp(ptr)
    expect_identical(.Random.seed,seed)
    expect_equal(fit$diagnostics$uniform_root_loglik,
      -.5*.er_reference_oracle(w,x,fit$rates),tolerance=1e-9)
    expect_true(fit$diagnostics$independence_deviance_gap>=0)
    # Independently scan more densely and require no materially better point.
    grid <- exp(seq(log(.001/scale),log(30/scale),length.out=101))
    oracle <- vapply(grid,function(q).er_reference_oracle(w,x,q),numeric(1))
    expect_lte(-2*fit$diagnostics$uniform_root_loglik,min(oracle)+1e-7)
    fast_ace_fit_release_cpp(ptr)
    expect_error(fast_ace_er_reference_cpp(ptr),"released")
  }
})

test_that("ER endpoints and invalid models are explicit", {
  tree <- ape::stree(4,"balanced"); tree$edge.length <- rep(1,6)
  w <- .fast_ace_workspace(.fast_ace_tree_workspace(tree),"ER",2L)
  for(x in list(rep(1L,4),c(1L,2L,1L,2L))) {
    ptr <- fast_ace_fit_prepare_cpp(w$edge,w$edge_length,x,w$rate$rate_index,1L)
    fit <- fast_ace_er_reference_cpp(ptr)
    expect_identical(fit$diagnostics$endpoint,if(all(x==1L))"zero" else "independent")
    expect_identical(fit$rates,if(all(x==1L))0 else Inf)
    fast_ace_fit_release_cpp(ptr)
  }
  ard <- .fast_ace_workspace(.fast_ace_tree_workspace(tree),"ARD",2L)
  ptr <- fast_ace_fit_prepare_cpp(ard$edge,ard$edge_length,rep(1:2,2),ard$rate$rate_index,1L)
  on.exit(fast_ace_fit_release_cpp(ptr))
  expect_error(fast_ace_er_reference_cpp(ptr),"common off-diagonal")
  expect_error(.heterogeneity_delta_control(list(reference_fit="profile")),"requires model='ER'")
  expect_error(.heterogeneity_delta_control(list(model="ER",reference_fit="other")),"arg")
})

test_that("independent ER endpoint simulation has analytic states and identical RNG", {
  set.seed(729); tree <- ape::rtree(16)
  tree$edge.length[1L] <- 0
  core <- compile_tree_cpp(tree$edge,tree$edge.length,16L)
  for(k in c(2L,4L,8L)) {
    index <- matrix(1L,k,k); diag(index) <- 0L
    ptr <- heterogeneity_mk_prepare_cpp(core,index,Inf)
    set.seed(957); states <- heterogeneity_mk_simulate_cpp(ptr,9L); seed <- .Random.seed
    set.seed(957); expected <- .hetero_er_states(core,k,Inf,9L)
    expect_identical(states,expected); expect_identical(.Random.seed,seed)
    heterogeneity_mk_release_cpp(ptr)
  }
  index <- .fast_ace_rate_index("ARD",3L)$rate_index
  expect_error(heterogeneity_mk_prepare_cpp(core,index,rep(Inf,6)),"invalid")
})

test_that("optional common profile fit leaves observed Delta and default path unchanged", {
  tree <- ape::stree(64,"balanced"); tree$edge.length <- rep(1,126)
  x <- setNames(rep(1:4,16),tree$tip.label)
  control <- list(model="ER",mcmc_sim=600L,thin=5L,burn=100L)
  run <- function(ctrl,test=TRUE) {
    set.seed(953)
    value <- suppressWarnings(fast_hetero(tree,x,method="Delta",test=test,overall=TRUE,
      nsim=2L,max_tips=16L,min_tips=12L,keep_null=TRUE,delta_control=ctrl,progress=FALSE))
    list(value=value,seed=.Random.seed)
  }
  default <- run(control)
  explicit <- run(c(control,list(reference_fit="nlminb")))
  # Optional control metadata is handled separately; numerical outputs match.
  expect_identical(default$value$clades,explicit$value$clades)
  expect_identical(default$value$overall,explicit$value$overall)
  expect_identical(default$value$heterogeneity,explicit$value$heterogeneity)
  expect_identical(default$value$null,explicit$value$null)
  expect_identical(default$seed,explicit$seed)
  profile <- run(c(control,list(reference_fit="profile")))
  expect_identical(profile$value$clades,default$value$clades)
  expect_identical(profile$value$overall,default$value$overall)
  expect_true(profile$value$heterogeneity$reference_model[[1]]$fit_diagnostics$endpoint %in%
    c("zero","interior","independent"))
  scan <- run(control,FALSE); profile_scan <- run(c(control,list(reference_fit="profile")),FALSE)
  expect_identical(scan$value$clades,profile_scan$value$clades)
  expect_identical(scan$value$overall,profile_scan$value$overall)
  expect_identical(scan$seed,profile_scan$seed)
})
