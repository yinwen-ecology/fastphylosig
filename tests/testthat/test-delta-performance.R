test_that("native V2 traversal preserves the complete frozen graph contract", {
  set.seed(419)
  for(type in c("random","balanced","left","star")) for(n in c(2L,8L,128L)) {
    tr <- if(type=="random") ape::rtree(n) else ape::stree(n,if(n==2L) "balanced" else type)
    tr$edge.length <- seq_len(nrow(tr$edge))*.01
    tr$edge.length[1] <- 0
    tr$tip.label <- paste0(rep(c("中","é","a","z"),length.out=n),seq_len(n))
    map <- c(seq_len(n),1000000L+sample(seq_len(tr$Nnode))*9L)
    tr$edge[] <- map[tr$edge]
    j <- sample(nrow(tr$edge)); tr$edge <- tr$edge[j,,drop=FALSE]; tr$edge.length <- tr$edge.length[j]
    expect_identical(.v2_build_graph(tr),.delta_opt_ref_graph(tr))
  }
  # Unary nodes and branch distances whose sum overflows are valid representations.
  tr <- structure(list(edge=matrix(c(9,8,8,1,8,2),3,2,byrow=TRUE),
    edge.length=c(1e308,1e308,0),tip.label=c("b","a"),Nnode=2L),class="phylo")
  expect_identical(.v2_build_graph(tr),.delta_opt_ref_graph(tr))
  keys <- lapply(c("z","aa","a","a","中","é"),charToRaw)
  expect_identical(.v2_order_bytes(keys),.delta_opt_ref_order(keys))
  expect_identical(.v2_order_bytes(list()),integer())
  set.seed(5); seed <- .Random.seed
  .v2_build_graph(tr)
  expect_identical(.Random.seed,seed)
})

test_that("native graph validation preserves structural failure messages", {
  tr <- ape::stree(8,"balanced"); tr$edge.length <- rep(.1,nrow(tr$edge))
  bad <- list(tr,tr,tr,tr,tr)
  bad[[1]]$edge[1,2] <- bad[[1]]$edge[1,1]
  bad[[2]]$edge[2,2] <- bad[[2]]$edge[1,2]
  bad[[3]]$edge[1,1] <- 1L
  bad[[4]]$edge[1,1] <- 100L
  bad[[5]]$edge[1,2] <- 100L
  message <- function(fun,tree) tryCatch({fun(tree);NULL},error=function(e)conditionMessage(e))
  for(tree in bad) expect_identical(message(.v2_build_graph,tree),message(.delta_opt_ref_graph,tree))
})

test_that("prepared ACE fits reuse the existing pruning kernel exactly", {
  set.seed(183)
  for(k in c(2L,4L,7L)) for(model in c("ER","ARD")) for(repeated in c(FALSE,TRUE)) {
    tr <- ape::rtree(80)
    if(repeated) tr$edge.length <- rep(c(.1,.1,.2),length.out=nrow(tr$edge))
    state <- rep(seq_len(k),length.out=80); state[c(1,9,31)] <- NA_integer_
    w <- .fast_ace_workspace(.fast_ace_tree_workspace(tr),model,k)
    ptr <- fast_ace_fit_prepare_cpp(w$edge,w$edge_length,state,w$rate$rate_index)
    params <- list(rep(.1,w$rate$np),rep(0,w$rate$np),runif(w$rate$np),
      rep(-1,w$rate$np),rep(Inf,w$rate$np),rep(.1,w$rate$np))
    for(p in params) {
      expected <- fast_ace_discrete_deviance_cpp(w$edge,w$edge_length,state,w$rate$rate_index,p,FALSE)
      expect_identical(fast_ace_discrete_deviance_cpp(w$edge,w$edge_length,state,w$rate$rate_index,p,TRUE),expected)
      expect_identical(fast_ace_fit_deviance_cpp(ptr,p),expected)
    }
    # nlminb must follow exactly the same path, including convergence/message.
    a <- suppressWarnings(.fast_ace_fit_workspace(state,as.character(seq_len(k)),w,.1,estimate_se=FALSE))
    b <- suppressWarnings(.fast_ace_fit_workspace(state,as.character(seq_len(k)),w,.1,estimate_se=FALSE,reuse_buffers=FALSE))
    expect_identical(a,b)
    seed <- .Random.seed
    fast_ace_fit_release_cpp(ptr)
    fast_ace_fit_release_cpp(ptr)
    expect_identical(.Random.seed,seed)
    expect_error(fast_ace_fit_deviance_cpp(ptr,params[[1]]),"has been released")
  }
})

test_that("corrected MCMC agrees with an independent raw-parameter posterior", {
  # The pre-0.3.2 frozen fixture pins a different target and estimator. Do not
  # regenerate it to bless a scientific change: compare the declared posterior.
  for (k in c(2L, 4L)) for (entropy in 1:3) {
    set.seed(19 + k)
    p <- matrix(runif(40 * k, .05, 1), 40, k)
    p <- p / rowSums(p) # no endpoint jitter, so initialization is the first RNG
    set.seed(714)
    fit <- delta_mcmc_cpp(p, .1, .5, 1000L, 5L, 101L, entropy, TRUE)
    seed <- .Random.seed
    set.seed(714)
    ref <- .delta_posterior_reference(fit$entropy, .1, .5, 1000L, 5L, 101L)
    expect_equal(fit$alpha_chain, ref$alpha, tolerance = 1e-9)
    expect_equal(fit$beta_chain, ref$beta, tolerance = 1e-9)
    expect_equal(fit$delta, mean(ref$beta / ref$alpha), tolerance = 1e-10)
    expect_identical(.Random.seed, seed)
    set.seed(714)
    compact <- delta_mcmc_cpp(p, .1, .5, 1000L, 5L, 101L, entropy, FALSE)
    expect_identical(compact$delta, fit$delta)
    expect_identical(.Random.seed, seed)
  }
})
