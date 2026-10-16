.nuisance_fixture <- function() {
  tree <- ape::stree(4,type="balanced")
  tree$edge.length <- c(.29,.51,.83,1.17,1.41,.67)
  core <- compile_tree_cpp(tree$edge,tree$edge.length,4L)
  x <- t(as.matrix(expand.grid(rep(list(1:3),4L))))
  storage.mode(x) <- "integer"
  list(tree=tree,core=core,x=x)
}
.nuisance_oracle <- function(z,Q) {
  # Independent base-R pruning and expm, not the native interval engine.
  transitions <- lapply(z$core$branch_length_by_node,function(t)
    if(t==0)diag(nrow(Q)) else expm::expm(Q*t))
  vapply(seq_len(ncol(z$x)),function(b) {
    a <- matrix(1,z$core$n_total,nrow(Q))
    for(i in seq_len(z$core$n_tip))if(!is.na(z$x[i,b])) {
      a[i,] <- 0;a[i,z$x[i,b]] <- 1
    }
    for(node in z$core$postorder) {
      p <- z$core$parent[node]
      if(p) a[p,] <- a[p,]*as.numeric(transitions[[node]] %*% a[node,])
    }
    log(mean(a[z$core$root,]))
  },numeric(1))
}
.nuisance_pointer <- function(z,model="ER",tail=rep(TRUE,ncol(z$x))) {
  g <- matrix(-4*log(3),ncol(z$x),2)
  index <- .fast_ace_rate_index(model,3L)$rate_index
  heterogeneity_nuisance_prepare_cpp(z$core,z$x,index,g,tail,model=="ER")
}

test_that("ER and ARD interval likelihoods contain independent pruning", {
  skip_if_not_installed("expm")
  set.seed(880)
  z <- .nuisance_fixture()
  seed <- .Random.seed
  er <- .nuisance_pointer(z); ard <- .nuisance_pointer(z,"ARD")
  on.exit({heterogeneity_nuisance_release_cpp(er);heterogeneity_nuisance_release_cpp(ard)})
  for(q in c(0,.0001,.2,3,100)) {
    fit <- heterogeneity_nuisance_box_cpp(er,q,q,.001,TRUE)
    Q <- matrix(q,3,3);diag(Q) <- -2*q
    expected <- .nuisance_oracle(z,Q)
    expect_true(all(fit$log_likelihood[,1]<=expected+1e-12))
    expect_true(all(fit$log_likelihood[,2]>=expected-1e-12))
    if(q>0)expect_lt(max(fit$log_likelihood[,2]-fit$log_likelihood[,1]),1e-8)
  }
  bounded <- heterogeneity_nuisance_box_cpp(er,.1,.8,.001,TRUE)$log_likelihood
  for(q in seq(.1,.8,length.out=7)) {
    Q <- matrix(q,3,3);diag(Q) <- -2*q
    expected <- .nuisance_oracle(z,Q)
    expect_true(all(bounded[,1]<=expected+1e-12))
    expect_true(all(bounded[,2]>=expected-1e-12))
  }
  index <- .fast_ace_rate_index("ARD",3L)$rate_index
  rates <- c(.03,.08,.12,.27,.4,.7)
  for(intensity in c(.01,1,100)) {
    fit <- heterogeneity_nuisance_box_cpp(ard,rates*intensity,rates*intensity,.001,TRUE)
    Q <- matrix(0,3,3);off <- row(Q)!=col(Q);Q[off] <- (rates*intensity)[index[off]]
    diag(Q) <- -rowSums(Q);expected <- .nuisance_oracle(z,Q)
    expect_true(all(fit$log_likelihood[,1]<=expected+1e-11))
    expect_true(all(fit$log_likelihood[,2]>=expected-1e-11))
    expect_lt(max(fit$log_likelihood[,2]-fit$log_likelihood[,1]),1e-7)
  }
  fit <- heterogeneity_nuisance_box_cpp(ard,rates*.7,rates*1.3,.001,TRUE)
  for(multiplier in c(.7,.9,1.1,1.3)) {
    Q <- matrix(0,3,3);off <- row(Q)!=col(Q);Q[off] <- (rates*multiplier)[index[off]]
    diag(Q) <- -rowSums(Q);expected <- .nuisance_oracle(z,Q)
    expect_true(all(fit$log_likelihood[,1]<=expected+1e-11))
    expect_true(all(fit$log_likelihood[,2]>=expected-1e-11))
  }
  expect_identical(.Random.seed,seed)
})

test_that("proposal density uses the sampler and missing-tip marginalization", {
  skip_if_not_installed("expm")
  set.seed(880)
  z <- .nuisance_fixture();z$x[2,seq(1,81,3)] <- NA_integer_
  index <- matrix(1L,3,3);diag(index) <- 0L
  proposals <- lapply(c(.1,1,Inf),function(q)heterogeneity_mk_prepare_cpp(z$core,index,q))
  on.exit(for(p in proposals)heterogeneity_mk_release_cpp(p))
  seed <- .Random.seed
  likelihood <- lapply(proposals,function(p)heterogeneity_mk_loglik_bounds_cpp(p,z$core,z$x))
  for(i in 1:2) {
    Q <- matrix(c(.1,1)[i],3,3);diag(Q) <- -rowSums(Q)+diag(Q)
    expected <- .nuisance_oracle(z,Q)
    expect_true(all(likelihood[[i]][,1]<=expected+1e-12))
    expect_true(all(likelihood[[i]][,2]>=expected-1e-12))
  }
  independent <- -colSums(!is.na(z$x))*log(3)
  expect_true(all(likelihood[[3]][,1]<=independent+1e-12))
  expect_true(all(likelihood[[3]][,2]>=independent-1e-12))
  mixture <- heterogeneity_nuisance_mixture_cpp(likelihood)
  expected <- log(rowMeans(sapply(likelihood,function(v)exp(rowMeans(v)))))
  expect_true(all(mixture[,1]<=expected+1e-12))
  expect_true(all(mixture[,2]>=expected-1e-12))
  expect_identical(.Random.seed,seed)
})

test_that("corrected importance ranks and continuous covers are conservative", {
  skip_if_not_installed("expm")
  z <- .nuisance_fixture();z$x <- z$x[,c(5,9,17,37,55,63,70),drop=FALSE]
  tail <- c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE)
  pointer <- .nuisance_pointer(z,tail=tail)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  q <- .21;Q <- matrix(q,3,3);diag(Q) <- -2*q
  w <- exp(.nuisance_oracle(z,Q)+4*log(3))
  expected <- sum(w[tail])/sum(w)
  point <- heterogeneity_nuisance_box_cpp(pointer,q,q,0)
  expect_lte(point$lower,expected+1e-12);expect_gte(point$upper,expected-1e-12)
  settings <- .heterogeneity_nuisance_control(list(max_boxes=12L,tolerance=0,gamma=0))
  domain <- list(bounds=matrix(c(.1,.7),1,2),model="ER")
  result <- .heterogeneity_nuisance_search(pointer,domain,settings,1/3)
  expect_false(result$converged)
  expect_equal(result$cells,12L)
  for(q in seq(.1,.7,length.out=21))
    expect_gte(result$upper,heterogeneity_nuisance_box_cpp(pointer,q,q,0)$upper-1e-10)
  expect_gte(result$upper,result$lower)
  expect_equal(heterogeneity_nuisance_add_gamma_cpp(1,1,.001),c(1,1))
})

test_that("likelihood compression preserves individual MCMC upper-tail flags", {
  skip_if_not_installed("expm")
  z <- .nuisance_fixture();z$x <- z$x[,c(5,5,9,9,5,17),drop=FALSE]
  tail <- c(TRUE,FALSE,FALSE,TRUE,FALSE,TRUE)
  pointer <- .nuisance_pointer(z,tail=tail)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  q <- .21;Q <- matrix(q,3,3);diag(Q) <- -2*q
  w <- exp(.nuisance_oracle(z,Q)+4*log(3))
  p <- heterogeneity_nuisance_box_cpp(pointer,q,q,0)
  expect_equal(p$n_unique_patterns,3L)
  expect_lte(p$lower,sum(w[tail])/sum(w)+1e-12)
  expect_gte(p$upper,sum(w[tail])/sum(w)-1e-12)
})

test_that("nuisance public controls declare the model without changing plugin defaults", {
  z <- .hetero_fixture();x <- setNames(rep(c("A","B","C"),length.out=64),z$tree$tip.label)
  args <- list(tree=z$tree,data=x,method="Delta",nodes=z$nodes,min_tips=20,
    test=TRUE,nsim=3,progress=FALSE,keep_null=TRUE)
  base <- list(model="ER",mcmc_sim=1200L,thin=5L,burn=100L)
  expect_error(do.call(fast_hetero,c(args,list(delta_control=c(base,list(reference_test="nuisance"))))),"state_space")
  expect_error(do.call(fast_hetero,c(args,list(delta_control=c(base,list(reference_test="nuisance",
    state_space=c("A","B")))))),"Observed states")
  set.seed(881)
  h <- suppressWarnings(do.call(fast_hetero,c(args,list(delta_control=c(base,list(reference_test="nuisance",
    state_space=c("A","B","C"),nuisance_control=list(max_boxes=3L)))))))
  expect_identical(h$heterogeneity$P_heterogeneity_kind,"nuisance_upper_bound")
  expect_gte(h$heterogeneity$P_heterogeneity,h$heterogeneity$P_nuisance_lower)
  expect_true(is.na(h$heterogeneity$MCSE_P_heterogeneity))
  expect_true(is.na(h$heterogeneity$P_heterogeneity_lower))
  expect_equal(length(h$null[[1]]),3L)
  expect_equal(nrow(h$heterogeneity$reference_model[[1]]$reference_patterns),64L)
  expect_equal(h$heterogeneity$nuisance_reference_attempted,3L)
  seed <- .Random.seed
  device <- tempfile(fileext=".pdf");grDevices::pdf(device)
  picture <- plot(h,tree=z$tree)
  grDevices::dev.off();unlink(device)
  expect_true(any(grepl("nuisance upper bound",picture$summary_text)))
  expect_identical(.Random.seed,seed)
  # A missing biological state remains part of the reference model.
  args$data <- setNames(rep(c("A","C"),32),z$tree$tip.label)
  set.seed(881)
  h <- suppressWarnings(do.call(fast_hetero,c(args,list(delta_control=c(base,list(reference_test="nuisance",
    state_space=c("A","B","C"),nuisance_control=list(max_boxes=2L)))))))
  expect_identical(h$heterogeneity$reference_model[[1]]$states,c("A","B","C"))
  # Full, unbounded ARD searches certified infinite tails without a rate cap.
  args$data <- x;base$model <- "ARD"
  h <- suppressWarnings(do.call(fast_hetero,c(args,list(delta_control=c(base,list(reference_test="nuisance",
    state_space=c("A","B","C"),nuisance_control=list(max_boxes=2L)))))))
  expect_true(h$heterogeneity$status %in% c("nuisance_bound","nuisance_budget_limited"))
  expect_lte(h$heterogeneity$P_heterogeneity,1)
  expect_equal(h$heterogeneity$nuisance_reference_attempted,3)
  expect_true(h$heterogeneity$reference_model[[1]]$search$boundary_covered)
  expect_true(any(is.infinite(h$heterogeneity$reference_model[[1]]$bounds[,2])))
  # A predefined finite domain runs the same bank and retains each diagnosis.
  h <- suppressWarnings(do.call(fast_hetero,c(args,list(delta_control=c(base,list(reference_test="nuisance",
    state_space=c("A","B","C"),nuisance_control=list(rate_bounds=matrix(c(rep(.03,6),rep(.7,6)),6,2),
    max_boxes=2L)))))))
  model <- h$heterogeneity$reference_model[[1]]
  expect_true(model$rate_bounds_supplied)
  expect_equal(dim(model$reference_diagnostics),c(3L,10L))
  expect_equal(ncol(model$reference_patterns),4L)
  expect_gte(h$heterogeneity$P_nuisance_upper,h$heterogeneity$P_nuisance_lower)
})

test_that("unbounded ARD forest and Euler envelopes include independent transitions", {
  skip_if_not_installed("expm")
  set.seed(161111);seed <- .Random.seed
  L <- matrix(.03,3,3);U <- matrix(.7,3,3);diag(L)<-diag(U)<-0
  positions <- which(row(L)!=col(L))
  for(mask in list(positions[1],positions[c(1,3)],positions[c(1,4,5)],positions)) {
    upper<-U;upper[mask]<-Inf
    fit<-heterogeneity_ard_boundary_cpp(L,upper,.7)
    expect_true(all(fit$lower<=fit$upper))
    expect_equal(fit$Euler_error,32/2^24)
    for(q in c(.04,.5,10,1e4)) {
      Q<-(L+U)/2;Q[mask]<-q;diag(Q)<- -rowSums(Q)
      P<-expm::expm(Q*.7);R<-solve(diag(3)-Q*.7/2^24)
      expect_true(all(fit$lower<=P+1e-8))
      expect_true(all(fit$upper>=P-1e-8))
      expect_true(all(fit$resolvent_lower<=R+1e-10))
      expect_true(all(fit$resolvent_upper>=R-1e-10))
    }
  }
  # A single arbitrarily fast exit has a useful bound, not [0,1].
  L[1,2]<-1e5;U[1,2]<-Inf
  fit<-heterogeneity_ard_boundary_cpp(L,U,1)
  expect_lt(max(fit$upper[,1]),1e-4)
  expect_identical(.Random.seed,seed)
})

test_that("stochastic-row pruning tightens likelihood bounds and can exclude infinite tails", {
  z<-.nuisance_fixture();pointer<-.nuisance_pointer(z,"ARD")
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  a<-heterogeneity_nuisance_box_cpp(pointer,rep(.1,6),rep(.7,6),0,TRUE,constrain_rows=FALSE)
  b<-heterogeneity_nuisance_box_cpp(pointer,rep(.1,6),rep(.7,6),0,TRUE)
  expect_true(all(b$log_likelihood[,1]>=a$log_likelihood[,1]-1e-10))
  expect_true(all(b$log_likelihood[,2]<=a$log_likelihood[,2]+1e-10))
  expect_true(any(b$log_likelihood[,2]<a$log_likelihood[,2]-1e-6))
  lo<-rep(.1,6);hi<-rep(.7,6);lo[3]<-1e5;hi[3]<-Inf
  tail<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,.001)
  expect_true(tail$excluded)
  expect_equal(tail$upper,0)
  full<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,.001,TRUE)
  expect_identical(tail[names(tail)],full[names(tail)])
})

test_that("ARD mixture components are fixed before seeing traits", {
  z<-.heterogeneity_delta_control(list(model="ARD",reference_test="nuisance",state_space=letters[1:3]))
  d<-.heterogeneity_nuisance_domain(z,letters[1:3]);seed<-.Random.seed
  p<-.heterogeneity_ard_proposals(d,z$nuisance_control,1)
  expect_length(p,6L)
  expect_true(all(vapply(p,function(r)length(r)==6L&&all(is.finite(r)&r>0),logical(1))))
  z$nuisance_control$proposal_family<-"ER"
  expect_length(.heterogeneity_ard_proposals(d,z$nuisance_control,1),0L)
  expect_identical(.Random.seed,seed)
  expect_error(.heterogeneity_delta_control(list(model="ER",reference_test="nuisance",
    state_space=letters[1:3],nuisance_control=list(proposal_family="ARD"))),"requires")
})
