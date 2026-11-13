test_that("lambda common branch units do not discard representable likelihoods", {
  set.seed(17)
  tree <- ape::rtree(20)
  x <- setNames(rnorm(20),tree$tip.label)
  make_core <- function(z) compile_tree_cpp(z$edge,z$edge.length,length(z$tip.label))
  core <- make_core(tree)
  grid <- c(0,.1,.25,.75,1)
  reference <- fast_lambda_tree_fixed_cpp(core,matrix(x),grid)
  public_ref <- fast_lambda(tree,x,test=TRUE,progress=FALSE)
  for (s in c(1e-309,1.9e-309,1e-300,2^-1020,2^1000)) {
    z <- tree; z$edge.length <- z$edge.length*s
    expect_true(all(is.finite(z$edge.length)&z$edge.length>0))
    got <- fast_lambda_tree_fixed_cpp(make_core(z),cbind(x,x*1e-200,x*1e200),grid)
    expect_true(all(got$valid))
    expect_equal(got$logLik[,1],reference$logLik[,1],tolerance=1e-10)
    expect_equal(got$logLik[,2],reference$logLik[,1]+20*200*log(10),tolerance=1e-10)
    expect_equal(got$logLik[,3],reference$logLik[,1]-20*200*log(10),tolerance=1e-10)
    expect_equal(got$gls_mean[,1],reference$gls_mean[,1],tolerance=1e-10)
    a <- suppressWarnings(fast_lambda(z,x,test=TRUE,progress=FALSE))
    expect_equal(a$lambda,public_ref$lambda,tolerance=1e-5)
    expect_equal(a$logL,public_ref$logL,tolerance=1e-10)
    expect_equal(a$P,public_ref$P,tolerance=1e-8)
    expect_true(is.finite(a$lambda)&is.finite(a$P))
    if(s==1e-309) expect_match(a$message,"sigma2 exceeds double range")
  }
})

test_that("joint branch and trait scaling preserves finite diffusion variance", {
  set.seed(27)
  tr <- ape::rtree(20); x <- rnorm(20)
  core <- function(z) compile_tree_cpp(z$edge,z$edge.length,20L)
  a <- fast_lambda_tree_fixed_cpp(core(tr),matrix(x),c(0,.4,1))
  z <- tr; z$edge.length <- z$edge.length*2^-1020
  b <- fast_lambda_tree_fixed_cpp(core(z),matrix(x*2^-510),c(0,.4,1))
  expect_true(all(b$valid))
  expect_equal(b$sigma2,a$sigma2,tolerance=1e-10)
  expect_equal(b$logLik,a$logLik+20*510*log(2),tolerance=1e-10)
})

test_that("pooled ML and REML and bootstrap are invariant to common branch units", {
  set.seed(917)
  trees <- list(ape::rtree(12),ape::rtree(14))
  cores <- function(ts) lapply(ts,function(z) compile_tree_cpp(z$edge,z$edge.length,length(z$tip.label)))
  values <- lapply(trees,function(z)rnorm(length(z$tip.label)))
  scaled <- lapply(trees,function(z){z$edge.length<-z$edge.length*2^-1020;z})
  for(reml in c(FALSE,TRUE)) {
    set.seed(1); a<-heterogeneity_lambda_calibrate_cpp(cores(trees),values,.1,3L,TRUE,reml,.4)
    set.seed(1); b<-heterogeneity_lambda_calibrate_cpp(cores(scaled),values,.1,3L,TRUE,reml,.4)
    expect_equal(b$null_logLik,a$null_logLik,tolerance=1e-10)
    expect_equal(b$lambda_common,a$lambda_common,tolerance=1e-10)
    expect_equal(b$P_heterogeneity,a$P_heterogeneity,tolerance=1e-10)
    expect_length(b$null,3L)
    expect_equal(b$null,a$null,tolerance=1e-5)
  }
})