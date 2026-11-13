test_that("lambda retains finite likelihood when trait variance is unrepresentable", {
  set.seed(17)
  tree <- ape::rtree(20)
  x <- setNames(rnorm(20), tree$tip.label)
  core <- compile_tree_cpp(tree$edge, tree$edge.length, 20L)
  base <- fast_lambda_tree_optimize_cpp(core, matrix(x), max_lambda=1, profile=TRUE, tol=1e-8)
  for (s in c(2^660, 2^-660, 1e200, 1e-200)) {
    got <- fast_lambda_tree_optimize_cpp(core, matrix(x*s), max_lambda=1, profile=TRUE, tol=1e-8)
    expect_true(got$valid[1])
    expect_equal(got$lambda, base$lambda, tolerance=2e-5)
    expect_equal(got$logLik, base$logLik-20*log(s), tolerance=1e-10)
    expect_equal(got$LR, base$LR, tolerance=1e-8)
    expect_equal(got$P, base$P, tolerance=1e-8)
    expect_equal(got$gls_mean/s, base$gls_mean, tolerance=2e-5)
    expect_match(got$status, "sigma2 .*double range")
    expect_true(if(s>1) is.infinite(got$sigma2) else got$sigma2==0)
    fixed <- fast_lambda_tree_fixed_cpp(core, cbind(x,x*s), c(0,.25,.75))
    expect_true(all(fixed$valid))
    expect_equal(fixed$logLik[,2], fixed$logLik[,1]-20*log(s), tolerance=1e-11)
    public <- fast_lambda(tree,x*s,test=TRUE,progress=FALSE)
    reference <- fast_lambda(tree,x,test=TRUE,progress=FALSE)
    expect_equal(public$lambda,reference$lambda,tolerance=2e-5)
    expect_equal(public$P,reference$P,tolerance=1e-7)
    expect_match(public$message,"sigma2 .*double range")
  }
  shifted <- x+1e8
  a <- fast_lambda_tree_fixed_cpp(core,matrix(shifted),c(0,.25,.75))
  b <- fast_lambda_tree_fixed_cpp(core,matrix(shifted*2^660),c(0,.25,.75))
  expect_equal(b$logLik,a$logLik-20*660*log(2),tolerance=1e-10)
  constant <- fast_lambda_tree_fixed_cpp(core,matrix(rep(1e200,20)),c(0,.5))
  expect_false(any(constant$valid))
})

test_that("pooled lambda ML and REML keep trait scale offsets", {
  set.seed(917)
  trees <- list(ape::rtree(12),ape::rtree(14))
  cores <- lapply(trees,function(z)compile_tree_cpp(z$edge,z$edge.length,length(z$tip.label)))
  values <- lapply(trees,function(z)rnorm(length(z$tip.label)))
  for (reml in c(FALSE,TRUE)) {
    set.seed(1)
    a <- heterogeneity_lambda_calibrate_cpp(cores,values,.1,3L,TRUE,reml,.4)
    set.seed(1)
    b <- heterogeneity_lambda_calibrate_cpp(cores,lapply(values,`*`,2^660),.1,3L,TRUE,reml,.4)
    expect_equal(b$lambda_common,a$lambda_common)
    expect_equal(b$null_logLik,a$null_logLik-(26-2*reml)*660*log(2),tolerance=1e-10)
    # Bootstrap is unit-free: fixed reference and random draws remain identical.
    b$null_logLik <- a$null_logLik
    expect_equal(b,a,tolerance=0)
  }
})