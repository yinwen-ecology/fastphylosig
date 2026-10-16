.hetero_dense_reml <- function(tree,x,lambda) {
  C <- ape::vcv.phylo(tree)[tree$tip.label,tree$tip.label]
  C <- lambda*(C-diag(diag(C)))+diag(diag(C))
  tryCatch({
    inv <- solve(C);precision <- sum(inv)
    mu <- sum(inv %*% x)/precision
    df <- length(x)-1L
    s2 <- as.numeric(crossprod(x-mu,inv %*% (x-mu)))/df
    -.5*(df*(log(2*pi)+log(s2)+1)+as.numeric(determinant(C,TRUE)$modulus)+log(precision))
  },error=function(e)-Inf)
}

test_that("shared REML reference matches dense covariance likelihood", {
  set.seed(901)
  trees <- lapply(1:3,function(i)ape::rtree(15+i))
  values <- lapply(trees,function(s)rnorm(length(s$tip.label)))
  cores <- lapply(trees,function(s)compile_tree_cpp(s$edge,s$edge.length,length(s$tip.label)))
  objective <- function(lambda)sum(vapply(seq_along(trees),function(i)
    .hetero_dense_reml(trees[[i]],values[[i]],lambda),numeric(1)))
  expected <- optimize(objective,c(0,1),maximum=TRUE,tol=1e-8)
  candidates <- c(0,expected$maximum,1)
  lambda <- candidates[which.max(vapply(candidates,objective,numeric(1)))]
  actual <- heterogeneity_lambda_calibrate_cpp(cores,values,.2,5L,FALSE,TRUE)
  expect_equal(actual$lambda_common,lambda,tolerance=2e-6)
  expect_equal(actual$null_logLik,objective(lambda),tolerance=1e-8)
  expect_identical(actual$reference_fit,"REML")
  fixed <- heterogeneity_lambda_calibrate_cpp(cores,values,.2,5L,FALSE,FALSE,.3)
  expect_equal(fixed$lambda_common,.3)
  expect_error(heterogeneity_lambda_calibrate_cpp(cores,values,.2,5L,FALSE,FALSE,-.1),"domain")
})

test_that("public heterogeneity uses REML reference without changing local lambda", {
  set.seed(902)
  tree <- ape::stree(64,"balanced");tree$edge.length <- rep(1,nrow(tree$edge))
  x <- setNames(rnorm(64),tree$tip.label)
  set.seed(903)
  reml <- fast_hetero(tree,x,method="lambda",min_tips=20,max_tips=32,
    test=TRUE,nsim=19,keep_null=TRUE,progress=FALSE)
  rng <- .Random.seed
  set.seed(903)
  ml <- fast_hetero(tree,x,method="lambda",min_tips=20,max_tips=32,
    test=TRUE,nsim=19,keep_null=TRUE,progress=FALSE,lambda_reference="ML")
  expect_identical(.Random.seed,rng)
  expect_identical(ml$clades$estimate,reml$clades$estimate)
  expect_identical(reml$heterogeneity$lambda_reference,"REML")
  expect_identical(ml$heterogeneity$lambda_reference,"ML")
})
