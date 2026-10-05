test_that('near-ultrametric lambda bounds agree with feasible edge constraints', {
  for (delta in c(0, 1e-9, 2e-8, 4e-8, .1)) {
    tree <- ape::read.tree(text = sprintf(
      '((a:1,b:1):1,(c:1,d:%.15g):1);', 1 + delta))
    ctx <- fastphylosig::prepare_tree(tree)
    tree <- ctx$tree
    group <- fastphylosig:::.prepared_tree_subset(ctx,need_matrix=FALSE)
    depth <- ape::node.depth.edgelength(tree)
    terminal <- tree$edge[,2] <= 4
    edge_bound <- min(depth[tree$edge[terminal,2]] / depth[tree$edge[terminal,1]])
    expected <- if (ape::is.ultrametric(tree)) edge_bound else 1
    upper <- fastphylosig:::.max_lambda(tree)
    expect_equal(upper,expected,tolerance=1e-14)
    expect_gte(min(depth[tree$edge[terminal,2]] -
                     upper * depth[tree$edge[terminal,1]]), -1e-14)
    x <- matrix(c(1,2,4,3),4)
    fixed <- fastphylosig:::fast_lambda_tree_fixed_cpp(group$compiled_tree,x,
      c(0,1,upper*(1-1e-7),upper,upper+1e-5))
    expect_equal(fixed$max_lambda,expected,tolerance=1e-14)
    expect_true(all(fixed$valid[1:3,]))
    expect_false(fixed$valid[5,1])
    expect_match(fixed$status[5,1],'outside')
    # An exact singular endpoint may be invalid; just inside it must match
    # independently constructed dense Pagel covariance likelihood.
    C <- ape::vcv.phylo(tree)
    for (i in 1:3) {
      lambda <- fixed$lambda[i]
      V <- C * lambda
      diag(V) <- diag(C)
      inverse <- solve(V)
      y <- x[,1]
      mu <- sum(inverse %*% y) / sum(inverse)
      variance <- as.numeric(crossprod(y-mu,inverse %*% (y-mu)))/4
      ll <- -.5*(4*(log(2*pi)+1+log(variance)) +
        as.numeric(determinant(V,logarithm=TRUE)$modulus))
      expect_equal(fixed$logLik[i,1],ll,tolerance=1e-7)
    }
    optimized <- fastphylosig:::fast_lambda_tree_optimize_cpp(
      group$compiled_tree,x,max_lambda=upper)
    expect_equal(optimized$max_lambda,upper,tolerance=1e-14)
    expect_true(optimized$valid[1])
    expect_lte(optimized$lambda[1],upper)
    oracle <- optimize(function(lambda) {
      V <- lambda*C; diag(V) <- diag(C)
      inverse <- solve(V)
      mu <- sum(inverse %*% x[,1])/sum(inverse)
      variance <- as.numeric(crossprod(x[,1]-mu,inverse %*% (x[,1]-mu)))/4
      -.5*(4*(log(2*pi)+1+log(variance)) +
        as.numeric(determinant(V,logarithm=TRUE)$modulus))
    }, c(0,upper), maximum=TRUE)
    expect_equal(optimized$logLik[1],oracle$objective,tolerance=1e-6)
    fit <- fastphylosig::fast_lambda(ctx,setNames(x[,1],tree$tip.label),
      test=FALSE,verbose=FALSE,progress=FALSE)
    expect_true(is.finite(fit$lambda))
    expect_lte(fit$lambda,upper)
  }
})

test_that('approximately ultrametric lambda upper bound never makes a terminal edge negative', {
  tree <- ape::read.tree(text='((a:1,b:1):1,(c:1,d:1.000000001):1);')
  upper <- fastphylosig:::.max_lambda(tree)
  expect_equal(upper,2)
  # phytools/ape classification is retained, but its max-height ratio can
  # exceed the exact feasible edge limit on slightly unequal tip heights.
  expect_true(ape::is.ultrametric(tree))
  depths <- ape::node.depth.edgelength(tree)
  legacy_upper <- max(depths[tree$edge[,2]]) / max(depths[tree$edge[,1]])
  expect_gt(legacy_upper,upper)
})
