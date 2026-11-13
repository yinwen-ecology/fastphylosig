test_that("cached search knots retain a narrow excluded LR gap near the boundary", {
  tr <- ape::read.tree(text = paste0(
    "(((t6:2.7223428872131088e-06,t1:1.2471560969957194e-05):0.010410820197158256,",
    "((t8:0.00043189953810771715,t2:40.374962856766381):18.96253731901621,",
    "t5:7.2376568355125297):902.55227937619657):1.3319443434518778,",
    "((t7:1.2576390145180813e-05,t3:123.06603191888971):4.3185557773381258e-06,",
    "t4:0.0032466927237239139):0.61025748467609153);"))
  x <- c(t6=.491690144941276,t1=.493796135038392,t8=6.52279701497337,
         t2=18.8043946827766,t5=17.7865960120158,t7=.222616579286159,
         t3=11.5361908541573,t4=.341041097705755)
  fit <- fast_lambda(tr, x, test=TRUE, progress=FALSE)
  dense <- fast_lambda(tr, x, test=TRUE, lambda_profile_points=5001, progress=FALSE)
  ints <- fit$lambda_CI_intervals
  expect_equal(nrow(ints), 2L)
  expect_true(fit$lambda_CI_disconnected)
  expect_equal(ints, dense$lambda_CI_intervals, tolerance=1e-9)
  # Independent dense Gaussian likelihood verifies all three finite crossings.
  V <- ape::vcv(tr); y <- x[rownames(V)]
  ell <- function(l) {
    C <- V*l; diag(C) <- diag(V); iv <- solve(C)
    mu <- sum(iv %*% y)/sum(iv); r <- y-mu
    s2 <- drop(t(r)%*%iv%*%r)/length(y)
    -.5*(length(y)*(log(2*pi)+1+log(s2)) +
           as.numeric(determinant(C, logarithm=TRUE)$modulus))
  }
  cut <- fit$lambda_CI_cutoff
  expect_lt(ell(.995), cut)
  expect_gt(ell(.9999), cut)
  roots <- vapply(list(c(.5,.8),c(.99,.995),c(.999,.9999)), function(b) {
    uniroot(function(z) ell(z)-cut, b, tol=1e-12)$root
  }, numeric(1))
  expect_equal(c(ints$lower[1], ints$upper[1], ints$lower[2]), roots, tolerance=1e-9)
  tab <- fast_lambda(prepare_tree(tr), cbind(x=x, twice=2*x),
                     test=TRUE, lambda_profile=TRUE, progress=FALSE)
  expect_equal(tab$lambda_CI_intervals_fast[[1]], ints, tolerance=1e-9)
  expect_equal(tab$lambda_CI_intervals_fast[[2]], ints, tolerance=1e-9)
  # Profiling enriches CI sampling, never the estimator or RNG consumption.
  set.seed(382)
  before <- .Random.seed
  no_profile <- fast_lambda(tr,x,test=TRUE,lambda_profile=FALSE,progress=FALSE)
  expect_identical(.Random.seed,before)
  expect_identical(fit$lambda,no_profile$lambda)
  expect_identical(fit$logL,no_profile$logL)
  grid <- fit$lambda_profile$lambda
  expect_true(all(c(2^-(7:40),1-2^-(7:40)) %in% grid))
  ref <- fastphylosig:::.lambda_profile_data_batch(
    1,fit$lambda,fit$logL,101,scalar_lik=function(z,y) ell(z),Y=matrix(x,ncol=1))[[1]]
  expect_equal(ref$lambda,grid)
  expect_equal(ref$logL,fit$lambda_profile$logL,tolerance=1e-8)
})

test_that("scalar profile fallback preserves grid-by-trait orientation", {
  fun <- fastphylosig:::.lambda_profile_data_batch
  for (p in 1:3) {
    Y <- matrix(seq_len(3*p),nrow=3)
    hats <- seq(.2,.8,length.out=p)
    vals <- colSums(Y)+hats
    scalar <- function(z,y) sum(y)+z
    batch <- function(z) outer(z,colSums(Y),`+`)
    a <- fun(1,hats,vals,5,scalar_lik=scalar,Y=Y)
    b <- fun(1,hats,vals,5,batch_lik=batch,scalar_lik=scalar,Y=Y)
    expect_equal(a,b)
    for (j in seq_len(p)) expect_equal(a[[j]]$logL,colSums(Y)[j]+a[[j]]$lambda)
  }
})
