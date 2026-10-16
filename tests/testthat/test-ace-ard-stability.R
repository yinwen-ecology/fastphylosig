test_that("ill-conditioned ARD objectives agree with independent matrix exponentials", {
  fixtures <- readRDS(test_path("fixtures","ace-ard-stability.rds"))
  set.seed(719); seed <- .Random.seed
  for(z in fixtures) {
    tw <- .fast_ace_tree_workspace(z$tree)
    ws <- .fast_ace_workspace(tw,"ARD",5L)
    y <- setNames(z$values,z$tree$tip.label)[tw$tip_label]
    for(reuse in c(FALSE,TRUE)) {
      dev <- fast_ace_discrete_deviance_cpp(ws$edge,ws$edge_length,y,
        ws$rate$rate_index,z$rates,reuse)
      expect_equal(as.numeric(dev),z$objective,tolerance=1e-8,info=z$key)
      a <- fast_ace_discrete_liks_cpp(ws$edge,ws$edge_length,y,
        ws$rate$rate_index,z$rates,FALSE,reuse)
      expect_true(all(is.finite(a$lik.anc) & a$lik.anc>=0 & a$lik.anc<=1))
      expect_lt(max(abs(rowSums(a$lik.anc)-1)),1e-14)
      expect_equal(a$deviance,z$objective,tolerance=1e-8)
    }
  }
  expect_identical(.Random.seed,seed)
})

test_that("interior entropy cannot conceal illegal probabilities", {
  # Each of these has at least one normalized entropy inside [0,1].
  for(p in list(c(-.1,.5,.6),c(.2,.2,.2),c(.6,.6,0)))
    for(type in 1:3) {
      set.seed(720);seed <- .Random.seed
      expect_error(delta_entropy_cpp(matrix(p,1L),type),"probabilities are invalid")
      expect_identical(.Random.seed,seed)
    }
})

test_that("equivalent converged multistart candidate replaces an unfinished minimum", {
  call <- 0L
  testthat::local_mocked_bindings(nlminb=function(...) {
    call <<- call+1L
    list(par=c(.2,.4),objective=if(call==1L)100 else 100+1e-8,
      convergence=if(call==1L)1L else 0L,iterations=1L,
      evaluations=c("function"=1L,"gradient"=1L),message="fixture")
  },.package="stats")
  set.seed(721);seed <- .Random.seed
  fit <- .fast_ace_multistart(function(p)sum(p),c(.1,.1),
    list(nl=3L,edge_length=c(.3,.7,1.4)))
  expect_identical(fit$convergence,0L)
  expect_identical(fit$diagnostics$selected_attempt,2L)
  expect_identical(fit$diagnostics$selection_recovery$minimum_attempt,1L)
  expect_identical(.Random.seed,seed)
})
