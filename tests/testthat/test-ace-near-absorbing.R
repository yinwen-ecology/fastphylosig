test_that("bounded recovery guards near-absorbing spectral likelihoods", {
  z <- readRDS(test_path("fixtures","ace-near-absorbing.rds"));w <- z$workspace
  set.seed(734);seed <- .Random.seed
  p <- fast_ace_fit_prepare_cpp(w$edge,w$edge_length,z$tip_state,w$rate$rate_index,1L)
  on.exit(fast_ace_fit_release_cpp(p),add=TRUE)
  # The independent constant uses exact IEEE inputs and 60/90-digit Decimal
  # Taylor/scaling/squaring, not the production eigen or expmat implementations.
  value <- fast_ace_fit_deviance_cpp(p,z$par,TRUE)
  expect_lt(abs(value-z$independent),1e-8)
  a <- fast_ace_fit_liks_cpp(p,z$par,FALSE,TRUE)
  expect_identical(a$deviance,value)
  expect_true(all(is.finite(a$lik.anc)&a$lik.anc>=0&a$lik.anc<=1))
  expect_lt(max(abs(rowSums(a$lik.anc)-1)),1e-12)
  expect_gt(fast_ace_fit_stability_cpp(p)[1L],0)
  expect_identical(.Random.seed,seed)
})

test_that("recovery retains the original objective when its guarded baseline changes", {
  testthat::local_mocked_bindings(nlminb=function(...)stop("unavailable"),
    optim=function(...)stop("unavailable"),.package="stats")
  first <- list(par=c(.1,.2),objective=12,convergence=1L,message="unfinished",
    iterations=150L,evaluations=c("function"=160L))
  a <- .fast_ace_recover_ard(function(p)12.001,first,recheck=TRUE)
  expect_identical(a$objective,12.001)
  expect_identical(a$convergence,1L)
  expect_identical(a$diagnostics$optimizer_retry$original,first)
  expect_identical(a$diagnostics$baseline_deviance,12.001)
  expect_equal(a$diagnostics$initial_numerical_correction,-.001)
  expect_error(.fast_ace_recover_ard(function(p)1e50,first,recheck=TRUE),
    "no numerically valid initial deviance")
})
