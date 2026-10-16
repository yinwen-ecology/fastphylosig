test_that("bounded ARD recovery agrees with independent likelihoods at matching rates", {
  skip_if_not_installed("expm")
  fixtures <- readRDS(test_path("fixtures","ace-ard-local-recovery.rds"))
  set.seed(731);seed <- .Random.seed
  for(z in fixtures) {
    a <- z$original;notes <- list()
    f <- withCallingHandlers(.fast_ace_fit_workspace(a$tip_state,a$lvls,a$workspace,.1,
      CI=TRUE,estimate_se=FALSE,recover_ard=TRUE),warning=function(w) {
        notes[[length(notes)+1L]] <<- w;invokeRestart("muffleWarning")
      })
    expect_identical(f$convergence,0L,info=a$key)
    # The saved oracle belongs to z$fit$par, not to a newly optimized fit.
    saved_reference <- .ace_recovery_reference_deviance(a, z$fit$par)
    saved_native <- fast_ace_discrete_deviance_cpp(a$workspace$edge,
      a$workspace$edge_length, a$tip_state, a$workspace$rate$rate_index,
      z$fit$par, TRUE)
    expect_equal(saved_reference, z$independent, tolerance=1e-8, info=a$key)
    expect_equal(as.numeric(saved_native), saved_reference,
      tolerance=1e-8, info=a$key)
    current_reference <- .ace_recovery_reference_deviance(a, f$rates)
    expect_equal(-2*f$loglik, current_reference, tolerance=1e-8, info=a$key)
    expect_true(all(is.finite(f$rates) & f$rates >= 0), info=a$key)
    expect_true(all(is.finite(f$lik.anc)&f$lik.anc>=0&f$lik.anc<=1))
    expect_lt(max(abs(rowSums(f$lik.anc)-1)),1e-12)
    expect_false(any(vapply(notes,inherits,logical(1),"fastphylosig_optimizer_warning")))
    # A platform may converge on its first attempt. Validate that answer above;
    # only a real retry has an unfinished baseline and a retry warning.
    retry <- f$diagnostics$optimizer_retry
    expect_identical(any(vapply(notes,inherits,logical(1),
      "fastphylosig_optimizer_retry_warning")), !is.null(retry))
    if(!is.null(retry)) {
      baseline_reference <- .ace_recovery_reference_deviance(a,retry$original$par)
      expect_lt(current_reference,baseline_reference)
      expect_identical(retry$original$convergence,1L)
      expect_true(isTRUE(retry$recovered))
    }
  }
  expect_identical(.Random.seed,seed)
})

test_that("captured unfinished ARD fits exercise recovery on every platform", {
  skip_if_not_installed("expm")
  fixtures <- readRDS(test_path("fixtures","ace-ard-local-recovery.rds"))
  recover <- function(a) {
    w <- a$workspace
    ptr <- fast_ace_fit_prepare_cpp(w$edge,w$edge_length,a$tip_state,
      w$rate$rate_index,1L)
    on.exit(fast_ace_fit_release_cpp(ptr))
    dev <- function(p) fast_ace_fit_deviance_cpp(ptr,p,TRUE)
    first <- a$fit
    expect_identical(first$convergence,1L)
    first$objective <- dev(first$par)
    fit <- .fast_ace_recover_ard(dev,first,recheck=TRUE)
    reference <- .ace_recovery_reference_deviance(a,fit$par)
    expect_identical(fit$convergence,0L,info=a$key)
    expect_equal(fit$objective,reference,tolerance=1e-8,info=a$key)
    expect_lt(reference,.ace_recovery_reference_deviance(a,first$par))
    expect_true(isTRUE(fit$diagnostics$optimizer_retry$recovered))
  }
  set.seed(733);seed <- .Random.seed
  for(z in fixtures) recover(z$original)
  expect_identical(.Random.seed,seed)
})

test_that("successful fits bypass recovery without another objective evaluation", {
  f <- list(par=c(.1,.2),objective=12,convergence=0L,message="converged")
  set.seed(732);seed <- .Random.seed
  expect_identical(.fast_ace_recover_ard(function(p)stop("unexpected evaluation"),f),f)
  expect_identical(.Random.seed,seed)
})

test_that("candidate selection retains a materially better unfinished fit", {
  a <- list(list(valid=TRUE,fit=list(objective=10,convergence=1L)),
    list(valid=TRUE,fit=list(objective=10.01,convergence=0L)))
  expect_identical(.fast_ace_select_candidate(a)$selected,1L)
  a[[2L]]$fit$objective <- 10+1e-8
  expect_identical(.fast_ace_select_candidate(a)$selected,2L)
})

test_that("failed recovery attempts retain the original unfinished fit and diagnostics", {
  testthat::local_mocked_bindings(nlminb=function(...)stop("nlminb failure"),
    optim=function(...)stop("optim failure"),.package="stats")
  f <- list(par=c(.1,.2),objective=12,convergence=1L,iterations=150L,
    evaluations=c("function"=160L),message="iteration limit")
  a <- .fast_ace_recover_ard(function(p)12,f)
  expect_identical(a$par,f$par)
  expect_identical(a$convergence,1L)
  expect_false(a$diagnostics$optimizer_retry$recovered)
  expect_identical(a$diagnostics$attempts$valid,c(TRUE,FALSE,FALSE))
})

test_that("optimizer retry has separate Delta and native summary diagnostics", {
  testthat::local_mocked_bindings(.fast_ace_lik_anc_workspace=function(...) {
    warning(structure(list(message="recovered",call=NULL,recovered=TRUE,
      diagnostics=list(method="fixture")),
      class=c("fastphylosig_optimizer_retry_warning","warning","condition")))
    matrix(c(.3,.7),3L,2L,byrow=TRUE)
  })
  tree <- list(tip.label=c("a","b","c","d"))
  y <- setNames(factor(c(1,1,2,2)),tree$tip.label)
  z <- suppressWarnings(.delta_one(tree,y,.1,.5,100L,10L,10L,1L,"ARD","fast",
    ace_workspace=list(),keep_chains=FALSE))$calculation_diagnostics
  expect_true(z$optimizer_retry)
  expect_true(z$optimizer_recovered)
  expect_false(z$optimizer_warning)
  expect_false(z$numerical_recovery)
  expect_false(z$other_warning)
  fit <- list(estimate=1,details=z,single_state=FALSE)
  flags <- heterogeneity_delta_flags_cpp(fit,list(fit,fit))
  expect_identical(unname(flags["nsim_optimizer_warning"]),0L)
  expect_identical(attr(flags,"optimizer_retry"),c(nsim_optimizer_retry=1L,nsim_optimizer_recovered=1L))
  healthy <- list(estimate=1,details=.delta_calculation_flags(),single_state=FALSE)
  expect_null(attr(heterogeneity_delta_flags_cpp(healthy,list(healthy)),"optimizer_retry"))
})
