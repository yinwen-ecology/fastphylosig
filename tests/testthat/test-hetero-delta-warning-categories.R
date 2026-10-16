test_that("Delta recovery, optimizer warnings and errors have distinct structured flags", {
  tr <- ape::stree(4, "balanced"); tr$edge.length <- rep(1, 6)
  y <- factor(c(1, 2, 1, 2)); names(y) <- tr$tip.label
  p <- rbind(c(.5,.5)-.Machine$double.eps, c(.5,.5), c(.8,.2))
  testthat::local_mocked_bindings(.fast_ace_lik_anc_workspace = function(...) p)
  set.seed(416)
  fit <- suppressWarnings(.delta_one(tr, y, .1, .5, 80L, 2L, 10L,
    3L, "ER", "fast", ace_workspace = list()))
  f <- fit$calculation_diagnostics
  expect_true(f$entropy_recovery)
  expect_identical(f$entropy_recovery_rows, 1L)
  expect_false(f$optimizer_warning)
  expect_false(f$other_warning)
  expect_false(f$calculation_failed)

  testthat::local_mocked_bindings(.fast_ace_lik_anc_workspace = function(...) {
    warning(structure(list(message="optimizer",call=NULL,convergence=1L),
      class=c("fastphylosig_optimizer_warning","warning","condition")))
    warning(structure(list(message="recovery",call=NULL,evaluations=3),
      class=c("fastphylosig_numerical_recovery_warning","warning","condition")))
    p
  })
  set.seed(416)
  again <- suppressWarnings(.delta_one(tr, y, .1, .5, 80L, 2L, 10L,
    3L, "ER", "fast", ace_workspace=list()))
  expect_identical(again$delta, fit$delta)
  expect_true(again$calculation_diagnostics$optimizer_warning)
  expect_true(again$calculation_diagnostics$numerical_recovery)
  expect_identical(again$calculation_diagnostics$numerical_recovery_evaluations, 3)
  expect_false(again$calculation_diagnostics$calculation_failed)

  testthat::local_mocked_bindings(.fast_ace_lik_anc_workspace = function(...) stop("ACE failure"))
  e <- tryCatch(.delta_one(tr,y,.1,.5,80L,2L,10L,3L,"ER","fast",
    ace_workspace=list()),error=identity)
  expect_match(conditionMessage(e),"ACE failure")
  expect_true(e$calculation_diagnostics$calculation_failed)
  expect_identical(e$calculation_diagnostics$failure_stage,"ACE")
  expect_false(e$calculation_diagnostics$optimizer_warning)
  testthat::local_mocked_bindings(.fast_ace_lik_anc_workspace = function(...)
    matrix(c(-.1,1.1),nrow=1L))
  set.seed(417);seed <- .Random.seed
  e <- tryCatch(.delta_one(tr,y,.1,.5,80L,2L,10L,3L,"ER","fast",
    ace_workspace=list()),error=identity)
  expect_match(conditionMessage(e),"cannot be repaired")
  expect_identical(e$calculation_diagnostics$failure_stage,"MCMC")
  expect_true(e$calculation_diagnostics$calculation_failed)
  expect_false(e$calculation_diagnostics$optimizer_warning)
  expect_identical(.Random.seed,seed)
})

test_that("Delta reference counts flagged replicates without discarding failures", {
  z <- .hetero_fixture()
  x <- setNames(rep(1:3,length.out=64),z$tree$tip.label)
  control <- .heterogeneity_delta_control(list(model="ER",mcmc_sim=80L,thin=2L,burn=10L))
  h <- suppressWarnings(fast_hetero(z$tree,x,method="Delta",nodes=z$nodes,
    min_tips=20,overall=TRUE,delta_control=control,progress=FALSE))
  X <- matrix(x,ncol=1,dimnames=list(z$tree$tip.label,"trait_1"))
  attr(X,"state_labels") <- list(as.character(1:3))
  index <- 0L
  testthat::local_mocked_bindings(.heterogeneity_delta_evaluate = function(...) {
    index <<- index+1L; b <- (index-1L)%/%3L+1L
    flags <- .delta_calculation_flags()
    if(b==1L) flags$entropy_recovery <- TRUE
    if(b==2L) flags$optimizer_warning <- TRUE
    if(b==3L) flags$numerical_recovery <- TRUE
    if(b==4L) { flags$calculation_failed <- TRUE; flags$failure_stage <- "MCMC" }
    list(estimate=if(b>=4L) NA_real_ else 1, fit_warning=flags$optimizer_warning,
      details=c(list(MCMC_warning=b==3L),flags))
  })
  core <- compile_tree_cpp(z$tree$edge,z$tree$edge.length,64L)
  out <- .heterogeneity_delta_tests(z$tree,core,rep(list(list(list())),2),X,
    z$nodes,h$clades,h$overall,control,20L,TRUE,5L,TRUE,FALSE,1L)$summary
  expect_identical(out$nsim_entropy_recovery,1L)
  expect_identical(out$nsim_numerical_recovery,1L)
  expect_identical(out$nsim_optimizer_warning,1L)
  expect_identical(out$nsim_fit_warning,out$nsim_optimizer_warning)
  expect_identical(out$nsim_MCMC_warning,1L)
  expect_identical(out$nsim_calculation_failure,1L)
  expect_identical(out$nsim_unavailable,2L)
  expect_identical(out$nsim_successful,3L)
  expect_true(is.na(out$P_heterogeneity))
  expect_identical(out$status,"simulation_failure")
})
