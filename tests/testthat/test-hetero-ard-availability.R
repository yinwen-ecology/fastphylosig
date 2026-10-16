test_that("reference tail bounds retain every requested draw and no RNG", {
  set.seed(432);seed <- .Random.seed
  for(draws in list(c(1,2,3),c(1,NA,3),c(NA,NA),c(Inf,NaN,2),numeric(),1)) {
    b <- heterogeneity_bootstrap_bounds_cpp(draws,2)
    if(length(draws)>=2L) {
      count <- sum(draws[is.finite(draws)]>=2);missing <- sum(!is.finite(draws))
      expect_identical(b$finite_exceedance_count,as.integer(count))
      expect_equal(b$P_heterogeneity_lower,(1+count)/(length(draws)+1))
      expect_equal(b$P_heterogeneity_upper,(1+count+missing)/(length(draws)+1))
    } else expect_true(is.na(b$P_heterogeneity_upper))
  }
  expect_true(is.na(heterogeneity_bootstrap_bounds_cpp(c(1,2),NA_real_)$P_heterogeneity_lower))
  expect_identical(.Random.seed,seed)
})

test_that("native reference flags reproduce isTRUE and finite-estimate semantics", {
  make <- function(value=1)list(estimate=value,single_state=FALSE,
    details=c(list(MCMC_warning=FALSE),.delta_calculation_flags()))
  fields <- c("MCMC_warning","optimizer_warning","numerical_recovery","entropy_recovery",
    "other_warning","calculation_failed")
  set.seed(435);seed <- .Random.seed
  for(field in fields) for(value in list(TRUE,FALSE,NA,1,NULL,c(TRUE,FALSE))) {
    f <- make();f$details[[field]] <- value
    flags <- heterogeneity_delta_flags_cpp(make(),list(f,make()))
    expect_identical(unname(flags[seq_along(fields)]),as.integer(vapply(fields,
      function(nm)isTRUE(f$details[[nm]]),logical(1))))
  }
  f <- make(NA_real_);f$single_state <- TRUE
  flags <- heterogeneity_delta_flags_cpp(make(),list(f))
  expect_identical(unname(flags[7:8]),c(1L,1L));expect_identical(.Random.seed,seed)
})

test_that("conservative single-state policy preserves draws and rejects actual errors", {
  z <- .hetero_fixture();x <- setNames(rep(1:3,length.out=64),z$tree$tip.label)
  control <- .heterogeneity_delta_control(list(model="ER",mcmc_sim=80L,thin=2L,burn=10L))
  h <- suppressWarnings(fast_hetero(z$tree,x,method="Delta",nodes=z$nodes,
    min_tips=20,overall=TRUE,delta_control=control,progress=FALSE))
  X <- matrix(x,ncol=1,dimnames=list(z$tree$tip.label,"trait_1"))
  attr(X,"state_labels") <- list(as.character(1:3))
  core <- compile_tree_cpp(z$tree$edge,z$tree$edge.length,64L)
  index <- 0L;actual_error <- FALSE
  testthat::local_mocked_bindings(.heterogeneity_delta_evaluate=function(...) {
    index <<- index+1L;b <- (index-1L)%/%3L+1L
    f <- .delta_calculation_flags();f$calculation_failed <- actual_error && b==3L
    list(estimate=if(b==3L)NA_real_ else 1,single_state=!actual_error&&b==3L,
      details=c(list(MCMC_warning=FALSE),f))
  })
  run <- function(policy) {
    index <<- 0L;set.seed(431);ctrl <- c(control,list(reference_unavailable=policy))
    out <- suppressWarnings(.heterogeneity_delta_tests(z$tree,core,rep(list(list(list())),2),X,
      z$nodes,h$clades,h$overall,ctrl,20L,TRUE,5L,TRUE,FALSE,1L))
    list(value=out,seed=.Random.seed)
  }
  strict <- run("strict");conservative <- run("conservative")
  expect_identical(strict$seed,conservative$seed)
  expect_identical(strict$value$null,conservative$value$null)
  a <- strict$value$summary;b <- conservative$value$summary
  expect_true(is.na(a$P_heterogeneity));expect_identical(b$status,"conservative_reference")
  expect_identical(b$P_heterogeneity,b$P_heterogeneity_upper)
  expect_equal(b$P_heterogeneity_upper-b$P_heterogeneity_lower,1/6)
  expect_identical(b$nsim_single_state,1L)
  figure <- h
  for(nm in names(b))figure$heterogeneity[[nm]] <- b[[nm]]
  figure$heterogeneity$nsim_requested <- 5L
  text <- .heterogeneity_summary_text(figure,"trait_1")
  expect_true(any(grepl("conservative upper-tail",text,fixed=TRUE)))
  expect_true(any(grepl("4/5 finite simulations",text,fixed=TRUE)))
  expect_true(any(grepl("Reference tail bounds",text,fixed=TRUE)))
  actual_error <- TRUE;failed <- run("conservative")$value$summary
  expect_true(is.na(failed$P_heterogeneity));expect_identical(failed$status,"simulation_failure")
  expect_identical(failed$nsim_calculation_failure,1L)
})

test_that("ARD multistart is deterministic, retains the best objective and reports nonconvergence", {
  workspace <- list(nl=3L,edge_length=c(.3,.7,1.4))
  dev <- function(p)sum((p-c(.2,.4))^2)
  set.seed(434);seed <- .Random.seed
  fit <- .fast_ace_multistart(dev,c(.1,.1),workspace)
  expect_lte(fit$objective,fit$diagnostics$baseline_deviance)
  expect_equal(fit$par,c(.2,.4),tolerance=1e-5)
  expect_identical(.Random.seed,seed)
  expect_equal(fit$objective,min(fit$diagnostics$attempts$objective))
  expect_identical(nrow(fit$diagnostics$attempts),4L)
  # A nonconverged best attempt must not be replaced by an inferior converged one.
  call <- 0L
  testthat::local_mocked_bindings(nlminb=function(...) {
    call <<- call+1L
    list(par=c(.2,.4),objective=as.numeric(call),convergence=if(call==1L)1L else 0L,
      iterations=1L,evaluations=c("function"=1L,"gradient"=1L),message="fixture")
  },.package="stats")
  best <- .fast_ace_multistart(dev,c(.1,.1),workspace)
  expect_identical(best$convergence,1L);expect_identical(best$objective,1)
  expect_identical(nrow(best$diagnostics$attempts),5L)
})

test_that("specialist controls stay out of observed Delta estimators", {
  expect_error(.heterogeneity_delta_control(list(model="ER",reference_fit="multistart")),"ARD")
  expect_error(.heterogeneity_delta_control(list(reference_unavailable="discard")),"arg")
  tr <- ape::stree(32,"balanced");tr$edge.length <- seq(.3,1.5,length.out=62)
  x <- setNames(rep(1:3,length.out=32),tr$tip.label)
  base <- list(model="ARD",mcmc_sim=80L,thin=2L,burn=10L)
  run <- function(ctrl) {
    set.seed(436);h <- suppressWarnings(fast_hetero(tr,x,method="Delta",min_tips=8L,
      max_tips=16L,overall=TRUE,delta_control=ctrl,progress=FALSE))
    list(clades=h$clades,overall=h$overall,seed=.Random.seed)
  }
  expect_identical(run(base),run(c(base,list(reference_fit="multistart",reference_unavailable="conservative"))))
})
