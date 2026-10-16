test_that("native scalar Delta warnings reproduce every R reason and threshold", {
  d <- list(ESS_alpha=1000,ESS_beta=1000,ESS_Delta=1000,
    split_Rhat_alpha=1,split_Rhat_beta=1,split_Rhat_Delta=1,MCSE_Delta=.001)
  for(nm in .delta_rank_fields())d[[nm]] <- if(startsWith(nm,"rank_Rhat"))1 else 1000
  d$relative_MCSE_Delta <- .001
  oracle <- function(d,delta=1) .delta_mcmc_warning_summary("reference",delta,
    d$ESS_alpha,d$ESS_beta,d$split_Rhat_alpha,d$split_Rhat_beta,d$MCSE_Delta,
    ess_delta=d$ESS_Delta,rhat_delta=d$split_Rhat_Delta,
    rank_diagnostics=d[.delta_rank_fields()])$by_trait$reason
  set.seed(914);seed <- .Random.seed
  for(delta in c(NA_real_,Inf,-Inf,NaN,0,-1,1))
    expect_identical(delta_warning_reason_cpp(delta,d),oracle(d,delta))
  for(nm in setdiff(names(d),"relative_MCSE_Delta")) {
    limit <- if(startsWith(nm,"rank_Rhat"))1.01 else if(startsWith(nm,"split_Rhat"))1.1 else
      if(startsWith(nm,"ESS_bulk")||startsWith(nm,"ESS_tail"))200 else if(nm=="MCSE_Delta").05 else 20
    for(value in c(NA_real_,NaN,Inf,-Inf,0,limit-.Machine$double.eps*max(1,limit),
        limit,limit+.Machine$double.eps*max(1,limit))) {
      z <- d;z[[nm]] <- value
      expect_identical(delta_warning_reason_cpp(1,z),oracle(z),info=paste(nm,value))
    }
  }
  for(value in c(NA_real_,NaN,Inf,-Inf,0,200,1.01)) {
    z <- lapply(d,function(x)value)
    expect_identical(delta_warning_reason_cpp(1,z),oracle(z))
  }
  expect_identical(.Random.seed,seed)
})

test_that("heterogeneity retains R warning text under custom numeric formatting", {
  tr <- ape::stree(8,"balanced");tr$edge.length <- rep(1,14)
  control <- .heterogeneity_delta_control(list(model="ER",mcmc_sim=80L,burn=10L,thin=2L))
  context <- .heterogeneity_delta_context(tr,tr$tip.label,tr$tip.label,control)
  testthat::local_mocked_bindings(delta_warning_reason_cpp=function(...)stop("native formatter used"))
  for(config in list(list(OutDec=","),list(scipen=-9))) {
    old <- options(config);set.seed(918)
    fit <- tryCatch(.heterogeneity_delta_values(context,rep(1:2,4)),finally=options(old))
    expect_identical(fit$status,"ok")
    expect_true(fit$details$MCMC_warning)
    expect_match(fit$message,"Delta MCMC:")
  }
})
