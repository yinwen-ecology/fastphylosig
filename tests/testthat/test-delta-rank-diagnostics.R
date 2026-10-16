test_that("native rank diagnostics match independent R covariance calculations", {
  set.seed(883)
  for (n in c(8L, 31L, 400L)) {
    a <- matrix(exp(rnorm(n*2)), n, 2)
    b <- matrix(exp(rnorm(n*2)), n, 2)
    variants <- list(list(a,b), list(a,b+rep(c(0,4),each=n)),
      list(a,round(b,1)+.1), list(a,2*a),
      list(a,exp(apply(log(b),2,cumsum)/4)))
    for (values in variants) {
      rng <- .Random.seed
      native <- fastphylosig:::delta_diagnostics_cpp(values[[1]],values[[2]])
      reference <- fastphylosig:::.delta_diagnostics_fallback(values[[1]],values[[2]])
      expect_identical(.Random.seed,rng)
      for (nm in fastphylosig:::.delta_rank_fields())
        expect_equal(native[[nm]],reference[[nm]],tolerance=1e-9,info=nm)
    }
  }
})

test_that("rank diagnostics match pinned upstream posterior oracle fixtures", {
  fixtures <- readRDS(test_path("fixtures","delta-rank-oracle.rds"))
  for (z in fixtures) {
    native <- fastphylosig:::delta_diagnostics_cpp(z$alpha,z$beta)
    for (nm in names(z$expected))
      expect_equal(native[[nm]],z$expected[[nm]],tolerance=1e-8,info=nm)
  }
})

test_that("rank diagnostics detect scale mismatch and constant traces", {
  set.seed(884)
  a <- matrix(exp(rnorm(4000)),2000,2)
  b <- cbind(exp(rnorm(2000,0,.05)),exp(rnorm(2000,0,2)))
  d <- fastphylosig:::delta_diagnostics_cpp(a,b)
  expect_gt(d$rank_Rhat_beta,1.1)
  fixed <- fastphylosig:::delta_diagnostics_cpp(a,2*a)
  expect_true(is.na(fixed$rank_Rhat_Delta))
  expect_true(is.na(fixed$ESS_bulk_Delta))
})

test_that("ratio diagnostics enter warnings even with acceptable parameter diagnostics", {
  fields <- fastphylosig:::.delta_rank_fields()
  modern <- setNames(as.list(rep(1000,length(fields))),fields)
  modern[grep("rank_Rhat",fields,value=TRUE)] <- list(1)
  modern$relative_MCSE_Delta <- .01
  modern$ESS_tail_Delta <- 50
  z <- fastphylosig:::.delta_mcmc_warning_summary("x",2,1000,1000,1,1,.02,
    ess_delta=1000,rhat_delta=1,rank_diagnostics=modern)
  expect_true(z$warning_emitted)
  expect_match(z$by_trait$reason,"ESS_tail_Delta")
  expect_equal(z$rank_rhat_maximum,1.01)
})

test_that("public Delta and heterogeneity tables expose modern diagnostics", {
  set.seed(885)
  tree <- ape::stree(32,"balanced");tree$edge.length <- rep(1,nrow(tree$edge))
  x <- setNames(rep(letters[1:3],length.out=32),tree$tip.label)
  d <- suppressWarnings(fast_delta(tree,x,mcmc_sim=100,thin=2,burn=20,progress=FALSE))
  h <- suppressWarnings(fast_hetero(tree,x,method="Delta",min_tips=8,max_tips=16,overall=TRUE,
    delta_control=list(model="ER",mcmc_sim=100,thin=2,burn=20),progress=FALSE))
  fields <- fastphylosig:::.delta_rank_fields()
  expect_true(all(fields %in% names(d)))
  expect_true(all(fields %in% names(h$clades)))
  expect_true(all(fields %in% names(h$overall)))
})
