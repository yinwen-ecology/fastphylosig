test_that('preflight rejects representations that V2 cannot prepare', {
  tree <- ape::read.tree(text='((a:1,b:1):1,c:2);')
  for (kind in c('factor','missing_Nnode','zero_Nnode')) {
    bad <- tree
    if(kind=='factor') bad$tip.label <- factor(tree$tip.label)
    if(kind=='missing_Nnode') bad$Nnode <- NULL
    if(kind=='zero_Nnode') bad$Nnode <- 0L
    snapshot <- serialize(bad,NULL)
    checked <- fastphylosig::check_tree(bad)
    expect_false(checked$ready)
    expect_false(any(checked$ready_by_signal))
    expect_true(any(checked$issues$severity=='ERROR'))
    expect_error(fastphylosig::prepare_tree(bad))
    expect_identical(serialize(bad,NULL),snapshot)
  }
  expect_true(fastphylosig::check_tree(tree,signal=c('K','lambda','D'))$ready)
})

test_that('explicit D probability columns survive differing aliases and overlays', {
  fit <- data.frame(trait='x',D_fast=.5,observed=2,mean_random=1,mean_brownian=3,
    P_random=.4,P_Brownian=.2,Pval1_fast=.7,Pval0_fast=.6,custom_P=.8)
  fit$random_fast <- I(list(c(0,1,2,2,3)))
  fit$brownian_fast <- I(list(c(1,2,2,3,4)))
  expect_equal(fastphylosig:::.signal_distribution_data(fit,p_col='Pval0_fast')$p_value,.6)
  expect_equal(fastphylosig:::.signal_distribution_data(fit,p_col='Pval1_fast')$p_value,.7)
  expect_equal(fastphylosig:::.signal_distribution_data(fit,p_col='custom_P')$p_value,.8)
  expect_equal(fastphylosig:::.signal_distribution_data(fit,null='brownian',
    p_col='Pval1_fast')$p_value,.2)
  expect_equal(fastphylosig:::.signal_distribution_data(fit,null='random',
    p_col='Pval0_fast')$p_value,.4)
  scalar <- structure(list(DEstimate=.5,binvar='x',P_random=.4,P_Brownian=.2,
    Pval1=.7,Pval0=.6,Parameters=list(Observed=2,MeanRandom=1,MeanBrownian=3),
    Permutations=list(random=c(0,1,2,2,3),brownian=c(1,2,2,3,4))),class='phylo.d')
  expect_equal(fastphylosig:::.signal_distribution_data(scalar,p_col='Pval0')$p_value,.6)
  expect_equal(fastphylosig:::.signal_distribution_data(scalar,p_col='Pval1')$p_value,.7)
  path <- tempfile(fileext='.pdf')
  grDevices::pdf(path)
  on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
  for (x in list(fit,scalar)) {
    alias <- if (is.data.frame(x)) 'Pval0_fast' else 'Pval0'
    plotted <- fastphylosig::plot_signal(x,p_col=alias)
    expect_equal(plotted$Pval0,.6)
    expect_equal(plotted$Pval1,.4)
    expect_equal(plotted$extreme_random,2)
    expect_equal(plotted$extreme_brownian,2)
    expect_equal(plotted$P_Brownian_display,fastphylosig:::.format_p_value(.6))
  }
})
