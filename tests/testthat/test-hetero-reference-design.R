.count_oracle <- function(y,k) {
  count<-rep(1,k);total<-k;answer<-0
  for(state in y[!is.na(y)]) {
    raw<-previous<-0;cdf<-numeric(k)
    for(j in seq_len(k)) {
      raw<-raw+count[j]/total
      previous<-max(previous,min(1,max(0,if(j==k)1 else raw)));cdf[j]<-previous
    }
    answer<-answer+log(cdf[state]-if(state>1)cdf[state-1]else 0)
    count[state]<-count[state]+1;total<-total+1
  }
  answer
}

test_that("count proposal densities normalize and match independent sequential kernels", {
  for(k in c(2L,3L,8L)) {
    x<-t(as.matrix(expand.grid(rep(list(seq_len(k)),3L))));storage.mode(x)<-'integer'
    x<-rbind(x[1,,drop=FALSE],NA_integer_,x[2:3,,drop=FALSE])
    set.seed(2020);seed<-.Random.seed
    fit<-heterogeneity_count_loglik_cpp(x,k)
    expected<-apply(x,2,.count_oracle,k=k)
    expect_true(all(fit[,1]<=expected+1e-13));expect_true(all(fit[,2]>=expected-1e-13))
    expect_lte(sum(exp(fit[,1])),1+1e-13);expect_gte(sum(exp(fit[,2])),1-1e-13)
    expect_identical(.Random.seed,seed)
    order<-rev(seq_len(ncol(x)))
    expect_identical(fit[order,,drop=FALSE],heterogeneity_count_loglik_cpp(x[,order,drop=FALSE],k))
    ideal<-apply(x,2,function(y)lgamma(k)+sum(lgamma(tabulate(y,k)+1))-lgamma(sum(!is.na(y))+k))
    expect_lt(max(abs(rowMeans(fit)-ideal)),1e-11)
  }
  x<-matrix(c(1L,NA_integer_,2L,2L),2)
  expect_error(heterogeneity_count_loglik_cpp(x,3L),"fixed missing mask")
  expect_error(heterogeneity_count_loglik_cpp(matrix(4L,2,2),3L),"invalid count proposal state")
})

test_that("count simulation has the stated RNG kernel and resets each draw", {
  mask<-rep(c(TRUE,FALSE,TRUE),7)
  for(k in c(2L,3L,8L)) {
    set.seed(2021);u<-runif(sum(mask));expected_seed<-.Random.seed
    count<-rep(1,k);answer<-rep(NA_integer_,length(mask));position<-0L
    for(i in which(mask)) {
      position<-position+1L;raw<-previous<-0;cdf<-numeric(k)
      for(j in seq_len(k)) {raw<-raw+count[j]/sum(count);previous<-max(previous,min(1,max(0,if(j==k)1 else raw)));cdf[j]<-previous}
      state<-which(u[position]<cdf)[1L];answer[i]<-state;count[state]<-count[state]+1
    }
    set.seed(2021);actual<-heterogeneity_count_simulate_cpp(k,mask)
    expect_identical(actual,answer);expect_identical(.Random.seed,expected_seed)
    set.seed(2021);expect_identical(heterogeneity_count_simulate_cpp(k,mask),answer)
  }
})

test_that("between-crown Delta reuses dispersion and leaves the global baseline out", {
  set.seed(2022);seed<-.Random.seed
  for(n in 2:12) {
    x<-exp(seq(-2,2,length.out=n))
    expect_equal(heterogeneity_delta_between_cpp(x),stats::sd(log(x)),tolerance=1e-14)
    expect_equal(heterogeneity_delta_between_cpp(x*17),heterogeneity_delta_between_cpp(x),tolerance=1e-14)
  }
  expect_true(is.na(heterogeneity_delta_between_cpp(c(1,NA_real_))))
  expect_identical(.Random.seed,seed)
})

test_that("public statistic changes retain estimators, reference bank and random stream", {
  tree<-ape::stree(16,type='balanced');tree$edge.length<-seq(.25,1.75,length.out=30)
  x<-setNames(rep(c('A','B','C'),length.out=16),tree$tip.label)
  nodes<-tree$edge[tree$edge[,1]==17,2]
  run<-function(statistic=NULL,counts=NULL) {
    set.seed(2023)
    control<-list(model='ER',reference_test='nuisance',state_space=c('A','B','C'),
      mcmc_sim=1200L,thin=5L,burn=100L,nuisance_control=list(max_boxes=2L))
    if(!is.null(statistic))control$statistic<-statistic
    if(!is.null(counts))control$nuisance_control$proposal_counts<-counts
    fit<-suppressWarnings(fast_hetero(tree,x,method='Delta',nodes=nodes,min_tips=5,
      test=TRUE,nsim=9,keep_null=TRUE,progress=FALSE,delta_control=control))
    list(fit=fit,seed=.Random.seed,uniform=runif(5),normal=rnorm(5))
  }
  old<-run();explicit<-run('local_global',FALSE);between<-run('between_clades')
  for(other in list(explicit,between)) {
    expect_identical(other$fit$clades,old$fit$clades)
    expect_identical(other$fit$overall,old$fit$overall)
    expect_identical(other$seed,old$seed);expect_identical(other$uniform,old$uniform);expect_identical(other$normal,old$normal)
    a<-old$fit$heterogeneity$reference_model[[1]];b<-other$fit$heterogeneity$reference_model[[1]]
    for(name in c('reference_patterns','reference_diagnostics','reference_unavailable_details','proposal_log_density_bounds'))
      expect_identical(a[[name]],b[[name]])
  }
  expect_identical(explicit$fit$null,old$fit$null)
  h<-between$fit$heterogeneity
  expect_equal(h$statistic,sd(log(between$fit$clades$estimate)),tolerance=1e-14)
  expect_identical(h$test_statistic,'SD(log(Delta_local))')
  count<-run('between_clades',TRUE)
  m<-count$fit$heterogeneity$reference_model[[1]]
  expect_identical(m$proposal_counts$observed_counts_used,FALSE)
  expect_equal(m$proposal_counts$weight,.5)
  expect_true(all(is.finite(m$proposal_log_density_bounds)))
  expect_equal(ncol(m$reference_patterns),10L)
  expect_error(.heterogeneity_delta_control(list(statistic='between_clades')),"require reference_test")
  expect_error(.heterogeneity_nuisance_control(list(proposal_counts=1)),"TRUE or FALSE")
})
