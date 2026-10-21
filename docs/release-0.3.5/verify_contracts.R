# Cross-install result/RNG contract and fixed-workload timing.
args <- commandArgs(TRUE)
out <- normalizePath(args[[1L]], winslash="/",mustWork=TRUE)
label <- args[[2L]]
if(length(args)>2L) .libPaths(c(normalizePath(args[[3L]],winslash="/",mustWork=TRUE),.libPaths()))
library(fastphylosig); library(ape)
ns <- asNamespace("fastphylosig")
numeric_result <- function(x) {
  if(is.numeric(x)||is.logical(x)||is.character(x))return(x)
  if(!is.list(x))return(NULL)
  x <- x[setdiff(names(x),c("timing","call"))]
  lapply(x,numeric_result)
}
capture <- function(fun,seed=70123L) {
  set.seed(seed); warnings <- character()
  value <- withCallingHandlers(fun(),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  })
  list(value=value,numeric=numeric_result(value),warnings=warnings,seed=.Random.seed)
}
set.seed(4001)
tree <- rtree(128); tree$edge.length <- exp(runif(nrow(tree$edge),-2,1))
ctx <- prepare_tree(tree)
X <- matrix(rnorm(128*3),128,3,dimnames=list(tree$tip.label,c("x","y","z")))
X[c(2,19,40),2] <- NA;X[c(3,20,90),3] <- NA
D <- X;D[] <- rep(c(0,1,0,0),length.out=length(D));D[is.na(X)] <- NA
Y <- X;Y[] <- rep(letters[1:3],length.out=length(Y));Y[is.na(X)] <- NA
y <- setNames(rep(letters[1:3],length.out=128),tree$tip.label)
nodes <- tree$edge[tree$edge[,1L]==129L,2L]
cases <- list(
  K_raw=function()fast_k(tree,X,test=TRUE,nsim=31,return_sim=TRUE,progress=FALSE),
  K_prepared=function()fast_k(ctx,X,test=TRUE,nsim=31,return_sim=TRUE,progress=FALSE),
  lambda=function()fast_lambda(tree,X,lambda_profile=TRUE,lambda_profile_points=31,progress=FALSE),
  D_raw=function()fast_d(tree,D,nsim=31,return_sim=TRUE,progress=FALSE),
  D_prepared=function()fast_d(ctx,D,nsim=31,return_sim=TRUE,progress=FALSE),
  Delta_ER_vector=function()fast_delta(tree,y,model="ER",progress=FALSE),
  Delta_ER_batch=function()fast_delta(tree,Y,model="ER",progress=FALSE),
  Delta_ER_prepared=function()fast_delta(ctx,Y,model="ER",progress=FALSE),
  Delta_ER_SE=function()fast_delta(tree,y,model="ER",entropy="SE",progress=FALSE),
  Delta_ER_GINI=function()fast_delta(tree,y,model="ER",entropy="GINI",progress=FALSE),
  Delta_ER_permutation=function()fast_delta(tree,y,model="ER",test=TRUE,nsim=9,
    return_sim=TRUE,progress=FALSE),
  hetero_K=function()fast_hetero(tree,X,method="K",nodes=nodes,min_tips=10,
    overall=TRUE,overall_test=TRUE,clade_test=TRUE,test=TRUE,nsim=9,keep_null=TRUE,progress=FALSE),
  hetero_lambda=function()fast_hetero(tree,X,method="lambda",nodes=nodes,min_tips=10,
    overall=TRUE,test=TRUE,nsim=9,keep_null=TRUE,progress=FALSE),
  hetero_D=function()fast_hetero(tree,D,method="D",nodes=nodes,min_tips=10,
    overall=TRUE,test=TRUE,nsim=5,d_nsim=31,keep_null=TRUE,progress=FALSE),
  hetero_Delta_ER=function()fast_hetero(tree,Y,method="Delta",nodes=nodes,min_tips=10,
    overall=TRUE,overall_test=TRUE,clade_test=TRUE,nsim=5,
    delta_control=list(model="ER"),progress=FALSE),
  ACE_ER_point=function(){f<-fast_ace(y,tree,model="ER",progress=FALSE);f$se<-NULL;f},
  ACE_ARD_point=function(){f<-fast_ace(y,tree,model="ARD",progress=FALSE);
    f[c("rates","loglik","convergence","message","index.matrix")]}
)

cases <- c(cases,list(
  Delta_ARD_vector=function()fast_delta(tree,y,model="ARD",progress=FALSE),
  Delta_ARD_batch=function()fast_delta(tree,Y,model="ARD",progress=FALSE),
  Delta_ARD_prepared=function()fast_delta(ctx,Y,model="ARD",progress=FALSE),
  Delta_ARD_SE=function()fast_delta(tree,y,model="ARD",entropy="SE",progress=FALSE),
  Delta_ARD_GINI=function()fast_delta(tree,y,model="ARD",entropy="GINI",progress=FALSE),
  Delta_ARD_permutation=function()fast_delta(tree,y,model="ARD",test=TRUE,nsim=5,
    return_sim=TRUE,progress=FALSE),
  hetero_Delta_ARD=function()fast_hetero(tree,Y,method="Delta",nodes=nodes,min_tips=10,
    overall=TRUE,overall_test=TRUE,clade_test=TRUE,nsim=5,
    delta_control=list(model="ARD"),progress=FALSE),
  ACE_ARD_full=function()suppressWarnings(fast_ace(y,tree,model="ARD",progress=FALSE))
))

results <- lapply(cases,function(fun)capture(fun))
saveRDS(results,file.path(out,paste0(label,"-contracts.rds")))
if(label!="baseline") {
  baseline <- readRDS(file.path(out,"baseline-contracts.rds"))
  comparison <- do.call(rbind,lapply(names(results),function(name)data.frame(case=name,
    result_identical=identical(results[[name]]$numeric,baseline[[name]]$numeric,num.eq=FALSE,single.NA=FALSE),
    rng_identical=identical(results[[name]]$seed,baseline[[name]]$seed),
    warnings_identical=identical(results[[name]]$warnings,baseline[[name]]$warnings))))
  print(comparison)
  write.csv(comparison,file.path(out,paste0(label,"-contracts.csv")),row.names=FALSE)
  stopifnot(all(comparison$result_identical),all(comparison$rng_identical),
            all(comparison$warnings_identical))
}
cat("Contract capture:",label,length(results),"cases\n")
