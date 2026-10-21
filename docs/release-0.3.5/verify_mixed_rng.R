args <- commandArgs(TRUE)
out <- args[1L];base_lib <- args[2L];candidate_lib <- args[3L]
.libPaths(c(base_lib,.libPaths()))
library(fastphylosig);library(ape)
ns <- asNamespace("fastphylosig")
candidate <- dyn.load(file.path(candidate_lib,"fastphylosig","libs","x64","fastphylosig.dll"))
symbols <- grep("^_fastphylosig_",ls(ns,all.names=TRUE),value=TRUE)
original <- mget(symbols,ns,inherits=FALSE)
replacement <- lapply(original,function(x)getNativeSymbolInfo(x$name,candidate,withRegistrationInfo=TRUE))
select <- function(label)for(nm in symbols) {
  unlockBinding(nm,ns);assign(nm,(if(label=="baseline")original else replacement)[[nm]],ns);lockBinding(nm,ns)
}
clean <- function(x) {
  if(is.numeric(x)||is.logical(x)||is.character(x))return(x)
  if(!is.list(x))return(NULL)
  lapply(x[setdiff(names(x),c("call","timing"))],clean)
}
set.seed(4104);tree <- rtree(128)
tree$edge.length <- exp(runif(nrow(tree$edge),-2,1))
y <- setNames(sample(letters[1:3],128,TRUE,prob=c(.6,.3,.1)),tree$tip.label)
workloads <- list(delta=function()fast_delta(tree,y,model="ARD",mcmc_sim=500,thin=5,burn=20,progress=FALSE),
  hetero=function()fast_hetero(tree,y,method="Delta",min_tips=32,max_tips=64,overall=TRUE,
    delta_control=list(model="ARD",mcmc_sim=500L,thin=5L,burn=20L),progress=FALSE))
capture <- function(fn,kind,normal) {
  RNGkind(kind,normal);set.seed(70123)
  warnings <- list()
  value <- withCallingHandlers(fn(),warning=function(w) {
    warnings[[length(warnings)+1L]] <<- list(message=conditionMessage(w),class=class(w))
    invokeRestart("muffleWarning")
  })
  list(value=clean(value),seed=.Random.seed,warnings=warnings,
    next_normal=rnorm(5),next_uniform=runif(5),next_sample=sample.int(100,5),future_seed=.Random.seed)
}
rows <- list()
for(kind in c("Mersenne-Twister","L'Ecuyer-CMRG"))for(normal in c("Inversion","Box-Muller"))
  for(name in names(workloads)) {
    select("baseline");a <- capture(workloads[[name]],kind,normal)
    select("candidate");b <- capture(workloads[[name]],kind,normal)
    equal <- identical(a,b,num.eq=FALSE,single.NA=FALSE)
    rows[[length(rows)+1L]] <- data.frame(workload=name,kind,normal,identical=equal)
    stopifnot(equal)
  }
select("baseline")
rows <- do.call(rbind,rows)
write.csv(rows,file.path(out,"mixed-rng.csv"),row.names=FALSE);print(rows)
