args <- commandArgs(TRUE)
options(fastphylosig.ace_threads=1)
.libPaths(c(normalizePath(args[1]),.libPaths()))
library(fastphylosig)
set.seed(902100)
tree <- ape::stree(64,"balanced")
tree$edge.length <- runif(126,.25,1.75)
nodes <- tree$edge[tree$edge[,1]==65,2]
continuous <- cbind(x=sin(1:64),missing=cos((1:64)*.3))
rownames(continuous) <- tree$tip.label
continuous[c(2,9,40),2] <- NA
binary <- setNames(rep(c("no","yes"),32),tree$tip.label)
multi <- setNames(rep(c("A","B","C"),length.out=64),tree$tip.label)
records <- list()
clean <- function(x) {
  if(is.list(x)) {
    x[intersect(names(x),c("timing","timing_seconds"))] <- NULL
    for(i in seq_along(x))x[i] <- list(clean(x[[i]]))
  }
  attr(x,"timing") <- NULL
  attr(x,"runtime") <- NULL
  x
}
run <- function(key,method,data,test,extra=list()) {
  set.seed(902101)
  fit <- suppressWarnings(do.call(fast_hetero,c(list(tree=tree,data=data,method=method,
    nodes=nodes,min_tips=20L,test=test,nsim=9L,d_nsim=39L,
    overall_test=TRUE,overall_nsim=3L,clade_test=TRUE,clade_nsim=3L,
    keep_null=TRUE,progress=FALSE),extra)))
  seed <- .Random.seed
  fit$call <- NULL
  # Only elapsed-time metadata are excluded; no numerical/diagnostic field is dropped.
  fit <- clean(fit)
  records[[key]] <<- list(fit=fit,seed=seed,uniform=runif(11),normal=rnorm(11))
}
for(method in c("K","lambda","D")) for(test in c(FALSE,TRUE)) {
  run(paste(method,test,sep="_"),method,if(method=="D")binary else continuous,test)
}
for(model in c("ER","ARD")) {
  control <- list(model=model,mcmc_sim=1200L,thin=5L,burn=100L)
  run(paste0("Delta_",model,"_description"),"Delta",multi,FALSE,list(delta_control=control))
  run(paste0("Delta_",model,"_plugin"),"Delta",multi,TRUE,list(delta_control=control))
  control <- c(control,list(reference_test="nuisance",state_space=c("A","B","C"),
    nuisance_control=list(max_boxes=4L)))
  run(paste0("Delta_",model,"_nuisance"),"Delta",multi,TRUE,list(delta_control=control))
}
saveRDS(records,args[2],version=2)
cat("Saved",length(records),"complete non-timing contracts; package",as.character(packageVersion("fastphylosig")),"\n")
