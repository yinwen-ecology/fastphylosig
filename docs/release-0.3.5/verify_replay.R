# Verify the retained production installation against the frozen corrected fit
# trajectory, reconstructed probabilities, complete saved chains and RNG.
args <- commandArgs(TRUE)
out <- args[1L]
if(length(args)>1L).libPaths(c(args[2L],.libPaths()))
label <- if(length(args)>2L)args[3L]else"installed"
library(fastphylosig)
stopifnot(packageVersion("fastphylosig")%in%c("0.3.4.9003","0.3.5"))
saved <- readRDS(file.path(out,"baseline-captured.rds"))
rows <- list()
for(name in names(saved)) {
  expected <- saved[[name]]
  evaluated <- 0L
  for(f in expected$fits) {
    ptr <- fastphylosig:::fast_ace_fit_prepare_cpp(f$edge,f$edge_length,f$tip_state,
                                                f$rate_index,f$n_threads)
    set.seed(901);seed <- .Random.seed
    for(q in f$queries) {
      z <- fastphylosig:::fast_ace_fit_deviance_cpp(ptr,q$par,q$guard)
      stopifnot(identical(z,q$value,num.eq=FALSE,single.NA=FALSE),identical(.Random.seed,seed))
      evaluated <- evaluated+1L
    }
    z <- fastphylosig:::fast_ace_fit_liks_cpp(ptr,f$par,FALSE,f$guard)
    stopifnot(identical(z,f$likelihood,num.eq=FALSE,single.NA=FALSE),identical(.Random.seed,seed))
    fastphylosig:::fast_ace_fit_release_cpp(ptr)
  }
  for(s in expected$samplers) {
    assign(".Random.seed",s$seed,.GlobalEnv)
    z <- do.call(fastphylosig:::delta_mcmc_cpp,c(list(s$probabilities),s$control))
    stopifnot(identical(z,s$value,num.eq=FALSE,single.NA=FALSE),identical(.Random.seed,s$seed_end))
  }
  rows[[name]] <- data.frame(workload=name,objective_values=evaluated,
    posteriors=length(expected$fits),samplers=length(expected$samplers),
    results_identical=TRUE,rng_identical=TRUE)
}
rows <- do.call(rbind,rows)
write.csv(rows,file.path(out,paste0(label,"-replay.csv")),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,paste0(label,"-session.txt")))
print(rows)
cat("Retained installed library:",find.package("fastphylosig"),"\n")
