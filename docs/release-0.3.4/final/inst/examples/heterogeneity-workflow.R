# Recommended descriptive Delta/D workflow. B controls ordinary whole-tree
# signal randomizations only; it does not enable Delta/D heterogeneity tests.
library(fastphylosig)
set.seed(1142)
tree <- ape::stree(256, 'balanced')
tree$edge.length <- rep(1, nrow(tree$edge))
multi <- ape::rTraitDisc(tree, model='ER', k=3, rate=.15,
  root.value=sample.int(3L,1L))
binary <- ape::rTraitDisc(tree, model='ER', k=2, rate=.15,
  root.value=sample.int(2L,1L))
B <- 199L
delta <- fast_hetero(tree, multi, method='Delta', min_tips=50L,
  max_tips=80L, test=FALSE, overall=TRUE, overall_test=TRUE, overall_nsim=B,
  delta_control=list(model='ER'), progress=FALSE)
d <- fast_hetero(tree, binary, method='D', nodes=delta$selection$nodes,
  min_tips=51L, test=FALSE, d_nsim=999L,
  overall_test=TRUE, progress=FALSE)
stopifnot(!delta$selection$overlap,delta$selection$coverage==1,
  length(delta$selection$nodes)==4L,
  all(delta$clades$n_species==64L),
  all(is.na(delta$heterogeneity$P_heterogeneity)),
  all(is.na(d$heterogeneity$P_heterogeneity)),
  all(delta$heterogeneity$nsim_requested==0L),
  identical(d$selection$nodes,delta$selection$nodes),
  all(c('ESS_Delta','split_Rhat_Delta','MCSE_Delta') %in% names(delta$clades)))
output <- Sys.getenv('FASTPHYLOSIG_EXAMPLE_OUT',tempdir())
dir.create(output,recursive=TRUE,showWarnings=FALSE)
for(method in c('delta','d')) {
  fit <- get(method); seed <- .Random.seed
  grDevices::png(file.path(output,paste0(method,'-heterogeneity-tree.png')),
    width=1600,height=1000,res=150)
  info <- plot(fit,tree=tree,show.tip.label=FALSE,node_labels=TRUE)
  grDevices::dev.off()
  stopifnot(identical(.Random.seed,seed),
    any(grepl('Overall',info$summary_text)),
    any(grepl('heterogeneity|signal-deviation',info$summary_text)))
  if(method=='delta')stopifnot(any(grepl('Descriptive Delta comparison',info$summary_text)))
  # State composition lists remain in the RDS; CSV is for scalar display fields.
  scalar <- fit$clades[!vapply(fit$clades,is.list,logical(1))]
  write.csv(scalar,file.path(output,paste0(method,'-clades.csv')),row.names=FALSE)
  saveRDS(fit,file.path(output,paste0(method,'-fit.rds')))
  print(fit$overall);print(fit$heterogeneity)
}
cat('Example artifacts:',output,'\n')
