# Reproducible teaching example for fastphylosig 0.3.4. Run at repository root.
# FASTPHYLOSIG_LIB may point to an isolated installation; optional for users.
lib <- Sys.getenv('FASTPHYLOSIG_LIB')
if(nzchar(lib)) .libPaths(c(lib,.libPaths()))
library(fastphylosig)
stopifnot(as.character(packageVersion('fastphylosig')) == '0.3.4')
options(fastphylosig.ace_threads=1L)
RNGkind('Mersenne-Twister','Inversion','Rejection')
set.seed(20261016)
out <- Sys.getenv('FASTPHYLOSIG_CASE_OUT','docs/assets/fast-hetero-case')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
tree <- ape::stree(512,'balanced')
tree$tip.label <- sprintf('sp%03d',seq_len(512))
tree$edge.length <- rep(1,nrow(tree$edge))
# A tree-only, predetermined cut at depth 3: eight disjoint 64-tip crowns.
depth <- ape::node.depth.edgelength(tree)
nodes <- which(seq_along(depth)>512 & depth==3)
stopifnot(length(nodes)==8L)
x <- setNames(rep(NA_real_,512),tree$tip.label)
binary <- multi <- setNames(rep(NA_character_,512),tree$tip.label)
group <- setNames(rep(NA_character_,512),tree$tip.label)
members <- vector('list',8)
for(i in seq_along(nodes)) {
  crown <- ape::extract.clade(tree,nodes[i])
  tips <- crown$tip.label; members[[i]] <- tips
  group[tips] <- paste0('C',i)
  if(i<=4L) {
    x[tips] <- ape::rTraitCont(crown,model='BM',sigma=1)[tips]
    binary[tips] <- as.character(ape::rTraitDisc(crown,model='ER',k=2,
      states=c('absent','present'),root.value=1,rate=.08)[tips])
    multi[tips] <- as.character(ape::rTraitDisc(crown,model='ER',k=3,
      states=c('A','B','C'),root.value=1,rate=.08)[tips])
  } else {
    x[tips] <- rnorm(64,sd=sqrt(6))
    binary[tips] <- sample(c('absent','present'),64,replace=TRUE)
    multi[tips] <- sample(c('A','B','C'),64,replace=TRUE)
  }
}
# Deliberate teaching cases, chosen before fitting: missing observations and
# one single-state categorical crown. They are not removed from outputs.
x[members[[1]][1:3]] <- NA_real_
binary[members[[8]]] <- 'absent'
multi[members[[8]]] <- 'A'
traits <- data.frame(species=tree$tip.label,crown=unname(group),
  continuous=unname(x),binary=unname(binary),multistate=unname(multi))
write.csv(traits,file.path(out,'traits.csv'),row.names=FALSE,na='NA')
ape::write.tree(tree,file=file.path(out,'tree.nwk'))
write.csv(data.frame(crown=paste0('C',seq_along(nodes)),node=nodes,
  n_tip_total=64,depth=3,height=6),file.path(out,'crowns.csv'),row.names=FALSE)
fits <- list(); warnings <- list(); plot_info <- list()
for(i in seq_along(c('K','lambda','D','Delta'))) {
  method <- c('K','lambda','D','Delta')[i]
  set.seed(20261016+i)
  notes <- character()
  args <- list(tree=tree,data=if(method=='D')binary else if(method=='Delta')multi else x,
    method=method,nodes=nodes,min_tips=50L,overall_test=TRUE,clade_test=FALSE,
    test=method %in% c('K','lambda'),nsim=999L,progress=FALSE)
  if(method=='D') {args$min_tips <- 51L;args$d_nsim <- 1000L}
  if(method=='lambda') args$lambda_reference <- 'REML'
  if(method=='Delta') {
    args$delta_control <- list(model='ER')
    args$overall_nsim <- 199L # ordinary whole-tree signal test, not heterogeneity
  }
  fit <- withCallingHandlers(do.call(fast_hetero,args),warning=function(w) {
    notes <<- c(notes,conditionMessage(w));invokeRestart('muffleWarning')
  })
  fits[[method]] <- fit; warnings[[method]] <- notes
  fields <- fit$clades[!vapply(fit$clades,is.list,logical(1))]
  fields$crown <- paste0('C',match(fields$node,nodes))
  for(nm in intersect(c('state_counts','state_frequencies'),names(fit$clades)))
    fields[[nm]] <- vapply(fit$clades[[nm]],function(v)
      paste(paste(names(v),v,sep='='),collapse=';'),character(1))
  write.csv(fields,file.path(out,paste0(method,'-clades.csv')),row.names=FALSE,na='NA')
  write.csv(fit$overall,file.path(out,paste0(method,'-overall.csv')),row.names=FALSE,na='NA')
  plain <- fit$heterogeneity[!vapply(fit$heterogeneity,is.list,logical(1))]
  write.csv(plain,file.path(out,paste0(method,'-heterogeneity.csv')),row.names=FALSE,na='NA')
  seed <- .Random.seed
  grDevices::png(file.path(out,paste0(method,'-tree.png')),width=1800,height=1400,res=170)
  plot_info[[method]] <- plot(fit,tree=tree,show.tip.label=FALSE,node_labels=TRUE,
    main=paste0(method,': eight prespecified crowns'),cex.main=1)
  grDevices::dev.off()
  stopifnot(identical(.Random.seed,seed))
  if(length(notes)) cat(method,'warnings retained:',paste(notes,collapse=' | '),'\n')
}
stopifnot(all(vapply(fits,function(z)identical(z$selection$nodes,as.integer(nodes)),logical(1))))
stopifnot(all(vapply(fits,function(z)z$selection$coverage==1,logical(1))))
stopifnot(all(is.na(fits$D$heterogeneity$P_heterogeneity)),
          all(is.na(fits$Delta$heterogeneity$P_heterogeneity)))
saveRDS(list(tree=tree,traits=traits,nodes=nodes,fits=fits,warnings=warnings),
  file.path(out,'case-results.rds'),version=2)
writeLines(jsonlite::toJSON(warnings,pretty=TRUE,auto_unbox=TRUE),
  file.path(out,'warnings.json'),useBytes=TRUE)
writeLines(jsonlite::toJSON(lapply(plot_info,function(z)z$summary_text),pretty=TRUE,auto_unbox=TRUE),
  file.path(out,'plot-annotations.json'),useBytes=TRUE)
capture.output(sessionInfo(),file=file.path(out,'session-info.txt'))
cat('Saved all four methods, warnings, diagnostics and native tree plots.\n')
