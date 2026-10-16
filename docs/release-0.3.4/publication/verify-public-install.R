root <- normalizePath(getwd(),winslash='/')
out <- file.path(root,'docs/release-0.3.4/publication')
archive <- file.path(out,'downloads/fastphylosig_0.3.4.tar.gz')
con <- file(archive,'rb');raw <- readBin(con,'raw',n=file.info(archive)$size);close(con)
stopifnot(toupper(digest::digest(raw,algo='sha256',serialize=FALSE)) ==
  '8B1A8940101997492D4F3E2AFE696717BFE63CC9905FE91E225A3534754B997D')
work <- tempfile('fastphylosig-public-034-',tmpdir=Sys.getenv('TEMP'))
dir.create(work);lib <- file.path(work,'library');dir.create(lib)
stopifnot(!grepl('[^ -~]',work),file.copy(archive,file.path(work,basename(archive))))
Sys.setenv(MAKEFLAGS='-j2',OMP_NUM_THREADS='2')
status <- system2(file.path(R.home('bin'),'R.exe'),c('CMD','INSTALL','--preclean','--no-multiarch',
  shQuote(paste0('--library=',lib)),shQuote(file.path(work,basename(archive)))),
  stdout=file.path(out,'public-install.log'),stderr=file.path(out,'public-install-stderr.log'))
stopifnot(status==0L)
.libPaths(c(lib,.libPaths()))
library(fastphylosig,lib.loc=lib)
stopifnot(as.character(packageVersion('fastphylosig'))=='0.3.4')
set.seed(3401)
tree <- ape::stree(64,'balanced');tree$edge.length <- rep(1,nrow(tree$edge))
continuous <- setNames(rnorm(64),tree$tip.label)
binary <- setNames(rep(c('absent','present'),32),tree$tip.label)
multi <- setNames(rep(c('A','B','C'),length.out=64),tree$tip.label)
nodes <- tree$edge[tree$edge[,1]==65,2]
rows <- list();notes <- list()
for(method in c('K','lambda','D','Delta')) {
  warnings <- character()
  fit <- withCallingHandlers(fast_hetero(tree,
    if(method=='D')binary else if(method=='Delta')multi else continuous,
    method=method,nodes=nodes,min_tips=20,overall=TRUE,test=FALSE,
    d_nsim=39,delta_control=list(model='ER'),progress=FALSE),warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')
    })
  stopifnot(nrow(fit$clades)==2L,all(is.finite(fit$clades$estimate)),is.finite(fit$overall$estimate))
  rows[[method]] <- data.frame(method=method,version='0.3.4',local_estimates=2L,
    overall_estimate=fit$overall$estimate,load_and_smoke_passed=TRUE,warnings=length(warnings))
  notes[[method]] <- warnings
}
write.csv(do.call(rbind,rows),file.path(out,'PUBLIC_INSTALL_SMOKE.csv'),row.names=FALSE)
writeLines(jsonlite::toJSON(list(library=lib,archive_sha256=
  '8B1A8940101997492D4F3E2AFE696717BFE63CC9905FE91E225A3534754B997D',
  install_exit=status,warnings=notes),pretty=TRUE,auto_unbox=TRUE),
  file.path(out,'PUBLIC_INSTALL_STATUS.json'))
cat('PASS: public source download installed in a new library; K/lambda/D/Delta loaded and returned finite estimates.\n')
