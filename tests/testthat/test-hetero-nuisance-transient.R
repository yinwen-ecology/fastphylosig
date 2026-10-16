.transient_fixture <- function(missing=FALSE,duplicate=FALSE) {
  tree<-ape::stree(4,type="balanced");tree$edge.length<-c(.29,.51,.83,1.17,1.41,.67)
  core<-compile_tree_cpp(tree$edge,tree$edge.length,4L)
  x<-t(as.matrix(expand.grid(rep(list(1:3),4L))));storage.mode(x)<-"integer"
  first<-match("1231",apply(x,2,paste0,collapse=""));x<-x[,c(first,setdiff(1:81,first))]
  if(missing)x[4,]<-NA_integer_
  tail<-apply(x,2,function(y)length(unique(y[!is.na(y)]))==3)
  if(duplicate) {x<-cbind(x,x[,1]);tail<-c(tail,FALSE)}
  index<-.fast_ace_rate_index("ARD",3L)$rate_index
  list(tree=tree,core=core,x=x,index=index,tail=tail,
    g=matrix(-sum(!is.na(x[,1]))*log(3),ncol(x),2))
}
.transient_oracle <- function(z,rates) {
  Q<-matrix(0,3,3);off<-z$index>0;Q[off]<-rates[z$index[off]];diag(Q)<- -rowSums(Q)
  P<-lapply(z$core$branch_length_by_node,function(t)if(t==0)diag(3) else expm::expm(Q*t))
  vapply(seq_len(ncol(z$x)),function(b) {
    a<-matrix(1,z$core$n_total,3)
    for(i in seq_len(z$core$n_tip))if(!is.na(z$x[i,b])) {a[i,]<-0;a[i,z$x[i,b]]<-1}
    for(node in z$core$postorder)if(z$core$parent[node])
      a[z$core$parent[node],]<-a[z$core$parent[node],]*as.numeric(P[[node]]%*%a[node,])
    log(mean(a[z$core$root,]))
  },numeric(1))
}
.transient_edges <- function() c(lapply(1:3,function(i) {
  j<-setdiff(1:3,i);rbind(c(i,j[1]),c(i,j[2]))
}),lapply(list(c(1,2),c(1,3),c(2,3)),function(p)rbind(p,rev(p))))

test_that("shared transient covers all multiclass families and their slow neighbours", {
  skip_if_not_installed("expm")
  set.seed(161017)
  for(missing in c(FALSE,TRUE)) {
    z<-.transient_fixture(missing,TRUE)
    pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,z$g,z$tail,FALSE)
    for(edges in .transient_edges())for(extra in c(0,.002)) {
      base<-z$index[edges];lo<-numeric(6);hi<-rep(extra,6);hi[base]<-Inf
      seed<-.Random.seed
      old<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,.001,TRUE)
      new<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,.001,TRUE,correlation=4L)
      expect_true(new$transient_used);expect_false(new$transient_fallback)
      expect_lte(new$transient_cells,8L);expect_identical(.Random.seed,seed)
      expect_true(all(new$log_likelihood[,1]>=old$log_likelihood[,1]))
      expect_true(all(new$log_likelihood[,2]<=old$log_likelihood[,2]))
      expect_lte(new$upper,old$upper)
      if(!new$excluded)expect_gte(new$lower,old$lower)
      for(draw in 1:4) {
        q<-runif(6,0,extra);q[base]<-if(draw==1)c(0,0) else if(draw==2)c(.2,0) else 10^runif(2,-3,3)
        ll<-.transient_oracle(z,q)
        expect_true(all(new$log_likelihood[,1]<=ll+1e-8))
        expect_true(all(new$log_likelihood[,2]>=ll-1e-8))
        weights<-exp(ll-z$g[,1])
        if(weights[1]>=.001) {
          p<-sum(weights[z$tail])/sum(weights)
          expect_false(new$excluded);expect_lte(new$lower,p+1e-8);expect_gte(new$upper,p-1e-8)
        }
      }
    }
    heterogeneity_nuisance_release_cpp(pointer)
  }
})

test_that("transient exclusion does not remove unconditional likelihoods", {
  z<-.transient_fixture()
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,z$g,z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  hi<-numeric(6);hi[z$index[cbind(c(1,2),c(2,1))]]<-Inf
  new<-heterogeneity_nuisance_box_cpp(pointer,numeric(6),hi,.001,TRUE,correlation=4L)
  expect_true(new$transient_used);expect_true(new$excluded);expect_identical(new$upper,0)
  expect_gt(max(new$log_likelihood[,2]),-2)
  zero<-heterogeneity_nuisance_box_cpp(pointer,numeric(6),hi,0,TRUE,correlation=4L)
  expect_identical(zero$log_likelihood,new$log_likelihood)
})

test_that("transient bounds keep the complete domain and fixed lower witnesses", {
  z<-.transient_fixture()
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,z$g,z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  lo<-numeric(6);hi<-rep(Inf,6)
  old<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,.001,TRUE)
  new<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,.001,TRUE,correlation=4L)
  expect_false(new$transient_used);expect_false(new$transient_fallback)
  expect_identical(new$log_likelihood,old$log_likelihood)
  expect_true(new$mass_used);expect_lt(new$upper,old$upper)
  control<-.heterogeneity_nuisance_control(list(max_boxes=1L,correlation="transient"))
  domain<-list(bounds=cbind(lo,hi),model="ARD")
  set.seed(161018);seed<-.Random.seed
  fit<-.heterogeneity_nuisance_search(pointer,domain,control,1)
  expect_gt(fit$witness_evaluations,0L);expect_identical(.Random.seed,seed)
  expect_true(fit$boundary_covered);expect_gt(fit$upper,.99)
  point<-heterogeneity_nuisance_box_cpp(pointer,fit$best_rates,fit$best_rates,control$gamma)
  expect_true(point$included)
  expect_lte(fit$lower,point$lower+control$gamma+1e-12)
})

test_that("automatic transient route skips large banks and preserves legacy requests", {
  z<-.transient_fixture();hi<-numeric(6);hi[z$index[rbind(c(1,2),c(1,3))]]<-Inf
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,z$g,z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  auto<-heterogeneity_nuisance_box_cpp(pointer,numeric(6),hi,.001,TRUE,correlation=1L)
  expect_true(auto$transient_used)
  expect_false(heterogeneity_nuisance_box_cpp(pointer,numeric(6),hi,.001,TRUE,correlation=3L)$forest_used)
  expect_identical(.heterogeneity_nuisance_control(list(correlation="transient"))$correlation,"transient")
  expect_error(.heterogeneity_nuisance_control(list(correlation="bad")),"arg")
  tree<-ape::stree(256,type="balanced");tree$edge.length<-rep(1,510)
  core<-compile_tree_cpp(tree$edge,tree$edge.length,256L)
  x<-matrix(rep(1:3,length.out=256*40),256,40)
  # Distinct columns keep the compressed bank above the cost safeguard.
  for(i in 1:40)x[1:5,i]<-as.integer(intToBits(i)[1:5])+1L
  storage.mode(x)<-"integer"
  large<-heterogeneity_nuisance_prepare_cpp(core,x,z$index,matrix(-256*log(3),40,2),rep(TRUE,40),FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(large),add=TRUE)
  skipped<-heterogeneity_nuisance_box_cpp(large,numeric(6),hi,0,TRUE,correlation=4L)
  expect_false(skipped$transient_used);expect_false(skipped$transient_fallback)
})

test_that("finite shared pruning retains nonnegative zero-boundary likelihoods", {
  z<-.transient_fixture()
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,z$g,z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  finite<-heterogeneity_nuisance_box_cpp(pointer,numeric(6),rep(.002,6),0,TRUE,correlation=1L)
  expect_true(finite$correlation_used)
  expect_false(anyNA(finite$log_likelihood))
  expect_true(all(finite$log_likelihood[,1]<=finite$log_likelihood[,2]))
})

test_that("whole-domain mass constraints tighten the unbounded root box", {
  z<-.transient_fixture()
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,z$g,z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  set.seed(181018);seed<-.Random.seed
  fit<-heterogeneity_nuisance_box_cpp(pointer,numeric(6),rep(Inf,6),.001,TRUE,correlation=1L)
  expect_true(fit$mass_used);expect_false(fit$mass_fallback)
  expect_lt(fit$upper,1-6^-6/2);expect_gte(fit$upper,1-6^-6-1e-12)
  expect_identical(.Random.seed,seed)
  # gamma is added only after maximization; this weak whole-domain improvement
  # need not lower the final, truncated upper bound below one.
  expect_identical(heterogeneity_nuisance_add_gamma_cpp(0,fit$upper,.001)[2],1)
})

test_that("subtree repeats avoid recomputing identical messages", {
  z<-.transient_fixture();hi<-numeric(6);hi[z$index[rbind(c(1,2),c(1,3))]]<-Inf
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,z$g,z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  fit<-heterogeneity_nuisance_box_cpp(pointer,numeric(6),hi,.001,TRUE,correlation=4L)
  expect_true(fit$transient_used)
  expect_equal(fit$transient_messages,8*(81+9+9))
  expect_equal(fit$transient_reuses,8*(3*81-99))
})

test_that("shared likelihoods use independent repeated-phenotype tail weights", {
  z<-.transient_fixture(duplicate=TRUE)
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,z$g,z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  fit<-heterogeneity_nuisance_box_cpp(pointer,numeric(6),rep(Inf,6),.1,TRUE,correlation=4L)
  expect_true(fit$mass_used);expect_lt(fit$upper,1)
  expect_identical(fit$log_likelihood[1,],fit$log_likelihood[ncol(z$x),])
})

test_that("a true whole-domain supremum of one is not artificially reduced", {
  z<-.transient_fixture()
  x<-cbind(c(1L,1L,2L,2L),rep(1L,4),c(3L,3L,1L,1L))
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,x,z$index,matrix(-4*log(3),3,2),c(TRUE,TRUE,FALSE),FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  q<-numeric(6);q[z$index[rbind(c(1,2),c(2,1))]]<-.2
  point<-heterogeneity_nuisance_box_cpp(pointer,q,q,.001,TRUE,correlation=4L)
  full<-heterogeneity_nuisance_box_cpp(pointer,numeric(6),rep(Inf,6),.001,TRUE,correlation=4L)
  expect_true(point$included);expect_gt(point$lower,1-1e-10)
  expect_identical(full$upper,1);expect_false(full$mass_used)
})

test_that("repeated-observation tail ratios are bounded on the complete domain", {
  z<-.transient_fixture();x<-matrix(rep(z$x[,1],4),4,4)
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,x,z$index,matrix(-4*log(3),4,2),c(TRUE,FALSE,FALSE,FALSE),FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  fit<-heterogeneity_nuisance_box_cpp(pointer,numeric(6),rep(Inf,6),.001,TRUE,correlation=4L)
  expect_true(fit$mass_used);expect_lte(fit$upper,.25+1e-12);expect_gte(fit$upper,.25-1e-12)
})

