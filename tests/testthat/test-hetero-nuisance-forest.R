.forest_fixture <- function(missing=FALSE) {
  tree<-ape::stree(4,type="balanced");tree$edge.length<-c(.29,.51,.83,1.17,1.41,.67)
  core<-compile_tree_cpp(tree$edge,tree$edge.length,4L)
  x<-t(as.matrix(expand.grid(rep(list(1:3),4L))));storage.mode(x)<-"integer"
  id<-match("1231",apply(x,2,paste0,collapse=""));x<-x[,c(id,setdiff(1:81,id))]
  if(missing)x[4,]<-NA_integer_
  tail<-apply(x,2,function(z)length(unique(z[!is.na(z)]))==3)
  index<-.fast_ace_rate_index("ARD",3L)$rate_index
  list(tree=tree,core=core,x=x,index=index,tail=tail,
    g=matrix(-sum(!is.na(x[,1]))*log(3),81,2))
}
.forest_oracle <- function(z,rates) {
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

test_that("shared maximum-forest bounds cover finite Q in unbounded boxes", {
  skip_if_not_installed("expm")
  set.seed(172101)
  for(missing in c(FALSE,TRUE)) {
    z<-.forest_fixture(missing)
    pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,z$g,z$tail,FALSE)
    for(hi in list(rep(Inf,6),rep(100,6),c(Inf,30,Inf,50,80,Inf))) {
      lo<-rep(10,6);seed<-.Random.seed
      old<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,.001,TRUE)
      new<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,.001,TRUE,correlation=3L)
      expect_true(new$forest_used);expect_false(new$forest_fallback)
      expect_lte(new$upper,old$upper);expect_gte(new$lower,old$lower)
      expect_true(all(new$log_likelihood[,1]>=old$log_likelihood[,1]))
      expect_true(all(new$log_likelihood[,2]<=old$log_likelihood[,2]))
      expect_identical(.Random.seed,seed)
      for(draw in 1:5) {
        q<-if(draw==1)lo else pmin(hi,lo*10^runif(6,0,3))
        ll<-.forest_oracle(z,q)
        expect_true(all(new$log_likelihood[,1]<=ll+1e-8))
        expect_true(all(new$log_likelihood[,2]>=ll-1e-8))
        weights<-exp(ll-z$g[,1])
        if(weights[1]>=.001) {
          p<-sum(weights[z$tail])/sum(weights)
          expect_lte(new$lower,p+1e-8);expect_gte(new$upper,p-1e-8)
        }
      }
    }
    heterogeneity_nuisance_release_cpp(pointer)
  }
})

test_that("automatic forest bounds tighten fast tails without assuming stationarity", {
  z<-.forest_fixture()
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,z$g,z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  old<-heterogeneity_nuisance_box_cpp(pointer,rep(10,6),rep(Inf,6),.001,TRUE)
  new<-heterogeneity_nuisance_box_cpp(pointer,rep(10,6),rep(Inf,6),.001,TRUE,correlation=1L)
  expect_true(new$forest_used);expect_lt(new$forest_tip_error,.001)
  expect_lte(new$forest_cells,8L);expect_gt(new$forest_cells,0L)
  expect_lt(new$upper,old$upper-.05)
  expect_lt(mean(new$log_likelihood[,2]),mean(old$log_likelihood[,2])-.1)
  zero<-heterogeneity_nuisance_box_cpp(pointer,rep(10,6),rep(Inf,6),0,TRUE,correlation=3L)
  expect_true(zero$forest_used);expect_identical(zero$forest_cells,0L)
  expect_lte(zero$upper,1);expect_gte(zero$lower,0)
})

test_that("nonmixing and unsupported forest requests preserve the old enclosure", {
  z<-.forest_fixture()
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,z$g,z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  for(lo in list(rep(0,6),c(10,0,0,0,0,0))) {
    old<-heterogeneity_nuisance_box_cpp(pointer,lo,rep(Inf,6),.001,TRUE)
    new<-heterogeneity_nuisance_box_cpp(pointer,lo,rep(Inf,6),.001,TRUE,correlation=3L)
    expect_false(new$forest_used);expect_false(new$forest_fallback)
    expect_identical(new[names(old)],old)
  }
  expect_identical(.heterogeneity_nuisance_control(list(correlation="forest"))$correlation,"forest")
})

test_that("large trees use linear-cost likelihood bounds without polynomial cells", {
  tree<-ape::stree(256,type="balanced");tree$edge.length<-rep(2,510)
  core<-compile_tree_cpp(tree$edge,tree$edge.length,256L)
  x<-cbind(rep(1:3,length.out=256),rep(1L,256),rep(3:1,length.out=256))
  storage.mode(x)<-"integer"
  index<-.fast_ace_rate_index("ARD",3L)$rate_index
  pointer<-heterogeneity_nuisance_prepare_cpp(core,x,index,matrix(-256*log(3),3,2),
    c(TRUE,FALSE,FALSE),FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  old<-heterogeneity_nuisance_box_cpp(pointer,rep(10,6),rep(Inf,6),.001,TRUE)
  new<-heterogeneity_nuisance_box_cpp(pointer,rep(10,6),rep(Inf,6),.001,TRUE,correlation=1L)
  expect_true(new$forest_used);expect_false(new$forest_fallback)
  expect_identical(new$forest_cells,0L)
  expect_lt(new$log_likelihood[1,2],-250)
  expect_true(all(new$log_likelihood[,2]<=old$log_likelihood[,2]))
  # At q=10 and t=2, ER's terminal transient is exp(-60); the exact
  # likelihood differs from the iid value by far less than this tolerance.
  expect_true(all(new$log_likelihood[,1]<=-256*log(3)+1e-8))
  expect_true(all(new$log_likelihood[,2]>=-256*log(3)-1e-8))
})

test_that("pairwise coupling covers cycles without a common reset column", {
  skip_if_not_installed("expm")
  z<-.forest_fixture()
  Q<-matrix(0,3,3);Q[cbind(1:3,c(2,3,1))]<-10
  off<-z$index>0;lo<-numeric(6);lo[z$index[off]]<-Q[off]
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,z$g,z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  new<-heterogeneity_nuisance_box_cpp(pointer,lo,rep(Inf,6),.001,TRUE,correlation=1L)
  expect_true(new$forest_used);expect_false(new$forest_fallback)
  expect_gte(new$forest_beta,10-1e-12)
  ll<-.forest_oracle(z,lo);weights<-exp(ll-z$g[,1])
  expect_true(all(new$log_likelihood[,1]<=ll+1e-8))
  expect_true(all(new$log_likelihood[,2]>=ll-1e-8))
  expect_lte(new$lower,sum(weights[z$tail])/sum(weights)+1e-8)
  expect_gte(new$upper,sum(weights[z$tail])/sum(weights)-1e-8)
})
