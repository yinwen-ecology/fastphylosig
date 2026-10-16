.correlated_fixture <- function() {
  tree<-ape::stree(4,type="balanced");tree$edge.length<-c(.29,.51,.83,1.17,1.41,.67)
  core<-compile_tree_cpp(tree$edge,tree$edge.length,4L)
  x<-t(as.matrix(expand.grid(rep(list(1:3),4L))));storage.mode(x)<-"integer"
  index<-.fast_ace_rate_index("ARD",3L)$rate_index
  list(tree=tree,core=core,x=x,index=index,tail=seq_len(81)%%3==1)
}
.correlated_oracle <- function(z,rates) {
  Q<-matrix(0,3,3);off<-z$index>0;Q[off]<-rates[z$index[off]];diag(Q)<- -rowSums(Q)
  P<-lapply(z$core$branch_length_by_node,function(t)if(t==0)diag(3) else expm::expm(Q*t))
  vapply(seq_len(ncol(z$x)),function(b) {
    a<-matrix(1,z$core$n_total,3)
    for(i in 1:4)if(!is.na(z$x[i,b])) {a[i,]<-0;a[i,z$x[i,b]]<-1}
    for(node in z$core$postorder)if(z$core$parent[node])
      a[z$core$parent[node],]<-a[z$core$parent[node],]*as.numeric(P[[node]]%*%a[node,])
    log(mean(a[z$core$root,]))
  },numeric(1))
}
test_that("shared-Q Taylor models enclose independent pruning and weighted ranks", {
  skip_if_not_installed("expm")
  z<-.correlated_fixture();g<-matrix(-4*log(3),81,2)
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,g,z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  set.seed(17101);seed<-.Random.seed
  center<-c(.03,.08,.12,.27,.4,.7)
  for(radius in c(.0005,.005,.02)) {
    lo<-center-radius;hi<-center+radius
    old<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,0,TRUE)
    new<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,0,TRUE,correlation=2L)
    expect_true(new$correlation_used);expect_false(new$correlation_fallback)
    expect_true(all(new$log_likelihood[,1]>=old$log_likelihood[,1]))
    expect_true(all(new$log_likelihood[,2]<=old$log_likelihood[,2]))
    expect_lte(new$upper,old$upper);expect_gte(new$lower,old$lower)
    for(multiplier in c(-1,-.5,0,.5,1)) {
      expected<-.correlated_oracle(z,center+radius*multiplier*c(1,-1,1,-1,1,-1))
      expect_true(all(new$log_likelihood[,1]<=expected+1e-10))
      expect_true(all(new$log_likelihood[,2]>=expected-1e-10))
      weights<-exp(expected-g[,1]);p<-sum(weights[z$tail])/sum(weights)
      expect_lte(new$lower,p+1e-10);expect_gte(new$upper,p-1e-10)
    }
  }
  expect_identical(.Random.seed,seed)
})
test_that("unbounded shared-Q requests retain the complete boundary envelope", {
  z<-.correlated_fixture()
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,
    matrix(-4*log(3),81,2),z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  lo<-rep(.1,6);hi<-rep(.2,6);hi[c(1,4)]<-Inf
  old<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,0,TRUE)
  new<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,0,TRUE,correlation=2L)
  expect_false(new$correlation_used);expect_false(new$correlation_fallback)
  expect_identical(new[names(old)],old)
  expect_error(.heterogeneity_nuisance_control(list(correlation="gradient")),"arg")
})

test_that("a constant subnormal rate remains an exact shared-Q coordinate", {
  z<-.correlated_fixture()
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,
    matrix(-4*log(3),81,2),z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  lo<-rep(.199,6);hi<-rep(.201,6)
  lo[1]<-hi[1]<-.Machine$double.xmin*.Machine$double.eps
  old<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,0,TRUE)
  new<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,0,TRUE,correlation=2L)
  expect_true(new$correlation_used);expect_false(new$correlation_fallback)
  expect_lte(new$upper,old$upper);expect_gte(new$lower,old$lower)
  expect_true(all(new$log_likelihood[,1]>=old$log_likelihood[,1]))
  expect_true(all(new$log_likelihood[,2]<=old$log_likelihood[,2]))
})

test_that("automatic shared-Q work is restricted by rate-box geometry", {
  z<-.correlated_fixture()
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,
    matrix(-4*log(3),81,2),z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  narrow<-heterogeneity_nuisance_box_cpp(pointer,rep(.199,6),rep(.201,6),
    0,correlation=1L)
  wide<-heterogeneity_nuisance_box_cpp(pointer,rep(.1,6),rep(.3,6),
    0,correlation=1L)
  original<-heterogeneity_nuisance_box_cpp(pointer,rep(.1,6),rep(.3,6),0)
  forced<-heterogeneity_nuisance_box_cpp(pointer,rep(.1,6),rep(.3,6),
    0,correlation=2L)
  expect_true(narrow$correlation_used)
  expect_false(wide$correlation_used);expect_false(wide$correlation_fallback)
  expect_identical(wide[names(original)],original)
  expect_true(forced$correlation_used)
  expect_lte(forced$upper,original$upper);expect_gte(forced$lower,original$lower)
})
test_that("confidence-set dual certificates retain every included shared Q", {
  skip_if_not_installed("expm")
  z<-.correlated_fixture()
  id<-match("1231",apply(z$x,2,paste0,collapse=""))
  z$x<-z$x[,c(id,setdiff(1:81,id))]
  z$tail<-seq_len(81)%%3==1
  center<-rep(.2,6);g<-matrix(-4*log(3),81,2)
  # Fixed synthetic fixture puts the confidence boundary through its center.
  gamma<-exp(.correlated_oracle(z,center)[1]-g[1,1])
  expect_gt(gamma,0);expect_lt(gamma,1)
  pointer<-heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,g,z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  lo<-center-.005;hi<-center+.005
  old<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,gamma,TRUE)
  new<-heterogeneity_nuisance_box_cpp(pointer,lo,hi,gamma,TRUE,correlation=2L)
  expect_true(new$correlation_used);expect_false(new$correlation_fallback)
  expect_false(new$excluded);expect_false(new$included)
  expect_lte(new$upper,old$upper);expect_gte(new$lower,old$lower)
  included<-0L
  for(offset in seq(-.005,.005,length.out=9)) {
    expected<-.correlated_oracle(z,center+offset)
    weights<-exp(expected-g[,1])
    if(weights[1]>=gamma) {
      included<-included+1L;p<-sum(weights[z$tail])/sum(weights)
      expect_lte(new$lower,p+1e-10);expect_gte(new$upper,p-1e-10)
    }
  }
  expect_gt(included,0L)
})
