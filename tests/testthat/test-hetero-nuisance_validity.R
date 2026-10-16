.validity_fixture <- function() {
  tree<-ape::stree(4,type="balanced");tree$edge.length<-c(.29,.51,.83,1.17,1.41,.67)
  x<-t(as.matrix(expand.grid(rep(list(1:3),4L))));storage.mode(x)<-"integer"
  first<-which(colSums(abs(x-c(1L,2L,3L,1L)))==0)
  x<-x[,c(first,setdiff(seq_len(81),first))]
  list(core=compile_tree_cpp(tree$edge,tree$edge.length,4L),x=x,
    index=.fast_ace_rate_index("ARD",3L)$rate_index,
    tail=apply(x,2,function(y)length(unique(y))==3L),g=matrix(-4*log(3),81,2))
}
.validity_pointer <- function(z,model="ER",tail=rep(TRUE,ncol(z$x)))
  heterogeneity_nuisance_prepare_cpp(z$core,z$x,.fast_ace_rate_index(model,3L)$rate_index,
    matrix(-4*log(3),ncol(z$x),2),tail,model=="ER")

test_that("support witnesses certify saturation without replacing the upper bound", {
  z <- .validity_fixture()
  z$x <- matrix(c(1L,2L,1L,2L,1L,3L,2L,3L,3L,1L,3L,2L),4L)
  pointer <- .validity_pointer(z,"ARD",c(TRUE,FALSE,FALSE))
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  bank <- attr(pointer,"bank_diagnostics")
  expect_equal(bank$observed_state_count,2L)
  expect_equal(bank$non_tail_count,2L)
  expect_true(3L %in% bank$support_isolated_states)
  domain <- list(bounds=cbind(numeric(6),rep(Inf,6)),model="ARD",
    index=.fast_ace_rate_index("ARD",3L)$rate_index)
  control <- .heterogeneity_nuisance_control(list(max_boxes=32L))
  set.seed(1901);seed <- .Random.seed
  fit <- .heterogeneity_nuisance_search(pointer,domain,control,1/3)
  expect_identical(.Random.seed,seed)
  expect_equal(c(fit$lower,fit$upper),c(1,1))
  expect_equal(fit$cells,1L)
  expect_identical(fit$P_one_reason,"certified_saturation")
  expect_identical(fit$stop_reason,"certified_P_one")
  expect_true(fit$boundary_covered)
})

test_that("an unresolved envelope is not labelled a certified maximum", {
  z <- .validity_fixture()
  pointer <- heterogeneity_nuisance_prepare_cpp(z$core,z$x,z$index,z$g,z$tail,FALSE)
  on.exit(heterogeneity_nuisance_release_cpp(pointer))
  domain <- list(bounds=cbind(numeric(6),rep(Inf,6)),model="ARD",index=z$index)
  fit <- .heterogeneity_nuisance_search(pointer,domain,
    .heterogeneity_nuisance_control(list(max_boxes=1L)),1/3)
  expect_equal(fit$upper,1)
  expect_lt(fit$lower,1)
  expect_identical(fit$P_one_reason,"unresolved_envelope")
  expect_false(fit$converged)
  expect_gt(fit$gap,.01)
})

test_that("all-tail, duplicate auxiliary marks and gamma ceilings stay distinct", {
  z <- .validity_fixture();z$x <- matrix(rep(c(1L,2L,1L,2L),2),4L)
  domain <- list(bounds=matrix(c(.5,.5),1,2),model="ER")
  for(tail in list(c(TRUE,TRUE),c(TRUE,FALSE))) {
    pointer <- .validity_pointer(z,tail=tail)
    fit <- .heterogeneity_nuisance_search(pointer,domain,
      .heterogeneity_nuisance_control(list(max_boxes=1L,gamma=if(all(tail)).6 else .5)),1/3)
    if(all(tail))expect_identical(fit$P_one_reason,"all_upper_tail") else {
      expect_identical(fit$P_one_reason,"gamma_ceiling")
      expect_equal(fit$raw_upper,.5,tolerance=1e-10)
      expect_equal(fit$bank_diagnostics$observed_duplicate_non_tail_count,1L)
    }
    expect_equal(fit$upper,1)
    heterogeneity_nuisance_release_cpp(pointer)
  }
})
