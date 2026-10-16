# Covariance oracles are independent of the native lambda evaluator.
.hetero_scale_crown <- function(shape = "balanced", height = 1) {
  tree <- if (shape == "irregular") ape::read.tree(text = paste0(
    "(((a:.5,b:.8):.4,(c:.4,d:.7):.6):.3,",
    "((e:.2,f:.6):.5,(g:.9,h:.4):.3):.8);")) else ape::stree(8, shape)
  if (is.null(tree$edge.length)) tree$edge.length <- rep(1, nrow(tree$edge))
  depth <- ape::node.depth.edgelength(tree)
  height0 <- max(depth) + .2
  terminal <- which(tree$edge[, 2] <= length(tree$tip.label))
  tree$edge.length[terminal] <- tree$edge.length[terminal] +
    height0 - depth[tree$edge[terminal, 2]]
  tree$edge.length <- tree$edge.length * height / height0
  tree
}
.hetero_scale_covariance <- function(C, lambda) lambda * C + (1-lambda) * diag(diag(C))
.hetero_scale_contrasts <- function(n) cbind(diag(n-1L), -1)

test_that("ultrametric crown contrasts have the derived local lambda scale", {
  for (shape in c("balanced", "left", "irregular")) {
    tree <- .hetero_scale_crown(shape)
    C <- ape::vcv.phylo(tree)
    L <- .hetero_scale_contrasts(nrow(C))
    for (s in c(0, .5, 3, 20)) for (lambda in c(0, .2, .6, .9, 1)) {
      a <- 1 + (1-lambda)*s
      local <- lambda/a
      global <- .hetero_scale_covariance(C+s, lambda)
      expected <- a * .hetero_scale_covariance(C, local)
      expect_equal(L %*% global %*% t(L), L %*% expected %*% t(L), tolerance=2e-12)
    }
  }
  expect_equal(.6/(1+.4*3), 3/11)
  expect_equal(.6*3/(3+.4), 9/17)
})

test_that("native crown extraction and NA root contraction preserve the scale relation", {
  a <- .hetero_scale_crown("balanced", 1)
  b <- .hetero_scale_crown("left", 3)
  a$tip.label <- paste0("A", seq_len(8)); b$tip.label <- paste0("B", seq_len(8))
  tree <- ape::read.tree(text=paste0("(",sub(";$","",ape::write.tree(a,digits=17)),":3,",
    sub(";$","",ape::write.tree(b,digits=17)),":1);"))
  core <- .heterogeneity_compile(tree,"lambda")
  plan <- heterogeneity_plan_cpp(core,tree$tip.label,NULL,3L,8L,100L)
  global <- .hetero_scale_covariance(ape::vcv.phylo(tree),.6)
  s_before <- numeric(2)
  for (i in 1:2) {
    crown <- plan$trees[[i]]
    # In the first crown, retain one daughter clade only, moving its root.
    sets <- list(crown$tip.label, crown$tip.label[1:4])
    for (labels in sets) {
      retained <- ape::keep.tip(crown,labels)
      C <- ape::vcv.phylo(retained)
      labels <- rownames(C)
      h <- diag(C)[1]; s <- 4-h
      if(length(labels)==8L)s_before[i] <- s
      L <- .hetero_scale_contrasts(length(labels))
      local <- .6*h/(h+.4*s); variance <- 1+.4*s/h
      expect_equal(L %*% global[labels,labels] %*% t(L),
        variance * L %*% .hetero_scale_covariance(C,local) %*% t(L),tolerance=2e-12)
      if(i==1 && length(labels)==4L) expect_gt(s,s_before[i])
    }
    compiled <- compile_tree_cpp(crown$edge,crown$edge.length,length(crown$tip.label))
    expect_equal(compiled$root_distance[seq_len(8)], rep(if(i==1)1 else 3,8))
  }
  expect_equal(s_before,c(3,1))
})

test_that("non-ultrametric restrictions need not have a scalar local lambda", {
  tree <- ape::stree(8,"balanced")
  tree$edge.length <- seq(.2,1.7,length.out=nrow(tree$edge))
  C <- ape::vcv.phylo(tree); L <- .hetero_scale_contrasts(8)
  A <- L %*% C %*% t(L)
  B <- L %*% diag(diag(C)) %*% t(L)
  target <- L %*% .hetero_scale_covariance(C+2,.6) %*% t(L)
  design <- cbind(as.vector(A),as.vector(B))
  residual <- qr.resid(qr(design),as.vector(target))
  expect_gt(max(abs(residual)),.01)
  # The BM boundary is unchanged after removing the ancestral constant.
  bm <- L %*% .hetero_scale_covariance(C+2,1) %*% t(L)
  expect_equal(bm,A)
})

test_that("native local lambda remains invariant to nuisance shifts and scale", {
  set.seed(9070)
  tree <- .hetero_scale_crown("irregular")
  x <- setNames(sin(seq_len(8)*1.7)+cos(seq_len(8)*.4),tree$tip.label)
  seed <- .Random.seed
  a <- fast_lambda(tree,x,lambda_profile=FALSE,verbose=FALSE,progress=FALSE)
  b <- fast_lambda(tree,3*x+2,lambda_profile=FALSE,verbose=FALSE,progress=FALSE)
  expect_equal(a$estimate,b$estimate,tolerance=1e-7)
  expect_identical(.Random.seed,seed)
})

test_that("scope annotations distinguish description and explicit research tests", {
  z <- .hetero_fixture()
  expect_identical(formals(fast_hetero)$test,FALSE)
  for(method in c("lambda","D","Delta")) {
    x <- if(method=="lambda")z$X else setNames(rep(1:2,32),z$tree$tip.label)
    set.seed(9071)
    fit <- suppressWarnings(fast_hetero(z$tree,x,method=method,nodes=z$nodes,
      min_tips=20,overall=TRUE,d_nsim=19L,
      delta_control=list(model="ER",mcmc_sim=400L,thin=5L,burn=50L),progress=FALSE))
    expect_true(all(is.na(fit$heterogeneity$P_heterogeneity)))
    expect_identical(fit$heterogeneity$nsim_requested,0L)
    seed <- .Random.seed
    text <- capture.output(print(fit))
    expect_true(any(grepl(switch(method,lambda="common local crown lambda",
      D="experimental common ER",Delta="local/overall signal description"),text,fixed=TRUE)))
    if(method=="Delta") {
      expect_true(any(grepl("Descriptive Delta comparison",.heterogeneity_summary_text(fit,"x"))))
      # Only replace metadata to exercise labelling; do not re-run a long test.
      fit$heterogeneity$nsim_requested <- 19L
      expect_true(any(grepl("not recommended",.heterogeneity_summary_text(fit,"x"))))
    }
    expect_identical(.Random.seed,seed)
  }
})
