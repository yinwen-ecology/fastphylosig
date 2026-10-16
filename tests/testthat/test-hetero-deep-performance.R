test_that("interval jumps preserve maximal disjoint selection and exact coverage", {
  set.seed(672)
  for (type in c("random", "balanced", "left", "star")) {
    tree <- if (type == "random") ape::rtree(256) else ape::stree(256, type)
    tree$edge.length <- rep(.2, nrow(tree$edge))
    # Source tip numbers need not follow traversal order.
    map <- c(sample.int(256), 256L + rev(seq_len(tree$Nnode)))
    tree$edge[] <- map[tree$edge]
    core <- compile_tree_cpp(tree$edge, tree$edge.length, 256L)
    ix <- .hetero_index_reference(core)
    for (bounds in list(c(3L,40L), c(20L,80L), c(50L,256L))) {
      selected <- integer(); last <- 0L
      for (node in core$preorder) if (node > 256L && ix$start[node] > last &&
          ix$n_tip[node] >= bounds[1] && ix$n_tip[node] <= bounds[2]) {
        selected <- c(selected,node); last <- ix$end[node]
      }
      if (!length(selected)) {
        expect_error(heterogeneity_plan_cpp(core, tree$tip.label, NULL, bounds[1],bounds[2],1e6), "No clades")
        next
      }
      plan <- heterogeneity_plan_cpp(core, tree$tip.label, NULL,bounds[1],bounds[2],1e6)
      expect_identical(plan$nodes,selected)
      expect_false(plan$overlap)
      expect_equal(sum(plan$covered),sum(plan$n_tip))
    }
  }
})

test_that("cached local K tests match public kernels, exclusions and RNG", {
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind,as.list(old_kind)),add = TRUE)
  set.seed(621)
  tree <- ape::rtree(180)
  X <- matrix(rnorm(180*6),180,6,dimnames = list(tree$tip.label,paste0("x",1:6)))
  X[seq(3,180,7),2:3] <- NA
  X[seq(4,180,5),4] <- NA
  X[,5] <- 1
  X[1:179,6] <- NA
  core <- .heterogeneity_compile(tree,"K")
  plan <- heterogeneity_plan_cpp(core,tree$tip.label,NULL,10L,60L,1e6)
  scan <- heterogeneity_scan_cpp(plan,X,"K",10L,TRUE,NULL,TRUE)
  eligible <- scan$status == "ok"
  for (kind in c("Mersenne-Twister","L'Ecuyer-CMRG")) {
    RNGkind(kind)
    set.seed(305)
    ref <- rep(NA_real_,length(eligible)); n <- integer(length(eligible))
    for (i in seq_along(plan$nodes)) {
      ids <- (i-1L)*ncol(X)+seq_len(ncol(X))
      j <- which(eligible[ids]); if (!length(j)) next
      result <- fast_k(plan$trees[[i]],X[plan$tips[[i]],j,drop=FALSE],test=TRUE,
        nsim=29,simulation_chunk=3,verbose=FALSE,progress=FALSE)
      ref[ids[j]] <- result$P_fast; n[ids[j]] <- result$nsim_successful
    }
    seed <- .Random.seed
    for (chunk in c(1L,3L,32L)) {
      set.seed(305)
      got <- heterogeneity_k_signal_cpp(scan$cores,scan$data,eligible,29L,chunk)
      expect_identical(got$P,ref)
      expect_identical(got$nsim_successful,n)
      expect_identical(.Random.seed,seed)
    }
  }
})

test_that("native ACE preflight retains binary/root diagnostics without RNG", {
  tree <- ape::stree(64,"balanced"); tree$edge.length <- rep(1,nrow(tree$edge))
  set.seed(54); seed <- .Random.seed
  expect_identical(.fast_ace_topology_guard(tree),63L)
  expect_identical(.Random.seed,seed)
  root <- 65L; map <- seq_len(127); map[c(root,66)] <- c(66L,root)
  malformed <- tree; malformed$edge[] <- map[malformed$edge]
  expect_error(.fast_ace_topology_guard(malformed),"structural root")
  star <- ape::stree(3,"star"); star$edge.length <- rep(1,3); star$Nnode <- 2L
  expect_error(.fast_ace_topology_guard(star),"binary tree|children")
})
