test_that('long mask keys are indexed by exact identity and survive serialization', {
  cache <- new.env(parent=emptyenv())
  for (n in c(37568L,37569L,50000L,100001L)) {
    masks <- list(seq_len(n),seq_len(n-1L),c(1L,n),integer())
    present <- matrix(FALSE,n,length(masks))
    for (j in seq_along(masks)) present[masks[[j]],j] <- TRUE
    cpp <- fastphylosig:::group_na_masks_cpp(present)
    fallback <- fastphylosig:::.analysis_na_groups_r(present)
    expect_equal(cpp$keep,masks)
    expect_equal(cpp$group_id,1:4)
    for (raw_keys in list(cpp$key,fallback$key)) {
      keys <- vapply(raw_keys,fastphylosig:::.mask_cache_key,character(1),cache=cache)
      expect_lte(max(nchar(keys,type='bytes')),10000)
      expect_length(unique(keys),4)
      for(j in 1:4) assign(keys[j],masks[[j]],cache)
      for(j in 1:4) expect_identical(get(keys[j],cache),masks[[j]])
      restored <- unserialize(serialize(cache,NULL))
      expect_identical(vapply(raw_keys,fastphylosig:::.mask_cache_key,
        character(1),cache=restored),keys)
    }
    if (n==37568L) expect_identical(
      fastphylosig:::.mask_cache_key(cpp$key[1],cache),cpp$key[1])
    if(n==37569L) expect_gt(nchar(cpp$key[1],type='bytes'),10000)
  }
  a <- fastphylosig:::.tree_mask_key(c(1,65),50000)
  b <- fastphylosig:::.tree_mask_key(c(1,65),50001)
  expect_false(identical(fastphylosig:::.mask_cache_key(a,cache),
                         fastphylosig:::.mask_cache_key(b,cache)))
  expect_identical(fastphylosig:::.tree_mask_key(c(65,1),50000),a)
})

test_that('50K preparation, subset caches and checked subsets accept long masks', {
  set.seed(251)
  tree <- ape::rtree(50000L)
  ctx <- fastphylosig::prepare_tree(tree,cache_budget=0)
  fingerprint <- ctx$fingerprint
  expect_equal(fastphylosig::cache_info(ctx)$n_structural_entries,1)
  full <- fastphylosig:::.prepared_tree_subset(ctx,need_matrix=FALSE)
  expect_equal(full$n_tip,50000)
  keep <- c(1:24999,25001:50000)
  subset <- fastphylosig:::.prepared_tree_subset(ctx,keep,need_matrix=FALSE)
  again <- fastphylosig:::.prepared_tree_subset(ctx,rev(keep),need_matrix=FALSE)
  expect_identical(subset$subtree_key,again$subtree_key)
  expect_identical(subset$retained_parent_indices,sort(as.integer(keep)))
  expect_identical(subset$tree$tip.label,ctx$tree$tip.label[keep])
  expect_equal(fastphylosig::cache_info(ctx)$n_structural_entries,2)
  restored <- unserialize(serialize(ctx,NULL,version=3L))
  expect_identical(fastphylosig:::.prepared_tree_subset(restored,keep,
    need_matrix=FALSE)$subtree_key,subset$subtree_key)
  expect_identical(restored$fingerprint,fingerprint)
  X <- matrix(seq_len(100000),50000,2,dimnames=list(ctx$tree$tip.label,c('a','b')))
  X[25000,2] <- NA
  prepared <- fastphylosig:::.prepare_analysis(ctx,X,'K','continuous',verbose=FALSE)
  expect_equal(lengths(prepared$na_patterns$keep),c(50000L,49999L))
  expect_true(all(vapply(prepared$groups,function(g)g$kernel_ready,logical(1))))
  testthat::local_mocked_bindings(.analysis_na_groups=fastphylosig:::.analysis_na_groups_r,
                                 .package='fastphylosig')
  fallback <- fastphylosig:::.prepare_analysis(ctx,X,'K','continuous',verbose=FALSE)
  expect_equal(lapply(fallback$na_patterns$keep,unname),prepared$na_patterns$keep)
  expect_equal(vapply(fallback$groups,function(g)g$group$n_tip,integer(1)),c(50000L,49999L))
  fit <- fastphylosig::fast_k(ctx,X[,1],test=FALSE,verbose=FALSE,progress=FALSE)
  expect_true(is.finite(as.numeric(fit)))
  expect_identical(ctx$fingerprint,fingerprint)
  expect_equal(fastphylosig::cache_info(ctx)$bytes_used,0)
})

test_that('numerical payload caches share the exact long-mask index', {
  ctx <- fastphylosig::prepare_tree(ape::read.tree(text='((a:1,b:1):1,c:2);'))
  group <- fastphylosig:::.prepared_tree_subset(ctx,need_matrix=FALSE)
  raw <- fastphylosig:::.tree_mask_key(seq_len(50000),50000)
  key <- fastphylosig:::.mask_cache_key(raw,ctx$structural_cache)
  # A small tree payload exercises the numeric cache without a 50K VCV.
  assign(key,group,ctx$structural_cache)
  first <- fastphylosig:::.numerical_payload(ctx,key,group,need_lambda=TRUE)
  again <- fastphylosig:::.numerical_payload(ctx,key,group,need_lambda=TRUE)
  expect_equal(first$C,ape::vcv.phylo(ctx$tree))
  expect_equal(first$C,again$C)
  expect_true(exists(key,ctx$numerical_cache,inherits=FALSE))
  info <- fastphylosig::cache_info(ctx)
  expect_true(info$entries$has_vcv[info$entries$entry==key])
  restored <- unserialize(serialize(ctx,NULL,version=3L))
  expect_identical(fastphylosig:::.mask_cache_key(raw,restored$structural_cache),key)
  expect_equal(fastphylosig:::.numerical_payload(restored,key,group)$C,first$C)
})
