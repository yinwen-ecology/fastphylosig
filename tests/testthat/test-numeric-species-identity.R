test_that('explicit numeric species identities survive both preparation layers', {
  tree <- ape::read.tree(text = '((3:1,1:1):1,2:2);')
  values <- setNames(c(10, 20, 30), as.character(1:3))
  inputs <- list(values, matrix(values, 3, dimnames = list(names(values), 'x')),
                 data.frame(x = unname(values), row.names = names(values)))
  for (x in inputs) {
    table <- fastphylosig:::.as_named_trait_table(x, tree, verbose = FALSE)
    mat <- fastphylosig:::.as_trait_matrix(x, tree, verbose = FALSE)
    expect_identical(rownames(table), names(values))
    expect_identical(rownames(mat), names(values))
    expect_equal(as.numeric(table[tree$tip.label, 1]), c(30, 10, 20))
    expect_equal(as.numeric(mat[tree$tip.label, 1]), c(30, 10, 20))
    # Downstream conversion cannot reinterpret explicitly assigned rownames.
    expect_identical(rownames(fastphylosig:::.as_trait_matrix(table, tree,
      verbose = FALSE)), names(values))
    subset <- if (is.null(dim(x))) x[1:2] else x[1:2, , drop = FALSE]
    expect_identical(rownames(fastphylosig:::.as_named_trait_table(subset, tree,
      verbose = FALSE)), c('1', '2'))
    expect_identical(rownames(fastphylosig:::.as_trait_matrix(subset, tree,
      verbose = FALSE)), c('1', '2'))
  }
  supplied <- fastphylosig::fast_k(tree, values, test = FALSE,
                                  verbose = FALSE, progress = FALSE)
  positional <- fastphylosig::fast_k(tree, unname(values[tree$tip.label]),
                                    test = FALSE, verbose = FALSE, progress = FALSE)
  expect_equal(as.numeric(supplied), as.numeric(positional))
})

test_that('only absent or automatic species rownames use positional matching', {
  tree <- ape::read.tree(text = '((3:1,1:1):1,2:2);')
  for (x in list(c(10,20,30), matrix(c(10,20,30),3), data.frame(x=c(10,20,30)))) {
    expect_identical(rownames(fastphylosig:::.as_named_trait_table(x,tree,
      verbose=FALSE)), tree$tip.label)
    expect_identical(rownames(fastphylosig:::.as_trait_matrix(x,tree,
      verbose=FALSE)), tree$tip.label)
    subset <- if (is.null(dim(x))) x[1:2] else x[1:2,,drop=FALSE]
    # Data-frame subsetting may create explicit rownames; construct an
    # actually automatic short frame for the missing-identity assertion.
    if (is.data.frame(x)) subset <- data.frame(x=c(10,20))
    expect_error(fastphylosig:::.as_named_trait_table(subset,tree,
      verbose=FALSE), 'species names')
    expect_error(fastphylosig:::.as_trait_matrix(subset,tree,
      verbose=FALSE), 'species names')
  }
  for (nm in list(c('1','1','2'), c('1','','2'), c('1',NA,'2'))) {
    for (x in list(setNames(c(10,20,30),nm),
                  matrix(c(10,20,30),3,dimnames=list(nm,'x')))) {
      expect_error(fastphylosig:::.as_named_trait_table(x,tree,verbose=FALSE))
      expect_error(fastphylosig:::.as_trait_matrix(x,tree,verbose=FALSE))
    }
  }
})
