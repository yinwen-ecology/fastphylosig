# Match a phylogeny and trait table ------------------------------------------

match_phylo_data <- function(tree, X, prune = TRUE, verbose = TRUE) {
  if (!inherits(tree, "phylo")) {
    stop("tree should be an object of class \"phylo\".", call. = FALSE)
  }
  X0 <- .as_named_match_table(X, tree)

  tree_tips <- tree$tip.label
  data_names <- rownames(X0)
  matched <- tree_tips[tree_tips %in% data_names]
  removed_tree <- setdiff(tree_tips, data_names)
  removed_data <- setdiff(data_names, tree_tips)
  if (length(matched) < 2) {
    stop("Fewer than two species are shared by tree and X.", call. = FALSE)
  }

  matched_tree <- if (isTRUE(prune) && length(removed_tree)) {
    ape::drop.tip(tree, removed_tree)
  } else {
    tree
  }
  matched_X <- X0[matched, , drop = FALSE]

  report <- data.frame(
    original_tree_tips = length(tree_tips),
    input_rows = nrow(X0),
    matched_species = length(matched),
    removed_tree_tips = length(removed_tree),
    removed_data_rows = length(removed_data),
    stringsAsFactors = FALSE
  )

  if (isTRUE(verbose)) {
    message(sprintf(
      "Matched %d species; removed %d tree tips and %d data rows.",
      length(matched), length(removed_tree), length(removed_data)
    ))
  }

  out <- list(
    tree = matched_tree,
    X = matched_X,
    report = report,
    matched_species = matched,
    tree_tips_removed = removed_tree,
    data_rows_removed = removed_data
  )
  class(out) <- "fastphylosig_match"
  out
}

.as_named_match_table <- function(X, tree) {
  if (is.data.frame(X)) {
    out <- X
  } else if (is.null(dim(X))) {
    out <- data.frame(x = I(X), stringsAsFactors = FALSE)
    rownames(out) <- names(X)
  } else {
    out <- as.data.frame(X, stringsAsFactors = FALSE)
  }

  default_rownames <- identical(rownames(out), as.character(seq_len(nrow(out))))
  if ((is.null(rownames(out)) || default_rownames) &&
      nrow(out) == ape::Ntip(tree)) {
    message("X has no species row names; assuming X is in tree$tip.label order")
    rownames(out) <- tree$tip.label
  }
  if (is.null(rownames(out))) {
    stop("X must have species names as names/rownames.", call. = FALSE)
  }
  if (anyDuplicated(rownames(out))) {
    stop("X species names must be unique.", call. = FALSE)
  }
  if (is.null(colnames(out))) {
    colnames(out) <- paste0("trait_", seq_len(ncol(out)))
  }
  out
}
