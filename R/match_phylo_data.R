# Match a phylogeny and trait table ------------------------------------------

#' Match a phylogeny and trait table
#'
#' Matches a phylogeny and trait table without silently discarding information.
#' The returned object contains the matched tree, the trait table reordered in
#' tree-tip order, and a report listing matched species, removed tree tips, and
#' removed data rows.
#'
#' @param tree An object of class \code{"phylo"} or a prepared tree context.
#' @param data Named vector, matrix, or data.frame. Species names must be
#'   supplied as names or row names unless rows are already in tree-tip order.
#' @param X Deprecated compatibility alias for `data`.
#' @param prune Whether to drop tree tips that are absent from \code{X}.
#' @param verbose Whether to print one short message when rows/tips are
#'   removed.  A complete match is silent.
#' @return A list with elements \code{tree}, \code{X}, \code{data},
#'   \code{report}, \code{matched_species}, \code{tree_tips_removed}, and
#'   \code{data_rows_removed}. \code{data} is an alias of \code{X}.
#' @export
match_tree_data <- function(tree, data = NULL, prune = TRUE, verbose = TRUE,
                            X = NULL) {
  if (is.null(data)) {
    data <- X
  } else if (!is.null(X)) {
    stop("supply data or its deprecated X alias, not both.", call. = FALSE)
  }
  .match_tree_data_core(tree = tree, data = data, prune = prune,
                        verbose = verbose)
}

.match_tree_data_core <- function(tree, data, prune = TRUE, verbose = TRUE,
                                  allow_insufficient = FALSE) {
  if (!is.logical(prune) || length(prune) != 1L || is.na(prune) ||
      !is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("prune and verbose must be TRUE or FALSE.", call. = FALSE)
  }
  ctx <- .as_prepared_tree(tree)
  tree <- ctx$tree
  X0 <- .as_named_trait_table(
    data, tree, verbose = FALSE, input_name = "X", vector_name = "x"
  )

  tree_tips <- as.character(tree$tip.label)
  data_names <- rownames(X0)
  matched <- tree_tips[tree_tips %in% data_names]
  removed_tree <- setdiff(tree_tips, data_names)
  removed_data <- setdiff(data_names, tree_tips)
  if (length(matched) < 2L && !isTRUE(allow_insufficient)) {
    stop("Fewer than two species are shared by tree and data.", call. = FALSE)
  }

  matched_tree <- if (isTRUE(prune) && length(matched) >= 2L) {
    .prepared_tree_subset(
      ctx, match(matched, tree_tips), need_matrix = FALSE
    )$tree
  } else {
    tree
  }
  matched_X <- if (length(matched)) X0[matched, , drop = FALSE] else
    X0[FALSE, , drop = FALSE]

  report <- data.frame(
    original_tree_tips = length(tree_tips),
    input_rows = nrow(X0),
    matched_species = length(matched),
    removed_tree_tips = length(removed_tree),
    removed_data_rows = length(removed_data),
    stringsAsFactors = FALSE
  )
  report$insufficient_retained <- length(matched) < 2L

  # A complete match is intentionally silent even when verbose = TRUE.  When
  # information is removed, emit exactly one message for this operation.
  if (isTRUE(verbose) && (length(removed_tree) || length(removed_data))) {
    message(sprintf(
      "Matched %d species; removed %d tree tips and %d data rows.",
      length(matched), length(removed_tree), length(removed_data)
    ))
  }

  out <- list(
    tree = matched_tree,
    X = matched_X,
    data = matched_X,
    report = report,
    matched_species = matched,
    tree_tips_removed = removed_tree,
    data_rows_removed = removed_data,
    base_keep = match(matched, tree_tips),
    insufficient_retained = length(matched) < 2L
  )
  class(out) <- "fastphylosig_match"
  out
}

#' @rdname match_tree_data
#' @export
match_phylo_data <- function(tree, X, prune = TRUE, verbose = TRUE) {
  .Deprecated("match_tree_data")
  match_tree_data(tree = tree, data = X, prune = prune, verbose = verbose)
}
