\name{check_tree}
\alias{check_tree}
\title{Check a phylogenetic tree before signal calculations}
\description{
  Inspect a phylogenetic tree without reordering, pruning, compiling, or
  otherwise changing it.  Readiness is reported separately for each selected
  signal method.  K, lambda, and D can use compatibility paths for polytomies;
  Delta follows \code{fast_ace()} and requires a rooted,
  fully-dichotomous tree with strictly positive branch lengths.  The summary
  reports both the conventional degree-two rooted flag and the structural root
  used by the compatibility kernels, so a root polytomy or an unrooted
  encoding is visible rather than silently treated as a binary rooted tree.
}
\usage{
check_tree(tree, signal = c("K", "lambda", "D", "Delta"))
}
\arguments{
  \item{tree}{A \code{phylo} object or an object returned by
    \code{prepare_tree()}.}
  \item{signal}{Character vector selecting one or more methods.  The default
    checks \code{"K"}, \code{"lambda"}, \code{"D"}, and \code{"Delta"}.}
}
\value{
  An object of class \code{fastphylosig_tree_check}.  Its \code{ready} field
  is a logical scalar and \code{ready_by_signal} is a named logical vector.
  \code{tree_summary} contains tip, topology, rootedness, ultrametricity,
  branch-length, polytomy, unary-node, and zero/near-zero branch diagnostics.
  Zero internal branches are conditional for lambda; zero terminal branches
  require user review at the lambda upper boundary.  The \code{issues} field
  is a data frame with \code{code}, \code{severity}, \code{signal},
  \code{message}, \code{check}, \code{problem}, \code{action}, and
  \code{auto_fixable} columns.  Severity is one of \code{INFO},
  \code{WARNING}, \code{ERROR}, or \code{USER_ACTION_REQUIRED}.
  Representation-only canonical node-ID diagnostics may set
  \code{auto_fixable = TRUE}; biological rooting, topology, and branch-length
  repairs are never automatic.  The input tree is never modified.
}
\examples{
if (requireNamespace("ape", quietly = TRUE)) {
  tree <- ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
  check_tree(tree)
  check_tree(tree, signal = c("K", "lambda"))
}
}
