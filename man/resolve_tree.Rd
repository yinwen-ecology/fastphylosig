\name{resolve_tree}
\alias{resolve_tree}
\title{Resolve a phylogenetic tree conservatively}
\usage{
resolve_tree(tree, signal, outgroup = NULL)
}
\arguments{
  \item{tree}{An object of class \code{"phylo"}.}
  \item{signal}{One of \code{"K"}, \code{"lambda"}, \code{"D"}, or
  \code{"Delta"}; matching is case-insensitive and exactly one value is
  required.}
  \item{outgroup}{Optional character vector of one or more existing tip
  labels. It is used only when the selected signal requires a rooted tree.}
}
\value{
A \code{"phylo"} object with class \code{"fastphylosig_resolved_tree"} and
an attribute named \code{"fastphylosig_resolution"}. The attribute records
the signal, original and final readiness statuses, representation changes,
and a logical \code{ready} flag.
}
\description{
Checks a tree with \code{check_tree()} and makes only safe representation
normalizations (edge ordering and internal node numbering). An unrooted tree
is rooted only when the selected signal requires it and an explicit valid
outgroup is supplied. The function does not midpoint-root, infer an outgroup,
resolve polytomies, collapse unary nodes, modify branch lengths, jitter zero
branches, or drop tips. Remaining blocking conditions produce an explicit
\code{USER_ACTION_REQUIRED} error.
}
\details{
After every representation or rooting change the tree is checked again. A
successful result is idempotent: resolving it again does not change its
topology, tip labels, or branch-length paths.
}
\seealso{\code{check_tree}, \code{prepare_tree}}
\examples{
\dontrun{
tree <- ape::rtree(8)
tree <- resolve_tree(tree, signal = "K")
attr(tree, "fastphylosig_resolution")
}
}
