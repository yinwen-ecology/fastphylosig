\name{prepare_tree}
\alias{prepare_tree}
\title{Prepare a phylogenetic tree for repeated signal calculations}
\description{
Builds a reusable tree context containing linear-size topology and traversal
data.
Dense covariance matrices, Cholesky factors, lambda spectral data, and
trait-specific numerical subsets are cached lazily when first requested.
Every external reuse validates the V2 context schema and exact protected
snapshot before cached or numerical work.
}
\usage{
prepare_tree(tree, max_cached_subsets = 16L,
  cache_budget = 512 * 1024^2)
}
\arguments{
\item{tree}{A rooted object of class \code{"phylo"} with finite,
non-negative branch lengths.}
\item{max_cached_subsets}{Maximum number of full or NA-pruned subset entries
retained in the structural cache, including the full-tree entry. Retained for
compatibility; dense memory is controlled separately by \code{cache_budget}.}
\item{cache_budget}{Maximum number of bytes retained in the numerical LRU
cache. A value of zero permits transient numerical calculations without
retaining their dense payloads.}
}
\value{
An object of class \code{c("fastphylosig_tree", "list")}. Pass it as the tree
argument to \code{fast_signal()}, \code{fast_k()}, \code{fast_lambda()},
\code{fast_d()}, and \code{fast_delta()}. For ACE, use
\code{fast_ace(x, prepared = ctx)}. The unified \code{fast_signal()} entry
point still requires an explicit \code{method}.
}
\details{
The tree context does not alter the tree or raw trait data. Initial preparation
stores a snapshot; later changes to the original \code{phylo} object do not
update it. On reuse, a versioned exact protected snapshot covers the cached
tree's tip labels, edge matrix, branch lengths, and \code{Nnode}. Changes to
those fields inside the context are rejected before cache or numerical work,
while \code{node.label} and other non-computational attributes are outside the
integrity contract; \code{node.label} may be reordered when internal nodes are
renumbered. V2 contexts can be saved with \code{saveRDS()} and reused
after \code{readRDS()}; missing, legacy, or unknown context schemas must be
rebuilt with \code{prepare_tree()}. Initial preparation stores only topology,
tip mappings, and traversal data. Public production K and lambda use the
tree-pruning engines and do not request dense covariance, Cholesky, or spectral
payloads; D, Delta, and ACE do not request them either. The bounded numerical
cache remains available to internal validation paths. Results retain the statistical
formulas of the corresponding reference functions. See \code{\link{cache_info}}
for cache diagnostics.
}
\examples{
set.seed(1)
tree <- ape::rtree(20)
ctx <- prepare_tree(tree)
X <- matrix(stats::rnorm(40), nrow = 20,
            dimnames = list(tree$tip.label, c("trait_1", "trait_2")))
fast_signal(ctx, data = X, method = "K", verbose = FALSE)
}
