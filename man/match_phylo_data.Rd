\name{match_phylo_data}
\alias{match_phylo_data}
\title{Deprecated alias for \code{match_tree_data}}
\usage{
match_phylo_data(tree, X, prune = TRUE, verbose = TRUE)
}
\arguments{
\item{tree}{An object of class \code{"phylo"} or a prepared tree context.}
\item{X}{Named vector, matrix, or data.frame. Species names must be supplied
as names or row names unless rows are already in \code{tree$tip.label} order.}
\item{prune}{Whether to drop tree tips that are absent from \code{X}.}
\item{verbose}{Whether to print a short matching message.}
}
\description{
This compatibility alias is deprecated. Use \code{match_tree_data()} instead.
The function returns the same object as \code{match_tree_data()} and emits a
deprecation warning when called.
}
\value{
See \code{\link{match_tree_data}}.
}
\examples{
set.seed(6)
tree <- ape::rtree(10)
X <- data.frame(trait = stats::rnorm(10), row.names = tree$tip.label)
\dontrun{
matched <- match_phylo_data(tree, X, verbose = FALSE)
matched$report
}
}
