\name{match_tree_data}
\alias{match_tree_data}
\title{Match a phylogeny and trait table}
\usage{
match_tree_data(tree, data = NULL, prune = TRUE, verbose = TRUE, X = NULL)
}
\arguments{
\item{tree}{An object of class \code{"phylo"} or a prepared tree context.}
\item{data}{Named vector, matrix, or data.frame. Species names must be
supplied as names or row names unless rows are already in
\code{tree$tip.label} order.}
\item{X}{Deprecated compatibility alias for \code{data}. Supply only one of
\code{data} and \code{X}.}
\item{prune}{Whether to drop tree tips that are absent from \code{X}.}
\item{verbose}{Whether to print a short matching message.}
}
\description{
Matches a phylogeny and trait table without silently discarding information.
The returned object contains the matched tree, the trait table reordered in
tree-tip order, and a report listing matched species, removed tree tips, and
removed data rows.
}
\value{
A list with elements \code{tree}, \code{X}, \code{data}, \code{report},
\code{matched_species}, \code{tree_tips_removed}, and
\code{data_rows_removed}. \code{data} is an alias of \code{X} for the
recommended check-resolve-match-calculate workflow.
}
\examples{
set.seed(6)
tree <- ape::rtree(10)
X <- data.frame(trait = stats::rnorm(10), row.names = tree$tip.label)
matched <- match_tree_data(tree, data = X, verbose = FALSE)
matched$report
}
