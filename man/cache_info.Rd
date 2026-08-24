\name{cache_info}
\alias{cache_info}
\title{Inspect a prepared-tree cache}
\usage{
cache_info(ctx)
}
\arguments{
\item{ctx}{An object returned by \code{prepare_tree()}.}
}
\description{
Reports structural and numerical cache use without modifying the prepared tree.
Dense numerical entries are managed by a byte-budget least-recently-used cache.
}
\value{
A list containing \code{budget}, \code{bytes_used}, numbers of structural and
numerical entries, hit/miss/eviction counters, and an \code{entries} data.frame.
The table reports retained tip count, structural, numerical and total bytes,
VCV, Cholesky and eigensystem availability, and last-access order for each key.
}
\examples{
set.seed(1)
tree <- ape::rtree(20)
ctx <- prepare_tree(tree, cache_budget = 64 * 1024^2)
cache_info(ctx)
}
