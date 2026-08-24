\name{fast_lambda}
\alias{fast_lambda}
\title{Fast Pagel's lambda phylogenetic signal}
\usage{
fast_lambda(tree, x = NULL, test = FALSE, se = NULL,
  verbose = TRUE, ncores = 1, X = NULL,
  lambda_profile = NULL, lambda_profile_points = 101,
  prepared = NULL, trait_chunk = 64, progress = interactive(), data = NULL)
}
\arguments{
\item{tree}{A connected \code{"phylo"} object with a unique structural root
and finite, non-negative branch lengths, or a context returned by
\code{prepare_tree()}. Zero terminal branches are not supported; zero internal
branches are allowed with a numerical warning. A conventional degree-two root,
binary topology, and ultrametricity are not required.}
\item{x, X, data}{Named continuous trait vector or trait matrix/data.frame.
The \code{data} spelling is an alias for \code{x}; supply only one spelling.}
\item{test}{Whether to compute the likelihood-ratio test against lambda zero.}
\item{se}{Currently only \code{NULL} is supported.}
\item{verbose}{Whether to report species matching.}
\item{ncores}{Number of OpenMP threads used by lambda kernels.}
\item{lambda_profile}{Whether to store a lambda likelihood profile. The
default stores it for vector input and skips it for batch input.}
\item{lambda_profile_points}{Number of regular profile grid points (at least
5). The fitted lambda, and lambda 1 when allowed by the tree, are inserted in
addition to the regular grid.}
\item{prepared}{Optional object returned by \code{prepare_tree()}; when
supplied it is used in place of \code{tree}.}
\item{trait_chunk}{Positive integer number of traits processed per C++ block.}
\item{progress}{Whether to report method-level progress and elapsed runtime.}
}
\description{
Computes Pagel's lambda using the validated tree-pruning production engine.
This is the fixed-method specialist counterpart of
\code{fast_signal(..., method = "lambda")}. The result class and shape are
kept compatible with the former continuous \code{fast_signal()} API. When a
profile is requested, the grid runs from zero to the tree-specific allowed
maximum. The nominal 95 percent profile interval uses
\code{max(logL) - 0.5 * qchisq(0.95, df = 1)} as its cutoff; bounds are found by
linear interpolation at threshold crossings and fall back to grid endpoints
when a crossing is absent. The interval is therefore an approximate
grid/profile interval, not an exact confidence interval.
}
\value{
For vector input, a \code{phylosig}-compatible result; matrix/data.frame input
returns one row per trait. Results include timing metadata and a small
\code{workflow} attribute. Per-trait \code{status} and
\code{message}/\code{note} fields distinguish valid fits from invalid,
insufficient, or undefined trait-level results.
}
\seealso{\code{fast_signal()}, \code{fast_k()}}
\examples{
tree <- ape::rtree(20)
x <- stats::rnorm(20)
names(x) <- tree$tip.label
fast_lambda(tree, x, verbose = FALSE)
}
