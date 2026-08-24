\name{fast_k}
\alias{fast_k}
\title{Fast Blomberg's K phylogenetic signal}
\usage{
fast_k(tree, x = NULL, test = FALSE, nsim = 1000, se = NULL,
  permutations = NULL, return_sim = FALSE, verbose = TRUE, ncores = 1,
  X = NULL, prepared = NULL, trait_chunk = 64, simulation_chunk = 128,
  keep_null = return_sim, progress = interactive(), data = NULL)
}
\arguments{
\item{tree}{A connected \code{"phylo"} object with a unique structural root
and finite, strictly positive branch lengths, or a context returned by
\code{prepare_tree()}. A conventional degree-two root, binary topology, and
ultrametricity are not required.}
\item{x, X, data}{Named continuous trait vector or trait matrix/data.frame.
The \code{data} spelling is an alias for \code{x}; supply only one spelling.}
\item{test}{Whether to compute the permutation significance test.}
\item{nsim}{Number of randomizations when \code{test = TRUE}.}
\item{se}{Currently only \code{NULL} is supported.}
\item{permutations}{Optional integer permutation matrix with exactly
\code{nsim} rows and one column per retained tip.}
\item{return_sim}{Whether to return simulated K values.}
\item{verbose}{Whether to report species matching.}
\item{ncores}{Number of OpenMP threads for randomization.}
\item{prepared}{Optional object returned by \code{prepare_tree()}; when
supplied it is used in place of \code{tree}.}
\item{trait_chunk}{Positive integer number of traits processed per C++ block.}
\item{simulation_chunk}{Positive integer number of randomizations processed
per streaming block.}
\item{keep_null}{Whether to retain the full K null distribution.}
\item{progress}{Whether to report method-level progress and elapsed runtime.}
}
\description{
Computes Blomberg's K using the validated tree-pruning production engine.
This is the fixed-method specialist counterpart of
\code{fast_signal(..., method = "K")}. The result class and shape are kept
compatible with the former continuous \code{fast_signal()} API. Matching,
trait-wise missing-value handling, and tree preparation follow the same
validated path as the unified dispatcher.
}
\value{
For vector input, a \code{phylosig}-compatible result (an atomic numeric when
\code{test = FALSE}); matrix/data.frame input returns one row per trait.
Results include timing metadata and a small \code{workflow} attribute.
Per-trait \code{status} and \code{message}/\code{note} fields identify
invalid, insufficient, or undefined traits; an atomic no-test result exposes
the same information as attributes.
}
\seealso{\code{fast_signal()}, \code{fast_lambda()}}
\examples{
tree <- ape::rtree(20)
x <- stats::rnorm(20)
names(x) <- tree$tip.label
fast_k(tree, x, verbose = FALSE)
}
