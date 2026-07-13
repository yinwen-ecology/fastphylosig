\name{fast_signal}
\alias{fast_signal}
\title{Fast phytools-compatible continuous-trait phylogenetic signal}
\usage{
fast_signal(tree, x = NULL, method = c("K", "lambda"), test = FALSE,
  nsim = 1000, se = NULL, start = NULL, control = list(),
  permutations = NULL, return_sim = test, verbose = TRUE, ncores = 1,
  X = NULL, lambda_profile = NULL, lambda_profile_points = 101)
}
\arguments{
\item{tree}{An object of class \code{"phylo"}.}
\item{x, X}{Named continuous trait vector or trait matrix/data.frame. Rows
must be species names for matrix/data.frame input.}
\item{method}{Either \code{"K"} or \code{"lambda"}.}
\item{test}{Whether to compute the significance test.}
\item{nsim}{Number of randomizations for K when \code{test = TRUE}.}
\item{se}{Currently only \code{NULL} is supported.}
\item{start, control}{Accepted for API similarity with
\code{phytools::phylosig()}, but unused when \code{se = NULL}.}
\item{permutations}{Optional permutation matrix for reproducible K
randomization. It must have exactly \code{nsim} rows; every row must be a
complete 1-based permutation of the retained species.}
\item{return_sim}{Whether to return simulated K values when
\code{test = TRUE}.}
\item{verbose}{Whether to report species matching.}
\item{ncores}{Number of OpenMP threads for C++ randomization work. Currently
used for \code{method = "K", test = TRUE}; deterministic K and lambda estimates
are unchanged and usually do not benefit from this argument.}
\item{lambda_profile}{Whether to store a likelihood profile grid for
\code{method = "lambda"}. The default \code{NULL} stores profiles for vector
input and skips them for matrix/data.frame input to avoid slowing large batch
jobs.}
\item{lambda_profile_points}{Number of grid points used for the stored lambda
profile. The optimized \code{lambda_hat}, and \code{lambda = 1} when inside the
optimization interval, are added to this grid.}
}
\description{
Computes Blomberg's K and Pagel's lambda with the same statistical formulas as
the inspected \code{phytools::phylosig()} implementation for \code{se = NULL}.
The speedup comes from reusing tree-level matrix work and moving repeated
trait arithmetic into Rcpp/RcppArmadillo. For K, the implementation evaluates
the same \code{C^{-1}} GLS mean and quadratic form with Cholesky triangular
solves rather than explicitly constructing \code{solve(C)}. Lambda likelihoods
also use Cholesky solves and log determinants instead of a full inverse at each
optimizer step.
}
\value{
For vector input, a compact \code{phylosig}-like object containing the main
estimate and requested test fields. For matrix/data.frame input, a data.frame
with one row per trait and columns containing fast estimates,
matching reports, and trait-wise missing-data information. Lambda vector input
also stores \code{lambda_profile} and an approximate 95 percent
profile-likelihood confidence interval for plotting.
}
\examples{
set.seed(2)
tree <- ape::rtree(30)
X <- matrix(stats::rnorm(30 * 3), 30, 3,
            dimnames = list(tree$tip.label, paste0("trait_", 1:3)))
fast_signal(tree, X, method = "K", test = FALSE, verbose = FALSE)
}
