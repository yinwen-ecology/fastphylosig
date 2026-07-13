\name{fast_d}
\alias{fast_d}
\title{Fast caper-compatible phylo.d calculation}
\usage{
fast_d(tree, x = NULL, test = TRUE, nsim = 1000, se = NULL,
  return_sim = test, verbose = TRUE, rnd.bias = NULL, permutations = NULL,
  random_states = NULL, brownian_states = NULL, ncores = 1, X = NULL,
  permut = NULL)
}
\arguments{
\item{tree}{An object of class \code{"phylo"}.}
\item{x, X}{Named binary trait vector or matrix/data.frame of binary traits.}
\item{test}{Whether to return randomization and Brownian P values.}
\item{nsim, permut}{Number of simulations. \code{permut} is accepted as a
\code{caper::phylo.d()}-style alias.}
\item{se}{Unused; binary states should be supplied directly in \code{x}.}
\item{return_sim}{Whether to return random and Brownian simulation vectors.}
\item{verbose}{Whether to report species matching.}
\item{rnd.bias}{Optional sampling weights for the random null model.}
\item{permutations}{Optional permutation matrix for the random null model.}
\item{random_states, brownian_states}{Optional state matrices for controlled
debugging.}
\item{ncores}{Number of OpenMP threads for the C++ contrast-sum simulation
kernel. Binary trees can use multiple threads; polytomies fall back to the
original serial vectorized path to preserve caper-compatible handling.}
}
\description{
Computes the Fritz and Purvis D statistic for binary traits using the same
null-model structure as \code{caper::phylo.d()}. The contrast-sum calculation
is moved into C++ for repeated simulation work. Brownian null states are
simulated directly along tree edges and then thresholded by the observed state
proportion, avoiding dense VCV construction for this step.
Binary labels are recoded to equally spaced 0/1 states, making D invariant to
the original labels.
}
\value{
For vector input, a compact \code{phylo.d}-like object without caper's
\code{NodalVals} or embedded \code{comparative.data}. For matrix/data.frame
input, a data.frame with one row per binary trait.
}
\examples{
set.seed(3)
tree <- ape::rtree(30)
z <- ape::rTraitCont(tree)
x <- as.integer(z > stats::median(z))
names(x) <- tree$tip.label
fast_d(tree, x, nsim = 20, return_sim = FALSE, verbose = FALSE)
}
