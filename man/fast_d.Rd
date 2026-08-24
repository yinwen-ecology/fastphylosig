\name{fast_d}
\alias{fast_d}
\title{Fast caper-compatible phylo.d calculation}
\usage{
fast_d(tree, x = NULL, test = TRUE, nsim = 1000, se = NULL,
  return_sim = test, verbose = TRUE, rnd.bias = NULL, permutations = NULL,
  random_states = NULL, brownian_states = NULL, ncores = 1, X = NULL,
  permut = NULL, prepared = NULL, chunk_size = 128,
  keep_null = return_sim, progress = interactive())
}
\arguments{
\item{tree}{An object of class \code{"phylo"} with finite, strictly positive
branch lengths and no single-child internal nodes. Safe representation-only
differences in root numbering and edge order are normalized automatically.
Ordinary binary trees and rooted polytomies are supported; ultrametricity is
not required, while unsafe topology changes remain errors.}
\item{x, X}{Named binary trait vector or matrix/data.frame of binary traits.
\code{X} is a compatibility alias for \code{x}; supply only one. Species
identifiers are matched by name and reordered to tree-tip order.}
\item{test}{Whether to return randomization and Brownian P values.}
\item{nsim, permut}{Number of simulations. \code{permut} is a
\code{caper::phylo.d()}-style compatibility alias for \code{nsim}; when
supplied, it takes precedence.}
\item{se}{Unused; binary states should be supplied directly in \code{x}.}
\item{return_sim}{Whether random and Brownian simulation vectors may be
returned. Its default is \code{test}; actual retention also requires
\code{keep_null = TRUE}.}
\item{verbose}{Whether to report species matching.}
\item{progress}{Whether to report method-level progress stages and elapsed
runtime. Defaults to \code{interactive()}. Setting this to \code{FALSE}
suppresses progress and species-matching messages; warnings and errors are
still signalled.}
\item{rnd.bias}{Optional sampling weights for the random null model.}
\item{permutations}{Optional permutation matrix for the random null model.}
\item{random_states, brownian_states}{Optional state matrices for controlled
debugging.}
\item{ncores}{Number of OpenMP threads for the C++ contrast-sum simulation
kernel. Binary trees can use multiple threads, although scaling is often modest
because traversal is memory-bandwidth bound. Polytomies fall back to the
serial vectorized path to preserve caper-compatible handling.}
\item{prepared}{An optional object returned by \code{prepare_tree()}. When
supplied, its cached tree calculations are reused and \code{tree} is ignored.}
\item{chunk_size}{Integer from 32 to 512 controlling the number of null states
generated and traversed per streaming block.}
\item{keep_null}{Whether full random and Brownian null-statistic vectors may
be retained. Both \code{return_sim = TRUE} and \code{keep_null = TRUE} are
required. Set either to \code{FALSE} to keep only online means, tail counts,
D, and P values.}
}
\description{
Computes the Fritz and Purvis D statistic for binary traits using the same
null-model structure as \code{caper::phylo.d()}. The contrast-sum calculation
is moved into C++ for repeated simulation work. Brownian null states are
simulated directly along tree edges and then thresholded by the observed state
proportion, avoiding dense VCV construction for this step.
Binary labels are recoded to equally spaced 0/1 states, making D invariant to
the original labels.
For ordinary binary trees without controlled state matrices, null generation,
thresholding, traversal, and accumulation are fused in chunks. Peak null-state
memory is therefore O(n times \code{chunk_size}) rather than O(n times
\code{nsim}). Controlled and polytomy cases use the compatibility path.
}
\value{
For vector input, a compact \code{phylo.d}-like object without caper's
\code{NodalVals} or embedded \code{comparative.data}. In addition to the
legacy \code{Pval1}/\code{Pval0} fields, it reports the null-specific aliases
\code{P_random}/\code{P_Brownian}, their Bernoulli Monte Carlo standard errors
\code{MCSE_P_random}/\code{MCSE_P_Brownian}, and
\code{nsim_requested}, \code{nsim_successful_random},
\code{nsim_successful_brownian}, \code{nsim_failed_random}, and
\code{nsim_failed_brownian}. The two P values are never combined. For
matrix/data.frame input, a data.frame with one row per binary trait and the
same fields (with the legacy \code{Pval1_fast}/\code{Pval0_fast} columns
retained for compatibility).

Rooted polytomies use the serial compatibility path. Unary internal nodes are
not supported, including unary nodes exposed after species matching or
trait-wise missing-value pruning.

Every trait also has a machine-readable \code{status} (\code{"ok"},
\code{"partial"}, \code{"invalid_trait"}, \code{"insufficient_data"}, or
\code{"failed"}) and a human-readable \code{message}/\code{note}. These fields
distinguish invalid traits and incomplete simulation sets from a valid result.

The MCSE uses
\code{sqrt(p * (1 - p) / nsim_successful)} separately for each null. Thus
the reported MCSE is exactly zero when the estimated P is 0 or 1, and is
\code{NA} when no finite simulation succeeded. Non-finite simulation
statistics are excluded from that null's denominator and explicitly reported
in the successful/failed counts and \code{note}; they are not silently counted
as non-extreme values. With \code{test = FALSE}, P and MCSE columns are
returned as \code{NA}, although simulation counts still describe the null
calibration used to compute D. Results expose a \code{timing} attribute with
at least \code{total_elapsed} (seconds); vector results also provide a
\code{timing} list field.
}
\seealso{\code{fast_signal()}, \code{check_tree()}}
\examples{
set.seed(3)
tree <- ape::rtree(30)
z <- ape::rTraitCont(tree)
x <- as.integer(z > stats::median(z))
names(x) <- tree$tip.label
fast_d(tree, x, nsim = 20, return_sim = FALSE, verbose = FALSE)
}
