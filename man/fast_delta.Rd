\name{fast_delta}
\alias{fast_delta}
\title{Fast Delta statistic for categorical traits}
\usage{
fast_delta(tree, x = NULL, test = FALSE, nsim = 1000, se = NULL,
  mcmc_sim = 10000, thin = 10, burn = 100, lambda0 = 0.1,
  proposal_sd = 0.5, entropy = c("LSE", "SE", "GINI"),
  model = "ARD", permutations = NULL, return_sim = test,
  verbose = TRUE, ncores = 1, X = NULL, prepared = NULL,
  progress = interactive())
}
\arguments{
\item{tree}{A conventionally rooted, fully binary \code{"phylo"} object with
finite, strictly positive branch lengths, or a compatible context returned by
\code{prepare_tree()}. The method does not require ultrametricity.}
\item{x, X}{Named categorical trait vector or matrix/data.frame. Matrix and
data.frame row names must be species names. \code{X} is a compatibility alias
for \code{x}; supply only one. Species identifiers are matched by name and
reordered to tree-tip order.}
\item{test}{Whether to compute a permutation P value.}
\item{nsim}{Number of trait permutations when \code{test = TRUE}.}
\item{se}{Unused; retained so unsupported measurement-error input fails
clearly. Use \code{proposal_sd} for the MCMC proposal.}
\item{mcmc_sim}{Number of MCMC iterations per chain and Delta calculation.}
\item{thin}{MCMC thinning interval.}
\item{burn}{First MCMC iteration eligible to be saved.}
\item{lambda0}{Rate of the exponential priors on alpha and beta.}
\item{proposal_sd}{Standard deviation of the log-normal random-walk proposal.}
\item{entropy}{Entropy transformation: \code{"LSE"}, normalized Shannon
entropy \code{"SE"}, or normalized Gini impurity \code{"GINI"}.}
\item{model}{Discrete transition model used by \code{fast_ace()}: \code{"ER"}
or \code{"ARD"}.}
\item{permutations}{Optional integer permutation matrix with exactly
\code{nsim} rows and one column per retained tip.}
\item{return_sim}{Whether to retain simulated Delta values.}
\item{verbose}{Whether to report species matching.}
\item{progress}{Whether to report method-level progress stages and elapsed
runtime. Defaults to \code{interactive()}. Setting this to \code{FALSE}
suppresses progress and species-matching messages; warnings and errors are
still signalled.}
\item{ncores}{Number of worker processes for the Delta permutation test. Small
jobs can be slower in parallel because process startup is not free. Replaying
stochastic results requires a fixed worker count; changing it repartitions
worker RNG streams and is not draw-for-draw reproducible.}
\item{prepared}{Optional object returned by \code{prepare_tree()}.}
}
\description{
Computes the single-tree Delta statistic for categorical traits. Ancestral
state probabilities are estimated with \code{fast_ace()}, node uncertainty is
transformed to entropy, and Delta is the posterior
\code{mean(beta) / mean(alpha)} from two chains. Ancestral reconstruction,
entropy, and MCMC calculations use Rcpp/RcppArmadillo. The ape comparison path
is retained internally for regression tests, not exposed as a user engine.
}
\value{
For vector input, a \code{phylo_delta}-compatible object; matrix/data.frame
input returns one row per trait. Results report Delta, alpha/beta chain means
and SDs, ESS, split R-hat, covariance-aware \code{MCSE_Delta}, requested and
successful iteration counts, and safe RNG metadata. Tests also report P,
\code{MCSE_P}, failed fits, and successful permutation counts. Inspect the
diagnostics before interpreting short runs. Results expose a
\code{timing} attribute with at least \code{total_elapsed} (seconds); vector
results also provide a \code{timing} list field. Per-trait \code{status} and
\code{message}/\code{note} fields distinguish successful, partial,
invalid-trait, insufficient-data, and failed calculations.

Fewer than two species shared globally by the tree and data is a call-level
error. After matching and trait-wise NA pruning, Delta returns a per-trait
insufficient-data result when a trait retains fewer than two species. A batch with one or more
traits retaining 2--19 species emits one warning and still returns estimates.
The metadata records all per-trait retained/removed counts and sample-size
status (including the fewer-than-two error), plus the affected trait details,
in
\code{attr(result, "analysis_metadata")$delta_retained_species}.

For a finite Delta estimate, non-finite chain diagnostics, ESS below 20, or
split R-hat above 1.1 emit one batch warning. This warning never reruns a
chain or alters requested iterations. Per-trait diagnostic warning details
are retained in
\code{attr(result, "analysis_metadata")$delta_mcmc_diagnostics}.

The public defaults include \code{test = FALSE}, \code{nsim = 1000},
\code{return_sim = test}, \code{mcmc_sim = 10000}, \code{thin = 10},
\code{burn = 100}, and \code{model = "ARD"}. Thus permutation P,
\code{MCSE_P}, and retained null values are absent unless testing is
enabled; null values additionally require \code{return_sim = TRUE}.
\code{MCSE_P} uses the finite successful permutation count. These are
computational defaults rather than a guarantee
of adequate mixing for every data set. Users should increase the chain when
ESS, split R-hat, or Monte Carlo error is inadequate.
}
\seealso{\code{fast_signal()}, \code{fast_ace()}}
\examples{
set.seed(4)
tree <- ape::rtree(25)
x <- sample(letters[1:3], ape::Ntip(tree), replace = TRUE)
names(x) <- tree$tip.label
fast_delta(tree, x, mcmc_sim = 10000, thin = 10, burn = 100,
           verbose = FALSE)
}
