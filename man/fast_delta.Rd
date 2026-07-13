\name{fast_delta}
\alias{fast_delta}
\title{Fast Delta statistic for categorical traits}
\usage{
fast_delta(tree, x = NULL, test = FALSE, nsim = 1000, se = NULL,
  mcmc_sim = 10000, thin = 10, burn = 100, lambda0 = 0.1,
  proposal_sd = 0.5, entropy = c("LSE", "SE", "GINI"),
  model = "ARD", ace_engine = c("fast", "ape"),
  permutations = NULL, return_sim = test, verbose = TRUE, ncores = 1,
  X = NULL)
}
\arguments{
\item{tree}{An object of class \code{"phylo"} with strictly positive branch
lengths.}
\item{x, X}{Named categorical trait vector or matrix/data.frame. Rows must be
species names for matrix/data.frame input.}
\item{test}{Whether to compute a permutation P value.}
\item{nsim}{Number of trait permutations when \code{test = TRUE}.}
\item{se}{Unused; retained only to make unsupported measurement-error input
fail clearly. Use \code{proposal_sd} for the MCMC proposal.}
\item{mcmc_sim}{Number of MCMC iterations for each Delta calculation.}
\item{thin}{MCMC thinning interval.}
\item{burn}{First MCMC iteration eligible to be saved.}
\item{lambda0}{Rate parameter of the exponential priors on alpha and beta.}
\item{proposal_sd}{Standard deviation of the log-normal random-walk proposal.}
\item{entropy}{Entropy transformation: \code{"LSE"} for the linear Shannon
entropy used in Borges et al. (2019), \code{"SE"} for normalized Shannon
entropy, or \code{"GINI"} for normalized Gini impurity.}
\item{model}{Discrete transition model passed to \code{\link{fast_ace}} or,
when \code{ace_engine = "ape"}, to \code{ape::ace()}. Defaults to
\code{"ARD"}.}
\item{ace_engine}{Use \code{"fast"} for \code{\link{fast_ace}} or
\code{"ape"} for the reference \code{ape::ace()} implementation.}
\item{permutations}{Optional permutation matrix for reproducible null
distributions. It must have exactly \code{nsim} rows; every row must be a
complete 1-based permutation of the retained species.}
\item{return_sim}{Whether to return simulated Delta values when
\code{test = TRUE}.}
\item{verbose}{Whether to report species matching.}
\item{ncores}{Number of worker processes for the Delta permutation test. This
is used only when \code{test = TRUE}; the observed single-tree Delta estimate
is unchanged.}
}
\description{
Computes the Delta statistic for categorical traits. The implementation follows
the single-tree workflow from Borges et al. (2019): ancestral state
probabilities are estimated with \code{\link{fast_ace}} by default, node uncertainties are
transformed to entropy values, and Delta is estimated as the posterior
\code{mean(beta) / mean(alpha)} from two Metropolis-Hastings/Gibbs chains. The
ancestral reconstruction, entropy, and MCMC steps are implemented in
Rcpp/RcppArmadillo when \code{ace_engine = "fast"}.
}
\value{
For vector input, a \code{phylo_delta} object. For matrix/data.frame input, a
data.frame with one row per trait. When \code{test = TRUE}, both forms include
\code{n_failed_sim}, which reports non-finite permutation fits; the P value
uses the remaining finite fits and is missing if all fits fail.
}
\examples{
set.seed(4)
tree <- ape::rtree(25)
x <- sample(letters[1:3], ape::Ntip(tree), replace = TRUE)
names(x) <- tree$tip.label
fast_delta(tree, x, mcmc_sim = 200, thin = 10, burn = 20,
           verbose = FALSE)
}
