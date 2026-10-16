\name{fast_ace}
\alias{fast_ace}
\title{Fast discrete ancestral-state reconstruction}
\usage{
fast_ace(x, phy = NULL, type = "discrete", method = "ML", CI = TRUE,
  model = c("ER", "ARD"), kappa = 1, ip = 0.1, marginal = FALSE,
  prepared = NULL, progress = interactive())
}
\arguments{
\item{x}{Categorical trait vector with one value per tip. Names, when present,
must be unique and have exactly the same species set as \code{phy$tip.label};
their order may differ and values are reordered to tree-tip order.}
\item{phy}{A rooted, fully dichotomous \code{"phylo"} object with finite,
strictly positive branch lengths.}
\item{type}{Currently only \code{"discrete"} is supported.}
\item{method}{Currently only \code{"ML"} is supported.}
\item{CI}{Whether to return ancestral likelihoods in \code{lik.anc}.}
\item{model}{Discrete transition model, \code{"ER"} or \code{"ARD"};
the default is \code{"ER"}.}
\item{kappa}{Optional branch-length exponent, matching \code{ape::ace()}.}
\item{ip}{Initial value for transition-rate optimization.}
\item{marginal}{Whether to return marginal likelihoods. The default
\code{FALSE} matches the default path used by \code{ape::ace()}.}
\item{prepared}{An optional object returned by \code{prepare_tree()}. When
supplied, it supplies a validated phylogenetic tree and takes precedence over
\code{phy}. V2 integrity is checked before cached evidence is reused. The ACE likelihood itself
uses the tree edges and does not require the dense covariance cache.}
\item{progress}{Whether to report method-level progress stages and elapsed
runtime. Defaults to \code{interactive()}. Setting this to \code{FALSE}
suppresses progress messages; warnings and errors are still signalled.}
}
\description{
Optimized subset of \code{ape::ace()} for discrete maximum-likelihood
ancestral-state reconstruction. The statistical target follows
\code{ape::ace(type = "discrete", method = "ML")} for \code{"ER"} and
\code{"ARD"} models, while the repeated pruning likelihood calculation is
implemented in Rcpp/RcppArmadillo.
}
\details{
Discrete fits retain the spectral pruning arithmetic for ordinary
eigensystems. Defective or ill-conditioned eigensystems use deterministic
matrix-exponential recovery; recovery counts are distinct from final optimizer
nonconvergence. Returned ancestral probability rows must be finite,
normalized and inside [0,1]. Only machine-scale boundary excursions may be
clipped and renormalized; material invalidity triggers a stable recomputation
or an explicit failure. These abnormal paths can change estimates; no random
starts or random numerical recovery are used.

Likelihood-only fits (\code{CI = FALSE}, including Delta's internal ACE fits)
reuse a private native workspace. Sufficiently large uncached branch-length
workloads evaluate independent complex exponent vectors in OpenMP batches
bounded by 2 MiB. Matrix products, pruning accumulation and optimizer calls
retain their sequential order. No random draws occur in the parallel region.
The \code{fastphylosig.ace_threads} option defaults to \code{0L} and selects at
most two threads. Set this option to \code{1L} for serial fitting,
or a positive integer to request more threads, capped by the OpenMP runtime
limit. Small fits, builds without OpenMP and nested OpenMP calls use the serial
path. Exactly repeated branch lengths use an exponent-vector cache bounded by
64 MiB per fitting workspace, with capacity determined by state count and
actual element size. This bounds exponent storage, not total memory;
mixed-length fits batch only the remaining uncached branches using a compact
worklist. The activation threshold counts uncached work, so a large tree with
very few uncached branches does not incur parallel startup costs.
The default \code{CI = TRUE} path is unaffected by this thread option.
}
\value{
An \code{ace}-like object containing \code{loglik}, \code{rates}, numerical
rate standard errors in \code{se}, \code{index.matrix}, optimizer status, and,
when \code{CI = TRUE}, \code{lik.anc}. A \code{timing} attribute and list
field report elapsed runtime in seconds.
}
\seealso{\code{fast_delta()}, \code{fast_signal()}}
\examples{
set.seed(5)
tree <- ape::rtree(30)
x <- ape::rTraitDisc(tree, model = "ER", k = 3, rate = 1,
                     states = letters[1:3])
fast_ace(x, tree, model = "ER")
}
