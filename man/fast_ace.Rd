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
\code{phy}. The ACE likelihood itself
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
