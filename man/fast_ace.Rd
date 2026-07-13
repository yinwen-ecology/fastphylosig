\name{fast_ace}
\alias{fast_ace}
\title{Fast discrete ancestral-state reconstruction}
\usage{
fast_ace(x, phy, type = "discrete", method = "ML", CI = TRUE,
  model = c("ER", "ARD"), kappa = 1, ip = 0.1, marginal = FALSE)
}
\arguments{
\item{x}{Categorical trait vector with one value per tip. Names, when present,
must be unique and match \code{phy$tip.label} exactly.}
\item{phy}{A rooted, fully dichotomous \code{"phylo"} object with branch
lengths.}
\item{type}{Currently only \code{"discrete"} is supported.}
\item{method}{Currently only \code{"ML"} is supported.}
\item{CI}{Whether to return ancestral likelihoods in \code{lik.anc}.}
\item{model}{Discrete transition model, \code{"ER"} or \code{"ARD"}.}
\item{kappa}{Optional branch-length exponent, matching \code{ape::ace()}.}
\item{ip}{Initial value for transition-rate optimization.}
\item{marginal}{Whether to return marginal likelihoods. The default
\code{FALSE} matches the default path used by \code{ape::ace()}.}
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
when \code{CI = TRUE}, \code{lik.anc}.
}
\examples{
set.seed(5)
tree <- ape::rtree(30)
x <- ape::rTraitDisc(tree, model = "ER", k = 3, rate = 1,
                     states = letters[1:3])
fast_ace(x, tree, model = "ER")
}
