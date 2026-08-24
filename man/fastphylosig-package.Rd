\name{fastphylosig-package}
\alias{fastphylosig-package}
\alias{fastphylosig}
\docType{package}
\title{Fast phylogenetic signal calculations}
\description{
One-command and method-specific workflows for selected phylogenetic signal
statistics. The method-explicit \code{fast_signal(tree, data, method, ...)}
entry point checks and matches the inputs, then runs K, lambda, D, or Delta.
Continuous-trait Blomberg's K and Pagel's lambda preserve the
\code{phytools::phylosig()} targets for \code{se = NULL}; the binary Fritz and
Purvis D statistic follows \code{caper::phylo.d()}; categorical Delta uses accelerated
discrete ancestral-state likelihoods and reports MCMC diagnostics.
}
\details{
Start with \code{\link{fast_signal}} for a standard analysis. Advanced users
can run tree and data preparation explicitly, then call \code{\link{fast_k}},
\code{\link{fast_lambda}}, \code{\link{fast_d}}, \code{\link{fast_delta}},
and \code{\link{fast_ace}}. Repeated analyses can reuse
\code{\link{prepare_tree}}. Utilities include \code{\link{cache_info}},
\code{\link{check_tree}}, \code{\link{resolve_tree}},
\code{\link{match_tree_data}}, and \code{\link{plot_signal}}.

Production K and lambda calculations use compiled tree-pruning algorithms. K
randomizations are evaluated in bounded chunks, and lambda optimization is
completed within one compiled batch call. Public utilities provide
method-specific tree checks, conservative representation normalization,
progress output, runtime reporting, and structured results without exposing
numerical-engine choices.

Missing values are handled trait by trait. Retained-tip patterns are grouped,
and species matching/removal counts are reported. D retains separate random
and Brownian null P values. K, D, and Delta report Monte Carlo precision where
applicable; Delta additionally reports ESS, split R-hat, and estimate MCSE.

\code{plot_signal()} provides statistic-specific diagnostics: a permutation
null plot for K, a profile-likelihood plot for lambda, a two-null Fritz and
Purvis D calibration plot, and a permutation null plot for Delta.
}
