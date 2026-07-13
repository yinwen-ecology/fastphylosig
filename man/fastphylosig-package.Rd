\name{fastphylosig-package}
\alias{fastphylosig-package}
\alias{fastphylosig}
\docType{package}
\title{Fast phylogenetic signal calculations}
\description{
Fast R and RcppArmadillo implementations of selected phylogenetic signal
statistics. Continuous-trait Blomberg's K and Pagel's lambda follow
\code{phytools::phylosig()} for \code{se = NULL}; binary Fritz and Purvis D
follows \code{caper::phylo.d()}; categorical Delta uses accelerated discrete
ancestral-state likelihoods. The package is designed for repeated calculations
on a shared tree, explicit tree/data matching, and trait-wise missing-data
handling.
}
\details{
The main user-facing functions are \code{\link{fast_signal}},
\code{\link{fast_d}}, \code{\link{fast_delta}}, and \code{\link{fast_ace}}.
Utility helpers include
\code{\link{match_phylo_data}} and \code{\link{plot_signal}}.
Continuous traits currently support Blomberg's K and Pagel's lambda through
\code{fast_signal()}.
Binary traits support a caper-compatible D statistic through \code{fast_d()}.
Categorical traits support the Delta statistic through \code{fast_delta()},
with discrete ancestral reconstruction accelerated by \code{fast_ace()}.

For continuous traits, K and lambda keep the \code{phytools::phylosig()}
statistical targets while using Cholesky solves to avoid unnecessary dense
matrix inverses.
\code{plot_signal()} provides method-specific diagnostics, including
permutation null plots for K and Delta, a profile-likelihood plot for lambda,
and a Fritz & Purvis D calibration plot for binary traits.
}
