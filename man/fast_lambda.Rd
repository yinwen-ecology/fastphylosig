\name{fast_lambda}
\alias{fast_lambda}
\title{Fast Pagel's lambda phylogenetic signal}
\usage{
fast_lambda(tree, x = NULL, test = FALSE, se = NULL,
  verbose = TRUE, ncores = 1, X = NULL,
  lambda_profile = NULL, lambda_profile_points = 101,
  prepared = NULL, trait_chunk = 64, progress = interactive(), data = NULL)
}
\arguments{
\item{tree}{A connected \code{"phylo"} object with a unique structural root
and finite, non-negative branch lengths, or a context returned by
\code{prepare_tree()}. Zero terminal branches are not supported; zero internal
branches are allowed with a numerical warning. A conventional degree-two root,
binary topology, and ultrametricity are not required.}
\item{x, X, data}{Named continuous trait vector or trait matrix/data.frame.
The \code{data} spelling is an alias for \code{x}; supply only one spelling.}
\item{test}{Whether to compute the likelihood-ratio test against lambda zero.}
\item{se}{Currently only \code{NULL} is supported.}
\item{verbose}{Whether to report species matching.}
\item{ncores}{Accepted for API compatibility. The current fixed-lambda and
optimization kernels run serially; increasing this value does not enable
parallel lambda evaluation.}
\item{lambda_profile}{Whether to store a lambda likelihood profile. The
default stores it for vector input and skips it for batch input.}
\item{lambda_profile_points}{Number of regular profile grid points (at least
5). The fitted lambda, and lambda 1 when allowed by the tree, are inserted in
addition to the regular grid. The optimizer search knots, including
geometrically spaced points near both boundaries, are also included.}
\item{prepared}{Optional object returned by \code{prepare_tree()}; when
supplied it is used in place of \code{tree} after V2 integrity validation.}
\item{trait_chunk}{Positive integer number of traits processed per C++ block.}
\item{progress}{Whether to report method-level progress and elapsed runtime.}
}
\description{
Computes Pagel's lambda using the validated tree-pruning production engine.
This is the fixed-method specialist counterpart of
\code{fast_signal(..., method = "lambda")}. The result class and shape are
kept compatible with the former continuous \code{fast_signal()} API. When a
profile is requested, the grid runs from zero to the tree-specific allowed
maximum. The nominal 95 percent profile interval uses
\code{max(logL) - 0.5 * qchisq(0.95, df = 1)} as its cutoff; bounds are found by
numerically refining threshold crossings bracketed by the grid using the
profile likelihood. If refinement fails, a warning is issued and linear
interpolation is used. Bounds fall back to grid endpoints when a crossing
is absent. This is an approximate likelihood-ratio interval, not an exact
finite-sample confidence interval.
}
\value{
For vector input, a \code{phylosig}-compatible result; matrix/data.frame input
returns one row per trait. Results include timing metadata and a small
\code{workflow} attribute. Per-trait \code{status} and
\code{message}/\code{note} fields distinguish valid fits from invalid,
insufficient, or undefined trait-level results.
For extreme finite trait scales, the fitted variance may overflow to
\code{Inf} or underflow to zero while lambda and log likelihood remain
valid. The \code{message}/\code{note} field records this variance-range
limitation; it does not indicate a constant trait.
}
\details{
The approximate likelihood-ratio confidence set can have disconnected
components. \code{lambda_CI_intervals} contains all grid-detected components
(as a data frame with \code{lower} and \code{upper} columns), and
\code{lambda_CI_disconnected} identifies a set with more than one component.
The legacy \code{lambda_CI} bounds give the envelope of this set; values in
between disconnected components are not included in the confidence set.
For table results the corresponding fields are
\code{lambda_CI_intervals_fast}, \code{lambda_CI_disconnected_fast}, and
\code{lambda_CI_lower_fast}/\code{lambda_CI_upper_fast}.
The set remains a profile-grid approximation: narrow components between grid
points may be missed, as may narrow excluded gaps between components.
Boundary-enriched search knots reduce this risk but do not eliminate it.
Non-finite evaluations are not interpolated across.
Increase \code{lambda_profile_points} to assess grid sensitivity.
}
\seealso{\code{fast_signal()}, \code{fast_k()}}
\examples{
tree <- ape::rtree(20)
x <- stats::rnorm(20)
names(x) <- tree$tip.label
fast_lambda(tree, x, verbose = FALSE)
}
