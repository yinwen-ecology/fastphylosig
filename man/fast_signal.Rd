\name{fast_signal}
\alias{fast_signal}
\title{Unified fast phylogenetic signal dispatcher}
\usage{
fast_signal(tree, data = NULL, method = NULL, ..., x = NULL, X = NULL,
  prepared = NULL, trait_chunk = 64)
}
\arguments{
\item{tree}{A \code{"phylo"} object or a context returned by
\code{prepare_tree()}.}
\item{x, data, X}{Named trait vector or trait matrix/data.frame. The
\code{data} spelling is accepted by the unified dispatcher; \code{x} and
\code{X} are retained as compatibility aliases. Matrix and data.frame row
names must be species names. Supply only one spelling.}
\item{method}{Exactly one method selector. Supported aliases are \code{K}/\code{k},
\code{lambda}/\code{Lambda}, the Unicode \code{U+03BB} alias, \code{D}/\code{d}, and
\code{Delta}/\code{delta}. The method is explicit and is never inferred from
the trait type; omitting it or supplying an unknown selector is an error.}
\item{prepared, trait_chunk}{Compatibility formals for prepared tree contexts
and C++ trait blocks; prepared contexts undergo V2 integrity validation, and
the formals are forwarded only to targets that support them.}
\item{...}{Arguments for the selected production function. K and lambda
accept the arguments documented by \code{fast_k()} and \code{fast_lambda()};
D and Delta accept those documented by \code{fast_d()} and
\code{fast_delta()}, respectively. Unknown or inapplicable arguments are
reported as errors. Legacy positional arguments after \code{method} remain
supported. When \code{data}, \code{x}, and \code{X} are supplied more than
once, the call fails rather than choosing one silently.}
}
\description{
Provides the one-command entry point for phylogenetic signal analysis and
dispatches to exactly one of the four public production functions:
\code{fast_k()}, \code{fast_lambda()}, \code{fast_d()}, or
\code{fast_delta()}. The selected function owns validation, progress output,
runtime timing, numerical calculation, and result construction. Method
aliases are normalized to \code{"K"}, \code{"lambda"}, \code{"D"}, or
\code{"Delta"}. Ordinary calls perform tree inspection, conservative
representation handling, explicit species matching, and trait-wise missing-
value grouping before the selected production kernel runs. No biological tree
repair is inferred from the method or trait type.
}
\details{
For a first analysis, supply a tree, named trait data, and one method selector
to this function. Advanced users can call \code{fast_k()},
\code{fast_lambda()}, \code{fast_d()}, or \code{fast_delta()} directly when
they need method-specific simulation, profile, null-retention, MCMC, or
diagnostic controls. \code{fast_ace()} is the specialist ancestral-state ML
function used by the Delta workflow.

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
points may be missed, and non-finite evaluations are not interpolated across.
Increase \code{lambda_profile_points} to assess grid sensitivity.
}
\value{
The original result class and shape returned by the selected production
function (for example \code{phylosig}, \code{phylo.d}, or
\code{phylo_delta}). A small \code{workflow} attribute records the canonical
method and production function. Lambda results may include an approximate
profile and nominal 95 percent interval: the profile uses a regular grid from
zero to the tree-specific allowed maximum, with the fitted value (and one when
allowed) inserted. Its cutoff is \code{max(logL) - 0.5 * qchisq(0.95, df = 1)};
limits are estimated by linear interpolation at profile threshold crossings,
or by grid endpoints when no crossing is available. This is a grid/profile
approximation rather than an exact confidence interval.
}
\seealso{
\code{fast_k()}, \code{fast_lambda()}, \code{fast_d()},
\code{fast_delta()}, \code{fast_ace()}, \code{check_tree()},
\code{match_tree_data()}, \code{prepare_tree()}
}
\examples{
set.seed(2)
tree <- ape::rtree(30)
X <- matrix(stats::rnorm(30 * 3), 30, 3,
            dimnames = list(tree$tip.label, paste0("trait_", 1:3)))
fit <- fast_signal(tree, data = X, method = "K", test = FALSE,
                   verbose = FALSE)
summary(fit)
fast_signal(tree, data = X[, 1], method = "lambda", verbose = FALSE)
}
