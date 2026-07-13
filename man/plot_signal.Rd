\name{plot_signal}
\alias{plot_signal}
\title{Plot fastphylosig results}
\usage{
plot_signal(result, p_col = NULL, alpha = 0.05, main = NULL,
  xlab = NULL, ylab = NULL,
  null = c("auto", "random", "brownian"),
  label = c("none", "significant", "top", "all"), top_n = 10,
  ridge_scale = 0.85, col = NULL, pch = 19, cex = 0.7, ...)
}
\arguments{
\item{result}{Result produced by \code{\link{fast_signal}},
\code{\link{fast_d}}, or \code{\link{fast_delta}}.}
\item{p_col}{Optional P-value column for data.frame results. For D results
the default is \code{"Pval1_fast"}.}
\item{alpha}{Significance level used for coloring observed signal lines.}
\item{null}{For D calibration plots, which null simulation to plot. \code{"auto"}
draws both the random and Brownian D null distributions for a single-trait
\code{fast_d()} result.}
\item{label}{Which observed traits to label in multi-trait ridge plots.}
\item{top_n}{Number of smallest-P traits to label when \code{label = "top"}.}
\item{ridge_scale}{Vertical scaling of each density ridge in multi-trait
plots.}
\item{main, xlab, ylab, col, pch, cex, ...}{Base plotting
arguments.}
}
\description{
Creates lightweight base R diagnostic plots for supported fastphylosig result
objects. For K and Delta, the x axis is the final phylogenetic-signal value and
the density/ridge shows simulated signal values from the permutation test; the
right-tail region \code{simulated >= observed} is shaded to show the P-value
tail. For single-trait D results, \code{null = "auto"} draws a Fritz & Purvis
D calibration plot: Brownian threshold and random association null
distributions are scaled onto the D axis and shown together as overlaid
mountain densities on the same baseline, the observed D is drawn as a vertical line, and reference lines mark
D = 0 (Brownian expectation) and D = 1 (random expectation). The internal
simulated contrast sums are converted back to D values before plotting.
\code{P_random} and \code{P_Brownian} are labelled separately; zero extreme
counts are displayed as \code{P < 1/nsim} rather than \code{P = 0}. For lambda, which
uses a deterministic likelihood-ratio test rather than random simulations, the
helper draws the stored profile log-likelihood curve when available, marks
\code{lambda_hat}, \code{lambda = 0}, \code{lambda = 1}, and an approximate
95 percent profile-likelihood confidence interval. The confidence interval is
the region with log-likelihood at least
\code{max(logL) - 0.5 * qchisq(0.95, df = 1)}. Old lambda results without a
stored profile fall back to a compact likelihood-ratio summary.
Observed signal values, P values, and simulation counts are printed on the
plots where applicable.
}
\value{
Invisibly returns the normalized plotting data. K and Delta plots include
\code{sim}, \code{null_model}, \code{n_sim}, and \code{tail_direction}.
Single-trait D overlay plots include \code{sim_random}, \code{sim_brownian},
\code{Pval1}, \code{Pval0}, \code{n_sim}, \code{plot_type = "D_calibration"},
\code{P_random_display}, and \code{P_Brownian_display}. Lambda plots include
\code{plot_type = "lambda_profile"} when a stored profile is available, plus
\code{profile}, \code{ci_lower}, \code{ci_upper}, and \code{LR}. Old lambda
results without a profile return \code{plot_type = "lambda_lr_summary"}.
}
\examples{
set.seed(7)
tree <- ape::rtree(20)
x <- stats::rnorm(20)
names(x) <- tree$tip.label
k <- fast_signal(tree, x, method = "K", test = TRUE, nsim = 20,
                 verbose = FALSE)
plot_signal(k)
}
