\name{result_methods}
\alias{print.fastphylosig_result}
\alias{summary.fastphylosig_result}
\alias{as.data.frame.fastphylosig_result}
\alias{print.fastphylosig_table}
\alias{summary.fastphylosig_table}
\alias{as.data.frame.fastphylosig_table}
\title{Stable result classes and display methods}
\description{
The public calculation functions retain their historical fields and legacy
classes while adding a common result class. Scalar results inherit from
\code{fastphylosig_result} after a statistic-specific class
(\code{fastphylosig_signal}, \code{fastphylosig_d}, or
\code{fastphylosig_delta}) and before the corresponding legacy class. Matrix
and data-frame results inherit from
\code{c("fastphylosig_table", "fastphylosig_result", "data.frame")}.
}
\usage{
\method{print}{fastphylosig_result}(x, ..., digits = getOption("digits"))
\method{summary}{fastphylosig_result}(object, ...)
\method{as.data.frame}{fastphylosig_result}(x, row.names = NULL,
  optional = FALSE, ...)
}
\arguments{
\item{x}{A scalar or table result returned by one of the fastphylosig
calculation functions.}
\item{object}{A result supplied to \code{summary()}.}
\item{row.names}{Optional row names for the data-frame representation.}
\item{optional}{Logical flag accepted for compatibility with the base
\code{as.data.frame()} generic.}
\item{digits}{Number of significant digits used by \code{print()}.}
\item{...}{Additional arguments passed to the underlying data-frame printer.}
}
\details{
Canonical aliases are additive. K results expose \code{estimate},
\code{p_value}, \code{p_mcse}, \code{n_sim_requested}, and
\code{n_sim_successful}; lambda results expose \code{estimate},
\code{logLik}, \code{LR}, and \code{p_value}. D results keep the random and
Brownian P values and Monte Carlo errors separate (\code{P_random},
\code{P_Brownian}, \code{MCSE_P_random}, and
\code{MCSE_P_Brownian}); no combined \code{p_value} is created. Delta
results expose \code{estimate}, \code{p_value}, \code{p_mcse},
\code{estimate_mcse}, chain-specific ESS/R-hat aliases, and simulation or
iteration count aliases.

The no-test K scalar remains an atomic numeric with its legacy
\code{phylosig} class. The display methods derive its estimate directly and
do not coerce it to a list. Printing, summarising, and coercing results to a
data frame use a compact canonical schema and omit retained simulation
vectors.
}
\value{
\code{print()} returns the input invisibly. \code{summary()} and
\code{as.data.frame()} return a compact data frame; the latter has one row
for a scalar result and preserves all non-simulation columns for a result
table.
}
\seealso{\code{fast_signal}, \code{fast_k}, \code{fast_lambda},
\code{fast_d}, \code{fast_delta}}
