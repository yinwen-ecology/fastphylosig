#ifndef FASTPHYLOSIG_HETEROGENEITY_BOOTSTRAP_H
#define FASTPHYLOSIG_HETEROGENEITY_BOOTSTRAP_H
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "numeric_utils.h"

// Shared summaries for model-reference upper-tail dispersion bootstraps.
// Any failed replicate invalidates P; successful draws are never substituted
// for the requested reference distribution.
namespace fastphylosig {
inline double hetero_sd(const std::vector<double>& x) {
  if (x.size() < 2) return NA_REAL;
  long double mean = 0, m2 = 0;
  int n = 0;
  for (double value : x) {
    if (!std::isfinite(value)) return NA_REAL;
    const long double delta = value - mean;
    mean += delta / ++n; m2 += delta * (value - mean);
  }
  return std::sqrt(std::max(0.0L, m2 / (n - 1)));
}
inline Rcpp::List hetero_bootstrap_summary(const std::vector<double>& draws,
    double observed, bool keep_null) {
  std::vector<double> good;
  int count = 0;
  for (double value : draws) if (std::isfinite(value)) {
    good.push_back(value); if (value >= observed) ++count;
  }
  const int nsim = draws.size(), successful = good.size();
  double p = NA_REAL, mcse = NA_REAL, mean = NA_REAL, upper = NA_REAL;
  if (nsim >= 2 && successful == nsim && std::isfinite(observed)) {
    p = (1.0 + count) / (nsim + 1.0);
    const double q = static_cast<double>(count) / nsim;
    mcse = std::sqrt(nsim * q * (1 - q)) / (nsim + 1.0);
    long double sum = 0;
    for (double value : good) sum += value;
    mean = sum / nsim;
    std::sort(good.begin(), good.end());
    upper = quantile_type7_sorted(good, 0.95);
  }
  return Rcpp::List::create(Rcpp::_["P_heterogeneity"] = p,
    Rcpp::_["MCSE_P_heterogeneity"] = mcse,
    Rcpp::_["null_mean"] = mean, Rcpp::_["null_upper"] = upper,
    Rcpp::_["nsim_successful"] = successful,
    Rcpp::_["exceedance_count"] = successful == nsim && nsim >= 2 ? count : NA_INTEGER,
    Rcpp::_["null"] = keep_null ? Rcpp::RObject(Rcpp::wrap(draws)) : Rcpp::RObject(R_NilValue));
}
}
#endif
