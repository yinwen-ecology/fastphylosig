// Diagnostics for the two-chain Delta sampler.
//
// This translation unit is deliberately separate from the Delta estimator.
// The estimator owns the sampler and returns the saved alpha/beta traces;
// this kernel only summarizes those traces.  In particular, no random
// numbers are drawn here, so collecting diagnostics cannot alter an estimate
// or its reproducibility.

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

double nan_value() {
  return std::numeric_limits<double>::quiet_NaN();
}

double chain_ess(const arma::vec& x) {
  const arma::uword n = x.n_elem;
  if (n <= 1) {
    return static_cast<double>(n);
  }

  const double mean = arma::mean(x);
  const arma::vec centered = x - mean;
  const double denominator = arma::dot(centered, centered);
  if (!std::isfinite(denominator) || denominator <= 0.0) {
    return static_cast<double>(n);
  }

  // Geyer's initial-positive-sequence estimate.  A finite lag cap keeps the
  // diagnostic bounded for very long traces while retaining the short-lag
  // dependence that dominates the variance estimate.
  const arma::uword max_lag = std::min<arma::uword>(n - 1, 1000);
  double tau = 1.0;
  for (arma::uword lag = 1; lag <= max_lag; lag += 2) {
    const arma::uword m = n - lag;
    const double acov1 = arma::dot(
      centered.head(m), centered.tail(m)
    ) / denominator;
    double pair = acov1;
    if (lag + 1 <= max_lag) {
      const arma::uword m2 = n - lag - 1;
      const double acov2 = arma::dot(
        centered.head(m2), centered.tail(m2)
      ) / denominator;
      pair += acov2;
    }
    if (!std::isfinite(pair) || pair <= 0.0) {
      break;
    }
    tau += 2.0 * pair;
  }
  if (!std::isfinite(tau) || tau < 1.0) {
    tau = 1.0;
  }
  const double out = static_cast<double>(n) / tau;
  return std::max(1.0, std::min(static_cast<double>(n), out));
}

arma::mat sample_covariance(const arma::mat& values) {
  if (values.n_rows < 2) {
    return arma::mat(values.n_cols, values.n_cols, arma::fill::zeros);
  }
  const arma::rowvec means = arma::mean(values, 0);
  arma::mat centered = values;
  centered.each_row() -= means;
  return (centered.t() * centered) /
    static_cast<double>(values.n_rows - 1);
}

// Estimate the covariance matrix of a chain mean with batch means.  The
// fallback is the usual IID covariance divided by n for traces too short to
// form two batches; callers expose the trace length so such estimates remain
// auditable.
arma::mat chain_mean_covariance(const arma::vec& alpha,
                                const arma::vec& beta) {
  const arma::uword n = alpha.n_elem;
  arma::mat values(n, 2);
  values.col(0) = alpha;
  values.col(1) = beta;
  if (n < 2) {
    return arma::mat(2, 2, arma::fill::zeros);
  }

  const arma::uword batch = std::max<arma::uword>(1, std::floor(std::sqrt(
    static_cast<double>(n)
  )));
  const arma::uword n_batch = n / batch;
  if (n_batch < 2) {
    return sample_covariance(values) / static_cast<double>(n);
  }

  arma::mat means(n_batch, 2, arma::fill::zeros);
  for (arma::uword b = 0; b < n_batch; ++b) {
    const arma::uword first = b * batch;
    const arma::uword last = first + batch;
    means.row(b) = arma::mean(values.rows(first, last - 1), 0);
  }
  const arma::mat between = sample_covariance(means);
  const arma::uword used = n_batch * batch;
  return between * (static_cast<double>(batch) /
                    static_cast<double>(used));
}

double split_rhat(const arma::mat& values) {
  const arma::uword n = values.n_rows;
  const arma::uword chains = values.n_cols;
  if (chains < 2 || n < 4) {
    return nan_value();
  }
  const arma::uword half = n / 2;
  if (half < 2) {
    return nan_value();
  }

  arma::vec means(2 * chains, arma::fill::zeros);
  double within_ss = 0.0;
  for (arma::uword c = 0; c < chains; ++c) {
    for (arma::uword part = 0; part < 2; ++part) {
      const arma::uword first = part == 0 ? 0 : n - half;
      const arma::vec block = values.col(c).subvec(first, first + half - 1);
      const double m = arma::mean(block);
      means(part * chains + c) = m;
      const arma::vec centered = block - m;
      within_ss += arma::dot(centered, centered);
    }
  }

  const double m_total = arma::mean(means);
  const arma::vec centered_means = means - m_total;
  const double between = static_cast<double>(half) * arma::dot(
    centered_means, centered_means
  ) / static_cast<double>(means.n_elem - 1);
  const double within = within_ss /
    static_cast<double>(means.n_elem * (half - 1));
  if (!std::isfinite(within) || within < 0.0 ||
      !std::isfinite(between) || between < 0.0) {
    return nan_value();
  }
  if (within == 0.0) {
    return between == 0.0 ? 1.0 : nan_value();
  }
  const double variance =
    (static_cast<double>(half - 1) / static_cast<double>(half)) * within +
    between / static_cast<double>(half);
  if (!std::isfinite(variance) || variance < 0.0) {
    return nan_value();
  }
  return std::sqrt(std::max(0.0, variance / within));
}

} // namespace

// [[Rcpp::export]]
Rcpp::List delta_diagnostics_cpp(const arma::mat& alpha_chain,
                                 const arma::mat& beta_chain) {
  if (alpha_chain.n_rows < 1 || alpha_chain.n_cols < 1 ||
      beta_chain.n_rows != alpha_chain.n_rows ||
      beta_chain.n_cols != alpha_chain.n_cols) {
    Rcpp::stop("alpha_chain and beta_chain must have matching non-empty dimensions.");
  }
  if (!alpha_chain.is_finite() || !beta_chain.is_finite()) {
    Rcpp::stop("alpha_chain and beta_chain must be finite.");
  }

  const arma::uword n = alpha_chain.n_rows;
  const arma::uword chains = alpha_chain.n_cols;
  const arma::vec alpha_all = arma::vectorise(alpha_chain);
  const arma::vec beta_all = arma::vectorise(beta_chain);
  const double alpha_mean = arma::mean(alpha_all);
  const double beta_mean = arma::mean(beta_all);

  arma::mat all(n * chains, 2);
  all.col(0) = alpha_all;
  all.col(1) = beta_all;
  const arma::mat sample_cov = sample_covariance(all);

  arma::vec chain_ess_alpha(chains, arma::fill::zeros);
  arma::vec chain_ess_beta(chains, arma::fill::zeros);
  arma::vec chain_alpha_mean(chains, arma::fill::zeros);
  arma::vec chain_beta_mean(chains, arma::fill::zeros);
  arma::vec chain_alpha_sd(chains, arma::fill::zeros);
  arma::vec chain_beta_sd(chains, arma::fill::zeros);
  arma::mat mean_cov(2, 2, arma::fill::zeros);

  for (arma::uword c = 0; c < chains; ++c) {
    const arma::vec a = alpha_chain.col(c);
    const arma::vec b = beta_chain.col(c);
    chain_ess_alpha(c) = chain_ess(a);
    chain_ess_beta(c) = chain_ess(b);
    chain_alpha_mean(c) = arma::mean(a);
    chain_beta_mean(c) = arma::mean(b);
    chain_alpha_sd(c) = n > 1 ? arma::stddev(a) : 0.0;
    chain_beta_sd(c) = n > 1 ? arma::stddev(b) : 0.0;
    mean_cov += chain_mean_covariance(a, b) /
      static_cast<double>(chains * chains);
  }

  const double ess_alpha = std::max(
    1.0, std::min(static_cast<double>(n * chains), arma::accu(chain_ess_alpha))
  );
  const double ess_beta = std::max(
    1.0, std::min(static_cast<double>(n * chains), arma::accu(chain_ess_beta))
  );
  double mcse_delta = nan_value();
  if (std::isfinite(alpha_mean) && alpha_mean != 0.0 && mean_cov.is_finite()) {
    const double g_alpha = -beta_mean / (alpha_mean * alpha_mean);
    const double g_beta = 1.0 / alpha_mean;
    const double variance = g_alpha * g_alpha * mean_cov(0, 0) +
      2.0 * g_alpha * g_beta * mean_cov(0, 1) +
      g_beta * g_beta * mean_cov(1, 1);
    if (std::isfinite(variance) && variance >= 0.0) {
      mcse_delta = std::sqrt(variance);
    }
  }

  const double rhat_alpha = split_rhat(alpha_chain);
  const double rhat_beta = split_rhat(beta_chain);
  Rcpp::List out;
  out["alpha_mean"] = alpha_mean;
  out["beta_mean"] = beta_mean;
  out["alpha_sd"] = n * chains > 1 ? arma::stddev(alpha_all) : 0.0;
  out["beta_sd"] = n * chains > 1 ? arma::stddev(beta_all) : 0.0;
  out["chain_alpha_mean"] = chain_alpha_mean;
  out["chain_beta_mean"] = chain_beta_mean;
  out["chain_alpha_sd"] = chain_alpha_sd;
  out["chain_beta_sd"] = chain_beta_sd;
  out["alpha_chain_mean"] = chain_alpha_mean;
  out["beta_chain_mean"] = chain_beta_mean;
  out["alpha_chain_sd"] = chain_alpha_sd;
  out["beta_chain_sd"] = chain_beta_sd;
  out["ESS_alpha"] = ess_alpha;
  out["ESS_beta"] = ess_beta;
  out["split_Rhat_alpha"] = rhat_alpha;
  out["split_Rhat_beta"] = rhat_beta;
  out["Rhat_alpha"] = rhat_alpha;
  out["Rhat_beta"] = rhat_beta;
  out["alpha_beta_cov"] = mean_cov(0, 1);
  out["mean_covariance"] = mean_cov;
  out["sample_alpha_beta_cov"] = sample_cov(0, 1);
  out["MCSE_Delta"] = mcse_delta;
  out["n_saved"] = static_cast<double>(n * chains);
  out["n_chains"] = static_cast<int>(chains);
  out["n_iter"] = static_cast<int>(n);
  return out;
}
