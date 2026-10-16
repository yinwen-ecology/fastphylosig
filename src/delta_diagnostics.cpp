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
#include <numeric>
#include <vector>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

double nan_value() {
  return std::numeric_limits<double>::quiet_NaN();
}

} // namespace

// Fixed production heterogeneity thresholds from .delta_mcmc_warning_summary.
// Only a scalar reason is needed per reference fit; avoid constructing an R
// data frame and affected-traits subset for every local/overall replicate.
// Custom thresholds continue to use the public R summary. No RNG calls here.
// [[Rcpp::export(rng = false)]]
Rcpp::CharacterVector delta_warning_reason_cpp(const double delta,
                                               const Rcpp::List& d) {
  if (!std::isfinite(delta)) return Rcpp::CharacterVector::create(NA_STRING);
  std::string reason;
  auto add = [&](const std::string& message) {
    if (!reason.empty()) reason += "; ";
    reason += message;
  };
  auto value = [&](const char* field) { return Rcpp::as<double>(d[field]); };
  for (const char* field : {"ESS_alpha", "ESS_beta"}) {
    const double x = value(field);
    if (!std::isfinite(x)) add(std::string(field)+" is non-finite");
    else if (x < 20) add(std::string(field)+" < 20");
  }
  for (const char* field : {"split_Rhat_alpha", "split_Rhat_beta"}) {
    const double x = value(field);
    if (!std::isfinite(x)) add(std::string(field)+" is non-finite");
    else if (x > 1.1) add(std::string(field)+" > 1.1");
  }
  const double mcse = value("MCSE_Delta");
  if (!std::isfinite(mcse)) add("MCSE_Delta is non-finite");
  const double ess = value("ESS_Delta"), rhat = value("split_Rhat_Delta");
  if (!std::isfinite(ess) || ess < 20) add("ESS_Delta is non-finite or < 20");
  if (!std::isfinite(rhat) || rhat > 1.1) add("split_Rhat_Delta is non-finite or > 1.1");
  for (const char* quantity : {"alpha", "beta", "Delta"}) {
    for (const char* prefix : {"rank_Rhat_", "ESS_bulk_", "ESS_tail_"}) {
      const std::string field = std::string(prefix)+quantity;
      const double x = value(field.c_str());
      const bool is_rhat = prefix[0]=='r';
      if (!std::isfinite(x) || (is_rhat ? x > 1.01 : x < 200))
        add(field+(is_rhat ? " is non-finite or > 1.01" : " is non-finite or < 200"));
    }
  }
  const double relative = mcse / std::abs(delta);
  if (std::isfinite(relative) && relative > .05) add("relative_MCSE_Delta > 0.05");
  if (reason.empty()) return Rcpp::CharacterVector::create(NA_STRING);
  return Rcpp::CharacterVector::create(reason);
}

namespace {

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
arma::mat chain_mean_covariance(const arma::mat& values) {
  const arma::uword n = values.n_rows;
  if (n < 2) {
    return arma::mat(values.n_cols, values.n_cols, arma::fill::zeros);
  }

  const arma::uword batch = std::max<arma::uword>(1, std::floor(std::sqrt(
    static_cast<double>(n)
  )));
  const arma::uword n_batch = n / batch;
  if (n_batch < 2) {
    return sample_covariance(values) / static_cast<double>(n);
  }

  arma::mat means(n_batch, values.n_cols, arma::fill::zeros);
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

// Vehtari et al. (2021): pooled mid-ranks, folded rank R-hat and
// split-chain bulk/tail ESS. All operations below are deterministic.
arma::mat split_chains(const arma::mat& x) {
  const arma::uword half = x.n_rows / 2;
  arma::mat out(half, 2 * x.n_cols);
  for (arma::uword c = 0; c < x.n_cols; ++c) {
    out.col(2*c) = x.col(c).head(half);
    out.col(2*c+1) = x.col(c).tail(half);
  }
  return out;
}

bool diagnostic_input(const arma::mat& x) {
  return x.n_elem > 1 && x.is_finite() &&
    x.max() - x.min() > std::numeric_limits<double>::epsilon();
}

struct NormalScores {
  arma::uword n;
  std::vector<double> values;
  explicit NormalScores(arma::uword count) : n(count) {
    if (count <= 500000) values.assign(2*count+1,nan_value());
  }
  double get(arma::uword twice_rank) {
    if (!values.empty() && std::isfinite(values[twice_rank])) return values[twice_rank];
    const double rank = static_cast<double>(twice_rank)/2;
    const double value = R::qnorm((rank-.375)/(n+.25),0,1,true,false);
    if (!values.empty()) values[twice_rank] = value;
    return value;
  }
};

std::vector<arma::uword> rank_order(const arma::mat& x) {
  std::vector<arma::uword> order(x.n_elem);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](arma::uword a, arma::uword b) {
    return x[a] < x[b];
  });
  return order;
}

arma::mat rank_normalize(const arma::mat& x,
                         const std::vector<arma::uword>& order,
                         NormalScores& scores) {
  arma::mat out(x.n_rows, x.n_cols);
  for (arma::uword first = 0; first < x.n_elem;) {
    arma::uword last = first + 1;
    while (last < x.n_elem && x[order[last]] == x[order[first]]) ++last;
    const double z = scores.get(first+last+1);
    for (arma::uword i = first; i < last; ++i) out[order[i]] = z;
    first = last;
  }
  return out;
}

// Splitting an odd-length chain omits its middle observation. Merge those
// observations back into the sorted split trace to retain the original full
// trace type-7 quantiles, including ties and middle outliers.
arma::vec sorted_full_trace(const arma::mat& x, const arma::mat& split,
                            const std::vector<arma::uword>& order) {
  arma::vec middle;
  if (x.n_rows % 2) middle = arma::sort(x.row(x.n_rows/2).t());
  arma::vec sorted(x.n_elem);
  arma::uword from_split = 0, from_middle = 0;
  for (arma::uword i=0; i<sorted.n_elem; ++i) {
    if (from_middle>=middle.n_elem ||
        (from_split<order.size() && split[order[from_split]]<=middle[from_middle]))
      sorted[i]=split[order[from_split++]];
    else sorted[i]=middle[from_middle++];
  }
  return sorted;
}

// Absolute distance is monotone on either side of the median. Merging the
// reversed lower half and forward upper half avoids another complete sort.
// Equal distances still receive their original pooled midrank.
std::vector<arma::uword> folded_rank_order(const arma::mat& split,
                                          const arma::mat& folded,
                                          const double median,
                                          const std::vector<arma::uword>& order) {
  arma::uword right=0;
  while (right<order.size() && split[order[right]]<median) ++right;
  arma::uword left=right;
  std::vector<arma::uword> folded_order(order.size());
  for (arma::uword i=0; i<order.size(); ++i) {
    if (right>=order.size() ||
        (left>0 && folded[order[left-1]]<=folded[order[right]]))
      folded_order[i]=order[--left];
    else folded_order[i]=order[right++];
  }
  return folded_order;
}

double basic_rhat(const arma::mat& x) {
  if (x.n_rows < 2 || x.n_cols < 2 || !diagnostic_input(x)) return nan_value();
  const arma::rowvec means = arma::mean(x, 0);
  const double within = arma::mean(arma::var(x, 0, 0));
  if (within <= 0) return nan_value();
  const double between = x.n_rows * arma::var(means);
  return std::sqrt((between / within + x.n_rows - 1) / x.n_rows);
}

// FFT autocovariances avoid the O(iterations * lag) cost of long traces.
// Scaling by n (rather than n-lag) gives Geyer's biased covariance estimate.
double multichain_ess(const arma::mat& x) {
  const arma::uword n = x.n_rows;
  if (n < 3 || !diagnostic_input(x)) return nan_value();
  arma::uword padded = 1;
  while (padded < 2*n) padded *= 2;
  arma::vec autocov(n, arma::fill::zeros);
#if !defined(ARMA_USE_FFTW3)
  // Matrix FFT reuses the same KissFFT preparation across columns. Keep the
  // historical reduction order; dominant matrix buffers stay within 64 MiB.
  // A configured FFTW backend retains scalar calls (its plan policy can differ).
  if (padded <= (64*1024*1024/64) / x.n_cols) {
    arma::mat buffers(padded,x.n_cols,arma::fill::zeros);
    for (arma::uword c=0;c<x.n_cols;++c)
      buffers.col(c).head(n)=x.col(c)-arma::mean(x.col(c));
    const arma::cx_mat spectrum=arma::fft(buffers);
    const arma::cx_mat power=arma::conv_to<arma::cx_mat>::from(arma::square(arma::abs(spectrum)));
    const arma::cx_mat inverse=arma::ifft(power);
    for (arma::uword c=0;c<x.n_cols;++c)
      autocov+=arma::real(inverse.col(c).head(n))/static_cast<double>(n*x.n_cols);
  } else
#endif
  {
    arma::vec buffer(padded, arma::fill::zeros);
    for (arma::uword c = 0; c < x.n_cols; ++c) {
      buffer.zeros();
      buffer.head(n) = x.col(c) - arma::mean(x.col(c));
      const arma::cx_vec spectrum = arma::fft(buffer);
      const arma::cx_vec power = arma::conv_to<arma::cx_vec>::from(arma::square(arma::abs(spectrum)));
      const arma::cx_vec inverse = arma::ifft(power);
      autocov += arma::real(inverse.head(n)) / static_cast<double>(n*x.n_cols);
    }
  }
  const double within = autocov[0] * n / (n-1);
  double variance = autocov[0];
  if (x.n_cols > 1) variance += arma::var(arma::mean(x, 0));
  if (!std::isfinite(variance) || variance <= 0) return nan_value();
  arma::vec rho(n, arma::fill::zeros);
  arma::uword lag = 0;
  double even = 1, odd = 1 - (within-autocov[1]) / variance;
  rho[0] = even; rho[1] = odd;
  while (lag + 5 < n && std::isfinite(even+odd) && even+odd > 0) {
    lag += 2;
    even = 1 - (within-autocov[lag]) / variance;
    odd = 1 - (within-autocov[lag+1]) / variance;
    if (even+odd >= 0) { rho[lag] = even; rho[lag+1] = odd; }
  }
  const arma::uword last = lag;
  if (even > 0) rho[last] = even;
  for (arma::uword t = 2; t+2 <= last; t += 2) {
    const double previous = rho[t-2] + rho[t-1];
    if (rho[t]+rho[t+1] > previous) rho[t] = rho[t+1] = previous/2;
  }
  const double total = static_cast<double>(x.n_elem);
  const double sum = last ? arma::accu(rho.head(last)) : rho[0];
  const double tau = std::max(1/std::log10(total), -1 + 2*sum + rho[last]);
  return total / tau;
}

double quantile_type7(const arma::vec& sorted, const double p) {
  const double index = (sorted.n_elem-1) * p;
  const arma::uword lower = static_cast<arma::uword>(std::floor(index));
  const arma::uword upper = std::min(lower+1, sorted.n_elem-1);
  return sorted[lower] + (index-lower)*(sorted[upper]-sorted[lower]);
}

struct RankDiagnostics { double rhat, bulk, tail; };
RankDiagnostics rank_diagnostics(const arma::mat& x, NormalScores& scores) {
  RankDiagnostics out{nan_value(), nan_value(), nan_value()};
  if (x.n_rows < 4 || !diagnostic_input(x)) return out;
  const arma::mat split = split_chains(x);
  const std::vector<arma::uword> order=rank_order(split);
  const arma::mat ranks = rank_normalize(split,order,scores);
  const arma::vec sorted = sorted_full_trace(x,split,order);
  const double median=quantile_type7(sorted, .5);
  const arma::mat folded = arma::abs(split - median);
  const double bulk_rhat = basic_rhat(ranks);
  const double folded_rhat = basic_rhat(rank_normalize(
    folded,folded_rank_order(split,folded,median,order),scores));
  if (std::isfinite(bulk_rhat) && std::isfinite(folded_rhat))
    out.rhat = std::max(bulk_rhat, folded_rhat);
  out.bulk = multichain_ess(ranks);
  const arma::mat low = arma::conv_to<arma::mat>::from(x <= quantile_type7(sorted, .05));
  const arma::mat high = arma::conv_to<arma::mat>::from(x <= quantile_type7(sorted, .95));
  const double low_ess = multichain_ess(split_chains(low));
  const double high_ess = multichain_ess(split_chains(high));
  if (std::isfinite(low_ess) && std::isfinite(high_ess)) out.tail = std::min(low_ess, high_ess);
  return out;
}

} // namespace

// Retain the historical entry RNG scope, including absent-seed initialization.
// The diagnostic itself makes no random draws.
// [[Rcpp::export]]
Rcpp::List delta_diagnostics_cpp(const arma::mat& alpha_chain,
                                 const arma::mat& beta_chain) {
  if (alpha_chain.n_rows < 1 || alpha_chain.n_cols < 1 ||
      beta_chain.n_rows != alpha_chain.n_rows ||
      beta_chain.n_cols != alpha_chain.n_cols) {
    Rcpp::stop("alpha_chain and beta_chain must have matching non-empty dimensions.");
  }
  if (!alpha_chain.is_finite() || !beta_chain.is_finite() ||
      arma::any(arma::vectorise(alpha_chain) <= 0.0) ||
      arma::any(arma::vectorise(beta_chain) <= 0.0)) {
    Rcpp::stop("alpha_chain and beta_chain must be finite and positive.");
  }

  const arma::uword n = alpha_chain.n_rows;
  const arma::uword chains = alpha_chain.n_cols;
  const arma::vec alpha_all = arma::vectorise(alpha_chain);
  const arma::vec beta_all = arma::vectorise(beta_chain);
  const double alpha_mean = arma::mean(alpha_all);
  const double beta_mean = arma::mean(beta_all);
  const arma::mat ratio = beta_chain / alpha_chain;

  arma::mat all(n * chains, 2);
  all.col(0) = alpha_all;
  all.col(1) = beta_all;
  const arma::mat sample_cov = sample_covariance(all);

  arma::vec chain_ess_alpha(chains, arma::fill::zeros);
  arma::vec chain_ess_beta(chains, arma::fill::zeros);
  arma::vec chain_ess_delta(chains, arma::fill::zeros);
  arma::vec chain_alpha_mean(chains, arma::fill::zeros);
  arma::vec chain_beta_mean(chains, arma::fill::zeros);
  arma::vec chain_alpha_sd(chains, arma::fill::zeros);
  arma::vec chain_beta_sd(chains, arma::fill::zeros);
  arma::mat joint_mean_cov(3, 3, arma::fill::zeros);

  for (arma::uword c = 0; c < chains; ++c) {
    const arma::vec a = alpha_chain.col(c);
    const arma::vec b = beta_chain.col(c);
    chain_ess_alpha(c) = chain_ess(a);
    chain_ess_beta(c) = chain_ess(b);
    chain_ess_delta(c) = chain_ess(ratio.col(c));
    chain_alpha_mean(c) = arma::mean(a);
    chain_beta_mean(c) = arma::mean(b);
    chain_alpha_sd(c) = n > 1 ? arma::stddev(a) : 0.0;
    chain_beta_sd(c) = n > 1 ? arma::stddev(b) : 0.0;
    arma::mat values(n, 3);
    values.col(0) = a;
    values.col(1) = b;
    values.col(2) = ratio.col(c);
    joint_mean_cov += chain_mean_covariance(values) /
      static_cast<double>(chains * chains);
  }

  const double ess_alpha = std::max(
    1.0, std::min(static_cast<double>(n * chains), arma::accu(chain_ess_alpha))
  );
  const double ess_beta = std::max(
    1.0, std::min(static_cast<double>(n * chains), arma::accu(chain_ess_beta))
  );
  const arma::mat mean_cov = joint_mean_cov.submat(0, 0, 1, 1);
  const double ratio_variance = joint_mean_cov(2, 2);
  const double mcse_delta = std::isfinite(ratio_variance) && ratio_variance >= 0.0 ?
    std::sqrt(ratio_variance) : nan_value();

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
  out["ESS_Delta"] = std::max(1.0, std::min(static_cast<double>(n * chains),
    arma::accu(chain_ess_delta)));
  out["split_Rhat_Delta"] = split_rhat(ratio);
  out["split_Rhat_alpha"] = rhat_alpha;
  out["split_Rhat_beta"] = rhat_beta;
  out["Rhat_alpha"] = rhat_alpha;
  out["Rhat_beta"] = rhat_beta;
  out["alpha_beta_cov"] = mean_cov(0, 1);
  out["mean_covariance"] = mean_cov;
  out["sample_alpha_beta_cov"] = sample_cov(0, 1);
  out["MCSE_Delta"] = mcse_delta;
  const arma::mat* quantities[] = {&alpha_chain, &beta_chain, &ratio};
  const char* labels[] = {"alpha", "beta", "Delta"};
  NormalScores scores(2*chains*(n/2));
  for (int i = 0; i < 3; ++i) {
    const RankDiagnostics modern = rank_diagnostics(*quantities[i],scores);
    out[std::string("rank_Rhat_")+labels[i]] = modern.rhat;
    out[std::string("ESS_bulk_")+labels[i]] = modern.bulk;
    out[std::string("ESS_tail_")+labels[i]] = modern.tail;
  }
  out["relative_MCSE_Delta"] = mcse_delta / arma::mean(arma::vectorise(ratio));
  out["n_saved"] = static_cast<double>(n * chains);
  out["n_chains"] = static_cast<int>(chains);
  out["n_iter"] = static_cast<int>(n);
  return out;
}
