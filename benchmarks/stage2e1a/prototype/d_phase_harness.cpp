// Private Stage 2E1A D phase-decomposition harness.
//
// Audit-only source.  The production streaming translation unit is included
// into this translation unit so the measured Brownian generation and D
// contrast arithmetic are the same implementation used by fast_d().  No
// production binding is changed by this file.

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

#define STAGE2E1A_PRODUCTION_CPP "../../../src/d_stream.cpp"
#include STAGE2E1A_PRODUCTION_CPP
#undef STAGE2E1A_PRODUCTION_CPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace stage2e1a {

typedef std::chrono::steady_clock stage_clock;

double elapsed_seconds(const stage_clock::time_point start,
                       const stage_clock::time_point end) {
  return std::chrono::duration_cast<std::chrono::duration<double> >(
    end - start
  ).count();
}

std::uint64_t mix_bits(const std::uint64_t current, const double value) {
  std::uint64_t bits = 0;
  if (std::isfinite(value)) {
    // This is an audit checksum, not a numerical transformation.  Quantising
    // to 1e12 keeps it stable across Rcpp scalar wrapping while retaining
    // enough information to detect a changed state or D value.
    bits = static_cast<std::uint64_t>(std::llround(
      std::fabs(value) * 1e12
    ));
  } else {
    bits = 0xffffffffffffffffULL;
  }
  return current * 1469598103934665603ULL + bits;
}

std::uint64_t mix_integer(const std::uint64_t current, const int value) {
  return current * 1099511628211ULL +
    static_cast<std::uint64_t>(value + 3);
}

int effective_threads_audit(const int requested) {
  if (requested < 1) Rcpp::stop("n_threads must be positive.");
#ifdef _OPENMP
  return std::max(1, std::min(requested, omp_get_max_threads()));
#else
  return 1;
#endif
}

bool openmp_compiled_audit() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}

int openmp_max_threads_audit() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

void validate_inputs(const Rcpp::IntegerMatrix& edge,
                     const arma::vec& edge_length,
                     const arma::vec& observed,
                     const int n_tip,
                     const int nsim,
                     const double prop_state1,
                     const int simulation_chunk) {
  if (n_tip < 2) Rcpp::stop("n_tip must be at least two.");
  if (nsim < 1) Rcpp::stop("nsim must be positive.");
  if (edge.ncol() != 2 ||
      edge_length.n_elem != static_cast<arma::uword>(edge.nrow())) {
    Rcpp::stop("edge must be a two-column matrix matching edge_length.");
  }
  if (observed.n_elem != static_cast<arma::uword>(n_tip)) {
    Rcpp::stop("observed must have one value per tip.");
  }
  if (!observed.is_finite()) Rcpp::stop("observed must be finite.");
  if (!std::isfinite(prop_state1) || prop_state1 < 0.0 ||
      prop_state1 > 1.0) {
    Rcpp::stop("prop_state1 must be in [0, 1].");
  }
  if (simulation_chunk < 1) Rcpp::stop("simulation_chunk must be positive.");
}

struct BrownianPayload {
  double elapsed_s;
  std::uint64_t checksum;
  int rng_draws;
  Rcpp::NumericMatrix states;
};

// This is the production brownian_tips() routine's call order, kept in the
// same translation unit through d_stream.cpp.  It exposes continuous states
// only for small correctness fixtures; formal timings use checksum mode.
BrownianPayload generate_brownian(const BinaryTree& tree,
                                  const arma::vec& edge_length,
                                  const int nsim,
                                  const bool return_states) {
  if (return_states &&
      static_cast<double>(tree.n_tip) * static_cast<double>(nsim) > 2e6) {
    Rcpp::stop("return_states is limited to small audit fixtures.");
  }
  Rcpp::NumericMatrix states(return_states ? tree.n_tip : 0,
                             return_states ? nsim : 0);
  std::vector<double> node_value(
    static_cast<std::size_t>(tree.max_node + 1), 0.0
  );
  std::vector<double> tip_values(static_cast<std::size_t>(tree.n_tip));
  std::uint64_t checksum = 0;
  int draws = 0;
  Rcpp::RNGScope rng_scope;
  const stage_clock::time_point started = stage_clock::now();
  for (int sim = 0; sim < nsim; ++sim) {
    brownian_tips(tree, edge_length, node_value, tip_values);
    draws += tree.n_edge;
    for (int tip = 0; tip < tree.n_tip; ++tip) {
      const double value = tip_values[static_cast<std::size_t>(tip)];
      checksum = mix_bits(checksum, value);
      if (return_states) states(tip, sim) = value;
    }
  }
  const stage_clock::time_point finished = stage_clock::now();
  BrownianPayload out;
  out.elapsed_s = elapsed_seconds(started, finished);
  out.checksum = checksum;
  out.rng_draws = draws;
  out.states = states;
  return out;
}

struct SortPayload {
  double elapsed_s;
  std::uint64_t checksum;
  std::vector<double> thresholds;
  Rcpp::IntegerMatrix binary;
};

void check_samples(const Rcpp::NumericMatrix& samples, const int n_tip,
                   const int nsim) {
  if (samples.nrow() != n_tip || samples.ncol() != nsim) {
    Rcpp::stop("samples must be n_tip x nsim.");
  }
  for (R_xlen_t i = 0; i < samples.size(); ++i) {
    if (!std::isfinite(samples[i])) {
      Rcpp::stop("samples must contain only finite values.");
    }
  }
}

SortPayload sort_and_threshold(const Rcpp::NumericMatrix& samples,
                               const double prop_state1,
                               const bool threshold,
                               const bool return_binary) {
  const int n_tip = samples.nrow();
  const int nsim = samples.ncol();
  std::vector<double> sorted(static_cast<std::size_t>(n_tip));
  std::vector<double> thresholds(static_cast<std::size_t>(nsim), NA_REAL);
  Rcpp::IntegerMatrix binary(return_binary ? n_tip : 0,
                             return_binary ? nsim : 0);
  std::uint64_t checksum = 0;
  Rcpp::RNGScope rng_scope;
  const stage_clock::time_point started = stage_clock::now();
  for (int sim = 0; sim < nsim; ++sim) {
    for (int tip = 0; tip < n_tip; ++tip) {
      sorted[static_cast<std::size_t>(tip)] = samples(tip, sim);
    }
    std::sort(sorted.begin(), sorted.end());
    for (int tip = 0; tip < n_tip; ++tip) {
      checksum = mix_bits(checksum, sorted[static_cast<std::size_t>(tip)]);
    }
    if (threshold) {
      const double q = fastphylosig::quantile_type7_sorted(
        sorted, prop_state1
      );
      thresholds[static_cast<std::size_t>(sim)] = q;
      checksum = mix_bits(checksum, q);
      if (return_binary) {
        for (int tip = 0; tip < n_tip; ++tip) {
          const int state = samples(tip, sim) < q ? 1 : 0;
          binary(tip, sim) = state;
          checksum = mix_integer(checksum, state);
        }
      }
    }
  }
  const stage_clock::time_point finished = stage_clock::now();
  SortPayload out;
  out.elapsed_s = elapsed_seconds(started, finished);
  out.checksum = checksum;
  out.thresholds = thresholds;
  out.binary = binary;
  return out;
}

Rcpp::IntegerMatrix as_binary(const Rcpp::NumericMatrix& samples,
                              const std::vector<double>& thresholds) {
  const int n_tip = samples.nrow();
  const int nsim = samples.ncol();
  if (static_cast<int>(thresholds.size()) != nsim) {
    Rcpp::stop("thresholds must contain one value per simulation.");
  }
  Rcpp::IntegerMatrix out(n_tip, nsim);
  for (int sim = 0; sim < nsim; ++sim) {
    if (!std::isfinite(thresholds[static_cast<std::size_t>(sim)])) {
      Rcpp::stop("thresholds must be finite.");
    }
    for (int tip = 0; tip < n_tip; ++tip) {
      out(tip, sim) = samples(tip, sim) <
        thresholds[static_cast<std::size_t>(sim)] ? 1 : 0;
    }
  }
  return out;
}

struct Summary {
  double mean_random;
  double mean_brownian;
  double p_random;
  double p_brownian;
  double mcse_random;
  double mcse_brownian;
  int n_successful_random;
  int n_successful_brownian;
  int n_failed_random;
  int n_failed_brownian;
};

Summary summarize(const arma::vec& random_values,
                  const arma::vec& brownian_values,
                  const double observed) {
  Summary out;
  const int n_random = static_cast<int>(random_values.n_elem);
  const int n_brownian = static_cast<int>(brownian_values.n_elem);
  out.n_successful_random = 0;
  out.n_successful_brownian = 0;
  out.n_failed_random = 0;
  out.n_failed_brownian = 0;
  double random_total = 0.0;
  double brownian_total = 0.0;
  int random_less = 0;
  int brownian_greater = 0;
  for (int i = 0; i < n_random; ++i) {
    const double value = random_values(i);
    if (!std::isfinite(value)) {
      ++out.n_failed_random;
    } else {
      ++out.n_successful_random;
      random_total += value;
      if (value < observed) ++random_less;
    }
  }
  for (int i = 0; i < n_brownian; ++i) {
    const double value = brownian_values(i);
    if (!std::isfinite(value)) {
      ++out.n_failed_brownian;
    } else {
      ++out.n_successful_brownian;
      brownian_total += value;
      if (value > observed) ++brownian_greater;
    }
  }
  out.mean_random = out.n_successful_random > 0
    ? random_total / static_cast<double>(out.n_successful_random) : NA_REAL;
  out.mean_brownian = out.n_successful_brownian > 0
    ? brownian_total / static_cast<double>(out.n_successful_brownian) : NA_REAL;
  out.p_random = out.n_successful_random > 0
    ? static_cast<double>(random_less) /
      static_cast<double>(out.n_successful_random) : NA_REAL;
  out.p_brownian = out.n_successful_brownian > 0
    ? static_cast<double>(brownian_greater) /
      static_cast<double>(out.n_successful_brownian) : NA_REAL;
  out.mcse_random = std::isfinite(out.p_random)
    ? std::sqrt(out.p_random * (1.0 - out.p_random) /
                static_cast<double>(out.n_successful_random)) : NA_REAL;
  out.mcse_brownian = std::isfinite(out.p_brownian)
    ? std::sqrt(out.p_brownian * (1.0 - out.p_brownian) /
                static_cast<double>(out.n_successful_brownian)) : NA_REAL;
  (void)n_random;
  (void)n_brownian;
  return out;
}

struct CompletePayload {
  double elapsed_s;
  double random_s;
  double brownian_s;
  double accounting_s;
  std::uint64_t random_checksum;
  std::uint64_t brownian_checksum;
  arma::vec random_values;
  arma::vec brownian_values;
  Summary summary;
  double observed_value;
  int rng_draws;
};

// Actual implementation with the edge matrix retained.  Keeping this
// separate makes the phase-specific wrappers explicit and avoids any hidden
// R-level work in the measured null loop.
CompletePayload complete_pipeline(const BinaryTree& tree,
                                  const Rcpp::IntegerMatrix& edge,
                                  const arma::vec& edge_length,
                                  const arma::vec& observed,
                                  const double prop_state1,
                                  const int nsim,
                                  const int n_threads,
                                  const int simulation_chunk,
                                  const bool return_sim) {
  const int n = tree.n_tip;
  const int chunk_limit = std::min(simulation_chunk, nsim);
  const stage_clock::time_point full_started = stage_clock::now();
  arma::mat observed_matrix(n, 1);
  for (int i = 0; i < n; ++i) observed_matrix(i, 0) = observed(i);
  const arma::vec observed_sum = contrast_sum_batch(
    observed_matrix, tree, edge_length, n_threads
  );
  const double observed_value = observed_sum(0);
  arma::vec random_values(return_sim ? nsim : 0);
  arma::vec brownian_values(return_sim ? nsim : 0);
  std::uint64_t random_checksum = 0;
  std::uint64_t brownian_checksum = 0;
  int rng_draws = 0;
  double random_s = 0.0;
  double brownian_s = 0.0;
  double random_total = 0.0;
  double brownian_total = 0.0;
  int random_success = 0;
  int brownian_success = 0;
  int random_less = 0;
  int brownian_greater = 0;
  std::vector<int> permutation(static_cast<std::size_t>(n));
  // Reuse production-equivalent per-tip workspaces across simulation chunks;
  // this keeps phase timings focused on the same allocation boundary as the
  // streaming implementation.
  std::vector<double> node_value(
    static_cast<std::size_t>(tree.max_node + 1), 0.0
  );
  std::vector<double> tip_values(static_cast<std::size_t>(n));
  std::vector<double> sorted_values(static_cast<std::size_t>(n));
  for (int first = 0; first < nsim; first += chunk_limit) {
    const int count = std::min(nsim, first + chunk_limit) - first;
    const stage_clock::time_point random_phase_started = stage_clock::now();
    arma::mat random_states(n, count, arma::fill::zeros);
    arma::mat brownian_states(n, count, arma::fill::zeros);
    for (int c = 0; c < count; ++c) {
      for (int tip = 0; tip < n; ++tip) permutation[
        static_cast<std::size_t>(tip)
      ] = tip;
      for (int tip = n - 1; tip > 0; --tip) {
        const double u = R::runif(0.0, static_cast<double>(tip + 1));
        int pick = static_cast<int>(std::floor(u));
        if (pick < 0) pick = 0;
        if (pick > tip) pick = tip;
        std::swap(permutation[static_cast<std::size_t>(tip)],
                  permutation[static_cast<std::size_t>(pick)]);
      }
      rng_draws += n - 1;
      for (int tip = 0; tip < n; ++tip) random_states(tip, c) =
        observed(permutation[static_cast<std::size_t>(tip)]);
    }
    const arma::vec random_chunk = contrast_sum_batch(
      random_states, tree, edge_length, n_threads
    );
    const stage_clock::time_point random_finished = stage_clock::now();
    random_s += elapsed_seconds(random_phase_started, random_finished);
    for (int c = 0; c < count; ++c) {
      const int sim = first + c;
      const double value = random_chunk(c);
      random_checksum = mix_bits(random_checksum, value);
      if (std::isfinite(value)) {
        ++random_success;
        random_total += value;
        if (value < observed_value) ++random_less;
      }
      if (return_sim) random_values(sim) = value;
    }

    const stage_clock::time_point brownian_started = stage_clock::now();
    for (int c = 0; c < count; ++c) {
      brownian_tips(tree, edge_length, node_value, tip_values);
      rng_draws += tree.n_edge;
      for (int tip = 0; tip < n; ++tip) {
        brownian_states(tip, c) = tip_values[static_cast<std::size_t>(tip)];
        sorted_values[static_cast<std::size_t>(tip)] =
          tip_values[static_cast<std::size_t>(tip)];
      }
      std::sort(sorted_values.begin(), sorted_values.end());
      const double threshold = fastphylosig::quantile_type7_sorted(
        sorted_values, prop_state1
      );
      for (int tip = 0; tip < n; ++tip) brownian_states(tip, c) =
        brownian_states(tip, c) < threshold ? 1.0 : 0.0;
    }
    const stage_clock::time_point brownian_generated = stage_clock::now();
    const arma::vec brownian_chunk = contrast_sum_batch(
      brownian_states, tree, edge_length, n_threads
    );
    const stage_clock::time_point brownian_finished = stage_clock::now();
    brownian_s += elapsed_seconds(brownian_started, brownian_generated) +
      elapsed_seconds(brownian_generated, brownian_finished);
    for (int c = 0; c < count; ++c) {
      const int sim = first + c;
      const double value = brownian_chunk(c);
      brownian_checksum = mix_bits(brownian_checksum, value);
      if (std::isfinite(value)) {
        ++brownian_success;
        brownian_total += value;
        if (value > observed_value) ++brownian_greater;
      }
      if (return_sim) brownian_values(sim) = value;
    }
  }
  const stage_clock::time_point accounting_started = stage_clock::now();
  Summary summary;
  summary.n_successful_random = random_success;
  summary.n_successful_brownian = brownian_success;
  summary.n_failed_random = nsim - random_success;
  summary.n_failed_brownian = nsim - brownian_success;
  summary.mean_random = random_success > 0
    ? random_total / static_cast<double>(random_success) : NA_REAL;
  summary.mean_brownian = brownian_success > 0
    ? brownian_total / static_cast<double>(brownian_success) : NA_REAL;
  summary.p_random = random_success > 0
    ? static_cast<double>(random_less) / static_cast<double>(random_success)
    : NA_REAL;
  summary.p_brownian = brownian_success > 0
    ? static_cast<double>(brownian_greater) /
      static_cast<double>(brownian_success) : NA_REAL;
  summary.mcse_random = std::isfinite(summary.p_random)
    ? std::sqrt(summary.p_random * (1.0 - summary.p_random) /
                static_cast<double>(random_success)) : NA_REAL;
  summary.mcse_brownian = std::isfinite(summary.p_brownian)
    ? std::sqrt(summary.p_brownian * (1.0 - summary.p_brownian) /
                static_cast<double>(brownian_success)) : NA_REAL;
  const stage_clock::time_point accounting_finished = stage_clock::now();
  CompletePayload out;
  out.elapsed_s = elapsed_seconds(full_started, accounting_finished);
  out.random_s = random_s;
  out.brownian_s = brownian_s;
  out.accounting_s = elapsed_seconds(accounting_started, accounting_finished);
  out.random_checksum = random_checksum;
  out.brownian_checksum = brownian_checksum;
  out.random_values = random_values;
  out.brownian_values = brownian_values;
  out.summary = summary;
  out.observed_value = observed_value;
  out.rng_draws = rng_draws;
  return out;
}

Rcpp::List base_metadata(const BinaryTree& tree, const int n_threads,
                         const int nsim, const std::string& phase) {
  return Rcpp::List::create(
    Rcpp::Named("phase") = phase,
    Rcpp::Named("n") = tree.n_tip,
    Rcpp::Named("n_edge") = tree.n_edge,
    Rcpp::Named("nsim") = nsim,
    Rcpp::Named("n_threads_requested") = n_threads,
    Rcpp::Named("n_threads_effective") = effective_threads_audit(n_threads),
    Rcpp::Named("openmp_compiled") = openmp_compiled_audit(),
    Rcpp::Named("openmp_max_threads") = openmp_max_threads_audit()
  );
}

} // namespace stage2e1a

// [[Rcpp::export]]
Rcpp::List stage2e1a_brownian_generation(
    const Rcpp::IntegerMatrix& edge,
    const arma::vec& edge_length,
    const int n_tip,
    const int nsim,
    const double prop_state1,
    const bool return_states = false,
    const int n_threads = 1) {
  arma::vec observed(static_cast<arma::uword>(n_tip), arma::fill::zeros);
  stage2e1a::validate_inputs(
    edge, edge_length, observed, n_tip, nsim, prop_state1, 1
  );
  const BinaryTree tree = validate_binary_tree(edge, edge_length, n_tip);
  stage2e1a::BrownianPayload payload = stage2e1a::generate_brownian(
    tree, edge_length, nsim, return_states
  );
  Rcpp::List out = stage2e1a::base_metadata(
    tree, n_threads, nsim, "brownian_generation"
  );
  out["elapsed_s"] = payload.elapsed_s;
  out["checksum"] = Rcpp::NumericVector::create(
    static_cast<double>(payload.checksum)
  );
  out["rng_draws"] = payload.rng_draws;
  if (return_states) out["states"] = payload.states;
  return out;
}

// [[Rcpp::export]]
Rcpp::List stage2e1a_brownian_phase(
    const Rcpp::NumericMatrix& samples,
    const double prop_state1,
    const std::string phase,
    const Rcpp::Nullable<Rcpp::NumericVector> thresholds,
    const Rcpp::IntegerMatrix& edge,
    const arma::vec& edge_length,
    const int n_tip,
    const int n_threads = 1,
    const bool return_payload = false) {
  const int n = n_tip > 0 ? n_tip : samples.nrow();
  const int nsim = samples.ncol();
  if (n < 2 || samples.nrow() != n) Rcpp::stop("invalid sample dimensions.");
  if (!std::isfinite(prop_state1) || prop_state1 < 0.0 ||
      prop_state1 > 1.0) Rcpp::stop("prop_state1 must be in [0, 1].");
  stage2e1a::check_samples(samples, n, nsim);
  std::vector<double> threshold_values;
  if (thresholds.isNotNull()) {
    Rcpp::NumericVector supplied(thresholds);
    if (supplied.size() != nsim) Rcpp::stop("thresholds must match nsim.");
    threshold_values.assign(supplied.begin(), supplied.end());
  }
  Rcpp::List out;
  if (phase == "sort") {
    stage2e1a::SortPayload sorted = stage2e1a::sort_and_threshold(
      samples, prop_state1, false, false
    );
    out = Rcpp::List::create(
      Rcpp::Named("phase") = phase,
      Rcpp::Named("elapsed_s") = sorted.elapsed_s,
      Rcpp::Named("checksum") = static_cast<double>(sorted.checksum),
      Rcpp::Named("n") = n, Rcpp::Named("nsim") = nsim
    );
  } else if (phase == "sort_threshold") {
    stage2e1a::SortPayload sorted = stage2e1a::sort_and_threshold(
      samples, prop_state1, true, false
    );
    Rcpp::NumericVector q(nsim);
    for (int i = 0; i < nsim; ++i) q(i) = sorted.thresholds[
      static_cast<std::size_t>(i)
    ];
    out = Rcpp::List::create(
      Rcpp::Named("phase") = phase,
      Rcpp::Named("elapsed_s") = sorted.elapsed_s,
      Rcpp::Named("checksum") = static_cast<double>(sorted.checksum),
      Rcpp::Named("thresholds") = q,
      Rcpp::Named("n") = n, Rcpp::Named("nsim") = nsim
    );
  } else if (phase == "binary_d") {
    if (threshold_values.empty()) {
      Rcpp::stop("binary_d requires thresholds from sort_threshold.");
    }
    const Rcpp::IntegerMatrix binary = stage2e1a::as_binary(
      samples, threshold_values
    );
    if (edge.ncol() != 2 || edge_length.n_elem !=
        static_cast<arma::uword>(edge.nrow())) {
      Rcpp::stop("binary_d requires edge and edge_length.");
    }
    const BinaryTree d_tree = validate_binary_tree(edge, edge_length, n);
    const stage2e1a::stage_clock::time_point started =
      stage2e1a::stage_clock::now();
    arma::mat converted(n, nsim);
    for (int j = 0; j < nsim; ++j) {
      for (int i = 0; i < n; ++i) converted(i, j) = binary(i, j);
    }
    const arma::vec values = contrast_sum_batch(
      converted, d_tree, edge_length, n_threads
    );
    const stage2e1a::stage_clock::time_point finished =
      stage2e1a::stage_clock::now();
    std::uint64_t checksum = 0;
    for (arma::uword i = 0; i < values.n_elem; ++i) {
      checksum = stage2e1a::mix_bits(checksum, values(i));
    }
    out = Rcpp::List::create(
      Rcpp::Named("phase") = phase,
      Rcpp::Named("elapsed_s") = stage2e1a::elapsed_seconds(started, finished),
      Rcpp::Named("checksum") = static_cast<double>(checksum),
      Rcpp::Named("n") = n, Rcpp::Named("nsim") = nsim
    );
    if (return_payload) {
      out["binary"] = binary;
      out["D"] = values;
    }
  } else {
    Rcpp::stop("phase must be sort, sort_threshold, or binary_d.");
  }
  out["n_threads_requested"] = n_threads;
  out["n_threads_effective"] = stage2e1a::effective_threads_audit(n_threads);
  out["openmp_compiled"] = stage2e1a::openmp_compiled_audit();
  out["openmp_max_threads"] = stage2e1a::openmp_max_threads_audit();
  return out;
}

// [[Rcpp::export]]
Rcpp::List stage2e1a_random_null(
    const Rcpp::IntegerMatrix& edge,
    const arma::vec& edge_length,
    const arma::vec& observed,
    const int n_tip,
    const int nsim,
    const int n_threads = 1,
    const int simulation_chunk = 128,
    const bool return_sim = false) {
  stage2e1a::validate_inputs(
    edge, edge_length, observed, n_tip, nsim, 0.5, simulation_chunk
  );
  const BinaryTree tree = validate_binary_tree(edge, edge_length, n_tip);
  const int chunk_limit = std::min(simulation_chunk, nsim);
  arma::vec values(return_sim ? nsim : 0);
  std::uint64_t checksum = 0;
  std::uint64_t draws = 0;
  std::vector<int> permutation(static_cast<std::size_t>(n_tip));
  Rcpp::RNGScope rng_scope;
  const stage2e1a::stage_clock::time_point started =
    stage2e1a::stage_clock::now();
  for (int first = 0; first < nsim; first += chunk_limit) {
    const int count = std::min(nsim, first + chunk_limit) - first;
    arma::mat states(n_tip, count, arma::fill::zeros);
    for (int c = 0; c < count; ++c) {
      for (int tip = 0; tip < n_tip; ++tip) permutation[
        static_cast<std::size_t>(tip)
      ] = tip;
      for (int tip = n_tip - 1; tip > 0; --tip) {
        const double u = R::runif(0.0, static_cast<double>(tip + 1));
        int pick = static_cast<int>(std::floor(u));
        if (pick < 0) pick = 0;
        if (pick > tip) pick = tip;
        std::swap(permutation[static_cast<std::size_t>(tip)],
                  permutation[static_cast<std::size_t>(pick)]);
      }
      draws += static_cast<std::uint64_t>(n_tip - 1);
      for (int tip = 0; tip < n_tip; ++tip) states(tip, c) =
        observed(permutation[static_cast<std::size_t>(tip)]);
    }
    const arma::vec chunk = contrast_sum_batch(
      states, tree, edge_length, n_threads
    );
    for (int c = 0; c < count; ++c) {
      const int sim = first + c;
      checksum = stage2e1a::mix_bits(checksum, chunk(c));
      if (return_sim) values(sim) = chunk(c);
    }
  }
  const stage2e1a::stage_clock::time_point finished =
    stage2e1a::stage_clock::now();
  Rcpp::List out = stage2e1a::base_metadata(
    tree, n_threads, nsim, "random_null"
  );
  out["elapsed_s"] = stage2e1a::elapsed_seconds(started, finished);
  out["checksum"] = static_cast<double>(checksum);
  out["rng_draws"] = static_cast<double>(draws);
  if (return_sim) out["D"] = values;
  return out;
}

// [[Rcpp::export]]
Rcpp::List stage2e1a_complete_d(
    const Rcpp::IntegerMatrix& edge,
    const arma::vec& edge_length,
    const arma::vec& observed,
    const int n_tip,
    const int nsim,
    const double prop_state1,
    const int n_threads = 1,
    const int simulation_chunk = 128,
    const bool return_sim = false) {
  stage2e1a::validate_inputs(
    edge, edge_length, observed, n_tip, nsim, prop_state1, simulation_chunk
  );
  const BinaryTree tree = validate_binary_tree(edge, edge_length, n_tip);
  Rcpp::RNGScope rng_scope;
  stage2e1a::CompletePayload payload = stage2e1a::complete_pipeline(
    tree, edge, edge_length, observed, prop_state1, nsim, n_threads,
    simulation_chunk, return_sim
  );
  Rcpp::List out = stage2e1a::base_metadata(
    tree, n_threads, nsim, "complete_d"
  );
  out["elapsed_s"] = payload.elapsed_s;
  out["random_s"] = payload.random_s;
  out["brownian_s"] = payload.brownian_s;
  out["accounting_s"] = payload.accounting_s;
  out["random_checksum"] = static_cast<double>(payload.random_checksum);
  out["brownian_checksum"] = static_cast<double>(payload.brownian_checksum);
  out["rng_draws"] = payload.rng_draws;
  out["mean_random"] = payload.summary.mean_random;
  out["mean_brownian"] = payload.summary.mean_brownian;
  out["p_random"] = payload.summary.p_random;
  out["p_brownian"] = payload.summary.p_brownian;
  out["mcse_random"] = payload.summary.mcse_random;
  out["mcse_brownian"] = payload.summary.mcse_brownian;
  out["nsim_requested"] = nsim;
  out["nsim_successful_random"] = payload.summary.n_successful_random;
  out["nsim_successful_brownian"] = payload.summary.n_successful_brownian;
  out["nsim_failed_random"] = payload.summary.n_failed_random;
  out["nsim_failed_brownian"] = payload.summary.n_failed_brownian;
  out["observed"] = payload.observed_value;
  if (return_sim) {
    out["random_D"] = payload.random_values;
    out["brownian_D"] = payload.brownian_values;
  }
  return out;
}
