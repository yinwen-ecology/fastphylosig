// Private Stage 2D2 Candidate C prototype for the K permutation/null engine.
//
// The production translation unit is included so the oracle and Candidate C
// use the same kperm::Tree, kperm::Cache, tip accessor, and Fisher-Yates
// implementation.  This file is audit-only: it is not registered in the
// package DLL and does not change R/, src/, tests/, or the public API.

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

// Hide the literal .cpp include from Rcpp::sourceCpp.  Otherwise sourceCpp
// compiles the included file as a second object and the same-TU oracle cannot
// be linked because every production symbol is defined twice.
#define STAGE2D2_PRODUCTION_CPP "../../../src/k_permutation.cpp"
#include STAGE2D2_PRODUCTION_CPP
#undef STAGE2D2_PRODUCTION_CPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace stage2d2_candidate_c {

// Candidate C's only numerical change is to hoist this exact production
// expression out of compute_one's two downward passes.  Tips retain alpha=1.
std::vector<double> precompute_alpha(const kperm::Tree& tree,
                                     const kperm::Cache& cache) {
  std::vector<double> alpha(static_cast<std::size_t>(tree.n_total), 1.0);
  for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
    const int node = tree.preorder[ii];
    alpha[static_cast<std::size_t>(node)] = node < tree.n_tip ? 1.0 :
      tree.branch[static_cast<std::size_t>(node)] *
      cache.outgoing[static_cast<std::size_t>(node)];
  }
  return alpha;
}

// This is a line-for-line numerical copy of kperm::compute_one(), with only
// the repeated alpha expression replaced by the call-level alpha vector.
void compute_one_candidate_c(const kperm::Tree& tree,
                             const kperm::Cache& cache,
                             const std::vector<double>& alpha,
                             const double* x, const int ncol,
                             const std::vector<int>& perm,
                             const int trait_chunk,
                             std::vector<double>& out) {
  const int n = tree.n_tip;
  const int p = ncol;
  out.assign(static_cast<std::size_t>(p),
             std::numeric_limits<double>::quiet_NaN());
  const int chunk_limit = std::max(1, std::min(trait_chunk, p));
  std::vector<double> message;
  std::vector<double> state;
  std::vector<double> baseline;
  // GLS offsets are kept relative to the per-trait baseline.  Never rebuild
  // the absolute mean b + delta for the residual pass: doing so would lose
  // the low-order trait variation when b is around 1e12--1e15.
  std::vector<double> delta;
  for (int col0 = 0; col0 < p; col0 += chunk_limit) {
    const int chunk = std::min(chunk_limit, p - col0);
    const std::size_t cells = static_cast<std::size_t>(tree.n_total) *
      static_cast<std::size_t>(chunk);
    message.assign(cells, 0.0);
    state.assign(cells, 0.0);
    baseline.assign(static_cast<std::size_t>(chunk), 0.0);
    delta.assign(static_cast<std::size_t>(chunk), 0.0);

    for (int j = 0; j < chunk; ++j) {
      baseline[static_cast<std::size_t>(j)] =
        kperm::tip_value(x, n, perm, 0, col0 + j);
    }
    // Upward Gaussian messages for x - baseline.
    for (std::size_t ii = 0; ii < tree.postorder.size(); ++ii) {
      const int node = tree.postorder[ii];
      const std::size_t base = static_cast<std::size_t>(node) *
        static_cast<std::size_t>(chunk);
      if (node < n) {
        for (int j = 0; j < chunk; ++j) {
          message[base + static_cast<std::size_t>(j)] =
            kperm::tip_value(x, n, perm, node, col0 + j) -
            baseline[static_cast<std::size_t>(j)];
        }
      } else {
        const int begin = tree.child_ptr[static_cast<std::size_t>(node)];
        const int end = tree.child_ptr[static_cast<std::size_t>(node + 1)];
        const double s = cache.aggregate[static_cast<std::size_t>(node)];
        for (int j = 0; j < chunk; ++j) {
          long double weighted = 0.0L;
          for (int k = begin; k < end; ++k) {
            const int child = tree.children[static_cast<std::size_t>(k)];
            weighted += static_cast<long double>(cache.outgoing[
              static_cast<std::size_t>(child)]) *
              static_cast<long double>(message[
                static_cast<std::size_t>(child) *
                static_cast<std::size_t>(chunk) + static_cast<std::size_t>(j)
              ]);
          }
          message[base + static_cast<std::size_t>(j)] =
            static_cast<double>(weighted / static_cast<long double>(s));
        }
      }
    }
    // Downward conditional states and baseline-relative GLS offset delta.
    for (int j = 0; j < chunk; ++j) {
      state[static_cast<std::size_t>(tree.root) *
            static_cast<std::size_t>(chunk) + static_cast<std::size_t>(j)] = 0.0;
    }
    for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
      const int node = tree.preorder[ii];
      const int par = tree.parent[static_cast<std::size_t>(node)];
      const double alpha_value = alpha[static_cast<std::size_t>(node)];
      const std::size_t base = static_cast<std::size_t>(node) *
        static_cast<std::size_t>(chunk);
      const std::size_t pbase = static_cast<std::size_t>(par) *
        static_cast<std::size_t>(chunk);
      for (int j = 0; j < chunk; ++j) {
        const double parent_state = pbase == base ? 0.0 :
          state[pbase + static_cast<std::size_t>(j)];
        state[base + static_cast<std::size_t>(j)] = parent_state +
          alpha_value * (message[base + static_cast<std::size_t>(j)] - parent_state);
      }
    }
    for (int j = 0; j < chunk; ++j) {
      long double qlinear = 0.0L;
      for (int tip = 0; tip < n; ++tip) {
        const int par = tree.parent[static_cast<std::size_t>(tip)];
        const double parent_state = state[
          static_cast<std::size_t>(par) * static_cast<std::size_t>(chunk) +
          static_cast<std::size_t>(j)
        ];
        qlinear += (static_cast<long double>(
          kperm::tip_value(x, n, perm, tip, col0 + j) -
          baseline[static_cast<std::size_t>(j)] - parent_state
        )) / static_cast<long double>(tree.branch[static_cast<std::size_t>(tip)]);
      }
      delta[static_cast<std::size_t>(j)] =
        static_cast<double>(qlinear / static_cast<long double>(cache.sum_inv));
    }

    // Repeat on (x - baseline) - delta to obtain the numerator and precision
    // energy.  This is algebraically x - GLS_mean, but is stable for traits
    // carrying a very large common offset.
    for (std::size_t ii = 0; ii < tree.postorder.size(); ++ii) {
      const int node = tree.postorder[ii];
      const std::size_t base = static_cast<std::size_t>(node) *
        static_cast<std::size_t>(chunk);
      if (node < n) {
        for (int j = 0; j < chunk; ++j) {
          message[base + static_cast<std::size_t>(j)] =
            kperm::tip_value(x, n, perm, node, col0 + j) -
            baseline[static_cast<std::size_t>(j)] -
            delta[static_cast<std::size_t>(j)];
        }
      } else {
        const int begin = tree.child_ptr[static_cast<std::size_t>(node)];
        const int end = tree.child_ptr[static_cast<std::size_t>(node + 1)];
        const double s = cache.aggregate[static_cast<std::size_t>(node)];
        for (int j = 0; j < chunk; ++j) {
          long double weighted = 0.0L;
          for (int k = begin; k < end; ++k) {
            const int child = tree.children[static_cast<std::size_t>(k)];
            weighted += static_cast<long double>(cache.outgoing[
              static_cast<std::size_t>(child)]) *
              static_cast<long double>(message[
                static_cast<std::size_t>(child) *
                static_cast<std::size_t>(chunk) + static_cast<std::size_t>(j)
              ]);
          }
          message[base + static_cast<std::size_t>(j)] =
            static_cast<double>(weighted / static_cast<long double>(s));
        }
      }
    }
    for (int j = 0; j < chunk; ++j) {
      state[static_cast<std::size_t>(tree.root) *
            static_cast<std::size_t>(chunk) + static_cast<std::size_t>(j)] = 0.0;
    }
    for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
      const int node = tree.preorder[ii];
      const int par = tree.parent[static_cast<std::size_t>(node)];
      const double alpha_value = alpha[static_cast<std::size_t>(node)];
      const std::size_t base = static_cast<std::size_t>(node) *
        static_cast<std::size_t>(chunk);
      const std::size_t pbase = static_cast<std::size_t>(par) *
        static_cast<std::size_t>(chunk);
      for (int j = 0; j < chunk; ++j) {
        const double parent_state = state[pbase + static_cast<std::size_t>(j)];
        state[base + static_cast<std::size_t>(j)] = parent_state +
          alpha_value * (message[base + static_cast<std::size_t>(j)] - parent_state);
      }
    }
    for (int j = 0; j < chunk; ++j) {
      long double numerator = 0.0L;
      long double denominator = 0.0L;
      for (int tip = 0; tip < n; ++tip) {
        const long double y = static_cast<long double>(
          kperm::tip_value(x, n, perm, tip, col0 + j) -
          baseline[static_cast<std::size_t>(j)] -
          delta[static_cast<std::size_t>(j)]
        );
        numerator += y * y;
      }
      for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
        const int node = tree.preorder[ii];
        const int par = tree.parent[static_cast<std::size_t>(node)];
        const long double parent_state = state[
          static_cast<std::size_t>(par) * static_cast<std::size_t>(chunk) +
          static_cast<std::size_t>(j)
        ];
        const long double child_state = node < n
          ? static_cast<long double>(
              kperm::tip_value(x, n, perm, node, col0 + j) -
              baseline[static_cast<std::size_t>(j)] -
              delta[static_cast<std::size_t>(j)]
            )
          : static_cast<long double>(state[
              static_cast<std::size_t>(node) * static_cast<std::size_t>(chunk) +
              static_cast<std::size_t>(j)
            ]);
        const long double diff = child_state - parent_state;
        denominator += diff * diff /
          static_cast<long double>(tree.branch[static_cast<std::size_t>(node)]);
      }
      const double den = static_cast<double>(denominator);
      const double num = static_cast<double>(numerator);
      out[static_cast<std::size_t>(col0 + j)] =
        den > 0.0 && std::isfinite(den) && std::isfinite(num)
        ? (num / den) / cache.normalization
        : std::numeric_limits<double>::quiet_NaN();
    }
  }
}

void validate_supplied(const Rcpp::IntegerMatrix& perm_matrix,
                       const int nsim, const int n) {
  if (perm_matrix.nrow() != nsim || perm_matrix.ncol() != n) {
    Rcpp::stop("permutations must have nsim rows and one column per tip.");
  }
  std::vector<int> seen(static_cast<std::size_t>(n), -1);
  for (int i = 0; i < nsim; ++i) {
    for (int r = 0; r < n; ++r) {
      const int idx = perm_matrix(i, r);
      if (idx < 1 || idx > n || seen[static_cast<std::size_t>(idx - 1)] == i) {
        Rcpp::stop("each permutation row must contain 1:n exactly once.");
      }
      seen[static_cast<std::size_t>(idx - 1)] = i;
    }
  }
}

struct EngineResult {
  Rcpp::NumericVector K;
  Rcpp::NumericVector P;
  Rcpp::NumericVector exceedance;
  Rcpp::NumericVector MCSE;
  Rcpp::IntegerVector successful;
  Rcpp::IntegerVector failed;
  Rcpp::LogicalVector success;
  Rcpp::Nullable<Rcpp::NumericMatrix> sim;
  std::size_t alpha_nodes;
  std::size_t alpha_bytes;
};

template <typename Compute>
EngineResult run_engine(const Rcpp::List& compiled_tree,
                        const Rcpp::NumericMatrix& X, const int nsim,
                        SEXP permutations, const int trait_chunk,
                        const bool store_sim, const bool include_observed,
                        const int n_threads, const int simulation_chunk,
                        Compute compute, const bool use_alpha) {
  if (nsim < 1) Rcpp::stop("nsim must be a positive integer.");
  if (trait_chunk < 1) Rcpp::stop("trait_chunk must be positive.");
  if (simulation_chunk < 1) {
    Rcpp::stop("simulation_chunk must be positive.");
  }
  const kperm::Tree tree = kperm::parse_tree(compiled_tree);
  const int n = tree.n_tip;
  const int p = X.ncol();
  if (X.nrow() != n || p < 1) {
    Rcpp::stop("X must have one row per compiled tree tip and at least one trait.");
  }
  const double* x = REAL(X);
  for (R_xlen_t i = 0; i < X.size(); ++i) {
    if (!std::isfinite(x[i])) Rcpp::stop("X must contain only finite values.");
  }
  const kperm::Cache cache = kperm::build_cache(tree);
  std::vector<double> alpha;
  if (use_alpha) alpha = precompute_alpha(tree, cache);

  const bool supplied = permutations != R_NilValue;
  Rcpp::IntegerMatrix perm_matrix;
  if (supplied) {
    perm_matrix = Rcpp::IntegerMatrix(permutations);
    validate_supplied(perm_matrix, nsim, n);
  }

  Rcpp::NumericVector observed(p);
  std::vector<int> identity(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) identity[static_cast<std::size_t>(i)] = i;
  std::vector<double> observed_vec;
  compute(tree, cache, alpha, x, p, identity, trait_chunk, observed_vec);
  for (int j = 0; j < p; ++j) observed[j] = observed_vec[static_cast<std::size_t>(j)];
  std::vector<char> valid_observed(static_cast<std::size_t>(p), 1);
  for (int j = 0; j < p; ++j) {
    valid_observed[static_cast<std::size_t>(j)] =
      std::isfinite(observed[j]) ? 1 : 0;
  }

  Rcpp::NumericVector exceedance(p, 0.0);
  Rcpp::NumericMatrix sim;
  if (store_sim) sim = Rcpp::NumericMatrix(nsim, p);
  double* sim_ptr = store_sim ? REAL(sim) : NULL;
  const int threads_eff = kperm::effective_threads(n_threads);
  const int chunk_limit = std::min(simulation_chunk, nsim);

  Rcpp::RNGScope rng_scope;
  if (supplied && threads_eff > 1) {
#ifdef _OPENMP
    std::vector<std::vector<double> > local(
      static_cast<std::size_t>(threads_eff),
      std::vector<double>(static_cast<std::size_t>(p), 0.0)
    );
#pragma omp parallel num_threads(threads_eff)
    {
      const int tid = omp_get_thread_num();
      std::vector<int> perm(static_cast<std::size_t>(n));
      std::vector<double> kval;
#pragma omp for schedule(static)
      for (int i = 0; i < nsim; ++i) {
        for (int r = 0; r < n; ++r) perm[static_cast<std::size_t>(r)] =
          perm_matrix(i, r) - 1;
        compute(tree, cache, alpha, x, p, perm, trait_chunk, kval);
        for (int j = 0; j < p; ++j) {
          const double value = kval[static_cast<std::size_t>(j)];
          if (store_sim) sim_ptr[static_cast<std::size_t>(i) +
                                  static_cast<std::size_t>(nsim) *
                                  static_cast<std::size_t>(j)] = value;
          if (valid_observed[static_cast<std::size_t>(j)] &&
              fastphylosig::inclusive_upper_tail(value, observed[j]))
            local[static_cast<std::size_t>(tid)][static_cast<std::size_t>(j)] += 1.0;
        }
      }
    }
    for (int t = 0; t < threads_eff; ++t) {
      for (int j = 0; j < p; ++j) exceedance[j] +=
        local[static_cast<std::size_t>(t)][static_cast<std::size_t>(j)];
    }
#else
    (void)threads_eff;
#endif
  } else if (supplied) {
    std::vector<int> perm(static_cast<std::size_t>(n));
    std::vector<double> kval;
    for (int i = 0; i < nsim; ++i) {
      for (int r = 0; r < n; ++r) perm[static_cast<std::size_t>(r)] =
        perm_matrix(i, r) - 1;
      compute(tree, cache, alpha, x, p, perm, trait_chunk, kval);
      for (int j = 0; j < p; ++j) {
        const double value = kval[static_cast<std::size_t>(j)];
        if (store_sim) sim_ptr[static_cast<std::size_t>(i) +
                                static_cast<std::size_t>(nsim) *
                                static_cast<std::size_t>(j)] = value;
        if (valid_observed[static_cast<std::size_t>(j)] &&
            fastphylosig::inclusive_upper_tail(value, observed[j]))
          exceedance[j] += 1.0;
      }
    }
  } else {
    std::vector<int> chunk_perms(
      static_cast<std::size_t>(chunk_limit) * static_cast<std::size_t>(n)
    );
    for (int first = 0; first < nsim; first += chunk_limit) {
      const int count = std::min(chunk_limit, nsim - first);
      for (int local_i = 0; local_i < count; ++local_i) {
        std::vector<int> perm(static_cast<std::size_t>(n));
        if (include_observed && first + local_i == 0) {
          perm = identity;
        } else {
          perm = identity;
          kperm::fisher_yates(perm);
        }
        std::copy(perm.begin(), perm.end(), chunk_perms.begin() +
                  static_cast<std::size_t>(local_i) *
                  static_cast<std::size_t>(n));
      }
      if (threads_eff > 1) {
#ifdef _OPENMP
        std::vector<std::vector<double> > local_counts(
          static_cast<std::size_t>(threads_eff),
          std::vector<double>(static_cast<std::size_t>(p), 0.0)
        );
#pragma omp parallel num_threads(threads_eff)
        {
          const int tid = omp_get_thread_num();
          std::vector<int> perm(static_cast<std::size_t>(n));
          std::vector<double> kval;
#pragma omp for schedule(static)
          for (int local_i = 0; local_i < count; ++local_i) {
            std::copy(
              chunk_perms.begin() + static_cast<std::size_t>(local_i) *
                static_cast<std::size_t>(n),
              chunk_perms.begin() + static_cast<std::size_t>(local_i + 1) *
                static_cast<std::size_t>(n), perm.begin()
            );
            compute(tree, cache, alpha, x, p, perm, trait_chunk, kval);
            const int global_i = first + local_i;
            for (int j = 0; j < p; ++j) {
              const double value = kval[static_cast<std::size_t>(j)];
              if (store_sim) sim_ptr[static_cast<std::size_t>(global_i) +
                                      static_cast<std::size_t>(nsim) *
                                      static_cast<std::size_t>(j)] = value;
              if (valid_observed[static_cast<std::size_t>(j)] &&
                  fastphylosig::inclusive_upper_tail(value, observed[j]))
                local_counts[static_cast<std::size_t>(tid)][
                  static_cast<std::size_t>(j)] += 1.0;
            }
          }
        }
        for (int t = 0; t < threads_eff; ++t) {
          for (int j = 0; j < p; ++j) exceedance[j] += local_counts[
            static_cast<std::size_t>(t)][static_cast<std::size_t>(j)];
        }
#else
        (void)threads_eff;
#endif
      } else {
        std::vector<int> perm(static_cast<std::size_t>(n));
        std::vector<double> kval;
        for (int local_i = 0; local_i < count; ++local_i) {
          std::copy(
            chunk_perms.begin() + static_cast<std::size_t>(local_i) *
              static_cast<std::size_t>(n),
            chunk_perms.begin() + static_cast<std::size_t>(local_i + 1) *
              static_cast<std::size_t>(n), perm.begin()
          );
          compute(tree, cache, alpha, x, p, perm, trait_chunk, kval);
          const int global_i = first + local_i;
          for (int j = 0; j < p; ++j) {
            const double value = kval[static_cast<std::size_t>(j)];
            if (store_sim) sim_ptr[static_cast<std::size_t>(global_i) +
                                    static_cast<std::size_t>(nsim) *
                                    static_cast<std::size_t>(j)] = value;
            if (valid_observed[static_cast<std::size_t>(j)] &&
                fastphylosig::inclusive_upper_tail(value, observed[j]))
              exceedance[j] += 1.0;
          }
        }
      }
    }
  }

  Rcpp::NumericVector p_value(p), mcse(p);
  Rcpp::IntegerVector nsim_successful(p, nsim);
  Rcpp::IntegerVector nsim_failed(p, 0);
  const int n_random = supplied ? nsim :
    (include_observed ? std::max(0, nsim - 1) : nsim);
  for (int j = 0; j < p; ++j) {
    if (!valid_observed[static_cast<std::size_t>(j)]) {
      p_value[j] = NA_REAL;
      mcse[j] = NA_REAL;
      nsim_successful[j] = 0;
      nsim_failed[j] = nsim;
      continue;
    }
    p_value[j] = exceedance[j] / static_cast<double>(nsim);
    if (!supplied && include_observed && n_random > 0) {
      const double q = std::max(0.0, std::min(1.0,
        (exceedance[j] - 1.0) / static_cast<double>(n_random)));
      mcse[j] = std::sqrt(q * (1.0 - q) * static_cast<double>(n_random)) /
        static_cast<double>(nsim);
    } else if (n_random > 0) {
      const double q = p_value[j];
      mcse[j] = std::sqrt(q * (1.0 - q) /
                          static_cast<double>(n_random));
    } else {
      mcse[j] = NA_REAL;
    }
  }

  Rcpp::Nullable<Rcpp::NumericMatrix> sim_value(R_NilValue);
  if (store_sim) sim_value = Rcpp::Nullable<Rcpp::NumericMatrix>(sim);
  EngineResult result = {
    observed, p_value, exceedance, mcse, nsim_successful,
    nsim_failed, Rcpp::LogicalVector(nsim, true), sim_value,
    use_alpha ? static_cast<std::size_t>(alpha.size()) : 0,
    use_alpha ? static_cast<std::size_t>(alpha.size()) * sizeof(double) : 0
  };
  return result;
}

struct OracleCompute {
  void operator()(const kperm::Tree& tree, const kperm::Cache& cache,
                  const std::vector<double>&,
                  const double* x, const int ncol,
                  const std::vector<int>& perm, const int trait_chunk,
                  std::vector<double>& out) const {
    kperm::compute_one(tree, cache, x, ncol, perm, trait_chunk, out);
  }
};

struct CandidateCCompute {
  void operator()(const kperm::Tree& tree, const kperm::Cache& cache,
                  const std::vector<double>& alpha,
                  const double* x, const int ncol,
                  const std::vector<int>& perm, const int trait_chunk,
                  std::vector<double>& out) const {
    compute_one_candidate_c(tree, cache, alpha, x, ncol, perm, trait_chunk, out);
  }
};

}  // namespace stage2d2_candidate_c

// [[Rcpp::export]]
Rcpp::List stage2d2_k_permutation_prototype(
    const Rcpp::List& compiled_tree,
    const Rcpp::NumericMatrix& X,
    const int nsim = 999,
    SEXP permutations = R_NilValue,
    const int trait_chunk = 64,
    const bool return_sim = true,
    const bool include_observed = true,
    const int n_threads = 1,
    const int simulation_chunk = 128,
    const std::string mode = "oracle",
    const bool return_ordered = true) {
  using namespace stage2d2_candidate_c;
  if (mode != "oracle" && mode != "C") {
    Rcpp::stop("Stage 2D2 prototype supports only mode='oracle' or mode='C'.");
  }
  const bool use_alpha = mode == "C";
  const bool store_sim = return_sim || return_ordered;
  EngineResult result = use_alpha
    ? run_engine(compiled_tree, X, nsim, permutations, trait_chunk, store_sim,
                 include_observed, n_threads, simulation_chunk,
                 CandidateCCompute(), true)
    : run_engine(compiled_tree, X, nsim, permutations, trait_chunk, store_sim,
                 include_observed, n_threads, simulation_chunk,
                 OracleCompute(), false);
  const bool supplied = permutations != R_NilValue;
  const int n_random = supplied ? nsim :
    (include_observed ? std::max(0, nsim - 1) : nsim);
  Rcpp::IntegerVector nsim_failed(result.failed);
  Rcpp::CharacterVector failure_reason(nsim, "");
  Rcpp::List memory = Rcpp::List::create(
    Rcpp::Named("bounded") = true,
    Rcpp::Named("working_memory_model") = use_alpha
      ? "production_compute_one_plus_call_level_alpha"
      : "production_compute_one",
    Rcpp::Named("alpha_precomputed") = use_alpha,
    Rcpp::Named("alpha_nodes") = static_cast<double>(result.alpha_nodes),
    Rcpp::Named("alpha_bytes") = static_cast<double>(result.alpha_bytes),
    Rcpp::Named("requested_sim_output_bytes") = store_sim
      ? static_cast<double>(static_cast<std::size_t>(nsim) *
                            static_cast<std::size_t>(X.ncol()) * sizeof(double))
      : 0.0,
    Rcpp::Named("note") = use_alpha
      ? "alpha is precomputed once after tree/cache construction and reused in both downward passes"
      : "same-translation-unit production compute_one oracle"
  );
  Rcpp::List accounting = Rcpp::List::create(
    Rcpp::Named("nsim_requested") = nsim,
    Rcpp::Named("n_randomizations") = n_random,
    Rcpp::Named("nsim_successful") = result.successful,
    Rcpp::Named("nsim_failed") = nsim_failed,
    Rcpp::Named("failure_reason") = failure_reason
  );
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("K") = result.K,
    Rcpp::Named("P") = result.P,
    Rcpp::Named("exceedance_count") = result.exceedance,
    Rcpp::Named("nsim_requested") = nsim,
    Rcpp::Named("nsim_successful") = result.successful,
    Rcpp::Named("nsim_failed") = nsim_failed,
    Rcpp::Named("status") = "ok",
    Rcpp::Named("failure_reason") = failure_reason,
    Rcpp::Named("n_randomizations") = n_random,
    Rcpp::Named("MCSE_P") = result.MCSE,
    Rcpp::Named("success") = result.success,
    Rcpp::Named("replicate_status") = result.success,
    Rcpp::Named("permutation_mode") = supplied ? "controlled" : "internal_rng",
    Rcpp::Named("include_observed") = !supplied && include_observed,
    Rcpp::Named("simulation_chunk") = simulation_chunk,
    Rcpp::Named("return_ordered") = return_ordered,
    Rcpp::Named("stage2d2_mode") = mode,
    Rcpp::Named("stage2d2_same_translation_unit") = true,
    Rcpp::Named("stage2d2_alpha_precomputed") = use_alpha,
    Rcpp::Named("stage2d2_alpha_nodes") = static_cast<double>(result.alpha_nodes),
    Rcpp::Named("stage2d2_alpha_bytes") = static_cast<double>(result.alpha_bytes),
    Rcpp::Named("memory") = memory,
    Rcpp::Named("accounting") = accounting,
    Rcpp::Named("candidate_metadata") = Rcpp::List::create(
      Rcpp::Named("candidate") = mode,
      Rcpp::Named("alpha_precomputed") = use_alpha,
      Rcpp::Named("alpha_reused_in_downward_passes") = use_alpha,
      Rcpp::Named("alpha_nodes") = static_cast<double>(result.alpha_nodes),
      Rcpp::Named("alpha_bytes") = static_cast<double>(result.alpha_bytes),
      Rcpp::Named("compute_one_equivalents") = static_cast<double>(nsim + 1),
      Rcpp::Named("n_randomizations") = n_random)
  );
  if (store_sim) out["sim_K"] = result.sim;
  return out;
}
