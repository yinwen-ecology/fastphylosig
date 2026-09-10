// Private Stage 2D2 Candidate A prototype for the K permutation/null engine.
//
// The production translation unit is included so the oracle and Candidate A
// use the same kperm::Tree, kperm::Cache, tip accessor, and Fisher-Yates
// implementation.  This file is audit-only: it is not registered in the
// package DLL and does not change R/, src/, tests/, or the public API.

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

#define STAGE2D2_PRODUCTION_CPP "../../../src/k_permutation.cpp"
#include STAGE2D2_PRODUCTION_CPP
#undef STAGE2D2_PRODUCTION_CPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace stage2d2 {

struct NumericWorkspace {
  std::vector<double> message;
  std::vector<double> state;
  std::vector<double> baseline;
  std::vector<double> delta;

  bool debug;
  std::uint64_t ensure_calls;
  std::uint64_t allocation_events;
  std::uint64_t reuse_calls;
  std::uint64_t debug_checks;
  std::uint64_t canary_checks;
  std::uint64_t stale_state_failures;
  std::uint64_t canary_failures;
  std::size_t peak_bytes;
  std::uint64_t canary_before;
  std::uint64_t canary_after;

  static double sentinel() {
    // The production inputs are finite.  Debug mode uses a NaN marker so a
    // stale message/state read becomes visible without changing normal mode.
    return std::numeric_limits<double>::quiet_NaN();
  }

  explicit NumericWorkspace(const bool debug_mode = false)
    : debug(debug_mode), ensure_calls(0), allocation_events(0),
      reuse_calls(0), debug_checks(0), canary_checks(0),
      stale_state_failures(0), canary_failures(0), peak_bytes(0),
      canary_before(0x5354324432414241ULL),
      canary_after(0x535432443241454EULL) {}

  std::size_t bytes() const {
    return (message.capacity() + state.capacity() +
            baseline.capacity() + delta.capacity()) * sizeof(double);
  }

  void ensure(const std::size_t cells, const std::size_t chunk) {
    ++ensure_calls;
    if (message.capacity() < cells) {
      message.resize(cells);
      ++allocation_events;
    }
    if (state.capacity() < cells) {
      state.resize(cells);
      ++allocation_events;
    }
    if (baseline.capacity() < chunk) {
      baseline.resize(chunk);
      ++allocation_events;
    }
    if (delta.capacity() < chunk) {
      delta.resize(chunk);
      ++allocation_events;
    }
    peak_bytes = std::max(peak_bytes, bytes());
  }

  void begin_chunk(const std::size_t cells, const std::size_t chunk) {
    ++reuse_calls;
    ensure(cells, chunk);
    if (!debug) return;
    std::fill(message.begin(), message.begin() +
              static_cast<std::ptrdiff_t>(cells), sentinel());
    std::fill(state.begin(), state.begin() +
              static_cast<std::ptrdiff_t>(cells), sentinel());
    std::fill(baseline.begin(), baseline.begin() +
              static_cast<std::ptrdiff_t>(chunk), sentinel());
    std::fill(delta.begin(), delta.begin() +
              static_cast<std::ptrdiff_t>(chunk), sentinel());
  }

  void check_canaries() {
    // These scalar canaries sit outside all reusable payloads.  They detect
    // accidental mutation of the workspace guard state in debug runs.
    ++canary_checks;
    if (debug && (ensure_calls == 0 || bytes() == 0 ||
                  canary_before != 0x5354324432414241ULL ||
                  canary_after != 0x535432443241454EULL)) {
      ++canary_failures;
      Rcpp::stop("Stage 2D2 workspace canary check failed.");
    }
  }

  void check_written(const std::vector<double>& values,
                     const std::size_t length, const char* label) {
    if (!debug) return;
    ++debug_checks;
    for (std::size_t i = 0; i < length; ++i) {
      if (std::isnan(values[i])) {
        ++stale_state_failures;
        Rcpp::stop(std::string("Stage 2D2 stale workspace value in ") + label);
      }
    }
  }

  void check_written_scalar(const std::vector<double>& values,
                            const std::size_t length, const char* label) {
    check_written(values, length, label);
  }
};

// This is a line-for-line numerical copy of kperm::compute_one(), with only
// the four local vectors replaced by a caller-owned reusable workspace.  The
// production oracle below calls kperm::compute_one() unchanged.
void compute_one_reusable(const kperm::Tree& tree, const kperm::Cache& cache,
                          const double* x, const int ncol,
                          const std::vector<int>& perm,
                          const int trait_chunk, NumericWorkspace& workspace,
                          std::vector<double>& out) {
  const int n = tree.n_tip;
  const int p = ncol;
  out.assign(static_cast<std::size_t>(p),
             std::numeric_limits<double>::quiet_NaN());
  const int chunk_limit = std::max(1, std::min(trait_chunk, p));

  for (int col0 = 0; col0 < p; col0 += chunk_limit) {
    const int chunk = std::min(chunk_limit, p - col0);
    const std::size_t cells = static_cast<std::size_t>(tree.n_total) *
      static_cast<std::size_t>(chunk);
    workspace.begin_chunk(cells, static_cast<std::size_t>(chunk));
    std::vector<double>& message = workspace.message;
    std::vector<double>& state = workspace.state;
    std::vector<double>& baseline = workspace.baseline;
    std::vector<double>& delta = workspace.delta;

    for (int j = 0; j < chunk; ++j) {
      baseline[static_cast<std::size_t>(j)] =
        kperm::tip_value(x, n, perm, 0, col0 + j);
    }
    workspace.check_written_scalar(baseline, static_cast<std::size_t>(chunk),
                                   "baseline");

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
    workspace.check_written(message, cells, "message after upward pass 1");

    // Downward conditional states and baseline-relative GLS offset delta.
    for (int j = 0; j < chunk; ++j) {
      state[static_cast<std::size_t>(tree.root) *
            static_cast<std::size_t>(chunk) + static_cast<std::size_t>(j)] = 0.0;
    }
    for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
      const int node = tree.preorder[ii];
      const int par = tree.parent[static_cast<std::size_t>(node)];
      const double alpha = node < n ? 1.0 :
        tree.branch[static_cast<std::size_t>(node)] *
        cache.outgoing[static_cast<std::size_t>(node)];
      const std::size_t base = static_cast<std::size_t>(node) *
        static_cast<std::size_t>(chunk);
      const std::size_t pbase = static_cast<std::size_t>(par) *
        static_cast<std::size_t>(chunk);
      for (int j = 0; j < chunk; ++j) {
        const double parent_state = pbase == base ? 0.0 :
          state[pbase + static_cast<std::size_t>(j)];
        state[base + static_cast<std::size_t>(j)] = parent_state +
          alpha * (message[base + static_cast<std::size_t>(j)] - parent_state);
      }
    }
    workspace.check_written(state, cells, "state after downward pass 1");

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
    workspace.check_written_scalar(delta, static_cast<std::size_t>(chunk),
                                   "delta");

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
    workspace.check_written(message, cells, "message after upward pass 2");

    for (int j = 0; j < chunk; ++j) {
      state[static_cast<std::size_t>(tree.root) *
            static_cast<std::size_t>(chunk) + static_cast<std::size_t>(j)] = 0.0;
    }
    for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
      const int node = tree.preorder[ii];
      const int par = tree.parent[static_cast<std::size_t>(node)];
      const double alpha = node < n ? 1.0 :
        tree.branch[static_cast<std::size_t>(node)] *
        cache.outgoing[static_cast<std::size_t>(node)];
      const std::size_t base = static_cast<std::size_t>(node) *
        static_cast<std::size_t>(chunk);
      const std::size_t pbase = static_cast<std::size_t>(par) *
        static_cast<std::size_t>(chunk);
      for (int j = 0; j < chunk; ++j) {
        const double parent_state = state[pbase + static_cast<std::size_t>(j)];
        state[base + static_cast<std::size_t>(j)] = parent_state +
          alpha * (message[base + static_cast<std::size_t>(j)] - parent_state);
      }
    }
    workspace.check_written(state, cells, "state after downward pass 2");

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
    workspace.check_canaries();
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
  std::size_t workspace_bytes;
  std::size_t workspace_peak_bytes;
  std::uint64_t workspace_allocations;
  std::uint64_t workspace_reuses;
  std::uint64_t debug_checks;
  std::uint64_t canary_checks;
  std::uint64_t stale_state_failures;
  std::uint64_t canary_failures;
};

void update_workspace_totals(const std::vector<NumericWorkspace>& workspaces,
                             std::size_t& bytes, std::size_t& peak_bytes,
                             std::uint64_t& allocations,
                             std::uint64_t& reuses,
                             std::uint64_t& debug_checks,
                             std::uint64_t& canary_checks,
                             std::uint64_t& stale_state_failures,
                             std::uint64_t& canary_failures) {
  bytes = 0;
  peak_bytes = 0;
  allocations = reuses = debug_checks = canary_checks = 0;
  stale_state_failures = canary_failures = 0;
  for (std::size_t i = 0; i < workspaces.size(); ++i) {
    bytes += workspaces[i].bytes();
    peak_bytes += workspaces[i].peak_bytes;
    allocations += workspaces[i].allocation_events;
    reuses += workspaces[i].reuse_calls;
    debug_checks += workspaces[i].debug_checks;
    canary_checks += workspaces[i].canary_checks;
    stale_state_failures += workspaces[i].stale_state_failures;
    canary_failures += workspaces[i].canary_failures;
  }
}

template <typename Compute>
EngineResult run_engine(const Rcpp::List& compiled_tree,
                        const Rcpp::NumericMatrix& X, const int nsim,
                        SEXP permutations, const int trait_chunk,
                        const bool store_sim, const bool include_observed,
                        const int n_threads, const int simulation_chunk,
                        Compute compute, const bool debug_workspace) {
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
  NumericWorkspace observed_workspace(debug_workspace);
  compute(tree, cache, x, p, identity, trait_chunk, observed_workspace,
          observed_vec);
  // The oracle callback uses a workspace-shaped signature but ignores it;
  // Candidate A uses the same call shape with the reusable implementation.
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
  std::vector<NumericWorkspace> workspaces(
    static_cast<std::size_t>(std::max(1, threads_eff)),
    NumericWorkspace(debug_workspace)
  );

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
        compute(tree, cache, x, p, perm, trait_chunk,
                workspaces[static_cast<std::size_t>(tid)], kval);
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
      compute(tree, cache, x, p, perm, trait_chunk, workspaces[0], kval);
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
            compute(tree, cache, x, p, perm, trait_chunk,
                    workspaces[static_cast<std::size_t>(tid)], kval);
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
          compute(tree, cache, x, p, perm, trait_chunk, workspaces[0], kval);
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

  std::size_t workspace_bytes = observed_workspace.bytes();
  std::size_t workspace_peak_bytes = observed_workspace.peak_bytes;
  std::uint64_t workspace_allocations = observed_workspace.allocation_events;
  std::uint64_t workspace_reuses = observed_workspace.reuse_calls;
  std::uint64_t debug_checks = observed_workspace.debug_checks;
  std::uint64_t canary_checks = observed_workspace.canary_checks;
  std::uint64_t stale_state_failures = observed_workspace.stale_state_failures;
  std::uint64_t canary_failures = observed_workspace.canary_failures;
  std::size_t worker_bytes = 0;
  std::size_t worker_peak_bytes = 0;
  std::uint64_t worker_allocations = 0;
  std::uint64_t worker_reuses = 0;
  std::uint64_t worker_debug_checks = 0;
  std::uint64_t worker_canary_checks = 0;
  std::uint64_t worker_stale_failures = 0;
  std::uint64_t worker_canary_failures = 0;
  update_workspace_totals(workspaces, worker_bytes, worker_peak_bytes,
                          worker_allocations, worker_reuses,
                          worker_debug_checks, worker_canary_checks,
                          worker_stale_failures, worker_canary_failures);
  workspace_bytes += worker_bytes;
  workspace_peak_bytes = std::max(workspace_peak_bytes,
                                  worker_peak_bytes + observed_workspace.peak_bytes);
  workspace_allocations += worker_allocations;
  workspace_reuses += worker_reuses;
  debug_checks += worker_debug_checks;
  canary_checks += worker_canary_checks;
  stale_state_failures += worker_stale_failures;
  canary_failures += worker_canary_failures;

  Rcpp::Nullable<Rcpp::NumericMatrix> sim_value(R_NilValue);
  if (store_sim) sim_value = Rcpp::Nullable<Rcpp::NumericMatrix>(sim);
  EngineResult result = {
    observed, p_value, exceedance, mcse, nsim_successful, nsim_failed,
    Rcpp::LogicalVector(nsim, true),
    sim_value,
    workspace_bytes, workspace_peak_bytes, workspace_allocations,
    workspace_reuses, debug_checks, canary_checks, stale_state_failures,
    canary_failures
  };
  return result;
}

struct OracleCompute {
  void operator()(const kperm::Tree& tree, const kperm::Cache& cache,
                  const double* x, const int ncol,
                  const std::vector<int>& perm, const int trait_chunk,
                  NumericWorkspace&, std::vector<double>& out) const {
    kperm::compute_one(tree, cache, x, ncol, perm, trait_chunk, out);
  }
};

struct CandidateACompute {
  void operator()(const kperm::Tree& tree, const kperm::Cache& cache,
                  const double* x, const int ncol,
                  const std::vector<int>& perm, const int trait_chunk,
                  NumericWorkspace& workspace, std::vector<double>& out) const {
    compute_one_reusable(tree, cache, x, ncol, perm, trait_chunk, workspace, out);
  }
};

} // namespace stage2d2

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
    const bool return_ordered = true,
    const bool debug_workspace = false) {
  using namespace stage2d2;
  if (mode != "oracle" && mode != "A") {
    Rcpp::stop("Stage 2D2 prototype supports only mode='oracle' or mode='A'.");
  }
  const bool store_sim = return_sim || return_ordered;
  EngineResult result = mode == "oracle"
    ? run_engine(compiled_tree, X, nsim, permutations, trait_chunk, store_sim,
                 include_observed, n_threads, simulation_chunk, OracleCompute(),
                 debug_workspace)
    : run_engine(compiled_tree, X, nsim, permutations, trait_chunk, store_sim,
                 include_observed, n_threads, simulation_chunk, CandidateACompute(),
                 debug_workspace);

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("K") = result.K,
    Rcpp::Named("P") = result.P,
    Rcpp::Named("exceedance_count") = result.exceedance,
    Rcpp::Named("nsim_requested") = nsim,
    Rcpp::Named("nsim_successful") = result.successful,
    Rcpp::Named("nsim_failed") = result.failed,
    Rcpp::Named("status") = "ok",
    Rcpp::Named("failure_reason") = Rcpp::CharacterVector(nsim, ""),
    Rcpp::Named("n_randomizations") =
      permutations == R_NilValue
        ? (include_observed ? std::max(0, nsim - 1) : nsim)
        : nsim,
    Rcpp::Named("MCSE_P") = result.MCSE,
    Rcpp::Named("success") = result.success,
    Rcpp::Named("permutation_mode") =
      permutations == R_NilValue ? "internal_rng" : "controlled",
    Rcpp::Named("include_observed") =
      permutations == R_NilValue && include_observed,
    Rcpp::Named("simulation_chunk") = simulation_chunk,
    Rcpp::Named("stage2d2_mode") = mode,
    Rcpp::Named("stage2d2_same_translation_unit") = true,
    Rcpp::Named("stage2d2_workspace_bytes") =
      static_cast<double>(result.workspace_bytes),
    Rcpp::Named("stage2d2_workspace_peak_bytes") =
      static_cast<double>(result.workspace_peak_bytes),
    Rcpp::Named("stage2d2_workspace_allocations") =
      static_cast<double>(result.workspace_allocations),
    Rcpp::Named("stage2d2_workspace_reuses") =
      static_cast<double>(result.workspace_reuses),
    Rcpp::Named("stage2d2_debug_checks") =
      static_cast<double>(result.debug_checks),
    Rcpp::Named("stage2d2_canary_checks") =
      static_cast<double>(result.canary_checks),
    Rcpp::Named("stage2d2_stale_state_failures") =
      static_cast<double>(result.stale_state_failures),
    Rcpp::Named("stage2d2_canary_failures") =
      static_cast<double>(result.canary_failures)
  );
  if (store_sim) out["sim_K"] = result.sim;
  return out;
}
