// Stage 2D2B Candidate B remediation: true blocked K null evaluator.
//
// This file is audit-only.  It includes the frozen production translation
// unit so that the exported oracle can be compared in the same translation
// unit, but the Candidate B path below never calls compute_one(),
// fast_k_tree_permutation_cpp(), or the oracle.  The candidate implements
// the scalar compute_one arithmetic directly over a bounded block of
// permutations.

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

// Hide the literal .cpp include from Rcpp::sourceCpp.  The production source
// is compiled once as part of this translation unit and supplies the frozen
// Tree/Cache parser, oracle, and Fisher-Yates reference implementation.
#define STAGE2D2B_PRODUCTION_CPP "../../../src/k_permutation.cpp"
#include STAGE2D2B_PRODUCTION_CPP
#undef STAGE2D2B_PRODUCTION_CPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace stage2d2b {

struct Counters {
  std::uint64_t blocks;
  std::uint64_t generated_permutations;
  std::uint64_t fisher_yates_swaps;
  std::uint64_t candidate_compute_count;
  std::uint64_t candidate_old_compute_one_calls;
  std::uint64_t oracle_call_count;
  std::uint64_t workspace_allocations;
  std::uint64_t workspace_reuses;
  std::uint64_t debug_checks;
  std::uint64_t canary_checks;
  std::uint64_t stale_state_failures;
  std::uint64_t canary_failures;

  Counters()
    : blocks(0), generated_permutations(0), fisher_yates_swaps(0),
      candidate_compute_count(0), candidate_old_compute_one_calls(0),
      oracle_call_count(0), workspace_allocations(0), workspace_reuses(0),
      debug_checks(0), canary_checks(0), stale_state_failures(0),
      canary_failures(0) {}
};

inline bool valid_block_size(const int block_size) {
  return block_size == 2 || block_size == 4 || block_size == 8 ||
         block_size == 16;
}

// Node-major/replicate-minor storage keeps the structural traversal shared
// by all members of a block.  The last dimension is the trait chunk, so a
// replicate's arithmetic remains independent of all other replicates.
inline std::size_t workspace_index(const int node, const int replicate,
                                   const int total, const int block_count,
                                   const int trait_in_chunk,
                                   const int chunk) {
  return (static_cast<std::size_t>(node) *
          static_cast<std::size_t>(block_count) +
          static_cast<std::size_t>(replicate)) *
         static_cast<std::size_t>(chunk) +
         static_cast<std::size_t>(trait_in_chunk);
}

inline std::size_t block_scalar_index(const int replicate,
                                      const int trait_in_chunk,
                                      const int chunk) {
  return static_cast<std::size_t>(replicate) *
           static_cast<std::size_t>(chunk) +
         static_cast<std::size_t>(trait_in_chunk);
}

inline double block_tip_value(const double* x, const int n_tip,
                              const std::vector<int>& permutations,
                              const int replicate, const int tip,
                              const int trait) {
  const int source_tip = permutations[
    static_cast<std::size_t>(replicate) * static_cast<std::size_t>(n_tip) +
    static_cast<std::size_t>(tip)
  ];
  return x[static_cast<std::size_t>(source_tip) +
           static_cast<std::size_t>(n_tip) * static_cast<std::size_t>(trait)];
}

// One bounded workspace is reused for every block in the serial evaluator.
// Debug mode poisons only the active cells before each chunk and checks that
// every cell read by a later pass was written in the current pass.
struct BlockWorkspace {
  std::vector<double> message;
  std::vector<double> state;
  std::vector<double> baseline;
  std::vector<double> delta;
  bool debug;
  std::uint64_t allocation_events;
  std::uint64_t reuse_calls;
  std::uint64_t debug_checks;
  std::uint64_t canary_checks;
  std::uint64_t stale_state_failures;
  std::uint64_t canary_failures;
  std::size_t peak_bytes;
  std::uint64_t canary_before;
  std::uint64_t canary_after;

  explicit BlockWorkspace(const bool debug_mode = false)
    : debug(debug_mode), allocation_events(0), reuse_calls(0),
      debug_checks(0), canary_checks(0), stale_state_failures(0),
      canary_failures(0), peak_bytes(0),
      canary_before(0x5354324432424241ULL),
      canary_after(0x535432443242454EULL) {}

  std::size_t bytes() const {
    return (message.capacity() + state.capacity() +
            baseline.capacity() + delta.capacity()) * sizeof(double);
  }

  void ensure(const std::size_t cells, const std::size_t scalar_cells) {
    if (message.capacity() < cells) {
      message.resize(cells);
      ++allocation_events;
    } else {
      message.resize(cells);
    }
    if (state.capacity() < cells) {
      state.resize(cells);
      ++allocation_events;
    } else {
      state.resize(cells);
    }
    if (baseline.capacity() < scalar_cells) {
      baseline.resize(scalar_cells);
      ++allocation_events;
    } else {
      baseline.resize(scalar_cells);
    }
    if (delta.capacity() < scalar_cells) {
      delta.resize(scalar_cells);
      ++allocation_events;
    } else {
      delta.resize(scalar_cells);
    }
    peak_bytes = std::max(peak_bytes, bytes());
  }

  void begin_chunk(const std::size_t cells, const std::size_t scalar_cells) {
    ++reuse_calls;
    ensure(cells, scalar_cells);
    if (!debug) return;
    const double marker = std::numeric_limits<double>::quiet_NaN();
    std::fill(message.begin(), message.end(), marker);
    std::fill(state.begin(), state.end(), marker);
    std::fill(baseline.begin(), baseline.end(), marker);
    std::fill(delta.begin(), delta.end(), marker);
  }

  void check_written(const std::vector<double>& values,
                     const char* label) {
    if (!debug) return;
    ++debug_checks;
    for (std::size_t i = 0; i < values.size(); ++i) {
      if (std::isnan(values[i])) {
        ++stale_state_failures;
        Rcpp::stop(std::string("Stage 2D2B stale workspace value in ") +
                   label);
      }
    }
  }

  void check_canaries() {
    ++canary_checks;
    if (debug && (canary_before != 0x5354324432424241ULL ||
                  canary_after != 0x535432443242454EULL)) {
      ++canary_failures;
      Rcpp::stop("Stage 2D2B workspace canary check failed.");
    }
  }
};

// This is the frozen compute_one() arithmetic lifted over a block.  For
// every fixed (replicate, trait), node/child traversal, operation order, and
// long-double reduction order are identical to the scalar implementation.
// Only independent replicates are interleaved at the outer loop level.
void compute_block(const kperm::Tree& tree, const kperm::Cache& cache,
                   const double* x, const int ncol,
                   const std::vector<int>& permutations,
                   const int block_count, const int trait_chunk,
                   BlockWorkspace& workspace, std::vector<double>& block_k,
                   Counters& counters) {
  const int n = tree.n_tip;
  const int total = tree.n_total;
  const int p = ncol;
  const int chunk_limit = std::max(1, std::min(trait_chunk, p));
  const double nan = std::numeric_limits<double>::quiet_NaN();

  if (block_count < 1 || block_count > 16) {
    Rcpp::stop("block_count must be between 1 and 16.");
  }
  if (static_cast<std::size_t>(block_count) * static_cast<std::size_t>(n) !=
      permutations.size()) {
    Rcpp::stop("permutation block has incompatible dimensions.");
  }
  block_k.assign(static_cast<std::size_t>(block_count) *
                 static_cast<std::size_t>(p), nan);
  counters.candidate_compute_count +=
    static_cast<std::uint64_t>(block_count);

  for (int col0 = 0; col0 < p; col0 += chunk_limit) {
    const int chunk = std::min(chunk_limit, p - col0);
    const std::size_t cells = static_cast<std::size_t>(total) *
      static_cast<std::size_t>(block_count) * static_cast<std::size_t>(chunk);
    const std::size_t scalar_cells = static_cast<std::size_t>(block_count) *
      static_cast<std::size_t>(chunk);
    workspace.begin_chunk(cells, scalar_cells);
    std::vector<double>& message = workspace.message;
    std::vector<double>& state = workspace.state;
    std::vector<double>& baseline = workspace.baseline;
    std::vector<double>& delta = workspace.delta;

    for (int r = 0; r < block_count; ++r) {
      for (int j = 0; j < chunk; ++j) {
        baseline[block_scalar_index(r, j, chunk)] =
          block_tip_value(x, n, permutations, r, 0, col0 + j);
      }
    }
    workspace.check_written(baseline, "baseline");

    // Upward Gaussian messages for x - baseline.
    for (std::size_t ii = 0; ii < tree.postorder.size(); ++ii) {
      const int node = tree.postorder[ii];
      const bool is_tip = node < n;
      const std::size_t begin = static_cast<std::size_t>(
        tree.child_ptr[static_cast<std::size_t>(node)]);
      const std::size_t end = static_cast<std::size_t>(
        tree.child_ptr[static_cast<std::size_t>(node + 1)]);
      const double s = is_tip ? 0.0 :
        cache.aggregate[static_cast<std::size_t>(node)];
      for (int r = 0; r < block_count; ++r) {
        for (int j = 0; j < chunk; ++j) {
          const std::size_t here = workspace_index(
            node, r, total, block_count, j, chunk);
          if (is_tip) {
            message[here] = block_tip_value(
              x, n, permutations, r, node, col0 + j) -
              baseline[block_scalar_index(r, j, chunk)];
          } else {
            long double weighted = 0.0L;
            for (std::size_t k = begin; k < end; ++k) {
              const int child = tree.children[k];
              weighted += static_cast<long double>(cache.outgoing[
                static_cast<std::size_t>(child)]) *
                static_cast<long double>(message[workspace_index(
                  child, r, total, block_count, j, chunk)]);
            }
            message[here] = static_cast<double>(
              weighted / static_cast<long double>(s));
          }
        }
      }
    }
    workspace.check_written(message, "message after upward pass 1");

    // Downward conditional states and baseline-relative GLS offset delta.
    for (int r = 0; r < block_count; ++r) {
      for (int j = 0; j < chunk; ++j) {
        state[workspace_index(tree.root, r, total, block_count, j, chunk)] =
          0.0;
      }
    }
    for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
      const int node = tree.preorder[ii];
      const int par = tree.parent[static_cast<std::size_t>(node)];
      const double alpha = node < n ? 1.0 :
        tree.branch[static_cast<std::size_t>(node)] *
        cache.outgoing[static_cast<std::size_t>(node)];
      for (int r = 0; r < block_count; ++r) {
        for (int j = 0; j < chunk; ++j) {
          const std::size_t here = workspace_index(
            node, r, total, block_count, j, chunk);
          const std::size_t parent_index = workspace_index(
            par, r, total, block_count, j, chunk);
          const double parent_state = par == node ? 0.0 : state[parent_index];
          state[here] = parent_state +
            alpha * (message[here] - parent_state);
        }
      }
    }
    workspace.check_written(state, "state after downward pass 1");

    for (int r = 0; r < block_count; ++r) {
      for (int j = 0; j < chunk; ++j) {
        long double qlinear = 0.0L;
        for (int tip = 0; tip < n; ++tip) {
          const int par = tree.parent[static_cast<std::size_t>(tip)];
          const double parent_state = state[workspace_index(
            par, r, total, block_count, j, chunk)];
          qlinear += static_cast<long double>(
            block_tip_value(x, n, permutations, r, tip, col0 + j) -
            baseline[block_scalar_index(r, j, chunk)] - parent_state) /
            static_cast<long double>(tree.branch[
              static_cast<std::size_t>(tip)]);
        }
        delta[block_scalar_index(r, j, chunk)] = static_cast<double>(
          qlinear / static_cast<long double>(cache.sum_inv));
      }
    }
    workspace.check_written(delta, "delta");

    // Repeat on (x - baseline) - delta for numerator and precision energy.
    for (std::size_t ii = 0; ii < tree.postorder.size(); ++ii) {
      const int node = tree.postorder[ii];
      const bool is_tip = node < n;
      const std::size_t begin = static_cast<std::size_t>(
        tree.child_ptr[static_cast<std::size_t>(node)]);
      const std::size_t end = static_cast<std::size_t>(
        tree.child_ptr[static_cast<std::size_t>(node + 1)]);
      const double s = is_tip ? 0.0 :
        cache.aggregate[static_cast<std::size_t>(node)];
      for (int r = 0; r < block_count; ++r) {
        for (int j = 0; j < chunk; ++j) {
          const std::size_t here = workspace_index(
            node, r, total, block_count, j, chunk);
          if (is_tip) {
            message[here] = block_tip_value(
              x, n, permutations, r, node, col0 + j) -
              baseline[block_scalar_index(r, j, chunk)] -
              delta[block_scalar_index(r, j, chunk)];
          } else {
            long double weighted = 0.0L;
            for (std::size_t k = begin; k < end; ++k) {
              const int child = tree.children[k];
              weighted += static_cast<long double>(cache.outgoing[
                static_cast<std::size_t>(child)]) *
                static_cast<long double>(message[workspace_index(
                  child, r, total, block_count, j, chunk)]);
            }
            message[here] = static_cast<double>(
              weighted / static_cast<long double>(s));
          }
        }
      }
    }
    workspace.check_written(message, "message after upward pass 2");

    for (int r = 0; r < block_count; ++r) {
      for (int j = 0; j < chunk; ++j) {
        state[workspace_index(tree.root, r, total, block_count, j, chunk)] =
          0.0;
      }
    }
    for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
      const int node = tree.preorder[ii];
      const int par = tree.parent[static_cast<std::size_t>(node)];
      const double alpha = node < n ? 1.0 :
        tree.branch[static_cast<std::size_t>(node)] *
        cache.outgoing[static_cast<std::size_t>(node)];
      for (int r = 0; r < block_count; ++r) {
        for (int j = 0; j < chunk; ++j) {
          const std::size_t here = workspace_index(
            node, r, total, block_count, j, chunk);
          const std::size_t parent_index = workspace_index(
            par, r, total, block_count, j, chunk);
          const double parent_state = state[parent_index];
          state[here] = parent_state +
            alpha * (message[here] - parent_state);
        }
      }
    }
    workspace.check_written(state, "state after downward pass 2");

    for (int r = 0; r < block_count; ++r) {
      for (int j = 0; j < chunk; ++j) {
        long double numerator = 0.0L;
        long double denominator = 0.0L;
        for (int tip = 0; tip < n; ++tip) {
          const long double y = static_cast<long double>(
            block_tip_value(x, n, permutations, r, tip, col0 + j) -
            baseline[block_scalar_index(r, j, chunk)] -
            delta[block_scalar_index(r, j, chunk)]);
          numerator += y * y;
        }
        for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
          const int node = tree.preorder[ii];
          const int par = tree.parent[static_cast<std::size_t>(node)];
          const long double parent_state = static_cast<long double>(state[
            workspace_index(par, r, total, block_count, j, chunk)]);
          const long double child_state = node < n
            ? static_cast<long double>(
                block_tip_value(x, n, permutations, r, node, col0 + j) -
                baseline[block_scalar_index(r, j, chunk)] -
                delta[block_scalar_index(r, j, chunk)])
            : static_cast<long double>(state[workspace_index(
                node, r, total, block_count, j, chunk)]);
          const long double diff = child_state - parent_state;
          denominator += diff * diff /
            static_cast<long double>(tree.branch[
              static_cast<std::size_t>(node)]);
        }
        const double den = static_cast<double>(denominator);
        const double num = static_cast<double>(numerator);
        block_k[static_cast<std::size_t>(r) * static_cast<std::size_t>(p) +
                static_cast<std::size_t>(col0 + j)] =
          den > 0.0 && std::isfinite(den) && std::isfinite(num)
          ? (num / den) / cache.normalization : nan;
      }
    }
    workspace.check_canaries();
  }
}

inline void fisher_yates_serial(std::vector<int>& permutation,
                                Counters& counters) {
  for (int i = static_cast<int>(permutation.size()) - 1; i > 0; --i) {
    const int j = static_cast<int>(std::floor(
      R::runif(0.0, static_cast<double>(i + 1))));
    std::swap(permutation[static_cast<std::size_t>(i)],
              permutation[static_cast<std::size_t>(j)]);
    ++counters.fisher_yates_swaps;
  }
}

inline void validate_supplied(const Rcpp::IntegerMatrix& permutations,
                              const int nsim, const int n) {
  if (permutations.nrow() != nsim || permutations.ncol() != n) {
    Rcpp::stop("permutations must have nsim rows and one column per tip.");
  }
  std::vector<int> seen(static_cast<std::size_t>(n), -1);
  for (int i = 0; i < nsim; ++i) {
    for (int r = 0; r < n; ++r) {
      const int idx = permutations(i, r);
      if (idx < 1 || idx > n ||
          seen[static_cast<std::size_t>(idx - 1)] == i) {
        Rcpp::stop("each permutation row must contain 1:n exactly once.");
      }
      seen[static_cast<std::size_t>(idx - 1)] = i;
    }
  }
}

Rcpp::List add_audit_fields(Rcpp::List out, const int nsim, const int p,
                             const int n_randomizations,
                             const bool include_observed,
                             const int block_size,
                             const std::size_t peak_index_bytes,
                             const std::size_t peak_workspace_bytes,
                             const std::size_t requested_sim_output_bytes,
                             const BlockWorkspace& workspace,
                             const Counters& counters,
                             const Rcpp::NumericVector& observed,
                             const std::string& candidate_name,
                             const bool blocked) {
  Rcpp::IntegerVector successful(p, nsim), failed(p, 0);
  for (int j = 0; j < p; ++j) {
    if (j >= observed.size() || !std::isfinite(observed[j])) {
      successful[j] = 0;
      failed[j] = nsim;
    }
  }
  out["nsim_successful"] = successful;
  Rcpp::LogicalVector replicate_status(nsim, true);
  Rcpp::CharacterVector failure_reason(nsim, "");
  out["nsim_failed"] = failed;
  out["status"] = "ok";
  out["replicate_status"] = replicate_status;
  out["failure_reason"] = failure_reason;
  out["memory"] = Rcpp::List::create(
    Rcpp::Named("bounded") = blocked,
    Rcpp::Named("working_memory_model") = "n_times_block",
    Rcpp::Named("block_size") = block_size,
    Rcpp::Named("peak_index_bytes") =
      static_cast<double>(peak_index_bytes),
    Rcpp::Named("peak_workspace_bytes") =
      static_cast<double>(peak_workspace_bytes),
    Rcpp::Named("requested_sim_output_bytes") =
      static_cast<double>(requested_sim_output_bytes),
    Rcpp::Named("note") =
      "candidate workspace is bounded by n x block; sim_K output is separate"
  );
  out["candidate_metadata"] = Rcpp::List::create(
    Rcpp::Named("candidate") = candidate_name,
    Rcpp::Named("blocked") = blocked,
    Rcpp::Named("serial_permutation_generation") = true,
    Rcpp::Named("blocks_processed") =
      static_cast<double>(counters.blocks),
    Rcpp::Named("generated_permutations") =
      static_cast<double>(counters.generated_permutations),
    Rcpp::Named("fisher_yates_swaps") =
      static_cast<double>(counters.fisher_yates_swaps),
    Rcpp::Named("compute_one_equivalents") =
      static_cast<double>(counters.candidate_compute_count),
    Rcpp::Named("candidate_compute_count") =
      static_cast<double>(counters.candidate_compute_count),
    Rcpp::Named("candidate_old_compute_one_calls") =
      static_cast<double>(counters.candidate_old_compute_one_calls),
    Rcpp::Named("oracle_call_count") =
      static_cast<double>(counters.oracle_call_count),
    Rcpp::Named("n_randomizations") = n_randomizations,
    Rcpp::Named("include_observed") = include_observed,
    Rcpp::Named("workspace_allocations") =
      static_cast<double>(workspace.allocation_events),
    Rcpp::Named("workspace_reuses") =
      static_cast<double>(workspace.reuse_calls),
    Rcpp::Named("debug_checks") =
      static_cast<double>(workspace.debug_checks),
    Rcpp::Named("canary_checks") =
      static_cast<double>(workspace.canary_checks),
    Rcpp::Named("stale_state_failures") =
      static_cast<double>(workspace.stale_state_failures),
    Rcpp::Named("canary_failures") =
      static_cast<double>(workspace.canary_failures)
  );
  out["candidate_old_compute_one_calls"] = static_cast<double>(
    counters.candidate_old_compute_one_calls);
  out["oracle_call_count"] = static_cast<double>(counters.oracle_call_count);
  out["candidate_compute_count"] = static_cast<double>(
    counters.candidate_compute_count);
  return out;
}

}  // namespace stage2d2b

// Same-TU frozen oracle.  Only this function calls the production entry
// point; Candidate B below has no production/oracle fallback.
// [[Rcpp::export]]
Rcpp::List stage2d2_k_permutation_oracle(
    const Rcpp::List& compiled_tree, const Rcpp::NumericMatrix& X,
    const int nsim = 1000, SEXP permutations = R_NilValue,
    const int trait_chunk = 64, const bool return_sim = true,
    const bool include_observed = true, const int n_threads = 1,
    const int simulation_chunk = 128, const int block_size = 4,
    const bool return_ordered = true, const bool collect_counters = false) {
  (void)block_size;
  (void)return_ordered;
  (void)collect_counters;
  stage2d2b::Counters counters;
  const bool supplied = permutations != R_NilValue;
  const int p = X.ncol();
  const int n_randomizations = supplied ? nsim :
    (include_observed ? std::max(0, nsim - 1) : nsim);
  Rcpp::List out = fast_k_tree_permutation_cpp(
    compiled_tree, X, nsim, permutations, trait_chunk, return_sim,
    include_observed, n_threads, simulation_chunk);
  // The frozen production call evaluates observed plus every null row.
  counters.oracle_call_count = static_cast<std::uint64_t>(nsim) + 1ULL;
  Rcpp::NumericVector oracle_observed = out["K"];
  return stage2d2b::add_audit_fields(
    out, nsim, p, n_randomizations, !supplied && include_observed,
    0, 0, 0, 0, stage2d2b::BlockWorkspace(false), counters,
    oracle_observed, "oracle", false
  );
}

// True Candidate B remediation.  All observed and null K values are
// generated by stage2d2b::compute_block().  The candidate path intentionally
// does not call the production evaluator or the oracle.
Rcpp::List run_candidate_b(
    const Rcpp::List& compiled_tree, const Rcpp::NumericMatrix& X,
    const int nsim = 1000, SEXP permutations = R_NilValue,
    const int trait_chunk = 64, const bool return_sim = true,
    const bool include_observed = true, const int n_threads = 1,
    const int simulation_chunk = 128, const int block_size = 4,
    const bool return_ordered = true, const bool collect_counters = false,
    const bool debug_workspace = false) {
  using namespace stage2d2b;
  (void)n_threads;
  (void)collect_counters;
  if (nsim < 1) Rcpp::stop("nsim must be a positive integer.");
  if (trait_chunk < 1) Rcpp::stop("trait_chunk must be positive.");
  if (simulation_chunk < 1) {
    Rcpp::stop("simulation_chunk must be a positive integer.");
  }
  if (!valid_block_size(block_size)) {
    Rcpp::stop("block_size must be one of 2, 4, 8, or 16.");
  }

  // Match the frozen production validation order before any candidate work.
  const kperm::Tree tree = kperm::parse_tree(compiled_tree);
  const int n = tree.n_tip;
  const int p = X.ncol();
  if (X.nrow() != n || p < 1) {
    Rcpp::stop("X must have one row per compiled tree tip and at least one trait.");
  }
  const double* x = REAL(X);
  for (R_xlen_t i = 0; i < X.size(); ++i) {
    if (!std::isfinite(x[i])) {
      Rcpp::stop("X must contain only finite values.");
    }
  }
  const kperm::Cache cache = kperm::build_cache(tree);

  const bool supplied = permutations != R_NilValue;
  Rcpp::IntegerMatrix supplied_permutations;
  if (supplied) {
    supplied_permutations = Rcpp::IntegerMatrix(permutations);
    validate_supplied(supplied_permutations, nsim, n);
  }

  Counters counters;
  Rcpp::RNGScope rng_scope;
  std::vector<int> identity(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) identity[static_cast<std::size_t>(i)] = i;

  // Observed K is independently evaluated as a one-member block.  It is not
  // copied from the oracle or from the production evaluator.
  BlockWorkspace workspace(debug_workspace);
  std::vector<double> observed_block;
  compute_block(tree, cache, x, p, identity, 1, trait_chunk, workspace,
                observed_block, counters);
  Rcpp::NumericVector observed(p);
  for (int j = 0; j < p; ++j) {
    observed[j] = observed_block[static_cast<std::size_t>(j)];
  }
  std::vector<char> valid_observed(static_cast<std::size_t>(p), 1);
  for (int j = 0; j < p; ++j) {
    valid_observed[static_cast<std::size_t>(j)] =
      std::isfinite(observed[j]) ? 1 : 0;
  }

  const bool result_include_observed = !supplied && include_observed;
  const int n_randomizations = supplied ? nsim :
    (result_include_observed ? std::max(0, nsim - 1) : nsim);
  const bool store_sim = return_sim || return_ordered;
  Rcpp::NumericMatrix sim;
  if (store_sim) sim = Rcpp::NumericMatrix(nsim, p);
  Rcpp::NumericVector exceedance(p, 0.0);

  std::vector<int> block_permutations(
    static_cast<std::size_t>(block_size) * static_cast<std::size_t>(n));
  std::vector<int> one_permutation(static_cast<std::size_t>(n));
  std::vector<double> block_k;
  std::size_t peak_index_bytes = 0;
  int processed = 0;
  while (processed < nsim) {
    const int count = std::min(block_size, nsim - processed);
    ++counters.blocks;
    const std::size_t block_cells = static_cast<std::size_t>(count) *
      static_cast<std::size_t>(n);
    if (block_permutations.size() != block_cells) {
      block_permutations.resize(block_cells);
    }
    peak_index_bytes = std::max(peak_index_bytes,
      block_permutations.size() * sizeof(int));

    for (int r = 0; r < count; ++r) {
      const int global = processed + r;
      if (supplied) {
        for (int tip = 0; tip < n; ++tip) {
          block_permutations[static_cast<std::size_t>(r) *
                             static_cast<std::size_t>(n) +
                             static_cast<std::size_t>(tip)] =
            supplied_permutations(global, tip) - 1;
        }
      } else {
        one_permutation = identity;
        if (!(result_include_observed && global == 0)) {
          fisher_yates_serial(one_permutation, counters);
          ++counters.generated_permutations;
        }
        for (int tip = 0; tip < n; ++tip) {
          block_permutations[static_cast<std::size_t>(r) *
                             static_cast<std::size_t>(n) +
                             static_cast<std::size_t>(tip)] =
            one_permutation[static_cast<std::size_t>(tip)];
        }
      }
    }

    // The result is consumed directly from the candidate block output.  No
    // scalar production call is made after this computation.
    compute_block(tree, cache, x, p, block_permutations, count,
                  trait_chunk, workspace, block_k, counters);
    for (int r = 0; r < count; ++r) {
      const int global = processed + r;
      for (int j = 0; j < p; ++j) {
        const double value = block_k[
          static_cast<std::size_t>(r) * static_cast<std::size_t>(p) +
          static_cast<std::size_t>(j)];
        if (store_sim) sim(global, j) = value;
        if (valid_observed[static_cast<std::size_t>(j)] &&
            fastphylosig::inclusive_upper_tail(
              value, observed[j])) {
          exceedance[j] += 1.0;
        }
      }
    }
    processed += count;
  }

  Rcpp::NumericVector p_value(p), mcse(p);
  Rcpp::IntegerVector successful(p, nsim), failed(p, 0);
  for (int j = 0; j < p; ++j) {
    if (!valid_observed[static_cast<std::size_t>(j)]) {
      p_value[j] = NA_REAL;
      mcse[j] = NA_REAL;
      successful[j] = 0;
      failed[j] = nsim;
      continue;
    }
    p_value[j] = exceedance[j] / static_cast<double>(nsim);
    if (!supplied && result_include_observed && n_randomizations > 0) {
      const double q = std::max(0.0, std::min(1.0,
        (exceedance[j] - 1.0) / static_cast<double>(n_randomizations)));
      mcse[j] = std::sqrt(q * (1.0 - q) *
                           static_cast<double>(n_randomizations)) /
        static_cast<double>(nsim);
    } else if (n_randomizations > 0) {
      const double q = p_value[j];
      mcse[j] = std::sqrt(q * (1.0 - q) /
                           static_cast<double>(n_randomizations));
    } else {
      mcse[j] = NA_REAL;
    }
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("K") = observed,
    Rcpp::Named("P") = p_value,
    Rcpp::Named("exceedance_count") = exceedance,
    Rcpp::Named("nsim_requested") = nsim,
    Rcpp::Named("nsim_successful") = successful,
    Rcpp::Named("n_randomizations") = n_randomizations,
    Rcpp::Named("MCSE_P") = mcse,
    Rcpp::Named("permutation_mode") = supplied ? "controlled" : "internal_rng",
    Rcpp::Named("include_observed") = result_include_observed,
    Rcpp::Named("simulation_chunk") = simulation_chunk,
    Rcpp::Named("sum_inv") = cache.sum_inv,
    Rcpp::Named("normalization") = cache.normalization,
    Rcpp::Named("stage2d2_mode") = "B-remediation",
    Rcpp::Named("stage2d2_same_translation_unit") = true,
    Rcpp::Named("n_threads_requested") = n_threads
  );
  if (store_sim) out["sim_K"] = sim;

  const std::size_t requested_sim_output_bytes = store_sim
    ? static_cast<std::size_t>(nsim) * static_cast<std::size_t>(p) *
      sizeof(double) : 0;
  const std::size_t peak_workspace_bytes = workspace.peak_bytes;
  out = add_audit_fields(
    out, nsim, p, n_randomizations, result_include_observed, block_size,
    peak_index_bytes, peak_workspace_bytes, requested_sim_output_bytes,
    workspace, counters, observed, "B-remediation", true
  );
  return out;
}

// Compatibility export for runners that use the historical Candidate B
// symbol.  It is an alias to the same independent blocked implementation.
// [[Rcpp::export]]
Rcpp::List stage2d2_k_permutation_candidate_b(
    const Rcpp::List& compiled_tree, const Rcpp::NumericMatrix& X,
    const int nsim = 1000, SEXP permutations = R_NilValue,
    const int trait_chunk = 64, const bool return_sim = true,
    const bool include_observed = true, const int n_threads = 1,
    const int simulation_chunk = 128, const int block_size = 4,
    const bool return_ordered = true, const bool collect_counters = false,
    const bool debug_workspace = false) {
  return run_candidate_b(
    compiled_tree, X, nsim, permutations, trait_chunk, return_sim,
    include_observed, n_threads, simulation_chunk, block_size,
    return_ordered, collect_counters, debug_workspace
  );
}
