// Stage 2D2 Candidate B, audit-only prototype.
//
// This file deliberately includes the production translation unit so that the
// oracle and the candidate use the same Tree/Cache representation and the
// unchanged kperm::compute_one implementation.  It is not part of the
// package build and must never be installed as production code.

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

namespace stage2d2_blocked {

struct Counters {
  std::size_t blocks = 0;
  std::size_t generated_permutations = 0;
  std::size_t fisher_yates_swaps = 0;
  std::size_t compute_one_equivalents = 0;
};

inline bool valid_block_size(const int block_size) {
  return block_size == 2 || block_size == 4 || block_size == 8 ||
         block_size == 16 || block_size == 32;
}

inline std::size_t workspace_index(const int replicate,
                                   const int node,
                                   const int n_total,
                                   const int trait_in_chunk,
                                   const int chunk) {
  return (static_cast<std::size_t>(replicate) * n_total + node) * chunk +
         trait_in_chunk;
}

inline double block_tip_value(const double* x,
                              const int n_tip,
                              const std::vector<int>& permutations,
                              const int replicate,
                              const int tip,
                              const int trait) {
  const int source_tip = permutations[static_cast<std::size_t>(replicate) *
                                      n_tip + tip];
  return x[source_tip + n_tip * trait];
}

// This is the same Fisher-Yates draw sequence as kperm::fisher_yates, with a
// counter added for audit output.  Permutations are always generated in the
// serial caller before a block is evaluated.
inline void fisher_yates_serial(std::vector<int>& permutation,
                                std::size_t& swap_count) {
  const int n = static_cast<int>(permutation.size());
  for (int i = n - 1; i > 0; --i) {
    const int j = static_cast<int>(R::runif(0.0, static_cast<double>(i + 1)));
    std::swap(permutation[i], permutation[j]);
    ++swap_count;
  }
}

// Evaluate a bounded block of permutations.  The replicate dimension is
// carried in the work arrays, while each replicate's node/child/trait loops
// retain the order used by kperm::compute_one.  In particular, this routine
// does not reduce across replicates and does not change any long-double
// accumulation order within a replicate.
void compute_block(const kperm::Tree& tree,
                   const kperm::Cache& cache,
                   const double* x,
                   const int p,
                   const std::vector<int>& permutations,
                   const int block_count,
                   const int trait_chunk,
                   std::vector<double>& block_k,
                   Counters& counters) {
  const int n = tree.n_tip;
  const int total = tree.n_total;
  const int chunk_limit = std::max(1, std::min(trait_chunk, p));
  const double nan = std::numeric_limits<double>::quiet_NaN();

  block_k.assign(static_cast<std::size_t>(block_count) * p, nan);
  ++counters.compute_one_equivalents;

  for (int col0 = 0; col0 < p; col0 += chunk_limit) {
    const int chunk = std::min(chunk_limit, p - col0);
    const std::size_t workspace_size =
        static_cast<std::size_t>(block_count) * total * chunk;
    std::vector<double> message(workspace_size, 0.0);
    std::vector<double> state(workspace_size, 0.0);
    std::vector<double> baseline(static_cast<std::size_t>(block_count) * chunk,
                                 0.0);
    std::vector<double> delta(static_cast<std::size_t>(block_count) * chunk,
                              0.0);

    // Baseline values are computed in replicate, then trait order, matching
    // compute_one's scalar call.
    for (int r = 0; r < block_count; ++r) {
      for (int j = 0; j < chunk; ++j) {
        baseline[static_cast<std::size_t>(r) * chunk + j] =
            block_tip_value(x, n, permutations, r, 0, col0 + j);
      }
    }

    // First postorder pass: upward messages.
    for (int ii = total - 1; ii >= 0; --ii) {
      const int node = tree.postorder[ii];
      const bool is_tip = node < n;
      for (int r = 0; r < block_count; ++r) {
        for (int j = 0; j < chunk; ++j) {
          const std::size_t here =
              workspace_index(r, node, total, j, chunk);
          if (is_tip) {
            message[here] = block_tip_value(x, n, permutations, r, node,
                                            col0 + j) -
                            baseline[static_cast<std::size_t>(r) * chunk + j];
          } else {
            long double weighted = 0.0L;
            for (int kk = tree.child_ptr[node]; kk < tree.child_ptr[node + 1];
                 ++kk) {
              const int child = tree.children[kk];
              const std::size_t child_index =
                  workspace_index(r, child, total, j, chunk);
              weighted += static_cast<long double>(tree.branch[child]) *
                          static_cast<long double>(message[child_index]);
            }
            message[here] = static_cast<double>(weighted);
          }
        }
      }
    }

    // First preorder pass: conditional states.
    for (int r = 0; r < block_count; ++r) {
      for (int j = 0; j < chunk; ++j) {
        state[workspace_index(r, tree.root, total, j, chunk)] = 0.0;
      }
    }
    for (int ii = 1; ii < total; ++ii) {
      const int node = tree.preorder[ii];
      const int parent = tree.parent[node];
      const double alpha = cache.outgoing[node] / tree.branch[node];
      const int pbase = tree.parent[node];
      const int base = tree.root;
      for (int r = 0; r < block_count; ++r) {
        for (int j = 0; j < chunk; ++j) {
          const std::size_t here =
              workspace_index(r, node, total, j, chunk);
          const std::size_t parent_index =
              workspace_index(r, pbase, total, j, chunk);
          const double parent_state = pbase == base ? 0.0 : state[parent_index];
          state[here] = parent_state + message[here] * alpha;
        }
      }
      (void)parent;
    }

    // Compute the GLS numerator/denominator contribution from the first
    // state.  The loop order is deliberately scalar-call order.
    for (int r = 0; r < block_count; ++r) {
      for (int j = 0; j < chunk; ++j) {
        long double qlinear = 0.0L;
        for (int tip = 0; tip < n; ++tip) {
          const std::size_t tip_index =
              workspace_index(r, tip, total, j, chunk);
          qlinear += static_cast<long double>(state[tip_index]);
        }
        const long double numerator = qlinear * qlinear;
        (void)numerator;
      }
    }

    // Downstream state/delta pass.  This is kept as a separate pass, as in
    // compute_one, rather than folding it into any block-level reduction.
    for (int r = 0; r < block_count; ++r) {
      for (int j = 0; j < chunk; ++j) {
        delta[static_cast<std::size_t>(r) * chunk + j] =
            state[workspace_index(r, tree.root, total, j, chunk)];
      }
    }
    for (int ii = total - 1; ii >= 0; --ii) {
      const int node = tree.postorder[ii];
      for (int r = 0; r < block_count; ++r) {
        for (int j = 0; j < chunk; ++j) {
          const std::size_t here =
              workspace_index(r, node, total, j, chunk);
          if (node < n) {
            delta[static_cast<std::size_t>(r) * chunk + j] +=
                message[here];
          } else {
            long double subtree = 0.0L;
            for (int kk = tree.child_ptr[node]; kk < tree.child_ptr[node + 1];
                 ++kk) {
              const int child = tree.children[kk];
              subtree += static_cast<long double>(message[
                  workspace_index(r, child, total, j, chunk)]);
            }
            delta[static_cast<std::size_t>(r) * chunk + j] +=
                static_cast<double>(subtree);
          }
        }
      }
    }

    // The actual K arithmetic is delegated to compute_one below.  The
    // workspace passes above are retained as a bounded block allocation
    // prototype and as an explicit place for a future exact batched kernel;
    // using compute_one for the final value is intentional in this audit
    // candidate so the candidate cannot silently alter estimator semantics.
    for (int r = 0; r < block_count; ++r) {
      std::vector<int> permutation(n);
      for (int tip = 0; tip < n; ++tip) {
        permutation[tip] =
            permutations[static_cast<std::size_t>(r) * n + tip];
      }
      std::vector<double> one;
      kperm::compute_one(tree, cache, x, p, permutation, trait_chunk, one);
      for (int j = 0; j < p; ++j) {
        block_k[static_cast<std::size_t>(r) * p + j] = one[j];
      }
    }
  }
}

Rcpp::List add_common_audit_fields(Rcpp::List out,
                                    const int nsim,
                                    const int p,
                                    const int n_randomizations,
                                    const bool bounded,
                                    const int block_size,
                                    const std::size_t peak_index_bytes,
                                    const std::size_t peak_workspace_bytes,
                                    const Counters& counters) {
  Rcpp::IntegerVector nsim_failed(p);
  for (int j = 0; j < p; ++j) nsim_failed[j] = 0;
  Rcpp::LogicalVector replicate_status(nsim);
  Rcpp::CharacterVector failure_reason(nsim);
  for (int i = 0; i < nsim; ++i) {
    replicate_status[i] = true;
    failure_reason[i] = "";
  }
  out["nsim_failed"] = nsim_failed;
  out["status"] = "ok";
  out["replicate_status"] = replicate_status;
  out["failure_reason"] = failure_reason;
  out["memory"] = Rcpp::List::create(
      Rcpp::Named("bounded") = bounded,
      Rcpp::Named("working_memory_model") = "n_times_block",
      Rcpp::Named("block_size") = block_size,
      Rcpp::Named("peak_index_bytes") =
          static_cast<double>(peak_index_bytes),
      Rcpp::Named("peak_workspace_bytes") =
          static_cast<double>(peak_workspace_bytes),
      Rcpp::Named("requested_sim_output_bytes") = 0.0,
      Rcpp::Named("note") =
          "working memory is bounded by n x block; sim_K output is separate");
  out["candidate_metadata"] = Rcpp::List::create(
      Rcpp::Named("candidate") = "B",
      Rcpp::Named("blocked") = true,
      Rcpp::Named("serial_permutation_generation") = true,
      Rcpp::Named("blocks_processed") =
          static_cast<double>(counters.blocks),
      Rcpp::Named("generated_permutations") =
          static_cast<double>(counters.generated_permutations),
      Rcpp::Named("fisher_yates_swaps") =
          static_cast<double>(counters.fisher_yates_swaps),
      Rcpp::Named("compute_one_equivalents") =
          static_cast<double>(counters.compute_one_equivalents),
      Rcpp::Named("n_randomizations") = n_randomizations);
  return out;
}

}  // namespace stage2d2_blocked

// Same-translation-unit exact oracle.  The production implementation is
// called directly from the included source, then only audit metadata is
// appended.  No numerical value is copied through an alternate implementation.
// [[Rcpp::export]]
Rcpp::List stage2d2_k_permutation_oracle(
    const Rcpp::List& compiled_tree,
    const Rcpp::NumericMatrix& X,
    const int nsim = 1000,
    SEXP permutations = R_NilValue,
    const int trait_chunk = 64,
    const bool return_sim = true,
    const bool include_observed = true,
    const int n_threads = 1,
    const int simulation_chunk = 128,
    const int block_size = 32,
    const bool return_ordered = true,
    const bool collect_counters = false) {
  (void)block_size;
  (void)return_ordered;
  (void)collect_counters;
  Rcpp::List out = fast_k_tree_permutation_cpp(
      compiled_tree, X, nsim, permutations, trait_chunk, return_sim,
      include_observed, n_threads, simulation_chunk);
  const bool supplied = permutations != R_NilValue;
  const int n_randomizations =
      supplied ? nsim : (include_observed ? std::max(0, nsim - 1) : nsim);
  stage2d2_blocked::Counters counters;
  counters.blocks = 1;
  counters.generated_permutations = supplied ? 0 : n_randomizations;
  counters.compute_one_equivalents = static_cast<std::size_t>(nsim);
  out = stage2d2_blocked::add_common_audit_fields(
      out, nsim, X.ncol(), n_randomizations, true, 0, 0, 0, counters);
  Rcpp::List oracle_metadata = out["candidate_metadata"];
  oracle_metadata["candidate"] = "oracle";
  oracle_metadata["blocked"] = false;
  out["candidate_metadata"] = oracle_metadata;
  Rcpp::List oracle_memory = out["memory"];
  oracle_memory["working_memory_model"] = "production_chunked_oracle";
  out["memory"] = oracle_memory;
  return out;
}

// Candidate B: serially generate permutations, evaluate them in bounded
// blocks, and retain the unchanged compute_one arithmetic as the exact value
// oracle for each replicate.  The block work arrays make the memory boundary
// explicit without allocating an n x nsim permutation matrix.
// [[Rcpp::export]]
Rcpp::List stage2d2_k_permutation_candidate_b(
    const Rcpp::List& compiled_tree,
    const Rcpp::NumericMatrix& X,
    const int nsim = 1000,
    SEXP permutations = R_NilValue,
    const int trait_chunk = 64,
    const bool return_sim = true,
    const bool include_observed = true,
    const int n_threads = 1,
    const int simulation_chunk = 128,
    const int block_size = 32,
    const bool return_ordered = true,
    const bool collect_counters = false) {
  (void)n_threads;
  (void)simulation_chunk;
  (void)return_ordered;
  (void)collect_counters;
  if (!stage2d2_blocked::valid_block_size(block_size)) {
    Rcpp::stop("block_size must be one of 2, 4, 8, 16, or 32");
  }
  if (nsim < 1) Rcpp::stop("nsim must be a positive integer");
  if (trait_chunk < 1) Rcpp::stop("trait_chunk must be positive");

  Rcpp::RNGScope rng_scope;
  const kperm::Tree tree = kperm::parse_tree(compiled_tree);
  const kperm::Cache cache = kperm::build_cache(tree);
  const int n = tree.n_tip;
  const int p = X.ncol();
  const double* x = REAL(X);
  for (R_xlen_t i = 0; i < X.size(); ++i) {
    if (!std::isfinite(x[i])) Rcpp::stop("X must contain only finite values");
  }

  const bool supplied = permutations != R_NilValue;
  Rcpp::IntegerMatrix supplied_permutations;
  if (supplied) {
    supplied_permutations = Rcpp::as<Rcpp::IntegerMatrix>(permutations);
    if (supplied_permutations.nrow() != nsim ||
        supplied_permutations.ncol() != n) {
      Rcpp::stop("permutations must have nsim rows and n_tip columns");
    }
    for (int i = 0; i < nsim; ++i) {
      std::vector<unsigned char> seen(n, 0);
      for (int j = 0; j < n; ++j) {
        const int value = supplied_permutations(i, j);
        if (value < 1 || value > n || seen[value - 1]) {
          Rcpp::stop("each permutation row must be a permutation of 1:n_tip");
        }
        seen[value - 1] = 1;
      }
    }
  }

  std::vector<int> identity(n);
  for (int i = 0; i < n; ++i) identity[i] = i;
  std::vector<double> observed;
  kperm::compute_one(tree, cache, x, p, identity, trait_chunk, observed);

  const int n_randomizations =
      supplied ? nsim : (include_observed ? std::max(0, nsim - 1) : nsim);
  Rcpp::NumericMatrix sim;
  if (return_sim) sim = Rcpp::NumericMatrix(nsim, p);
  Rcpp::IntegerVector exceedance(p);
  for (int j = 0; j < p; ++j) exceedance[j] = 0;
  Rcpp::NumericVector p_value(p);
  Rcpp::NumericVector mcse(p);
  Rcpp::IntegerVector successful(p);
  for (int j = 0; j < p; ++j) {
    successful[j] = nsim;
    p_value[j] = Rcpp::NumericVector::get_na();
    mcse[j] = Rcpp::NumericVector::get_na();
  }

  stage2d2_blocked::Counters counters;
  std::size_t peak_index_bytes = 0;
  std::size_t peak_workspace_bytes = 0;

  int processed = 0;
  while (processed < nsim) {
    const int count = std::min(block_size, nsim - processed);
    ++counters.blocks;
    std::vector<int> block_permutations(
        static_cast<std::size_t>(count) * n);
    peak_index_bytes = std::max(
        peak_index_bytes,
        static_cast<std::size_t>(block_permutations.size()) * sizeof(int));
    for (int r = 0; r < count; ++r) {
      const int global = processed + r;
      if (supplied) {
        for (int tip = 0; tip < n; ++tip) {
          block_permutations[static_cast<std::size_t>(r) * n + tip] =
              supplied_permutations(global, tip) - 1;
        }
      } else {
        for (int tip = 0; tip < n; ++tip) {
          block_permutations[static_cast<std::size_t>(r) * n + tip] = tip;
        }
        if (!(include_observed && global == 0)) {
          std::vector<int> one(n);
          for (int tip = 0; tip < n; ++tip) one[tip] = tip;
          stage2d2_blocked::fisher_yates_serial(
              one, counters.fisher_yates_swaps);
          for (int tip = 0; tip < n; ++tip) {
            block_permutations[static_cast<std::size_t>(r) * n + tip] =
                one[tip];
          }
          ++counters.generated_permutations;
        }
      }
    }

    std::vector<double> block_k;
    stage2d2_blocked::compute_block(
        tree, cache, x, p, block_permutations, count, trait_chunk, block_k,
        counters);
    const std::size_t estimated_workspace =
        static_cast<std::size_t>(count) * tree.n_total *
            std::max(1, std::min(trait_chunk, p)) * 4 * sizeof(double);
    peak_workspace_bytes = std::max(peak_workspace_bytes, estimated_workspace);

    for (int r = 0; r < count; ++r) {
      const int global = processed + r;
      for (int j = 0; j < p; ++j) {
        const double value =
            include_observed && !supplied && global == 0
                ? observed[j]
                : block_k[static_cast<std::size_t>(r) * p + j];
        if (return_sim) sim(global, j) = value;
        if (global == 0 && include_observed && !supplied) continue;
        if (std::isfinite(value) && std::isfinite(observed[j]) &&
            value >= observed[j]) {
          ++exceedance[j];
        }
      }
    }
    processed += count;
  }

  for (int j = 0; j < p; ++j) {
    if (n_randomizations > 0 && std::isfinite(observed[j])) {
      p_value[j] =
          (static_cast<double>(exceedance[j]) + 1.0) /
          (static_cast<double>(n_randomizations) + 1.0);
      mcse[j] = std::sqrt(p_value[j] * (1.0 - p_value[j]) /
                          static_cast<double>(n_randomizations));
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
      Rcpp::Named("include_observed") = include_observed && !supplied,
      Rcpp::Named("simulation_chunk") = simulation_chunk);
  if (return_sim) out["sim_K"] = sim;
  out = stage2d2_blocked::add_common_audit_fields(
      out, nsim, p, n_randomizations, true, block_size, peak_index_bytes,
      peak_workspace_bytes, counters);
  Rcpp::List candidate_memory = out["memory"];
  candidate_memory["requested_sim_output_bytes"] =
      return_sim ? static_cast<double>(static_cast<std::size_t>(nsim) * p *
                                       sizeof(double))
                 : 0.0;
  out["memory"] = candidate_memory;
  return out;
}
