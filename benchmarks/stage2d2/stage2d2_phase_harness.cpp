// Private Stage 2D2 phase-decomposition harness for the K null engine.
//
// This file is audit-only.  It includes the production translation unit so
// that the measured evaluator is exactly kperm::compute_one() and the
// generation path is exactly kperm::fisher_yates().  Nothing in this file is
// linked into the package DLL or exposed through the package namespace.

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

// Keep the source path behind a macro.  Rcpp::sourceCpp otherwise treats a
// literal .cpp include as a second compilation unit, which would duplicate
// every production symbol at link time.  The compiler still expands this
// include in the current translation unit.
#define STAGE2D2_PRODUCTION_CPP "../../src/k_permutation.cpp"
#include STAGE2D2_PRODUCTION_CPP
#undef STAGE2D2_PRODUCTION_CPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#ifdef _OPENMP
# include <omp.h>
#endif

namespace stage2d2 {

typedef std::chrono::steady_clock stage_clock;

std::string checksum_text(const std::uint64_t value) {
  const char* digits = "0123456789abcdef";
  std::string out(16, '0');
  std::uint64_t current = value;
  for (int i = 15; i >= 0; --i) {
    out[static_cast<std::size_t>(i)] =
      digits[static_cast<std::size_t>(current & 0xfULL)];
    current >>= 4U;
  }
  return out;
}

std::uint64_t mix_double(const std::uint64_t current, const double value) {
  const std::uint64_t scaled = std::isfinite(value)
    ? static_cast<std::uint64_t>(std::llround(std::fabs(value) * 1e9))
    : 0xffffffffffffffffULL;
  return current * 1469598103934665603ULL + scaled;
}

double elapsed_seconds(const stage_clock::time_point start,
                       const stage_clock::time_point end) {
  return std::chrono::duration_cast<std::chrono::duration<double> >(
    end - start
  ).count();
}

int effective_threads(const int requested) {
  if (requested < 1) Rcpp::stop("n_threads must be positive.");
#ifdef _OPENMP
  return std::max(1, std::min(requested, omp_get_max_threads()));
#else
  return 1;
#endif
}

bool openmp_compiled() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}

int openmp_max_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

struct TimingResult {
  double elapsed_s;
  std::uint64_t checksum;
  std::uint64_t rng_draws;
  std::uint64_t fisher_yates_swaps;
  std::vector<int> trace;
  int trace_rows;
};

// Keep a small trace so the R harness can replay the first exact generated
// permutations without materializing the full n x nsim null matrix.
void append_trace(std::vector<int>& trace, const int trace_rows,
                  const int n, const int row,
                  const std::vector<int>& perm) {
  if (row >= trace_rows) return;
  const std::size_t begin = static_cast<std::size_t>(row) *
    static_cast<std::size_t>(n);
  std::copy(perm.begin(), perm.end(), trace.begin() + begin);
}

TimingResult run_generation(const int n, const int nsim,
                            const bool include_observed,
                            const int trace_rows_requested) {
  if (n < 2) Rcpp::stop("n must be at least two.");
  if (nsim < 1) Rcpp::stop("nsim must be positive.");
  const int trace_rows = std::max(0, std::min(trace_rows_requested, nsim));
  std::vector<int> identity(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) identity[static_cast<std::size_t>(i)] = i;
  std::vector<int> trace(static_cast<std::size_t>(trace_rows) *
                         static_cast<std::size_t>(n), 0);
  std::uint64_t checksum = 0;
  std::uint64_t draws = 0;
  std::uint64_t swaps = 0;

  // Rcpp's RNGScope is deliberately inside this private exported call.  The
  // generated permutations therefore consume the same R RNG stream as the
  // production entry point, including the identity-first convention.
  Rcpp::RNGScope rng_scope;
  const stage_clock::time_point started = stage_clock::now();
  for (int i = 0; i < nsim; ++i) {
    // The production engine constructs a fresh n-element vector for every
    // replicate, resets it to identity, and then shuffles in place.
    std::vector<int> perm(static_cast<std::size_t>(n));
    std::copy(identity.begin(), identity.end(), perm.begin());
    const bool identity_row = include_observed && i == 0;
    if (!identity_row) {
      // Call the production helper itself.  Its current contract is one
      // R::runif(0, pos + 1) draw and one swap for every pos = n-1,...,1.
      kperm::fisher_yates(perm);
      draws += static_cast<std::uint64_t>(n - 1);
      swaps += static_cast<std::uint64_t>(n - 1);
    }
    append_trace(trace, trace_rows, n, i, perm);
    checksum = checksum * 1469598103934665603ULL +
      static_cast<std::uint64_t>(perm[static_cast<std::size_t>(i % n)] + 1);
  }
  const stage_clock::time_point finished = stage_clock::now();
  TimingResult out;
  out.elapsed_s = elapsed_seconds(started, finished);
  out.checksum = checksum;
  out.rng_draws = draws;
  out.fisher_yates_swaps = swaps;
  out.trace = trace;
  out.trace_rows = trace_rows;
  return out;
}

// Build a bounded, deterministic family of valid controlled permutations.
// The family is intentionally independent of R's RNG; it is used only to
// time the null-K evaluator after permutation preparation has been excluded.
std::vector<std::vector<int> > fixed_controlled_permutations(
    const int n, const int requested) {
  const int count = std::max(1, std::min(requested, 8));
  std::vector<std::vector<int> > out(static_cast<std::size_t>(count));
  for (int pattern = 0; pattern < count; ++pattern) {
    std::vector<int>& perm = out[static_cast<std::size_t>(pattern)];
    perm.resize(static_cast<std::size_t>(n));
    const int shift = pattern / 2;
    const bool reversed = pattern % 2 == 1;
    for (int r = 0; r < n; ++r) {
      const int base = reversed ? (n - 1 - r) : r;
      perm[static_cast<std::size_t>(r)] = (base + shift) % n;
    }
  }
  return out;
}

struct EvaluationResult {
  double elapsed_s;
  std::uint64_t checksum;
  std::uint64_t trait_reads;
};

EvaluationResult run_controlled_evaluation(
    const kperm::Tree& tree, const kperm::Cache& cache,
    const double* x, const int ncol, const int nsim,
    const int trait_chunk, const int requested_threads,
    const int controlled_count) {
  const int threads = effective_threads(requested_threads);
  const std::vector<std::vector<int> > controlled =
    fixed_controlled_permutations(tree.n_tip, controlled_count);
  const int pattern_count = static_cast<int>(controlled.size());
  std::vector<std::uint64_t> checksums(static_cast<std::size_t>(threads), 0);
  std::vector<std::uint64_t> reads(static_cast<std::size_t>(threads), 0);
  const stage_clock::time_point started = stage_clock::now();

  if (threads > 1) {
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
    {
      const int tid = omp_get_thread_num();
      std::vector<double> kval;
#pragma omp for schedule(static)
      for (int i = 0; i < nsim; ++i) {
        const std::vector<int>& perm = controlled[
          static_cast<std::size_t>(i % pattern_count)
        ];
        compute_one(tree, cache, x, ncol, perm, trait_chunk, kval);
        const double value = kval[static_cast<std::size_t>(i % ncol)];
        checksums[static_cast<std::size_t>(tid)] =
          checksums[static_cast<std::size_t>(tid)] *
          1469598103934665603ULL +
          static_cast<std::uint64_t>(std::llround(
            std::fabs(value) * 1000000000.0
          ));
        reads[static_cast<std::size_t>(tid)] +=
          static_cast<std::uint64_t>(tree.n_tip) *
          static_cast<std::uint64_t>(5) *
          static_cast<std::uint64_t>(ncol) +
          static_cast<std::uint64_t>(ncol);
      }
    }
#else
    (void)threads;
#endif
  } else {
    std::vector<double> kval;
    for (int i = 0; i < nsim; ++i) {
      const std::vector<int>& perm = controlled[
        static_cast<std::size_t>(i % pattern_count)
      ];
      compute_one(tree, cache, x, ncol, perm, trait_chunk, kval);
      const double value = kval[static_cast<std::size_t>(i % ncol)];
      checksums[0] = checksums[0] * 1469598103934665603ULL +
        static_cast<std::uint64_t>(std::llround(
          std::fabs(value) * 1000000000.0
        ));
      reads[0] += static_cast<std::uint64_t>(tree.n_tip) *
        static_cast<std::uint64_t>(5) * static_cast<std::uint64_t>(ncol) +
        static_cast<std::uint64_t>(ncol);
    }
  }
  const stage_clock::time_point finished = stage_clock::now();
  EvaluationResult out;
  out.elapsed_s = elapsed_seconds(started, finished);
  out.checksum = 0;
  out.trait_reads = 0;
  for (int t = 0; t < threads; ++t) {
    out.checksum ^= checksums[static_cast<std::size_t>(t)] +
      static_cast<std::uint64_t>(t + 1);
    out.trait_reads += reads[static_cast<std::size_t>(t)];
  }
  return out;
}

struct FullResult {
  double elapsed_s;
  double observed_s;
  double generation_and_chunk_s;
  double evaluation_s;
  double accounting_s;
  std::uint64_t generation_checksum;
  std::uint64_t evaluation_checksum;
  std::uint64_t rng_draws;
  std::uint64_t fisher_yates_swaps;
  std::vector<int> trace;
  int trace_rows;
  std::size_t chunk_bytes;
  std::size_t generation_bytes;
  std::size_t worker_bytes;
  std::size_t peak_index_bytes;
};

FullResult run_full_pipeline(
    const kperm::Tree& tree, const kperm::Cache& cache,
    const double* x, const int ncol, const int nsim,
    const int trait_chunk, const bool include_observed,
    const int requested_threads, const int simulation_chunk,
    const int trace_rows_requested) {
  const int threads = effective_threads(requested_threads);
  const int chunk_limit = std::min(simulation_chunk, nsim);
  const int n = tree.n_tip;
  const int trace_rows = std::max(0, std::min(trace_rows_requested, nsim));
  const stage_clock::time_point full_started = stage_clock::now();
  std::vector<int> identity(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) identity[static_cast<std::size_t>(i)] = i;
  std::vector<int> trace(static_cast<std::size_t>(trace_rows) *
                         static_cast<std::size_t>(n), 0);
  std::vector<int> chunk_perms(
    static_cast<std::size_t>(chunk_limit) * static_cast<std::size_t>(n)
  );
  std::uint64_t checksum = 0;
  std::uint64_t evaluation_checksum = 0;
  std::uint64_t draws = 0;
  std::uint64_t swaps = 0;

  // Production computes the observed identity once before entering the null
  // loop.  It is part of this full-pipeline total, but gets its own field so
  // the null-work budget is transparent.
  std::vector<double> observed;
  const stage_clock::time_point observed_started = stage_clock::now();
  compute_one(tree, cache, x, ncol, identity, trait_chunk, observed);
  const stage_clock::time_point observed_finished = stage_clock::now();

  std::vector<char> valid_observed(static_cast<std::size_t>(ncol), 1);
  for (int j = 0; j < ncol; ++j) valid_observed[
    static_cast<std::size_t>(j)
  ] = std::isfinite(observed[static_cast<std::size_t>(j)]) ? 1 : 0;
  std::vector<double> exceedance(static_cast<std::size_t>(ncol), 0.0);

  Rcpp::RNGScope rng_scope;
  double generation_and_chunk_s = 0.0;
  double evaluation_s = 0.0;
  for (int first = 0; first < nsim; first += chunk_limit) {
    const int count = std::min(chunk_limit, nsim - first);
    const stage_clock::time_point generation_started = stage_clock::now();
    for (int local_i = 0; local_i < count; ++local_i) {
      std::vector<int> perm(static_cast<std::size_t>(n));
      std::copy(identity.begin(), identity.end(), perm.begin());
      const int global_i = first + local_i;
      const bool identity_row = include_observed && global_i == 0;
      if (!identity_row) {
        kperm::fisher_yates(perm);
        draws += static_cast<std::uint64_t>(n - 1);
        swaps += static_cast<std::uint64_t>(n - 1);
      }
      append_trace(trace, trace_rows, n, global_i, perm);
      checksum = checksum * 1469598103934665603ULL +
        static_cast<std::uint64_t>(perm[static_cast<std::size_t>(global_i % n)] + 1);
      const std::size_t begin = static_cast<std::size_t>(local_i) *
        static_cast<std::size_t>(n);
      std::copy(perm.begin(), perm.end(), chunk_perms.begin() + begin);
    }
    generation_and_chunk_s += elapsed_seconds(
      generation_started, stage_clock::now()
    );

    const stage_clock::time_point evaluation_started = stage_clock::now();
    if (threads > 1) {
#ifdef _OPENMP
      std::vector<std::uint64_t> local_checksums(
        static_cast<std::size_t>(threads), 0
      );
      std::vector<std::vector<double> > local_exceedance(
        static_cast<std::size_t>(threads),
        std::vector<double>(static_cast<std::size_t>(ncol), 0.0)
      );
#pragma omp parallel num_threads(threads)
      {
        const int tid = omp_get_thread_num();
        std::vector<int> perm(static_cast<std::size_t>(n));
        std::vector<double> kval;
#pragma omp for schedule(static)
        for (int local_i = 0; local_i < count; ++local_i) {
          const std::size_t begin = static_cast<std::size_t>(local_i) *
            static_cast<std::size_t>(n);
          std::copy(chunk_perms.begin() + begin,
                    chunk_perms.begin() + begin + static_cast<std::size_t>(n),
                    perm.begin());
          compute_one(tree, cache, x, ncol, perm, trait_chunk, kval);
          for (int j = 0; j < ncol; ++j) {
            const double value = kval[static_cast<std::size_t>(j)];
            local_checksums[static_cast<std::size_t>(tid)] = mix_double(
              local_checksums[static_cast<std::size_t>(tid)], value
            );
            if (valid_observed[static_cast<std::size_t>(j)] &&
                fastphylosig::inclusive_upper_tail(
                  value, observed[static_cast<std::size_t>(j)]
                )) local_exceedance[static_cast<std::size_t>(tid)][
                  static_cast<std::size_t>(j)
                ] += 1.0;
          }
        }
      }
      for (int tid = 0; tid < threads; ++tid) {
        evaluation_checksum ^= local_checksums[static_cast<std::size_t>(tid)] +
          static_cast<std::uint64_t>(tid + 1);
        for (int j = 0; j < ncol; ++j) exceedance[
          static_cast<std::size_t>(j)
        ] += local_exceedance[static_cast<std::size_t>(tid)][
          static_cast<std::size_t>(j)
        ];
      }
#else
      (void)threads;
#endif
    } else {
      std::vector<int> perm(static_cast<std::size_t>(n));
      std::vector<double> kval;
      for (int local_i = 0; local_i < count; ++local_i) {
        const std::size_t begin = static_cast<std::size_t>(local_i) *
          static_cast<std::size_t>(n);
        std::copy(chunk_perms.begin() + begin,
                  chunk_perms.begin() + begin + static_cast<std::size_t>(n),
                  perm.begin());
        compute_one(tree, cache, x, ncol, perm, trait_chunk, kval);
        for (int j = 0; j < ncol; ++j) {
          const double value = kval[static_cast<std::size_t>(j)];
          evaluation_checksum = mix_double(evaluation_checksum, value);
          if (valid_observed[static_cast<std::size_t>(j)] &&
              fastphylosig::inclusive_upper_tail(
                value, observed[static_cast<std::size_t>(j)]
              )) exceedance[static_cast<std::size_t>(j)] += 1.0;
        }
      }
    }
    evaluation_s += elapsed_seconds(evaluation_started, stage_clock::now());
  }

  // The production engine then constructs P/MCSE/accounting vectors.  This
  // accounting is intentionally kept outside the null evaluator and included
  // only in the full-pipeline total.
  const stage_clock::time_point accounting_started = stage_clock::now();
  std::vector<double> p_value(static_cast<std::size_t>(ncol), 0.0);
  std::vector<double> mcse(static_cast<std::size_t>(ncol), 0.0);
  const int n_random = include_observed ? std::max(0, nsim - 1) : nsim;
  for (int j = 0; j < ncol; ++j) {
    p_value[static_cast<std::size_t>(j)] = exceedance[
      static_cast<std::size_t>(j)
    ] / static_cast<double>(nsim);
    if (include_observed && n_random > 0) {
      const double q = std::max(0.0, std::min(1.0,
        (exceedance[static_cast<std::size_t>(j)] - 1.0) /
        static_cast<double>(n_random))
      );
      mcse[static_cast<std::size_t>(j)] = std::sqrt(
        q * (1.0 - q) * static_cast<double>(n_random)
      ) / static_cast<double>(nsim);
    } else if (n_random > 0) {
      const double q = p_value[static_cast<std::size_t>(j)];
      mcse[static_cast<std::size_t>(j)] = std::sqrt(
        q * (1.0 - q) / static_cast<double>(n_random)
      );
    } else {
      mcse[static_cast<std::size_t>(j)] = NA_REAL;
    }
  }
  const stage_clock::time_point accounting_finished = stage_clock::now();

  FullResult out;
  out.elapsed_s = elapsed_seconds(full_started, accounting_finished);
  out.observed_s = elapsed_seconds(observed_started, observed_finished);
  out.generation_and_chunk_s = generation_and_chunk_s;
  out.evaluation_s = evaluation_s;
  out.accounting_s = elapsed_seconds(accounting_started, accounting_finished);
  out.generation_checksum = checksum;
  out.evaluation_checksum = evaluation_checksum;
  out.rng_draws = draws;
  out.fisher_yates_swaps = swaps;
  out.trace = trace;
  out.trace_rows = trace_rows;
  out.chunk_bytes = chunk_perms.size() * sizeof(int);
  out.generation_bytes = static_cast<std::size_t>(n) * sizeof(int);
  out.worker_bytes = static_cast<std::size_t>(threads) *
    static_cast<std::size_t>(n) * sizeof(int);
  out.peak_index_bytes = out.chunk_bytes + out.generation_bytes +
    out.worker_bytes;
  return out;
}

Rcpp::IntegerMatrix as_trace_matrix(const std::vector<int>& trace,
                                    const int rows, const int n) {
  Rcpp::IntegerMatrix out(rows, n);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < n; ++j) out(i, j) = trace[
      static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
      static_cast<std::size_t>(j)
    ] + 1;
  }
  return out;
}

} // namespace stage2d2

// [[Rcpp::export]]
Rcpp::List stage2d2_phase_profile(
    const Rcpp::List& compiled_tree,
    const Rcpp::NumericMatrix& X,
    const int nsim = 999,
    const int trait_chunk = 1,
    const bool include_observed = true,
    const int n_threads = 1,
    const int simulation_chunk = 128,
    const std::string phase = "full",
    const int controlled_count = 8,
    const int trace_rows = 3) {
  using namespace stage2d2;
  if (nsim < 1) Rcpp::stop("nsim must be positive.");
  if (trait_chunk < 1) Rcpp::stop("trait_chunk must be positive.");
  if (simulation_chunk < 1) {
    Rcpp::stop("simulation_chunk must be positive.");
  }
  if (controlled_count < 1) Rcpp::stop("controlled_count must be positive.");
  if (trace_rows < 0) Rcpp::stop("trace_rows must be non-negative.");
  if (phase != "generation" && phase != "evaluation" && phase != "full") {
    Rcpp::stop("phase must be generation, evaluation, or full.");
  }

  // Parse and cache once before the phase timer.  This harness measures the
  // prepared null engine; tree preparation is a separate production budget.
  const kperm::Tree tree = kperm::parse_tree(compiled_tree);
  if (X.nrow() != tree.n_tip || X.ncol() < 1) {
    Rcpp::stop("X must have one row per compiled tree tip and one trait.");
  }
  const kperm::Cache cache = kperm::build_cache(tree);
  const double* x = REAL(X);
  for (R_xlen_t i = 0; i < X.size(); ++i) {
    if (!std::isfinite(x[i])) Rcpp::stop("X must contain finite values.");
  }

  const int threads_eff = effective_threads(n_threads);
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("phase") = phase,
    Rcpp::Named("n") = tree.n_tip,
    Rcpp::Named("n_total") = tree.n_total,
    Rcpp::Named("traits") = X.ncol(),
    Rcpp::Named("nsim") = nsim,
    Rcpp::Named("include_observed") = include_observed,
    Rcpp::Named("n_threads_requested") = n_threads,
    Rcpp::Named("n_threads_effective") = threads_eff,
    Rcpp::Named("openmp_compiled") = openmp_compiled(),
    Rcpp::Named("openmp_max_threads") = openmp_max_threads(),
    Rcpp::Named("simulation_chunk") = simulation_chunk,
    Rcpp::Named("controlled_count") = std::max(1, std::min(controlled_count, 8)),
    Rcpp::Named("trace_rows_requested") = trace_rows,
    Rcpp::Named("evaluator_source") = "src/k_permutation.cpp::kperm::compute_one",
    Rcpp::Named("rng_source") = "src/k_permutation.cpp::R::runif + Fisher-Yates"
  );

  if (phase == "generation") {
    const TimingResult result = run_generation(
      tree.n_tip, nsim, include_observed, trace_rows
    );
    out["elapsed_s"] = result.elapsed_s;
    out["generation_s"] = result.elapsed_s;
    out["checksum"] = checksum_text(result.checksum);
    out["rng_draws"] = static_cast<double>(result.rng_draws);
    out["fisher_yates_swaps"] = static_cast<double>(result.fisher_yates_swaps);
    out["trace"] = as_trace_matrix(result.trace, result.trace_rows, tree.n_tip);
    out["trace_rows"] = result.trace_rows;
    out["peak_index_bytes"] = static_cast<double>(tree.n_tip * sizeof(int));
    return out;
  }

  if (phase == "evaluation") {
    const EvaluationResult result = run_controlled_evaluation(
      tree, cache, x, X.ncol(), nsim, trait_chunk, n_threads, controlled_count
    );
    out["elapsed_s"] = result.elapsed_s;
    out["evaluation_s"] = result.elapsed_s;
    out["checksum"] = checksum_text(result.checksum);
    out["trait_reads"] = static_cast<double>(result.trait_reads);
    out["peak_index_bytes"] = static_cast<double>(
      std::max(1, std::min(controlled_count, 8)) * tree.n_tip * sizeof(int)
    );
    return out;
  }

  const FullResult result = run_full_pipeline(
    tree, cache, x, X.ncol(), nsim, trait_chunk, include_observed,
    n_threads, simulation_chunk, trace_rows
  );
  out["elapsed_s"] = result.elapsed_s;
  out["full_pipeline_s"] = result.elapsed_s;
  out["observed_s"] = result.observed_s;
  out["generation_and_chunk_s"] = result.generation_and_chunk_s;
  out["evaluation_s"] = result.evaluation_s;
  out["accounting_s"] = result.accounting_s;
  out["generation_checksum"] = checksum_text(result.generation_checksum);
  out["evaluation_checksum"] = checksum_text(result.evaluation_checksum);
  out["rng_draws"] = static_cast<double>(result.rng_draws);
  out["fisher_yates_swaps"] = static_cast<double>(result.fisher_yates_swaps);
  out["trace"] = as_trace_matrix(result.trace, result.trace_rows, tree.n_tip);
  out["trace_rows"] = result.trace_rows;
  out["chunk_bytes"] = static_cast<double>(result.chunk_bytes);
  out["generation_bytes"] = static_cast<double>(result.generation_bytes);
  out["worker_bytes"] = static_cast<double>(result.worker_bytes);
  out["peak_index_bytes"] = static_cast<double>(result.peak_index_bytes);
  return out;
}
