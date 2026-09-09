// Private Stage 2D1 prototype for the K permutation/null engine.
//
// This file is intentionally outside src/.  It is compiled with
// Rcpp::sourceCpp() by the audit harness and is not a package binding.  The
// oracle path mirrors src/k_permutation.cpp. A and B change only the bounded
// permutation-index transport described in the Stage 2D1 protocol. The tree
// parse, cache construction, K arithmetic, RNG calls, and tail rule are kept
// in the same order as production. Candidate C is intentionally absent:
// production already reads X through permutation indices and has no complete
// permuted-trait gather to fuse.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#ifdef _OPENMP
# include <omp.h>
#endif

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

namespace stage2d1 {

struct Counters {
  std::uint64_t null_replicates = 0;
  std::uint64_t index_buffer_allocations = 0;
  std::uint64_t index_initialization_elements = 0;
  std::uint64_t fisher_yates_swaps = 0;
  std::uint64_t generation_to_chunk_elements = 0;
  std::uint64_t chunk_to_worker_elements = 0;
  std::uint64_t controlled_to_worker_elements = 0;
  std::uint64_t index_reads = 0;
  std::uint64_t trait_reads = 0;
  std::uint64_t compute_workspace_allocations = 0;
  std::uint64_t compute_workspace_bytes = 0;
  std::uint64_t trait_gather_allocations = 0;
  std::size_t compute_workspace_peak_bytes = 0;

  void add(const Counters& rhs) {
    null_replicates += rhs.null_replicates;
    index_buffer_allocations += rhs.index_buffer_allocations;
    index_initialization_elements += rhs.index_initialization_elements;
    fisher_yates_swaps += rhs.fisher_yates_swaps;
    generation_to_chunk_elements += rhs.generation_to_chunk_elements;
    chunk_to_worker_elements += rhs.chunk_to_worker_elements;
    controlled_to_worker_elements += rhs.controlled_to_worker_elements;
    index_reads += rhs.index_reads;
    trait_reads += rhs.trait_reads;
    compute_workspace_allocations += rhs.compute_workspace_allocations;
    compute_workspace_bytes += rhs.compute_workspace_bytes;
    trait_gather_allocations += rhs.trait_gather_allocations;
    compute_workspace_peak_bytes = std::max(
      compute_workspace_peak_bytes, rhs.compute_workspace_peak_bytes
    );
  }
};

struct Tree {
  int n_tip;
  int n_total;
  int root;
  std::vector<int> parent;
  std::vector<int> child_ptr;
  std::vector<int> children;
  std::vector<int> preorder;
  std::vector<int> postorder;
  std::vector<double> branch;
  std::vector<double> root_distance;
};

struct Cache {
  std::vector<double> outgoing;
  std::vector<double> aggregate;
  double sum_inv;
  double normalization;
};

bool has(const Rcpp::List& x, const char* name) {
  return x.containsElementNamed(name) && x[name] != R_NilValue;
}

int scalar_int(SEXP x, const std::string& label, const int min_value) {
  if (x == R_NilValue || Rf_xlength(x) != 1 ||
      (TYPEOF(x) != INTSXP && TYPEOF(x) != REALSXP)) {
    Rcpp::stop(label + " must be a scalar integer.");
  }
  const double value = TYPEOF(x) == INTSXP
    ? static_cast<double>(INTEGER(x)[0]) : REAL(x)[0];
  if (!std::isfinite(value) || std::floor(value) != value ||
      value < static_cast<double>(min_value) ||
      value > static_cast<double>(std::numeric_limits<int>::max())) {
    Rcpp::stop(label + " must be a finite integer.");
  }
  return static_cast<int>(value);
}

std::vector<int> ints(SEXP x, const R_xlen_t expected,
                      const std::string& label) {
  if (x == R_NilValue || Rf_xlength(x) != expected ||
      (TYPEOF(x) != INTSXP && TYPEOF(x) != REALSXP)) {
    Rcpp::stop(label + " has an incompatible length or type.");
  }
  std::vector<int> out(static_cast<std::size_t>(expected));
  for (R_xlen_t i = 0; i < expected; ++i) {
    const double v = TYPEOF(x) == INTSXP
      ? static_cast<double>(INTEGER(x)[i]) : REAL(x)[i];
    if (!std::isfinite(v) || std::floor(v) != v ||
        v < static_cast<double>(std::numeric_limits<int>::min()) ||
        v > static_cast<double>(std::numeric_limits<int>::max())) {
      Rcpp::stop(label + " must contain finite integer values.");
    }
    out[static_cast<std::size_t>(i)] = static_cast<int>(v);
  }
  return out;
}

std::vector<double> doubles(SEXP x, const R_xlen_t expected,
                            const std::string& label) {
  if (x == R_NilValue || Rf_xlength(x) != expected ||
      (TYPEOF(x) != INTSXP && TYPEOF(x) != REALSXP)) {
    Rcpp::stop(label + " has an incompatible length or type.");
  }
  std::vector<double> out(static_cast<std::size_t>(expected));
  for (R_xlen_t i = 0; i < expected; ++i) {
    const double v = TYPEOF(x) == INTSXP
      ? static_cast<double>(INTEGER(x)[i]) : REAL(x)[i];
    if (!std::isfinite(v)) Rcpp::stop(label + " must be finite.");
    out[static_cast<std::size_t>(i)] = v;
  }
  return out;
}

Tree parse_tree(const Rcpp::List& compiled) {
  if (!has(compiled, "n_tip") || !has(compiled, "n_total") ||
      !has(compiled, "root") || !has(compiled, "parent") ||
      !has(compiled, "child_ptr") || !has(compiled, "children") ||
      !has(compiled, "preorder") || !has(compiled, "postorder") ||
      !has(compiled, "branch_length_by_node")) {
    Rcpp::stop("compiled_tree is missing canonical tree-core fields.");
  }
  const int n_tip = scalar_int(compiled["n_tip"], "compiled_tree$n_tip", 2);
  const int n_total = scalar_int(compiled["n_total"],
                                 "compiled_tree$n_total", n_tip + 1);
  const int root_raw = scalar_int(compiled["root"], "compiled_tree$root", 1);
  const int root = root_raw - 1;
  if (root < n_tip || root >= n_total) {
    Rcpp::stop("compiled_tree$root must be an internal node.");
  }

  const std::vector<int> parent_raw = ints(
    compiled["parent"], n_total, "compiled_tree$parent"
  );
  const std::vector<int> ptr_raw = ints(
    compiled["child_ptr"], n_total + 1, "compiled_tree$child_ptr"
  );
  const std::vector<int> child_raw = ints(
    compiled["children"], n_total - 1, "compiled_tree$children"
  );
  const std::vector<int> preorder_raw = ints(
    compiled["preorder"], n_total, "compiled_tree$preorder"
  );
  const std::vector<int> postorder_raw = ints(
    compiled["postorder"], n_total, "compiled_tree$postorder"
  );
  const std::vector<double> branch = doubles(
    compiled["branch_length_by_node"], n_total,
    "compiled_tree$branch_length_by_node"
  );

  Tree tree;
  tree.n_tip = n_tip;
  tree.n_total = n_total;
  tree.root = root;
  tree.parent.resize(static_cast<std::size_t>(n_total));
  tree.child_ptr.resize(static_cast<std::size_t>(n_total + 1));
  tree.children.resize(static_cast<std::size_t>(n_total - 1));
  tree.preorder.resize(static_cast<std::size_t>(n_total));
  tree.postorder.resize(static_cast<std::size_t>(n_total));
  tree.branch = branch;
  for (int i = 0; i < n_total; ++i) {
    const int par = parent_raw[static_cast<std::size_t>(i)];
    tree.parent[static_cast<std::size_t>(i)] = par == 0 ? -1 : par - 1;
    const int pre = preorder_raw[static_cast<std::size_t>(i)];
    const int post = postorder_raw[static_cast<std::size_t>(i)];
    if (pre < 1 || pre > n_total || post < 1 || post > n_total) {
      Rcpp::stop("compiled_tree traversal arrays contain invalid node ids.");
    }
    tree.preorder[static_cast<std::size_t>(i)] = pre - 1;
    tree.postorder[static_cast<std::size_t>(i)] = post - 1;
  }
  if (ptr_raw.front() != 1 || ptr_raw.back() != n_total) {
    Rcpp::stop("compiled_tree$child_ptr must use the canonical 1-based offsets.");
  }
  for (int i = 0; i <= n_total; ++i) {
    const int begin = ptr_raw[static_cast<std::size_t>(i)] - 1;
    if (begin < 0 || begin > n_total - 1 ||
        (i > 0 && begin < tree.child_ptr[static_cast<std::size_t>(i - 1)])) {
      Rcpp::stop("compiled_tree$child_ptr is not monotone.");
    }
    tree.child_ptr[static_cast<std::size_t>(i)] = begin;
  }
  if (tree.child_ptr.back() != n_total - 1) {
    Rcpp::stop("compiled_tree$child_ptr does not span all edges.");
  }
  for (int i = 0; i < n_total - 1; ++i) {
    const int child = child_raw[static_cast<std::size_t>(i)] - 1;
    if (child < 0 || child >= n_total) {
      Rcpp::stop("compiled_tree$children contains an invalid node id.");
    }
    tree.children[static_cast<std::size_t>(i)] = child;
  }
  for (int node = 0; node < n_total; ++node) {
    if (node != root && tree.branch[static_cast<std::size_t>(node)] <= 0.0) {
      Rcpp::stop("K tree permutation kernel requires positive non-root branches.");
    }
  }

  tree.root_distance.assign(static_cast<std::size_t>(n_total), 0.0);
  for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
    const int node = tree.preorder[ii];
    const int par = tree.parent[static_cast<std::size_t>(node)];
    if (par < 0) Rcpp::stop("compiled_tree has an invalid parent traversal.");
    tree.root_distance[static_cast<std::size_t>(node)] =
      tree.root_distance[static_cast<std::size_t>(par)] +
      tree.branch[static_cast<std::size_t>(node)];
  }
  for (int tip = 0; tip < n_tip; ++tip) {
    if (!std::isfinite(tree.root_distance[static_cast<std::size_t>(tip)])) {
      Rcpp::stop("compiled_tree tip root distances are not finite.");
    }
  }
  return tree;
}

Cache build_cache(const Tree& tree) {
  Cache cache;
  const int n = tree.n_total;
  cache.aggregate.assign(static_cast<std::size_t>(n), 0.0);
  cache.outgoing.assign(static_cast<std::size_t>(n), 0.0);
  for (std::size_t ii = 0; ii < tree.postorder.size(); ++ii) {
    const int node = tree.postorder[ii];
    long double sum = 0.0L;
    if (node < tree.n_tip) {
      sum = 1.0L / static_cast<long double>(
        tree.branch[static_cast<std::size_t>(node)]
      );
    } else {
      const int begin = tree.child_ptr[static_cast<std::size_t>(node)];
      const int end = tree.child_ptr[static_cast<std::size_t>(node + 1)];
      for (int k = begin; k < end; ++k) {
        sum += static_cast<long double>(cache.outgoing[
          static_cast<std::size_t>(tree.children[static_cast<std::size_t>(k)])
        ]);
      }
    }
    if (!(sum > 0.0L) || !std::isfinite(static_cast<double>(sum))) {
      Rcpp::stop("tree Schur-complement precision is not finite and positive.");
    }
    cache.aggregate[static_cast<std::size_t>(node)] = static_cast<double>(sum);
    if (node == tree.root || node < tree.n_tip) {
      cache.outgoing[static_cast<std::size_t>(node)] = static_cast<double>(sum);
    } else {
      const long double l = static_cast<long double>(
        tree.branch[static_cast<std::size_t>(node)]
      );
      const long double out = 1.0L / (l + 1.0L / sum);
      if (!(out > 0.0L) || !std::isfinite(static_cast<double>(out))) {
        Rcpp::stop("tree outgoing precision is not finite and positive.");
      }
      cache.outgoing[static_cast<std::size_t>(node)] = static_cast<double>(out);
    }
  }

  std::vector<double> state(static_cast<std::size_t>(n), 0.0);
  for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
    const int node = tree.preorder[ii];
    const int par = tree.parent[static_cast<std::size_t>(node)];
    const double alpha = node < tree.n_tip ? 1.0 :
      tree.branch[static_cast<std::size_t>(node)] *
      cache.outgoing[static_cast<std::size_t>(node)];
    state[static_cast<std::size_t>(node)] =
      state[static_cast<std::size_t>(par)] +
      alpha * (1.0 - state[static_cast<std::size_t>(par)]);
  }
  long double sum_inv = 0.0L;
  for (int tip = 0; tip < tree.n_tip; ++tip) {
    const int par = tree.parent[static_cast<std::size_t>(tip)];
    sum_inv += (1.0L - static_cast<long double>(
      state[static_cast<std::size_t>(par)]
    )) / static_cast<long double>(tree.branch[static_cast<std::size_t>(tip)]);
  }
  cache.sum_inv = static_cast<double>(sum_inv);
  long double trace = 0.0L;
  for (int tip = 0; tip < tree.n_tip; ++tip) {
    trace += static_cast<long double>(
      tree.root_distance[static_cast<std::size_t>(tip)]
    );
  }
  cache.normalization = static_cast<double>(
    (trace - static_cast<long double>(tree.n_tip) / sum_inv) /
    static_cast<long double>(tree.n_tip - 1)
  );
  if (!(cache.sum_inv > 0.0) || !std::isfinite(cache.sum_inv) ||
      !(cache.normalization > 0.0) || !std::isfinite(cache.normalization)) {
    Rcpp::stop("K tree normalization is not finite and positive.");
  }
  return cache;
}

struct PermView {
  const int* zero_based;
  const int* one_based_matrix;
  int row;
  int nrow;
  bool from_matrix;

  inline int at(const int tip) const {
    return from_matrix
      ? one_based_matrix[row + nrow * tip] - 1
      : zero_based[tip];
  }
};

inline double tip_value(const double* x, const int n_tip,
                        const PermView& perm, const int tip, const int col,
                        Counters* counters) {
  if (counters != NULL) {
    ++counters->index_reads;
    ++counters->trait_reads;
  }
  return x[static_cast<std::size_t>(perm.at(tip)) +
           static_cast<std::size_t>(n_tip) * static_cast<std::size_t>(col)];
}

void compute_one(const Tree& tree, const Cache& cache,
                 const double* x, const int ncol,
                 const PermView& perm, const int trait_chunk,
                 std::vector<double>& out, Counters* counters) {
  const int n = tree.n_tip;
  const int p = ncol;
  out.assign(static_cast<std::size_t>(p),
             std::numeric_limits<double>::quiet_NaN());
  const int chunk_limit = std::max(1, std::min(trait_chunk, p));
  std::vector<double> message;
  std::vector<double> state;
  std::vector<double> baseline;
  std::vector<double> delta;
  for (int col0 = 0; col0 < p; col0 += chunk_limit) {
    const int chunk = std::min(chunk_limit, p - col0);
    const std::size_t cells = static_cast<std::size_t>(tree.n_total) *
      static_cast<std::size_t>(chunk);
    message.assign(cells, 0.0);
    state.assign(cells, 0.0);
    baseline.assign(static_cast<std::size_t>(chunk), 0.0);
    delta.assign(static_cast<std::size_t>(chunk), 0.0);
    if (counters != NULL) {
      counters->compute_workspace_allocations += 4;
      const std::size_t bytes = 2U * cells * sizeof(double) +
        2U * static_cast<std::size_t>(chunk) * sizeof(double);
      counters->compute_workspace_bytes += bytes;
      counters->compute_workspace_peak_bytes = std::max(
        counters->compute_workspace_peak_bytes, bytes
      );
    }

    for (int j = 0; j < chunk; ++j) {
      baseline[static_cast<std::size_t>(j)] =
        tip_value(x, n, perm, 0, col0 + j, counters);
    }
    for (std::size_t ii = 0; ii < tree.postorder.size(); ++ii) {
      const int node = tree.postorder[ii];
      const std::size_t base = static_cast<std::size_t>(node) *
        static_cast<std::size_t>(chunk);
      if (node < n) {
        for (int j = 0; j < chunk; ++j) {
          message[base + static_cast<std::size_t>(j)] =
            tip_value(x, n, perm, node, col0 + j, counters) -
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
    for (int j = 0; j < chunk; ++j) {
      long double qlinear = 0.0L;
      for (int tip = 0; tip < n; ++tip) {
        const int par = tree.parent[static_cast<std::size_t>(tip)];
        const double parent_state = state[
          static_cast<std::size_t>(par) * static_cast<std::size_t>(chunk) +
          static_cast<std::size_t>(j)
        ];
        qlinear += (static_cast<long double>(
          tip_value(x, n, perm, tip, col0 + j, counters) -
          baseline[static_cast<std::size_t>(j)] - parent_state
        )) / static_cast<long double>(tree.branch[static_cast<std::size_t>(tip)]);
      }
      delta[static_cast<std::size_t>(j)] =
        static_cast<double>(qlinear / static_cast<long double>(cache.sum_inv));
    }
    for (std::size_t ii = 0; ii < tree.postorder.size(); ++ii) {
      const int node = tree.postorder[ii];
      const std::size_t base = static_cast<std::size_t>(node) *
        static_cast<std::size_t>(chunk);
      if (node < n) {
        for (int j = 0; j < chunk; ++j) {
          message[base + static_cast<std::size_t>(j)] =
            tip_value(x, n, perm, node, col0 + j, counters) -
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
    for (int j = 0; j < chunk; ++j) {
      long double numerator = 0.0L;
      long double denominator = 0.0L;
      for (int tip = 0; tip < n; ++tip) {
        const long double y = static_cast<long double>(
          tip_value(x, n, perm, tip, col0 + j, counters) -
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
              tip_value(x, n, perm, node, col0 + j, counters) -
              baseline[static_cast<std::size_t>(j)] -
              delta[static_cast<std::size_t>(j)]
            )
          : static_cast<long double>(state[
              static_cast<std::size_t>(node) *
              static_cast<std::size_t>(chunk) + static_cast<std::size_t>(j)
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

void fisher_yates(std::vector<int>& perm, Counters* counters) {
  for (int i = static_cast<int>(perm.size()) - 1; i > 0; --i) {
    const int j = static_cast<int>(std::floor(R::runif(0.0, i + 1.0)));
    std::swap(perm[static_cast<std::size_t>(i)],
              perm[static_cast<std::size_t>(j)]);
    if (counters != NULL) ++counters->fisher_yates_swaps;
  }
}

inline bool inclusive_upper_tail(const double value, const double observed) {
  if (value >= observed) return true;
  if (!std::isfinite(value) || !std::isfinite(observed)) return false;
  const double scale = std::max(1.0, std::max(std::fabs(value),
                                               std::fabs(observed)));
  const double tie = 8.0 * std::numeric_limits<double>::epsilon() * scale;
  return observed - value <= tie;
}

enum class Mode { oracle, candidate_a, candidate_b };

Mode parse_mode(const std::string& value) {
  if (value == "oracle" || value == "production") return Mode::oracle;
  if (value == "A" || value == "a") return Mode::candidate_a;
  if (value == "B" || value == "b") return Mode::candidate_b;
  Rcpp::stop("mode must be one of oracle, A, or B; C is already present in production.");
  return Mode::oracle;
}

int effective_threads(const int requested) {
  if (requested < 1) Rcpp::stop("n_threads must be positive.");
#ifdef _OPENMP
  return std::max(1, std::min(requested, omp_get_max_threads()));
#else
  return 1;
#endif
}

Rcpp::List counters_list(const Counters& c, const std::size_t chunk_bytes,
                         const std::size_t generation_bytes,
                         const std::size_t worker_bytes,
                         const std::size_t output_bytes,
                         const int effective_n_threads) {
  const std::size_t peak = chunk_bytes + generation_bytes +
    static_cast<std::size_t>(effective_n_threads) *
      (worker_bytes + c.compute_workspace_peak_bytes);
  return Rcpp::List::create(
    Rcpp::Named("null_replicates") =
      static_cast<double>(c.null_replicates),
    Rcpp::Named("index_buffer_allocations") =
      static_cast<double>(c.index_buffer_allocations),
    Rcpp::Named("index_initialization_elements") =
      static_cast<double>(c.index_initialization_elements),
    Rcpp::Named("fisher_yates_swaps") =
      static_cast<double>(c.fisher_yates_swaps),
    Rcpp::Named("generation_to_chunk_elements") =
      static_cast<double>(c.generation_to_chunk_elements),
    Rcpp::Named("chunk_to_worker_elements") =
      static_cast<double>(c.chunk_to_worker_elements),
    Rcpp::Named("controlled_to_worker_elements") =
      static_cast<double>(c.controlled_to_worker_elements),
    Rcpp::Named("index_reads") = static_cast<double>(c.index_reads),
    Rcpp::Named("trait_reads") = static_cast<double>(c.trait_reads),
    Rcpp::Named("trait_gather_allocations") =
      static_cast<double>(c.trait_gather_allocations),
    Rcpp::Named("compute_workspace_allocations") =
      static_cast<double>(c.compute_workspace_allocations),
    Rcpp::Named("compute_workspace_bytes") =
      static_cast<double>(c.compute_workspace_bytes),
    Rcpp::Named("compute_workspace_peak_bytes") =
      static_cast<double>(c.compute_workspace_peak_bytes),
    Rcpp::Named("permutation_chunk_bytes") =
      static_cast<double>(chunk_bytes),
    Rcpp::Named("reusable_generation_bytes") =
      static_cast<double>(generation_bytes),
    Rcpp::Named("worker_index_bytes") = static_cast<double>(worker_bytes),
    Rcpp::Named("return_sim_bytes") = static_cast<double>(output_bytes),
    Rcpp::Named("peak_workspace_proxy_bytes") = static_cast<double>(peak),
    Rcpp::Named("effective_n_threads") = effective_n_threads
  );
}

} // namespace stage2d1

// [[Rcpp::export]]
Rcpp::List stage2d1_k_permutation_prototype(
    const Rcpp::List& compiled_tree,
    const Rcpp::NumericMatrix& X,
    const int nsim = 1000,
    SEXP permutations = R_NilValue,
    const int trait_chunk = 64,
    const bool return_sim = false,
    const bool include_observed = true,
    const int n_threads = 1,
    const int simulation_chunk = 128,
    const std::string mode = "oracle",
    const bool return_ordered = false,
    const bool collect_counters = true) {
  using namespace stage2d1;
  if (nsim < 1) Rcpp::stop("nsim must be a positive integer.");
  if (trait_chunk < 1) Rcpp::stop("trait_chunk must be positive.");
  if (simulation_chunk < 1) {
    Rcpp::stop("simulation_chunk must be a positive integer.");
  }
  const Mode candidate = parse_mode(mode);
  const Tree tree = parse_tree(compiled_tree);
  const int n = tree.n_tip;
  const int p = X.ncol();
  if (X.nrow() != n || p < 1) {
    Rcpp::stop("X must have one row per compiled tree tip and at least one trait.");
  }
  const double* x = REAL(X);
  for (R_xlen_t i = 0; i < X.size(); ++i) {
    if (!std::isfinite(x[i])) Rcpp::stop("X must contain only finite values.");
  }
  const Cache cache = build_cache(tree);

  const bool supplied = permutations != R_NilValue;
  Rcpp::IntegerMatrix perm_matrix;
  const int* permutation_data = NULL;
  if (supplied) {
    perm_matrix = Rcpp::IntegerMatrix(permutations);
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
    permutation_data = INTEGER(perm_matrix);
  }

  Counters counters;
  Counters* active_counters = collect_counters ? &counters : NULL;
  Rcpp::NumericVector observed(p);
  std::vector<int> identity(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) identity[static_cast<std::size_t>(i)] = i;
  std::vector<double> observed_vec;
  const PermView identity_view = { identity.data(), NULL, 0, 0, false };
  compute_one(tree, cache, x, p, identity_view, trait_chunk,
              observed_vec, active_counters);
  std::vector<char> valid_observed(static_cast<std::size_t>(p), 1);
  for (int j = 0; j < p; ++j) {
    observed[j] = observed_vec[static_cast<std::size_t>(j)];
    valid_observed[static_cast<std::size_t>(j)] =
      std::isfinite(observed[j]) ? 1 : 0;
  }

  Rcpp::NumericVector exceedance(p, 0.0);
  Rcpp::NumericMatrix sim;
  const bool keep_ordered = return_sim || return_ordered;
  if (keep_ordered) sim = Rcpp::NumericMatrix(nsim, p);
  double* sim_ptr = keep_ordered ? REAL(sim) : NULL;
  const int threads_eff = effective_threads(n_threads);

  if (supplied) {
    if (threads_eff > 1) {
#ifdef _OPENMP
      std::vector<Counters> local_counters(
        static_cast<std::size_t>(threads_eff)
      );
      std::vector<std::vector<double> > local_counts(
        static_cast<std::size_t>(threads_eff),
        std::vector<double>(static_cast<std::size_t>(p), 0.0)
      );
#pragma omp parallel num_threads(threads_eff)
      {
        const int tid = omp_get_thread_num();
        Counters& local = local_counters[static_cast<std::size_t>(tid)];
        std::vector<int> worker_perm;
        if (candidate == Mode::oracle || candidate == Mode::candidate_a) {
          worker_perm.resize(static_cast<std::size_t>(n));
          if (collect_counters) ++local.index_buffer_allocations;
        }
        std::vector<double> kval;
#pragma omp for schedule(static)
        for (int i = 0; i < nsim; ++i) {
          PermView view;
          if (candidate == Mode::oracle || candidate == Mode::candidate_a) {
            for (int r = 0; r < n; ++r) {
              worker_perm[static_cast<std::size_t>(r)] =
                permutation_data[i + nsim * r] - 1;
            }
            if (collect_counters) {
              local.controlled_to_worker_elements +=
                static_cast<std::uint64_t>(n);
            }
            view = { worker_perm.data(), NULL, 0, 0, false };
          } else {
            view = { NULL, permutation_data, i, nsim, true };
          }
          compute_one(tree, cache, x, p, view, trait_chunk, kval,
                      collect_counters ? &local : NULL);
          if (collect_counters) ++local.null_replicates;
          for (int j = 0; j < p; ++j) {
            const double value = kval[static_cast<std::size_t>(j)];
            if (keep_ordered) sim_ptr[static_cast<std::size_t>(i) +
                                      static_cast<std::size_t>(nsim) *
                                      static_cast<std::size_t>(j)] = value;
            if (valid_observed[static_cast<std::size_t>(j)] &&
                inclusive_upper_tail(value, observed[j])) {
              local_counts[static_cast<std::size_t>(tid)]
                [static_cast<std::size_t>(j)] += 1.0;
            }
          }
        }
      }
      for (int t = 0; t < threads_eff; ++t) {
        if (collect_counters) {
          counters.add(local_counters[static_cast<std::size_t>(t)]);
        }
        for (int j = 0; j < p; ++j) {
          exceedance[j] += local_counts[static_cast<std::size_t>(t)]
            [static_cast<std::size_t>(j)];
        }
      }
#else
      (void)threads_eff;
#endif
    } else {
      std::vector<int> worker_perm;
      if (candidate == Mode::oracle || candidate == Mode::candidate_a) {
        worker_perm.resize(static_cast<std::size_t>(n));
        if (collect_counters) ++counters.index_buffer_allocations;
      }
      std::vector<double> kval;
      for (int i = 0; i < nsim; ++i) {
        PermView view;
        if (candidate == Mode::oracle || candidate == Mode::candidate_a) {
          for (int r = 0; r < n; ++r) {
            worker_perm[static_cast<std::size_t>(r)] =
              permutation_data[i + nsim * r] - 1;
          }
          if (collect_counters) {
            counters.controlled_to_worker_elements +=
              static_cast<std::uint64_t>(n);
          }
          view = { worker_perm.data(), NULL, 0, 0, false };
        } else {
          view = { NULL, permutation_data, i, nsim, true };
        }
        compute_one(tree, cache, x, p, view, trait_chunk, kval,
                    active_counters);
        if (collect_counters) ++counters.null_replicates;
        for (int j = 0; j < p; ++j) {
          const double value = kval[static_cast<std::size_t>(j)];
          if (keep_ordered) sim_ptr[static_cast<std::size_t>(i) +
                                    static_cast<std::size_t>(nsim) *
                                    static_cast<std::size_t>(j)] = value;
          if (valid_observed[static_cast<std::size_t>(j)] &&
              inclusive_upper_tail(value, observed[j])) exceedance[j] += 1.0;
        }
      }
    }
  } else {
    const int chunk_limit = std::min(simulation_chunk, nsim);
    std::vector<int> chunk_perms(
      static_cast<std::size_t>(chunk_limit) * static_cast<std::size_t>(n)
    );
    const std::size_t chunk_bytes = chunk_perms.size() * sizeof(int);
    const bool direct_evaluation = candidate == Mode::candidate_b;
    // A and B reuse one generation buffer for the complete public call.
    // The oracle intentionally allocates a fresh vector per replicate so the
    // allocation cost removed by A remains measurable.
    std::vector<int> generation_perm;
    if (candidate != Mode::oracle) {
      generation_perm.resize(static_cast<std::size_t>(n));
      if (collect_counters) ++counters.index_buffer_allocations;
    }
    for (int first = 0; first < nsim; first += chunk_limit) {
      const int count = std::min(chunk_limit, nsim - first);
      for (int local_i = 0; local_i < count; ++local_i) {
        int* row = chunk_perms.data() +
          static_cast<std::size_t>(local_i) * static_cast<std::size_t>(n);
        std::vector<int> per_replicate_perm;
        std::vector<int>* generation = &generation_perm;
        if (candidate == Mode::oracle) {
          per_replicate_perm.resize(static_cast<std::size_t>(n));
          if (collect_counters) ++counters.index_buffer_allocations;
          generation = &per_replicate_perm;
        }
        std::copy(identity.begin(), identity.end(), generation->begin());
        if (collect_counters) {
          counters.index_initialization_elements += static_cast<std::uint64_t>(n);
        }
        if (!(include_observed && first + local_i == 0)) {
          fisher_yates(*generation, active_counters);
        }
        std::copy(generation->begin(), generation->end(), row);
        if (collect_counters) {
          counters.generation_to_chunk_elements += static_cast<std::uint64_t>(n);
        }
      }

      if (threads_eff > 1) {
#ifdef _OPENMP
        std::vector<Counters> local_counters(
          static_cast<std::size_t>(threads_eff)
        );
        std::vector<std::vector<double> > local_counts(
          static_cast<std::size_t>(threads_eff),
          std::vector<double>(static_cast<std::size_t>(p), 0.0)
        );
#pragma omp parallel num_threads(threads_eff)
        {
          const int tid = omp_get_thread_num();
          Counters& local = local_counters[static_cast<std::size_t>(tid)];
          std::vector<int> worker_perm;
          if (!direct_evaluation) {
            worker_perm.resize(static_cast<std::size_t>(n));
            if (collect_counters) ++local.index_buffer_allocations;
          }
          std::vector<double> kval;
#pragma omp for schedule(static)
          for (int local_i = 0; local_i < count; ++local_i) {
            const int* row = chunk_perms.data() +
              static_cast<std::size_t>(local_i) * static_cast<std::size_t>(n);
            PermView view;
            if (direct_evaluation) {
              view = { row, NULL, 0, 0, false };
            } else {
              std::copy(row, row + n, worker_perm.begin());
              if (collect_counters) {
                local.chunk_to_worker_elements += static_cast<std::uint64_t>(n);
              }
              view = { worker_perm.data(), NULL, 0, 0, false };
            }
            compute_one(tree, cache, x, p, view, trait_chunk, kval,
                        collect_counters ? &local : NULL);
            if (collect_counters) ++local.null_replicates;
            const int global_i = first + local_i;
            for (int j = 0; j < p; ++j) {
              const double value = kval[static_cast<std::size_t>(j)];
              if (keep_ordered) sim_ptr[static_cast<std::size_t>(global_i) +
                                        static_cast<std::size_t>(nsim) *
                                        static_cast<std::size_t>(j)] = value;
              if (valid_observed[static_cast<std::size_t>(j)] &&
                  inclusive_upper_tail(value, observed[j])) {
                local_counts[static_cast<std::size_t>(tid)]
                  [static_cast<std::size_t>(j)] += 1.0;
              }
            }
          }
        }
        for (int t = 0; t < threads_eff; ++t) {
          if (collect_counters) {
            counters.add(local_counters[static_cast<std::size_t>(t)]);
          }
          for (int j = 0; j < p; ++j) {
            exceedance[j] += local_counts[static_cast<std::size_t>(t)]
              [static_cast<std::size_t>(j)];
          }
        }
#else
        (void)threads_eff;
#endif
      } else {
        std::vector<int> worker_perm;
        if (!direct_evaluation) {
          worker_perm.resize(static_cast<std::size_t>(n));
          if (collect_counters) ++counters.index_buffer_allocations;
        }
        std::vector<double> kval;
        for (int local_i = 0; local_i < count; ++local_i) {
          const int* row = chunk_perms.data() +
            static_cast<std::size_t>(local_i) * static_cast<std::size_t>(n);
          PermView view;
          if (direct_evaluation) {
            view = { row, NULL, 0, 0, false };
          } else {
            std::copy(row, row + n, worker_perm.begin());
            if (collect_counters) {
              counters.chunk_to_worker_elements += static_cast<std::uint64_t>(n);
            }
            view = { worker_perm.data(), NULL, 0, 0, false };
          }
          compute_one(tree, cache, x, p, view, trait_chunk, kval,
                      active_counters);
          if (collect_counters) ++counters.null_replicates;
          const int global_i = first + local_i;
          for (int j = 0; j < p; ++j) {
            const double value = kval[static_cast<std::size_t>(j)];
            if (keep_ordered) sim_ptr[static_cast<std::size_t>(global_i) +
                                      static_cast<std::size_t>(nsim) *
                                      static_cast<std::size_t>(j)] = value;
            if (valid_observed[static_cast<std::size_t>(j)] &&
                inclusive_upper_tail(value, observed[j])) exceedance[j] += 1.0;
          }
        }
      }
    }

    const std::size_t generation_bytes =
      static_cast<std::size_t>(n) * sizeof(int);
    const std::size_t worker_bytes = direct_evaluation
      ? 0U : static_cast<std::size_t>(n) * sizeof(int);
    const std::size_t output_bytes = keep_ordered
      ? static_cast<std::size_t>(nsim) * static_cast<std::size_t>(p) *
        sizeof(double) : 0U;
    Rcpp::NumericVector p_value(p), mcse(p);
    Rcpp::IntegerVector nsim_successful(p, nsim);
    const int n_random = include_observed ? std::max(0, nsim - 1) : nsim;
    for (int j = 0; j < p; ++j) {
      if (!std::isfinite(observed[j])) {
        p_value[j] = NA_REAL;
        mcse[j] = NA_REAL;
        nsim_successful[j] = 0;
        continue;
      }
      p_value[j] = exceedance[j] / static_cast<double>(nsim);
      if (include_observed && n_random > 0) {
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
    const std::size_t peak_chunk = chunk_bytes + generation_bytes +
      static_cast<std::size_t>(threads_eff) *
        (worker_bytes + counters.compute_workspace_peak_bytes);
    Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("K") = observed,
      Rcpp::Named("P") = p_value,
      Rcpp::Named("exceedance_count") = exceedance,
      Rcpp::Named("nsim_requested") = nsim,
      Rcpp::Named("nsim_successful") = nsim_successful,
      Rcpp::Named("n_randomizations") = n_random,
      Rcpp::Named("MCSE_P") = mcse,
      Rcpp::Named("permutation_mode") = "internal_rng",
      Rcpp::Named("include_observed") = include_observed,
      Rcpp::Named("simulation_chunk") = simulation_chunk,
      Rcpp::Named("mode") = mode,
      Rcpp::Named("sum_inv") = cache.sum_inv,
      Rcpp::Named("normalization") = cache.normalization,
      Rcpp::Named("peak_workspace_proxy_bytes") =
        static_cast<double>(peak_chunk),
      Rcpp::Named("counters") = counters_list(
        counters, chunk_bytes, generation_bytes, worker_bytes, output_bytes,
        threads_eff
      )
    );
    if (return_sim) out["sim_K"] = sim;
    if (return_ordered) out["ordered_null_K"] = sim;
    return out;
  }

  // The controlled path reaches this common result section after processing.
  const std::size_t output_bytes = keep_ordered
    ? static_cast<std::size_t>(nsim) * static_cast<std::size_t>(p) * sizeof(double)
    : 0U;
  Rcpp::NumericVector p_value(p), mcse(p);
  Rcpp::IntegerVector nsim_successful(p, nsim);
  const int n_random = nsim;
  for (int j = 0; j < p; ++j) {
    if (!std::isfinite(observed[j])) {
      p_value[j] = NA_REAL;
      mcse[j] = NA_REAL;
      nsim_successful[j] = 0;
      continue;
    }
    p_value[j] = exceedance[j] / static_cast<double>(nsim);
    if (n_random > 0) {
      const double q = p_value[j];
      mcse[j] = std::sqrt(q * (1.0 - q) / static_cast<double>(n_random));
    } else {
      mcse[j] = NA_REAL;
    }
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("K") = observed,
    Rcpp::Named("P") = p_value,
    Rcpp::Named("exceedance_count") = exceedance,
    Rcpp::Named("nsim_requested") = nsim,
    Rcpp::Named("nsim_successful") = nsim_successful,
    Rcpp::Named("n_randomizations") = n_random,
    Rcpp::Named("MCSE_P") = mcse,
    Rcpp::Named("permutation_mode") = "controlled",
    Rcpp::Named("include_observed") = false,
    Rcpp::Named("simulation_chunk") = simulation_chunk,
    Rcpp::Named("mode") = mode,
    Rcpp::Named("sum_inv") = cache.sum_inv,
    Rcpp::Named("normalization") = cache.normalization,
    Rcpp::Named("peak_workspace_proxy_bytes") =
      static_cast<double>(threads_eff * counters.compute_workspace_peak_bytes),
    Rcpp::Named("counters") = counters_list(
      counters, 0U, 0U,
      (candidate == Mode::oracle || candidate == Mode::candidate_a)
        ? static_cast<std::size_t>(n) * sizeof(int) : 0U,
      output_bytes, threads_eff
    )
  );
  if (return_sim) out["sim_K"] = sim;
  if (return_ordered) out["ordered_null_K"] = sim;
  return out;
}
