// Bounded-memory permutation kernel for Blomberg's K.
//
// This module deliberately consumes the compiled-tree representation used by
// fast_k_tree_batch_cpp().  The Schur-complement messages below are the same
// Brownian-tree precision calculation as the production K tree engine; only
// the permutation loop is new.  In particular, the implementation never
// constructs an n_tip x nsim state matrix.  A trait chunk and one permutation
// are processed at a time, while an optional sim_K matrix is retained only
// when the caller explicitly requests it.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#ifdef _OPENMP
# include <omp.h>
#endif

// [[Rcpp::plugins(cpp11)]]

namespace kperm {

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

// Parse the canonical list returned by compile_tree_cpp().  The list is
// intentionally parsed into compact zero-based arrays once per call; this is
// O(n_tip) and does not depend on the number of requested permutations.
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

  // 1'Q1 can be obtained by the same downward interpolation used by the
  // ordinary trait pass, applied to an all-ones vector.
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

inline double tip_value(const double* x, const int n_tip,
                        const std::vector<int>& perm,
                        const int tip, const int col) {
  return x[static_cast<std::size_t>(perm[static_cast<std::size_t>(tip)]) +
           static_cast<std::size_t>(n_tip) * static_cast<std::size_t>(col)];
}

void compute_one(const Tree& tree, const Cache& cache,
                 const double* x, const int ncol,
                 const std::vector<int>& perm, const int trait_chunk,
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
        tip_value(x, n, perm, 0, col0 + j);
    }
    // Upward Gaussian messages for x - baseline.
    for (std::size_t ii = 0; ii < tree.postorder.size(); ++ii) {
      const int node = tree.postorder[ii];
      const std::size_t base = static_cast<std::size_t>(node) *
        static_cast<std::size_t>(chunk);
      if (node < n) {
        for (int j = 0; j < chunk; ++j) {
          message[base + static_cast<std::size_t>(j)] =
            tip_value(x, n, perm, node, col0 + j) -
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
          tip_value(x, n, perm, tip, col0 + j) -
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
            tip_value(x, n, perm, node, col0 + j) -
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
          tip_value(x, n, perm, tip, col0 + j) -
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
              tip_value(x, n, perm, node, col0 + j) -
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

void fisher_yates(std::vector<int>& perm) {
  for (int i = static_cast<int>(perm.size()) - 1; i > 0; --i) {
    const int j = static_cast<int>(std::floor(R::runif(0.0, i + 1.0)));
    std::swap(perm[static_cast<std::size_t>(i)],
              perm[static_cast<std::size_t>(j)]);
  }
}

int effective_threads(const int requested) {
  if (requested < 1) Rcpp::stop("n_threads must be positive.");
#ifdef _OPENMP
  return std::max(1, std::min(requested, omp_get_max_threads()));
#else
  return 1;
#endif
}

} // namespace kperm

// [[Rcpp::export]]
Rcpp::List fast_k_tree_permutation_cpp(
    const Rcpp::List& compiled_tree,
    const Rcpp::NumericMatrix& X,
    const int nsim = 1000,
    SEXP permutations = R_NilValue,
    const int trait_chunk = 64,
    const bool return_sim = false,
    const bool include_observed = true,
    const int n_threads = 1,
    const int simulation_chunk = 128) {
  using namespace kperm;
  if (nsim < 1) Rcpp::stop("nsim must be a positive integer.");
  if (trait_chunk < 1) Rcpp::stop("trait_chunk must be positive.");
  if (simulation_chunk < 1) {
    Rcpp::stop("simulation_chunk must be a positive integer.");
  }
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
  }

  Rcpp::NumericVector observed(p);
  std::vector<int> identity(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) identity[static_cast<std::size_t>(i)] = i;
  std::vector<double> observed_vec;
  compute_one(tree, cache, x, p, identity, trait_chunk, observed_vec);
  std::vector<char> valid_observed(static_cast<std::size_t>(p), 1);
  for (int j = 0; j < p; ++j) {
    observed[j] = observed_vec[static_cast<std::size_t>(j)];
    // phytools::phylosig() cannot form a meaningful K randomization test
    // when the observed K is non-finite (for example, a constant trait).
    // Mark that trait invalid rather than reporting an artificial P = 0.
    valid_observed[static_cast<std::size_t>(j)] =
      std::isfinite(observed[j]) ? 1 : 0;
  }

  Rcpp::NumericVector exceedance(p, 0.0);
  Rcpp::NumericMatrix sim;
  if (return_sim) sim = Rcpp::NumericMatrix(nsim, p);
  // Use the raw column-major payload inside optional OpenMP loops; this keeps
  // the threaded path free of Rcpp proxy operations.
  double* sim_ptr = return_sim ? REAL(sim) : NULL;
  const int threads_eff = effective_threads(n_threads);

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
        compute_one(tree, cache, x, p, perm, trait_chunk, kval);
        for (int j = 0; j < p; ++j) {
          const double value = kval[static_cast<std::size_t>(j)];
          if (return_sim) sim_ptr[static_cast<std::size_t>(i) +
                                  static_cast<std::size_t>(nsim) *
                                  static_cast<std::size_t>(j)] = value;
          if (valid_observed[static_cast<std::size_t>(j)] &&
              value >= observed[j]) local[static_cast<std::size_t>(tid)]
            [static_cast<std::size_t>(j)] += 1.0;
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
    // Controlled matrices are already materialized by the caller.  Keep the
    // existing serial fallback when OpenMP is unavailable or n_threads=1.
    std::vector<int> perm(static_cast<std::size_t>(n));
    std::vector<double> kval;
    for (int i = 0; i < nsim; ++i) {
      for (int r = 0; r < n; ++r) perm[static_cast<std::size_t>(r)] =
        perm_matrix(i, r) - 1;
      compute_one(tree, cache, x, p, perm, trait_chunk, kval);
      for (int j = 0; j < p; ++j) {
        const double value = kval[static_cast<std::size_t>(j)];
        if (return_sim) sim_ptr[static_cast<std::size_t>(i) +
                                static_cast<std::size_t>(nsim) *
                                static_cast<std::size_t>(j)] = value;
        if (valid_observed[static_cast<std::size_t>(j)] &&
            value >= observed[j]) exceedance[j] += 1.0;
      }
    }
  } else {
    // Internal RNG mode is genuinely chunked: only simulation_chunk
    // permutations are materialized at any time.  Fisher-Yates draws are
    // generated serially under R's RNG, so the stream is independent of both
    // OpenMP thread count and the chosen chunk size.
    const int chunk_limit = std::min(simulation_chunk, nsim);
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
          fisher_yates(perm);
        }
        std::copy(perm.begin(), perm.end(),
                  chunk_perms.begin() +
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
            compute_one(tree, cache, x, p, perm, trait_chunk, kval);
            const int global_i = first + local_i;
            for (int j = 0; j < p; ++j) {
              const double value = kval[static_cast<std::size_t>(j)];
              if (return_sim) sim_ptr[static_cast<std::size_t>(global_i) +
                                      static_cast<std::size_t>(nsim) *
                                      static_cast<std::size_t>(j)] = value;
              if (valid_observed[static_cast<std::size_t>(j)] &&
                  value >= observed[j]) local_counts[
                    static_cast<std::size_t>(tid)
                  ][static_cast<std::size_t>(j)] += 1.0;
            }
          }
        }
        for (int t = 0; t < threads_eff; ++t) {
          for (int j = 0; j < p; ++j) exceedance[j] += local_counts[
            static_cast<std::size_t>(t)
          ][static_cast<std::size_t>(j)];
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
          compute_one(tree, cache, x, p, perm, trait_chunk, kval);
          const int global_i = first + local_i;
          for (int j = 0; j < p; ++j) {
            const double value = kval[static_cast<std::size_t>(j)];
            if (return_sim) sim_ptr[static_cast<std::size_t>(global_i) +
                                    static_cast<std::size_t>(nsim) *
                                    static_cast<std::size_t>(j)] = value;
            if (valid_observed[static_cast<std::size_t>(j)] &&
                value >= observed[j]) exceedance[j] += 1.0;
          }
        }
      }
    }
  }

  Rcpp::NumericVector p_value(p), mcse(p);
  Rcpp::IntegerVector nsim_successful(p, nsim);
  const int n_random = supplied ? nsim :
    (include_observed ? std::max(0, nsim - 1) : nsim);
  for (int j = 0; j < p; ++j) {
    if (!valid_observed[static_cast<std::size_t>(j)]) {
      p_value[j] = NA_REAL;
      mcse[j] = NA_REAL;
      nsim_successful[j] = 0;
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

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("K") = observed,
    Rcpp::Named("P") = p_value,
    Rcpp::Named("exceedance_count") = exceedance,
    Rcpp::Named("nsim_requested") = nsim,
    Rcpp::Named("nsim_successful") = nsim_successful,
    Rcpp::Named("n_randomizations") = n_random,
    Rcpp::Named("MCSE_P") = mcse,
    Rcpp::Named("permutation_mode") = supplied ? "controlled" : "internal_rng",
    Rcpp::Named("include_observed") = (!supplied && include_observed),
    Rcpp::Named("simulation_chunk") = simulation_chunk,
    Rcpp::Named("sum_inv") = cache.sum_inv,
    Rcpp::Named("normalization") = cache.normalization
  );
  if (return_sim) out["sim_K"] = sim;
  return out;
}
