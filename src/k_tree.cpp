// Exact Brownian-tree Schur-complement kernel for Blomberg's K.
//
// The covariance of the tip values under ape::vcv.phylo() is the covariance
// obtained by putting a fixed (zero) state at the root and assigning an
// independent Brownian increment of variance l_e to every edge e.  Eliminating
// the internal states therefore gives the precision quadratic x' Q x by a
// Gaussian message pass on the tree.  The implementation below keeps the
// edge precisions (which depend only on the tree) once, and processes trait
// columns in bounded chunks.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

namespace {

struct TreeArrays {
  int n_tip;
  int n_total;
  int root; // zero based
  std::vector<int> parent;
  std::vector<int> child_ptr; // zero based CSR offsets, length n_total + 1
  std::vector<int> children;  // zero based node ids
  std::vector<int> preorder;
  std::vector<int> postorder;
  std::vector<double> branch; // branch[node] is edge parent -> node
  std::vector<double> root_distance;
};

bool has_field(const Rcpp::List& x, const char* name) {
  return x.containsElementNamed(name) && x[name] != R_NilValue;
}

double scalar_number(SEXP value, const std::string& label) {
  if (value == R_NilValue || Rf_xlength(value) != 1 ||
      (TYPEOF(value) != INTSXP && TYPEOF(value) != REALSXP)) {
    Rcpp::stop(label + " must be a single numeric value.");
  }
  const double out = TYPEOF(value) == INTSXP
    ? static_cast<double>(INTEGER(value)[0])
    : REAL(value)[0];
  if (!std::isfinite(out)) {
    Rcpp::stop(label + " must be finite.");
  }
  return out;
}

int scalar_integer(SEXP value, const std::string& label, const bool positive) {
  const double raw = scalar_number(value, label);
  if (raw < static_cast<double>(std::numeric_limits<int>::min()) ||
      raw > static_cast<double>(std::numeric_limits<int>::max()) ||
      std::floor(raw) != raw || (positive && raw < 1.0)) {
    Rcpp::stop(label + " must be an integer" +
               std::string(positive ? " >= 1." : "."));
  }
  return static_cast<int>(raw);
}

SEXP first_field(const Rcpp::List& x,
                 const char* primary,
                 const char* secondary = NULL,
                 const char* tertiary = NULL) {
  if (has_field(x, primary)) return x[primary];
  if (secondary != NULL && has_field(x, secondary)) return x[secondary];
  if (tertiary != NULL && has_field(x, tertiary)) return x[tertiary];
  return R_NilValue;
}

std::vector<double> numeric_vector(SEXP value,
                                   const R_xlen_t expected,
                                   const std::string& label) {
  if (value == R_NilValue ||
      (TYPEOF(value) != INTSXP && TYPEOF(value) != REALSXP) ||
      Rf_xlength(value) != expected) {
    Rcpp::stop(label + " has an incompatible length or type.");
  }
  std::vector<double> out(static_cast<std::size_t>(expected));
  for (R_xlen_t i = 0; i < expected; ++i) {
    const double v = TYPEOF(value) == INTSXP
      ? static_cast<double>(INTEGER(value)[i])
      : REAL(value)[i];
    if (!std::isfinite(v)) {
      Rcpp::stop(label + " must contain only finite values.");
    }
    out[static_cast<std::size_t>(i)] = v;
  }
  return out;
}

std::vector<int> integer_vector(SEXP value,
                                const R_xlen_t expected,
                                const std::string& label) {
  if (value == R_NilValue ||
      (TYPEOF(value) != INTSXP && TYPEOF(value) != REALSXP) ||
      Rf_xlength(value) != expected) {
    Rcpp::stop(label + " has an incompatible length or type.");
  }
  std::vector<int> out(static_cast<std::size_t>(expected));
  for (R_xlen_t i = 0; i < expected; ++i) {
    const double v = TYPEOF(value) == INTSXP
      ? static_cast<double>(INTEGER(value)[i])
      : REAL(value)[i];
    if (!std::isfinite(v) || std::floor(v) != v ||
        v < static_cast<double>(std::numeric_limits<int>::min()) ||
        v > static_cast<double>(std::numeric_limits<int>::max())) {
      Rcpp::stop(label + " must contain integer values.");
    }
    out[static_cast<std::size_t>(i)] = static_cast<int>(v);
  }
  return out;
}

// Read a compiled tree.  The canonical representation is the list produced
// by compile_tree_cpp(), but aliases used by older prepared-tree objects are
// accepted as well.  Traversal orders are rebuilt from CSR so that an
// explicit list need not carry postorder/preorder fields.
TreeArrays parse_tree(const Rcpp::List& compiled) {
  if (!has_field(compiled, "n_tip") || !has_field(compiled, "n_total")) {
    Rcpp::stop("compiled_tree must contain n_tip and n_total.");
  }
  const int n_tip = scalar_integer(compiled["n_tip"], "compiled_tree$n_tip", true);
  const int n_total = scalar_integer(
    compiled["n_total"], "compiled_tree$n_total", true
  );
  if (n_tip < 2 || n_total <= n_tip) {
    Rcpp::stop("compiled_tree must contain at least two tips and one internal node.");
  }

  SEXP root_sexp = first_field(compiled, "root");
  int root = 0;
  if (root_sexp != R_NilValue) {
    const int root_raw = scalar_integer(root_sexp, "compiled_tree$root", false);
    // Canonical compiled trees use one-based node ids; an explicit CSR list
    // may use zero-based ids (root == 0 is unambiguous).
    root = root_raw == 0 ? 0 : root_raw - 1;
  }

  SEXP ptr_sexp = first_field(
    compiled, "child_ptr", "children_ptr", "children_offset"
  );
  SEXP child_sexp = first_field(compiled, "children");
  if (ptr_sexp == R_NilValue || child_sexp == R_NilValue) {
    Rcpp::stop("compiled_tree must contain CSR child_ptr and children arrays.");
  }
  std::vector<int> raw_ptr = integer_vector(
    ptr_sexp, static_cast<R_xlen_t>(n_total + 1),
    "compiled_tree child_ptr"
  );
  std::vector<int> raw_children = integer_vector(
    child_sexp, static_cast<R_xlen_t>(n_total - 1),
    "compiled_tree children"
  );

  // compile_tree_cpp uses one-based offsets and node ids.  Supporting a
  // zero-based explicit array too costs little and makes the low-level entry
  // point useful in tests and serialized contexts.
  const bool ptr_zero_based = raw_ptr.front() == 0 &&
    raw_ptr.back() == n_total - 1;
  const bool ptr_one_based = raw_ptr.front() == 1 &&
    raw_ptr.back() == n_total;
  if (!ptr_zero_based && !ptr_one_based) {
    Rcpp::stop("compiled_tree child_ptr must be a contiguous CSR offset array.");
  }
  std::vector<int> child_ptr(static_cast<std::size_t>(n_total + 1));
  for (int i = 0; i <= n_total; ++i) {
    child_ptr[static_cast<std::size_t>(i)] =
      ptr_zero_based ? raw_ptr[static_cast<std::size_t>(i)]
                     : raw_ptr[static_cast<std::size_t>(i)] - 1;
    if (child_ptr[static_cast<std::size_t>(i)] < 0 ||
        child_ptr[static_cast<std::size_t>(i)] > n_total - 1 ||
        (i > 0 && child_ptr[static_cast<std::size_t>(i)] <
                    child_ptr[static_cast<std::size_t>(i - 1)])) {
      Rcpp::stop("compiled_tree child_ptr is not monotone and in range.");
    }
  }
  if (child_ptr.back() != n_total - 1) {
    Rcpp::stop("compiled_tree child_ptr does not span all edges.");
  }

  bool children_zero_based = false;
  for (std::size_t i = 0; i < raw_children.size(); ++i) {
    if (raw_children[i] == 0) {
      children_zero_based = true;
      break;
    }
  }
  std::vector<int> children(raw_children.size());
  for (std::size_t i = 0; i < raw_children.size(); ++i) {
    const int child = children_zero_based
      ? raw_children[i]
      : raw_children[i] - 1;
    if (child < 0 || child >= n_total) {
      Rcpp::stop("compiled_tree children contain an out-of-range node.");
    }
    children[i] = child;
  }

  SEXP branch_sexp = first_field(
    compiled, "branch_length_by_node", "branch_length"
  );
  if (branch_sexp == R_NilValue) {
    Rcpp::stop("compiled_tree must contain branch_length_by_node.");
  }
  std::vector<double> branch = numeric_vector(
    branch_sexp, static_cast<R_xlen_t>(n_total),
    "compiled_tree branch_length_by_node"
  );

  if (root_sexp == R_NilValue) {
    // Infer the unique root from the child incidence structure.
    std::vector<int> indegree(static_cast<std::size_t>(n_total), 0);
    for (std::size_t i = 0; i < children.size(); ++i) {
      ++indegree[static_cast<std::size_t>(children[i])];
    }
    int candidate = -1;
    for (int node = 0; node < n_total; ++node) {
      if (indegree[static_cast<std::size_t>(node)] == 0) {
        if (candidate >= 0) {
          Rcpp::stop("compiled_tree must have exactly one root.");
        }
        candidate = node;
      }
    }
    if (candidate < 0) Rcpp::stop("compiled_tree has no root.");
    root = candidate;
  }
  if (root < n_tip || root >= n_total) {
    Rcpp::stop("compiled_tree root must be an internal node.");
  }

  std::vector<int> parent(static_cast<std::size_t>(n_total), -1);
  for (int node = 0; node < n_total; ++node) {
    const int begin = child_ptr[static_cast<std::size_t>(node)];
    const int end = child_ptr[static_cast<std::size_t>(node + 1)];
    for (int k = begin; k < end; ++k) {
      const int child = children[static_cast<std::size_t>(k)];
      int& par = parent[static_cast<std::size_t>(child)];
      if (par >= 0) {
        Rcpp::stop("compiled_tree has a child with more than one parent.");
      }
      par = node;
    }
  }
  if (parent[static_cast<std::size_t>(root)] >= 0) {
    Rcpp::stop("compiled_tree root must have no parent.");
  }
  for (int node = 0; node < n_total; ++node) {
    const int degree = child_ptr[static_cast<std::size_t>(node + 1)] -
      child_ptr[static_cast<std::size_t>(node)];
    if (node < n_tip && degree != 0) {
      Rcpp::stop("tip nodes must not have children in compiled_tree.");
    }
    if (node >= n_tip && degree == 0) {
      Rcpp::stop("internal nodes must have at least one child.");
    }
    if (node != root && parent[static_cast<std::size_t>(node)] < 0) {
      Rcpp::stop("compiled_tree is disconnected from its root.");
    }
  }

  for (int node = 0; node < n_total; ++node) {
    if (!std::isfinite(branch[static_cast<std::size_t>(node)])) {
      Rcpp::stop("compiled_tree branch lengths must be finite.");
    }
    if (node != root && branch[static_cast<std::size_t>(node)] <= 0.0) {
      Rcpp::stop(
        "fast_k_tree_batch_cpp currently requires strictly positive branch lengths."
      );
    }
  }

  // Iterative DFS gives deterministic traversal and detects cycles without
  // recursion depth depending on the number of tips.
  struct Frame { int node; int next; };
  std::vector<char> seen(static_cast<std::size_t>(n_total), 0);
  std::vector<Frame> stack;
  std::vector<int> preorder;
  std::vector<int> postorder;
  stack.push_back(Frame{root, child_ptr[static_cast<std::size_t>(root)]});
  seen[static_cast<std::size_t>(root)] = 1;
  preorder.push_back(root);
  while (!stack.empty()) {
    Frame& frame = stack.back();
    const int end = child_ptr[static_cast<std::size_t>(frame.node + 1)];
    if (frame.next < end) {
      const int child = children[static_cast<std::size_t>(frame.next++)];
      if (seen[static_cast<std::size_t>(child)]) {
        Rcpp::stop("compiled_tree contains a cycle or repeated node.");
      }
      seen[static_cast<std::size_t>(child)] = 1;
      preorder.push_back(child);
      stack.push_back(Frame{child, child_ptr[static_cast<std::size_t>(child)]});
    } else {
      postorder.push_back(frame.node);
      stack.pop_back();
    }
  }
  for (int node = 0; node < n_total; ++node) {
    if (!seen[static_cast<std::size_t>(node)]) {
      Rcpp::stop("compiled_tree is disconnected from its root.");
    }
  }

  std::vector<double> root_distance(static_cast<std::size_t>(n_total), 0.0);
  for (std::size_t i = 1; i < preorder.size(); ++i) {
    const int node = preorder[i];
    const int par = parent[static_cast<std::size_t>(node)];
    const long double distance =
      static_cast<long double>(root_distance[static_cast<std::size_t>(par)]) +
      static_cast<long double>(branch[static_cast<std::size_t>(node)]);
    if (!std::isfinite(static_cast<double>(distance))) {
      Rcpp::stop("compiled_tree root distances are not finite.");
    }
    root_distance[static_cast<std::size_t>(node)] =
      static_cast<double>(distance);
  }

  TreeArrays out;
  out.n_tip = n_tip;
  out.n_total = n_total;
  out.root = root;
  out.parent.swap(parent);
  out.child_ptr.swap(child_ptr);
  out.children.swap(children);
  out.preorder.swap(preorder);
  out.postorder.swap(postorder);
  out.branch.swap(branch);
  out.root_distance.swap(root_distance);
  return out;
}

// Compute one scalar edge precision for every node.  For an internal node v,
// s_v is the parallel sum of its child precisions and w_v is the result after
// adding the branch above v.  This is exactly the scalar Schur complement of
// the latent state at v.
void compute_edge_weights(const TreeArrays& tree,
                          std::vector<double>& aggregate,
                          std::vector<double>& outgoing) {
  const int n = tree.n_total;
  aggregate.assign(static_cast<std::size_t>(n), 0.0);
  outgoing.assign(static_cast<std::size_t>(n), 0.0);
  for (std::size_t ii = 0; ii < tree.postorder.size(); ++ii) {
    const int node = tree.postorder[ii];
    long double s = 0.0L;
    if (node < tree.n_tip) {
      const double l = tree.branch[static_cast<std::size_t>(node)];
      s = 1.0L / static_cast<long double>(l);
    } else {
      const int begin = tree.child_ptr[static_cast<std::size_t>(node)];
      const int end = tree.child_ptr[static_cast<std::size_t>(node + 1)];
      for (int k = begin; k < end; ++k) {
        s += static_cast<long double>(outgoing[
          static_cast<std::size_t>(tree.children[static_cast<std::size_t>(k)])
        ]);
      }
    }
    if (!(s > 0.0L) || !std::isfinite(static_cast<double>(s))) {
      Rcpp::stop("tree Schur-complement precision is not finite and positive.");
    }
    aggregate[static_cast<std::size_t>(node)] = static_cast<double>(s);
    if (node == tree.root) {
      outgoing[static_cast<std::size_t>(node)] = static_cast<double>(s);
    } else if (node < tree.n_tip) {
      outgoing[static_cast<std::size_t>(node)] = static_cast<double>(s);
    } else {
      // 1/(l + 1/s) avoids forming l*s, which can overflow even though the
      // resulting Schur complement is perfectly finite.
      const long double l = static_cast<long double>(
        tree.branch[static_cast<std::size_t>(node)]
      );
      const long double w = 1.0L / (l + 1.0L / s);
      if (!(w > 0.0L) || !std::isfinite(static_cast<double>(w))) {
        Rcpp::stop("tree Schur-complement precision is not finite and positive.");
      }
      outgoing[static_cast<std::size_t>(node)] = static_cast<double>(w);
    }
  }
}

double edge_alpha(const TreeArrays& tree,
                  const int node,
                  const double outgoing) {
  if (node < tree.n_tip) return 1.0;
  const long double alpha =
    static_cast<long double>(tree.branch[static_cast<std::size_t>(node)]) *
    static_cast<long double>(outgoing);
  if (!std::isfinite(static_cast<double>(alpha))) {
    Rcpp::stop("tree state interpolation is not finite.");
  }
  // Roundoff can put a mathematically [0,1] coefficient a few ulps outside
  // its interval.  Clamping only that roundoff keeps tip states exact without
  // changing valid Schur-complement values.
  if (alpha <= 0.0L) return 0.0;
  if (alpha >= 1.0L) return 1.0;
  return static_cast<double>(alpha);
}

// Upward messages for one trait chunk.  Traits are represented as
// (X - baseline) - delta throughout the two passes.  Keeping the baseline and
// centered GLS delta separate is important for large common offsets: forming
// X - (baseline + delta) in one operation can discard the low bits of the
// residuals before the Schur complement sees them.
double tree_branch_precision(const TreeArrays& tree,
                             const int node,
                             const std::vector<double>& outgoing);

void upward_messages(const TreeArrays& tree,
                     const Rcpp::NumericMatrix& X,
                     const int col0,
                     const int chunk,
                     const std::vector<double>& outgoing,
                     const std::vector<double>& baseline,
                     const std::vector<double>& delta,
                     std::vector<double>& message) {
  const int n = tree.n_total;
  for (std::size_t ii = 0; ii < tree.postorder.size(); ++ii) {
    const int node = tree.postorder[ii];
    const std::size_t base = static_cast<std::size_t>(node) *
      static_cast<std::size_t>(chunk);
    if (node < tree.n_tip) {
      for (int j = 0; j < chunk; ++j) {
        const double value =
          (X(node, col0 + j) - baseline[static_cast<std::size_t>(j)]) -
          delta[static_cast<std::size_t>(j)];
        if (!std::isfinite(value)) {
          Rcpp::stop("X must contain only finite trait values.");
        }
        message[base + static_cast<std::size_t>(j)] = value;
      }
      continue;
    }

    const int begin = tree.child_ptr[static_cast<std::size_t>(node)];
    const int end = tree.child_ptr[static_cast<std::size_t>(node + 1)];
    const double s = tree_branch_precision(tree, node, outgoing);
    if (!(s > 0.0) || !std::isfinite(s)) {
      Rcpp::stop("tree Schur-complement precision is not finite and positive.");
    }
    for (int j = 0; j < chunk; ++j) {
      long double weighted = 0.0L;
      for (int k = begin; k < end; ++k) {
        const int child = tree.children[static_cast<std::size_t>(k)];
        weighted += static_cast<long double>(outgoing[
          static_cast<std::size_t>(child)
        ]) * static_cast<long double>(message[
          static_cast<std::size_t>(child) * static_cast<std::size_t>(chunk) +
          static_cast<std::size_t>(j)
        ]);
      }
      const double mean = static_cast<double>(weighted / static_cast<long double>(s));
      if (!std::isfinite(mean)) {
        Rcpp::stop("tree Gaussian messages are not finite.");
      }
      message[base + static_cast<std::size_t>(j)] = mean;
    }
  }
  (void)n;
}

// The aggregate precision used by upward_messages.  For tips outgoing is
// already the terminal precision; for an internal node the aggregate is the
// parallel sum of child outgoing precisions.  Recompute it from the CSR list
// here to keep the helper independent of a second mutable array.
double tree_branch_precision(const TreeArrays& tree,
                             const int node,
                             const std::vector<double>& outgoing) {
  if (node < tree.n_tip) return outgoing[static_cast<std::size_t>(node)];
  long double s = 0.0L;
  const int begin = tree.child_ptr[static_cast<std::size_t>(node)];
  const int end = tree.child_ptr[static_cast<std::size_t>(node + 1)];
  for (int k = begin; k < end; ++k) {
    s += static_cast<long double>(outgoing[
      static_cast<std::size_t>(tree.children[static_cast<std::size_t>(k)])
    ]);
  }
  return static_cast<double>(s);
}

void downward_states(const TreeArrays& tree,
                     const int col0,
                     const int chunk,
                     const std::vector<double>& outgoing,
                     const std::vector<double>& message,
                     std::vector<double>& state) {
  const std::size_t stride = static_cast<std::size_t>(chunk);
  const std::size_t root_base = static_cast<std::size_t>(tree.root) * stride;
  for (int j = 0; j < chunk; ++j) {
    state[root_base + static_cast<std::size_t>(j)] = 0.0;
  }
  for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
    const int node = tree.preorder[ii];
    const int parent = tree.parent[static_cast<std::size_t>(node)];
    const double alpha = edge_alpha(
      tree, node, outgoing[static_cast<std::size_t>(node)]
    );
    const std::size_t base = static_cast<std::size_t>(node) * stride;
    const std::size_t pbase = static_cast<std::size_t>(parent) * stride;
    for (int j = 0; j < chunk; ++j) {
      const double zp = state[pbase + static_cast<std::size_t>(j)];
      const double m = message[base + static_cast<std::size_t>(j)];
      const double value = zp + alpha * (m - zp);
      if (!std::isfinite(value)) {
        Rcpp::stop("tree Gaussian states are not finite.");
      }
      state[base + static_cast<std::size_t>(j)] = value;
    }
  }
  (void)col0;
}

double sum_inv_one(const TreeArrays& tree,
                   const std::vector<double>& outgoing) {
  const std::size_t n = static_cast<std::size_t>(tree.n_total);
  std::vector<double> state(n, 0.0);
  std::vector<double> message(n, 1.0);
  std::vector<double> one_offset(1, 0.0);
  // For an all-ones trait, every upward mean is exactly one.  Passing a
  // one-filled message vector through the ordinary downward interpolation is
  // sufficient and avoids forming an R matrix just for Q1.
  (void)one_offset;
  const std::size_t root_base = static_cast<std::size_t>(tree.root);
  state[root_base] = 0.0;
  for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
    const int node = tree.preorder[ii];
    const int parent = tree.parent[static_cast<std::size_t>(node)];
    const double alpha = edge_alpha(
      tree, node, outgoing[static_cast<std::size_t>(node)]
    );
    state[static_cast<std::size_t>(node)] =
      state[static_cast<std::size_t>(parent)] +
      alpha * (1.0 - state[static_cast<std::size_t>(parent)]);
  }
  long double sum = 0.0L;
  for (int tip = 0; tip < tree.n_tip; ++tip) {
    const int parent = tree.parent[static_cast<std::size_t>(tip)];
    const long double current =
      (1.0L - static_cast<long double>(state[static_cast<std::size_t>(parent)])) /
      static_cast<long double>(tree.branch[static_cast<std::size_t>(tip)]);
    sum += current;
  }
  const double out = static_cast<double>(sum);
  if (!std::isfinite(out) || out <= 0.0) {
    Rcpp::stop("1'Q1 is not finite and positive for compiled_tree.");
  }
  return out;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List fast_k_tree_batch_cpp(const Rcpp::List& compiled_tree,
                                 const Rcpp::NumericMatrix& X,
                                 const int trait_chunk = 64) {
  const TreeArrays tree = parse_tree(compiled_tree);
  const R_xlen_t nrow = X.nrow();
  const R_xlen_t ncol = X.ncol();
  if (nrow != static_cast<R_xlen_t>(tree.n_tip)) {
    Rcpp::stop("X must have one row per compiled_tree tip (in tip-id order).");
  }
  if (ncol < 1) {
    Rcpp::stop("X must contain at least one trait column.");
  }
  if (ncol > static_cast<R_xlen_t>(std::numeric_limits<int>::max())) {
    Rcpp::stop("X has too many trait columns for this kernel.");
  }
  if (trait_chunk < 1) {
    Rcpp::stop("trait_chunk must be a positive integer.");
  }

  const int p = static_cast<int>(ncol);
  std::vector<double> aggregate;
  std::vector<double> outgoing;
  compute_edge_weights(tree, aggregate, outgoing);
  const double sum_inv = sum_inv_one(tree, outgoing);

  long double trace_ld = 0.0L;
  for (int tip = 0; tip < tree.n_tip; ++tip) {
    trace_ld += static_cast<long double>(
      tree.root_distance[static_cast<std::size_t>(tip)]
    );
  }
  const double trace = static_cast<double>(trace_ld);
  const double norm_const =
    (trace - static_cast<double>(tree.n_tip) / sum_inv) /
    static_cast<double>(tree.n_tip - 1);
  if (!std::isfinite(norm_const) || norm_const <= 0.0) {
    Rcpp::stop("the Brownian-tree K normalization is not finite and positive.");
  }

  Rcpp::NumericVector means(p);
  Rcpp::NumericVector numerator(p);
  Rcpp::NumericVector denominator(p);
  Rcpp::NumericVector kval(p);

  const int chunk_limit = std::max(1, std::min(trait_chunk, p));
  for (int col0 = 0; col0 < p; col0 += chunk_limit) {
    const int chunk = std::min(chunk_limit, p - col0);
    const std::size_t cells = static_cast<std::size_t>(tree.n_total) *
      static_cast<std::size_t>(chunk);
    std::vector<double> message(cells, 0.0);
    std::vector<double> state(cells, 0.0);
    // Translate each trait by one observed tip before the first solve.  The
    // Brownian precision is not translation invariant (the root is fixed),
    // but q(x) = q(x-c1) + c q(1) is exact.  This keeps the state and current
    // arithmetic well-scaled for traits carrying a large common offset.
    std::vector<double> baseline(static_cast<std::size_t>(chunk), 0.0);
    std::vector<double> delta(static_cast<std::size_t>(chunk), 0.0);
    for (int j = 0; j < chunk; ++j) {
      baseline[static_cast<std::size_t>(j)] = X(0, col0 + j);
      if (!std::isfinite(baseline[static_cast<std::size_t>(j)])) {
        Rcpp::stop("X must contain only finite trait values.");
      }
    }

    // Pass 1: solve Qx through the tree and obtain 1'Qx from terminal edge
    // currents.  This is a direct bilinear calculation; in particular it
    // never obtains the GLS mean by subtracting two large energies.
    upward_messages(
      tree, X, col0, chunk, outgoing, baseline, delta, message
    );
    downward_states(tree, col0, chunk, outgoing, message, state);
    for (int j = 0; j < chunk; ++j) {
      long double qlinear = 0.0L;
      for (int tip = 0; tip < tree.n_tip; ++tip) {
        const int parent = tree.parent[static_cast<std::size_t>(tip)];
        const double x = X(tip, col0 + j);
        if (!std::isfinite(x)) {
          Rcpp::stop("X must contain only finite trait values.");
        }
        const double centered_tip = x -
          baseline[static_cast<std::size_t>(j)];
        const double zp = state[static_cast<std::size_t>(parent) *
          static_cast<std::size_t>(chunk) + static_cast<std::size_t>(j)];
        qlinear += (static_cast<long double>(centered_tip) -
                    static_cast<long double>(zp)) /
          static_cast<long double>(tree.branch[static_cast<std::size_t>(tip)]);
      }
      delta[static_cast<std::size_t>(j)] =
        static_cast<double>(qlinear) / sum_inv;
      const double a = baseline[static_cast<std::size_t>(j)] +
        delta[static_cast<std::size_t>(j)];
      if (!std::isfinite(a)) {
        Rcpp::stop("GLS means are not finite for X.");
      }
      means[col0 + j] = a;
    }

    // Pass 2: run the same Schur complement on residual traits.  The edge
    // energy at the resulting conditional states is (x-a)'Q(x-a), while the
    // unweighted tip residual sum is the numerator used by phytools.
    // Reuse the centered delta from pass 1.  Do not replace it with the full
    // GLS mean: X - (baseline + delta) loses low-order residual bits for
    // offsets around 1e12--1e15.
    upward_messages(
      tree, X, col0, chunk, outgoing, baseline, delta, message
    );
    downward_states(tree, col0, chunk, outgoing, message, state);
    for (int j = 0; j < chunk; ++j) {
      long double num = 0.0L;
      long double den = 0.0L;
      for (int tip = 0; tip < tree.n_tip; ++tip) {
        const double centered_tip = X(tip, col0 + j) -
          baseline[static_cast<std::size_t>(j)];
        const long double y = static_cast<long double>(centered_tip) -
          static_cast<long double>(delta[static_cast<std::size_t>(j)]);
        num += y * y;
      }
      for (std::size_t ii = 1; ii < tree.preorder.size(); ++ii) {
        const int node = tree.preorder[ii];
        const int parent = tree.parent[static_cast<std::size_t>(node)];
        const double zp = state[static_cast<std::size_t>(parent) *
          static_cast<std::size_t>(chunk) + static_cast<std::size_t>(j)];
        const long double child_state = node < tree.n_tip
          ? static_cast<long double>(
              X(node, col0 + j) - baseline[static_cast<std::size_t>(j)]
            ) - static_cast<long double>(delta[static_cast<std::size_t>(j)])
          : static_cast<long double>(state[static_cast<std::size_t>(node) *
              static_cast<std::size_t>(chunk) + static_cast<std::size_t>(j)]);
        const long double d = child_state - static_cast<long double>(zp);
        den += d * d /
          static_cast<long double>(tree.branch[static_cast<std::size_t>(node)]);
      }
      const double num_d = static_cast<double>(num);
      const double den_d = static_cast<double>(den);
      numerator[col0 + j] = num_d;
      denominator[col0 + j] = den_d;
      kval[col0 + j] = (den_d > 0.0 && std::isfinite(den_d))
        ? (num_d / den_d) / norm_const
        : std::numeric_limits<double>::quiet_NaN();
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("K") = kval,
    Rcpp::Named("gls_mean") = means,
    Rcpp::Named("numerator") = numerator,
    Rcpp::Named("denominator") = denominator,
    Rcpp::Named("sum_inv") = sum_inv,
    Rcpp::Named("normalization") = norm_const
  );
}
