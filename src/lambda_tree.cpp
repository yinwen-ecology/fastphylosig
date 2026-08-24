// Fixed-lambda Gaussian pruning kernel.
//
// For a compiled rooted tree this file evaluates the likelihood surface
// implied by phytools::phylosig(method = "lambda", se = NULL):
//
//   C_lambda = lambda * (C - diag(diag(C))) + diag(diag(C)).
//
// The same covariance is obtained by replacing every non-terminal child edge
// by lambda * l and every terminal edge by
// lambda * l + (1 - lambda) * root_distance_tip.  Gaussian messages then
// evaluate the GLS mean, ML sigma2 (with n in the denominator), and logLik in
// linear time in the tree size for each lambda and trait column.  The fixed
// evaluator remains the numerical oracle for the bounded Brent optimizer
// exposed at the end of this file.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

namespace {

typedef long double ld;

struct TreeArrays {
  int n_tip;
  int n_total;
  int root; // zero based
  std::vector<int> parent;
  std::vector<int> child_ptr; // zero based CSR offsets
  std::vector<int> children;  // zero based node ids
  std::vector<int> preorder;
  std::vector<int> postorder;
  std::vector<double> branch;
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
  if (!std::isfinite(out)) Rcpp::stop(label + " must be finite.");
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
                                   const std::string& label,
                                   const bool nonnegative = false) {
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
    if (!std::isfinite(v) || (nonnegative && v < 0.0)) {
      Rcpp::stop(label + " must contain finite" +
                 std::string(nonnegative ? " non-negative" : "") +
                 " values.");
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

// Read the canonical compiled-tree list emitted by compile_tree_cpp().  The
// zero-based CSR variant is accepted as well because it is useful in direct
// low-level tests and serialized contexts.
TreeArrays parse_tree(const Rcpp::List& compiled) {
  if (!has_field(compiled, "n_tip") || !has_field(compiled, "n_total")) {
    Rcpp::stop("compiled_tree must contain n_tip and n_total.");
  }
  const int n_tip = scalar_integer(compiled["n_tip"],
                                   "compiled_tree$n_tip", true);
  const int n_total = scalar_integer(compiled["n_total"],
                                     "compiled_tree$n_total", true);
  if (n_tip < 2 || n_total <= n_tip) {
    Rcpp::stop("compiled_tree must contain at least two tips and one internal node.");
  }

  SEXP root_sexp = first_field(compiled, "root");
  int root = 0;
  if (root_sexp != R_NilValue) {
    const int root_raw = scalar_integer(root_sexp, "compiled_tree$root", false);
    root = root_raw == 0 ? 0 : root_raw - 1;
  }

  SEXP ptr_sexp = first_field(compiled, "child_ptr", "children_ptr",
                              "children_offset");
  SEXP child_sexp = first_field(compiled, "children");
  if (ptr_sexp == R_NilValue || child_sexp == R_NilValue) {
    Rcpp::stop("compiled_tree must contain CSR child_ptr and children arrays.");
  }
  const std::vector<int> raw_ptr = integer_vector(
    ptr_sexp, static_cast<R_xlen_t>(n_total + 1),
    "compiled_tree child_ptr");
  const std::vector<int> raw_children = integer_vector(
    child_sexp, static_cast<R_xlen_t>(n_total - 1),
    "compiled_tree children");

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
    const int child = children_zero_based ? raw_children[i]
                                          : raw_children[i] - 1;
    if (child < 0 || child >= n_total) {
      Rcpp::stop("compiled_tree children contain an out-of-range node.");
    }
    children[i] = child;
  }

  SEXP branch_sexp = first_field(compiled, "branch_length_by_node",
                                 "branch_length");
  if (branch_sexp == R_NilValue) {
    Rcpp::stop("compiled_tree must contain branch_length_by_node.");
  }
  const std::vector<double> branch = numeric_vector(
    branch_sexp, static_cast<R_xlen_t>(n_total),
    "compiled_tree branch_length_by_node", true);

  if (root_sexp == R_NilValue) {
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

  struct Frame { int node; int next; };
  std::vector<char> seen(static_cast<std::size_t>(n_total), 0);
  std::vector<int> preorder;
  std::vector<int> postorder;
  preorder.reserve(static_cast<std::size_t>(n_total));
  postorder.reserve(static_cast<std::size_t>(n_total));
  std::vector<Frame> stack;
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
    const ld distance = static_cast<ld>(root_distance[static_cast<std::size_t>(par)]) +
      static_cast<ld>(branch[static_cast<std::size_t>(node)]);
    if (!std::isfinite(distance) || distance < 0.0L ||
        distance > static_cast<ld>(std::numeric_limits<double>::max())) {
      Rcpp::stop("compiled_tree root distances are not finite.");
    }
    root_distance[static_cast<std::size_t>(node)] = static_cast<double>(distance);
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
  out.branch = branch;
  out.root_distance.swap(root_distance);
  return out;
}

double phytools_max_lambda(const TreeArrays& tree) {
  const double h = tree.root_distance[0];
  if (!std::isfinite(h) || h <= 0.0) return 1.0;
  bool ultrametric = true;
  const double tol = 256.0 * std::numeric_limits<double>::epsilon() *
    std::max(1.0, std::abs(h));
  for (int tip = 1; tip < tree.n_tip; ++tip) {
    if (std::abs(tree.root_distance[static_cast<std::size_t>(tip)] - h) > tol) {
      ultrametric = false;
      break;
    }
  }
  if (!ultrametric) return 1.0;
  double max_parent_height = 0.0;
  for (int node = 0; node < tree.n_total; ++node) {
    if (node == tree.root) continue;
    const int par = tree.parent[static_cast<std::size_t>(node)];
    max_parent_height = std::max(
      max_parent_height, tree.root_distance[static_cast<std::size_t>(par)]
    );
  }
  if (!std::isfinite(max_parent_height) || max_parent_height <= 0.0) {
    return 1.0;
  }
  const double out = h / max_parent_height;
  return std::isfinite(out) && out > 0.0 ? out : 1.0;
}

struct LambdaWeights {
  bool ok;
  std::string status;
  std::vector<ld> edge;
  std::vector<ld> aggregate_a;
  std::vector<ld> aggregate_b1;
  std::vector<ld> aggregate_d1;
  std::vector<ld> a;
  std::vector<ld> b1;
  std::vector<ld> d1;
  std::vector<ld> factor;
  std::vector<ld> k;
  ld k_root;
  ld d1_root;
  ld logdet;
};

LambdaWeights make_weights(const TreeArrays& tree,
                           const double lambda,
                           const double max_lambda) {
  LambdaWeights out;
  out.ok = false;
  out.edge.assign(static_cast<std::size_t>(tree.n_total), 0.0L);
  out.aggregate_a.assign(static_cast<std::size_t>(tree.n_total), 0.0L);
  out.aggregate_b1.assign(static_cast<std::size_t>(tree.n_total), 0.0L);
  out.aggregate_d1.assign(static_cast<std::size_t>(tree.n_total), 0.0L);
  out.a.assign(static_cast<std::size_t>(tree.n_total), 0.0L);
  out.b1.assign(static_cast<std::size_t>(tree.n_total), 0.0L);
  out.d1.assign(static_cast<std::size_t>(tree.n_total), 0.0L);
  out.factor.assign(static_cast<std::size_t>(tree.n_total), 1.0L);
  out.k.assign(static_cast<std::size_t>(tree.n_total), 0.0L);
  if (!std::isfinite(lambda)) {
    out.status = "lambda is not finite";
    return out;
  }
  const double bound_tol = 256.0 * std::numeric_limits<double>::epsilon() *
    std::max(1.0, std::abs(max_lambda));
  if (lambda < 0.0 || lambda > max_lambda + bound_tol) {
    std::ostringstream msg;
    msg << "lambda is outside phytools domain [0, " << max_lambda << "]";
    out.status = msg.str();
    return out;
  }
  // Values rounded just beyond the upper bound are treated as the boundary,
  // not silently extrapolated into a different covariance model.
  const double lam_use = lambda > max_lambda ? max_lambda : lambda;
  const ld lam = static_cast<ld>(lam_use);
  const ld one_minus = 1.0L - lam;

  for (int node = 0; node < tree.n_total; ++node) {
    if (node < tree.n_tip) {
      out.edge[static_cast<std::size_t>(node)] =
        lam * static_cast<ld>(tree.branch[static_cast<std::size_t>(node)]) +
        one_minus * static_cast<ld>(
          tree.root_distance[static_cast<std::size_t>(node)]
        );
    } else if (node != tree.root) {
      out.edge[static_cast<std::size_t>(node)] =
        lam * static_cast<ld>(tree.branch[static_cast<std::size_t>(node)]);
    }
    const ld v = out.edge[static_cast<std::size_t>(node)];
    if (!std::isfinite(v) || v < 0.0L ||
        (node < tree.n_tip && v <= 0.0L)) {
      out.status = node < tree.n_tip
        ? "lambda produces a non-positive terminal variance"
        : "lambda produces a negative internal edge variance";
      return out;
    }
  }

  const ld log_two_pi = std::log(static_cast<ld>(
    6.2831853071795864769252867665590057683943387987502L
  ));
  for (std::size_t ii = 0; ii < tree.postorder.size(); ++ii) {
    const int node = tree.postorder[ii];
    ld a0 = 0.0L;
    ld b10 = 0.0L;
    ld d10 = 0.0L;
    ld k0 = 0.0L;
    if (node < tree.n_tip) {
      const ld v = out.edge[static_cast<std::size_t>(node)];
      if (!(v > 0.0L) || !std::isfinite(v)) {
        out.status = "terminal variance is not finite and positive";
        return out;
      }
      a0 = b10 = d10 = 1.0L / v;
      k0 = -0.5L * (log_two_pi + std::log(v));
      out.aggregate_a[static_cast<std::size_t>(node)] = a0;
      out.aggregate_b1[static_cast<std::size_t>(node)] = b10;
      out.aggregate_d1[static_cast<std::size_t>(node)] = d10;
      out.a[static_cast<std::size_t>(node)] = a0;
      out.b1[static_cast<std::size_t>(node)] = b10;
      out.d1[static_cast<std::size_t>(node)] = d10;
      out.k[static_cast<std::size_t>(node)] = k0;
      continue;
    }

    const int begin = tree.child_ptr[static_cast<std::size_t>(node)];
    const int end = tree.child_ptr[static_cast<std::size_t>(node + 1)];
    for (int kk = begin; kk < end; ++kk) {
      const int child = tree.children[static_cast<std::size_t>(kk)];
      a0 += out.a[static_cast<std::size_t>(child)];
      b10 += out.b1[static_cast<std::size_t>(child)];
      d10 += out.d1[static_cast<std::size_t>(child)];
      k0 += out.k[static_cast<std::size_t>(child)];
    }
    if (!(a0 > 0.0L) || !std::isfinite(a0) || !std::isfinite(b10) ||
        !std::isfinite(d10) || !std::isfinite(k0)) {
      out.status = "tree precision message is not finite";
      return out;
    }
    out.aggregate_a[static_cast<std::size_t>(node)] = a0;
    out.aggregate_b1[static_cast<std::size_t>(node)] = b10;
    out.aggregate_d1[static_cast<std::size_t>(node)] = d10;
    if (node == tree.root) {
      out.a[static_cast<std::size_t>(node)] = a0;
      out.b1[static_cast<std::size_t>(node)] = b10;
      out.d1[static_cast<std::size_t>(node)] = d10;
      out.k[static_cast<std::size_t>(node)] = k0;
      out.k_root = k0;
      out.d1_root = d10;
      continue;
    }
    const ld v = out.edge[static_cast<std::size_t>(node)];
    const ld denominator = 1.0L + a0 * v;
    if (!(denominator > 0.0L) || !std::isfinite(denominator)) {
      out.status = "tree message denominator is not finite and positive";
      return out;
    }
    const ld inv_denominator = 1.0L / denominator;
    ld d1out = d10 - v * b10 * b10 * inv_denominator;
    if (!std::isfinite(d1out)) {
      out.status = "tree GLS precision message is not finite";
      return out;
    }
    if (d1out < 0.0L) {
      const ld scale = std::max<ld>(1.0L, std::max(std::abs(d10),
                                                    std::abs(v * b10 * b10 *
                                                             inv_denominator)));
      if (std::abs(d1out) > 1024.0L * std::numeric_limits<ld>::epsilon() *
          scale) {
        out.status = "tree GLS precision message is negative";
        return out;
      }
      d1out = 0.0L;
    }
    out.factor[static_cast<std::size_t>(node)] = inv_denominator;
    out.a[static_cast<std::size_t>(node)] = a0 * inv_denominator;
    out.b1[static_cast<std::size_t>(node)] = b10 * inv_denominator;
    out.d1[static_cast<std::size_t>(node)] = d1out;
    out.k[static_cast<std::size_t>(node)] =
      k0 - 0.5L * std::log(denominator);
  }

  if (!(out.d1_root > 0.0L) || !std::isfinite(out.d1_root) ||
      !std::isfinite(out.k_root)) {
    out.status = "1'Q1 is not finite and positive";
    return out;
  }
  out.logdet = -2.0L * out.k_root -
    static_cast<ld>(tree.n_tip) * log_two_pi;
  if (!std::isfinite(out.logdet) ||
      out.logdet > static_cast<ld>(std::numeric_limits<double>::max()) ||
      out.logdet < -static_cast<ld>(std::numeric_limits<double>::max())) {
    out.status = "lambda covariance log determinant is not finite";
    return out;
  }
  out.ok = true;
  out.status = "ok";
  return out;
}

// First pass for one trait chunk.  It computes 1'Qz for z = X - baseline and
// propagates the mixed bilinear message z'Q1; keeping these quantities
// separate avoids subtracting a large common offset from the final GLS mean.
void mixed_pass(const TreeArrays& tree,
                const LambdaWeights& w,
                const Rcpp::NumericMatrix& X,
                const int col0,
                const int chunk,
                const std::vector<ld>& baseline,
                std::vector<ld>& by,
                std::vector<ld>& ey,
                ld& cross_root) {
  const std::size_t stride = static_cast<std::size_t>(chunk);
  for (std::size_t ii = 0; ii < tree.postorder.size(); ++ii) {
    const int node = tree.postorder[ii];
    const std::size_t base = static_cast<std::size_t>(node) * stride;
    if (node < tree.n_tip) {
      const ld v = w.edge[static_cast<std::size_t>(node)];
      for (int j = 0; j < chunk; ++j) {
        const ld z = static_cast<ld>(X(node, col0 + j)) -
          baseline[static_cast<std::size_t>(j)];
        const ld b = z / v;
        by[base + static_cast<std::size_t>(j)] = b;
        ey[base + static_cast<std::size_t>(j)] = b;
      }
      continue;
    }
    const int begin = tree.child_ptr[static_cast<std::size_t>(node)];
    const int end = tree.child_ptr[static_cast<std::size_t>(node + 1)];
    const bool is_root = node == tree.root;
    const ld f = is_root ? 1.0L : w.factor[static_cast<std::size_t>(node)];
    const ld v = is_root ? 0.0L : w.edge[static_cast<std::size_t>(node)];
    const ld b10 = w.aggregate_b1[static_cast<std::size_t>(node)];
    for (int j = 0; j < chunk; ++j) {
      ld by0 = 0.0L;
      ld ey0 = 0.0L;
      for (int kk = begin; kk < end; ++kk) {
        const int child = tree.children[static_cast<std::size_t>(kk)];
        const std::size_t child_base = static_cast<std::size_t>(child) * stride;
        by0 += by[child_base + static_cast<std::size_t>(j)];
        ey0 += ey[child_base + static_cast<std::size_t>(j)];
      }
      by[base + static_cast<std::size_t>(j)] = by0 * f;
      ey[base + static_cast<std::size_t>(j)] =
        is_root ? ey0 : ey0 - v * by0 * b10 * f;
      if (!std::isfinite(by[base + static_cast<std::size_t>(j)]) ||
          !std::isfinite(ey[base + static_cast<std::size_t>(j)])) {
        Rcpp::stop("lambda pruning messages are not finite.");
      }
    }
  }
  cross_root = 0.0L;
  const std::size_t root_base = static_cast<std::size_t>(tree.root) * stride;
  // The root stores one value per trait.  The caller reads ey directly after
  // this routine; this assignment is kept for a scalar diagnostic API.
  (void)cross_root;
  (void)root_base;
}

// Second pass on residual tips.  Unlike qz - cross^2/q1, this direct residual
// quadratic avoids cancellation for near-constant traits and large offsets.
void residual_pass(const TreeArrays& tree,
                   const LambdaWeights& w,
                   const Rcpp::NumericMatrix& X,
                   const int col0,
                   const int chunk,
                   const std::vector<ld>& baseline,
                   const std::vector<ld>& delta,
                   std::vector<ld>& by,
                   std::vector<ld>& dy,
                   std::vector<ld>& q_root) {
  const std::size_t stride = static_cast<std::size_t>(chunk);
  for (std::size_t ii = 0; ii < tree.postorder.size(); ++ii) {
    const int node = tree.postorder[ii];
    const std::size_t base = static_cast<std::size_t>(node) * stride;
    if (node < tree.n_tip) {
      const ld v = w.edge[static_cast<std::size_t>(node)];
      for (int j = 0; j < chunk; ++j) {
        const ld z = (static_cast<ld>(X(node, col0 + j)) -
                      baseline[static_cast<std::size_t>(j)]) -
          delta[static_cast<std::size_t>(j)];
        by[base + static_cast<std::size_t>(j)] = z / v;
        dy[base + static_cast<std::size_t>(j)] = z * z / v;
      }
      continue;
    }
    const int begin = tree.child_ptr[static_cast<std::size_t>(node)];
    const int end = tree.child_ptr[static_cast<std::size_t>(node + 1)];
    const bool is_root = node == tree.root;
    const ld f = is_root ? 1.0L : w.factor[static_cast<std::size_t>(node)];
    const ld v = is_root ? 0.0L : w.edge[static_cast<std::size_t>(node)];
    const ld a0 = w.aggregate_a[static_cast<std::size_t>(node)];
    for (int j = 0; j < chunk; ++j) {
      ld by0 = 0.0L;
      ld dy0 = 0.0L;
      for (int kk = begin; kk < end; ++kk) {
        const int child = tree.children[static_cast<std::size_t>(kk)];
        const std::size_t child_base = static_cast<std::size_t>(child) * stride;
        by0 += by[child_base + static_cast<std::size_t>(j)];
        dy0 += dy[child_base + static_cast<std::size_t>(j)];
      }
      by[base + static_cast<std::size_t>(j)] = by0 * f;
      dy[base + static_cast<std::size_t>(j)] =
        is_root ? dy0 : dy0 - v * by0 * by0 * f;
      if (!std::isfinite(by[base + static_cast<std::size_t>(j)]) ||
          !std::isfinite(dy[base + static_cast<std::size_t>(j)]) ||
          (dy[base + static_cast<std::size_t>(j)] < 0.0L &&
           std::abs(dy[base + static_cast<std::size_t>(j)]) >
             1024.0L * std::numeric_limits<ld>::epsilon() *
             std::max<ld>(1.0L, std::abs(dy0)))) {
        Rcpp::stop("lambda residual quadratic is not finite.");
      }
      // Roundoff can leave a mathematically zero Schur complement a tiny bit
      // negative.  It is diagnosed at the root, where zero means a constant
      // trait, but harmless local negatives are clipped here.
      if (dy[base + static_cast<std::size_t>(j)] < 0.0L) {
        dy[base + static_cast<std::size_t>(j)] = 0.0L;
      }
    }
    (void)a0;
  }
  q_root.assign(static_cast<std::size_t>(chunk), 0.0L);
  const std::size_t root_base = static_cast<std::size_t>(tree.root) * stride;
  for (int j = 0; j < chunk; ++j) {
    q_root[static_cast<std::size_t>(j)] =
      dy[root_base + static_cast<std::size_t>(j)];
  }
}

struct LambdaTraitPoint {
  double gls_mean;
  double sigma2;
  double log_lik;
  bool valid;
  std::string status;
};

// Scalar fixed-surface evaluator used by the optimizer.  Unlike the public
// batch kernel, this deliberately evaluates one trait at a time so that Brent
// searches for p traits cost O(n * p * evaluations), rather than evaluating
// all p traits at every trait-specific lambda candidate.  The pruning passes
// are exactly the same mixed/residual passes used by the fixed batch entry
// point above.
LambdaTraitPoint evaluate_lambda_trait(const TreeArrays& tree,
                                      const Rcpp::NumericMatrix& X,
                                      const int trait,
                                      const double lambda,
                                      const double max_lambda) {
  LambdaTraitPoint out;
  out.gls_mean = NA_REAL;
  out.sigma2 = NA_REAL;
  out.log_lik = R_NegInf;
  out.valid = false;
  out.status = "undefined";

  const LambdaWeights w = make_weights(tree, lambda, max_lambda);
  if (!w.ok) {
    out.status = w.status;
    return out;
  }
  const int chunk = 1;
  const std::size_t cells = static_cast<std::size_t>(tree.n_total);
  std::vector<ld> by(cells, 0.0L);
  std::vector<ld> ey(cells, 0.0L);
  std::vector<ld> baseline(1, static_cast<ld>(X(0, trait)));
  ld unused_cross = 0.0L;
  mixed_pass(tree, w, X, trait, chunk, baseline, by, ey, unused_cross);

  const std::size_t root = static_cast<std::size_t>(tree.root);
  const ld delta = ey[root] / w.d1_root;
  const ld mean = baseline[0] + delta;
  if (!std::isfinite(delta) || !std::isfinite(mean) ||
      mean > static_cast<ld>(std::numeric_limits<double>::max()) ||
      mean < -static_cast<ld>(std::numeric_limits<double>::max())) {
    out.status = "GLS mean is not finite";
    return out;
  }
  out.gls_mean = static_cast<double>(mean);

  std::vector<ld> dy(cells, 0.0L);
  std::vector<ld> delta_vec(1, delta);
  std::vector<ld> q_root;
  residual_pass(tree, w, X, trait, chunk, baseline, delta_vec, by, dy, q_root);
  const ld q = q_root[0];
  const ld sig = q / static_cast<ld>(tree.n_tip);
  if (!std::isfinite(q) || q <= 0.0L || !std::isfinite(sig) || sig <= 0.0L) {
    out.status = "non-positive residual quadratic (constant trait)";
    return out;
  }
  const ld log_two_pi = std::log(static_cast<ld>(
    6.2831853071795864769252867665590057683943387987502L
  ));
  const ld ll = -0.5L * q / sig -
    0.5L * static_cast<ld>(tree.n_tip) * log_two_pi -
    0.5L * (static_cast<ld>(tree.n_tip) * std::log(sig) + w.logdet);
  if (!std::isfinite(ll) ||
      sig > static_cast<ld>(std::numeric_limits<double>::max()) ||
      ll < -static_cast<ld>(std::numeric_limits<double>::max()) ||
      ll > static_cast<ld>(std::numeric_limits<double>::max())) {
    out.status = "sigma2 or logLik is not finite";
    return out;
  }
  out.sigma2 = static_cast<double>(sig);
  out.log_lik = static_cast<double>(ll);
  out.valid = true;
  out.status = "ok";
  return out;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List fast_lambda_tree_fixed_cpp(
    const Rcpp::List& compiled_tree,
    const Rcpp::NumericMatrix& X,
    const Rcpp::NumericVector& lambda,
    const int trait_chunk = 64,
    const int n_threads = 1) {
  const TreeArrays tree = parse_tree(compiled_tree);
  if (X.nrow() != tree.n_tip) {
    Rcpp::stop("X must have one row per compiled_tree tip (in tip-id order).");
  }
  if (X.ncol() < 1) Rcpp::stop("X must contain at least one trait column.");
  if (lambda.size() < 1) Rcpp::stop("lambda must contain at least one value.");
  if (trait_chunk < 1) Rcpp::stop("trait_chunk must be a positive integer.");
  if (n_threads < 1) Rcpp::stop("n_threads must be a positive integer.");
  (void)n_threads; // The first fixed-surface kernel is deterministic and serial.

  const int p = X.ncol();
  const int n_lambda = lambda.size();
  for (int i = 0; i < tree.n_tip; ++i) {
    for (int j = 0; j < p; ++j) {
      if (!std::isfinite(X(i, j))) {
        Rcpp::stop("X must contain only finite trait values.");
      }
    }
  }

  const double max_lambda = phytools_max_lambda(tree);
  Rcpp::NumericMatrix means(n_lambda, p);
  Rcpp::NumericMatrix sigma2(n_lambda, p);
  Rcpp::NumericMatrix log_lik(n_lambda, p);
  Rcpp::LogicalMatrix valid(n_lambda, p);
  Rcpp::CharacterMatrix status(n_lambda, p);
  Rcpp::LogicalVector lambda_valid(n_lambda);
  Rcpp::CharacterVector lambda_status(n_lambda);
  const double na = NA_REAL;

  for (int li = 0; li < n_lambda; ++li) {
    const double lam = lambda[li];
    const LambdaWeights w = make_weights(tree, lam, max_lambda);
    if (!w.ok) {
      lambda_valid[li] = false;
      lambda_status[li] = w.status;
      for (int j = 0; j < p; ++j) {
        means(li, j) = na;
        sigma2(li, j) = na;
        log_lik(li, j) = R_NegInf;
        valid(li, j) = false;
        status(li, j) = w.status;
      }
      continue;
    }

    bool row_ok = true;
    std::string row_status = "ok";
    const int chunk_limit = std::max(1, std::min(trait_chunk, p));
    for (int col0 = 0; col0 < p; col0 += chunk_limit) {
      const int chunk = std::min(chunk_limit, p - col0);
      const std::size_t cells = static_cast<std::size_t>(tree.n_total) *
        static_cast<std::size_t>(chunk);
      std::vector<ld> by(cells, 0.0L);
      std::vector<ld> ey(cells, 0.0L);
      std::vector<ld> baseline(static_cast<std::size_t>(chunk), 0.0L);
      for (int j = 0; j < chunk; ++j) {
        baseline[static_cast<std::size_t>(j)] =
          static_cast<ld>(X(0, col0 + j));
      }
      ld unused_cross = 0.0L;
      mixed_pass(tree, w, X, col0, chunk, baseline, by, ey, unused_cross);

      std::vector<ld> delta(static_cast<std::size_t>(chunk), 0.0L);
      for (int j = 0; j < chunk; ++j) {
        const std::size_t root_base = static_cast<std::size_t>(tree.root) *
          static_cast<std::size_t>(chunk);
        const ld cross = ey[root_base + static_cast<std::size_t>(j)];
        const ld d = w.d1_root;
        const ld delta_j = cross / d;
        delta[static_cast<std::size_t>(j)] = delta_j;
        const ld mean = baseline[static_cast<std::size_t>(j)] + delta_j;
        if (!std::isfinite(delta_j) || !std::isfinite(mean) ||
            mean > static_cast<ld>(std::numeric_limits<double>::max()) ||
            mean < -static_cast<ld>(std::numeric_limits<double>::max())) {
          means(li, col0 + j) = na;
          sigma2(li, col0 + j) = na;
          log_lik(li, col0 + j) = R_NegInf;
          valid(li, col0 + j) = false;
          status(li, col0 + j) = "GLS mean is not finite";
          row_ok = false;
          row_status = "one or more traits are undefined";
        } else {
          means(li, col0 + j) = static_cast<double>(mean);
        }
      }

      std::vector<ld> dy(cells, 0.0L);
      std::vector<ld> q_root;
      residual_pass(tree, w, X, col0, chunk, baseline, delta, by, dy, q_root);
      const ld log_two_pi = std::log(static_cast<ld>(
        6.2831853071795864769252867665590057683943387987502L
      ));
      for (int j = 0; j < chunk; ++j) {
        if (!std::isfinite(delta[static_cast<std::size_t>(j)])) continue;
        const ld q = q_root[static_cast<std::size_t>(j)];
        const ld sig = q / static_cast<ld>(tree.n_tip);
        const bool q_ok = std::isfinite(q) && q > 0.0L &&
          std::isfinite(sig) && sig > 0.0L;
        if (!q_ok) {
          sigma2(li, col0 + j) = na;
          log_lik(li, col0 + j) = R_NegInf;
          valid(li, col0 + j) = false;
          status(li, col0 + j) = "non-positive residual quadratic (constant trait)";
          row_ok = false;
          row_status = "one or more traits are undefined";
          continue;
        }
        const ld log_sig = std::log(sig);
        const ld ll = -0.5L * q / sig -
          0.5L * static_cast<ld>(tree.n_tip) * log_two_pi -
          0.5L * (static_cast<ld>(tree.n_tip) * log_sig + w.logdet);
        if (!std::isfinite(ll) ||
            sig > static_cast<ld>(std::numeric_limits<double>::max()) ||
            ll < -static_cast<ld>(std::numeric_limits<double>::max()) ||
            ll > static_cast<ld>(std::numeric_limits<double>::max())) {
          sigma2(li, col0 + j) = na;
          log_lik(li, col0 + j) = R_NegInf;
          valid(li, col0 + j) = false;
          status(li, col0 + j) = "sigma2 or logLik is not finite";
          row_ok = false;
          row_status = "one or more traits are undefined";
          continue;
        }
        sigma2(li, col0 + j) = static_cast<double>(sig);
        log_lik(li, col0 + j) = static_cast<double>(ll);
        valid(li, col0 + j) = true;
        status(li, col0 + j) = "ok";
      }
    }
    lambda_valid[li] = row_ok;
    lambda_status[li] = row_status;
  }

  return Rcpp::List::create(
    Rcpp::Named("lambda") = lambda,
    Rcpp::Named("gls_mean") = means,
    Rcpp::Named("sigma2") = sigma2,
    Rcpp::Named("logLik") = log_lik,
    Rcpp::Named("valid") = valid,
    Rcpp::Named("status") = status,
    Rcpp::Named("lambda_valid") = lambda_valid,
    Rcpp::Named("lambda_status") = lambda_status,
    Rcpp::Named("max_lambda") = max_lambda
  );
}

// Cache scalar fixed-surface evaluations by (lambda, trait).  Brent searches
// each trait independently, so evaluating all traits at every candidate would
// introduce an avoidable O(p^2) factor.
class LambdaTraitCache {
 public:
  LambdaTraitCache(const TreeArrays& tree,
                     const Rcpp::NumericMatrix& X,
                     const double max_lambda,
                     const int trait_chunk,
                     const int n_threads)
    : tree_(tree), X_(X), max_lambda_(max_lambda),
      trait_chunk_(trait_chunk), n_threads_(n_threads) {
    // Scalar optimizer evaluations do not use chunking or OpenMP yet; retain
    // the controls in the cache so the public signature stays aligned with
    // the fixed batch kernel while keeping this first optimizer deterministic.
  }

  const LambdaTraitPoint& eval(const double lambda, const int trait) {
    const std::pair<double, int> key(lambda, trait);
    std::map<std::pair<double, int>, LambdaTraitPoint>::const_iterator hit =
      cache_.find(key);
    if (hit != cache_.end()) return hit->second;

    const LambdaTraitPoint point = evaluate_lambda_trait(
      tree_, X_, trait, lambda, max_lambda_
    );
    const std::pair<
      std::map<std::pair<double, int>, LambdaTraitPoint>::iterator, bool
    > put = cache_.insert(std::make_pair(key, point));
    return put.first->second;
  }

 private:
  const TreeArrays& tree_;
  const Rcpp::NumericMatrix& X_;
  double max_lambda_;
  int trait_chunk_;
  int n_threads_;
  std::map<std::pair<double, int>, LambdaTraitPoint> cache_;
};

// Brent's bounded maximizer, written in the same tolerance regime as
// stats::optimize().  The implementation minimizes -logLik; +Inf denotes an
// invalid fixed-surface evaluation and is never selected as an optimum.
double brent_minimize(const double lower,
                      const double upper,
                      const double tol,
                      const std::function<double(double)>& objective) {
  const double cgold = 0.3819660112501051518;
  const double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());
  const int max_iter = 1000;
  double a = lower;
  double b = upper;
  double x = a + cgold * (b - a);
  double w = x;
  double v = x;
  double fx = objective(x);
  double fw = fx;
  double fv = fx;
  double d = 0.0;
  double e = 0.0;

  for (int iter = 0; iter < max_iter; ++iter) {
    const double xm = 0.5 * (a + b);
    // This is R's Brent_fmin stopping scale (src/library/stats/src/optimize.c):
    // relative sqrt(eps) term plus one third of the user tolerance.  Using
    // tol*abs(x) here is subtly different near a boundary and shifts lambda
    // estimates by several 1e-5 for phytools' default tolerance.
    const double tol1 = sqrt_eps * std::abs(x) + tol / 3.0;
    const double tol2 = 2.0 * tol1;
    if (std::abs(x - xm) <= tol2 - 0.5 * (b - a)) break;

    bool parabolic = false;
    if (std::abs(e) > tol1 && std::isfinite(fx) &&
        std::isfinite(fw) && std::isfinite(fv)) {
      const double r = (x - w) * (fx - fv);
      const double q = (x - v) * (fx - fw);
      double p = (x - v) * q - (x - w) * r;
      double q2 = 2.0 * (q - r);
      if (q2 > 0.0) p = -p;
      q2 = std::abs(q2);
      const double etemp = e;
      e = d;
      if (q2 > 0.0 && std::abs(p) <
          std::abs(0.5 * q2 * etemp) &&
          p > q2 * (a - x) && p < q2 * (b - x)) {
        d = p / q2;
        const double u = x + d;
        if (u - a < tol2 || b - u < tol2) {
          d = (xm >= x ? tol1 : -tol1);
        }
        parabolic = true;
      }
    }
    if (!parabolic) {
      e = (x >= xm) ? (a - x) : (b - x);
      d = cgold * e;
    }
    const double u = x + (std::abs(d) >= tol1 ? d :
                           (d >= 0.0 ? tol1 : -tol1));
    const double fu = objective(u);
    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }
  return x;
}

double chi_square1_upper_tail(const double statistic) {
  if (!std::isfinite(statistic) || statistic < 0.0) return NA_REAL;
  return std::erfc(std::sqrt(0.5 * statistic));
}

// [[Rcpp::export]]
Rcpp::List fast_lambda_tree_optimize_cpp(
    const Rcpp::List& compiled_tree,
    const Rcpp::NumericMatrix& X,
    SEXP max_lambda = R_NilValue,
    const bool test = true,
    const bool profile = false,
    const int profile_points = 101,
    const int trait_chunk = 64,
    const int n_threads = 1,
    const double tol = 0.0001220703125) {
  const TreeArrays tree = parse_tree(compiled_tree);
  if (X.nrow() != tree.n_tip) {
    Rcpp::stop("X must have one row per compiled_tree tip (in tip-id order).");
  }
  if (X.ncol() < 1) Rcpp::stop("X must contain at least one trait column.");
  if (trait_chunk < 1) Rcpp::stop("trait_chunk must be a positive integer.");
  if (n_threads < 1) Rcpp::stop("n_threads must be a positive integer.");
  if (!std::isfinite(tol) || tol <= 0.0) {
    Rcpp::stop("tol must be finite and positive.");
  }
  if (profile && profile_points < 2) {
    Rcpp::stop("profile_points must be at least two when profile is TRUE.");
  }
  for (int i = 0; i < tree.n_tip; ++i) {
    for (int j = 0; j < X.ncol(); ++j) {
      if (!std::isfinite(X(i, j))) {
        Rcpp::stop("X must contain only finite trait values.");
      }
    }
  }

  const double tree_max_lambda = phytools_max_lambda(tree);
  double upper = tree_max_lambda;
  if (max_lambda != R_NilValue) {
    if (Rf_xlength(max_lambda) != 1 ||
        (TYPEOF(max_lambda) != INTSXP && TYPEOF(max_lambda) != REALSXP)) {
      Rcpp::stop("max_lambda must be NULL or a single numeric value.");
    }
    const double requested = TYPEOF(max_lambda) == INTSXP
      ? static_cast<double>(INTEGER(max_lambda)[0])
      : REAL(max_lambda)[0];
    if (!std::isfinite(requested) || requested <= 0.0) {
      Rcpp::stop("max_lambda must be finite and positive.");
    }
    const double boundary_tol = 256.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, std::abs(tree_max_lambda));
    if (requested > tree_max_lambda + boundary_tol) {
      Rcpp::stop("max_lambda exceeds the phytools tree boundary.");
    }
    upper = std::min(requested, tree_max_lambda);
  }
  if (!(upper > 0.0) || !std::isfinite(upper)) {
    Rcpp::stop("the lambda optimization interval is not finite and positive.");
  }

  const int p = X.ncol();
  LambdaTraitCache cache(tree, X, upper, trait_chunk, n_threads);

  Rcpp::NumericVector lambda_hat(p, NA_REAL);
  Rcpp::NumericVector log_lik_hat(p, R_NegInf);
  Rcpp::NumericVector means(p, NA_REAL);
  Rcpp::NumericVector scales(p, NA_REAL);
  Rcpp::NumericVector log_lik_zero(p, R_NegInf);
  Rcpp::NumericVector lr(p, NA_REAL);
  Rcpp::NumericVector p_value(p, NA_REAL);
  Rcpp::LogicalVector valid(p, false);
  Rcpp::CharacterVector status(p);

  for (int j = 0; j < p; ++j) {
    const LambdaTraitPoint& at_zero = cache.eval(0.0, j);
    log_lik_zero[j] = at_zero.log_lik;
    const std::function<double(double)> objective =
      [&cache, j](const double x) {
        const LambdaTraitPoint& point = cache.eval(x, j);
        const double value = point.log_lik;
        return std::isfinite(value) && point.valid ? -value :
          std::numeric_limits<double>::infinity();
      };

    double candidate = brent_minimize(0.0, upper, tol, objective);
    candidate = std::max(0.0, std::min(upper, candidate));

    // Keep the Brent interior candidate as stats::optimize() does.  The
    // lambda=0 endpoint is evaluated separately for logLik0/LR and is never
    // substituted as an exact estimate; an invalid upper boundary remains
    // +Inf in the minimizer and therefore cannot be selected.
    if (!std::isfinite(candidate)) {
      status[j] = "lambda likelihood is undefined on this trait/tree";
      continue;
    }
    const LambdaTraitPoint& best = cache.eval(candidate, j);
    means[j] = best.gls_mean;
    if (!best.valid || !std::isfinite(best.log_lik)) {
      status[j] = "lambda likelihood is undefined on this trait/tree";
      continue;
    }
    lambda_hat[j] = candidate;
    log_lik_hat[j] = best.log_lik;
    scales[j] = best.sigma2;
    valid[j] = true;
    status[j] = best.status;
    if (test && std::isfinite(log_lik_zero[j])) {
      const double statistic = std::max(0.0,
        2.0 * (log_lik_hat[j] - log_lik_zero[j]));
      lr[j] = statistic;
      p_value[j] = chi_square1_upper_tail(statistic);
    }
  }

  Rcpp::RObject profile_obj = R_NilValue;
  Rcpp::RObject profile_lambda_obj = R_NilValue;
  Rcpp::RObject profile_log_lik_obj = R_NilValue;
  if (profile) {
    std::vector<double> grid;
    grid.reserve(static_cast<std::size_t>(profile_points) +
                static_cast<std::size_t>(p) + 3U);
    for (int i = 0; i < profile_points; ++i) {
      grid.push_back(upper * static_cast<double>(i) /
                     static_cast<double>(profile_points - 1));
    }
    if (upper >= 1.0) grid.push_back(1.0);
    grid.push_back(0.0);
    grid.push_back(upper);
    for (int j = 0; j < p; ++j) {
      if (std::isfinite(lambda_hat[j])) grid.push_back(lambda_hat[j]);
    }
    std::sort(grid.begin(), grid.end());
    std::vector<double> unique_grid;
    unique_grid.reserve(grid.size());
    for (std::size_t i = 0; i < grid.size(); ++i) {
      if (unique_grid.empty() ||
          grid[i] != unique_grid.back()) unique_grid.push_back(grid[i]);
    }
    Rcpp::NumericVector profile_lambda(unique_grid.size());
    Rcpp::NumericMatrix profile_log_lik(
      static_cast<int>(unique_grid.size()), p
    );
    for (std::size_t i = 0; i < unique_grid.size(); ++i) {
      profile_lambda[static_cast<R_xlen_t>(i)] = unique_grid[i];
      for (int j = 0; j < p; ++j) {
        profile_log_lik(static_cast<int>(i), j) =
          cache.eval(unique_grid[i], j).log_lik;
      }
    }
    Rcpp::List profile_list = Rcpp::List::create(
      Rcpp::Named("lambda") = profile_lambda,
      Rcpp::Named("logLik") = profile_log_lik
    );
    profile_obj = profile_list;
    profile_lambda_obj = profile_lambda;
    profile_log_lik_obj = profile_log_lik;
  }

  return Rcpp::List::create(
    Rcpp::Named("lambda") = lambda_hat,
    Rcpp::Named("lambda_hat") = lambda_hat,
    Rcpp::Named("logLik") = log_lik_hat,
    Rcpp::Named("logL") = log_lik_hat,
    Rcpp::Named("gls_mean") = means,
    Rcpp::Named("sigma2") = scales,
    Rcpp::Named("logLik0") = log_lik_zero,
    Rcpp::Named("logL0") = log_lik_zero,
    Rcpp::Named("LR") = lr,
    Rcpp::Named("P") = p_value,
    Rcpp::Named("valid") = valid,
    Rcpp::Named("status") = status,
    Rcpp::Named("max_lambda") = upper,
    Rcpp::Named("profile") = profile_obj,
    Rcpp::Named("profile_lambda") = profile_lambda_obj,
    Rcpp::Named("profile_logLik") = profile_log_lik_obj
  );
}
