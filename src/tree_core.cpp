// Compiled representation and validation for rooted phylogenies.
//
// The representation is intentionally a plain R list.  It contains the
// topology in CSR form and the traversal/root-distance arrays needed by later
// kernels, while remaining safe to copy and serialize from prepare_tree().

#include "tree_core.h"

#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

namespace {

struct ParsedTree {
  int n_tip;
  int n_total;
  int n_node;
  int root;
  std::vector<int> parent;
  std::vector<std::vector<int> > children;
  std::vector<double> edge_length;
  std::vector<double> branch_length_by_node;
  std::vector<int> edge_index_by_node;
  std::vector<int> preorder;
  std::vector<int> postorder;
  std::vector<double> root_distance;
  std::vector<int> tip_order;
};

bool finite_integer_value(const double value) {
  return std::isfinite(value) &&
    value >= 1.0 &&
    value <= static_cast<double>(std::numeric_limits<int>::max()) &&
    std::floor(value) == value;
}

int parse_n_tip(SEXP n_tip) {
  if (Rf_xlength(n_tip) != 1 ||
      (TYPEOF(n_tip) != INTSXP && TYPEOF(n_tip) != REALSXP)) {
    Rcpp::stop("n_tip must be a single positive integer.");
  }

  const double value = TYPEOF(n_tip) == INTSXP
    ? static_cast<double>(INTEGER(n_tip)[0])
    : REAL(n_tip)[0];
  if (!finite_integer_value(value)) {
    Rcpp::stop("n_tip must be a single positive integer.");
  }
  return static_cast<int>(value);
}

std::vector<int> parse_edge(SEXP edge, int* n_edge) {
  if (!Rf_isMatrix(edge) ||
      (TYPEOF(edge) != INTSXP && TYPEOF(edge) != REALSXP)) {
    Rcpp::stop("edge must be a two-column integer matrix.");
  }

  SEXP dim = Rf_getAttrib(edge, R_DimSymbol);
  if (dim == R_NilValue || Rf_xlength(dim) != 2) {
    Rcpp::stop("edge must be a two-column integer matrix.");
  }
  const R_xlen_t nr = INTEGER(dim)[0];
  const R_xlen_t nc = INTEGER(dim)[1];
  if (nc != 2) {
    Rcpp::stop("edge must have exactly two columns (parent, child).");
  }
  if (nr < 1) {
    Rcpp::stop("edge must contain at least one row.");
  }
  if (nr > static_cast<R_xlen_t>(std::numeric_limits<int>::max())) {
    Rcpp::stop("edge has too many rows for an integer representation.");
  }

  std::vector<int> values(static_cast<std::size_t>(nr) * 2U);
  for (R_xlen_t row = 0; row < nr; ++row) {
    for (R_xlen_t col = 0; col < 2; ++col) {
      const R_xlen_t offset = row + nr * col;
      double value = 0.0;
      if (TYPEOF(edge) == INTSXP) {
        const int raw = INTEGER(edge)[offset];
        if (raw == NA_INTEGER) {
          Rcpp::stop("edge contains missing node numbers.");
        }
        value = static_cast<double>(raw);
      } else {
        value = REAL(edge)[offset];
      }
      if (!finite_integer_value(value)) {
        Rcpp::stop("edge node numbers must be finite positive integers.");
      }
      values[static_cast<std::size_t>(offset)] = static_cast<int>(value);
    }
  }
  *n_edge = static_cast<int>(nr);
  return values;
}

std::vector<double> parse_edge_length(SEXP edge_length, int n_edge) {
  if ((TYPEOF(edge_length) != REALSXP && TYPEOF(edge_length) != INTSXP) ||
      Rf_xlength(edge_length) != static_cast<R_xlen_t>(n_edge)) {
    Rcpp::stop("edge_length must be a numeric vector with one value per edge.");
  }

  std::vector<double> out(static_cast<std::size_t>(n_edge));
  for (int i = 0; i < n_edge; ++i) {
    const double value = TYPEOF(edge_length) == REALSXP
      ? REAL(edge_length)[i]
      : static_cast<double>(INTEGER(edge_length)[i]);
    if (!std::isfinite(value) || value < 0.0) {
      Rcpp::stop("edge_length must contain finite non-negative values.");
    }
    out[static_cast<std::size_t>(i)] = value;
  }
  return out;
}

ParsedTree parse_and_compile(SEXP edge, SEXP edge_length, SEXP n_tip_sexp) {
  const int n_tip = parse_n_tip(n_tip_sexp);

  int n_edge = 0;
  const std::vector<int> edge_values = parse_edge(edge, &n_edge);
  const std::vector<double> lengths = parse_edge_length(edge_length, n_edge);

  int max_node = 0;
  for (std::size_t i = 0; i < edge_values.size(); ++i) {
    max_node = std::max(max_node, edge_values[i]);
  }
  if (max_node <= n_tip) {
    Rcpp::stop(
      "edge must include at least one internal node (node ids n_tip + 1 onward)."
    );
  }
  if (max_node < 1) {
    Rcpp::stop("edge node numbers must be positive.");
  }
  const int n_total = max_node;
  if (n_total > std::numeric_limits<int>::max() - 1) {
    Rcpp::stop("tree has too many nodes for an integer representation.");
  }
  const int n_node = n_total - n_tip;
  if (n_node < 1) {
    Rcpp::stop("tree must contain at least one internal node.");
  }
  if (n_edge != n_total - 1) {
    Rcpp::stop(
      "a rooted tree must have exactly n_total - 1 edges; edge topology is invalid."
    );
  }

  std::vector<int> indegree(static_cast<std::size_t>(n_total + 1), 0);
  std::vector<int> outdegree(static_cast<std::size_t>(n_total + 1), 0);
  std::vector<char> seen(static_cast<std::size_t>(n_total + 1), 0);
  std::vector<int> parent(static_cast<std::size_t>(n_total + 1), 0);
  std::vector<std::vector<int> > children(
    static_cast<std::size_t>(n_total + 1)
  );
  std::vector<double> branch_length_by_node(
    static_cast<std::size_t>(n_total + 1), 0.0
  );
  std::vector<int> edge_index_by_node(
    static_cast<std::size_t>(n_total + 1), 0
  );

  for (int i = 0; i < n_edge; ++i) {
    const int p = edge_values[static_cast<std::size_t>(i)];
    const int c = edge_values[static_cast<std::size_t>(i + n_edge)];
    if (p < 1 || p > n_total || c < 1 || c > n_total) {
      Rcpp::stop("edge node numbers are outside the contiguous node range.");
    }
    if (p == c) {
      Rcpp::stop("edge contains a self-loop; tree topology is invalid.");
    }
    if (indegree[static_cast<std::size_t>(c)] != 0) {
      Rcpp::stop(
        "each node must have at most one parent; duplicate child in edge."
      );
    }
    indegree[static_cast<std::size_t>(c)] = 1;
    outdegree[static_cast<std::size_t>(p)] += 1;
    parent[static_cast<std::size_t>(c)] = p;
    children[static_cast<std::size_t>(p)].push_back(c);
    seen[static_cast<std::size_t>(p)] = 1;
    seen[static_cast<std::size_t>(c)] = 1;
    branch_length_by_node[static_cast<std::size_t>(c)] =
      lengths[static_cast<std::size_t>(i)];
    edge_index_by_node[static_cast<std::size_t>(c)] = i + 1;
  }

  // Standard phylo node ids are contiguous.  Checking every id catches a
  // missing internal id that would otherwise be indistinguishable from an
  // isolated node after max_node is inferred.
  for (int node = 1; node <= n_total; ++node) {
    if (!seen[static_cast<std::size_t>(node)]) {
      Rcpp::stop(
        "edge node numbers must be contiguous from 1 through the largest node id."
      );
    }
  }

  int root = 0;
  int root_count = 0;
  for (int node = 1; node <= n_total; ++node) {
    const int in = indegree[static_cast<std::size_t>(node)];
    const int out = outdegree[static_cast<std::size_t>(node)];
    if (node <= n_tip) {
      if (in != 1 || out != 0) {
        Rcpp::stop(
          "tip nodes (1:n_tip) must occur exactly once as children and never as parents."
        );
      }
    } else {
      if (in == 0) {
        root = node;
        root_count += 1;
      } else if (in != 1) {
        Rcpp::stop("internal nodes must have exactly one parent except the root.");
      }
      if (out == 0) {
        Rcpp::stop("internal nodes must have at least one child.");
      }
    }
  }
  if (root_count != 1 || root <= n_tip) {
    Rcpp::stop("tree topology must contain exactly one internal root.");
  }

  // Validate connectedness and produce traversal order in one iterative DFS.
  // Input edge order is retained within each parent's child list, preserving
  // the deterministic ordering supplied by ape while supporting polytomies.
  std::vector<int> preorder;
  std::vector<int> postorder;
  preorder.reserve(static_cast<std::size_t>(n_total));
  postorder.reserve(static_cast<std::size_t>(n_total));
  std::vector<char> visited(static_cast<std::size_t>(n_total + 1), 0);
  std::vector<double> root_distance(
    static_cast<std::size_t>(n_total + 1), 0.0
  );

  struct Frame {
    int node;
    std::size_t next_child;
  };
  std::vector<Frame> stack;
  stack.reserve(static_cast<std::size_t>(n_total));
  visited[static_cast<std::size_t>(root)] = 1;
  preorder.push_back(root);
  stack.push_back(Frame{root, 0});

  while (!stack.empty()) {
    Frame& frame = stack.back();
    const std::vector<int>& child_list =
      children[static_cast<std::size_t>(frame.node)];
    if (frame.next_child < child_list.size()) {
      const int child = child_list[frame.next_child++];
      if (visited[static_cast<std::size_t>(child)]) {
        Rcpp::stop("edge topology contains a cycle or repeated node.");
      }
      visited[static_cast<std::size_t>(child)] = 1;
      const double distance =
        root_distance[static_cast<std::size_t>(frame.node)] +
        branch_length_by_node[static_cast<std::size_t>(child)];
      if (!std::isfinite(distance)) {
        Rcpp::stop("root distances overflow; branch lengths are too large.");
      }
      root_distance[static_cast<std::size_t>(child)] = distance;
      preorder.push_back(child);
      stack.push_back(Frame{child, 0});
    } else {
      postorder.push_back(frame.node);
      stack.pop_back();
    }
  }

  for (int node = 1; node <= n_total; ++node) {
    if (!visited[static_cast<std::size_t>(node)]) {
      Rcpp::stop("edge topology is disconnected from its root.");
    }
  }
  if (static_cast<int>(preorder.size()) != n_total ||
      static_cast<int>(postorder.size()) != n_total) {
    Rcpp::stop("edge topology does not form a rooted tree.");
  }

  std::vector<int> tip_order;
  tip_order.reserve(static_cast<std::size_t>(n_tip));
  for (std::size_t i = 0; i < preorder.size(); ++i) {
    const int node = preorder[i];
    if (node <= n_tip) {
      tip_order.push_back(node);
    }
  }
  if (static_cast<int>(tip_order.size()) != n_tip) {
    Rcpp::stop("tree topology does not reach every tip.");
  }

  ParsedTree out;
  out.n_tip = n_tip;
  out.n_total = n_total;
  out.n_node = n_node;
  out.root = root;
  out.parent.assign(parent.begin() + 1, parent.end());
  out.children.swap(children);
  out.edge_length = lengths;
  out.branch_length_by_node.assign(
    branch_length_by_node.begin() + 1, branch_length_by_node.end()
  );
  out.edge_index_by_node.assign(
    edge_index_by_node.begin() + 1, edge_index_by_node.end()
  );
  out.preorder = preorder;
  out.postorder = postorder;
  out.root_distance.assign(root_distance.begin() + 1, root_distance.end());
  out.tip_order = tip_order;
  return out;
}

Rcpp::IntegerVector as_integer_vector(const std::vector<int>& values) {
  Rcpp::IntegerVector out(static_cast<R_xlen_t>(values.size()));
  std::copy(values.begin(), values.end(), out.begin());
  return out;
}

Rcpp::NumericVector as_numeric_vector(const std::vector<double>& values) {
  Rcpp::NumericVector out(static_cast<R_xlen_t>(values.size()));
  std::copy(values.begin(), values.end(), out.begin());
  return out;
}

Rcpp::List compiled_list(const ParsedTree& tree) {
  const int n_total = tree.n_total;
  const int n_edge = n_total - 1;

  // CSR offsets are returned with an R-friendly one-based base: children for
  // node p occupy children[child_ptr[p] : child_ptr[p + 1] - 1].  Node ids and
  // traversal arrays are likewise one-based; parent[root] is the 0 sentinel.
  Rcpp::IntegerVector child_ptr(static_cast<R_xlen_t>(n_total + 1));
  Rcpp::IntegerVector children(static_cast<R_xlen_t>(n_edge));
  int cursor = 0;
  child_ptr[0] = 1;
  for (int node = 1; node <= n_total; ++node) {
    const std::vector<int>& child_list =
      tree.children[static_cast<std::size_t>(node)];
    for (std::size_t j = 0; j < child_list.size(); ++j) {
      children[cursor++] = child_list[j];
    }
    child_ptr[node] = cursor + 1;
  }

  Rcpp::IntegerVector tip_nodes(tree.n_tip);
  Rcpp::IntegerVector tip_index(n_total);
  for (int tip = 1; tip <= tree.n_tip; ++tip) {
    tip_nodes[tip - 1] = tip;
    tip_index[tip - 1] = tip;
  }
  for (int node = tree.n_tip + 1; node <= n_total; ++node) {
    tip_index[node - 1] = 0;
  }

  const Rcpp::IntegerVector parent = as_integer_vector(tree.parent);
  const Rcpp::NumericVector branch_length_by_node =
    as_numeric_vector(tree.branch_length_by_node);
  const Rcpp::IntegerVector edge_index =
    as_integer_vector(tree.edge_index_by_node);
  const Rcpp::IntegerVector preorder = as_integer_vector(tree.preorder);
  const Rcpp::IntegerVector postorder = as_integer_vector(tree.postorder);
  const Rcpp::NumericVector root_distance =
    as_numeric_vector(tree.root_distance);
  const Rcpp::NumericVector edge_length = as_numeric_vector(tree.edge_length);
  const Rcpp::IntegerVector tip_order = as_integer_vector(tree.tip_order);

  // This is an explicit lower-bound estimate for the payload arrays, useful
  // for cache diagnostics.  It intentionally excludes R object headers.
  const std::uint64_t bytes_int = sizeof(int);
  const std::uint64_t bytes_double = sizeof(double);
  const std::uint64_t payload_bytes =
    bytes_int * static_cast<std::uint64_t>(n_total) +       // parent
    bytes_int * static_cast<std::uint64_t>(n_total + 1) +   // child_ptr
    bytes_int * static_cast<std::uint64_t>(n_edge) +        // children
    bytes_int * static_cast<std::uint64_t>(n_total) +       // edge_index
    bytes_int * static_cast<std::uint64_t>(2 * n_total) +   // traversals
    bytes_int * static_cast<std::uint64_t>(tree.n_tip) +    // tip maps
    bytes_double * static_cast<std::uint64_t>(n_total) +    // branch lengths
    bytes_double * static_cast<std::uint64_t>(n_total);     // root distances

  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["n_tip"] = tree.n_tip,
    Rcpp::_["n_node"] = tree.n_node,
    Rcpp::_["n_total"] = tree.n_total,
    Rcpp::_["n_edge"] = n_edge,
    Rcpp::_["root"] = tree.root,
    Rcpp::_["parent"] = parent,
    Rcpp::_["child_ptr"] = child_ptr,
    Rcpp::_["children"] = children,
    Rcpp::_["branch_length_by_node"] = branch_length_by_node,
    Rcpp::_["edge_length"] = edge_length,
    Rcpp::_["edge_index"] = edge_index,
    Rcpp::_["preorder"] = preorder,
    Rcpp::_["postorder"] = postorder,
    Rcpp::_["root_distance"] = root_distance,
    Rcpp::_["tip_nodes"] = tip_nodes,
    Rcpp::_["tip_index"] = tip_index,
    Rcpp::_["tip_order"] = tip_order,
    Rcpp::_["topology_bytes"] = static_cast<double>(payload_bytes)
  );

  // Descriptive aliases make the list convenient for both C++ callers and R
  // code while preserving one canonical storage vector for each quantity.
  out["children_ptr"] = child_ptr;
  out["children_offset"] = child_ptr;
  out["branch_length"] = branch_length_by_node;
  out["node_to_tip"] = tip_index;
  out.attr("class") = Rcpp::CharacterVector::create(
    "fastphylosig_compiled_tree", "list"
  );
  out.attr("index_base") = 1;
  out.attr("root_branch_length") = 0.0;
  return out;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List compile_tree_cpp(SEXP edge, SEXP edge_length, SEXP n_tip) {
  return compiled_list(parse_and_compile(edge, edge_length, n_tip));
}

// [[Rcpp::export]]
Rcpp::List tree_core_validate_cpp(SEXP edge, SEXP edge_length, SEXP n_tip) {
  try {
    const ParsedTree tree = parse_and_compile(edge, edge_length, n_tip);
    return Rcpp::List::create(
      Rcpp::_["valid"] = true,
      Rcpp::_["message"] = "",
      Rcpp::_["n_tip"] = tree.n_tip,
      Rcpp::_["n_node"] = tree.n_node,
      Rcpp::_["n_total"] = tree.n_total,
      Rcpp::_["root"] = tree.root
    );
  } catch (const std::exception& ex) {
    return Rcpp::List::create(
      Rcpp::_["valid"] = false,
      Rcpp::_["message"] = std::string(ex.what())
    );
  }
}
