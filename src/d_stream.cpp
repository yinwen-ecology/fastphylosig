// Streaming simulation kernel for the binary-tree Fritz--Purvis D statistic.
//
// This file deliberately keeps the simulation and contrast traversal in one
// place.  Unlike the legacy oracle path, it never constructs an n_tip x nsim
// matrix: a tip vector is generated, thresholded (for the Brownian null), and
// consumed by the scalar contrast traversal before the next simulation is
// started.  The chunk loop bounds the amount of temporary work and is exposed
// as a compatibility/performance knob.

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "numeric_utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

int effective_threads(const int n_threads) {
#ifdef _OPENMP
  return n_threads > 1 ? n_threads : 1;
#else
  (void) n_threads;
  return 1;
#endif
}

struct BinaryTree {
  int n_tip;
  int n_edge;
  int max_node;
  int root;
  std::vector<int> child0;
  std::vector<int> child1;
  std::vector<int> edge0;
  std::vector<int> edge1;
  std::vector<int> child_to_edge;
  std::vector< std::vector<int> > children;
  std::vector< std::vector<double> > branch_lengths;
  std::vector<int> parent_order;
};

BinaryTree validate_binary_tree(const Rcpp::IntegerMatrix& edge,
                                const arma::vec& edge_length,
                                const int n_tip) {
  if (n_tip < 2) {
    Rcpp::stop("n_tip must be at least two for phylo.d.");
  }
  if (edge.ncol() != 2 ||
      edge_length.n_elem != static_cast<arma::uword>(edge.nrow())) {
    Rcpp::stop("edge must be a two-column matrix matching edge_length.");
  }

  BinaryTree tree;
  tree.n_tip = n_tip;
  tree.n_edge = edge.nrow();
  tree.max_node = n_tip;
  for (int i = 0; i < tree.n_edge; ++i) {
    tree.max_node = std::max(tree.max_node, edge(i, 0));
    tree.max_node = std::max(tree.max_node, edge(i, 1));
    if (edge_length(i) < 0.0 || !std::isfinite(edge_length(i))) {
      Rcpp::stop("edge_length must be finite and non-negative.");
    }
    if (edge(i, 0) < 1 || edge(i, 1) < 1) {
      Rcpp::stop("edge contains node numbers outside the expected range.");
    }
  }

  tree.root = n_tip + 1;
  if (tree.root > tree.max_node) {
    Rcpp::stop("Could not identify the root node.");
  }

  tree.child0.assign(tree.max_node + 1, -1);
  tree.child1.assign(tree.max_node + 1, -1);
  tree.edge0.assign(tree.max_node + 1, -1);
  tree.edge1.assign(tree.max_node + 1, -1);
  tree.child_to_edge.assign(tree.max_node + 1, -1);
  tree.children.assign(tree.max_node + 1, std::vector<int>());
  tree.branch_lengths.assign(tree.max_node + 1, std::vector<double>());

  std::vector<int> parent_count(tree.max_node + 1, 0);
  std::vector<int> is_child(tree.max_node + 1, 0);
  for (int i = 0; i < tree.n_edge; ++i) {
    const int parent = edge(i, 0);
    const int child = edge(i, 1);
    if (parent > tree.max_node || child > tree.max_node) {
      Rcpp::stop("edge contains node numbers outside the expected range.");
    }
    ++parent_count[parent];
    is_child[child] = 1;
    tree.children[parent].push_back(child);
    tree.branch_lengths[parent].push_back(edge_length(i));
    tree.child_to_edge[child] = i;
  }

  // A normal ape binary tree has exactly two children per internal node.  The
  // R wrapper routes polytomies to the oracle path; this guard also prevents a
  // direct call from silently changing that contract.
  for (int node = n_tip + 1; node <= tree.max_node; ++node) {
    if (parent_count[node] != 2) {
      Rcpp::stop("phylo_d_stream_cpp requires an ordinary binary tree.");
    }
  }
  if (is_child[tree.root]) {
    Rcpp::stop("Could not identify the root node.");
  }

  // Preserve the row/group traversal used by phylo_d_sum_binary_batch().
  int row = 0;
  while (row < tree.n_edge) {
    const int parent = edge(row, 0);
    tree.parent_order.push_back(parent);
    int n_child = 0;
    while (row < tree.n_edge && edge(row, 0) == parent) {
      if (n_child == 0) {
        tree.child0[parent] = edge(row, 1);
        tree.edge0[parent] = row;
      } else if (n_child == 1) {
        tree.child1[parent] = edge(row, 1);
        tree.edge1[parent] = row;
      }
      ++n_child;
      ++row;
    }
    if (n_child != 2 || tree.child0[parent] < 1 || tree.child1[parent] < 1) {
      Rcpp::stop("phylo_d_stream_cpp requires an ordinary binary tree.");
    }
  }

  return tree;
}

double contrast_sum(const std::vector<double>& tip_values,
                    const BinaryTree& tree,
                    const arma::vec& edge_length) {
  if (static_cast<int>(tip_values.size()) != tree.n_tip) {
    Rcpp::stop("tip state vector has an incompatible length.");
  }

  // Use one scalar per node and edge.  This is the same arithmetic and edge
  // adjustment order as phylo_d_sum_binary_batch(), but without a simulation
  // column dimension.
  std::vector<double> node_value(tree.max_node + 1, 0.0);
  for (int tip = 1; tip <= tree.n_tip; ++tip) {
    node_value[tip] = tip_values[static_cast<std::size_t>(tip - 1)];
  }
  std::vector<double> edge_len(static_cast<std::size_t>(tree.n_edge));
  for (int i = 0; i < tree.n_edge; ++i) {
    edge_len[static_cast<std::size_t>(i)] = edge_length(i);
  }

  double sum = 0.0;
  for (std::size_t parent_index = 0;
       parent_index < tree.parent_order.size(); ++parent_index) {
    const int parent_id = tree.parent_order[parent_index];
    const int child0 = tree.child0[parent_id];
    const int child1 = tree.child1[parent_id];
    const int edge0 = tree.edge0[parent_id];
    const int edge1 = tree.edge1[parent_id];

    const double bl0 = edge_len[static_cast<std::size_t>(edge0)];
    const double bl1 = edge_len[static_cast<std::size_t>(edge1)];
    if (bl0 <= 0.0 || bl1 <= 0.0 ||
        !std::isfinite(bl0) || !std::isfinite(bl1)) {
      Rcpp::stop("phylo_d requires strictly positive branch lengths.");
    }
    const double inv0 = 1.0 / bl0;
    const double inv1 = 1.0 / bl1;
    const double denom = inv0 + inv1;
    const double curr_nv =
      (node_value[child0] * inv0 + node_value[child1] * inv1) / denom;
    const double curr_contr =
      std::abs(node_value[child0] - curr_nv) +
      std::abs(node_value[child1] - curr_nv);
    node_value[parent_id] = curr_nv;
    sum += curr_contr;

    const double curr_bl_adj = (bl0 * bl1) / (bl0 + bl1);
    if (parent_id != tree.root) {
      const int parent_edge = tree.child_to_edge[parent_id];
      if (parent_edge >= 0) {
        edge_len[static_cast<std::size_t>(parent_edge)] += curr_bl_adj;
      }
    }
  }

  return sum;
}

arma::vec contrast_sum_batch(const arma::mat& states,
                            const BinaryTree& tree,
                            const arma::vec& edge_length,
                            const int n_threads) {
  const int n_col = static_cast<int>(states.n_cols);
  if (states.n_rows != static_cast<arma::uword>(tree.n_tip)) {
    Rcpp::stop("states must have one row per tip.");
  }
  for (int i = 0; i < tree.n_edge; ++i) {
    if (edge_length(i) <= 0.0 || !std::isfinite(edge_length(i))) {
      Rcpp::stop("phylo_d requires strictly positive branch lengths.");
    }
  }
  arma::mat node_value(tree.max_node, n_col, arma::fill::zeros);
  node_value.rows(0, tree.n_tip - 1) = states;
  arma::mat edge_len(tree.n_edge, n_col, arma::fill::zeros);
  for (int i = 0; i < tree.n_edge; ++i) {
    edge_len.row(i).fill(edge_length(i));
  }
  arma::vec sums(n_col, arma::fill::zeros);
  const int n_threads_eff = effective_threads(n_threads);

  // Each simulation column is independent once tip states are available.  A
  // bounded chunk therefore permits OpenMP traversal without sharing RNG
  // state, while preserving the legacy column-wise arithmetic order.
  for (std::size_t parent_index = 0;
       parent_index < tree.parent_order.size(); ++parent_index) {
    const int parent_id = tree.parent_order[parent_index];
    const int child0 = tree.child0[parent_id];
    const int child1 = tree.child1[parent_id];
    const int edge0 = tree.edge0[parent_id];
    const int edge1 = tree.edge1[parent_id];

#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads_eff) \
  if (n_threads_eff > 1 && n_col > 1) schedule(static)
#endif
    for (int c = 0; c < n_col; ++c) {
      const double bl0 = edge_len(edge0, c);
      const double bl1 = edge_len(edge1, c);
      const double inv0 = 1.0 / bl0;
      const double inv1 = 1.0 / bl1;
      const double denom = inv0 + inv1;
      const double curr_nv =
        (node_value(child0 - 1, c) * inv0 +
         node_value(child1 - 1, c) * inv1) / denom;
      const double curr_contr =
        std::abs(node_value(child0 - 1, c) - curr_nv) +
        std::abs(node_value(child1 - 1, c) - curr_nv);
      node_value(parent_id - 1, c) = curr_nv;
      sums(c) += curr_contr;
      const double curr_bl_adj = (bl0 * bl1) / (bl0 + bl1);
      if (parent_id != tree.root) {
        const int parent_edge = tree.child_to_edge[parent_id];
        if (parent_edge >= 0) edge_len(parent_edge, c) += curr_bl_adj;
      }
    }
  }
  return sums;
}

void brownian_tips(const BinaryTree& tree,
                   const arma::vec& edge_length,
                   std::vector<double>& node_value,
                   std::vector<double>& tip_values) {
  std::fill(node_value.begin(), node_value.end(), 0.0);
  std::vector<int> stack;
  stack.reserve(static_cast<std::size_t>(tree.max_node));
  stack.push_back(tree.root);

  while (!stack.empty()) {
    const int parent = stack.back();
    stack.pop_back();
    const std::vector<int>& children = tree.children[parent];
    const std::vector<double>& lengths = tree.branch_lengths[parent];
    for (std::size_t k = 0; k < children.size(); ++k) {
      const int child = children[k];
      const double bl = lengths[k];
      // Keep the RNG draw for zero branches, as in brownian_tree_threshold_cpp.
      node_value[child] = node_value[parent] + std::sqrt(bl) * R::rnorm(0.0, 1.0);
      if (!tree.children[child].empty()) stack.push_back(child);
    }
  }

  for (int tip = 1; tip <= tree.n_tip; ++tip) {
    tip_values[static_cast<std::size_t>(tip - 1)] = node_value[tip];
  }
}

} // namespace

// [[Rcpp::export]]
Rcpp::List phylo_d_stream_cpp(const arma::vec& observed,
                              const Rcpp::IntegerMatrix& edge,
                              const arma::vec& edge_length,
                              const int n_tip,
                              const int nsim,
                              const double prop_state1,
                              const int chunk_size = 128,
                              const bool return_sim = true,
                              const int n_threads = 1) {
  (void) n_threads; // RNG draws are serial; traversal itself is O(n_tip).
  if (nsim < 1) {
    Rcpp::stop("nsim must be positive.");
  }
  if (prop_state1 < 0.0 || prop_state1 > 1.0 ||
      !std::isfinite(prop_state1)) {
    Rcpp::stop("prop_state1 must be in [0, 1].");
  }
  if (chunk_size < 32 || chunk_size > 512) {
    Rcpp::stop("chunk_size must be an integer between 32 and 512.");
  }
  if (observed.n_elem != static_cast<arma::uword>(n_tip)) {
    Rcpp::stop("observed must have one value per tip.");
  }
  if (!observed.is_finite()) {
    Rcpp::stop("observed must contain only finite values.");
  }

  const BinaryTree tree = validate_binary_tree(edge, edge_length, n_tip);
  std::vector<double> observed_values(static_cast<std::size_t>(n_tip));
  for (int tip = 0; tip < n_tip; ++tip) {
    observed_values[static_cast<std::size_t>(tip)] = observed(tip);
  }
  const double observed_sum = contrast_sum(observed_values, tree, edge_length);

  const bool keep_null = return_sim;
  Rcpp::NumericVector random_sums(keep_null ? nsim : 0);
  Rcpp::NumericVector brownian_sums(keep_null ? nsim : 0);
  std::vector<double> node_value(static_cast<std::size_t>(tree.max_node + 1));

  double random_total = 0.0;
  double brownian_total = 0.0;
  double random_less = 0.0;
  double brownian_greater = 0.0;

  // Simulations are processed in bounded chunks.  Each state vector is
  // immediately traversed and discarded, so no n_tip x nsim state matrix is
  // ever allocated.  A chunk matrix is intentionally bounded by chunk_size;
  // this is the unit of OpenMP traversal and the only simulation-state
  // storage retained at any point.
  for (int first = 0; first < nsim; first += chunk_size) {
    const int last = std::min(nsim, first + chunk_size);
    const int n_chunk = last - first;
    arma::mat random_states(n_tip, n_chunk, arma::fill::zeros);
    arma::mat brownian_states(n_tip, n_chunk, arma::fill::zeros);
    std::vector<int> permutation(static_cast<std::size_t>(n_tip));

    // Generate unweighted random-association states (one permutation per
    // column).  This preserves observed prevalence exactly, as sample(ds)
    // does in the legacy R path.
    for (int c = 0; c < n_chunk; ++c) {
      // sample(ds) in the legacy R path is an unweighted permutation (not
      // independent sampling with replacement).  Fisher--Yates preserves the
      // observed state prevalence in every random-association replicate.
      for (int tip = 0; tip < n_tip; ++tip) {
        permutation[static_cast<std::size_t>(tip)] = tip;
      }
      for (int tip = n_tip - 1; tip > 0; --tip) {
        const double u = R::runif(0.0, static_cast<double>(tip + 1));
        int pick = static_cast<int>(std::floor(u));
        if (pick < 0) pick = 0;
        if (pick > tip) pick = tip;
        std::swap(permutation[static_cast<std::size_t>(tip)],
                  permutation[static_cast<std::size_t>(pick)]);
      }
      for (int tip = 0; tip < n_tip; ++tip) {
        random_states(tip, c) =
          observed(permutation[static_cast<std::size_t>(tip)]);
      }
    }

    // Brownian simulation, type-7 thresholding, and state construction are
    // fused per column.  The raw matrix exists only for this bounded chunk so
    // each column can be thresholded before its contrast traversal.
    std::vector<double> sorted_values(static_cast<std::size_t>(n_tip));
    std::vector<double> brownian_tip_values(static_cast<std::size_t>(n_tip));
    for (int c = 0; c < n_chunk; ++c) {
      brownian_tips(tree, edge_length, node_value, brownian_tip_values);
      for (int tip = 0; tip < n_tip; ++tip) {
        sorted_values[static_cast<std::size_t>(tip)] =
          brownian_tip_values[static_cast<std::size_t>(tip)];
      }
      std::sort(sorted_values.begin(), sorted_values.end());
      const double threshold =
        fastphylosig::quantile_type7_sorted(sorted_values, prop_state1);
      for (int tip = 0; tip < n_tip; ++tip) {
        brownian_states(tip, c) =
          brownian_tip_values[static_cast<std::size_t>(tip)] < threshold ?
          1.0 : 0.0;
      }
    }

    const arma::vec random_chunk = contrast_sum_batch(
      random_states, tree, edge_length, n_threads
    );
    const arma::vec brownian_chunk = contrast_sum_batch(
      brownian_states, tree, edge_length, n_threads
    );
    for (int c = 0; c < n_chunk; ++c) {
      const int sim = first + c;
      const double random_sum = random_chunk(c);
      const double brownian_sum = brownian_chunk(c);
      random_total += random_sum;
      brownian_total += brownian_sum;
      if (random_sum < observed_sum) random_less += 1.0;
      if (brownian_sum > observed_sum) brownian_greater += 1.0;
      if (keep_null) {
        random_sums[sim] = random_sum;
        brownian_sums[sim] = brownian_sum;
      }
    }
  }

  const double mean_random = random_total / static_cast<double>(nsim);
  const double mean_brownian = brownian_total / static_cast<double>(nsim);
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("observed") = observed_sum,
    Rcpp::Named("mean_random") = mean_random,
    Rcpp::Named("mean_brownian") = mean_brownian,
    Rcpp::Named("p_random") = random_less / static_cast<double>(nsim),
    Rcpp::Named("p_brownian") = brownian_greater / static_cast<double>(nsim)
  );
  if (keep_null) {
    out["random"] = random_sums;
    out["brownian"] = brownian_sums;
  }
  return out;
}
