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
#include "d_kernels.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

struct BinaryTree {
  int n_tip;
  int n_edge;
  int max_node;
  int root;
  std::vector< std::vector<int> > children;
  std::vector< std::vector<double> > branch_lengths;
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
    int n_child = 0;
    while (row < tree.n_edge && edge(row, 0) == parent) {
      ++n_child;
      ++row;
    }
    if (n_child != 2) {
      Rcpp::stop("phylo_d_stream_cpp requires an ordinary binary tree.");
    }
  }

  return tree;
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
  const fastphylosig::DContrasts contrasts =
    fastphylosig::DContrasts::from_edges(edge, edge_length.memptr(), n_tip);
  std::vector<double> contrast_work(contrasts.n_total());
  const double observed_sum = contrasts.sum(observed_values.data(), contrast_work.data());

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

    const std::vector<double> random_chunk = contrasts.batch(random_states.memptr(), n_chunk, n_threads);
    const std::vector<double> brownian_chunk = contrasts.batch(brownian_states.memptr(), n_chunk, n_threads);
    for (int c = 0; c < n_chunk; ++c) {
      const int sim = first + c;
      const double random_sum = random_chunk[c];
      const double brownian_sum = brownian_chunk[c];
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
