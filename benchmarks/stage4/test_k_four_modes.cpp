// Test-only bridge into the production K evaluator.  Including the production
// translation unit keeps these calls source-identical to k_permutation.cpp;
// no force-mode entry point is added to the installed package API.
//
// Load with Rcpp::sourceCpp() from the repository root.  The inclusion path is
// relative to this file: benchmarks/stage4/../../src/k_permutation.cpp.
// [[Rcpp::plugins(cpp17)]]
// Hide the path behind a macro so sourceCpp does not also compile the included
// .cpp as a separate translation unit (the test bridge needs one TU).
#define FASTPHYLOSIG_STAGE4_PRODUCTION_SOURCE "../../src/k_permutation.cpp"
#include FASTPHYLOSIG_STAGE4_PRODUCTION_SOURCE
#undef FASTPHYLOSIG_STAGE4_PRODUCTION_SOURCE

// [[Rcpp::export]]
Rcpp::List stage4_test_k_four_modes(
    const Rcpp::List& compiled_tree,
    const Rcpp::NumericMatrix& X,
    const Rcpp::IntegerMatrix& permutations) {
  using namespace kperm;

  const Tree tree = parse_tree(compiled_tree);
  const Cache cache = build_cache(tree);
  const int n = tree.n_tip;
  const int p = X.ncol();
  const int nsim = permutations.nrow();
  if (X.nrow() != n || p < 1) {
    Rcpp::stop("X must have one row per compiled tree tip and at least one trait.");
  }
  if (permutations.ncol() != n || nsim < 1) {
    Rcpp::stop("permutations must have nsim rows and one column per tip.");
  }

  // The test caller supplies valid 1-based permutations.  Copy them once into
  // the production evaluator's row-major, zero-based layout.
  std::vector<int> perm_flat(static_cast<std::size_t>(nsim) * n);
  for (int i = 0; i < nsim; ++i) {
    std::vector<unsigned char> seen(static_cast<std::size_t>(n), 0);
    for (int r = 0; r < n; ++r) {
      const int value = permutations(i, r);
      if (value < 1 || value > n || seen[static_cast<std::size_t>(value - 1)]) {
        Rcpp::stop("each permutation row must contain 1:n exactly once.");
      }
      seen[static_cast<std::size_t>(value - 1)] = 1;
      perm_flat[static_cast<std::size_t>(i) * n + r] = value - 1;
    }
  }

  const double* x = REAL(X);
  const std::vector<double> q = kfused::build_q(tree, cache);

  const auto evaluate = [&](const bool hybrid, const int block_cap) {
    std::vector<int> identity(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) identity[static_cast<std::size_t>(i)] = i;

    std::vector<double> observed;
    if (hybrid) {
      kfused::compute_one_hybrid(tree, cache, q, x, p, identity, p, observed);
    } else {
      compute_one(tree, cache, x, p, identity, p, observed);
    }

    Rcpp::NumericMatrix null(nsim, p);
    std::vector<double> message, state, baseline, delta, block_out;
    for (int first = 0; first < nsim; first += block_cap) {
      const int count = std::min(block_cap, nsim - first);
      kfused::compute_block(
        tree, cache, q, x, p, 0, p,
        perm_flat.data() + static_cast<std::size_t>(first) * n,
        count, block_cap, hybrid,
        message, state, baseline, delta, block_out
      );
      for (int i = 0; i < count; ++i) {
        for (int j = 0; j < p; ++j) {
          null(first + i, j) = block_out[static_cast<std::size_t>(i) * p + j];
        }
      }
    }

    Rcpp::IntegerVector exceedance(p);
    Rcpp::NumericVector p_value(p);
    Rcpp::LogicalMatrix tail(nsim, p);
    for (int j = 0; j < p; ++j) {
      int hits = 0;
      for (int i = 0; i < nsim; ++i) {
        const bool is_upper_tail =
          fastphylosig::inclusive_upper_tail(null(i, j), observed[j]);
        tail(i, j) = is_upper_tail;
        if (is_upper_tail) ++hits;
      }
      exceedance[j] = hits;
      p_value[j] = static_cast<double>(hits) / static_cast<double>(nsim);
    }

    return Rcpp::List::create(
      Rcpp::Named("K_obs") = Rcpp::wrap(observed),
      Rcpp::Named("K_null") = null,
      Rcpp::Named("P") = p_value,
      Rcpp::Named("exceedance") = exceedance,
      Rcpp::Named("tail_decision") = tail
    );
  };

  return Rcpp::List::create(
    Rcpp::Named("legacy_block1") = evaluate(false, 1),
    Rcpp::Named("hybrid_block1") = evaluate(true, 1),
    Rcpp::Named("legacy_block4") = evaluate(false, 4),
    Rcpp::Named("hybrid_block4") = evaluate(true, 4),
    Rcpp::Named("permutations") = permutations
  );
}
