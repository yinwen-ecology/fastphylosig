#ifndef FASTPHYLOSIG_CONTINUOUS_KERNELS_H
#define FASTPHYLOSIG_CONTINUOUS_KERNELS_H
#include <Rcpp.h>

// Shared production kernels; callers inside C++ bypass generated .Call wrappers.
Rcpp::List fast_k_tree_batch_cpp(const Rcpp::List& compiled_tree,
                               const Rcpp::NumericMatrix& X, int trait_chunk);
Rcpp::List fast_k_tree_permutation_cpp(const Rcpp::List& compiled_tree,
    const Rcpp::NumericMatrix& X, int nsim, SEXP permutations, int trait_chunk,
    bool return_sim, bool include_observed, int n_threads, int simulation_chunk);
Rcpp::List fast_lambda_tree_optimize_cpp(const Rcpp::List& compiled_tree,
    const Rcpp::NumericMatrix& X, SEXP max_lambda, bool test, bool profile,
    int profile_points, int trait_chunk, int n_threads, double tol);
Rcpp::NumericVector heterogeneity_bm_k_cpp(const Rcpp::List& compiled_tree,
                                         int nsim, int simulation_chunk);
#endif
