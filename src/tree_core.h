#ifndef FASTPHYLOSIG_TREE_CORE_H
#define FASTPHYLOSIG_TREE_CORE_H

// This header deliberately contains only the small public declarations used
// by tree_core.cpp.  The compiled representation returned to R is a regular
// list (rather than an external pointer), so it is safe to copy and serialize
// in a prepared-tree context.

#include <Rcpp.h>

Rcpp::List compile_tree_cpp(SEXP edge, SEXP edge_length, SEXP n_tip);

// Return a non-throwing validation diagnostic.  Invalid input is represented
// by valid = FALSE and a human-readable message; compile_tree_cpp() itself
// throws for invalid input.
Rcpp::List tree_core_validate_cpp(SEXP edge, SEXP edge_length, SEXP n_tip);

#endif  // FASTPHYLOSIG_TREE_CORE_H
