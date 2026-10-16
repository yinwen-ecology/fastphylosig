#ifndef FASTPHYLOSIG_BOUNDED_OPTIMIZE_H
#define FASTPHYLOSIG_BOUNDED_OPTIMIZE_H
#include <functional>

// Implemented in lambda_tree.cpp; shared deterministic Brent search.
double brent_minimize(double lower, double upper, double tol,
                      const std::function<double(double)>& objective);
#endif
