#ifndef FASTPHYLOSIG_K_SCALE_H
#define FASTPHYLOSIG_K_SCALE_H
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace fastphylosig {
// Use exact binary scaling only outside the ordinary arithmetic range. This
// must work when long double has the same exponent range as double (e.g. ARM
// macOS). Scaling both axes prevents squared energies and currents overflowing.
inline int k_extreme_exponent(const double magnitude) {
  if (!(magnitude > 0) || !std::isfinite(magnitude)) return 0;
  const int exponent = std::ilogb(magnitude);
  return exponent > 128 || exponent < -128 ? exponent : 0;
}
template<class Tree>
int normalize_k_tree(Tree& tree) {
  double largest = 0;
  for (double value : tree.branch) largest = std::max(largest, std::fabs(value));
  const int exponent = k_extreme_exponent(largest);
  if (exponent == 0) return 0;
  // Do not silently erase a genuine positive edge in an exceptionally wide
  // dynamic-range tree. Such a tree retains the original checked path.
  for (double value : tree.branch)
    if (value > 0 && !(std::scalbn(value, -exponent) > 0)) return 0;
  for (double& value : tree.branch) value = std::scalbn(value, -exponent);
  for (double& value : tree.root_distance) value = std::scalbn(value, -exponent);
  return exponent;
}
inline Rcpp::NumericMatrix normalize_k_traits(const Rcpp::NumericMatrix& input,
                                             std::vector<int>& exponents) {
  exponents.assign(input.ncol(), 0);
  bool changed = false;
  for (int j = 0; j < input.ncol(); ++j) {
    double largest = 0;
    for (int i = 0; i < input.nrow(); ++i) {
      const double value = input(i,j);
      if (!std::isfinite(value)) Rcpp::stop("X must contain only finite trait values.");
      largest = std::max(largest, std::fabs(value));
    }
    exponents[j] = k_extreme_exponent(largest);
    changed = changed || exponents[j] != 0;
  }
  if (!changed) return input;
  Rcpp::NumericMatrix output = Rcpp::clone(input);
  for (int j = 0; j < input.ncol(); ++j) if (exponents[j] != 0)
    for (int i = 0; i < input.nrow(); ++i)
      output(i,j) = std::scalbn(input(i,j), -exponents[j]);
  return output;
}
}
#endif
