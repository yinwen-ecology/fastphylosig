#ifndef FASTPHYLOSIG_NUMERIC_UTILS_H
#define FASTPHYLOSIG_NUMERIC_UTILS_H

#include <cmath>
#include <limits>
#include <vector>

namespace fastphylosig {

inline double quantile_type7_sorted(const std::vector<double>& values,
                                    const double probability) {
  const int n = static_cast<int>(values.size());
  if (n < 1) return std::numeric_limits<double>::quiet_NaN();

  const double h = 1.0 + static_cast<double>(n - 1) * probability;
  const int lo = static_cast<int>(std::floor(h));
  const double gamma = h - static_cast<double>(lo);
  if (lo <= 1) return values.front();
  if (lo >= n) return values.back();
  return (1.0 - gamma) * values[static_cast<std::size_t>(lo - 1)] +
    gamma * values[static_cast<std::size_t>(lo)];
}

} // namespace fastphylosig

#endif
