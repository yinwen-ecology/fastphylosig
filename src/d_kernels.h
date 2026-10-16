#ifndef FASTPHYLOSIG_D_KERNELS_H
#define FASTPHYLOSIG_D_KERNELS_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace fastphylosig {

// Trait-independent PIC factors shared by fast_d's scalar/batch paths and
// heterogeneity. Retain each path's established floating-point convention:
// reciprocal arithmetic for public D; scaled weights/extended sum for hetero.
class DContrasts {
 public:
  enum Arithmetic { reciprocal, scaled };
  DContrasts() : nt_(0), nn_(0), arithmetic_(reciprocal) {}
  DContrasts(int nt, int nn, const std::vector<int>& post,
      const std::vector<int>& left, const std::vector<int>& right,
      const std::vector<double>& branch, Arithmetic arithmetic = reciprocal)
      : nt_(nt), nn_(nn), arithmetic_(arithmetic) {
    std::vector<double> variance = branch;
    for (int node : post) if (node >= nt) {
      const int a = left[node], b = right[node];
      nodes_.push_back(node); left_.push_back(a); right_.push_back(b);
      if (b < 0) { w0_.push_back(0); w1_.push_back(0); denom_.push_back(1); continue; }
      const double va = variance[a], vb = variance[b];
      if (!(va > 0) || !(vb > 0) || !std::isfinite(va) || !std::isfinite(vb))
        Rcpp::stop("phylo_d requires strictly positive branch lengths.");
      if (arithmetic == scaled) {
        const double scale = std::max(va, vb), x = va/scale, y = vb/scale;
        w0_.push_back(y/(x+y)); w1_.push_back(1-w0_.back()); denom_.push_back(1);
        variance[node] += scale*x*y/(x+y);
      } else {
        w0_.push_back(1/va); w1_.push_back(1/vb); denom_.push_back(w0_.back()+w1_.back());
        variance[node] += (va*vb)/(va+vb);
      }
    }
  }

  // Edge groups are already in the public D postorder. Unary nodes retain the
  // legacy copy-through behavior; the validated streaming path is binary.
  static DContrasts from_edges(const Rcpp::IntegerMatrix& edge,
      const double* lengths, int nt) {
    int nn = nt;
    for (int i = 0; i < edge.nrow(); ++i) nn = std::max(nn, std::max(edge(i,0), edge(i,1)));
    std::vector<int> post, left(nn, -1), right(nn, -1);
    std::vector<double> branch(nn, 0);
    for (int i = 0; i < edge.nrow(); ++i) {
      const int node = edge(i,0)-1, child = edge(i,1)-1;
      if (node < nt || child < 0) Rcpp::stop("invalid D edge node");
      branch[child] = lengths[i];
      if (i == 0 || edge(i,0) != edge(i-1,0)) post.push_back(node);
      if (left[node] < 0) left[node] = child;
      else if (right[node] < 0) right[node] = child;
      else Rcpp::stop("binary D kernel received a polytomy");
    }
    return DContrasts(nt, nn, post, left, right, branch);
  }

  int n_total() const { return nn_; }
  template<class Sum = double>
  double sum(const double* tips, double* state) const {
    std::copy(tips, tips+nt_, state);
    Sum total = 0;
    for (size_t k = 0; k < nodes_.size(); ++k) {
      const int node = nodes_[k], a = left_[k], b = right_[k];
      if (b < 0) { state[node] = state[a]; continue; }
      const double mean = arithmetic_ == scaled ?
        state[a]*w0_[k]+state[b]*w1_[k] :
        (state[a]*w0_[k]+state[b]*w1_[k])/denom_[k];
      state[node] = mean;
      total += std::abs(state[a]-mean)+std::abs(state[b]-mean);
    }
    return static_cast<double>(total);
  }

  std::vector<double> batch(const double* tips, int ncol, int n_threads) const {
    std::vector<double> work(static_cast<size_t>(nn_)*ncol), sums(ncol);
#ifdef _OPENMP
    const int threads = std::max(1, n_threads);
#pragma omp parallel for num_threads(threads) if (threads > 1 && ncol > 1) schedule(static)
#else
    (void)n_threads;
#endif
    for (int j = 0; j < ncol; ++j)
      sums[j] = sum(tips+static_cast<size_t>(nt_)*j, work.data()+static_cast<size_t>(nn_)*j);
    return sums;
  }

 private:
  int nt_, nn_;
  Arithmetic arithmetic_;
  std::vector<int> nodes_, left_, right_;
  std::vector<double> w0_, w1_, denom_;
};
} // namespace fastphylosig
#endif
