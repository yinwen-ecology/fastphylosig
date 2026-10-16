// Binary D crowns and a homogeneous two-state ER reference conditioned on
// each crown's observed state count. Polynomial pruning fits this conditional
// likelihood; a downward pass samples exact counts without rejection. Tree
// factors, count probabilities and D normalization are reused across draws.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "heterogeneity_bootstrap.h"
#include "d_kernels.h"

namespace {
struct DTree {
  int nt, nn, root;
  std::vector<int> parent, pre, post, left, right, count;
  std::vector<double> branch, branch_sd;
  mutable std::vector<int> simulated_state, simulated_count;
  mutable std::vector<double> simulated_tips;
  fastphylosig::DContrasts contrasts;
  mutable std::vector<double> contrast_work;
  double total = 0;
  explicit DTree(const Rcpp::List& core) {
    nt = Rcpp::as<int>(core["n_tip"]); nn = Rcpp::as<int>(core["n_total"]);
    root = Rcpp::as<int>(core["root"]) - 1;
    parent = Rcpp::as<std::vector<int> >(core["parent"]);
    pre = Rcpp::as<std::vector<int> >(core["preorder"]);
    post = Rcpp::as<std::vector<int> >(core["postorder"]);
    branch = Rcpp::as<std::vector<double> >(core["branch_length_by_node"]);
    const std::vector<int> ptr = Rcpp::as<std::vector<int> >(core["child_ptr"]);
    const std::vector<int> child = Rcpp::as<std::vector<int> >(core["children"]);
    if (nt < 3 || nn != 2*nt-1 || root < nt || root >= nn ||
        parent.size() != static_cast<size_t>(nn) || pre.size() != parent.size() ||
        post.size() != parent.size() || branch.size() != parent.size() ||
        ptr.size() != parent.size()+1 || child.size()+1 != parent.size() ||
        ptr.front() != 1 || ptr.back() != nn) Rcpp::stop("invalid binary crown core");
    left.assign(nn, -1); right.assign(nn, -1); count.assign(nn, 0);
    for (int node = 0; node < nn; ++node) {
      --parent[node]; --pre[node]; --post[node];
      if (node != root && (parent[node] < nt || parent[node] >= nn ||
          !std::isfinite(branch[node]) || branch[node] <= 0))
        Rcpp::stop("D crowns require strictly positive branch lengths");
      if (node != root) total += branch[node];
      if (ptr[node] < 1 || ptr[node+1] < ptr[node] || ptr[node+1] > nn ||
          ptr[node+1]-ptr[node] != (node < nt ? 0 : 2)) Rcpp::stop("D crowns must be bifurcating");
      if (node >= nt) {
        left[node] = child[ptr[node]-1]-1; right[node] = child[ptr[node]]-1;
        if (left[node] < 0 || right[node] < 0 || left[node] >= nn || right[node] >= nn ||
            left[node] == right[node]) Rcpp::stop("invalid binary crown children");
      }
    }
    std::vector<bool> seen(nn, false);
    for (int node : pre) {
      if (node < 0 || node >= nn || seen[node] ||
          (node != root && !seen[parent[node]])) Rcpp::stop("invalid binary crown preorder");
      seen[node] = true;
    }
    std::fill(seen.begin(), seen.end(), false);
    for (int node : post) {
      if (node < 0 || node >= nn || seen[node]) Rcpp::stop("invalid binary crown postorder");
      if (node < nt) count[node] = 1;
      else {
        const int a = left[node], b = right[node];
        if (!seen[a] || !seen[b] || parent[a] != node || parent[b] != node)
          Rcpp::stop("invalid binary crown topology");
        count[node] = count[a] + count[b];
      }
      seen[node] = true;
    }
    if (!std::isfinite(total) || total <= 0 || count[root] != nt)
      Rcpp::stop("invalid binary crown branch scale");
    contrasts = fastphylosig::DContrasts(nt, nn, post, left, right, branch,
      fastphylosig::DContrasts::scaled);
    contrast_work.resize(nn); branch_sd.resize(nn);
    for (int node : pre) if (node != root) branch_sd[node] = std::sqrt(branch[node]);
  }
  double contrast(const std::vector<double>& tips) const {
    return contrasts.sum<long double>(tips.data(), contrast_work.data());
  }
  struct CountWeights {
    std::vector<std::vector<long double> > p0, p1;
    std::vector<double> change, likelihood0, likelihood1;
    std::vector<long double> mixed0, mixed1;
  };
  mutable CountWeights fit_work;
  void count_weights(double q, int m, CountWeights& w) const {
    w.p0.resize(nn); w.p1.resize(nn); w.change.resize(nn);
    for (int node : pre) if (node != root) w.change[node] = -.5*std::expm1(-2*q*branch[node]);
    for (int node : post) {
      w.p0[node].assign(std::min(count[node], m)+1, 0);
      w.p1[node].assign(w.p0[node].size(), 0);
      if (node < nt) { w.p0[node][0] = 1; w.p1[node][1] = 1; continue; }
      const int a = left[node], b = right[node];
      const long double ca = w.change[a], cb = w.change[b];
      // Reuse child transition mixtures without changing convolution order.
      w.mixed0.resize(w.p0[b].size()); w.mixed1.resize(w.p0[b].size());
      for (size_t j = 0; j < w.p0[b].size(); ++j) {
        w.mixed0[j] = (1-cb)*w.p0[b][j]+cb*w.p1[b][j];
        w.mixed1[j] = cb*w.p0[b][j]+(1-cb)*w.p1[b][j];
      }
      for (size_t i = 0; i < w.p0[a].size(); ++i) {
        const long double left0 = (1-ca)*w.p0[a][i]+ca*w.p1[a][i];
        const long double left1 = ca*w.p0[a][i]+(1-ca)*w.p1[a][i];
        for (size_t j = 0; j < w.p0[b].size() && i+j < w.p0[node].size(); ++j) {
          w.p0[node][i+j] += left0 * w.mixed0[j];
          w.p1[node][i+j] += left1 * w.mixed1[j];
        }
      }
    }
  }
  double log_likelihood(const std::vector<double>& tips, CountWeights& w) const {
    w.likelihood0.resize(nn); w.likelihood1.resize(nn);
    std::vector<double>& p0 = w.likelihood0;
    std::vector<double>& p1 = w.likelihood1;
    for (int tip = 0; tip < nt; ++tip) { p0[tip] = tips[tip] == 0; p1[tip] = tips[tip] == 1; }
    long double log_scale = 0;
    for (int node : post) if (node >= nt) {
      double a0 = 1, a1 = 1;
      for (int child : {left[node], right[node]}) {
        const double change = w.change[child];
        a0 *= (1-change)*p0[child] + change*p1[child];
        a1 *= change*p0[child] + (1-change)*p1[child];
      }
      const double scale = std::max(a0, a1);
      if (scale <= 0) return R_NegInf;
      p0[node] = a0/scale; p1[node] = a1/scale; log_scale += std::log(scale);
    }
    return log_scale + std::log(0.5*(p0[root]+p1[root]));
  }
  double conditional_log_likelihood(const std::vector<double>& tips, double q) const {
    const int m = std::accumulate(tips.begin(), tips.end(), 0);
    count_weights(q, m, fit_work);
    const long double mass = .5L*(fit_work.p0[root][m]+fit_work.p1[root][m]);
    return mass > 0 ? log_likelihood(tips, fit_work)-std::log(mass) : R_NegInf;
  }
  const std::vector<double>& simulate_count(int m, const CountWeights& w) const {
    simulated_state.resize(nn); simulated_count.resize(nn); simulated_tips.resize(nt);
    std::vector<int>& state = simulated_state;
    std::vector<int>& counts = simulated_count;
    counts[root] = m;
    const long double root_mass = w.p0[root][m]+w.p1[root][m];
    if (!(root_mass > 0)) Rcpp::stop("binary count reference has zero probability");
    state[root] = R::runif(0, 1)*root_mass >= w.p0[root][m];
    for (int node : pre) if (node >= nt) {
      const int a = left[node], b = right[node], total_count = counts[node];
      const long double ca = w.change[a], cb = w.change[b];
      const auto message = [&](int child, int k, long double change) {
        const long double s = state[node] == 0 ? 1-change : change;
        return s*w.p0[child][k]+(1-s)*w.p1[child][k];
      };
      const int lo = std::max(0, total_count-static_cast<int>(w.p0[b].size())+1);
      const int hi = std::min(total_count, static_cast<int>(w.p0[a].size())-1);
      long double mass = 0;
      for (int k = lo; k <= hi; ++k) mass += message(a, k, ca)*message(b, total_count-k, cb);
      if (!(mass > 0)) Rcpp::stop("binary count reference split has zero probability");
      const long double u = R::runif(0, 1)*mass;
      long double cumulative = 0; int split = hi;
      for (int k = lo; k <= hi; ++k) {
        cumulative += message(a, k, ca)*message(b, total_count-k, cb);
        if (u < cumulative) { split = k; break; }
      }
      counts[a] = split; counts[b] = total_count-split;
      for (int child : {a, b}) {
        const long double change = child == a ? ca : cb;
        const long double s = state[node] == 0 ? 1-change : change;
        const int k = counts[child];
        const long double p0 = s*w.p0[child][k], p1 = (1-s)*w.p1[child][k];
        if (!(p0+p1 > 0)) Rcpp::stop("binary count reference state has zero probability");
        state[child] = R::runif(0, 1)*(p0+p1) >= p0;
      }
    }
    std::copy(state.begin(), state.begin()+nt, simulated_tips.begin());
    if (std::accumulate(simulated_tips.begin(), simulated_tips.end(), 0) != m) Rcpp::stop("binary reference did not preserve state count");
    return simulated_tips;
  }
};

// Golden-section refinement of all coarse-grid local maxima (endpoints
// included). A dimensionless lower bound 1e-8 avoids the undefined q=0
// conditional event. Infinite q is the exact independent-state endpoint.
double fit_rate(const std::vector<DTree>& trees,
    const std::vector<std::vector<double> >& data, const std::vector<int>& selected, double scale) {
  const auto objective = [&](double r) {
    const double q = r == 1 ? R_PosInf : std::max(1e-8, r/(1-r))/scale;
    long double ll = 0;
    for (int i : selected) ll += trees[i].conditional_log_likelihood(data[i], q);
    return static_cast<double>(ll);
  };
  std::vector<double> grid(65);
  double best = R_NegInf, best_r = NA_REAL;
  for (int k = 0; k <= 64; ++k) {
    grid[k] = objective(k/64.0);
    if (grid[k] > best) { best = grid[k]; best_r = k/64.0; }
  }
  for (int k = 1; k < 64; ++k) if (grid[k] >= grid[k-1] && grid[k] >= grid[k+1] &&
      (grid[k] > grid[k-1] || grid[k] > grid[k+1])) {
    double lo = (k-1)/64.0, hi = (k+1)/64.0;
    double a = hi-(hi-lo)*0.6180339887498949, b = lo+(hi-lo)*0.6180339887498949;
    double fa = objective(a), fb = objective(b);
    for (int iter = 0; iter < 70 && hi-lo > 1e-8; ++iter) {
      if (fa > fb) { hi = b; b = a; fb = fa; a = hi-(hi-lo)*0.6180339887498949; fa = objective(a); }
      else { lo = a; a = b; fa = fb; b = lo+(hi-lo)*0.6180339887498949; fb = objective(b); }
    }
    const double r = fa > fb ? a : b, value = std::max(fa, fb);
    if (value > best) { best = value; best_r = r; }
  }
  return best_r == 1 ? R_PosInf : std::max(1e-8, best_r/(1-best_r))/scale;
}
}

// [[Rcpp::export]]
Rcpp::List heterogeneity_d_calibrate_cpp(const Rcpp::List& compiled,
    const Rcpp::List& data, const int d_nsim, const int nsim,
    const bool test, const bool keep_null) {
  const int n = compiled.size();
  if (n < 1 || data.size() != n || d_nsim < 2 || (test && nsim < 2))
    Rcpp::stop("invalid D heterogeneity dimensions");
  std::vector<DTree> trees;
  std::vector<std::vector<double> > values;
  std::vector<int> state_counts;
  Rcpp::NumericVector estimates(n, NA_REAL), contrasts(n), mr(n), mb(n);
  Rcpp::NumericVector p_random(n, NA_REAL), p_brownian(n, NA_REAL);
  Rcpp::CharacterVector status(n, "ok");
  std::vector<int> selected;
  for (int i = 0; i < n; ++i) {
    Rcpp::checkUserInterrupt(); trees.emplace_back(Rcpp::List(compiled[i]));
    const DTree& tree = trees.back();
    values.push_back(Rcpp::as<std::vector<double> >(data[i]));
    if (values.back().size() != static_cast<size_t>(tree.nt)) Rcpp::stop("incompatible D crown data");
    int n1 = 0;
    for (double value : values.back()) {
      if (value != 0 && value != 1) Rcpp::stop("D crown data require 0/1 states");
      n1 += value == 1;
    }
    if (n1 == 0 || n1 == tree.nt) Rcpp::stop("D crown data must contain both states");
    // Use minority states and a fixed orientation to preserve relabelling
    // invariance even for finite Monte Carlo normalization and tied statistics.
    if (n1 > tree.nt-n1 || (n1 == tree.nt-n1 && values.back()[0] == 1)) {
      for (double& value : values.back()) value = 1-value;
      n1 = tree.nt-n1;
    }
    state_counts.push_back(n1);
    if (test && tree.nt > 1000) Rcpp::stop("D count-conditioned testing supports at most 1000 observed tips per crown; use smaller crowns");
    contrasts[i] = tree.contrast(values.back());
    long double rs = 0, bs = 0;
    int random_less = 0, brownian_greater = 0;
    std::vector<int> order(tree.nt);
    std::vector<double> bm(tree.nn, 0), random_state(tree.nt, 0), bm_state(tree.nt, 0);
    for (int b = 0; b < d_nsim; ++b) {
      if (b % 128 == 0) Rcpp::checkUserInterrupt();
      std::iota(order.begin(), order.end(), 0);
      for (int k = tree.nt-1; k > 0; --k) std::swap(order[k], order[static_cast<int>(R::runif(0, k+1))]);
      std::fill(random_state.begin(), random_state.end(), 0);
      for (int k = 0; k < n1; ++k) random_state[order[k]] = 1;
      const double random_contrast = tree.contrast(random_state);
      rs += random_contrast;
      random_less += random_contrast < contrasts[i];
      bm[tree.root] = 0;
      for (int node : tree.pre) if (node != tree.root)
        bm[node] = bm[tree.parent[node]] + R::rnorm(0, 1)*tree.branch_sd[node];
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](int a, int b) { return bm[a] < bm[b]; });
      std::fill(bm_state.begin(), bm_state.end(), 0);
      for (int k = 0; k < n1; ++k) bm_state[order[k]] = 1;
      const double bm_contrast = tree.contrast(bm_state);
      bs += bm_contrast;
      brownian_greater += bm_contrast > contrasts[i];
    }
    mr[i] = rs/d_nsim; mb[i] = bs/d_nsim;
    p_random[i] = static_cast<double>(random_less)/d_nsim;
    p_brownian[i] = static_cast<double>(brownian_greater)/d_nsim;
    const double denom = mr[i]-mb[i];
    if (std::isfinite(denom) && std::abs(denom) > 64*std::numeric_limits<double>::epsilon()*std::max(1.0, std::abs(mr[i]))) {
      estimates[i] = (contrasts[i]-mb[i])/denom; selected.push_back(i);
    } else status[i] = "undefined_normalization";
  }
  std::vector<double> observed;
  double scale = 0; int edges = 0;
  for (int i : selected) {
    observed.push_back(estimates[i]);
    scale += trees[i].total; edges += trees[i].nn-1;
  }
  const double h = fastphylosig::hetero_sd(observed);
  double rate = NA_REAL;
  std::vector<double> draws(test ? nsim : 0, NA_REAL), simulated(selected.size());
  if (test && selected.size() >= 2) {
    rate = fit_rate(trees, values, selected, scale/edges);
    if (!std::isnan(rate)) for (int i : selected)
      trees[i].count_weights(rate, state_counts[i], trees[i].fit_work);
    if (!std::isnan(rate)) for (int b = 0; b < nsim; ++b) {
      if (b % 32 == 0) Rcpp::checkUserInterrupt();
      for (size_t k = 0; k < selected.size(); ++k) {
        const int i = selected[k];
        const std::vector<double>& x = trees[i].simulate_count(state_counts[i], trees[i].fit_work);
        simulated[k] = (trees[i].contrast(x)-mb[i])/(mr[i]-mb[i]);
      }
      draws[b] = fastphylosig::hetero_sd(simulated);
    }
  }
  Rcpp::List out = fastphylosig::hetero_bootstrap_summary(draws, h, keep_null);
  out["estimate"] = estimates; out["status"] = status;
  out["contrast"] = contrasts; out["D_mean_random"] = mr; out["D_mean_Brownian"] = mb;
  out["dispersion"] = h; out["q_common"] = rate;
  out["P_random"] = p_random; out["P_Brownian"] = p_brownian;
  return out;
}
