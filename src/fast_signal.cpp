// Fast numerical kernels for phytools-compatible phylosig calculations.
//
// K and lambda keep the same statistical targets as phytools::phylosig().
// Dense inverse matrices are avoided where possible: Cholesky triangular solves
// evaluate the same C^{-1} GLS means and quadratic forms more cheaply.

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

int effective_threads(const int n_threads) {
#ifdef _OPENMP
  return n_threads > 1 ? n_threads : 1;
#else
  return 1;
#endif
}

arma::vec chol_solve_upper(const arma::mat& cholC, const arma::vec& b) {
  const arma::vec y = arma::solve(arma::trimatl(cholC.t()), b);
  return arma::solve(arma::trimatu(cholC), y);
}

double chol_quad_upper(const arma::mat& cholC, const arma::vec& x) {
  const arma::vec z = arma::solve(arma::trimatl(cholC.t()), x);
  return arma::dot(z, z);
}

double k_stat_chol_one(const arma::vec& x,
                       const arma::mat& cholC,
                       const arma::vec& inv_one,
                       const double norm_const,
                       const double sum_invC) {
  const double a = arma::dot(inv_one, x) / sum_invC;
  const arma::vec centered = x - a;
  const double numerator = arma::dot(centered, centered);
  const double denominator = chol_quad_upper(cholC, centered);
  return (numerator / denominator) / norm_const;
}

arma::mat lambda_transform(const arma::mat& C, const double lambda) {
  arma::mat Cl = lambda * C;
  Cl.diag() = C.diag();
  return Cl;
}

double lambda_loglik_internal(const double lambda,
                              const arma::mat& C,
                              const arma::vec& y) {
  const unsigned int n = C.n_rows;
  arma::mat Cl = lambda_transform(C, lambda);

  arma::mat cholCl;
  const bool chol_ok = arma::chol(cholCl, Cl);
  if (!chol_ok || !cholCl.is_finite()) {
    return -std::numeric_limits<double>::infinity();
  }

  const arma::vec ones(n, arma::fill::ones);
  const arma::vec inv_one = chol_solve_upper(cholCl, ones);
  if (!inv_one.is_finite()) {
    return -std::numeric_limits<double>::infinity();
  }
  const double sum_invCl = arma::accu(inv_one);
  const double a = arma::dot(inv_one, y) / sum_invCl;
  const arma::vec centered = y - a;
  const double quad = chol_quad_upper(cholCl, centered);
  const double sig2 = quad / static_cast<double>(n);
  if (!std::isfinite(sig2) || sig2 <= 0.0) {
    return -std::numeric_limits<double>::infinity();
  }

  const double logdetCl = 2.0 * arma::accu(arma::log(cholCl.diag()));
  if (!std::isfinite(logdetCl)) {
    return -std::numeric_limits<double>::infinity();
  }

  const double log_two_pi = std::log(6.283185307179586476925286766559);
  const double logdet_sig2_Cl =
    static_cast<double>(n) * std::log(sig2) + logdetCl;

  return -0.5 * quad / sig2 -
    0.5 * static_cast<double>(n) * log_two_pi -
    0.5 * logdet_sig2_Cl;
}

arma::mat ace_build_q(const Rcpp::IntegerMatrix& rate_index,
                      const arma::vec& par) {
  const int k = rate_index.nrow();
  arma::mat Q(k, k, arma::fill::zeros);
  for (int j = 0; j < k; ++j) {
    for (int i = 0; i < k; ++i) {
      if (i == j) continue;
      const int idx = rate_index(i, j);
      if (idx > 0) {
        if (idx > static_cast<int>(par.n_elem)) {
          Rcpp::stop("rate_index refers to a missing parameter.");
        }
        Q(i, j) = par(idx - 1);
      }
    }
  }
  for (int i = 0; i < k; ++i) {
    Q(i, i) = -arma::accu(Q.row(i));
  }
  return Q;
}

bool ace_eigen_decomp(const arma::mat& Q,
                      arma::cx_vec& lambda,
                      arma::cx_mat& gamma,
                      arma::cx_mat& inv_gamma) {
  const bool ok = arma::eig_gen(lambda, gamma, Q);
  if (!ok || !lambda.is_finite() || !gamma.is_finite()) {
    return false;
  }
  inv_gamma = arma::inv(gamma);
  return inv_gamma.is_finite();
}

arma::vec ace_transition_multiply(const arma::cx_vec& lambda,
                                  const arma::cx_mat& gamma,
                                  const arma::cx_mat& inv_gamma,
                                  const arma::vec& child_lik,
                                  const double edge_length) {
  arma::cx_vec z = inv_gamma * arma::conv_to<arma::cx_vec>::from(child_lik);
  z %= arma::exp(lambda * edge_length);
  return arma::real(gamma * z);
}

arma::mat ace_transition_matrix(const arma::cx_vec& lambda,
                                const arma::cx_mat& gamma,
                                const arma::cx_mat& inv_gamma,
                                const double edge_length) {
  const arma::cx_mat P =
    gamma * arma::diagmat(arma::exp(lambda * edge_length)) * inv_gamma;
  return arma::real(P);
}

double ace_discrete_pruning_internal(const Rcpp::IntegerMatrix& edge,
                                     const arma::vec& edge_length,
                                     const Rcpp::IntegerVector& tip_state,
                                     const Rcpp::IntegerMatrix& rate_index,
                                     const arma::vec& par,
                                     arma::mat* out_liks = NULL,
                                     const bool marginal = false) {
  const int n_tip = tip_state.size();
  const int n_node = n_tip - 1;
  const int n_total = n_tip + n_node;
  const int n_edge = edge.nrow();
  const int k = rate_index.nrow();

  if (edge.ncol() != 2 || n_edge != 2 * n_node) {
    Rcpp::stop("edge must be a postorder binary-tree edge matrix.");
  }
  if (edge_length.n_elem != static_cast<arma::uword>(n_edge)) {
    Rcpp::stop("edge_length must have one value per edge.");
  }
  if (rate_index.ncol() != k) {
    Rcpp::stop("rate_index must be square.");
  }

  arma::mat liks(n_total, k, arma::fill::zeros);
  for (int i = 0; i < n_tip; ++i) {
    if (Rcpp::IntegerVector::is_na(tip_state[i])) {
      liks.row(i).ones();
    } else {
      const int state = tip_state[i];
      if (state < 1 || state > k) {
        Rcpp::stop("tip_state values must be 1-based state indices.");
      }
      liks(i, state - 1) = 1.0;
    }
  }

  const arma::mat Q = ace_build_q(rate_index, par);
  arma::cx_vec lambda;
  arma::cx_mat gamma;
  arma::cx_mat inv_gamma;
  if (!ace_eigen_decomp(Q, lambda, gamma, inv_gamma)) {
    return std::numeric_limits<double>::infinity();
  }

  double log_comp = 0.0;
  for (int e = 0; e < n_edge; e += 2) {
    const int e_right = e + 1;
    const int anc = edge(e, 0) - 1;
    const int des_left = edge(e, 1) - 1;
    const int des_right = edge(e_right, 1) - 1;

    if (anc < n_tip || anc >= n_total ||
        des_left < 0 || des_left >= n_total ||
        des_right < 0 || des_right >= n_total) {
      Rcpp::stop("edge contains node numbers outside the expected range.");
    }

    const arma::vec v_left = ace_transition_multiply(
      lambda, gamma, inv_gamma, liks.row(des_left).t(), edge_length(e)
    );
    const arma::vec v_right = ace_transition_multiply(
      lambda, gamma, inv_gamma, liks.row(des_right).t(), edge_length(e_right)
    );
    const arma::vec v = v_left % v_right;
    const double comp = arma::accu(v);
    if (!std::isfinite(comp) || comp <= 0.0) {
      return std::numeric_limits<double>::infinity();
    }
    liks.row(anc) = (v / comp).t();
    log_comp += std::log(comp);
  }

  if (out_liks != NULL) {
    arma::mat lik_anc = liks.rows(n_tip, n_total - 1);
    if (marginal) {
      *out_liks = lik_anc;
    } else {
      for (int e = n_edge - 2; e >= 0; e -= 2) {
        const int e_right = e + 1;
        const int anc = edge(e, 0) - n_tip - 1;

        const int des_left = edge(e, 1) - n_tip - 1;
        if (des_left >= 0) {
          const arma::mat P = ace_transition_matrix(
            lambda, gamma, inv_gamma, edge_length(e)
          );
          const arma::rowvec denom = lik_anc.row(des_left) * P;
          if (arma::any(denom == 0.0)) {
            return std::numeric_limits<double>::infinity();
          }
          const arma::rowvec tmp = lik_anc.row(anc) / denom;
          lik_anc.row(des_left) = (tmp * P) % lik_anc.row(des_left);
        }

        const int des_right = edge(e_right, 1) - n_tip - 1;
        if (des_right >= 0) {
          const arma::mat P = ace_transition_matrix(
            lambda, gamma, inv_gamma, edge_length(e_right)
          );
          const arma::rowvec denom = lik_anc.row(des_right) * P;
          if (arma::any(denom == 0.0)) {
            return std::numeric_limits<double>::infinity();
          }
          const arma::rowvec tmp = lik_anc.row(anc) / denom;
          lik_anc.row(des_right) = (tmp * P) % lik_anc.row(des_right);
        }

        const arma::vec rowsums = arma::sum(lik_anc, 1);
        for (arma::uword r = 0; r < lik_anc.n_rows; ++r) {
          if (rowsums(r) != 0.0) {
            lik_anc.row(r) /= rowsums(r);
          }
        }
      }
      *out_liks = lik_anc;
    }
  }

  return -2.0 * log_comp;
}

arma::vec delta_entropy_internal(const arma::mat& probabilities,
                                 const int entropy_type) {
  const int n = static_cast<int>(probabilities.n_rows);
  const int k = static_cast<int>(probabilities.n_cols);
  if (k < 2) {
    Rcpp::stop("Delta requires at least two states.");
  }

  arma::vec out(n, arma::fill::zeros);
  const double kd = static_cast<double>(k);
  const double uniform = 1.0 / kd;

  for (int i = 0; i < n; ++i) {
    if (entropy_type == 1) {
      double total = 0.0;
      for (int j = 0; j < k; ++j) {
        double p = probabilities(i, j);
        if (p > uniform) {
          p = p / (1.0 - kd) - 1.0 / (1.0 - kd);
        }
        total += p;
      }
      out(i) = total;
    } else if (entropy_type == 2) {
      double total = 0.0;
      for (int j = 0; j < k; ++j) {
        const double p = probabilities(i, j);
        if (p > 0.0) {
          total -= p * std::log(p);
        }
      }
      out(i) = total / std::log(kd);
    } else {
      double sumsq = 0.0;
      for (int j = 0; j < k; ++j) {
        const double p = probabilities(i, j);
        sumsq += p * p;
      }
      out(i) = ((1.0 - sumsq) * kd) / (kd - 1.0);
    }
  }

  bool has_zero = false;
  bool has_one = false;
  for (int i = 0; i < n; ++i) {
    if (out(i) == 0.0) has_zero = true;
    if (out(i) == 1.0) has_one = true;
  }
  if (has_zero) {
    const double eps = R::runif(0.0, 1.0) / 10000.0;
    for (int i = 0; i < n; ++i) {
      if (out(i) == 0.0) out(i) += eps;
    }
  }
  if (has_one) {
    const double eps = R::runif(0.0, 1.0) / 10000.0;
    for (int i = 0; i < n; ++i) {
      if (out(i) == 1.0) out(i) -= eps;
    }
  }
  return out;
}

double delta_lp_alpha(const double a,
                      const double b,
                      const arma::vec& x,
                      const double lambda0,
                      const double sum_log_x) {
  const double n = static_cast<double>(x.n_elem);
  return n * (std::lgamma(a + b) - std::lgamma(a)) -
    a * (lambda0 - sum_log_x);
}

double delta_lp_beta(const double a,
                     const double b,
                     const arma::vec& x,
                     const double lambda0,
                     const double sum_log_1_minus_x) {
  const double n = static_cast<double>(x.n_elem);
  return n * (std::lgamma(a + b) - std::lgamma(b)) -
    b * (lambda0 - sum_log_1_minus_x);
}

double delta_mh_alpha(const double a,
                      const double b,
                      const arma::vec& x,
                      const double lambda0,
                      const double proposal_sd,
                      const double sum_log_x) {
  const double current = delta_lp_alpha(a, b, x, lambda0, sum_log_x);
  for (int attempt = 0; attempt < 10000; ++attempt) {
    const double candidate = std::exp(R::rnorm(std::log(a), proposal_sd));
    const double proposed = delta_lp_alpha(candidate, b, x, lambda0, sum_log_x);
    const double log_ratio = proposed - current;
    if (std::isfinite(log_ratio)) {
      if (log_ratio >= 0.0 || std::log(R::runif(0.0, 1.0)) < log_ratio) {
        return candidate;
      }
      return a;
    }
  }
  return a;
}

double delta_mh_beta(const double a,
                     const double b,
                     const arma::vec& x,
                     const double lambda0,
                     const double proposal_sd,
                     const double sum_log_1_minus_x) {
  const double current =
    delta_lp_beta(a, b, x, lambda0, sum_log_1_minus_x);
  for (int attempt = 0; attempt < 10000; ++attempt) {
    const double candidate = std::exp(R::rnorm(std::log(b), proposal_sd));
    const double proposed =
      delta_lp_beta(a, candidate, x, lambda0, sum_log_1_minus_x);
    const double log_ratio = proposed - current;
    if (std::isfinite(log_ratio)) {
      if (log_ratio >= 0.0 || std::log(R::runif(0.0, 1.0)) < log_ratio) {
        return candidate;
      }
      return b;
    }
  }
  return b;
}

arma::vec delta_chain_sums(double alpha,
                           double beta,
                           const arma::vec& x,
                           const double lambda0,
                           const double proposal_sd,
                           const int sim,
                           const int thin,
                           const int burn,
                           const double sum_log_x,
                           const double sum_log_1_minus_x) {
  double alpha_sum = 0.0;
  double beta_sum = 0.0;
  double saved = 0.0;

  for (int iter = 1; iter <= sim; ++iter) {
    alpha = delta_mh_alpha(
      alpha, beta, x, lambda0, proposal_sd, sum_log_x
    );
    beta = delta_mh_beta(
      alpha, beta, x, lambda0, proposal_sd, sum_log_1_minus_x
    );

    if (iter >= burn && ((iter - burn) % thin == 0)) {
      alpha_sum += alpha;
      beta_sum += beta;
      saved += 1.0;
    }
  }

  return arma::vec({alpha_sum, beta_sum, saved});
}

bool phylo_d_has_polytomy(const Rcpp::IntegerMatrix& edge,
                          const int n_tip) {
  int max_node = n_tip;
  const int n_edge = edge.nrow();
  for (int i = 0; i < n_edge; ++i) {
    max_node = std::max(max_node, edge(i, 0));
  }

  std::vector<int> child_count(max_node + 1, 0);
  for (int i = 0; i < n_edge; ++i) {
    const int parent = edge(i, 0);
    if (parent >= 0 && parent <= max_node) {
      ++child_count[parent];
      if (child_count[parent] > 2) {
        return true;
      }
    }
  }
  return false;
}

double phylo_d_sum_one_binary(const arma::vec& state,
                              const Rcpp::IntegerMatrix& edge,
                              const arma::vec& edge_length,
                              const int n_tip) {
  const int n_edge = edge.nrow();
  int max_node = n_tip;
  for (int i = 0; i < n_edge; ++i) {
    max_node = std::max(max_node, edge(i, 0));
    max_node = std::max(max_node, edge(i, 1));
  }
  const int root = n_tip + 1;

  std::vector<double> node_value(max_node + 1, 0.0);
  for (int i = 0; i < n_tip; ++i) {
    node_value[i + 1] = state(i);
  }

  arma::vec edge_len = edge_length;
  std::vector<int> child_to_edge(max_node + 1, -1);
  for (int i = 0; i < n_edge; ++i) {
    child_to_edge[edge(i, 1)] = i;
  }

  double sum = 0.0;
  int row = 0;
  while (row < n_edge) {
    const int parent = edge(row, 0);
    int child0 = -1;
    int child1 = -1;
    double bl0 = 0.0;
    double bl1 = 0.0;
    int n_child = 0;

    while (row < n_edge && edge(row, 0) == parent) {
      if (n_child == 0) {
        child0 = edge(row, 1);
        bl0 = edge_len(row);
      } else if (n_child == 1) {
        child1 = edge(row, 1);
        bl1 = edge_len(row);
      }
      ++n_child;
      ++row;
    }

    if (n_child < 2) {
      if (n_child == 1) {
        node_value[parent] = node_value[child0];
      }
      continue;
    }

    const double inv0 = 1.0 / bl0;
    const double inv1 = 1.0 / bl1;
    const double denom = inv0 + inv1;
    const double curr_nv =
      (node_value[child0] * inv0 + node_value[child1] * inv1) / denom;
    const double curr_contr =
      std::abs(node_value[child0] - curr_nv) +
      std::abs(node_value[child1] - curr_nv);
    const double curr_bl_adj = (bl0 * bl1) / (bl0 + bl1);

    node_value[parent] = curr_nv;
    sum += curr_contr;

    if (parent != root) {
      const int parent_edge = child_to_edge[parent];
      if (parent_edge >= 0) {
        edge_len(parent_edge) += curr_bl_adj;
      }
    }
  }

  return sum;
}

} // namespace

// [[Rcpp::export]]
arma::vec fast_k_chol_batch_cpp(const arma::mat& X,
                                const arma::mat& cholC,
                                const double traceC) {
  if (X.n_rows != cholC.n_rows || cholC.n_rows != cholC.n_cols) {
    Rcpp::stop("X rows must match the square cholC matrix.");
  }

  const arma::vec ones(cholC.n_rows, arma::fill::ones);
  const arma::vec inv_one = chol_solve_upper(cholC, ones);
  const double sum_invC = arma::accu(inv_one);
  const int n = static_cast<int>(cholC.n_rows);
  const double norm_const =
    (traceC - static_cast<double>(n) / sum_invC) /
    static_cast<double>(n - 1);

  arma::vec out(X.n_cols);
  for (arma::uword j = 0; j < X.n_cols; ++j) {
    out(j) = k_stat_chol_one(X.col(j), cholC, inv_one,
                             norm_const, sum_invC);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List fast_k_chol_permutation_cpp(
    const arma::mat& X,
    const arma::mat& cholC,
    const double traceC,
    const Rcpp::IntegerMatrix& permutations,
    const int n_threads = 1) {
  if (X.n_rows != cholC.n_rows || cholC.n_rows != cholC.n_cols) {
    Rcpp::stop("X rows must match the square cholC matrix.");
  }
  if (permutations.ncol() != static_cast<int>(X.n_rows)) {
    Rcpp::stop("permutations must have one column per species.");
  }

  const int nsim = permutations.nrow();
  const int n = permutations.ncol();
  const int p = static_cast<int>(X.n_cols);
  const arma::imat perms = Rcpp::as<arma::imat>(permutations);
  for (int i = 0; i < nsim; ++i) {
    for (int r = 0; r < n; ++r) {
      const int idx = perms(i, r);
      if (idx < 1 || idx > n) {
        Rcpp::stop("permutations must contain 1-based indices in 1:n.");
      }
    }
  }

  const arma::vec ones(cholC.n_rows, arma::fill::ones);
  const arma::vec inv_one = chol_solve_upper(cholC, ones);
  const double sum_invC = arma::accu(inv_one);
  const double norm_const =
    (traceC - static_cast<double>(n) / sum_invC) /
    static_cast<double>(n - 1);
  const int n_threads_eff = effective_threads(n_threads);

  arma::vec observed(p);
  arma::vec p_value(p);
  arma::mat sim_k(nsim, p);

  for (int j = 0; j < p; ++j) {
    const arma::vec x = X.col(j);
    observed(j) = k_stat_chol_one(x, cholC, inv_one,
                                  norm_const, sum_invC);
    double ge_count = 0.0;

#ifdef _OPENMP
    if (n_threads_eff > 1) {
#pragma omp parallel num_threads(n_threads_eff) reduction(+:ge_count)
      {
        arma::vec permuted(n);
#pragma omp for schedule(static)
        for (int i = 0; i < nsim; ++i) {
          for (int r = 0; r < n; ++r) {
            permuted(r) = x(perms(i, r) - 1);
          }

          const double k = k_stat_chol_one(permuted, cholC, inv_one,
                                           norm_const, sum_invC);
          sim_k(i, j) = k;
          if (k >= observed(j)) {
            ge_count += 1.0;
          }
        }
      }
    } else
#endif
    {
      arma::vec permuted(n);
      for (int i = 0; i < nsim; ++i) {
        for (int r = 0; r < n; ++r) {
          permuted(r) = x(perms(i, r) - 1);
        }

        const double k = k_stat_chol_one(permuted, cholC, inv_one,
                                         norm_const, sum_invC);
        sim_k(i, j) = k;
        if (k >= observed(j)) {
          ge_count += 1.0;
        }
      }
    }

    p_value(j) = ge_count / static_cast<double>(nsim);
  }

  return Rcpp::List::create(
    Rcpp::Named("K") = observed,
    Rcpp::Named("P") = p_value,
    Rcpp::Named("sim_K") = sim_k
  );
}

// [[Rcpp::export]]
Rcpp::List fast_k_chol_permutation_p_cpp(
    const arma::mat& X,
    const arma::mat& cholC,
    const double traceC,
    const Rcpp::IntegerMatrix& permutations,
    const int n_threads = 1) {
  if (X.n_rows != cholC.n_rows || cholC.n_rows != cholC.n_cols) {
    Rcpp::stop("X rows must match the square cholC matrix.");
  }
  if (permutations.ncol() != static_cast<int>(X.n_rows)) {
    Rcpp::stop("permutations must have one column per species.");
  }

  const int nsim = permutations.nrow();
  const int n = permutations.ncol();
  const int p = static_cast<int>(X.n_cols);
  const arma::imat perms = Rcpp::as<arma::imat>(permutations);
  for (int i = 0; i < nsim; ++i) {
    for (int r = 0; r < n; ++r) {
      const int idx = perms(i, r);
      if (idx < 1 || idx > n) {
        Rcpp::stop("permutations must contain 1-based indices in 1:n.");
      }
    }
  }

  const arma::vec ones(cholC.n_rows, arma::fill::ones);
  const arma::vec inv_one = chol_solve_upper(cholC, ones);
  const double sum_invC = arma::accu(inv_one);
  const double norm_const =
    (traceC - static_cast<double>(n) / sum_invC) /
    static_cast<double>(n - 1);
  const int n_threads_eff = effective_threads(n_threads);

  arma::vec observed(p);
  arma::vec p_value(p);

  for (int j = 0; j < p; ++j) {
    const arma::vec x = X.col(j);
    observed(j) = k_stat_chol_one(x, cholC, inv_one,
                                  norm_const, sum_invC);
    double ge_count = 0.0;

#ifdef _OPENMP
    if (n_threads_eff > 1) {
#pragma omp parallel num_threads(n_threads_eff) reduction(+:ge_count)
      {
        arma::vec permuted(n);
#pragma omp for schedule(static)
        for (int i = 0; i < nsim; ++i) {
          for (int r = 0; r < n; ++r) {
            permuted(r) = x(perms(i, r) - 1);
          }

          const double k = k_stat_chol_one(permuted, cholC, inv_one,
                                           norm_const, sum_invC);
          if (k >= observed(j)) {
            ge_count += 1.0;
          }
        }
      }
    } else
#endif
    {
      arma::vec permuted(n);
      for (int i = 0; i < nsim; ++i) {
        for (int r = 0; r < n; ++r) {
          permuted(r) = x(perms(i, r) - 1);
        }

        const double k = k_stat_chol_one(permuted, cholC, inv_one,
                                         norm_const, sum_invC);
        if (k >= observed(j)) {
          ge_count += 1.0;
        }
      }
    }

    p_value(j) = ge_count / static_cast<double>(nsim);
  }

  return Rcpp::List::create(
    Rcpp::Named("K") = observed,
    Rcpp::Named("P") = p_value
  );
}

// [[Rcpp::export]]
double lambda_loglik_cpp(const double lambda,
                         const arma::mat& C,
                         const arma::vec& y) {
  if (C.n_rows != C.n_cols || C.n_rows != y.n_elem) {
    Rcpp::stop("C must be square and match y length.");
  }
  return lambda_loglik_internal(lambda, C, y);
}

// [[Rcpp::export]]
arma::vec phylo_d_sums_cpp(const arma::mat& states,
                           const Rcpp::IntegerMatrix& edge,
                           const arma::vec& edge_length,
                           const int n_tip,
                           const int n_threads = 1) {
  const int n_edge = edge.nrow();
  const int n_col = static_cast<int>(states.n_cols);
  if (states.n_rows != static_cast<arma::uword>(n_tip)) {
    Rcpp::stop("states must have one row per tip.");
  }
  if (edge.ncol() != 2 || edge_length.n_elem != static_cast<arma::uword>(n_edge)) {
    Rcpp::stop("edge must be a two-column matrix matching edge_length.");
  }

  const int n_threads_eff = effective_threads(n_threads);
  if (n_threads_eff > 1 && !phylo_d_has_polytomy(edge, n_tip)) {
    arma::vec sums(n_col, arma::fill::zeros);
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads_eff) schedule(static)
#endif
    for (int j = 0; j < n_col; ++j) {
      sums(j) = phylo_d_sum_one_binary(
        states.col(j), edge, edge_length, n_tip
      );
    }
    return sums;
  }

  int max_node = n_tip;
  for (int i = 0; i < n_edge; ++i) {
    max_node = std::max(max_node, edge(i, 0));
    max_node = std::max(max_node, edge(i, 1));
  }
  const int root = n_tip + 1;

  arma::mat node_value(max_node, n_col, arma::fill::zeros);
  node_value.rows(0, n_tip - 1) = states;
  arma::vec sums(n_col, arma::fill::zeros);
  arma::vec edge_len = edge_length;

  std::vector<int> child_to_edge(max_node + 1, -1);
  for (int i = 0; i < n_edge; ++i) {
    child_to_edge[edge(i, 1)] = i;
  }

  int row = 0;
  while (row < n_edge) {
    const int parent = edge(row, 0);
    std::vector<int> children;
    std::vector<double> bl;
    while (row < n_edge && edge(row, 0) == parent) {
      children.push_back(edge(row, 1));
      bl.push_back(edge_len(row));
      ++row;
    }

    const int n_child = static_cast<int>(children.size());
    if (n_child < 2) {
      if (n_child == 1) {
        node_value.row(parent - 1) = node_value.row(children[0] - 1);
      }
      continue;
    }

    arma::rowvec curr_nv(n_col, arma::fill::zeros);
    arma::rowvec curr_contr(n_col, arma::fill::zeros);
    double curr_bl_adj = 0.0;

    if (n_child == 2) {
      const double inv0 = 1.0 / bl[0];
      const double inv1 = 1.0 / bl[1];
      const double denom = inv0 + inv1;
      curr_nv = (node_value.row(children[0] - 1) * inv0 +
                 node_value.row(children[1] - 1) * inv1) / denom;
      curr_contr =
        arma::abs(node_value.row(children[0] - 1) - curr_nv) +
        arma::abs(node_value.row(children[1] - 1) - curr_nv);
      curr_bl_adj = (bl[0] * bl[1]) / (bl[0] + bl[1]);
    } else {
      const int ref_col = (n_col > 1) ? 1 : 0; // caper uses ref.var = "V1".
      std::vector<double> rv(n_child);
      double mean_rv = 0.0;
      for (int i = 0; i < n_child; ++i) {
        rv[i] = node_value(children[i] - 1, ref_col);
        mean_rv += rv[i];
      }
      mean_rv /= static_cast<double>(n_child);

      double var_rv = 0.0;
      for (int i = 0; i < n_child; ++i) {
        const double diff = rv[i] - mean_rv;
        var_rv += diff * diff;
      }
      if (n_child > 1) {
        var_rv /= static_cast<double>(n_child - 1);
      }

      std::vector<int> group(n_child, 0);
      if (var_rv < std::numeric_limits<double>::epsilon()) {
        group[0] = 1;
      } else {
        for (int i = 0; i < n_child; ++i) {
          group[i] = (rv[i] > mean_rv) ? 1 : 0;
        }
      }

      arma::rowvec sum0(n_col, arma::fill::zeros);
      arma::rowvec sum1(n_col, arma::fill::zeros);
      double w0 = 0.0;
      double w1 = 0.0;
      for (int i = 0; i < n_child; ++i) {
        const double inv = 1.0 / bl[i];
        if (group[i]) {
          sum1 += node_value.row(children[i] - 1) * inv;
          w1 += inv;
        } else {
          sum0 += node_value.row(children[i] - 1) * inv;
          w0 += inv;
        }
      }
      if (w0 <= 0.0 || w1 <= 0.0) {
        Rcpp::stop("Invalid polytomy split in phylo.d calculation.");
      }

      const arma::rowvec sub0 = sum0 / w0;
      const arma::rowvec sub1 = sum1 / w1;
      curr_nv = (sub0 * w0 + sub1 * w1) / (w0 + w1);
      for (int i = 0; i < n_child; ++i) {
        curr_contr += arma::abs(node_value.row(children[i] - 1) - curr_nv);
      }
      curr_bl_adj = 1.0 / (w0 + w1);
    }

    node_value.row(parent - 1) = curr_nv;
    sums += curr_contr.t();

    if (parent != root) {
      const int parent_edge = child_to_edge[parent];
      if (parent_edge >= 0) {
        edge_len(parent_edge) += curr_bl_adj;
      }
    }
  }

  return sums;
}

double quantile_type7_sorted(const std::vector<double>& values,
                             const double probability) {
  const int n = static_cast<int>(values.size());
  const double h = 1.0 + static_cast<double>(n - 1) * probability;
  const int lo = static_cast<int>(std::floor(h));
  const double gamma = h - static_cast<double>(lo);

  if (lo <= 1) {
    return values[0];
  }
  if (lo >= n) {
    return values[n - 1];
  }
  return (1.0 - gamma) * values[lo - 1] + gamma * values[lo];
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix brownian_threshold_cpp(const arma::mat& samples,
                                           const double prop_state1) {
  const int n = static_cast<int>(samples.n_rows);
  const int p = static_cast<int>(samples.n_cols);
  Rcpp::IntegerMatrix out(n, p);

  if (prop_state1 < 0.0 || prop_state1 > 1.0) {
    Rcpp::stop("prop_state1 must be in [0, 1].");
  }

  std::vector<double> values(n);
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < n; ++i) {
      values[i] = samples(i, j);
    }
    std::sort(values.begin(), values.end());
    const double threshold = quantile_type7_sorted(values, prop_state1);

    for (int i = 0; i < n; ++i) {
      out(i, j) = samples(i, j) < threshold ? 1 : 0;
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix brownian_tree_threshold_cpp(
    const Rcpp::IntegerMatrix& edge,
    const arma::vec& edge_length,
    const int n_tip,
    const int nsim,
    const double prop_state1) {
  const int n_edge = edge.nrow();
  if (edge.ncol() != 2 ||
      edge_length.n_elem != static_cast<arma::uword>(n_edge)) {
    Rcpp::stop("edge must be a two-column matrix matching edge_length.");
  }
  if (n_tip < 1 || nsim < 1) {
    Rcpp::stop("n_tip and nsim must be positive.");
  }
  if (prop_state1 < 0.0 || prop_state1 > 1.0) {
    Rcpp::stop("prop_state1 must be in [0, 1].");
  }

  int max_node = n_tip;
  for (int i = 0; i < n_edge; ++i) {
    max_node = std::max(max_node, edge(i, 0));
    max_node = std::max(max_node, edge(i, 1));
    if (edge_length(i) < 0.0 || !std::isfinite(edge_length(i))) {
      Rcpp::stop("edge_length must be finite and non-negative.");
    }
  }

  std::vector< std::vector<int> > children(max_node + 1);
  std::vector< std::vector<double> > branch_lengths(max_node + 1);
  std::vector<int> is_child(max_node + 1, 0);
  for (int i = 0; i < n_edge; ++i) {
    const int parent = edge(i, 0);
    const int child = edge(i, 1);
    if (parent < 1 || child < 1 || parent > max_node || child > max_node) {
      Rcpp::stop("edge contains node numbers outside the expected range.");
    }
    children[parent].push_back(child);
    branch_lengths[parent].push_back(edge_length(i));
    is_child[child] = 1;
  }

  int root = n_tip + 1;
  if (root > max_node || children[root].empty() || is_child[root]) {
    root = -1;
    for (int node = n_tip + 1; node <= max_node; ++node) {
      if (!children[node].empty() && !is_child[node]) {
        root = node;
        break;
      }
    }
  }
  if (root < 1) {
    Rcpp::stop("Could not identify the root node.");
  }

  Rcpp::IntegerMatrix out(n_tip, nsim);
  std::vector<double> node_value(max_node + 1, 0.0);
  std::vector<double> tip_values(n_tip);
  std::vector<double> sorted_values(n_tip);
  std::vector<int> stack;
  stack.reserve(max_node);

  for (int sim = 0; sim < nsim; ++sim) {
    std::fill(node_value.begin(), node_value.end(), 0.0);
    stack.clear();
    stack.push_back(root);

    while (!stack.empty()) {
      const int parent = stack.back();
      stack.pop_back();
      for (std::size_t k = 0; k < children[parent].size(); ++k) {
        const int child = children[parent][k];
        const double bl = branch_lengths[parent][k];
        node_value[child] =
          node_value[parent] + std::sqrt(bl) * R::rnorm(0.0, 1.0);
        if (!children[child].empty()) {
          stack.push_back(child);
        }
      }
    }

    for (int tip = 0; tip < n_tip; ++tip) {
      tip_values[tip] = node_value[tip + 1];
      sorted_values[tip] = tip_values[tip];
    }
    std::sort(sorted_values.begin(), sorted_values.end());
    const double threshold =
      quantile_type7_sorted(sorted_values, prop_state1);

    for (int tip = 0; tip < n_tip; ++tip) {
      out(tip, sim) = tip_values[tip] < threshold ? 1 : 0;
    }
  }

  return out;
}

// [[Rcpp::export]]
double fast_ace_discrete_deviance_cpp(const Rcpp::IntegerMatrix& edge,
                                      const arma::vec& edge_length,
                                      const Rcpp::IntegerVector& tip_state,
                                      const Rcpp::IntegerMatrix& rate_index,
                                      const arma::vec& par) {
  if (arma::any(par < 0.0) || !par.is_finite()) {
    return 1e50;
  }
  const double dev = ace_discrete_pruning_internal(
    edge, edge_length, tip_state, rate_index, par
  );
  if (!std::isfinite(dev)) return 1e50;
  return dev;
}

// [[Rcpp::export]]
Rcpp::List fast_ace_discrete_liks_cpp(const Rcpp::IntegerMatrix& edge,
                                      const arma::vec& edge_length,
                                      const Rcpp::IntegerVector& tip_state,
                                      const Rcpp::IntegerMatrix& rate_index,
                                      const arma::vec& par,
                                      const bool marginal = false) {
  if (arma::any(par < 0.0) || !par.is_finite()) {
    Rcpp::stop("par must be finite and non-negative.");
  }
  arma::mat lik_anc;
  const double dev = ace_discrete_pruning_internal(
    edge, edge_length, tip_state, rate_index, par, &lik_anc, marginal
  );
  if (!std::isfinite(dev)) {
    Rcpp::stop("ancestral likelihood calculation failed.");
  }
  return Rcpp::List::create(
    Rcpp::Named("deviance") = dev,
    Rcpp::Named("loglik") = -0.5 * dev,
    Rcpp::Named("lik.anc") = lik_anc
  );
}

// [[Rcpp::export]]
arma::vec delta_entropy_cpp(const arma::mat& probabilities,
                            const int entropy_type) {
  if (!probabilities.is_finite()) {
    Rcpp::stop("probabilities must be finite.");
  }
  return delta_entropy_internal(probabilities, entropy_type);
}

// [[Rcpp::export]]
Rcpp::List delta_mcmc_cpp(const arma::mat& probabilities,
                          const double lambda0,
                          const double proposal_sd,
                          const int sim,
                          const int thin,
                          const int burn,
                          const int entropy_type) {
  if (sim < 1 || thin < 1 || burn < 1 || burn > sim) {
    Rcpp::stop("sim, thin, and burn must satisfy 1 <= burn <= sim.");
  }
  if (lambda0 <= 0.0 || proposal_sd <= 0.0) {
    Rcpp::stop("lambda0 and proposal_sd must be positive.");
  }

  arma::vec x = delta_entropy_internal(probabilities, entropy_type);
  if (arma::any(x <= 0.0) || arma::any(x >= 1.0) || !x.is_finite()) {
    Rcpp::stop("Delta entropies must be finite and strictly inside (0, 1).");
  }

  const double sum_log_x = arma::accu(arma::log(x));
  const double sum_log_1_minus_x = arma::accu(arma::log(1.0 - x));

  const arma::vec chain1 = delta_chain_sums(
    R::rexp(1.0), R::rexp(1.0), x, lambda0, proposal_sd, sim, thin, burn,
    sum_log_x, sum_log_1_minus_x
  );
  const arma::vec chain2 = delta_chain_sums(
    R::rexp(1.0), R::rexp(1.0), x, lambda0, proposal_sd, sim, thin, burn,
    sum_log_x, sum_log_1_minus_x
  );

  const double saved = chain1(2) + chain2(2);
  if (saved <= 0.0) {
    Rcpp::stop("No MCMC samples were saved.");
  }
  const double alpha_mean = (chain1(0) + chain2(0)) / saved;
  const double beta_mean = (chain1(1) + chain2(1)) / saved;
  const double delta = beta_mean / alpha_mean;

  return Rcpp::List::create(
    Rcpp::Named("delta") = delta,
    Rcpp::Named("alpha_mean") = alpha_mean,
    Rcpp::Named("beta_mean") = beta_mean,
    Rcpp::Named("n_saved") = saved,
    Rcpp::Named("entropy") = x
  );
}
