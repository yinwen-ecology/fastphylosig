// Fast numerical kernels for phytools-compatible phylosig calculations.
//
// K and lambda keep the same statistical targets as phytools::phylosig().
// Dense inverse matrices are avoided where possible: Cholesky triangular solves
// evaluate the same C^{-1} GLS means and quadratic forms more cheaply.

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cfenv>
#include <limits>
#include <stdexcept>
#include <vector>
#include <unordered_map>

#include "numeric_utils.h"
#include "ace_exponent_cache.h"
#include "d_kernels.h"
#include "bounded_optimize.h"
#include "heterogeneity_bounds.h"

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

arma::vec k_stats_chol_matrix(const arma::mat& X,
                              const arma::mat& cholC,
                              const arma::vec& inv_one,
                              const double norm_const,
                              const double sum_invC) {
  const arma::rowvec means = (inv_one.t() * X) / sum_invC;
  arma::mat centered = X;
  centered.each_row() -= means;
  const arma::mat z = arma::solve(
    arma::trimatl(cholC.t()), centered
  );
  const arma::vec numerator = arma::sum(arma::square(centered), 0).t();
  const arma::vec denominator = arma::sum(arma::square(z), 0).t();
  arma::vec out = numerator;
  out /= denominator;
  out /= norm_const;
  return out;
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

double lambda_loglik_spectral_internal(
    const double lambda,
    const arma::vec& eigval,
    const arma::mat& eigvec,
    const arma::vec& inv_sqrt_diag,
    const double log_diag,
    const arma::vec& y) {
  const arma::uword n = eigval.n_elem;
  if (eigvec.n_rows != n || eigvec.n_cols != n ||
      inv_sqrt_diag.n_elem != n || y.n_elem != n) {
    Rcpp::stop("spectral lambda inputs have incompatible dimensions.");
  }
  const arma::vec m = (1.0 - lambda) + lambda * eigval;
  if (!m.is_finite() || arma::any(m <= 0.0)) {
    return -std::numeric_limits<double>::infinity();
  }
  const arma::vec inv_m = 1.0 / m;
  const arma::vec w = inv_sqrt_diag;
  const arma::vec z = inv_sqrt_diag % y;
  const arma::vec wproj = eigvec.t() * w;
  const arma::vec zproj = eigvec.t() * z;
  const double sum_inv = arma::dot(arma::square(wproj), inv_m);
  if (!std::isfinite(sum_inv) || sum_inv <= 0.0) {
    return -std::numeric_limits<double>::infinity();
  }
  const double a = arma::dot(zproj % wproj, inv_m) / sum_inv;
  const arma::vec centered = zproj - a * wproj;
  const double quad = arma::dot(arma::square(centered), inv_m);
  const double sig2 = quad / static_cast<double>(n);
  if (!std::isfinite(sig2) || sig2 <= 0.0) {
    return -std::numeric_limits<double>::infinity();
  }
  const double logdet = log_diag + arma::accu(arma::log(m));
  if (!std::isfinite(logdet)) {
    return -std::numeric_limits<double>::infinity();
  }
  const double log_two_pi = std::log(6.283185307179586476925286766559);
  const double logdet_sig2 = static_cast<double>(n) * std::log(sig2) + logdet;
  return -0.5 * quad / sig2 -
    0.5 * static_cast<double>(n) * log_two_pi -
    0.5 * logdet_sig2;
}

double lambda_loglik_spectral_projected_internal(
    const double lambda,
    const arma::vec& eigval,
    const arma::vec& wproj,
    const arma::vec& zproj,
    const double log_diag) {
  const arma::uword n = eigval.n_elem;
  if (wproj.n_elem != n || zproj.n_elem != n) {
    Rcpp::stop("projected spectral lambda inputs have incompatible dimensions.");
  }
  const arma::vec m = (1.0 - lambda) + lambda * eigval;
  if (!m.is_finite() || arma::any(m <= 0.0)) {
    return -std::numeric_limits<double>::infinity();
  }
  const arma::vec inv_m = 1.0 / m;
  const double sum_inv = arma::dot(arma::square(wproj), inv_m);
  if (!std::isfinite(sum_inv) || sum_inv <= 0.0) {
    return -std::numeric_limits<double>::infinity();
  }
  const double a = arma::dot(zproj % wproj, inv_m) / sum_inv;
  const arma::vec centered = zproj - a * wproj;
  const double quad = arma::dot(arma::square(centered), inv_m);
  const double sig2 = quad / static_cast<double>(n);
  if (!std::isfinite(sig2) || sig2 <= 0.0) {
    return -std::numeric_limits<double>::infinity();
  }
  const double logdet = log_diag + arma::accu(arma::log(m));
  if (!std::isfinite(logdet)) {
    return -std::numeric_limits<double>::infinity();
  }
  const double log_two_pi = std::log(6.283185307179586476925286766559);
  return -0.5 * quad / sig2 -
    0.5 * static_cast<double>(n) * log_two_pi -
    0.5 * (static_cast<double>(n) * std::log(sig2) + logdet);
}

arma::mat lambda_loglik_spectral_projected_batch_internal(
    const arma::vec& lambda,
    const arma::vec& eigval,
    const arma::vec& wproj,
    const arma::mat& zproj,
    const double log_diag,
    const int n_threads) {
  const arma::uword n = eigval.n_elem;
  if (wproj.n_elem != n || zproj.n_rows != n) {
    Rcpp::stop("projected spectral lambda batch inputs have incompatible dimensions.");
  }
  arma::mat out(lambda.n_elem, zproj.n_cols);
  const int n_threads_eff = effective_threads(n_threads);

#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads_eff) schedule(static)
#endif
  for (int li = 0; li < static_cast<int>(lambda.n_elem); ++li) {
    const double lam = lambda(li);
    const arma::vec m = (1.0 - lam) + lam * eigval;
    if (!m.is_finite() || arma::any(m <= 0.0)) {
      out.row(li).fill(-std::numeric_limits<double>::infinity());
      continue;
    }
    const arma::vec inv_m = 1.0 / m;
    const double logdet = log_diag + arma::accu(arma::log(m));
    const double sum_inv = arma::dot(arma::square(wproj), inv_m);
    if (!std::isfinite(sum_inv) || sum_inv <= 0.0 || !std::isfinite(logdet)) {
      out.row(li).fill(-std::numeric_limits<double>::infinity());
      continue;
    }
    for (arma::uword j = 0; j < zproj.n_cols; ++j) {
      const arma::vec centered = zproj.col(j) -
        (arma::dot(zproj.col(j) % wproj, inv_m) / sum_inv) * wproj;
      const double quad = arma::dot(arma::square(centered), inv_m);
      const double sig2 = quad / static_cast<double>(n);
      if (!std::isfinite(sig2) || sig2 <= 0.0) {
        out(li, j) = -std::numeric_limits<double>::infinity();
      } else {
        const double log_two_pi =
          std::log(6.283185307179586476925286766559);
        out(li, j) = -0.5 * quad / sig2 -
          0.5 * static_cast<double>(n) * log_two_pi -
          0.5 * (static_cast<double>(n) * std::log(sig2) + logdet);
      }
    }
  }
  return out;
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

// Delta reuses per-edge vector storage while keeping the same matrix products
// and element order as the ordinary ACE transition path.
struct AceTransitionWorkspace {
  arma::cx_vec child;
  arma::cx_vec transformed;
  arma::cx_vec exp_lambda;
  arma::cx_vec gamma_transformed;

  void reset(const arma::uword k) {
    child.set_size(k);
    transformed.set_size(k);
    exp_lambda.set_size(k);
    gamma_transformed.set_size(k);
  }
};

// Owned by one fit or analysis cache; never stored in a serializable prepared
// tree or shared between processes. Internal rows are overwritten before
// they are read, so only the fixed tip rows need initialization.
struct AcePruningWorkspace {
  arma::mat liks;
  arma::cx_vec lambda;
  arma::cx_mat gamma, inv_gamma, tip_transformed;
  std::vector<bool> transformed_ready;
  std::vector<int> edge_exp_slot;
  arma::cx_mat exponents;
  std::vector<bool> exp_ready;
  // A bounded batch is used for sufficiently large uncached-length workloads.
  // Exponents are independent; all pruning/products/log sums remain serial.
  arma::cx_mat exponent_batch;
  std::vector<arma::cx_vec> exponent_workers;
  // Only mixed-length batches need a compact worklist. Their uncached slots
  // encode -(work index + 1); distinct-only fits use edge rows directly.
  std::vector<int> exponent_edges;
  int batch_begin = 0, batch_end = 0;
  AceTransitionWorkspace transition;
  AceExponentCache eigen_exponents;
  arma::vec left, right, product;
  bool initialized = false;
  arma::uword er_recoveries = 0;
  arma::uword stable_recoveries = 0, probability_recoveries = 0;
};

struct AceFitWorkspace {
  Rcpp::IntegerMatrix edge, rate_index;
  Rcpp::IntegerVector tip_state;
  arma::vec edge_length;
  AcePruningWorkspace scratch;
  AceFitWorkspace(const Rcpp::IntegerMatrix& e, const arma::vec& lengths,
                  const Rcpp::IntegerVector& states, const Rcpp::IntegerMatrix& rates,
                  const int requested_threads)
    : edge(Rcpp::clone(e)), rate_index(Rcpp::clone(rates)),
      tip_state(Rcpp::clone(states)), edge_length(lengths) {
    // Exact repeated branch lengths only. Cache exponent vectors rather than
    // transition matrices, retaining the original multiplication/rounding.
    // Bound storage as in the Mk simulation cache. The slots and uncached
    // worklist are fixed for this immutable fit input.
    std::unordered_map<double,int> count, slot;
    count.reserve(lengths.n_elem);
    for (double length : lengths) if (std::isfinite(length)) ++count[length];
    scratch.edge_exp_slot.assign(lengths.n_elem,-1);
    const int k = rates.nrow();
    // Exponent vectors use O(k), not O(k^2), bytes per repeated length. Let the
    // actual byte budget determine capacity instead of an unrelated 4096 cap.
    const int limit = k>0 ? static_cast<int>(64u*1024u*1024u/(sizeof(std::complex<double>)*k)) : 0;
    for (arma::uword i=0; i<lengths.n_elem; ++i) {
      if (!std::isfinite(lengths[i]) || count[lengths[i]]<2) continue;
      auto found = slot.find(lengths[i]);
      if (found==slot.end() && static_cast<int>(slot.size())<limit)
        found=slot.emplace(lengths[i],static_cast<int>(slot.size())).first;
      if (found!=slot.end()) scratch.edge_exp_slot[i]=found->second;
    }
    scratch.exponents.set_size(k,slot.size());
    scratch.exp_ready.resize(slot.size());
    const arma::uword uncached = std::count(
      scratch.edge_exp_slot.begin(),scratch.edge_exp_slot.end(),-1);
    if (slot.empty()) scratch.edge_exp_slot.clear();
    int threads = 1;
#ifdef _OPENMP
    if (!omp_in_parallel()) threads = std::max(1,std::min(
      requested_threads == 0 ? 2 : requested_threads, omp_get_max_threads()));
#endif
    // Count actual uncached work, not all edges: repeated lengths must not
    // trigger small/empty parallel batches. No global thread settings change.
    const arma::uword cap = k>0 ? (2u*1024u*1024u)/(sizeof(std::complex<double>)*k) : 0;
    // Work items are independent: an odd width is valid and avoids a second
    // OpenMP region containing just one leftover branch in mixed-length fits.
    const arma::uword width = std::min(uncached,cap);
    if (threads>1 && width*static_cast<arma::uword>(k)>=8192) {
      if (!slot.empty()) {
        scratch.exponent_edges.reserve(uncached);
        for (arma::uword row=0; row<lengths.n_elem; ++row) {
          if (scratch.edge_exp_slot[row]>=0) continue;
          scratch.edge_exp_slot[row] = -1-static_cast<int>(scratch.exponent_edges.size());
          scratch.exponent_edges.push_back(static_cast<int>(row));
        }
      }
      scratch.exponent_batch.set_size(k,width);
      scratch.exponent_workers.resize(threads);
      for (auto& worker : scratch.exponent_workers) worker.set_size(k);
    }
    // Small distinct-length ARD fits are below the exponent batch threshold.
    // Their exact eigenvalues can recur across different Q optimizer trials.
    if(k==3 && rates.ncol()==3 && lengths.n_elem>=64 && lengths.n_elem<=8191u/3u &&
        scratch.edge_exp_slot.empty() && scratch.exponent_workers.empty() &&
        std::any_of(rates.begin(),rates.end(),[](int value){return value>1;}) &&
        lengths.is_finite() && arma::all(lengths>0))
      scratch.eigen_exponents.configure(lengths.n_elem);
  }
};

bool ace_is_er(const Rcpp::IntegerMatrix& index) {
  const int k = index.nrow();
  if (k < 2 || index.ncol() != k) return false;
  for (int i=0; i<k; ++i) for (int j=0; j<k; ++j)
    if (index(i,j) != (i==j ? 0 : 1)) return false;
  return true;
}

void ace_er_transition_terms(const int k, const double q, const double length,
                             double& decay, double& mix) {
  // Zero branches stay identities even if k*q overflows at a finite rate.
  const double exponent = length==0 ? 0 : -(k*q)*length;
  decay=std::exp(exponent); mix=-std::expm1(exponent)/k;
}

// Shared by the optional ER reference objective and the numerical fallback.
// Exact-length slots retain two real coefficients in existing complex storage.
void ace_er_transition_message(const int child, const int edge, const int k,
                                const double q, const double length,
                                AcePruningWorkspace& s, arma::vec& out) {
  const int slot=s.edge_exp_slot.empty() ? -1 : s.edge_exp_slot[edge];
  double decay, mix;
  if (slot>=0 && s.exp_ready[slot]) {
    decay=s.exponents(0,slot).real(); mix=s.exponents(1,slot).real();
  } else {
    ace_er_transition_terms(k,q,length,decay,mix);
    if (slot>=0) {
      s.exponents(0,slot)=decay; s.exponents(1,slot)=mix; s.exp_ready[slot]=true;
    }
  }
  double total=0;
  for (int j=0; j<k; ++j) total+=s.liks(child,j);
  for (int j=0; j<k; ++j) out[j]=decay*s.liks(child,j)+mix*total;
}

arma::mat ace_er_transition_matrix(const int k, const double q, const double length) {
  double decay, mix;
  ace_er_transition_terms(k,q,length,decay,mix);
  arma::mat P(k,k); P.fill(mix); P.diag()+=decay;
  return P;
}

// Reference-only ER objective. Reuses ACE topology, tip rows and exact-length
// slots, with P(t)v = exp(-k*q*t)v + (1-exp(-k*q*t))*sum(v)/k.
void ace_er_reference_check(AceFitWorkspace& w) {
  const int nt = w.tip_state.size(), k = w.rate_index.nrow();
  if (nt < 2 || k < 2 || w.rate_index.ncol() != k || w.edge.ncol() != 2 ||
      w.edge.nrow() != 2*(nt-1) || w.edge_length.n_elem != static_cast<arma::uword>(w.edge.nrow()))
    Rcpp::stop("invalid prepared ER tree");
  for (int i=0; i<k; ++i) for (int j=0; j<k; ++j)
    if (w.rate_index(i,j) != (i==j ? 0 : 1))
      Rcpp::stop("ER reference fitting requires one common off-diagonal rate");
  for (int state : w.tip_state)
    if (!Rcpp::IntegerVector::is_na(state) && (state<1 || state>k))
      Rcpp::stop("tip_state values must be 1-based state indices");
  for (int e=0; e<w.edge.nrow(); ++e) {
    if (!std::isfinite(w.edge_length[e]) || w.edge_length[e]<=0)
      Rcpp::stop("ER reference fitting requires positive finite branch lengths");
    if (w.edge(e,0)<=nt || w.edge(e,0)>2*nt-1 || w.edge(e,1)<1 || w.edge(e,1)>2*nt-1 ||
        (e%2 && w.edge(e,0)!=w.edge(e-1,0))) Rcpp::stop("invalid prepared ER edges");
  }
}

double ace_er_reference_deviance(AceFitWorkspace& w, double q) {
  if (std::isnan(q) || q<0) return std::numeric_limits<double>::infinity();
  const int nt = w.tip_state.size(), k = w.rate_index.nrow();
  AcePruningWorkspace& s = w.scratch;
  if (!s.initialized) {
    s.liks.zeros(2*nt-1,k);
    for (int i=0; i<nt; ++i) {
      if (Rcpp::IntegerVector::is_na(w.tip_state[i])) s.liks.row(i).ones();
      else s.liks(i,w.tip_state[i]-1)=1;
    }
    s.initialized=true;
  }
  std::fill(s.exp_ready.begin(),s.exp_ready.end(),false);
  s.left.set_size(k); s.right.set_size(k); s.product.set_size(k);
  const auto message = [&](int edge, arma::vec& out) {
    ace_er_transition_message(w.edge(edge,1)-1,edge,k,q,w.edge_length[edge],s,out);
  };
  double log_comp=0;
  for (int e=0; e<w.edge.nrow(); e+=2) {
    message(e,s.left); message(e+1,s.right);
    s.product=s.left%s.right;
    const double scale=arma::accu(s.product);
    if (!(scale>0) || !std::isfinite(scale)) return std::numeric_limits<double>::infinity();
    s.liks.row(w.edge(e,0)-1)=(s.product/scale).t();
    log_comp+=std::log(scale);
  }
  // Same root-sum constant as ACE; subtract log(k) for a uniform-root loglik.
  return -2*log_comp;
}

// The Armadillo expression promotes the branch scalar to complex before
// multiplication. Keep that operation (including signed zeros) exactly.
// Within the normal real range, cexp uses the same exp result;
// overflowing/underflowing, nonfinite and complex inputs retain std::exp.
inline std::complex<double> ace_exponent_value(const std::complex<double>& eigenvalue,
                                               const double length) {
  const std::complex<double> factor(length,0.0);
  const std::complex<double> z=eigenvalue*factor;
  return z.imag()==0 && z.real()>=-700 && z.real()<=700 ?
    std::complex<double>(std::exp(z.real()),z.imag()) : std::exp(z);
}

void ace_exponents_reuse(const arma::cx_vec& lambda,const double length,arma::cx_vec& out) {
  out.set_size(lambda.n_elem);
  for(arma::uword j=0;j<lambda.n_elem;++j)out[j]=ace_exponent_value(lambda[j],length);
}

// The reuse workspace is already sized and distinct from the operands.
// Delegate to the same Armadillo gemv dispatch used by glue_times, retaining
// its BLAS/ATLAS/emulation and scalar-argument choices exactly.
inline void ace_gemv_reuse(arma::cx_vec& out,const arma::cx_mat& A,
                           const arma::cx_vec& x) {
  arma::gemv<false,false,false>::apply(out.memptr(),A,x.memptr());
}

void ace_transition_multiply_reuse(const arma::cx_vec& lambda,
                                    const arma::cx_mat& gamma,
                                    const arma::cx_mat& inv_gamma,
                                    const arma::mat& liks, const int child,
                                    const double edge_length,
                                    AceTransitionWorkspace& ws,
                                    arma::vec& out,
                                    const arma::cx_mat* tip_transformed = NULL,
                                    const int tip_column = 0,
                                    const arma::cx_mat* exponents = NULL,
                                    const int exp_column = 0,
                                    const AceExponentCache* eigen_exponents = NULL,
                                    const int edge_row = 0) {
  if (tip_transformed != NULL) {
    std::copy_n(tip_transformed->colptr(tip_column),ws.transformed.n_elem,ws.transformed.memptr());
  } else {
    for (arma::uword i = 0; i < liks.n_cols; ++i) {
      ws.child(i) = std::complex<double>(liks(child,i), 0.0);
    }
    ace_gemv_reuse(ws.transformed,inv_gamma,ws.child);
  }
  if(eigen_exponents != NULL) {
    // Same component-wise complex products and order as %=, without copying
    // cached values through exp_lambda. This cache is only enabled for k=3.
    for(arma::uword j=0;j<3;++j)ws.transformed[j] *= eigen_exponents->value(j,edge_row);
  } else {
    if(exponents != NULL)std::copy_n(exponents->colptr(exp_column),ws.exp_lambda.n_elem,ws.exp_lambda.memptr());
    else ace_exponents_reuse(lambda,edge_length,ws.exp_lambda);
    ws.transformed %= ws.exp_lambda;
  }
  ace_gemv_reuse(ws.gamma_transformed,gamma,ws.transformed);
  for(arma::uword i=0;i<out.n_elem;++i)out[i]=std::real(ws.gamma_transformed[i]);
}

arma::vec ace_transition_multiply(const arma::cx_vec& lambda,
                                  const arma::cx_mat& gamma,
                                  const arma::cx_mat& inv_gamma,
                                  const arma::vec& child_lik,
                                  const double edge_length) {
  arma::cx_vec z = inv_gamma * arma::conv_to<arma::cx_vec>::from(child_lik);
  z %= arma::exp(lambda * edge_length);
  arma::vec out = arma::real(gamma * z);
  return out;
}

arma::mat ace_transition_matrix(const arma::cx_vec& lambda,
                                const arma::cx_mat& gamma,
                                const arma::cx_mat& inv_gamma,
                                const double edge_length) {
  const arma::cx_mat P =
    gamma * arma::diagmat(arma::exp(lambda * edge_length)) * inv_gamma;
  return arma::real(P);
}

// A cheap upper bound for the eigenvector condition number; no SVD or changes
// to the ordinary spectral arithmetic. Beyond sqrt(eps) amplification, use the
// same scaling-and-squaring matrix exponential as the Mk simulation fallback.
bool ace_stable_eigensystem(const arma::cx_mat& gamma, const arma::cx_mat& inverse) {
  return arma::norm(gamma,"fro") * arma::norm(inverse,"fro") <=
    1.0/std::sqrt(std::numeric_limits<double>::epsilon());
}

// A tiny positive exit rate can lose relative accuracy in the unbalanced
// spectral decomposition even when its eigenvectors pass the condition check.
// Exact absorbing rows are handled well by the existing path. Enable this
// additional guard only during an unfinished Delta ARD fit's bounded recovery.
bool ace_near_absorbing_generator(const arma::mat& Q) {
  double largest=0, smallest=std::numeric_limits<double>::infinity();
  for (arma::uword i=0; i<Q.n_rows; ++i) {
    const double exit=-Q(i,i);
    largest=std::max(largest,exit);
    if (exit>0) smallest=std::min(smallest,exit);
  }
  return smallest < std::sqrt(std::numeric_limits<double>::epsilon())*largest;
}

bool ace_valid_transition(const arma::mat& P, const double tolerance = 1e-8) {
  return P.is_finite() && P.min() >= -tolerance &&
    arma::max(arma::abs(arma::sum(P,1)-1.0)) <= tolerance;
}

arma::mat ace_stable_transition(const arma::mat& Q, const double length) {
  arma::mat P;
  const double tolerance=std::min(1e-12,8.0*Q.n_rows*std::numeric_limits<double>::epsilon());
  bool ok=false;
  try { ok=arma::expmat(P,Q*length) && ace_valid_transition(P,1e-10) &&
    P.min()>=-tolerance && P.max()<=1+tolerance; }
  catch (const std::runtime_error&) { ok=false; }
  if (!ok) {
    // An invalid optimizer trial is a rejected objective, not an exception
    // that aborts the whole deterministic search.
    P.set_size(Q.n_rows,Q.n_cols);
    P.fill(std::numeric_limits<double>::quiet_NaN());
    return P;
  }
  // Only the recovery path changes arithmetic. Remove roundoff before it can
  // propagate through products and divisions in the existing joint algorithm.
  for (arma::uword i=0; i<P.n_rows; ++i) {
    bool outside=false;
    for (arma::uword j=0; j<P.n_cols; ++j) if (P(i,j)<0 || P(i,j)>1) outside=true;
    if (outside) {
      P.row(i)=arma::clamp(P.row(i),0.0,1.0);
      P.row(i)/=arma::accu(P.row(i));
    }
  }
  return P;
}

double ace_discrete_pruning_internal(const Rcpp::IntegerMatrix& edge,
                                     const arma::vec& edge_length,
                                     const Rcpp::IntegerVector& tip_state,
                                     const Rcpp::IntegerMatrix& rate_index,
                                     const arma::vec& par,
                                     arma::mat* out_liks = NULL,
                                     const bool marginal = false,
                                     const bool reuse_buffers = false,
                                     AcePruningWorkspace* scratch = NULL,
                                     const bool force_stable = false,
                                     const bool recovery_guard = false) {
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

  AcePruningWorkspace local;
  AcePruningWorkspace& ws = scratch == NULL ? local : *scratch;
  arma::mat& liks = ws.liks;
  if (!ws.initialized) {
    liks.zeros(n_total, k);
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
    ws.initialized = true;
  }
  const arma::mat Q = ace_build_q(rate_index, par);
  arma::cx_vec& lambda = ws.lambda;
  arma::cx_mat& gamma = ws.gamma;
  arma::cx_mat& inv_gamma = ws.inv_gamma;
  bool analytic_er = false, stable = force_stable ||
    (recovery_guard && ace_near_absorbing_generator(Q));
  const auto recover_er = [&]() {
    return par.n_elem==1 && std::isfinite(par[0]) && par[0]>=0 && ace_is_er(rate_index);
  };
  try {
    if (!stable && !ace_eigen_decomp(Q, lambda, gamma, inv_gamma)) {
      analytic_er = recover_er(); stable = !analytic_er;
    }
  } catch (const std::runtime_error&) {
    analytic_er = recover_er(); stable = !analytic_er;
  }
  if (!stable && !analytic_er && !ace_stable_eigensystem(gamma,inv_gamma)) {
    analytic_er = recover_er(); stable = !analytic_er;
  }
  if (analytic_er) ++ws.er_recoveries;
  if (stable) ++ws.stable_recoveries;
  double log_comp = 0.0;
  AceTransitionWorkspace& transition_workspace = ws.transition;
  arma::vec& v_left_reuse = ws.left;
  arma::vec& v_right_reuse = ws.right;
  arma::vec& v_reuse = ws.product;
  if (reuse_buffers || analytic_er || stable) {
    transition_workspace.reset(static_cast<arma::uword>(k));
    v_left_reuse.set_size(static_cast<arma::uword>(k));
    v_right_reuse.set_size(static_cast<arma::uword>(k));
    v_reuse.set_size(static_cast<arma::uword>(k));
  }
  // Fixed tip vectors repeat throughout a trait. Transform each distinct state
  // once with the SAME matrix-vector multiplication; branch exponentials and
  // the following product retain their original evaluation order.
  arma::cx_mat& tip_transformed = ws.tip_transformed;
  std::vector<bool>& transformed_ready = ws.transformed_ready;
  if (reuse_buffers && !analytic_er && !stable) {
    tip_transformed.set_size(k,k+1);
    transformed_ready.assign(k+1,false);
  }
  std::fill(ws.exp_ready.begin(),ws.exp_ready.end(),false);
  ws.batch_begin = ws.batch_end = 0;
  const AceExponentCache* eigen_exponents = NULL;
  if(reuse_buffers && !analytic_er && !stable && ws.eigen_exponents.prepare(lambda,edge_length,
      [](const std::complex<double>& value,const double length) {
        return ace_exponent_value(value,length);
      }))
    eigen_exponents = &ws.eigen_exponents;
  auto prepare_exponent_batch = [&](const int begin) {
    const bool mixed = !ws.exponent_edges.empty();
    const int count = mixed ? static_cast<int>(ws.exponent_edges.size()) : n_edge;
    const int end = std::min(count,begin+static_cast<int>(ws.exponent_batch.n_cols));
    const int threads = static_cast<int>(ws.exponent_workers.size());
    const int rounding = std::fegetround();
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
    {
      const int tid = omp_get_thread_num();
      const int previous_rounding = std::fegetround();
      if (previous_rounding != rounding) std::fesetround(rounding);
      arma::cx_vec& values = ws.exponent_workers[tid];
#pragma omp for schedule(static)
      for (int i=begin; i<end; ++i) {
        const int row = mixed ? ws.exponent_edges[i] : i;
        // Same exponent evaluator and element order as the serial ACE path.
        ace_exponents_reuse(lambda,edge_length(row),values);
        ws.exponent_batch.col(i-begin) = values;
      }
      if (previous_rounding != rounding) std::fesetround(previous_rounding);
    }
#else
    (void)threads;
    (void)rounding;
    for (int i=begin; i<end; ++i) {
      const int row = mixed ? ws.exponent_edges[i] : i;
      arma::cx_vec& values=ws.exponent_workers[0];
      ace_exponents_reuse(lambda,edge_length(row),values);
      ws.exponent_batch.col(i-begin) = values;
    }
#endif
    ws.batch_begin = begin; ws.batch_end = end;
  };
  auto transition = [&](const int child, const int edge_row, arma::vec& out) {
    const double length = edge_length(edge_row);
    if (analytic_er) {
      ace_er_transition_message(child,edge_row,k,par[0],length,ws,out);
      return;
    }
    if (stable) {
      out = ace_stable_transition(Q,length)*liks.row(child).t();
      return;
    }
    int column = 0;
    const arma::cx_mat* cached = NULL;
    if (child < n_tip) {
      column = Rcpp::IntegerVector::is_na(tip_state[child]) ? k : tip_state[child]-1;
      if (!transformed_ready[column]) {
        for (int j=0; j<k; ++j)
          transition_workspace.child(j) = std::complex<double>(liks(child,j),0.0);
        tip_transformed.col(column) = inv_gamma * transition_workspace.child;
        transformed_ready[column] = true;
      }
      cached = &tip_transformed;
    }
    const int slot = ws.edge_exp_slot.empty() ? -1 : ws.edge_exp_slot[edge_row];
    if (slot >= 0 && !ws.exp_ready[slot]) {
      ace_exponents_reuse(lambda,length,transition_workspace.exp_lambda);
      ws.exponents.col(slot) = transition_workspace.exp_lambda;
      ws.exp_ready[slot] = true;
    }
    const arma::cx_mat* exponents = slot>=0 ? &ws.exponents : NULL;
    int exp_column = slot;
    if (slot<0 && !ws.exponent_workers.empty()) {
      const int work = ws.exponent_edges.empty() ? edge_row : -(slot+1);
      if (work>=ws.batch_end) prepare_exponent_batch(work);
      exponents = &ws.exponent_batch;
      exp_column = work-ws.batch_begin;
    }
    ace_transition_multiply_reuse(lambda,gamma,inv_gamma,liks,child,
                                  length,transition_workspace,out,cached,column,
                                  exponents,exp_column,eigen_exponents,edge_row);
  };
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

    double comp = 0.0;
    if (reuse_buffers || analytic_er || stable) {
      transition(des_left, e, v_left_reuse);
      transition(des_right, e_right, v_right_reuse);
      v_reuse = v_left_reuse % v_right_reuse;
      comp = arma::accu(v_reuse);
      if (!std::isfinite(comp) || comp <= 0.0) {
        return std::numeric_limits<double>::infinity();
      }
      liks.row(anc) = (v_reuse / comp).t();
    } else {
      const arma::vec v_left = ace_transition_multiply(
        lambda, gamma, inv_gamma, liks.row(des_left).t(), edge_length(e)
      );
      const arma::vec v_right = ace_transition_multiply(
        lambda, gamma, inv_gamma, liks.row(des_right).t(), edge_length(e_right)
      );
      const arma::vec v = v_left % v_right;
      comp = arma::accu(v);
      if (!std::isfinite(comp) || comp <= 0.0) {
        return std::numeric_limits<double>::infinity();
      }
      liks.row(anc) = (v / comp).t();
    }
    log_comp += std::log(comp);
  }

  if (out_liks != NULL) {
    arma::mat lik_anc = liks.rows(n_tip, n_total - 1);
    if (marginal) {
      *out_liks = lik_anc;
    } else {
      // Historically every row was renormalized after EVERY sibling pair.
      // Defer those independent row operations until a row is read/returned.
      // Preserve exact floating-point behavior: stop at a fixed point, or skip
      // complete two-cycles (some nearly unit rows alternate indefinitely).
      // Merely normalizing each modified row once would change Delta's RNG
      // stream when rounding changes entropy endpoint jitter.
      std::vector<int> normalized_at(n_node, 0);
      arma::rowvec previous(k), current(k);
      auto normalize_to = [&](const int row, const int step) {
        int remaining = step - normalized_at[row];
        bool have_previous = false;
        while (remaining > 0) {
          double sum = lik_anc(row, 0);
          for (int col = 1; col < k; ++col) sum += lik_anc(row, col);
          if (sum == 0.0 || sum == 1.0) break;
          current = lik_anc.row(row);
          lik_anc.row(row) /= sum;
          --remaining;
          bool fixed = true, cycle = have_previous;
          for (int col = 0; col < k; ++col) {
            fixed = fixed && lik_anc(row, col) == current[col];
            cycle = cycle && lik_anc(row, col) == previous[col];
          }
          if (fixed) break;
          if (cycle) {
            if (remaining % 2) lik_anc.row(row) = current;
            break;
          }
          previous = current;
          have_previous = true;
        }
        normalized_at[row] = step;
      };
      // Q has parent states in rows. The child contribution at each parent
      // state is P * L_child, represented here as L_child * P.t(). Retain the
      // historical ER arithmetic exactly: its P is mathematically symmetric,
      // and changing its floating-point direction can alter Delta endpoint RNG.
      const bool er_downpass = ace_is_er(rate_index);
      const auto propagate = [&](const int parent, const int child, const arma::mat& P) {
        const arma::rowvec denom = er_downpass ?
          arma::rowvec(lik_anc.row(child) * P) :
          arma::rowvec(lik_anc.row(child) * P.t());
        arma::rowvec tmp;
        if (arma::any(denom == 0.0)) {
          if (er_downpass) return false;
          // A reducible ARD Q can make a parent state impossible. Its posterior
          // and child contribution are both zero, so that state's message is
          // zero rather than 0/0. A positive posterior with zero contribution
          // remains a numerical failure; do not silently discard its mass.
          tmp.zeros(k);
          for (int j = 0; j < k; ++j) {
            if (denom[j] != 0.0) tmp[j] = lik_anc(parent, j) / denom[j];
            else if (lik_anc(parent, j) != 0.0) return false;
          }
        } else tmp = lik_anc.row(parent) / denom;
        lik_anc.row(child) = (tmp * P) % lik_anc.row(child);
        return true;
      };
      int step = 0;
      for (int e = n_edge - 2; e >= 0; e -= 2) {
        const int e_right = e + 1;
        const int anc = edge(e, 0) - n_tip - 1;
        normalize_to(anc, step);

        const int des_left = edge(e, 1) - n_tip - 1;
        if (des_left >= 0) {
          normalize_to(des_left, step);
          const arma::mat P = stable ? ace_stable_transition(Q,edge_length(e)) :
            analytic_er ? ace_er_transition_matrix(k,par[0],edge_length(e)) :
            ace_transition_matrix(lambda, gamma, inv_gamma, edge_length(e));
          if (!propagate(anc, des_left, P))
            return std::numeric_limits<double>::infinity();
        }

        const int des_right = edge(e_right, 1) - n_tip - 1;
        if (des_right >= 0) {
          normalize_to(des_right, step);
          const arma::mat P = stable ? ace_stable_transition(Q,edge_length(e_right)) :
            analytic_er ? ace_er_transition_matrix(k,par[0],edge_length(e_right)) :
            ace_transition_matrix(lambda, gamma, inv_gamma, edge_length(e_right));
          if (!propagate(anc, des_right, P))
            return std::numeric_limits<double>::infinity();
        }

        ++step;
      }
      for (int row = 0; row < n_node; ++row) normalize_to(row, step);
      *out_liks = lik_anc;
    }
    const double tolerance=std::min(1e-12,8.0*k*std::numeric_limits<double>::epsilon());
    for (arma::uword i=0; i<out_liks->n_rows; ++i) {
      double total=0; bool outside=false, invalid=false;
      for (int j=0; j<k; ++j) {
        const double p=(*out_liks)(i,j); total+=p;
        outside=outside || p<0 || p>1;
        invalid=invalid || !std::isfinite(p) || p < -tolerance || p > 1+tolerance;
      }
      invalid=invalid || std::abs(total-1)>tolerance;
      if (invalid) {
        if (!stable && !analytic_er) return ace_discrete_pruning_internal(edge,edge_length,
          tip_state,rate_index,par,out_liks,marginal,reuse_buffers,&ws,true,recovery_guard);
        return std::numeric_limits<double>::infinity();
      }
      if (outside) {
        out_liks->row(i)=arma::clamp(out_liks->row(i),0.0,1.0);
        out_liks->row(i)/=arma::accu(out_liks->row(i));
        ++ws.probability_recoveries;
      }
    }
  }

  return -2.0 * log_comp;
}

struct DeltaEntropyRecovery {
  int below = 0, above = 0;
  double maximum_correction = 0;
};

double delta_entropy_tolerance(const int k) {
  return std::min(1e-12,8.0*k*std::numeric_limits<double>::epsilon());
}

// Called only for a tiny out-of-range entropy. Do not change healthy rows or
// their arithmetic, and never normalize malformed probabilities into validity.
bool delta_entropy_repairable_row(const arma::mat& probabilities, const int row,
                                  const double tolerance) {
  double total = 0;
  for (arma::uword j=0; j<probabilities.n_cols; ++j) {
    const double p = probabilities(row,j);
    if (!std::isfinite(p) || p<0 || p>1) return false;
    total+=p;
  }
  return std::abs(total-1)<=tolerance;
}

arma::vec delta_entropy_internal(const arma::mat& probabilities,
                                 const int entropy_type,
                                 DeltaEntropyRecovery* recovered = NULL) {
  const int n = static_cast<int>(probabilities.n_rows);
  const int k = static_cast<int>(probabilities.n_cols);
  if (k < 2) {
    Rcpp::stop("Delta requires at least two states.");
  }

  arma::vec out(n, arma::fill::zeros);
  const double kd = static_cast<double>(k);
  const double uniform = 1.0 / kd;
  // Same log(k) value, shared across nodes; do not replace division with a
  // reciprocal multiply, which can alter endpoint rounding and Delta's RNG.
  const double log_k = entropy_type==2 ? std::log(kd) : 0;
  const double tolerance = delta_entropy_tolerance(k);
  DeltaEntropyRecovery recovery;

  for (int i = 0; i < n; ++i) {
    // Validate every row, including rows whose entropy happens to be interior.
    // This precedes ALL endpoint draws and does not renormalize healthy rows.
    if (!delta_entropy_repairable_row(probabilities,i,tolerance))
      Rcpp::stop("Delta probabilities are invalid and cannot be repaired as roundoff.");
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
      out(i) = total / log_k;
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
    const double value = out(i);
    // Most nodes are interior. Fuse validation and endpoint detection instead
    // of making every healthy node perform four comparisons.
    if (value>0 && value<1) continue;
    if (value==0) { has_zero=true; continue; }
    if (value==1) { has_one=true; continue; }
    // All three normalized entropy definitions have range [0,1]. Only repair
    // a machine-scale excursion supported by a valid probability row.
    const double endpoint = value<0 ? 0 : 1;
    const double correction = std::abs(value-endpoint);
    if (!std::isfinite(value) || correction>tolerance ||
        !delta_entropy_repairable_row(probabilities,i,tolerance))
      Rcpp::stop("Delta entropy is outside [0, 1] and cannot be repaired as roundoff.");
    if (endpoint==0) { ++recovery.below; has_zero=true; }
    else { ++recovery.above; has_one=true; }
    recovery.maximum_correction=std::max(recovery.maximum_correction,correction);
    out(i)=endpoint;
  }
  if (recovered!=NULL) *recovered=recovery;
  if (recovery.below+recovery.above>0)
    Rcpp::warning("Delta repaired entropy roundoff in %d row(s) before endpoint jitter.",
      recovery.below+recovery.above);
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

// Fixed, valid (0,1) acceptance bounds: reuse R's generator and endpoint
// rejection without generic runif bounds/scaling. Keep the conditional draw.
inline double delta_unif_open() {
  double u;
  do { u=::unif_rand(); } while (u<=0 || u>=1);
  return u;
}

// Lognormal random-walk MH for the Beta likelihood and independent exponential
// priors on the raw alpha/beta parameters. Cache accepted-state lgamma values.
double delta_mh_update(const double value, const double other,
                      double& gamma_value, double& gamma_sum,
                      double& log_value,
                      const double n, const double penalty,
                      const double proposal_sd, const bool alpha_update) {
  const double current = n * (gamma_sum - gamma_value) - value * penalty;
  const double log_candidate = R::rnorm(log_value, proposal_sd);
  const double candidate = std::exp(log_candidate);
  // Reject invalid proposals once: retrying would condition the proposal law.
  if (!std::isfinite(candidate) || candidate <= 0.0) return value;
  const double proposed_sum = std::lgamma(alpha_update ? candidate + other : other + candidate);
  const double proposed_value = std::lgamma(candidate);
  const double proposed = n * (proposed_sum - proposed_value) - candidate * penalty;
  const double log_ratio = proposed - current + log_candidate - log_value;
  if (std::isfinite(log_ratio)) {
    if (log_ratio >= 0.0 || std::log(delta_unif_open()) < log_ratio) {
      gamma_value = proposed_value;
      gamma_sum = proposed_sum;
      // Cache log(exp(log_candidate)), not the proposal's log_candidate:
      // round-tripping exactly preserves the next proposal argument.
      log_value = std::log(candidate);
      return candidate;
    }
    return value;
  }
  return value;
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
                           const double sum_log_1_minus_x,
                           arma::vec* alpha_trace = nullptr,
                           arma::vec* beta_trace = nullptr) {
  double alpha_sum = 0.0;
  double beta_sum = 0.0;
  double ratio_sum = 0.0;
  double saved = 0.0;
  double gamma_alpha = std::lgamma(alpha), gamma_beta = std::lgamma(beta);
  double gamma_sum = std::lgamma(alpha + beta);
  double log_alpha = std::log(alpha), log_beta = std::log(beta);
  const double n = static_cast<double>(x.n_elem);
  const double alpha_penalty = lambda0 - sum_log_x;
  const double beta_penalty = lambda0 - sum_log_1_minus_x;
  int until_save = burn;

  for (int iter = 0; iter < sim; ++iter) {
    alpha = delta_mh_update(
      alpha, beta, gamma_alpha, gamma_sum, log_alpha, n, alpha_penalty, proposal_sd, true
    );
    beta = delta_mh_update(
      beta, alpha, gamma_beta, gamma_sum, log_beta, n, beta_penalty, proposal_sd, false
    );

    // First save at burn, then every thin iterations. Avoid a division in
    // every sampler step without changing the saved floating-point sums.
    if (--until_save == 0) {
      until_save = thin;
      if (alpha_trace != nullptr && beta_trace != nullptr) {
        const arma::uword trace_index = static_cast<arma::uword>(saved);
        if (trace_index < alpha_trace->n_elem &&
            trace_index < beta_trace->n_elem) {
          (*alpha_trace)(trace_index) = alpha;
          (*beta_trace)(trace_index) = beta;
        }
      }
      alpha_sum += alpha;
      beta_sum += beta;
      ratio_sum += beta / alpha;
      saved += 1.0;
    }
  }

  return arma::vec({alpha_sum, beta_sum, saved, ratio_sum});
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

arma::vec phylo_d_sum_binary_batch(const arma::mat& states,
                                   const Rcpp::IntegerMatrix& edge,
                                   const arma::vec& edge_length,
                                   const int n_tip,
                                   const int n_threads) {
  const fastphylosig::DContrasts contrasts =
    fastphylosig::DContrasts::from_edges(edge, edge_length.memptr(), n_tip);
  return arma::vec(contrasts.batch(states.memptr(), states.n_cols, n_threads));
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

  return k_stats_chol_matrix(
    X, cholC, inv_one, norm_const, sum_invC
  );
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

  const arma::vec observed = k_stats_chol_matrix(
    X, cholC, inv_one, norm_const, sum_invC
  );
  arma::vec ge_count(p, arma::fill::zeros);
  arma::mat sim_k(nsim, p);

#ifdef _OPENMP
  if (n_threads_eff > 1) {
    std::vector<arma::vec> local_counts(n_threads_eff);
    for (int t = 0; t < n_threads_eff; ++t) local_counts[t].zeros(p);
#pragma omp parallel num_threads(n_threads_eff)
    {
      const int tid = omp_get_thread_num();
      arma::mat permuted(n, p);
      arma::vec kval;
#pragma omp for schedule(static)
      for (int i = 0; i < nsim; ++i) {
        for (int r = 0; r < n; ++r) {
          permuted.row(r) = X.row(perms(i, r) - 1);
        }
        kval = k_stats_chol_matrix(
          permuted, cholC, inv_one, norm_const, sum_invC
        );
        sim_k.row(i) = kval.t();
        for (int j = 0; j < p; ++j) {
          if (fastphylosig::inclusive_upper_tail(kval(j), observed(j))) {
            local_counts[tid](j) += 1.0;
          }
        }
      }
    }
    for (int t = 0; t < n_threads_eff; ++t) ge_count += local_counts[t];
  } else
#endif
  {
    arma::mat permuted(n, p);
    for (int i = 0; i < nsim; ++i) {
      for (int r = 0; r < n; ++r) {
        permuted.row(r) = X.row(perms(i, r) - 1);
      }
      const arma::vec kval = k_stats_chol_matrix(
        permuted, cholC, inv_one, norm_const, sum_invC
      );
      sim_k.row(i) = kval.t();
      for (int j = 0; j < p; ++j) {
        if (fastphylosig::inclusive_upper_tail(kval(j), observed(j))) {
          ge_count(j) += 1.0;
        }
      }
    }
  }
  arma::vec p_value = ge_count / static_cast<double>(nsim);

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

  const arma::vec observed = k_stats_chol_matrix(
    X, cholC, inv_one, norm_const, sum_invC
  );
  arma::vec ge_count(p, arma::fill::zeros);

#ifdef _OPENMP
  if (n_threads_eff > 1) {
    std::vector<arma::vec> local_counts(n_threads_eff);
    for (int t = 0; t < n_threads_eff; ++t) local_counts[t].zeros(p);
#pragma omp parallel num_threads(n_threads_eff)
    {
      const int tid = omp_get_thread_num();
      arma::mat permuted(n, p);
      arma::vec kval;
#pragma omp for schedule(static)
      for (int i = 0; i < nsim; ++i) {
        for (int r = 0; r < n; ++r) {
          permuted.row(r) = X.row(perms(i, r) - 1);
        }
        kval = k_stats_chol_matrix(
          permuted, cholC, inv_one, norm_const, sum_invC
        );
        for (int j = 0; j < p; ++j) {
          if (fastphylosig::inclusive_upper_tail(kval(j), observed(j))) {
            local_counts[tid](j) += 1.0;
          }
        }
      }
    }
    for (int t = 0; t < n_threads_eff; ++t) ge_count += local_counts[t];
  } else
#endif
  {
    arma::mat permuted(n, p);
    for (int i = 0; i < nsim; ++i) {
      for (int r = 0; r < n; ++r) {
        permuted.row(r) = X.row(perms(i, r) - 1);
      }
      const arma::vec kval = k_stats_chol_matrix(
        permuted, cholC, inv_one, norm_const, sum_invC
      );
      for (int j = 0; j < p; ++j) {
        if (fastphylosig::inclusive_upper_tail(kval(j), observed(j))) {
          ge_count(j) += 1.0;
        }
      }
    }
  }
  arma::vec p_value = ge_count / static_cast<double>(nsim);

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
double lambda_loglik_spectral_cpp(
    const double lambda,
    const arma::vec& eigval,
    const arma::mat& eigvec,
    const arma::vec& inv_sqrt_diag,
    const double log_diag,
    const arma::vec& y) {
  return lambda_loglik_spectral_internal(
    lambda, eigval, eigvec, inv_sqrt_diag, log_diag, y
  );
}

// [[Rcpp::export]]
double lambda_loglik_spectral_projected_cpp(
    const double lambda,
    const arma::vec& eigval,
    const arma::vec& wproj,
    const arma::vec& zproj,
    const double log_diag) {
  return lambda_loglik_spectral_projected_internal(
    lambda, eigval, wproj, zproj, log_diag
  );
}

// [[Rcpp::export]]
arma::mat lambda_loglik_spectral_projected_batch_cpp(
    const arma::vec& lambda,
    const arma::vec& eigval,
    const arma::vec& wproj,
    const arma::mat& zproj,
    const double log_diag,
    const int n_threads = 1) {
  return lambda_loglik_spectral_projected_batch_internal(
    lambda, eigval, wproj, zproj, log_diag, n_threads
  );
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

  if (!phylo_d_has_polytomy(edge, n_tip)) {
    // A single traversal shares topology and adjusted branch lengths across
    // all simulation columns.  The arithmetic order remains column-wise.
    return phylo_d_sum_binary_batch(
      states, edge, edge_length, n_tip, n_threads
    );
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
    const double threshold =
      fastphylosig::quantile_type7_sorted(values, prop_state1);

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
      fastphylosig::quantile_type7_sorted(sorted_values, prop_state1);

    for (int tip = 0; tip < n_tip; ++tip) {
      out(tip, sim) = tip_values[tip] < threshold ? 1 : 0;
    }
  }

  return out;
}

namespace {
// Uses the ACE rate-matrix and transition routines in this translation unit.
// The immutable transitions and mutable node states live only for one scan.
struct HeterogeneityMkWorkspace {
  int nt, nn, k, root;
  std::vector<int> pre, parent, state, transition;
  arma::cube cumulative;
  HeterogeneityMkWorkspace(const Rcpp::List& core,
      const Rcpp::IntegerMatrix& index, const arma::vec& rates) {
    nt = Rcpp::as<int>(core["n_tip"]); nn = Rcpp::as<int>(core["n_total"]);
    k = index.nrow(); root = Rcpp::as<int>(core["root"]) - 1;
    pre = Rcpp::as<std::vector<int>>(core["preorder"]);
    parent = Rcpp::as<std::vector<int>>(core["parent"]);
    const Rcpp::NumericVector branch = core["branch_length_by_node"];
    bool independent_er = rates.n_elem==1 && rates[0]==R_PosInf;
    if (independent_er) for (int i=0; i<k; ++i) for (int j=0; j<index.ncol(); ++j)
      if (index(i,j)!=(i==j ? 0 : 1)) independent_er=false;
    if (nt < 2 || nn <= nt || k < 2 || index.ncol() != k ||
        pre.size() != static_cast<size_t>(nn) || parent.size() != pre.size() ||
        branch.size() != nn || root < 0 || root >= nn || pre[0] != root + 1 ||
        parent[root] != 0 || (!rates.is_finite() && !independent_er) || arma::any(rates < 0))
      Rcpp::stop("invalid common Mk simulation workspace");
    // Equal branch lengths have exactly the same transition matrix for a
    // fixed Q. Pool them without rounding lengths or changing probabilities.
    // Sorting integer node IDs avoids one hash-node allocation per branch.
    transition.resize(nn, -1);
    std::vector<int> unique;
    {
      std::vector<bool> seen(nn, false);
      for (int pos = 0; pos < nn; ++pos) {
        if ((pos & 1023) == 0) Rcpp::checkUserInterrupt();
        const int node = pre[pos] - 1;
        if (node < 0 || node >= nn || seen[node]) Rcpp::stop("invalid Mk preorder");
        seen[node] = true;
        if (!pos) continue;
        const int p = parent[node] - 1;
        if (p < 0 || p >= nn || p == node || !seen[p] ||
            !std::isfinite(branch[node]) || branch[node] < 0)
          Rcpp::stop("invalid Mk branch/parent");
        unique.push_back(node);
      }
    }
    std::sort(unique.begin(), unique.end(), [&](int a, int b) { return branch[a] < branch[b]; });
    size_t count = 0;
    double previous = -1;
    for (int node : unique) {
      if (branch[node] != previous) {
        previous = branch[node]; unique[count++] = node;
      }
      transition[node] = count - 1;
    }
    unique.resize(count);
    if (static_cast<double>(unique.size()) * k * k * sizeof(double) > 64.0*1024*1024)
      Rcpp::stop("common Mk transition cache exceeds 64 MiB for distinct branch lengths; reduce tree size or state count");
    arma::mat Q;
    if (!independent_er) Q = ace_build_q(index, rates);
    arma::cx_vec lambda; arma::cx_mat gamma, inverse;
    bool eigen = false;
    if (!independent_er) try { eigen = ace_eigen_decomp(Q, lambda, gamma, inverse); } catch (...) { eigen = false; }
    cumulative.set_size(k, k, unique.size());
    for (size_t slice = 0; slice < unique.size(); ++slice) {
      if ((slice & 1023) == 0) Rcpp::checkUserInterrupt();
      const int node = unique[slice];
      arma::mat P;
      if (independent_er) {
        if (branch[node]==0) P.eye(k,k);
        else { P.ones(k,k); P/=k; }
      } else P = eigen ? ace_transition_matrix(lambda, gamma, inverse, branch[node]) :
        arma::expmat(Q * branch[node]);
      // A defective/ill-conditioned eigensystem can produce a bad transition.
      if (!ace_valid_transition(P))
        P = arma::expmat(Q * branch[node]);
      if (!ace_valid_transition(P))
        Rcpp::stop("common Mk transition probabilities are numerically invalid");
      for (int i = 0; i < k; ++i) {
        double total = 0;
        for (int j = 0; j < k; ++j) total += std::max(0.0, P(i,j));
        double cdf = 0;
        for (int j = 0; j < k; ++j) {
          cdf += std::max(0.0, P(i,j))/total;
          cumulative(i,j,slice) = j == k - 1 ? 1.0 : cdf;
        }
      }
    }
    state.resize(nn);
  }
};
}

// [[Rcpp::export(rng = false)]]
SEXP heterogeneity_mk_prepare_cpp(const Rcpp::List& core,
    const Rcpp::IntegerMatrix& rate_index, const arma::vec& rates) {
  Rcpp::XPtr<HeterogeneityMkWorkspace> pointer(new HeterogeneityMkWorkspace(core, rate_index, rates), true);
  pointer.attr("class") = "fastphylosig_mk_workspace";
  return pointer;
}

// [[Rcpp::export(rng = false)]]
void heterogeneity_mk_release_cpp(SEXP workspace) {
  if (TYPEOF(workspace) != EXTPTRSXP || !Rf_inherits(workspace,"fastphylosig_mk_workspace"))
    Rcpp::stop("invalid common Mk workspace");
  Rcpp::XPtr<HeterogeneityMkWorkspace> ptr(workspace);
  ptr.release();
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix heterogeneity_mk_simulate_cpp(SEXP workspace, const int nsim) {
  if (TYPEOF(workspace) != EXTPTRSXP || !Rf_inherits(workspace, "fastphylosig_mk_workspace") ||
      R_ExternalPtrAddr(workspace) == nullptr || nsim < 1)
    Rcpp::stop("invalid common Mk simulation request");
  Rcpp::XPtr<HeterogeneityMkWorkspace> ptr(workspace);
  HeterogeneityMkWorkspace& w = *ptr;
  if (static_cast<double>(w.nt)*nsim > std::numeric_limits<int>::max())
    Rcpp::stop("common Mk simulation output is too large");
  Rcpp::IntegerMatrix out(w.nt, nsim);
  for (int b = 0; b < nsim; ++b) {
    Rcpp::checkUserInterrupt();
    w.state[w.root] = std::min(w.k-1, static_cast<int>(R::runif(0, w.k)));
    for (int pos = 1; pos < w.nn; ++pos) {
      if ((pos & 8191) == 0) Rcpp::checkUserInterrupt();
      const int node = w.pre[pos] - 1, p = w.parent[node] - 1;
      const double u = R::runif(0, 1);
      int state = 0;
      while (state < w.k - 1 && u >= w.cumulative(w.state[p], state, w.transition[node])) ++state;
      w.state[node] = state;
    }
    for (int tip = 0; tip < w.nt; ++tip) out(tip,b) = w.state[tip]+1;
  }
  return out;
}

// The proposal density must describe the actual normalized CDF sampler, not
// a separate eigensystem evaluation. This read-only operation consumes no RNG.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix heterogeneity_mk_loglik_bounds_cpp(SEXP workspace,
    const Rcpp::List& core, const Rcpp::IntegerMatrix& states) {
  if(TYPEOF(workspace)!=EXTPTRSXP || !Rf_inherits(workspace,"fastphylosig_mk_workspace"))
    Rcpp::stop("invalid common Mk workspace");
  Rcpp::XPtr<HeterogeneityMkWorkspace> w(workspace);
  if(!w.get())Rcpp::stop("released common Mk workspace");
  hetero_bounds::Tree tree(core,w->k);
  if(tree.nt!=w->nt || tree.nn!=w->nn || tree.root!=w->root ||
     tree.pre!=w->pre || tree.parent!=w->parent)
    Rcpp::stop("proposal workspace and tree do not match");
  std::vector<std::vector<hetero_bounds::I>> p(tree.lengths.size());
  for(int node:tree.pre) {
    --node;if(node==tree.root)continue;
    auto& slice=p[tree.length_id[node]];if(!slice.empty())continue;
    slice.resize(tree.k*tree.k);const int source=w->transition[node];
    for(int i=0;i<tree.k;++i) {
      double previous=0;
      for(int j=0;j<tree.k;++j) {
        const double a=hetero_bounds::sampler_cdf(w->cumulative(i,j,source),previous);
        slice[i*tree.k+j]=hetero_bounds::prob(hetero_bounds::sub(hetero_bounds::I(a),hetero_bounds::I(previous)));
        previous=a;
      }
    }
  }
  return hetero_bounds::matrix(tree.likelihood(states,p));
}

// Enter the RNG scope once per fit, just as the original first objective call
// did (including initialization when .Random.seed is absent). No draws occur.
// [[Rcpp::export]]
SEXP fast_ace_fit_prepare_cpp(const Rcpp::IntegerMatrix& edge,
                              const arma::vec& edge_length,
                              const Rcpp::IntegerVector& tip_state,
                              const Rcpp::IntegerMatrix& rate_index,
                              const int n_threads = 0) {
  if (n_threads < 0) Rcpp::stop("n_threads must be non-negative (0 selects automatic threads).");
  return Rcpp::XPtr<AceFitWorkspace>(new AceFitWorkspace(
    edge,edge_length,tip_state,rate_index,n_threads),true);
}

// [[Rcpp::export(rng = false)]]
double fast_ace_fit_deviance_cpp(SEXP prepared, const arma::vec& par,
                               const bool recovery_guard = false) {
  Rcpp::XPtr<AceFitWorkspace> ws(prepared);
  if (ws.get() == NULL) Rcpp::stop("ACE fit workspace has been released.");
  if (arma::any(par < 0.0) || !par.is_finite()) return 1e50;
  const double dev = ace_discrete_pruning_internal(
    ws->edge, ws->edge_length, ws->tip_state, ws->rate_index, par,
    NULL, false, true, &ws->scratch, false, recovery_guard);
  return std::isfinite(dev) ? dev : 1e50;
}

// [[Rcpp::export(rng = false)]]
double fast_ace_er_deviance_cpp(SEXP prepared, double rate) {
  Rcpp::XPtr<AceFitWorkspace> w(prepared);
  if (w.get()==NULL) Rcpp::stop("ACE fit workspace has been released.");
  ace_er_reference_check(*w);
  return ace_er_reference_deviance(*w,rate);
}

// Explicit experimental common-Q fit, never used for observed Delta or ACE.
// [[Rcpp::export(rng = false)]]
Rcpp::List fast_ace_er_reference_cpp(SEXP prepared) {
  Rcpp::XPtr<AceFitWorkspace> w(prepared);
  if (w.get()==NULL) Rcpp::stop("ACE fit workspace has been released.");
  ace_er_reference_check(*w);
  const double eps=std::numeric_limits<double>::epsilon(), k=w->rate_index.nrow();
  const double lower=std::max(std::log(std::numeric_limits<double>::min()),
    std::log(eps)-std::log(k)-std::log(w->edge_length.max()));
  const double upper=std::min(std::log(std::numeric_limits<double>::max())-1,
    std::log(-std::log(eps))-std::log(k)-std::log(w->edge_length.min()));
  const int intervals=std::max(2,static_cast<int>(std::ceil((upper-lower)/std::log(2.0))));
  int evaluations=0;
  const auto objective=[&](double q) {
    if ((++evaluations & 31)==0) Rcpp::checkUserInterrupt();
    return ace_er_reference_deviance(*w,q);
  };
  const double zero=objective(0), independent=objective(R_PosInf);
  double best=std::min(zero,independent), rate=zero<=independent ? 0 : R_PosInf;
  const double tie=1e-10*(1+std::abs(independent));
  const auto select=[&](double q,double value) {
    if (value<best-tie) { best=value; rate=q; }
  };
  std::vector<double> grid(intervals+1), values(intervals+1);
  for (int i=0; i<=intervals; ++i) {
    grid[i]=lower+(upper-lower)*i/intervals;
    values[i]=objective(std::exp(grid[i]));
    select(std::exp(grid[i]),values[i]);
  }
  for (int i=1; i<intervals; ++i)
    if (values[i]<=values[i-1] && values[i]<=values[i+1] &&
        (values[i]<values[i-1]-tie || values[i]<values[i+1]-tie)) {
      const double point=brent_minimize(grid[i-1],grid[i+1],1e-8,
        [&](double x) { return objective(std::exp(x)); });
      const double q=std::exp(point); select(q,objective(q));
    }
  const std::string endpoint=rate==0 ? "zero" : (std::isinf(rate) ? "independent" : "interior");
  return Rcpp::List::create(Rcpp::Named("rates")=Rcpp::NumericVector::create(rate),
    Rcpp::Named("loglik")=-.5*best,
    Rcpp::Named("convergence")=std::isfinite(best) ? 0 : 1,
    Rcpp::Named("diagnostics")=Rcpp::List::create(
      Rcpp::Named("method")="ER log-grid + Brent (analytic transitions)",
      Rcpp::Named("endpoint")=endpoint,
      Rcpp::Named("independence_deviance_gap")=independent-best,
      Rcpp::Named("uniform_root_loglik")=-.5*best-std::log(k),
      Rcpp::Named("objective_evaluations")=evaluations,
      Rcpp::Named("grid_points")=intervals+1));
}

// Same RNG-scope boundary as preparation: states change between fits, while
// topology/rate mapping and exact exponent slots remain immutable.
// [[Rcpp::export]]
void fast_ace_fit_update_cpp(SEXP prepared, const Rcpp::IntegerVector& states) {
  Rcpp::XPtr<AceFitWorkspace> ws(prepared);
  if (ws.get() == NULL) Rcpp::stop("ACE fit workspace has been released.");
  if (states.size() != ws->tip_state.size()) Rcpp::stop("ACE state count does not match the prepared tree.");
  const int k = ws->rate_index.nrow();
  for (int state : states) if (!Rcpp::IntegerVector::is_na(state) && (state < 1 || state > k))
    Rcpp::stop("tip_state values must be 1-based state indices.");
  std::copy(states.begin(),states.end(),ws->tip_state.begin());
  ws->scratch.initialized = false;
  ws->scratch.er_recoveries = 0;
  ws->scratch.stable_recoveries = ws->scratch.probability_recoveries = 0;
}

// One read per fit, no RNG scope and no per-objective R callbacks.
// [[Rcpp::export(rng = false)]]
double fast_ace_fit_recoveries_cpp(SEXP prepared) {
  Rcpp::XPtr<AceFitWorkspace> ws(prepared);
  if (ws.get() == NULL) Rcpp::stop("ACE fit workspace has been released.");
  return ws->scratch.er_recoveries;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector fast_ace_fit_stability_cpp(SEXP prepared) {
  Rcpp::XPtr<AceFitWorkspace> ws(prepared);
  if (ws.get() == NULL) Rcpp::stop("ACE fit workspace has been released.");
  return Rcpp::NumericVector::create(ws->scratch.stable_recoveries,
    ws->scratch.probability_recoveries);
}

// [[Rcpp::export]]
Rcpp::List fast_ace_fit_liks_cpp(SEXP prepared, const arma::vec& par,
                                const bool marginal = false,
                                const bool recovery_guard = false) {
  Rcpp::XPtr<AceFitWorkspace> ws(prepared);
  if (ws.get() == NULL) Rcpp::stop("ACE fit workspace has been released.");
  if (arma::any(par < 0.0) || !par.is_finite()) Rcpp::stop("par must be finite and non-negative.");
  arma::mat likelihood;
  const double dev = ace_discrete_pruning_internal(ws->edge,ws->edge_length,
    ws->tip_state,ws->rate_index,par,&likelihood,marginal,true,&ws->scratch,false,recovery_guard);
  if (!std::isfinite(dev)) Rcpp::stop("ancestral likelihood calculation failed.");
  return Rcpp::List::create(Rcpp::Named("deviance")=dev,
    Rcpp::Named("loglik")=-.5*dev,Rcpp::Named("lik.anc")=likelihood);
}

// Storage budget counts owned buffers, not R data or allocator overhead.
// [[Rcpp::export(rng = false)]]
double fast_ace_fit_bytes_cpp(SEXP prepared) {
  Rcpp::XPtr<AceFitWorkspace> ptr(prepared);
  if (ptr.get() == NULL) return 0;
  const AceFitWorkspace& w = *ptr;
  const AcePruningWorkspace& s = w.scratch;
  double bytes = sizeof(w) + sizeof(int)*(w.edge.size()+w.rate_index.size()+w.tip_state.size())
    + sizeof(double)*(w.edge_length.n_elem+s.liks.n_elem+s.left.n_elem+s.right.n_elem+s.product.n_elem);
  bytes += sizeof(std::complex<double>)*(s.lambda.n_elem+s.gamma.n_elem+s.inv_gamma.n_elem+
    s.tip_transformed.n_elem+s.exponents.n_elem+s.exponent_batch.n_elem+
    s.transition.child.n_elem+s.transition.transformed.n_elem+
    s.transition.exp_lambda.n_elem+s.transition.gamma_transformed.n_elem);
  for (const auto& worker : s.exponent_workers) bytes += sizeof(std::complex<double>)*worker.n_elem;
  bytes += sizeof(int)*(s.edge_exp_slot.capacity()+s.exponent_edges.capacity())+
    (s.transformed_ready.capacity()+s.exp_ready.capacity()+7)/8;
  return bytes+s.eigen_exponents.bytes();
}

// [[Rcpp::export(rng = false)]]
void fast_ace_fit_release_cpp(SEXP prepared) {
  Rcpp::XPtr<AceFitWorkspace> ws(prepared);
  ws.release();
}

// [[Rcpp::export]]
Rcpp::NumericVector fast_ace_discrete_deviance_cpp(const Rcpp::IntegerMatrix& edge,
                                      const arma::vec& edge_length,
                                      const Rcpp::IntegerVector& tip_state,
                                       const Rcpp::IntegerMatrix& rate_index,
                                       const arma::vec& par,
                                       const bool reuse_buffers = false) {
  if (arma::any(par < 0.0) || !par.is_finite()) {
    return Rcpp::NumericVector::create(1e50);
  }
  AcePruningWorkspace scratch;
  const double dev = ace_discrete_pruning_internal(
    edge, edge_length, tip_state, rate_index, par, NULL, false, reuse_buffers, &scratch
  );
  Rcpp::NumericVector value = Rcpp::NumericVector::create(std::isfinite(dev) ? dev : 1e50);
  if (scratch.er_recoveries) value.attr("er_recoveries")=static_cast<double>(scratch.er_recoveries);
  if (scratch.stable_recoveries) value.attr("stable_recoveries")=static_cast<double>(scratch.stable_recoveries);
  return value;
}

// [[Rcpp::export]]
Rcpp::List fast_ace_discrete_liks_cpp(const Rcpp::IntegerMatrix& edge,
                                      const arma::vec& edge_length,
                                      const Rcpp::IntegerVector& tip_state,
                                      const Rcpp::IntegerMatrix& rate_index,
                                       const arma::vec& par,
                                       const bool marginal = false,
                                       const bool reuse_buffers = false) {
  if (arma::any(par < 0.0) || !par.is_finite()) {
    Rcpp::stop("par must be finite and non-negative.");
  }
  arma::mat lik_anc;
  AcePruningWorkspace scratch;
  const double dev = ace_discrete_pruning_internal(
    edge, edge_length, tip_state, rate_index, par, &lik_anc, marginal,
    reuse_buffers, &scratch
  );
  if (!std::isfinite(dev)) {
    Rcpp::stop("ancestral likelihood calculation failed.");
  }
  Rcpp::List value = Rcpp::List::create(
    Rcpp::Named("deviance") = dev,
    Rcpp::Named("loglik") = -0.5 * dev,
    Rcpp::Named("lik.anc") = lik_anc
  );
  if (scratch.er_recoveries) value["er_recoveries"]=static_cast<double>(scratch.er_recoveries);
  if (scratch.stable_recoveries) value["stable_recoveries"]=static_cast<double>(scratch.stable_recoveries);
  if (scratch.probability_recoveries) value["probability_recoveries"]=static_cast<double>(scratch.probability_recoveries);
  return value;
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
                          const int entropy_type,
                          const bool return_chains = true) {
  if (sim < 1 || thin < 1 || burn < 1 || burn > sim) {
    Rcpp::stop("sim, thin, and burn must satisfy 1 <= burn <= sim.");
  }
  if (!std::isfinite(lambda0) || !std::isfinite(proposal_sd) ||
      lambda0 <= 0.0 || proposal_sd <= 0.0) {
    Rcpp::stop("lambda0 and proposal_sd must be finite and positive.");
  }

  DeltaEntropyRecovery recovery;
  arma::vec x = delta_entropy_internal(probabilities, entropy_type, &recovery);
  if (arma::any(x <= 0.0) || arma::any(x >= 1.0) || !x.is_finite()) {
    Rcpp::stop("Delta entropies must be finite and strictly inside (0, 1).");
  }

  const double sum_log_x = arma::accu(arma::log(x));
  const double sum_log_1_minus_x = arma::accu(arma::log(1.0 - x));

  // The trace length is determined solely by the controls.  Keeping this
  // allocation and the trace writes local to the existing sampler preserves
  // its random-number call order and its Delta estimate exactly.
  const int saved_per_chain = 1 + (sim - burn) / thin;
  arma::vec alpha_trace1;
  arma::vec beta_trace1;
  arma::vec alpha_trace2;
  arma::vec beta_trace2;
  if (return_chains) {
    alpha_trace1.set_size(saved_per_chain);
    beta_trace1.set_size(saved_per_chain);
    alpha_trace2.set_size(saved_per_chain);
    beta_trace2.set_size(saved_per_chain);
  }

  // Explicit sequencing makes the initialization RNG contract portable.
  const double alpha1 = R::rexp(1.0), beta1 = R::rexp(1.0);
  const arma::vec chain1 = delta_chain_sums(
    alpha1, beta1, x, lambda0, proposal_sd, sim, thin, burn,
    sum_log_x, sum_log_1_minus_x,
    return_chains ? &alpha_trace1 : nullptr,
    return_chains ? &beta_trace1 : nullptr
  );
  const double alpha2 = R::rexp(1.0), beta2 = R::rexp(1.0);
  const arma::vec chain2 = delta_chain_sums(
    alpha2, beta2, x, lambda0, proposal_sd, sim, thin, burn,
    sum_log_x, sum_log_1_minus_x,
    return_chains ? &alpha_trace2 : nullptr,
    return_chains ? &beta_trace2 : nullptr
  );

  const double saved = chain1(2) + chain2(2);
  if (saved <= 0.0) {
    Rcpp::stop("No MCMC samples were saved.");
  }
  const double alpha_mean = (chain1(0) + chain2(0)) / saved;
  const double beta_mean = (chain1(1) + chain2(1)) / saved;
  const double delta = (chain1(3) + chain2(3)) / saved;

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("delta") = delta,
    Rcpp::Named("alpha_mean") = alpha_mean,
    Rcpp::Named("beta_mean") = beta_mean,
    Rcpp::Named("n_saved") = saved,
    Rcpp::Named("entropy") = x
  );
  if (recovery.below+recovery.above>0)
    out["entropy_recovery"]=Rcpp::List::create(
      Rcpp::Named("below")=recovery.below,Rcpp::Named("above")=recovery.above,
      Rcpp::Named("maximum_correction")=recovery.maximum_correction,
      Rcpp::Named("tolerance")=delta_entropy_tolerance(probabilities.n_cols));
  if (return_chains) {
    out["alpha_chain"] = arma::join_rows(alpha_trace1, alpha_trace2);
    out["beta_chain"] = arma::join_rows(beta_trace1, beta_trace2);
  } else {
    // Keep a stable result schema for permutation workers while avoiding any
    // trace allocation or copy when diagnostics are not requested.
    out["alpha_chain"] = arma::mat(0, 0);
    out["beta_chain"] = arma::mat(0, 0);
  }
  return out;
}
