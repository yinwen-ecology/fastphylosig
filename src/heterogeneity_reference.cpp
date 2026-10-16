#include <Rcpp.h>
#include <array>
#include "heterogeneity_bounds.h"
#include "heterogeneity_bootstrap.h"

namespace {
// Shared by simulation and density: these are the actual sequential sampling
// thresholds, including floating-point normalization and the final endpoint.
void count_cdf(const std::array<double,8>& count,int k,double total,
    std::array<double,8>& cdf) {
  double raw=0,previous=0;
  for(int j=0;j<k;++j) {
    raw+=count[j]/total;
    previous=hetero_bounds::sampler_cdf(j==k-1 ? 1.0 : raw,previous);
    cdf[j]=previous;
  }
}
void validate_states(int k) {
  if(k<2||k>8)Rcpp::stop("count proposal requires 2--8 states");
}
}

// [[Rcpp::export]]
Rcpp::IntegerVector heterogeneity_count_simulate_cpp(const int k,
    const Rcpp::LogicalVector& observed) {
  validate_states(k);
  for(int x:observed)if(x==NA_LOGICAL)Rcpp::stop("invalid count proposal mask");
  Rcpp::IntegerVector out(observed.size(),NA_INTEGER);
  std::array<double,8> count{},cdf{};count.fill(1);double total=k;
  for(int i=0;i<observed.size();++i)if(observed[i]) {
    if((i&1023)==0)Rcpp::checkUserInterrupt();
    count_cdf(count,k,total,cdf);
    double u=R::runif(0,1);int j=0;
    while(j<k-1&&u>=cdf[j])++j;
    out[i]=j+1;++count[j];++total;
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix heterogeneity_count_loglik_cpp(const Rcpp::IntegerMatrix& x,
    const int k) {
  validate_states(k);
  if(!x.ncol())Rcpp::stop("empty count proposal bank");
  std::vector<hetero_bounds::I> answer;answer.reserve(x.ncol());
  for(int b=0;b<x.ncol();++b) {
    if((b&15)==0)Rcpp::checkUserInterrupt();
    std::array<double,8> count{},cdf{};count.fill(1);double total=k;
    hetero_bounds::I log(0);
    for(int i=0;i<x.nrow();++i) {
      int state=x(i,b);
      if((state==NA_INTEGER)!=(x(i,0)==NA_INTEGER))Rcpp::stop("count proposal requires a fixed missing mask");
      if(state==NA_INTEGER)continue;
      if(state<1||state>k)Rcpp::stop("invalid count proposal state");
      count_cdf(count,k,total,cdf);int j=state-1;
      auto probability=hetero_bounds::prob(hetero_bounds::sub(
        hetero_bounds::I(cdf[j]),hetero_bounds::I(j ? cdf[j-1] : 0)));
      log=hetero_bounds::add(log,hetero_bounds::log_interval(probability));
      ++count[j];++total;
    }
    answer.push_back(log);
  }
  return hetero_bounds::matrix(answer);
}

// [[Rcpp::export(rng = false)]]
double heterogeneity_delta_between_cpp(const Rcpp::NumericVector& local) {
  std::vector<double> values;values.reserve(local.size());
  for(double x:local) {
    if(!(x>0)||!std::isfinite(x))return NA_REAL;
    values.push_back(std::log(x));
  }
  return fastphylosig::hetero_sd(values);
}
