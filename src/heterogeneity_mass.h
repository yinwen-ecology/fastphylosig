#ifndef FASTPHYLOSIG_HETEROGENEITY_MASS_H
#define FASTPHYLOSIG_HETEROGENEITY_MASS_H
#include "heterogeneity_taylor.h"
#include "heterogeneity_ard_boundary.h"

namespace hetero_mass {
using hetero_bounds::I;
// Whole-domain constraints on mutually exclusive phenotype probabilities.
// Three-state CTMCs share a stationary vector and satisfy the return bound
// P_jj(t)>=max(0,(5*pi_j-1)/4). With a uniform root this gives sum L_mono>=6^-E.
struct Bank {
  bool eligible=true;
  int observed_constant=-1;
  double floor=0;
  double group_upper=1;
  I observed_density;
  std::vector<I> a,b;
  std::array<int,3> constants{{-1,-1,-1}};
  Bank(const hetero_bounds::Tree& tree,const Rcpp::IntegerMatrix& x,
      const std::vector<int>& ids,const std::vector<I>& proposal,const std::vector<int>& tail):
      a(x.ncol()),b(x.ncol()) {
    int observed=0;
    for(int i=0;i<x.nrow();++i) {
      observed+=x(i,0)!=NA_INTEGER;
      for(int p=1;p<x.ncol();++p)if((x(i,p)==NA_INTEGER)!=(x(i,0)==NA_INTEGER))eligible=false;
    }
    if(!observed)eligible=false;
    if(!eligible)return;
    I f(1),sixth=hetero_bounds::divpos(I(1),I(6));int edges=tree.nn-1;
    while(edges) {if(edges&1)f=hetero_bounds::posmul(f,sixth);edges>>=1;if(edges)sixth=hetero_bounds::posmul(sixth,sixth);}
    floor=f.lo;
    for(int p=0;p<x.ncol();++p) {
      int state=0;
      for(int i=0;i<x.nrow();++i)if(x(i,p)!=NA_INTEGER) {
        if(state&&state!=x(i,p)) {state=0;break;}state=x(i,p);
      }
      if(state) {constants[state-1]=p;if(!p)observed_constant=state-1;}
    }
    observed_density=hetero_taylor::positive_exp(proposal[0]);
    double pivot=-hetero_bounds::inf;for(I g:proposal)pivot=std::max(pivot,-g.lo);
    for(size_t i=0;i<ids.size();++i) {
      I weight=hetero_bounds::exp_nonpositive(hetero_taylor::minus(I(-pivot),proposal[i]));
      I& target=tail[i] ? a[ids[i]] : b[ids[i]];
      target=hetero_bounds::exact_posadd(target,weight);
    }
    group_upper=0;
    for(size_t p=0;p<a.size();++p)if(a[p].hi>0||b[p].hi>0)
      group_upper=std::max(group_upper,hetero_bounds::fraction(a[p],b[p]).hi);
  }
  bool useful(double gamma) const {
    if(!eligible)return false;
    if(gamma>0&&group_upper<1)return true;
    if(gamma>0&&b[0].lo>0)return true;
    if(!(floor>0))return false;
    for(int p:constants)if(p<0||b[p].lo<=0)return false;
    return true;
  }
  I upper(double gamma,I original) const {
    if(!useful(gamma))return original;
    // The rank is a weighted mean of phenotype-specific tail ratios. C_gamma
    // with gamma>0 guarantees a positive denominator, even if h underflows.
    if(gamma>0)original.hi=std::min(original.hi,group_upper);
    double h=hetero_bounds::posmul(I(gamma),observed_density).lo,s=floor;
    if(h>1)return I(1,0);
    if(observed_constant>=0)s=std::max(0.0,hetero_bounds::dn(s-h));
    I allocated=hetero_bounds::posadd(I(h),I(s));
    if(allocated.lo>1)return I(1,0);
    // A rounding-overlap at total mass one must relax, never exclude.
    if(allocated.hi>1) {s=0;allocated=I(h);}
    const I remaining=hetero_bounds::prob(hetero_bounds::sub(I(1),allocated));
    auto certificate=[&](double u) {
      double largest=0,mono=-hetero_bounds::inf;I obs;
      std::array<I,3> mono_coeff{};
      for(size_t p=0;p<a.size();++p) {
        I value=hetero_taylor::minus(hetero_taylor::times(hetero_taylor::minus(I(1),I(u)),a[p]),
          hetero_taylor::times(I(u),b[p]));
        for(int j=0;j<3;++j)if(constants[j]==static_cast<int>(p))mono_coeff[j]=value;
        largest=std::max(largest,value.hi);if(!p)obs=value;
      }
      for(I c:mono_coeff)mono=std::max(mono,c.hi);
      I bound=hetero_taylor::plus(hetero_taylor::times(I(h),obs),
        hetero_taylor::plus(hetero_taylor::times(I(s),I(mono)),hetero_taylor::times(remaining,I(largest))));
      return bound.hi<=0;
    };
    double lo=0,hi=original.hi;
    for(int step=0;step<48;++step) {double mid=lo/2+hi/2;if(certificate(mid))hi=mid;else lo=mid;}
    return I(original.lo,hi);
  }
};
}
#endif
