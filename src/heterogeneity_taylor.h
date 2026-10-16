#ifndef FASTPHYLOSIG_HETEROGENEITY_TAYLOR_H
#define FASTPHYLOSIG_HETEROGENEITY_TAYLOR_H
#include "heterogeneity_bounds.h"
#include <array>
#include <stdexcept>

// First-order Taylor models with six shared rate coordinates. The remainder
// contains every nonlinear, truncation and rounding error; no local-linear
// approximation is ever substituted for a whole-box certificate.
namespace hetero_taylor {
using hetero_bounds::I;
using hetero_bounds::add;
using hetero_bounds::sub;
using hetero_bounds::mul;
inline bool zero(I a) { return a.lo==0&&a.hi==0; }
inline I plus(I a,I b) { return zero(a) ? b : zero(b) ? a : add(a,b); }
inline I times(I a,I b) { return zero(a)||zero(b) ? I(0) : mul(a,b); }
inline I minus(I a,I b) { return plus(a,I(-b.hi,-b.lo)); }
inline I meet(I a,I b) {
  I z(std::max(a.lo,b.lo),std::min(a.hi,b.hi));
  if(z.lo>z.hi)Rcpp::stop("inconsistent shared-Q enclosure");
  return z;
}
struct T {
  I c,r,box;
  std::array<I,6> a{};
  explicit T(I value=I(0)):c(value),r(0),box(value) {}
  I linear() const {
    double radius=0;
    for(I x:a)if(!zero(x))radius=hetero_bounds::up(radius+std::max(std::abs(x.lo),std::abs(x.hi)));
    return I(-radius,radius);
  }
  void finish() {
    I z=plus(plus(c,r),linear());
    if(!std::isfinite(z.lo)||!std::isfinite(z.hi))
      throw std::overflow_error("shared-Q Taylor range overflow");
    box=meet(box,z);
  }
  bool constant() const {
    if(!zero(r))return false;
    for(I x:a)if(!zero(x))return false;
    return true;
  }
};
inline T variable(I bounds,int coordinate) {
  if(bounds.lo==bounds.hi)return T(bounds);
  T z;const double center=bounds.lo/2+bounds.hi/2;
  z.c=I(center);z.box=bounds;
  z.a[coordinate]=I(std::max(hetero_bounds::up(center-bounds.lo),
    hetero_bounds::up(bounds.hi-center)));
  z.finish();return z;
}
inline T sum(const T& x,const T& y) {
  T z;z.c=plus(x.c,y.c);z.r=plus(x.r,y.r);z.box=plus(x.box,y.box);
  for(int j=0;j<6;++j)z.a[j]=plus(x.a[j],y.a[j]);
  z.finish();return z;
}
inline T scale(const T& x,I value) {
  if(zero(value))return T();
  T z;z.c=times(x.c,value);z.r=times(x.r,value);z.box=times(x.box,value);
  for(int j=0;j<6;++j)z.a[j]=times(x.a[j],value);
  z.finish();return z;
}
inline T difference(const T& x,const T& y) { return sum(x,scale(y,I(-1))); }
inline T product(const T& x,const T& y) {
  if(x.constant())return scale(y,meet(x.c,x.box));
  if(y.constant())return scale(x,meet(y.c,y.box));
  T z;z.c=times(x.c,y.c);z.box=times(x.box,y.box);
  for(int j=0;j<6;++j)z.a[j]=plus(times(x.c,y.a[j]),times(y.c,x.a[j]));
  I u=meet(plus(x.linear(),x.r),minus(x.box,x.c));
  I v=meet(plus(y.linear(),y.r),minus(y.box,y.c));
  z.r=plus(plus(times(x.c,y.r),times(y.c,x.r)),times(u,v));
  z.finish();return z;
}
inline void constrain(T& x,I valid) { x.box=meet(x.box,valid); }
inline void contract(T& x,const T& alternative) {
  I error=plus(alternative.r,minus(alternative.c,x.c));
  double radius=0;
  for(int j=0;j<6;++j) {
    I d=minus(alternative.a[j],x.a[j]);
    if(!zero(d))radius=hetero_bounds::up(radius+std::max(std::abs(d.lo),std::abs(d.hi)));
  }
  error=plus(error,I(-radius,radius));x.r=meet(x.r,error);
  x.box=meet(x.box,alternative.box);x.finish();
}
using Matrix=std::array<T,9>;
inline void stochastic(Matrix& p) {
  std::vector<I> bounds;for(const T& x:p)bounds.push_back(hetero_bounds::prob(x.box));
  hetero_bounds::tighten(bounds,3);
  for(int i=0;i<9;++i)constrain(p[i],bounds[i]);
  // Exact row identity provides a second affine representation of each entry.
  for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
    T other;for(int h=0;h<3;++h)if(h!=j)other=sum(other,p[i*3+h]);
    contract(p[i*3+j],difference(T(I(1)),other));
  }
}
inline Matrix identity() { Matrix p;for(int i=0;i<3;++i)p[i*3+i]=T(I(1));return p; }
inline Matrix matmul(const Matrix& a,const Matrix& b) {
  Matrix out;
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int h=0;h<3;++h)
    out[i*3+j]=sum(out[i*3+j],product(a[i*3+h],b[h*3+j]));
  stochastic(out);return out;
}
inline std::vector<Matrix> transitions(const hetero_bounds::Tree& tree,
    const std::vector<I>& rates,const Rcpp::IntegerMatrix& index,
    const std::vector<std::vector<I>>& original) {
  double lambda=0;
  for(int i=0;i<3;++i) {
    I total;for(int j=0;j<3;++j)if(i!=j)total=hetero_bounds::posadd(total,rates[i*3+j]);
    lambda=std::max(lambda,total.hi);
  }
  if(!(lambda>0)||!std::isfinite(lambda))throw std::overflow_error("shared-Q uniformization scale");
  Matrix jump;
  for(int i=0;i<3;++i) {
    T total;
    for(int j=0;j<3;++j)if(i!=j) {
      jump[i*3+j]=scale(variable(rates[i*3+j],index(i,j)-1),hetero_bounds::divpos(I(1),I(lambda)));
      total=sum(total,jump[i*3+j]);
    }
    jump[i*3+i]=difference(T(I(1)),total);
  }
  stochastic(jump);
  // J is shared by every branch. Compute its 18 powers once per rate box.
  std::array<Matrix,19> powers;powers[0]=identity();
  for(int n=1;n<=18;++n)powers[n]=matmul(powers[n-1],jump);
  std::vector<Matrix> answer;answer.reserve(tree.lengths.size());
  for(size_t l=0;l<tree.lengths.size();++l) {
    I mu=hetero_bounds::posmul(I(lambda),I(tree.lengths[l]));int scaling=0;
    while(mu.hi>.125) {
      mu=I(std::max(0.0,hetero_bounds::dn(std::ldexp(mu.lo,-1))),hetero_bounds::up(std::ldexp(mu.hi,-1)));
      if(++scaling>24)throw std::overflow_error("shared-Q stiff scaling budget");
    }
    Matrix p=powers[0];I coefficient(1);
    for(int n=1;n<=18;++n) {
      coefficient=hetero_bounds::divpos(hetero_bounds::posmul(coefficient,mu),I(n));
      for(int a=0;a<9;++a)p[a]=sum(p[a],scale(powers[n][a],coefficient));
    }
    I next=hetero_bounds::divpos(hetero_bounds::posmul(coefficient,mu),I(19));
    I tail=hetero_bounds::divpos(next,sub(I(1),hetero_bounds::divpos(mu,I(20))));
    I decay=hetero_bounds::negexp(mu);
    for(T& x:p) { x=scale(x,decay);x.r=plus(x.r,I(0,tail.hi));x.box=plus(x.box,I(0,tail.hi));x.finish(); }
    stochastic(p);
    for(int n=0;n<scaling;++n)p=matmul(p,p);
    for(int a=0;a<9;++a)constrain(p[a],original[l][a]);
    answer.push_back(p);
  }
  return answer;
}
struct Likelihood { T value;long long exponent;I log; };
// One binary-node kernel, shared by serial pruning and the subtree-repeat DAG.
// Child, state, multiply/add and scaling order is the original serial order.
inline long long message(T* target,const T* const child[2],
    const long long shifts[2],const Matrix* const p[2]) {
  for(int i=0;i<3;++i)target[i]=T(I(1));
  long long shift=0;
  for(int c=0;c<2;++c) {
    shift+=shifts[c];const auto& transition=*p[c];
    for(int i=0;i<3;++i) {
      T total;I row[3],values[3];
      for(int h=0;h<3;++h) {
        total=sum(total,product(transition[i*3+h],child[c][h]));
        row[h]=transition[i*3+h].box;values[h]=child[c][h].box;
      }
      constrain(total,hetero_bounds::stochastic_dot(row,values,3));
      target[i]=product(target[i],total);
      if(target[i].box.lo<0)constrain(target[i],I(0,hetero_bounds::inf));
    }
  }
  double maximum=0;for(int i=0;i<3;++i)maximum=std::max(maximum,target[i].box.hi);
  int scaling=0;if(maximum>0)std::frexp(maximum,&scaling);
  const double factor=std::ldexp(1.0,-scaling);
  if(!std::isfinite(factor))throw std::overflow_error("shared-Q message scaling");
  for(int i=0;i<3;++i)target[i]=scale(target[i],I(factor));
  return shift+scaling;
}
inline Likelihood root_likelihood(const T* root,long long exponent) {
  T total;for(int i=0;i<3;++i)total=sum(total,root[i]);
  total=scale(total,hetero_bounds::divpos(I(1),I(3)));
  if(total.box.lo<0)constrain(total,I(0,hetero_bounds::inf));
  I log=hetero_bounds::log_interval(total.box);
  I shift=times(I(static_cast<double>(exponent)),hetero_bounds::log2_interval());
  if(log.lo!=-hetero_bounds::inf)log.lo=add(I(log.lo),shift).lo;
  if(log.hi!=-hetero_bounds::inf)log.hi=std::min(0.0,add(I(log.hi),shift).hi);
  return {total,exponent,log};
}
inline std::vector<Likelihood> likelihood(const hetero_bounds::Tree& tree,
    const Rcpp::IntegerMatrix& x,const std::vector<Matrix>& p) {
  std::vector<T> scratch(tree.nn*3);std::vector<long long> exponent(tree.nn);
  std::vector<Likelihood> answer;answer.reserve(x.ncol());
  for(int b=0;b<x.ncol();++b) {
    if((b&15)==0)Rcpp::checkUserInterrupt();
    for(int pos=tree.nn-1;pos>=0;--pos) {
      int node=tree.pre[pos]-1;T* target=&scratch[node*3];
      if(node<tree.nt) {
        for(int i=0;i<3;++i)target[i]=T(I(x(node,b)==NA_INTEGER||x(node,b)==i+1 ? 1 : 0));
        exponent[node]=0;continue;
      }
      const int l=tree.children[node][0],r=tree.children[node][1];
      const T* child[]={&scratch[l*3],&scratch[r*3]};
      const long long shifts[]={exponent[l],exponent[r]};
      const Matrix* matrices[]={&p[tree.length_id[l]],&p[tree.length_id[r]]};
      exponent[node]=message(target,child,shifts,matrices);
    }
    answer.push_back(root_likelihood(&scratch[tree.root*3],exponent[tree.root]));
  }
  return answer;
}
inline I positive_exp(I x) {
  auto point=[](double z) {
    if(z<=0)return hetero_bounds::negexp(I(-z));
    I inverse=hetero_bounds::negexp(I(z));
    if(inverse.lo<=0)throw std::overflow_error("shared-Q weight scaling");
    return hetero_bounds::divpos(I(1),inverse);
  };
  return I(point(x.lo).lo,point(x.hi).hi);
}
inline bool certificate(const T& contrast,const T* constraint,bool upper) {
  auto valid=[&](const T& x) { return upper ? x.box.hi<=0 : x.box.lo>=0; };
  if(valid(contrast))return true;
  if(!constraint)return false;
  // Each nonnegative multiplier is a separate weak-duality certificate.
  // Approximate coefficient cancellation only selects candidates; it never
  // supplies the bound itself. Interval evaluation certifies every candidate.
  for(int j=0;j<6;++j) {
    double a=contrast.a[j].lo/2+contrast.a[j].hi/2;
    double b=constraint->a[j].lo/2+constraint->a[j].hi/2;
    if(b==0)continue;
    double eta=(upper ? -a : a)/b;
    if(!(eta>0)||!std::isfinite(eta)||eta>1e100)continue;
    T dual=sum(contrast,scale(*constraint,I(upper ? eta : -eta)));
    if(valid(dual))return true;
  }
  return false;
}
inline I probability(const T& tail,const T& other,I original,const T* constraint=nullptr) {
  // A zero denominator is left to the original conservative interval method.
  if(!constraint&&plus(tail.box,other.box).lo<=0)return original;
  I result=original;
  double lo=result.lo,hi=result.hi;
  for(int step=0;step<40;++step) {
    double u=lo/2+hi/2;
    T contrast=difference(scale(tail,minus(I(1),I(u))),scale(other,I(u)));
    if(certificate(contrast,constraint,true))hi=u;else lo=u;
  }
  result.hi=hi;lo=result.lo;hi=result.hi;
  for(int step=0;step<40;++step) {
    double u=lo/2+hi/2;
    T contrast=difference(scale(tail,minus(I(1),I(u))),scale(other,I(u)));
    if(certificate(contrast,constraint,false))lo=u;else hi=u;
  }
  result.lo=lo;return result;
}
// Shared by finite-rate and compact transient models. Preserve the operation
// order of the original finite-rate rank certificate (including scaling).
inline I rank(const std::vector<Likelihood>& shared,const std::vector<int>& id,
    const std::vector<I>& proposal,const std::vector<int>& upper_tail,
    double gamma,I original) {
  double pivot=-hetero_bounds::inf;
  for(size_t i=0;i<id.size();++i)
    pivot=std::max(pivot,hetero_bounds::sub(shared[id[i]].log,proposal[i]).hi);
  if(!std::isfinite(pivot))return original;
  T tail,other,observed_weight;
  for(size_t i=0;i<id.size();++i) {
    const auto& model=shared[id[i]];
    if(model.value.box.hi==0)continue;
    I shift=times(I(static_cast<double>(model.exponent)),hetero_bounds::log2_interval());
    I log_scale=hetero_bounds::sub(hetero_bounds::sub(shift,proposal[i]),I(pivot));
    auto weight=scale(model.value,positive_exp(log_scale));
    if(i==0)observed_weight=weight;
    if(upper_tail[i])tail=sum(tail,weight);else other=sum(other,weight);
  }
  if(gamma>0) {
    I threshold=times(I(gamma),positive_exp(I(-pivot)));
    auto constraint=difference(observed_weight,T(threshold));
    return probability(tail,other,original,&constraint);
  }
  return probability(tail,other,original);
}
}
#endif
