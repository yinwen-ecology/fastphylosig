#ifndef FASTPHYLOSIG_HETEROGENEITY_BOUNDS_H
#define FASTPHYLOSIG_HETEROGENEITY_BOUNDS_H
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Positive interval arithmetic for nuisance-domain certificates, deliberately
// separate from the unchanged ACE objective and Delta estimator. No fast-math.
namespace hetero_bounds {
const double inf = std::numeric_limits<double>::infinity();
inline double dn(double x) { return std::nextafter(x,-inf); }
inline double up(double x) { return std::nextafter(x,inf); }
// Sequential inverse-CDF sampling reaches category j only after every prior
// threshold was crossed. Clamp and cummax also cover a rounded CDF above one
// or a final threshold slightly below its predecessor.
inline double sampler_cdf(double raw,double previous) {
  return std::max(previous,std::min(1.0,std::max(0.0,raw)));
}
struct I {
  double lo, hi;
  I(double a=0) : lo(a),hi(a) {}
  I(double a,double b) : lo(a),hi(b) {}
};
inline I add(I a,I b) { return I(dn(a.lo+b.lo),up(a.hi+b.hi)); }
inline I sub(I a,I b) { return I(dn(a.lo-b.hi),up(a.hi-b.lo)); }
inline I mul(I a,I b) {
  double x[]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
  return I(dn(*std::min_element(x,x+4)),up(*std::max_element(x,x+4)));
}
inline I posadd(I a,I b) { I x=add(a,b);x.lo=std::max(0.0,x.lo);return x; }
inline I posmul(I a,I b) {
  if(a.hi==0 || b.hi==0)return I(0);
  return I(std::max(0.0,dn(a.lo*b.lo)),up(a.hi*b.hi));
}
inline I divpos(I a,I b) {
  if(b.lo<=0)Rcpp::stop("non-positive interval divisor");
  return posmul(a,I(std::max(0.0,dn(1/b.hi)),up(1/b.lo)));
}
inline I prob(I x) { return I(std::max(0.0,x.lo),std::min(1.0,x.hi)); }
inline I negexp_point(double x) {
  if(x==0)return I(1);
  if(x==inf)return I(0);
  if(!(x>=0))Rcpp::stop("negative exponential argument");
  int e=0;std::frexp(x,&e);
  const int s=std::max(0,e+3);
  const double z=std::ldexp(x,-s); // exact power-of-two scaling
  I term(1),sum(1);
  for(int n=1;n<=18;++n) {
    term=divpos(posmul(term,I(z)),I(n));sum=posadd(sum,term);
  }
  I next=divpos(posmul(term,I(z)),I(19));
  I tail=divpos(next,sub(I(1),divpos(I(z),I(20))));
  I out=divpos(I(1),I(sum.lo,up(sum.hi+tail.hi)));
  for(int j=0;j<s;++j) {
    out=prob(posmul(out,out));
    if(out.hi<=std::numeric_limits<double>::denorm_min())return I(0,out.hi);
  }
  return out;
}
inline I negexp(I x) {
  return I(negexp_point(x.hi).lo,negexp_point(std::max(0.0,x.lo)).hi);
}
inline I log_mantissa(double m) {
  // m in [1,2]. log(m)=2*atanh((m-1)/(m+1)); positive remainder.
  if(m==1)return I(0);
  I z=divpos(prob(sub(I(m),I(1))),add(I(m),I(1)));
  I zz=posmul(z,z),power=z,sum(0);
  for(int n=0;n<32;++n) {
    sum=posadd(sum,divpos(power,I(2*n+1)));power=posmul(power,zz);
  }
  I tail=divpos(posmul(I(2),power),mul(I(65),sub(I(1),zz)));
  I out=posmul(I(2),sum);out.hi=up(out.hi+tail.hi);return out;
}
inline I log2_interval() { static const I x=log_mantissa(2);return x; }
inline I logpoint(double x) {
  if(x==0)return I(-inf);
  if(x==inf)return I(inf);
  if(!(x>0))Rcpp::stop("non-positive logarithm argument");
  int e=0;double m=std::frexp(x,&e)*2;--e;
  return add(log_mantissa(m),mul(I(e),log2_interval()));
}
inline I log_interval(I x) {
  return I(x.lo==0 ? -inf : logpoint(x.lo).lo,
           x.hi==0 ? -inf : logpoint(x.hi).hi);
}
inline I exp_nonpositive(I x) {
  if(x.hi==-inf)return I(0);
  return negexp(I(std::max(0.0,-x.hi),x.lo==-inf ? inf : std::max(0.0,-x.lo)));
}
inline std::vector<I> matmul(const std::vector<I>& a,const std::vector<I>& b,int k) {
  std::vector<I> c(k*k);
  for(int i=0;i<k;++i)for(int j=0;j<k;++j) {
    I z(0);for(int h=0;h<k;++h)z=posadd(z,posmul(a[i*k+h],b[h*k+j]));
    c[i*k+j]=prob(z);
  }
  return c;
}
inline void tighten(std::vector<I>& a,int k) {
  for(int i=0;i<k;++i)for(int j=0;j<k;++j) {
    I other(0);for(int h=0;h<k;++h)if(h!=j)other=posadd(other,a[i*k+h]);
    I complement=prob(sub(I(1),other));
    a[i*k+j].lo=std::max(a[i*k+j].lo,complement.lo);
    a[i*k+j].hi=std::min(a[i*k+j].hi,complement.hi);
  }
}
// Exact box/simplex linear-program extrema, with outward arithmetic. Narrow
// point intervals retain the original inexpensive dot-product evaluation.
inline I stochastic_dot(const I* p,const I* x,int k,bool constrain=true) {
  I naive(0);double width=0;
  for(int j=0;j<k;++j) { naive=posadd(naive,posmul(p[j],x[j]));width=std::max(width,p[j].hi-p[j].lo); }
  if(!constrain||width<1e-9)return naive;
  I limits[2];
  for(int side=0;side<2;++side) {
    int order[8];for(int j=0;j<k;++j)order[j]=j;
    auto value=[&](int j){return side ? x[j].hi : x[j].lo;};
    for(int j=1;j<k;++j) {
      int id=order[j],h=j;
      while(h>0&&(side ? value(id)>value(order[h-1]) : value(id)<value(order[h-1]))) {
        order[h]=order[h-1];--h;
      }
      order[h]=id;
    }
    I total(0),used(0);
    for(int j=0;j<k;++j) {
      used=posadd(used,I(p[j].lo));total=posadd(total,posmul(I(p[j].lo),I(value(j))));
    }
    I remaining=prob(sub(I(1),used));
    for(int h=0;h<k;++h) {
      int j=order[h];I cap=prob(sub(I(p[j].hi),I(p[j].lo)));
      I take(std::min(remaining.lo,cap.lo),std::min(remaining.hi,cap.hi));
      total=posadd(total,posmul(take,I(value(j))));
      remaining=prob(sub(remaining,take));
    }
    limits[side]=total;
  }
  return I(std::max(naive.lo,limits[0].lo),std::min(naive.hi,limits[1].hi));
}
inline std::vector<I> er_transition(double t,int k,double lower,double upper) {
  I decay;
  if(lower==inf)decay=I(0);
  else {
    I x=posmul(I(k),posmul(I(t),I(lower,upper)));
    decay=negexp(x);
  }
  std::vector<I> out(k*k);
  I off=divpos(prob(sub(I(1),decay)),I(k));
  for(int i=0;i<k;++i)for(int j=0;j<k;++j)
    out[i*k+j]=i==j ? prob(posadd(off,decay)) : off;
  // off+decay has dependency; the equivalent diagonal formula is sharper.
  I diagonal=divpos(posadd(I(1),posmul(I(k-1),decay)),I(k));
  for(int i=0;i<k;++i)out[i*k+i]=prob(diagonal);
  tighten(out,k);return out;
}
inline std::vector<I> ard_transition(double t,int k,const std::vector<I>& rates) {
  std::vector<I> wide(k*k,I(0,1));
  bool all_zero=true;
  for(const I& rate:rates)if(rate.hi!=0)all_zero=false;
  if(all_zero) { for(int i=0;i<k*k;++i)wide[i]=I(i/k==i%k ? 1 : 0);return wide; }
  double lambda=0;
  for(int i=0;i<k;++i) {
    I total(0);for(int j=0;j<k;++j)if(i!=j)total=posadd(total,rates[i*k+j]);
    lambda=std::max(lambda,total.hi);
  }
  if(!std::isfinite(lambda))return wide;
  if(lambda==0) { for(int i=0;i<k*k;++i)wide[i]=I(i/k==i%k ? 1 : 0);return wide; }
  I mu=posmul(I(lambda),I(t));
  if(!std::isfinite(mu.hi))return wide;
  int scaling=0;
  while(mu.hi>.125) {
    // The upper/lower endpoints can have very different magnitudes. Retain
    // outward bounds even if halving the lower endpoint enters subnormals.
    mu=I(std::max(0.0,dn(std::ldexp(mu.lo,-1))),up(std::ldexp(mu.hi,-1)));
    ++scaling;
  }
  std::vector<I> jump(k*k),power(k*k),sum(k*k);
  for(int i=0;i<k;++i) {
    I total(0);
    for(int j=0;j<k;++j)if(i!=j) {
      jump[i*k+j]=prob(divpos(rates[i*k+j],I(lambda)));
      total=posadd(total,jump[i*k+j]);
    }
    jump[i*k+i]=prob(sub(I(1),total));
    power[i*k+i]=sum[i*k+i]=I(1);
  }
  tighten(jump,k);I coefficient(1);
  for(int n=1;n<=18;++n) {
    power=matmul(power,jump,k);
    coefficient=divpos(posmul(coefficient,mu),I(n));
    for(int a=0;a<k*k;++a)sum[a]=posadd(sum[a],posmul(coefficient,power[a]));
  }
  I next=divpos(posmul(coefficient,mu),I(19));
  I tail=divpos(next,sub(I(1),divpos(mu,I(20))));
  I decay=negexp(mu);
  for(int a=0;a<k*k;++a) {
    sum[a]=prob(posmul(decay,sum[a]));sum[a].hi=std::min(1.0,up(sum[a].hi+tail.hi));
  }
  tighten(sum,k);
  for(int j=0;j<scaling;++j) { sum=matmul(sum,sum,k);tighten(sum,k); }
  return sum;
}

struct Tree {
  int nt,nn,k,root;
  double total_length=0;
  std::vector<int> pre,parent,length_id;
  std::vector<std::vector<int>> children;
  std::vector<double> lengths;
  mutable std::vector<I> scratch;
  mutable std::vector<long long> exponent;
  Tree(const Rcpp::List& core,int states):k(states) {
    nt=Rcpp::as<int>(core["n_tip"]);nn=Rcpp::as<int>(core["n_total"]);
    root=Rcpp::as<int>(core["root"])-1;
    pre=Rcpp::as<std::vector<int>>(core["preorder"]);
    parent=Rcpp::as<std::vector<int>>(core["parent"]);
    Rcpp::NumericVector branch=core["branch_length_by_node"];
    if(k<2||k>8||nt<2||nn<=nt||root<0||root>=nn||int(pre.size())!=nn||
       int(parent.size())!=nn||branch.size()!=nn||pre[0]!=root+1)
      Rcpp::stop("invalid nuisance tree/state workspace");
    children.resize(nn);length_id.assign(nn,-1);std::vector<bool> seen(nn,false);
    std::vector<int> sorted;
    for(int pos=0;pos<nn;++pos) {
      int node=pre[pos]-1;
      if(node<0||node>=nn||seen[node])Rcpp::stop("invalid nuisance preorder");
      seen[node]=true;
      if(!pos) { if(parent[node]!=0)Rcpp::stop("invalid nuisance root");continue; }
      int p=parent[node]-1;
      if(p<0||p>=nn||!seen[p]||!(branch[node]>0)||!std::isfinite(branch[node]))
        Rcpp::stop("nuisance tests require strictly positive finite branches");
      total_length+=branch[node];
      children[p].push_back(node);sorted.push_back(node);
    }
    for(int node=0;node<nn;++node)
      if(children[node].size()!=static_cast<size_t>(node<nt ? 0 : 2))
        Rcpp::stop("nuisance tests require a binary tree");
    std::sort(sorted.begin(),sorted.end(),[&](int a,int b){return branch[a]<branch[b];});
    for(int node:sorted) {
      if(lengths.empty()||lengths.back()!=branch[node])lengths.push_back(branch[node]);
      length_id[node]=lengths.size()-1;
    }
    scratch.resize(nn*k);exponent.resize(nn);
  }
  void validate(const Rcpp::IntegerMatrix& x) const {
    if(x.nrow()!=nt)Rcpp::stop("nuisance patterns do not match the tree");
    for(int state:x)if(state!=NA_INTEGER&&(state<1||state>k))
      Rcpp::stop("nuisance states are outside the declared state space");
  }
  std::vector<I> likelihood(const Rcpp::IntegerMatrix& x,
      const std::vector<std::vector<I>>& transition,bool constrain=true) const {
    validate(x);std::vector<I> answer(x.ncol());
    for(int b=0;b<x.ncol();++b) {
      if((b&31)==0)Rcpp::checkUserInterrupt();
      for(int pos=nn-1;pos>=0;--pos) {
        int node=pre[pos]-1;I* target=&scratch[node*k];
        if(node<nt) {
          for(int i=0;i<k;++i)target[i]=I(x(node,b)==NA_INTEGER||x(node,b)==i+1 ? 1 : 0);
          exponent[node]=0;continue;
        }
        for(int i=0;i<k;++i)target[i]=I(1);
        long long shift=0;
        for(int child:children[node]) {
          shift+=exponent[child];const auto& p=transition[length_id[child]];
          for(int i=0;i<k;++i) {
            I total=stochastic_dot(&p[i*k],&scratch[child*k],k,constrain);
            target[i]=posmul(target[i],total);
          }
        }
        double maximum=0;for(int i=0;i<k;++i)maximum=std::max(maximum,target[i].hi);
        int scale=0;if(maximum>0)std::frexp(maximum,&scale);
        for(int i=0;i<k;++i)target[i]=I(std::max(0.0,dn(std::ldexp(target[i].lo,-scale))),
                                     up(std::ldexp(target[i].hi,-scale)));
        exponent[node]=shift+scale;
      }
      I total(0);for(int i=0;i<k;++i)total=posadd(total,scratch[root*k+i]);
      total=divpos(total,I(k));I value=log_interval(total);
      I shift=mul(I(static_cast<double>(exponent[root])),log2_interval());
      if(value.lo!=-inf)value.lo=add(I(value.lo),shift).lo;
      if(value.hi!=-inf)value.hi=add(I(value.hi),shift).hi;
      value.hi=std::min(0.0,value.hi);answer[b]=value;
    }
    return answer;
  }
};
inline Rcpp::NumericMatrix matrix(const std::vector<I>& values) {
  Rcpp::NumericMatrix out(values.size(),2);
  for(size_t i=0;i<values.size();++i) { out(i,0)=values[i].lo;out(i,1)=values[i].hi; }
  Rcpp::colnames(out)=Rcpp::CharacterVector::create("lower","upper");return out;
}
}
#endif
