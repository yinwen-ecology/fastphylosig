#ifndef FASTPHYLOSIG_ARD_BOUNDARY_H
#define FASTPHYLOSIG_ARD_BOUNDARY_H
#include "heterogeneity_bounds.h"
#include <array>

namespace hetero_bounds {
inline I exact_posadd(I a,I b) {
  if(a.hi==0)return b;
  if(b.hi==0)return a;
  return posadd(a,b);
}
struct Forest3 { int mask;std::array<int,3> root; };
inline const std::vector<Forest3>& forests3() {
  static const auto forests=[]() {
    const int from[]={0,0,1,1,2,2},to[]={1,2,0,2,0,1};
    std::vector<Forest3> out;
    for(int mask=0;mask<64;++mask) {
      int next[]={-1,-1,-1},count=0;bool valid=true;
      for(int e=0;e<6;++e)if(mask&(1<<e)) {
        if(next[from[e]]>=0)valid=false;
        next[from[e]]=to[e];++count;
      }
      if(!valid||count>2)continue;
      Forest3 f;f.mask=mask;
      for(int i=0;i<3;++i) {
        int node=i,steps=0;
        while(next[node]>=0&&steps<3) { node=next[node];++steps; }
        if(steps==3)valid=false;
        f.root[i]=node;
      }
      if(valid)out.push_back(f);
    }
    return out;
  }();
  return forests;
}
inline I fraction(I n,I other) {
  if(n.hi==0)return I(0);
  if(other.hi==0)return I(1);
  double dlo=std::max(0.0,dn(n.hi+other.lo)),dhi=up(n.lo+other.hi);
  return I(dhi>0 ? std::max(0.0,dn(n.lo/dhi)) : 0,
    dlo>0 ? std::min(1.0,up(n.hi/dlo)) : 1);
}
inline std::vector<I> resolvent3(double t,const std::vector<I>& rates,int power) {
  // Endpoint projective pairs for h*q. Each pair is independently normalized;
  // its rational interpolation still covers the complete endpoint interval.
  I a[6][2],b[6][2];int edge=0;
  double h=std::ldexp(t,-power);
  I step(std::max(0.0,dn(h)),up(h));
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)if(i!=j) {
    for(int side=0;side<2;++side) {
      double q=side ? rates[i*3+j].hi : rates[i*3+j].lo;
      if(q==0) { a[edge][side]=I(0);b[edge][side]=I(1); }
      else if(q==inf) { a[edge][side]=I(1);b[edge][side]=I(0); }
      else {
        I r=posmul(step,I(q));
        if(r.hi<=1) { a[edge][side]=r;b[edge][side]=I(1); }
        else if(r.lo>0) { a[edge][side]=I(1);b[edge][side]=divpos(I(1),r); }
        else { a[edge][side]=b[edge][side]=I(0,1); }
      }
    }
    ++edge;
  }
  std::vector<I> out(9,I(1,0));
  for(int corner=0;corner<64;++corner) {
    I numer[9];I total(0);
    for(const auto& forest:forests3()) {
      I weight(1);
      for(int e=0;e<6;++e) {
        int side=(corner>>e)&1;
        weight=posmul(weight,(forest.mask&(1<<e)) ? a[e][side] : b[e][side]);
      }
      total=exact_posadd(total,weight);
      for(int i=0;i<3;++i)numer[i*3+forest.root[i]]=exact_posadd(numer[i*3+forest.root[i]],weight);
    }
    if(total.hi==0)continue; // identically zero Bernstein denominator coefficient
    for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
      I other(0);for(int l=0;l<3;++l)if(l!=j)other=exact_posadd(other,numer[i*3+l]);
      I f=fraction(numer[i*3+j],other);
      out[i*3+j].lo=std::min(out[i*3+j].lo,f.lo);
      out[i*3+j].hi=std::max(out[i*3+j].hi,f.hi);
    }
  }
  tighten(out,3);return out;
}
inline std::vector<I> stochastic_square3(const std::vector<I>& a) {
  std::vector<I> out(9);I column[3];
  for(int j=0;j<3;++j) {
    for(int h=0;h<3;++h)column[h]=a[h*3+j];
    for(int i=0;i<3;++i)out[i*3+j]=prob(stochastic_dot(&a[i*3],column,3));
  }
  tighten(out,3);return out;
}
inline I occupation(double t,I exit,I incoming,int initial) {
  if(exit.hi==inf||incoming.hi==inf)return I(0,1);
  I sum=exact_posadd(exit,incoming);if(sum.hi==0)return I(initial);
  if(sum.lo<=0)return I(0,1);
  I decay=negexp(posmul(sum,I(t)));
  return prob(posadd(posmul(I(initial),decay),
    posmul(divpos(incoming,sum),prob(sub(I(1),decay)))));
}
inline void occupation_tighten3(std::vector<I>& out,double t,const std::vector<I>& rates) {
  for(int j=0;j<3;++j) {
    I exits(0);double binlo=inf,binhi=0;
    for(int h=0;h<3;++h)if(h!=j) {
      exits=exact_posadd(exits,rates[j*3+h]);
      binlo=std::min(binlo,rates[h*3+j].lo);binhi=std::max(binhi,rates[h*3+j].hi);
    }
    for(int i=0;i<3;++i) {
      I lower=exits.hi==inf ? I(0) : occupation(t,I(exits.hi),I(binlo),i==j);
      I upper=binhi==inf ? I(1) : occupation(t,I(exits.lo),I(binhi),i==j);
      out[i*3+j].lo=std::max(out[i*3+j].lo,lower.lo);
      out[i*3+j].hi=std::min(out[i*3+j].hi,upper.hi);
    }
  }
  for(int i=0;i<3;++i) {
    I exits(0);for(int j=0;j<3;++j)if(i!=j)exits=exact_posadd(exits,rates[i*3+j]);
    I stay=negexp(posmul(I(exits.hi),I(t)));
    out[i*3+i].lo=std::max(out[i*3+i].lo,stay.lo);
    for(int j=0;j<3;++j)if(i!=j)out[i*3+j].hi=std::min(out[i*3+j].hi,prob(sub(I(1),stay)).hi);
  }
  tighten(out,3);
}
inline std::vector<I> ard_box_transition(double t,const std::vector<I>& rates,int power=24) {
  bool unbounded=false;for(I r:rates)if(r.hi==inf)unbounded=true;
  if(!unbounded)return ard_transition(t,3,rates);
  bool unrestricted=true;
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)if(i!=j)
    if(rates[i*3+j].lo!=0||rates[i*3+j].hi!=inf)unrestricted=false;
  if(unrestricted)return std::vector<I>(9,I(0,1));
  auto out=resolvent3(t,rates,power);
  for(int s=0;s<power;++s) {
    bool wide=true;for(I x:out)if(x.lo!=0||x.hi!=1)wide=false;
    if(wide)break;
    out=stochastic_square3(out);
  }
  const double error=std::ldexp(32.0,-power);
  for(I& x:out)x=prob(I(dn(x.lo-error),up(x.hi+error)));
  occupation_tighten3(out,t,rates);return out;
}
}
#endif
