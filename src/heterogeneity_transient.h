#ifndef FASTPHYLOSIG_HETEROGENEITY_TRANSIENT_H
#define FASTPHYLOSIG_HETEROGENEITY_TRANSIENT_H
#include "heterogeneity_taylor.h"
#include "heterogeneity_ard_boundary.h"
#include "heterogeneity_repeats.h"

// Three-state multiclass boundaries and their bounded-extra-jump neighbours.
// The two shared coordinates retain the *whole* base transient across branches;
// the remaining four rates have an explicit coupling error, never set to zero.
namespace hetero_transient {
using hetero_bounds::I;
using hetero_taylor::T;
struct Context {
  int family=-1,a=0,b=0;
  double delta=hetero_bounds::inf,error=1;
  I u,rho;
  hetero_taylor::Matrix jump;
  Context(const hetero_bounds::Tree& tree,const std::vector<I>& q) {
    const int edges[6][2]={{1,2},{3,5},{6,7},{1,3},{2,6},{5,7}};
    for(int f=0;f<6;++f) {
      double d=0;
      for(int i=0;i<3;++i) {
        I row;
        for(int j=0;j<3;++j)if(i!=j&&i*3+j!=edges[f][0]&&i*3+j!=edges[f][1])
          row=hetero_bounds::exact_posadd(row,I(q[i*3+j].hi));
        d=std::max(d,row.hi);
      }
      if(d<delta) {delta=d;family=f;a=edges[f][0];b=edges[f][1];}
    }
    if(family<0||!std::isfinite(delta))return;
    I s=hetero_bounds::exact_posadd(q[a],q[b]);
    u=hetero_bounds::negexp(hetero_bounds::posmul(s,I(tree.lengths.front())));
    rho=hetero_bounds::fraction(q[a],q[b]);
    error=hetero_bounds::prob(hetero_bounds::sub(I(1),
      hetero_bounds::negexp(hetero_bounds::posmul(I(delta),I(tree.lengths.back()))))).hi;
    jump=hetero_taylor::identity();
    if(delta>0) {
      int coordinate=2;
      for(int i=0;i<3;++i) {
        T total;
        for(int j=0;j<3;++j)if(i!=j&&i*3+j!=a&&i*3+j!=b) {
          jump[i*3+j]=hetero_taylor::variable(hetero_bounds::prob(
            hetero_bounds::divpos(q[i*3+j],I(delta))),coordinate++);
          total=hetero_taylor::sum(total,jump[i*3+j]);
        }
        jump[i*3+i]=hetero_taylor::difference(T(I(1)),total);
      }
      hetero_taylor::stochastic(jump);
    }
  }
};
// alpha >= 0, x in [0,1]; 0^0=1 is an enclosing endpoint convention.
inline I power(I x,I alpha) {
  if(alpha.hi==0)return I(1);
  auto endpoint=[](double u,double a) {
    if(a==0||u==1)return I(1);
    if(u==0)return I(0);
    return hetero_bounds::exp_nonpositive(hetero_taylor::times(I(a),hetero_bounds::logpoint(u)));
  };
  return I(endpoint(x.lo,alpha.hi).lo,endpoint(x.hi,alpha.lo).hi);
}
inline T decay(const T& u,double t,double t0) {
  if(t==t0)return u;
  I alpha=hetero_bounds::divpos(I(t),I(t0));alpha.lo=std::max(1.0,alpha.lo);
  if(!std::isfinite(alpha.hi))throw std::overflow_error("transient branch ratio");
  const double c=u.c.lo;
  I derivative=hetero_bounds::posmul(alpha,power(u.box,
    I(std::max(0.0,hetero_bounds::dn(alpha.lo-1)),hetero_bounds::up(alpha.hi-1))));
  const double slope=derivative.lo/2+derivative.hi/2;
  T out=hetero_taylor::sum(T(power(I(c),alpha)),
    hetero_taylor::scale(hetero_taylor::difference(u,T(I(c))),I(slope)));
  I remainder=hetero_taylor::times(hetero_taylor::minus(derivative,I(slope)),
    hetero_taylor::minus(u.box,I(c)));
  out.r=hetero_taylor::plus(out.r,remainder);
  out.box=power(u.box,alpha);out.finish();return out;
}
inline I logmean_point(double v) {
  if(v==0)return I(0);if(v==1)return I(1);
  I log=hetero_bounds::logpoint(v);
  return hetero_bounds::prob(hetero_bounds::divpos(hetero_bounds::prob(hetero_bounds::sub(I(1),I(v))),I(-log.hi,-log.lo)));
}
inline I logmean_derivative(double v) {
  I log=hetero_bounds::logpoint(v),x(-log.hi,-log.lo);
  if(x.hi<=.125) {
    I term=hetero_bounds::divpos(I(1),I(2)),sum=term;
    for(int n=1;n<=18;++n) {term=hetero_bounds::divpos(hetero_bounds::posmul(term,x),I(n+2));sum=hetero_bounds::posadd(sum,term);}
    I next=hetero_bounds::divpos(hetero_bounds::posmul(term,x),I(21));
    I tail=hetero_bounds::divpos(next,hetero_bounds::sub(I(1),hetero_bounds::divpos(x,I(22))));
    sum.hi=hetero_bounds::up(sum.hi+tail.hi);return sum;
  }
  I denominator=hetero_bounds::posmul(I(v),hetero_bounds::posmul(x,x));
  if(denominator.lo<=0)return I(.5,hetero_bounds::inf);
  I numerator=hetero_bounds::prob(hetero_bounds::sub(I(1),hetero_bounds::posmul(I(v),hetero_bounds::posadd(I(1),x))));
  return hetero_bounds::divpos(numerator,denominator);
}
inline T logmean(const T& v) {
  const double l=std::max(0.0,v.box.lo),h=std::min(1.0,v.box.hi);
  const I zl=logmean_point(l),zh=logmean_point(h);
  if(l==h)return T(I(zl.lo,zh.hi));
  const I width=hetero_bounds::sub(I(h),I(l));
  if(width.lo<=0)return T(I(zl.lo,zh.hi));
  I slope=hetero_bounds::divpos(hetero_bounds::prob(hetero_bounds::sub(zh,zl)),width);
  T chord=hetero_taylor::sum(T(zl),hetero_taylor::scale(hetero_taylor::difference(v,T(I(l))),slope));
  const double c=l/2+h/2;I derivative=logmean_derivative(c);double remainder=1;
  if(std::isfinite(derivative.hi)) {
    I zc=logmean_point(c);
    I left=hetero_taylor::minus(hetero_taylor::plus(zc,hetero_taylor::times(derivative,hetero_bounds::sub(I(l),I(c)))),zl);
    I right=hetero_taylor::minus(hetero_taylor::plus(zc,hetero_taylor::times(derivative,hetero_bounds::sub(I(h),I(c)))),zh);
    remainder=std::min(1.0,std::max(0.0,std::max(left.hi,right.hi)));
  }
  chord.r=hetero_taylor::plus(chord.r,I(0,remainder));chord.box=I(zl.lo,zh.hi);chord.finish();return chord;
}
inline I poisson_tail2(I mu) {
  I term=hetero_bounds::divpos(hetero_bounds::posmul(mu,mu),I(2)),sum=term;
  for(int n=3;n<=18;++n) {term=hetero_bounds::divpos(hetero_bounds::posmul(term,mu),I(n));sum=hetero_bounds::posadd(sum,term);}
  I next=hetero_bounds::divpos(hetero_bounds::posmul(term,mu),I(19));
  sum=hetero_bounds::posadd(sum,hetero_bounds::divpos(next,hetero_bounds::sub(I(1),hetero_bounds::divpos(mu,I(20)))));
  return hetero_bounds::prob(hetero_bounds::posmul(hetero_bounds::negexp(mu),sum));
}
inline std::vector<hetero_taylor::Matrix> transitions(const hetero_bounds::Tree& tree,
    const Context& ctx,I U,I R) {
  const T u=hetero_taylor::variable(U,0),rho=hetero_taylor::variable(R,1);
  const T other=hetero_taylor::difference(T(I(1)),rho);
  auto E=hetero_taylor::identity();const int i=ctx.a/3,j=ctx.a%3,k=ctx.b%3;
  if(ctx.family<3) {E[i*3+i]=T();E[i*3+j]=rho;E[i*3+k]=other;}
  else {E[i*3+i]=E[j*3+i]=other;E[i*3+j]=E[j*3+j]=rho;}
  hetero_taylor::stochastic(E);
  hetero_taylor::Matrix EJ,JE,EJE;
  if(ctx.delta>0) {EJ=hetero_taylor::matmul(E,ctx.jump);JE=hetero_taylor::matmul(ctx.jump,E);EJE=hetero_taylor::matmul(EJ,E);}
  std::vector<hetero_taylor::Matrix> answer;answer.reserve(tree.lengths.size());
  for(double t:tree.lengths) {
    T v=decay(u,t,tree.lengths.front());
    T moved=hetero_taylor::difference(T(I(1)),v);
    T a=hetero_taylor::product(rho,moved),b=hetero_taylor::product(other,moved);
    auto p=hetero_taylor::identity();
    if(ctx.family<3) {
      p[i*3+i]=v;p[i*3+j]=a;p[i*3+k]=b;
    } else {
      p[i*3+i]=hetero_taylor::difference(T(I(1)),a);p[i*3+j]=a;
      p[j*3+i]=b;p[j*3+j]=hetero_taylor::difference(T(I(1)),b);
    }
    if(ctx.delta>0) {
      T z=logmean(v),d=hetero_taylor::difference(z,v),e=hetero_taylor::difference(
        hetero_taylor::sum(T(I(1)),v),hetero_taylor::scale(z,I(2)));
      hetero_taylor::constrain(d,I(0,1));hetero_taylor::constrain(e,I(0,1));
      I mu=hetero_bounds::posmul(I(ctx.delta),I(t)),decay_b=hetero_bounds::negexp(mu);
      I one=hetero_bounds::posmul(mu,decay_b),tail=poisson_tail2(mu);
      const double error=hetero_bounds::prob(hetero_bounds::sub(I(1),
        hetero_bounds::negexp(hetero_bounds::posmul(I(ctx.delta),I(t))))).hi;
      for(int a=0;a<9;++a) {
        T old=p[a];old.r=hetero_taylor::plus(old.r,I(-error,error));old.box=hetero_taylor::plus(old.box,I(-error,error));old.finish();
        T w=hetero_taylor::sum(hetero_taylor::product(v,ctx.jump[a]),hetero_taylor::sum(
          hetero_taylor::product(d,hetero_taylor::sum(EJ[a],JE[a])),hetero_taylor::product(e,EJE[a])));
        hetero_taylor::constrain(w,I(0,1));
        p[a]=hetero_taylor::sum(hetero_taylor::scale(p[a],decay_b),hetero_taylor::scale(w,one));
        p[a].r=hetero_taylor::plus(p[a].r,I(0,tail.hi));p[a].box=hetero_taylor::plus(p[a].box,I(0,tail.hi));p[a].finish();
        hetero_taylor::contract(p[a],old);
      }
    }
    hetero_taylor::stochastic(p);answer.push_back(p);
  }
  return answer;
}
struct Result { I probability;std::vector<I> likelihood;int cells=0;bool excluded=false; };
inline Result evaluate(const hetero_bounds::Tree& tree,const Rcpp::IntegerMatrix& x,
    const Context& ctx,const std::vector<int>& id,const std::vector<I>& proposal,
    const std::vector<int>& tail,double gamma,I original,const std::vector<I>& old,
    const hetero_taylor::RepeatPlan* repeats=nullptr) {
  struct Cell {I u,rho;};std::vector<Cell> cells{{ctx.u,ctx.rho}};
  while(cells.size()<8) {
    int which=0;double width=-1;
    for(size_t i=0;i<cells.size();++i) {
      double w=std::max(cells[i].u.hi-cells[i].u.lo,cells[i].rho.hi-cells[i].rho.lo);
      if(w>width) {width=w;which=i;}
    }
    Cell left=cells[which],right=left;
    I& l=left.u.hi-left.u.lo>=left.rho.hi-left.rho.lo ? left.u : left.rho;
    I& r=left.u.hi-left.u.lo>=left.rho.hi-left.rho.lo ? right.u : right.rho;
    const double mid=l.lo/2+l.hi/2;
    if(!(mid>l.lo&&mid<l.hi))break;
    l.hi=mid;r.lo=mid;cells[which]=left;cells.push_back(right);
  }
  Result out;out.probability=I(1,0);out.likelihood.assign(old.size(),I(hetero_bounds::inf,-hetero_bounds::inf));
  const I threshold=gamma>0 ? hetero_bounds::add(hetero_bounds::logpoint(gamma),proposal[0]) : I(-hetero_bounds::inf);
  bool physical=false,conditional=false;
  for(const Cell& cell:cells) {
    ++out.cells;
    auto matrices=transitions(tree,ctx,cell.u,cell.rho);
    auto shared=repeats ? repeats->likelihood(tree,matrices) : hetero_taylor::likelihood(tree,x,matrices);
    bool empty=false;
    for(size_t i=0;i<old.size();++i) {
      I ll(std::max(old[i].lo,shared[i].log.lo),std::min(old[i].hi,shared[i].log.hi));
      if(ll.lo>ll.hi) {empty=true;break;}
      shared[i].log=ll;
      I shift=hetero_taylor::times(I(static_cast<double>(shared[i].exponent)),hetero_bounds::log2_interval());
      hetero_taylor::constrain(shared[i].value,hetero_taylor::positive_exp(hetero_bounds::sub(ll,shift)));
    }
    if(empty)continue;
    physical=true;
    // Unconditional likelihood union MUST include cells outside C_gamma.
    for(size_t i=0;i<old.size();++i) {
      out.likelihood[i].lo=std::min(out.likelihood[i].lo,shared[i].log.lo);
      out.likelihood[i].hi=std::max(out.likelihood[i].hi,shared[i].log.hi);
    }
    if(gamma>0&&shared[0].log.hi<threshold.lo)continue;
    conditional=true;
    I p=hetero_taylor::rank(shared,id,proposal,tail,gamma,original);
    out.probability.lo=std::min(out.probability.lo,p.lo);
    out.probability.hi=std::max(out.probability.hi,p.hi);
  }
  if(!physical)Rcpp::stop("empty transient cover for a nonempty rate box");
  if(!conditional) {out.probability=I(0);out.excluded=true;}
  return out;
}
}
#endif
