#ifndef FASTPHYLOSIG_FOREST_SHARED_H
#define FASTPHYLOSIG_FOREST_SHARED_H
#include "heterogeneity_ard_boundary.h"
#include "heterogeneity_taylor.h"

// A shared maximum-forest projection, with a coupling certificate for every
// terminal branch. It never assumes that finite branches are at stationarity.
namespace hetero_forest {
using hetero_bounds::I;
using hetero_taylor::plus;
using hetero_taylor::minus;
using hetero_taylor::times;
using Poly=std::vector<I>;
inline I positive(I x) { return I(std::max(0.0,x.lo),std::max(0.0,x.hi)); }
inline bool eligible(const Rcpp::IntegerMatrix& x) {
  int n=0;
  for(int i=0;i<x.nrow();++i)if(x(i,0)!=NA_INTEGER)++n;
  if(n<2||double(x.nrow())*x.ncol()>2000000)return false;
  for(int b=1;b<x.ncol();++b)for(int i=0;i<x.nrow();++i)
    if((x(i,b)==NA_INTEGER)!=(x(i,0)==NA_INTEGER))return false;
  return true;
}
inline std::array<I,3> stationary(const std::vector<I>& rates) {
  I a[6][2],b[6][2];int e=0;
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)if(i!=j) {
    for(int side=0;side<2;++side) {
      const double q=side ? rates[i*3+j].hi : rates[i*3+j].lo;
      if(q==hetero_bounds::inf) { a[e][side]=I(1);b[e][side]=I(0); }
      else if(q<=1) { a[e][side]=I(q);b[e][side]=I(1); }
      else { a[e][side]=I(1);b[e][side]=hetero_bounds::divpos(I(1),I(q)); }
    }
    ++e;
  }
  std::array<I,3> out{{I(1,0),I(1,0),I(1,0)}};
  for(int corner=0;corner<64;++corner) {
    I numer[3];
    for(const auto& f:hetero_bounds::forests3()) {
      if(!f.mask||!(f.mask&(f.mask-1)))continue;
      I weight(1);
      for(int j=0;j<6;++j)weight=hetero_bounds::posmul(weight,
        f.mask&(1<<j) ? a[j][(corner>>j)&1] : b[j][(corner>>j)&1]);
      numer[f.root[0]]=hetero_bounds::exact_posadd(numer[f.root[0]],weight);
    }
    if(numer[0].hi==0&&numer[1].hi==0&&numer[2].hi==0)continue;
    for(int j=0;j<3;++j) {
      I other(0);for(int k=0;k<3;++k)if(k!=j)other=hetero_bounds::exact_posadd(other,numer[k]);
      I value=hetero_bounds::fraction(numer[j],other);
      out[j].lo=std::min(out[j].lo,value.lo);out[j].hi=std::max(out[j].hi,value.hi);
    }
  }
  for(int j=0;j<3;++j) {
    I other;for(int k=0;k<3;++k)if(k!=j)other=plus(other,out[k]);
    out[j]=hetero_taylor::meet(out[j],hetero_bounds::prob(minus(I(1),other)));
  }
  return out;
}
struct Context {
  double beta=0,eps=1;
  std::vector<I> errors;
  std::array<I,3> pi;
  Context(const hetero_bounds::Tree& tree,const std::vector<I>& rates,
      const Rcpp::IntegerMatrix& x) {
    I sum;
    for(int j=0;j<3;++j) {
      double incoming=hetero_bounds::inf;
      for(int i=0;i<3;++i)if(i!=j)incoming=std::min(incoming,rates[i*3+j].lo);
      sum=hetero_bounds::exact_posadd(sum,I(incoming));
    }
    double pairwise=hetero_bounds::inf;
    for(int i=0;i<3;++i)for(int j=i+1;j<3;++j) {
      I coalesce=hetero_bounds::exact_posadd(I(rates[i*3+j].lo),I(rates[j*3+i].lo));
      for(int k=0;k<3;++k)if(k!=i&&k!=j)coalesce=hetero_bounds::exact_posadd(coalesce,
        I(std::min(rates[i*3+k].lo,rates[j*3+k].lo)));
      pairwise=std::min(pairwise,coalesce.lo);
    }
    beta=std::max(sum.lo,pairwise);if(!(beta>0))return;
    errors.resize(tree.nt);eps=0;std::vector<I> decay(tree.lengths.size(),I(-1));
    for(int i=0;i<tree.nt;++i)if(x(i,0)!=NA_INTEGER) {
      int id=tree.length_id[i];
      if(decay[id].lo<0)decay[id]=hetero_bounds::negexp(hetero_bounds::posmul(I(beta),I(tree.lengths[id])));
      errors[i]=decay[id];
      eps=std::max(eps,errors[i].hi);
    }
    pi=stationary(rates);
  }
};
struct Grid {
  int n;
  std::vector<std::array<int,3>> alpha;
  std::vector<I> inverse_multinomial;
  std::vector<std::array<int,6>> neighbor;
  std::vector<std::array<I,6>> probability;
  explicit Grid(int degree):n(degree) {
    for(int a=0;a<=n;++a)for(int b=0;b<=n-a;++b)alpha.push_back({a,b,n-a-b});
    for(const auto& x:alpha) {
      I inverse(1);
      for(int k=1;k<=x[0];++k)inverse=hetero_bounds::posmul(inverse,
        hetero_bounds::divpos(I(k),I(n-x[0]+k)));
      for(int k=1;k<=x[1];++k)inverse=hetero_bounds::posmul(inverse,
        hetero_bounds::divpos(I(k),I(n-x[0]-x[1]+k)));
      inverse_multinomial.push_back(inverse);
      std::array<int,6> at;std::array<I,6> p;int edge=0;
      for(int i=0;i<3;++i)for(int j=0;j<3;++j)if(i!=j) {
        auto y=x;--y[i];++y[j];at[edge]=x[i]>0 ? id(y) : -1;
        p[edge++]=x[i]>0 ? hetero_bounds::divpos(I(x[i]),I(3*n)) : I(0);
      }
      neighbor.push_back(at);probability.push_back(p);
    }
  }
  int id(const std::array<int,3>& x) const { return x[0]*(n+1)-x[0]*(x[0]-1)/2+x[1]; }
  Poly jump(const Poly& b) const {
    double maximum=0,minimum=hetero_bounds::inf;
    for(I z:b) { maximum=std::max(maximum,z.hi);minimum=std::min(minimum,z.lo); }
    Poly out(b.size());const I third=hetero_bounds::divpos(I(1),I(3));
    for(size_t i=0;i<b.size();++i) {
      I z=hetero_bounds::posmul(b[i],third);
      for(int e=0;e<6;++e)if(neighbor[i][e]>=0)
        z=hetero_bounds::exact_posadd(z,hetero_bounds::posmul(probability[i][e],b[neighbor[i][e]]));
      out[i]=hetero_taylor::meet(z,I(minimum,maximum));
    }
    return out;
  }
  Poly upper(const Poly& b,double eps) const {
    I mu=positive(times(I(n),hetero_bounds::log_interval(
      hetero_bounds::posadd(I(1),hetero_bounds::posmul(I(3),I(eps))))));
    Poly sum=b,power=b;I coefficient(1);double maximum=0;
    for(I z:b)maximum=std::max(maximum,z.hi);
    if(maximum==0)return b;
    for(int k=0;k<256;++k) {
      I next=hetero_bounds::divpos(hetero_bounds::posmul(coefficient,mu),I(k+1));
      if(mu.hi<k+2) {
        I tail=hetero_bounds::divpos(next,minus(I(1),hetero_bounds::divpos(mu,I(k+2))));
        if(tail.hi<=1e-18) {
          I remainder=hetero_bounds::posmul(tail,I(maximum));
          for(I& z:sum)z=hetero_bounds::posadd(z,I(0,remainder.hi));
          return sum;
        }
      }
      coefficient=next;power=jump(power);
      for(size_t i=0;i<sum.size();++i)sum[i]=hetero_bounds::exact_posadd(sum[i],
        hetero_bounds::posmul(coefficient,power[i]));
    }
    throw std::overflow_error("forest polynomial uniformization budget");
  }
};
struct Summary {
  Rcpp::IntegerMatrix x;
  std::vector<std::array<int,3>> counts;
  std::vector<I> profile;
  int n=0;
  explicit Summary(const Rcpp::IntegerMatrix& patterns):x(patterns) {
    for(int i=0;i<x.nrow();++i)if(x(i,0)!=NA_INTEGER)++n;
    for(int p=0;p<x.ncol();++p) {
      std::array<int,3> c{};for(int i=0;i<x.nrow();++i)if(x(i,p)!=NA_INTEGER)++c[x(i,p)-1];
      counts.push_back(c);I value;
      for(int j=0;j<3;++j)if(c[j])value=plus(value,times(I(c[j]),
        hetero_bounds::log_interval(hetero_bounds::divpos(I(c[j]),I(n)))));
      profile.push_back(value);
    }
  }
  std::vector<I> likelihood(const Context& ctx) const {
    I growth=times(I(n),hetero_bounds::log_interval(
      hetero_bounds::posadd(I(1),hetero_bounds::posmul(I(3),I(ctx.eps)))));
    std::vector<std::array<I,3>> logp(x.nrow());
    for(int i=0;i<x.nrow();++i)if(x(i,0)!=NA_INTEGER)for(int j=0;j<3;++j)
      logp[i][j]=hetero_bounds::log_interval(hetero_bounds::prob(
        plus(ctx.pi[j],I(-ctx.errors[i].hi,ctx.errors[i].hi))));
    std::vector<I> out;out.reserve(x.ncol());
    for(int p=0;p<x.ncol();++p) {
      double lo=0,hi=0;
      for(int i=0;i<x.nrow();++i)if(x(i,p)!=NA_INTEGER) {
        I log=logp[i][x(i,p)-1];
        if(log.lo==-hetero_bounds::inf)lo=-hetero_bounds::inf;
        else if(std::isfinite(lo))lo=plus(I(lo),I(log.lo)).lo;
        if(log.hi==-hetero_bounds::inf)hi=-hetero_bounds::inf;
        else if(std::isfinite(hi))hi=plus(I(hi),I(log.hi)).hi;
      }
      out.emplace_back(lo,std::min(0.0,std::min(hi,plus(profile[p],growth).hi)));
    }
    return out;
  }
};
struct Bank {
  Grid grid;
  const Rcpp::IntegerMatrix& x;
  const std::vector<std::array<int,3>>& counts;
  std::vector<I> a_weight,b_weight;
  Poly a,b,observed;
  I gamma_scale;
  Bank(const Summary& summary,const std::vector<int>& ids,
      const std::vector<I>& proposal,const std::vector<int>& tail):
      grid(summary.n),x(summary.x),counts(summary.counts),
      a_weight(x.ncol()),b_weight(x.ncol()),a(grid.alpha.size()),b(a.size()),observed(a.size()) {
    double pivot=-hetero_bounds::inf;for(I p:proposal)pivot=std::max(pivot,-p.lo);
    gamma_scale=hetero_taylor::positive_exp(I(-pivot));I obs_weight;
    for(size_t i=0;i<ids.size();++i) {
      I weight=hetero_bounds::exp_nonpositive(minus(I(-pivot),proposal[i]));
      I& target=tail[i] ? a_weight[ids[i]] : b_weight[ids[i]];
      target=hetero_bounds::exact_posadd(target,weight);if(!i)obs_weight=weight;
    }
    for(int p=0;p<x.ncol();++p) {
      int id=grid.id(counts[p]);
      a[id]=hetero_bounds::exact_posadd(a[id],hetero_bounds::posmul(a_weight[p],grid.inverse_multinomial[id]));
      b[id]=hetero_bounds::exact_posadd(b[id],hetero_bounds::posmul(b_weight[p],grid.inverse_multinomial[id]));
      if(!p)observed[id]=hetero_bounds::posmul(obs_weight,grid.inverse_multinomial[id]);
    }
  }
  std::array<Poly,5> polynomials(const Context& ctx,double gamma) const {
    std::array<Poly,5> out{{grid.upper(a,ctx.eps),a,grid.upper(b,ctx.eps),b,grid.upper(observed,ctx.eps)}};
    for(int p=0;p<x.ncol();++p) {
      I errors[3];
      for(int i=0;i<x.nrow();++i)if(x(i,p)!=NA_INTEGER)
        errors[x(i,p)-1]=hetero_bounds::exact_posadd(errors[x(i,p)-1],I(ctx.errors[i].hi));
      for(int j=0;j<3;++j)if(counts[p][j])for(int k=0;k<3;++k) {
        auto c=counts[p];--c[j];++c[k];int id=grid.id(c);
        I term=hetero_bounds::posmul(errors[j],grid.inverse_multinomial[id]);
        out[1][id]=minus(out[1][id],hetero_bounds::posmul(a_weight[p],term));
        out[3][id]=minus(out[3][id],hetero_bounds::posmul(b_weight[p],term));
      }
    }
    I threshold=hetero_bounds::posmul(I(gamma),gamma_scale);
    for(I& c:out[4])c=minus(c,threshold);
    return out;
  }
};
// Feasibility of all coefficient inequalities c_i+eta*h_i<=0, eta>=0.
// Interval division only selects a constant. Every inequality is rechecked.
inline bool certificate(const Poly& c,const Poly& h) {
  double low=0,high=hetero_bounds::inf;bool already=true;
  for(size_t i=0;i<c.size();++i) {
    double x=c[i].hi,y=h[i].hi;already=already&&x<=0;
    if(y<0&&x>0)low=std::max(low,hetero_bounds::divpos(I(x),I(-y)).hi);
    if(y>0) {
      if(x>0)return false;
      high=std::min(high,hetero_bounds::divpos(I(-x),I(y)).lo);
    }
    if(y==0&&x>0)return false;
  }
  if(already)return true;
  if(!(low<=high)||!std::isfinite(low)||low>1e100)return false;
  double eta=std::isfinite(high) ? low/2+high/2 : low;
  if(!std::isfinite(eta)||eta>1e100)eta=low;
  for(size_t i=0;i<c.size();++i)if(plus(I(c[i].hi),times(I(eta),I(h[i].hi))).hi>0)return false;
  return true;
}
struct Cell {
  std::array<Poly,5> p;
  std::array<std::array<double,3>,3> vertices{{{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
  I bounds;
  bool excluded=false;
};
inline void evaluate(Cell& cell,const Context& ctx,I original) {
  for(int j=0;j<3;++j) {
    double lo=1,hi=0;for(const auto& v:cell.vertices) { lo=std::min(lo,v[j]);hi=std::max(hi,v[j]); }
    if(hi<ctx.pi[j].lo||lo>ctx.pi[j].hi)cell.excluded=true;
  }
  double hmax=-hetero_bounds::inf;for(I h:cell.p[4])hmax=std::max(hmax,h.hi);
  if(hmax<0)cell.excluded=true;
  if(cell.excluded) { cell.bounds=I(0);return; }
  Poly contrast(cell.p[0].size());double lo=original.lo,hi=original.hi;
  for(int step=0;step<32;++step) {
    double u=lo/2+hi/2;
    for(size_t i=0;i<contrast.size();++i)contrast[i]=minus(times(minus(I(1),I(u)),I(cell.p[0][i].hi)),times(I(u),I(cell.p[3][i].lo)));
    if(certificate(contrast,cell.p[4]))hi=u;else lo=u;
  }
  cell.bounds.hi=hi;lo=original.lo;
  for(int step=0;step<32;++step) {
    double u=lo/2+hi/2;
    for(size_t i=0;i<contrast.size();++i)contrast[i]=minus(times(I(u),I(cell.p[2][i].hi)),times(minus(I(1),I(u)),I(cell.p[1][i].lo)));
    if(certificate(contrast,cell.p[4]))lo=u;else hi=u;
  }
  cell.bounds.lo=lo;
  // Independent valid bounds may prove that this relaxed pi cell has no
  // corresponding Q in the original confidence set.
  if(cell.bounds.lo>cell.bounds.hi) { cell.excluded=true;cell.bounds=I(0); }
}
inline std::array<Cell,2> split(const Cell& parent,const Grid& grid) {
  int a=0,b=1;double longest=-1;
  for(int i=0;i<3;++i)for(int j=i+1;j<3;++j) {
    double d=0;for(int k=0;k<3;++k)d+=std::pow(parent.vertices[i][k]-parent.vertices[j][k],2);
    if(d>longest) { longest=d;a=i;b=j; }
  }
  const int c=3-a-b;
  std::array<Cell,2> out;
  for(Cell& child:out) { child.vertices=parent.vertices;for(Poly& p:child.p)p.resize(grid.alpha.size()); }
  for(int j=0;j<3;++j)out[0].vertices[b][j]=out[1].vertices[a][j]=parent.vertices[a][j]/2+parent.vertices[b][j]/2;
  for(int p=0;p<5;++p)for(int fixed=0;fixed<=grid.n;++fixed) {
    int degree=grid.n-fixed;Poly work(degree+1);
    auto id=[&](int k){std::array<int,3> x{};x[a]=degree-k;x[b]=k;x[c]=fixed;return grid.id(x);};
    for(int k=0;k<=degree;++k)work[k]=parent.p[p][id(k)];
    out[0].p[p][id(0)]=work[0];out[1].p[p][id(degree)]=work[degree];
    for(int level=1;level<=degree;++level) {
      for(int k=0;k<=degree-level;++k)work[k]=times(plus(work[k],work[k+1]),I(.5));
      out[0].p[p][id(level)]=work[0];out[1].p[p][id(degree-level)]=work[degree-level];
    }
  }
  return out;
}
struct Result { I bounds;bool excluded;int cells; };
inline Result probability(const Bank& bank,const Context& ctx,double gamma,I original) {
  if(gamma<=0)return {original,false,0};
  Cell root;root.p=bank.polynomials(ctx,gamma);evaluate(root,ctx,original);
  std::vector<Cell> cells;cells.push_back(std::move(root));
  while(cells.size()<8) {
    Rcpp::checkUserInterrupt();int at=-1;double upper=-1,lower=1;
    for(size_t i=0;i<cells.size();++i)if(!cells[i].excluded) {
      if(cells[i].bounds.hi>upper) { upper=cells[i].bounds.hi;at=i; }
      lower=std::min(lower,cells[i].bounds.lo);
    }
    if(at<0||upper-lower<=1e-4)break;
    auto children=split(cells[at],bank.grid);
    for(Cell& child:children)evaluate(child,ctx,original);
    cells[at]=std::move(children[0]);cells.push_back(std::move(children[1]));
  }
  I bounds(1,0);bool excluded=true;
  for(const Cell& cell:cells)if(!cell.excluded) {
    excluded=false;bounds.lo=std::min(bounds.lo,cell.bounds.lo);bounds.hi=std::max(bounds.hi,cell.bounds.hi);
  }
  return {excluded ? I(0) : hetero_taylor::meet(bounds,original),excluded,int(cells.size())};
}
}
#endif
