#include <Rcpp.h>
#include <cstdint>
#include <unordered_map>
#include <memory>
#include "heterogeneity_bounds.h"
#include "heterogeneity_ard_boundary.h"
#include "heterogeneity_taylor.h"
#include "heterogeneity_forest_shared.h"
#include "heterogeneity_transient.h"
#include "heterogeneity_mass.h"

namespace {
using hetero_bounds::I;
struct NuisanceWorkspace {
  hetero_bounds::Tree tree;
  Rcpp::IntegerMatrix patterns,index,observed;
  std::vector<I> proposal;
  std::vector<int> tail;
  std::vector<int> pattern_id;
  std::unique_ptr<hetero_forest::Summary> forest_summary;
  std::unique_ptr<hetero_forest::Bank> forest;
  std::unique_ptr<hetero_taylor::RepeatPlan> repeats;
  std::unique_ptr<hetero_mass::Bank> mass;
  int n_tail=0,observed_states=0,observed_duplicate_other=0;
  std::vector<int> support_isolated;
  bool er;
  NuisanceWorkspace(const Rcpp::List& core,const Rcpp::IntegerMatrix& x,
      const Rcpp::IntegerMatrix& rates,const Rcpp::NumericMatrix& g,
      const Rcpp::LogicalVector& upper_tail,bool equal):
      tree(core,rates.nrow()),index(Rcpp::clone(rates)),er(equal) {
    tree.validate(x);
    if(index.ncol()!=tree.k || g.nrow()!=x.ncol() || g.ncol()!=2 ||
       upper_tail.size()!=x.ncol() || !x.ncol() || upper_tail[0]!=TRUE)
      Rcpp::stop("invalid nuisance bank dimensions");
    if(!er && tree.k!=3)Rcpp::stop("nuisance ARD currently requires three states");
    // Repeated phenotypes have identical likelihoods even when independent
    // Delta MCMC gives different H. Compress patterns, never their tail flags.
    std::unordered_map<uint64_t,std::vector<int>> buckets;
    std::vector<int> representatives;
    for(int b=0;b<x.ncol();++b) {
      uint64_t hash=1469598103934665603ULL;
      for(int i=0;i<x.nrow();++i) { hash^=static_cast<uint32_t>(x(i,b));hash*=1099511628211ULL; }
      auto& bucket=buckets[hash];int found=-1;
      for(int candidate:bucket) {
        int source=representatives[candidate];bool same=true;
        for(int i=0;i<x.nrow();++i)if(x(i,b)!=x(i,source)) { same=false;break; }
        if(same) { found=candidate;break; }
      }
      if(found<0) { found=representatives.size();representatives.push_back(b);bucket.push_back(found); }
      pattern_id.push_back(found);
    }
    patterns=Rcpp::IntegerMatrix(x.nrow(),representatives.size());
    observed=Rcpp::IntegerMatrix(x.nrow(),1);
    std::copy(x.begin(),x.begin()+x.nrow(),observed.begin());
    for(size_t b=0;b<representatives.size();++b)
      std::copy(x.begin()+representatives[b]*x.nrow(),x.begin()+(representatives[b]+1)*x.nrow(),
        patterns.begin()+b*x.nrow());
    for(int i=0;i<tree.k;++i)for(int j=0;j<tree.k;++j)
      if(index(i,j)!=(i==j ? 0 : index(i,j)) || (i!=j &&
        (index(i,j)<1 || index(i,j)>(er ? 1 : 6))))Rcpp::stop("invalid nuisance rate index");
    for(int i=0;i<x.ncol();++i) {
      if(!std::isfinite(g(i,0))||!std::isfinite(g(i,1))||g(i,0)>g(i,1)||
         Rcpp::LogicalVector::is_na(upper_tail[i]))Rcpp::stop("invalid nuisance proposal bounds/tail");
      proposal.emplace_back(g(i,0),g(i,1));tail.push_back(upper_tail[i]==TRUE);
      n_tail+=tail.back();
      observed_duplicate_other+=!tail.back()&&pattern_id[i]==0;
    }
    // Support obstructions are diagnostic candidates, never upper-bound
    // certificates. Check their likelihood and C_gamma membership normally.
    std::vector<int> masks(patterns.ncol(),0);
    for(int b=0;b<patterns.ncol();++b)for(int i=0;i<patterns.nrow();++i)
      if(patterns(i,b)!=NA_INTEGER)masks[b]|=1<<(patterns(i,b)-1);
    for(int j=0;j<tree.k;++j)observed_states+=(masks[0]>>j)&1;
    if(!er)for(int isolated=0;isolated<3;++isolated) {
      const int bit=1<<isolated;
      auto compatible=[&](int mask){return !(mask&bit)||mask==bit;};
      if(!compatible(masks[0]))continue;
      bool possible=true;
      for(size_t i=0;i<tail.size();++i)
        if(!tail[i]&&compatible(masks[pattern_id[i]])) {possible=false;break;}
      if(possible)support_isolated.push_back(isolated+1);
    }
  }
};
NuisanceWorkspace& get_workspace(SEXP pointer) {
  if(TYPEOF(pointer)!=EXTPTRSXP || !Rf_inherits(pointer,"fastphylosig_nuisance_workspace"))
    Rcpp::stop("invalid nuisance workspace");
  Rcpp::XPtr<NuisanceWorkspace> p(pointer);
  if(!p.get())Rcpp::stop("released nuisance workspace");
  return *p;
}
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector heterogeneity_nuisance_add_gamma_cpp(double lower,double upper,double gamma) {
  if(!std::isfinite(lower)||!std::isfinite(upper)||lower<0||upper>1||lower>upper||
      !std::isfinite(gamma)||gamma<0||gamma>=1)Rcpp::stop("invalid nuisance bound addition");
  I value=hetero_bounds::add(I(lower,upper),I(gamma));
  return Rcpp::NumericVector::create(std::min(1.0,std::max(0.0,value.lo)),std::min(1.0,value.hi));
}

// [[Rcpp::export(rng = false)]]
SEXP heterogeneity_nuisance_prepare_cpp(const Rcpp::List& core,
    const Rcpp::IntegerMatrix& patterns,const Rcpp::IntegerMatrix& rate_index,
    const Rcpp::NumericMatrix& log_proposal,const Rcpp::LogicalVector& tail,
    const bool equal_rates) {
  Rcpp::XPtr<NuisanceWorkspace> p(new NuisanceWorkspace(core,patterns,rate_index,
    log_proposal,tail,equal_rates),true);
  p.attr("class")="fastphylosig_nuisance_workspace";
  p.attr("bank_diagnostics")=Rcpp::List::create(
    Rcpp::_["upper_tail_count"]=p->n_tail,
    Rcpp::_["non_tail_count"]=static_cast<int>(p->tail.size())-p->n_tail,
    Rcpp::_["observed_state_count"]=p->observed_states,
    Rcpp::_["observed_duplicate_non_tail_count"]=p->observed_duplicate_other,
    Rcpp::_["support_isolated_states"]=Rcpp::wrap(p->support_isolated));
  return p;
}

// [[Rcpp::export(rng = false)]]
void heterogeneity_nuisance_release_cpp(SEXP pointer) {
  get_workspace(pointer);Rcpp::XPtr<NuisanceWorkspace> p(pointer);p.release();
}

// Density mixtures retain outward bounds, including the log-space operations.
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix heterogeneity_nuisance_mixture_cpp(const Rcpp::List& components) {
  if(!components.size())Rcpp::stop("empty nuisance proposal mixture");
  Rcpp::NumericMatrix first=components[0];const int n=first.nrow();
  if(first.ncol()!=2)Rcpp::stop("invalid mixture component bounds");
  std::vector<I> answer(n);
  for(int i=0;i<n;++i) {
    double pivot=-hetero_bounds::inf;
    for(SEXP element:components) {
      Rcpp::NumericMatrix c(element);
      if(c.nrow()!=n||c.ncol()!=2)Rcpp::stop("inconsistent mixture bounds");
      if(std::isnan(c(i,0))||std::isnan(c(i,1))||c(i,0)>c(i,1))Rcpp::stop("invalid mixture bounds");
      pivot=std::max(pivot,c(i,1));
    }
    if(pivot==-hetero_bounds::inf) { answer[i]=I(pivot);continue; }
    I sum(0);
    for(SEXP element:components) {
      Rcpp::NumericMatrix c(element);
      I log=hetero_bounds::sub(I(c(i,0),c(i,1)),I(pivot));
      sum=hetero_bounds::posadd(sum,hetero_bounds::exp_nonpositive(log));
    }
    answer[i]=hetero_bounds::add(I(pivot),hetero_bounds::log_interval(
      hetero_bounds::divpos(sum,I(components.size()))));
  }
  return hetero_bounds::matrix(answer);
}

// No optimizer decides validity: every returned upper bound encloses an entire
// continuous rate box, including the complete unbounded ARD boundary.
// [[Rcpp::export(rng = false)]]
Rcpp::List heterogeneity_nuisance_box_cpp(SEXP pointer,const Rcpp::NumericVector& lower,
    const Rcpp::NumericVector& upper,const double gamma,const bool keep_likelihood=false,
    const int boundary_power=24,const bool constrain_rows=true,const int correlation=0) {
  NuisanceWorkspace& w=get_workspace(pointer);const int d=w.er ? 1 : 6;
  if(lower.size()!=d||upper.size()!=d||!std::isfinite(gamma)||gamma<0||gamma>=1)
    Rcpp::stop("invalid nuisance parameter box/gamma");
  if(boundary_power<8||boundary_power>30)Rcpp::stop("boundary_power must be in 8..30");
  if(correlation<0||correlation>4)Rcpp::stop("correlation must be in 0..4");
  for(int j=0;j<d;++j)if((!std::isfinite(lower[j])&&!(w.er&&lower[j]==R_PosInf))||lower[j]<0||
      std::isnan(upper[j])||upper[j]<lower[j])Rcpp::stop("invalid nuisance rate bounds");
  std::vector<std::vector<I>> transition;
  transition.reserve(w.tree.lengths.size());
  std::vector<I> rates(w.tree.k*w.tree.k);
  if(!w.er)for(int i=0;i<w.tree.k;++i)for(int j=0;j<w.tree.k;++j)if(i!=j) {
    int a=w.index(i,j)-1;rates[i*w.tree.k+j]=I(lower[a],upper[a]);
  }
  for(double t:w.tree.lengths)transition.push_back(w.er ?
    hetero_bounds::er_transition(t,w.tree.k,lower[0],upper[0]) :
    hetero_bounds::ard_box_transition(t,rates,boundary_power));
  // Most high-rate tails are discarded by the observed likelihood alone.
  // Do not prune a complete reference bank merely to discover this fact.
  if(gamma>0&&!keep_likelihood) {
    const I observed=w.tree.likelihood(w.observed,transition,constrain_rows)[0];
    const I threshold=hetero_bounds::add(hetero_bounds::logpoint(gamma),w.proposal[0]);
    if(observed.hi<threshold.lo)return Rcpp::List::create(Rcpp::_["lower"]=0.0,
      Rcpp::_["upper"]=0.0,Rcpp::_["excluded"]=true,Rcpp::_["included"]=false,
      Rcpp::_["observed_weight_upper"]=1.0,Rcpp::_["ESS_lower"]=0.0,
      Rcpp::_["n_unique_patterns"]=w.patterns.ncol());
  }
  auto unique=w.tree.likelihood(w.patterns,transition,constrain_rows);
  std::vector<hetero_taylor::Likelihood> shared;
  bool correlation_used=false,correlation_fallback=false;
  int likelihood_tightened=0;
  if((correlation==1||correlation==2)&&!w.er) {
    bool finite=true,varying=false;double radius=0,largest_exit=0;
    for(int i=0;i<3;++i) {
      double width=0,exit=0;
      for(int j=0;j<3;++j)if(i!=j) {
        const I q=rates[i*3+j];finite=finite&&std::isfinite(q.hi);
        varying=varying||q.hi>q.lo;width+=q.hi/2-q.lo/2;exit+=q.hi;
      }
      radius=std::max(radius,width);largest_exit=std::max(largest_exit,exit);
    }
    const double t=w.tree.lengths.back();
    const bool economical=radius*t<=.15&&radius*w.tree.total_length<=.5&&largest_exit*t<=256&&
      double(w.tree.nt)*w.patterns.ncol()<=200000;
    if(finite&&varying&&(correlation==2||economical)) {
      try {
        auto p=hetero_taylor::transitions(w.tree,rates,w.index,transition);
        shared=hetero_taylor::likelihood(w.tree,w.patterns,p);
        auto tightened=unique;
        for(size_t i=0;i<unique.size();++i) {
          tightened[i]=hetero_taylor::meet(unique[i],shared[i].log);
          if(tightened[i].lo>unique[i].lo||tightened[i].hi<unique[i].hi)++likelihood_tightened;
          I shift=hetero_taylor::times(I(static_cast<double>(shared[i].exponent)),hetero_bounds::log2_interval());
          hetero_taylor::constrain(shared[i].value,hetero_taylor::positive_exp(hetero_bounds::sub(tightened[i],shift)));
        }
        unique.swap(tightened);correlation_used=true;
      } catch(const std::overflow_error&) {
        shared.clear();correlation_fallback=true;likelihood_tightened=0;
      }
    }
  }
  std::unique_ptr<hetero_forest::Context> forest_context;
  bool forest_used=false,forest_fallback=false;int forest_cells=0;
  if(!w.er&&(correlation==3||(correlation==1&&!correlation_used))&&
      hetero_forest::eligible(w.patterns)) {
    bool varying=false;for(int j=0;j<d;++j)varying=varying||upper[j]>lower[j];
    if(varying)try {
      forest_context.reset(new hetero_forest::Context(w.tree,rates,w.patterns));
      if(forest_context->beta>0&&forest_context->eps<=(correlation==3 ? .25 : .05)) {
        if(!w.forest_summary)w.forest_summary.reset(new hetero_forest::Summary(w.patterns));
        if(w.forest_summary->n<=128&&double(w.tree.nt)*w.patterns.ncol()<=200000&&!w.forest)
          w.forest.reset(new hetero_forest::Bank(*w.forest_summary,w.pattern_id,w.proposal,w.tail));
        auto extra=w.forest_summary->likelihood(*forest_context),tightened=unique;
        for(size_t i=0;i<unique.size();++i)tightened[i]=hetero_taylor::meet(unique[i],extra[i]);
        unique.swap(tightened);forest_used=true;
      } else forest_context.reset();
    } catch(const std::overflow_error&) { forest_context.reset();forest_fallback=true; }
  }
  std::vector<I> ll(w.pattern_id.size());
  for(size_t i=0;i<ll.size();++i)ll[i]=unique[w.pattern_id[i]];
  bool excluded=false,included=gamma==0;
  if(gamma>0) {
    I threshold=hetero_bounds::add(hetero_bounds::logpoint(gamma),w.proposal[0]);
    excluded=ll[0].hi<threshold.lo;included=ll[0].lo>=threshold.hi;
  }
  I probability(0,1);double observed_upper=1,ess_lower=0;
  if(excluded)probability=I(0);
  else {
    std::vector<I> log_weight(ll.size());double pivot=-hetero_bounds::inf;
    for(size_t i=0;i<ll.size();++i) {
      log_weight[i]=hetero_bounds::sub(ll[i],w.proposal[i]);
      pivot=std::max(pivot,log_weight[i].hi);
    }
    if(pivot!=-hetero_bounds::inf) {
      I numerator(0),other(0),squares(0),observed(0);
      for(size_t i=0;i<ll.size();++i) {
        I weight=hetero_bounds::exp_nonpositive(hetero_bounds::sub(log_weight[i],I(pivot)));
        if(w.tail[i])numerator=hetero_bounds::posadd(numerator,weight);
        else other=hetero_bounds::posadd(other,weight);
        squares=hetero_bounds::posadd(squares,hetero_bounds::posmul(weight,weight));
        if(i==0)observed=weight;
      }
      const double lo_sum=hetero_bounds::up(numerator.lo+other.hi);
      const double hi_sum=std::max(0.0,hetero_bounds::dn(numerator.hi+other.lo));
      probability=I(lo_sum>0 ? std::max(0.0,hetero_bounds::dn(numerator.lo/lo_sum)) : 0,
        hi_sum>0 ? std::min(1.0,hetero_bounds::up(numerator.hi/hi_sum)) : 1);
      I total=hetero_bounds::posadd(numerator,other);
      observed_upper=total.lo>0 ? std::min(1.0,hetero_bounds::up(observed.hi/total.lo)) : 1;
      if(squares.hi>0)ess_lower=std::max(0.0,hetero_bounds::dn(total.lo*total.lo/squares.hi));
    }
  }
  const double before_separator=probability.hi;
  if(correlation_used&&!excluded) {
    try {
      probability=hetero_taylor::rank(shared,w.pattern_id,w.proposal,w.tail,gamma,probability);
    } catch(const std::overflow_error&) { correlation_fallback=true; }
  }
  const double before_forest=probability.hi;
  if(forest_used&&!excluded&&w.forest)try {
    auto fit=hetero_forest::probability(*w.forest,*forest_context,gamma,probability);
    probability=fit.bounds;forest_cells=fit.cells;
    if(fit.excluded) { excluded=true;included=false; }
  } catch(const std::overflow_error&) { forest_fallback=true; }
  bool transient_used=false,transient_fallback=false;int transient_cells=0;
  double transient_messages=0,transient_reuses=0;
  double transient_delta=hetero_bounds::inf,transient_error=1;
  const double before_transient=probability.hi;
  if(!w.er&&(correlation==4||(correlation==1&&!correlation_used&&!forest_used))&&
      double(w.tree.nt)*w.patterns.ncol()<=8192) {
    bool varying=false;for(int j=0;j<d;++j)varying=varying||upper[j]>lower[j];
    if(varying)try {
      hetero_transient::Context ctx(w.tree,rates);
      transient_delta=ctx.delta;transient_error=ctx.error;
      if(ctx.family>=0&&ctx.error<=(correlation==4 ? .5 : .15)) {
        if(!w.repeats)w.repeats.reset(new hetero_taylor::RepeatPlan(w.tree,w.patterns));
        auto fit=hetero_transient::evaluate(w.tree,w.patterns,ctx,w.pattern_id,
          w.proposal,w.tail,gamma,probability,unique,w.repeats.get());
        transient_used=true;transient_cells=fit.cells;
        transient_messages=double(fit.cells)*w.repeats->messages;
        transient_reuses=double(fit.cells)*(w.repeats->serial_messages-w.repeats->messages);
        unique=std::move(fit.likelihood);probability=fit.probability;
        for(size_t i=0;i<ll.size();++i)ll[i]=unique[w.pattern_id[i]];
        if(gamma>0) {
          I threshold=hetero_bounds::add(hetero_bounds::logpoint(gamma),w.proposal[0]);
          excluded=fit.excluded||ll[0].hi<threshold.lo;
          included=!excluded&&ll[0].lo>=threshold.hi;
        }
      }
    } catch(const std::overflow_error&) { transient_fallback=true; }
  }
  const double before_mass=probability.hi;
  bool mass_used=false,mass_fallback=false;
  if(!w.er&&(correlation==1||correlation==4)&&!excluded)try {
    if(!w.mass)w.mass.reset(new hetero_mass::Bank(w.tree,w.patterns,w.pattern_id,w.proposal,w.tail));
    mass_used=w.mass->useful(gamma);
    if(mass_used) {
      probability=w.mass->upper(gamma,probability);
      if(probability.lo>probability.hi) {probability=I(0);excluded=true;included=false;}
    }
  } catch(const std::overflow_error&) {mass_fallback=true;}
  Rcpp::List out=Rcpp::List::create(Rcpp::_["lower"]=probability.lo,
    Rcpp::_["upper"]=probability.hi,Rcpp::_["excluded"]=excluded,
    Rcpp::_["included"]=included,Rcpp::_["observed_weight_upper"]=observed_upper,
    Rcpp::_["ESS_lower"]=ess_lower,Rcpp::_["n_unique_patterns"]=w.patterns.ncol());
  if(keep_likelihood)out["log_likelihood"]=hetero_bounds::matrix(ll);
  if(correlation&&!w.er) {
    out["correlation_used"]=correlation_used;out["correlation_fallback"]=correlation_fallback;
    out["correlation_likelihood_tightened"]=likelihood_tightened;
    out["P_upper_before_separator"]=before_separator;
  }
  if((correlation==1||correlation==3)&&!w.er) {
    out["forest_used"]=forest_used;out["forest_fallback"]=forest_fallback;
    out["forest_cells"]=forest_cells;out["P_upper_before_forest"]=before_forest;
    out["forest_beta"]=forest_context ? forest_context->beta : 0;
    out["forest_tip_error"]=forest_context ? forest_context->eps : 1;
  }
  if((correlation==1||correlation==4)&&!w.er) {
    out["transient_used"]=transient_used;out["transient_fallback"]=transient_fallback;
    out["transient_cells"]=transient_cells;out["transient_delta"]=transient_delta;
    out["transient_error"]=transient_error;out["P_upper_before_transient"]=before_transient;
    out["transient_messages"]=transient_messages;out["transient_reuses"]=transient_reuses;
    out["mass_used"]=mass_used;out["mass_fallback"]=mass_fallback;
    out["mass_constant_floor"]=w.mass ? w.mass->floor : 0;
    out["P_upper_before_mass"]=before_mass;
  }
  return out;
}

// Independent transition-envelope audit; no simulation and no RNG scope.
// [[Rcpp::export(rng = false)]]
Rcpp::List heterogeneity_ard_boundary_cpp(const Rcpp::NumericMatrix& lower,
    const Rcpp::NumericMatrix& upper,double time=1,const int boundary_power=24) {
  if(lower.nrow()!=3||lower.ncol()!=3||upper.nrow()!=3||upper.ncol()!=3||
      !(time>0)||!std::isfinite(time)||boundary_power<8||boundary_power>30)
    Rcpp::stop("invalid three-state transition audit");
  std::vector<I> rates(9);
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)if(i!=j) {
    if(!std::isfinite(lower(i,j))||lower(i,j)<0||std::isnan(upper(i,j))||upper(i,j)<lower(i,j))
      Rcpp::stop("invalid transition rate box");
    rates[i*3+j]=I(lower(i,j),upper(i,j));
  }
  const auto resolvent=hetero_bounds::resolvent3(time,rates,boundary_power);
  const auto transition=hetero_bounds::ard_box_transition(time,rates,boundary_power);
  Rcpp::NumericMatrix rlo(3,3),rhi(3,3),lo(3,3),hi(3,3);
  for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
    rlo(i,j)=resolvent[i*3+j].lo;rhi(i,j)=resolvent[i*3+j].hi;
    lo(i,j)=transition[i*3+j].lo;hi(i,j)=transition[i*3+j].hi;
  }
  return Rcpp::List::create(Rcpp::_["resolvent_lower"]=rlo,Rcpp::_["resolvent_upper"]=rhi,
    Rcpp::_["lower"]=lo,Rcpp::_["upper"]=hi,Rcpp::_["Euler_error"]=std::ldexp(32.0,-boundary_power));
}
