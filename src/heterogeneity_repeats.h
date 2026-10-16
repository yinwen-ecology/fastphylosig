#ifndef FASTPHYLOSIG_HETEROGENEITY_REPEATS_H
#define FASTPHYLOSIG_HETEROGENEITY_REPEATS_H
#include "heterogeneity_taylor.h"
#include <unordered_map>
#include <cstdint>

namespace hetero_taylor {
// Parameter-independent repeats at the SAME tree node. Numeric messages are
// recomputed for every rectangle; only their index plan and scratch are reused.
struct RepeatPlan {
  struct Task {int node,left,right,state;};
  std::vector<Task> tasks;
  std::vector<int> root_id;
  mutable std::vector<T> scratch;
  mutable std::vector<long long> exponents;
  long long messages=0,serial_messages=0;
  RepeatPlan(const hetero_bounds::Tree& tree,const Rcpp::IntegerMatrix& x) {
    std::vector<std::vector<int>> id(tree.nn,std::vector<int>(x.ncol()));
    for(int pos=tree.nn-1;pos>=0;--pos) {
      int node=tree.pre[pos]-1;
      std::unordered_map<uint64_t,int> seen;
      for(int b=0;b<x.ncol();++b) {
        int l=-1,r=-1,state=0;uint64_t key;
        if(node<tree.nt) {state=x(node,b)==NA_INTEGER ? 0 : x(node,b);key=state;}
        else {
          l=id[tree.children[node][0]][b];r=id[tree.children[node][1]][b];
          key=(uint64_t(static_cast<uint32_t>(l))<<32)|static_cast<uint32_t>(r);
        }
        auto found=seen.find(key);
        if(found==seen.end()) {
          int at=tasks.size();tasks.push_back({node,l,r,state});seen.emplace(key,at);id[node][b]=at;
          if(node>=tree.nt)++messages;
        } else id[node][b]=found->second;
      }
    }
    root_id=std::move(id[tree.root]);serial_messages=static_cast<long long>(tree.nn-tree.nt)*x.ncol();
    scratch.resize(tasks.size()*3);exponents.resize(tasks.size());
  }
  std::vector<Likelihood> likelihood(const hetero_bounds::Tree& tree,const std::vector<Matrix>& p) const {
    for(size_t at=0;at<tasks.size();++at) {
      if((at&255)==0)Rcpp::checkUserInterrupt();
      const Task& task=tasks[at];T* target=&scratch[at*3];
      if(task.node<tree.nt) {
        for(int i=0;i<3;++i)target[i]=T(I(task.state==0||task.state==i+1 ? 1 : 0));
        exponents[at]=0;continue;
      }
      const T* child[]={&scratch[task.left*3],&scratch[task.right*3]};
      const long long shifts[]={exponents[task.left],exponents[task.right]};
      const auto& children=tree.children[task.node];
      const Matrix* matrices[]={&p[tree.length_id[children[0]]],&p[tree.length_id[children[1]]]};
      exponents[at]=message(target,child,shifts,matrices);
    }
    std::vector<Likelihood> out;out.reserve(root_id.size());
    for(int id:root_id)out.push_back(root_likelihood(&scratch[id*3],exponents[id]));
    return out;
  }
};
}
#endif
