#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace {
void v2_fail(const char* message) {
  Rcpp::stop(std::string("fastphylosig V2: ") + message);
}
}

// Stable unsigned UTF-8 byte ordering, independent of character collation.
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector canonical_order_bytes_cpp(const Rcpp::List& keys) {
  std::vector<Rcpp::RawVector> bytes;
  bytes.reserve(keys.size());
  for (R_xlen_t i = 0; i < keys.size(); ++i) bytes.emplace_back(keys[i]);
  Rcpp::IntegerVector order = Rcpp::seq_len(keys.size());
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
    return std::lexicographical_compare(bytes[a-1].begin(), bytes[a-1].end(),
                                        bytes[b-1].begin(), bytes[b-1].end());
  });
  return order;
}

// R validates source fields/UTF-8 and retains all source attributes. This
// representation-only traversal accepts sparse internal IDs and zero lengths;
// method-specific numerical readiness is checked separately.
// [[Rcpp::export(rng = false)]]
Rcpp::List canonical_tree_graph_cpp(const Rcpp::IntegerMatrix& edge,
                                    const Rcpp::NumericVector& edge_length,
                                    const int n_tip, const int n_node,
                                    const Rcpp::IntegerVector& tip_order) {
  const int total = n_tip + n_node, ne = edge.nrow();
  if (n_tip < 2 || n_node < 1 || edge.ncol() != 2 ||
      ne != total-1 || edge_length.size() != ne || tip_order.size() != n_tip)
    v2_fail("invalid canonical graph input dimensions");
  std::unordered_map<int,int> lookup;
  lookup.reserve(total);
  Rcpp::NumericVector node_ids(total);
  for (int i=0; i<n_tip; ++i) { node_ids[i]=i+1; lookup.emplace(i+1,i); }
  int next = n_tip;
  // unique(as.numeric(edge)) visits the complete parent column first.
  for (int col=0; col<2; ++col) for (int r=0; r<ne; ++r) {
    const int id=edge(r,col);
    if (id <= 0) v2_fail("edge endpoints include an unrecognised node ID");
    if (id>n_tip && lookup.find(id)==lookup.end()) {
      if (next == total) v2_fail("edge endpoints do not define exactly Nnode internal IDs");
      lookup.emplace(id,next); node_ids[next++]=id;
    }
  }
  if (next != total) v2_fail("edge endpoints do not define exactly Nnode internal IDs");
  Rcpp::IntegerVector parent_idx(ne), child_idx(ne), degrees(total), indegree(total);
  std::vector<std::vector<int>> kids(total), ordered(total);
  for (int r=0; r<ne; ++r) {
    const int p=lookup.at(edge(r,0)), c=lookup.at(edge(r,1));
    if (p==c) v2_fail("tree contains a self-loop");
    parent_idx[r]=p+1; child_idx[r]=c+1; ++degrees[p]; ++indegree[c];
  }
  for (int i=0; i<total; ++i) kids[i].reserve(degrees[i]);
  for (int r=0; r<ne; ++r) kids[parent_idx[r]-1].push_back(child_idx[r]-1);
  for (int i=0; i<n_tip; ++i) if (degrees[i]!=0) v2_fail("tip nodes cannot be parents");
  for (int i=0; i<n_tip; ++i) if (indegree[i]!=1) v2_fail("each tip must have exactly one parent");
  for (int i=n_tip; i<total; ++i) if (indegree[i]>1) v2_fail("internal nodes must not have multiple parents");
  int root=-1, roots=0;
  for (int i=n_tip; i<total; ++i) if (indegree[i]==0) { root=i; ++roots; }
  if (roots!=1) v2_fail("tree must have exactly one structural root");
  for (int i=n_tip; i<total; ++i) if (degrees[i]<1) v2_fail("internal nodes must have at least one child");
  std::vector<int> colour(total,0), post, stack;
  post.reserve(total); stack.reserve(2*total+4); stack.push_back(root+1);
  int stack_peak=1;
  while (!stack.empty()) {
    const int entry=stack.back(); stack.pop_back();
    const int node=std::abs(entry)-1;
    if (entry<0) { colour[node]=2; post.push_back(node); continue; }
    if (colour[node]==1) v2_fail("tree contains a directed cycle");
    if (colour[node]==2) continue;
    colour[node]=1; stack.push_back(-entry);
    for (auto it=kids[node].rbegin(); it!=kids[node].rend(); ++it) {
      if (colour[*it]==1) v2_fail("tree contains a directed cycle");
      if (colour[*it]!=2) stack.push_back(*it+1);
    }
    stack_peak=std::max(stack_peak,static_cast<int>(stack.size()));
  }
  if (static_cast<int>(post.size())!=total) v2_fail("tree is disconnected from its structural root");
  Rcpp::IntegerVector tip_rank(n_tip), min_rank(total,NA_INTEGER);
  for (int i=0; i<n_tip; ++i) {
    if (tip_order[i]<1 || tip_order[i]>n_tip || tip_rank[tip_order[i]-1]!=0)
      v2_fail("invalid canonical tip order");
    tip_rank[tip_order[i]-1]=i+1;
  }
  for (int i=0; i<n_tip; ++i) min_rank[i]=tip_rank[i];
  for (int node : post) if (node>=n_tip) {
    ordered[node]=kids[node];
    std::sort(ordered[node].begin(),ordered[node].end(),
              [&](int a,int b){return min_rank[a]<min_rank[b];});
    int prev=-1;
    for (int c : ordered[node]) {
      if (min_rank[c]==NA_INTEGER || min_rank[c]==prev)
        v2_fail("child subtrees do not have distinct canonical minimum tip ranks");
      prev=min_rank[c];
    }
    min_rank[node]=min_rank[ordered[node][0]];
  }
  Rcpp::IntegerVector canonical_id(total), preorder(total), canonical_postorder(total), postorder(total);
  for (int i=0; i<n_tip; ++i) canonical_id[i]=i+1;
  stack.clear(); stack.push_back(root+1); int count=0, internal_id=n_tip+1;
  while (!stack.empty()) {
    const int node=stack.back()-1; stack.pop_back(); preorder[count++]=node+1;
    if (node>=n_tip) canonical_id[node]=internal_id++;
    for (auto it=ordered[node].rbegin(); it!=ordered[node].rend(); ++it) stack.push_back(*it+1);
  }
  stack.push_back(root+1); count=0;
  while (!stack.empty()) {
    const int entry=stack.back(); stack.pop_back(); const int node=std::abs(entry)-1;
    if (entry<0) { canonical_postorder[count++]=node+1; continue; }
    stack.push_back(-entry);
    for (auto it=ordered[node].rbegin(); it!=ordered[node].rend(); ++it) stack.push_back(*it+1);
  }
  Rcpp::NumericVector lengths(total,NA_REAL);
  Rcpp::IntegerVector edge_rows(total), parents(total,NA_INTEGER), descriptor_id=Rcpp::clone(canonical_id);
  for (int i=0; i<n_tip; ++i) descriptor_id[i]=tip_rank[i];
  for (int r=0; r<ne; ++r) { int c=child_idx[r]-1; lengths[c]=edge_length[r]; edge_rows[c]=r+1; parents[c]=parent_idx[r]; }
  Rcpp::List children(total), ordered_children(total), descriptor(total);
  for (int i=0; i<total; ++i) {
    Rcpp::IntegerVector ch(kids[i].size());
    for (size_t j=0; j<kids[i].size(); ++j) ch[j]=kids[i][j]+1;
    children[i]=ch; postorder[i]=post[i]+1;
    if (i<n_tip) {
      descriptor[descriptor_id[i]-1]=Rcpp::List::create(
        Rcpp::_["kind"]="tip", Rcpp::_["label_rank"]=tip_rank[i],
        Rcpp::_["child_descriptor_ids"]=Rcpp::IntegerVector(0));
    } else {
      Rcpp::IntegerVector oc(ordered[i].size()), ranks(oc.size()), ids(oc.size());
      for (R_xlen_t j=0; j<oc.size(); ++j) {
        int c=ordered[i][j]; oc[j]=c+1; ranks[j]=min_rank[c]; ids[j]=descriptor_id[c];
      }
      ordered_children[i]=oc;
      descriptor[descriptor_id[i]-1]=Rcpp::List::create(
        Rcpp::_["kind"]="internal", Rcpp::_["child_count"]=static_cast<int>(oc.size()),
        Rcpp::_["child_min_tip_ranks"]=ranks, Rcpp::_["child_descriptor_ids"]=ids);
    }
  }
  return Rcpp::List::create(
    Rcpp::_["node_ids"]=node_ids, Rcpp::_["parent_idx"]=parent_idx,
    Rcpp::_["child_idx"]=child_idx, Rcpp::_["children"]=children,
    Rcpp::_["ordered_children"]=ordered_children, Rcpp::_["root_idx"]=root+1,
    Rcpp::_["tip_idx"]=Rcpp::seq_len(n_tip), Rcpp::_["internal_idx"]=Rcpp::seq(n_tip+1,total),
    Rcpp::_["tip_order"]=tip_order, Rcpp::_["tip_rank"]=tip_rank,
    Rcpp::_["min_tip_rank"]=min_rank, Rcpp::_["canonical_id"]=canonical_id,
    Rcpp::_["parent_for_child"]=parents, Rcpp::_["descriptor_id"]=descriptor_id,
    Rcpp::_["descriptor"]=descriptor, Rcpp::_["postorder"]=postorder,
    Rcpp::_["canonical_preorder"]=preorder, Rcpp::_["canonical_postorder"]=canonical_postorder,
    Rcpp::_["edge_len_by_child"]=lengths, Rcpp::_["edge_row_for_child"]=edge_rows,
    Rcpp::_["metrics"]=Rcpp::List::create(Rcpp::_["schema_version"]=2,
      Rcpp::_["full_edge_scans"]=1, Rcpp::_["descendant_label_vectors_materialized"]=0,
      Rcpp::_["child_tuple_entries"]=ne, Rcpp::_["descriptor_entries"]=total,
      Rcpp::_["child_sort_calls"]=n_node, Rcpp::_["explicit_stack_peak"]=stack_peak));
}
