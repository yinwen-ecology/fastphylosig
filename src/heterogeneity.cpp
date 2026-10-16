#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>
#include "tree_core.h"
#include "continuous_kernels.h"
#include "numeric_utils.h"
#include "heterogeneity_bootstrap.h"

namespace {
struct CladeIndex {
  Rcpp::IntegerVector count, start, end;
};
CladeIndex clade_index(const Rcpp::List& core) {
  const int n_tip = Rcpp::as<int>(core["n_tip"]), n = Rcpp::as<int>(core["n_total"]);
  Rcpp::IntegerVector parent = core["parent"], pre = core["preorder"], post = core["postorder"];
  if (n_tip < 2 || n <= n_tip || parent.size() != n || pre.size() != n || post.size() != n)
    Rcpp::stop("inconsistent compiled clade tree");
  CladeIndex x{Rcpp::IntegerVector(n), Rcpp::IntegerVector(n), Rcpp::IntegerVector(n)};
  std::vector<int> size(n, 1);
  for (int i = 0; i < n; ++i) {
    const int node = pre[i] - 1;
    if (node < 0 || node >= n || x.start[node]) Rcpp::stop("invalid clade preorder");
    x.start[node] = i + 1;
  }
  std::vector<char> seen(n, 0);
  for (int i = 0; i < n_tip; ++i) x.count[i] = 1;
  for (int i = 0; i < n; ++i) {
    const int node = post[i] - 1;
    if (node < 0 || node >= n || seen[node]) Rcpp::stop("invalid clade postorder");
    seen[node] = 1;
    const int p = parent[node] - 1;
    if (parent[node] < 0 || parent[node] > n ||
        (p < 0 && node != pre[0] - 1)) Rcpp::stop("invalid clade root/parent");
    if (p >= 0) {
      if (p >= n || x.start[p] >= x.start[node] || seen[p]) Rcpp::stop("invalid clade parent");
      x.count[p] += x.count[node];
      size[p] += size[node];
    }
  }
  for (int node = 0; node < n; ++node) x.end[node] = x.start[node] + size[node] - 1;
  return x;
}

double lambda_bound(const Rcpp::List& core) {
  const int n_tip = Rcpp::as<int>(core["n_tip"]);
  const Rcpp::NumericVector depth = core["root_distance"];
  const Rcpp::IntegerVector parent = core["parent"];
  double h = 0.0, low = std::numeric_limits<double>::infinity();
  for (int tip = 0; tip < n_tip; ++tip) { h = std::max(h, depth[tip]); low = std::min(low, depth[tip]); }
  if (!std::isfinite(h) || h <= 0.0 || (h - low) / h > std::sqrt(std::numeric_limits<double>::epsilon())) return 1.0;
  double bound = std::numeric_limits<double>::infinity();
  for (int tip = 0; tip < n_tip; ++tip) {
    const double pd = depth[parent[tip] - 1];
    if (pd > 0.0) bound = std::min(bound, depth[tip] / pd);
  }
  return std::isfinite(bound) && bound > 0.0 ? bound : 1.0;
}
double finite_mean(const std::vector<double>& x) {
  long double sum = 0.0L;
  for (double v : x) sum += v;
  long double mean = sum / x.size(), correction = 0.0L;
  for (double v : x) correction += static_cast<long double>(v) - mean;
  return static_cast<double>(mean + correction / x.size());
}
} // namespace

// Linear tree annotation, with source-node identity and exact crown membership
// checks. Canonical plot IDs accommodate ape's root = n_tip + 1 convention.
// [[Rcpp::export(rng = false)]]
Rcpp::List heterogeneity_tree_map_cpp(const Rcpp::List& core,
    const Rcpp::IntegerVector& nodes, const Rcpp::IntegerVector& tip_ids,
    const Rcpp::IntegerVector& tip_offsets) {
  const int n_tip = Rcpp::as<int>(core["n_tip"]), n = Rcpp::as<int>(core["n_total"]);
  const Rcpp::IntegerVector pre = core["preorder"], parent = core["parent"];
  const Rcpp::NumericVector branch = core["branch_length_by_node"];
  const CladeIndex ix = clade_index(core);
  const int nc = nodes.size();
  if (nc < 1 || tip_offsets.size() != nc + 1 || tip_offsets[0] != 0 ||
      tip_offsets[nc] != tip_ids.size() || branch.size() != n)
    Rcpp::stop("inconsistent tree annotation dimensions");
  std::vector<int> selected(n, 0), seen(n_tip, -1), owner(n, 0);
  for (int i = 0; i < nc; ++i) {
    Rcpp::checkUserInterrupt();
    const int id = nodes[i];
    if (id <= n_tip || id > n || selected[id - 1])
      Rcpp::stop("invalid or repeated annotation node");
    const int node = id - 1;
    selected[node] = i + 1;
    if (tip_offsets[i] < 0 || tip_offsets[i + 1] < tip_offsets[i] ||
        tip_offsets[i + 1] > tip_ids.size() ||
        tip_offsets[i + 1] - tip_offsets[i] != ix.count[node])
      Rcpp::stop("tree does not match the result's selected crowns");
    for (int k = tip_offsets[i]; k < tip_offsets[i + 1]; ++k) {
      const int id_tip = tip_ids[k];
      if (id_tip < 1 || id_tip > n_tip)
        Rcpp::stop("tree does not match the result's selected crowns");
      const int tip = id_tip - 1;
      if (seen[tip] == i ||
          ix.start[tip] < ix.start[node] || ix.start[tip] > ix.end[node])
        Rcpp::stop("tree does not match the result's selected crowns");
      seen[tip] = i;
    }
  }
  Rcpp::IntegerVector source_to_plot(n), source_internal(n - n_tip), node_owner(n);
  for (int tip = 0; tip < n_tip; ++tip) source_to_plot[tip] = tip + 1;
  int internal = 0;
  for (int node : pre) {
    if (node > n_tip) {
      source_internal[internal] = node;
      source_to_plot[node - 1] = n_tip + (++internal);
    }
  }
  Rcpp::IntegerMatrix edge(n - 1, 2);
  Rcpp::IntegerVector edge_owner(n - 1);
  Rcpp::NumericVector length(n - 1);
  for (int pos = 0; pos < n; ++pos) {
    if ((pos & 8191) == 0) Rcpp::checkUserInterrupt();
    const int node = pre[pos] - 1, p = parent[node] - 1;
    owner[node] = selected[node] ? selected[node] : (p >= 0 ? owner[p] : 0);
    node_owner[source_to_plot[node] - 1] = owner[node];
    if (pos) {
      edge(pos - 1, 0) = source_to_plot[p];
      edge(pos - 1, 1) = source_to_plot[node];
      length[pos - 1] = branch[node];
      // The incoming stem is outside a selected crown.
      edge_owner[pos - 1] = selected[node] ? owner[p] : owner[node];
    }
  }
  return Rcpp::List::create(Rcpp::_["edge"] = edge, Rcpp::_["edge_length"] = length,
    Rcpp::_["edge_owner"] = edge_owner, Rcpp::_["node_owner"] = node_owner,
    Rcpp::_["source_to_plot"] = source_to_plot, Rcpp::_["source_internal"] = source_internal);
}

// Iterative selection/extraction keeps source IDs and avoids a whole-tree
// lookup allocation per crown. The mapping workspace is shared by all crowns.
// [[Rcpp::export(rng = false)]]
Rcpp::List heterogeneity_plan_cpp(const Rcpp::List& core,
    const Rcpp::CharacterVector& labels, SEXP requested_nodes,
    const int min_tips, const int max_tips, const double max_total_tips) {
  const int n_tip = Rcpp::as<int>(core["n_tip"]), n = Rcpp::as<int>(core["n_total"]);
  const Rcpp::IntegerVector pre = core["preorder"], parent = core["parent"];
  const Rcpp::NumericVector branch = core["branch_length_by_node"];
  if (labels.size() != n_tip || branch.size() != n || min_tips < 3 || max_tips < 3 ||
      (requested_nodes == R_NilValue && max_tips < min_tips))
    Rcpp::stop("invalid clade plan dimensions");
  const CladeIndex ix = clade_index(core);
  std::vector<int> nodes;
  if (requested_nodes != R_NilValue) nodes = Rcpp::as<std::vector<int> >(requested_nodes);
  else {
    for (int pos = 0; pos < n;) {
      const int node = pre[pos];
      if (node > n_tip &&
          ix.count[node - 1] >= min_tips && ix.count[node - 1] <= max_tips) {
        nodes.push_back(node); pos = ix.end[node - 1];
      } else if (ix.count[node - 1] < min_tips) pos = ix.end[node - 1];
      else ++pos;
    }
  }
  if (nodes.empty()) Rcpp::stop("No clades meet the size bounds; lower min_tips or increase max_tips.");
  std::vector<char> seen(n, 0);
  double total = 0.0;
  for (int node : nodes) {
    if (node <= n_tip || node > n || seen[node - 1]) Rcpp::stop("invalid or repeated clade node");
    seen[node - 1] = 1; total += ix.count[node - 1];
  }
  if (total > max_total_tips) Rcpp::stop("Selected clades exceed max_total_tips; use fewer nodes or disjoint clades.");
  std::vector<int> ordered(nodes);
  std::sort(ordered.begin(), ordered.end(), [&](int a, int b) { return ix.start[a - 1] < ix.start[b - 1]; });
  bool overlap = false; int last = 0;
  for (int node : ordered) { if (ix.start[node - 1] <= last) overlap = true; last = std::max(last, ix.end[node - 1]); }
  Rcpp::List trees(nodes.size()), tips(nodes.size());
  Rcpp::LogicalVector covered(n_tip, false);
  Rcpp::IntegerVector counts(nodes.size());
  std::vector<int> mapping(n, -1);
  for (std::size_t k = 0; k < nodes.size(); ++k) {
    Rcpp::checkUserInterrupt();
    const int first = ix.start[nodes[k] - 1] - 1, end = ix.end[nodes[k] - 1];
    std::vector<int> tip_ids, internal;
    tip_ids.reserve(ix.count[nodes[k] - 1]);
    internal.reserve(end - first - ix.count[nodes[k] - 1]);
    for (int pos = first; pos < end; ++pos) {
      if (pre[pos] <= n_tip) tip_ids.push_back(pre[pos]); else internal.push_back(pre[pos]);
    }
    std::sort(tip_ids.begin(), tip_ids.end());
    const int nt = tip_ids.size(); counts[k] = nt;
    Rcpp::CharacterVector local_labels(nt);
    for (int i = 0; i < nt; ++i) {
      mapping[tip_ids[i] - 1] = i + 1; covered[tip_ids[i] - 1] = true;
      local_labels[i] = labels[tip_ids[i] - 1];
    }
    for (std::size_t i = 0; i < internal.size(); ++i) mapping[internal[i] - 1] = nt + i + 1;
    Rcpp::IntegerMatrix edge(end - first - 1, 2);
    Rcpp::NumericVector length(edge.nrow());
    for (int pos = first + 1; pos < end; ++pos) {
      const int child = pre[pos] - 1, row = pos - first - 1;
      edge(row, 0) = mapping[parent[child] - 1]; edge(row, 1) = mapping[child]; length[row] = branch[child];
    }
    Rcpp::List tree = Rcpp::List::create(Rcpp::_["edge"] = edge, Rcpp::_["edge.length"] = length,
      Rcpp::_["tip.label"] = local_labels, Rcpp::_["Nnode"] = static_cast<int>(internal.size()));
    tree.attr("class") = "phylo";
    trees[k] = tree; tips[k] = Rcpp::wrap(tip_ids);
  }
  return Rcpp::List::create(Rcpp::_["nodes"] = Rcpp::wrap(nodes), Rcpp::_["trees"] = trees,
    Rcpp::_["tips"] = tips, Rcpp::_["n_tip"] = counts, Rcpp::_["covered"] = covered,
    Rcpp::_["overlap"] = overlap, Rcpp::_["total"] = total);
}

// One C++ boundary covers all crowns and traits. Missingness groups use packed
// exact keys: hash collisions cannot merge masks, as in na_mask.cpp. Ape owns
// biological NA pruning (including its unary/root-collapse conventions).
// No RNGScope here: lazy-loading an R dependency inside a native RNGScope can
// overwrite its restored .Random.seed with stale native state on return.
// [[Rcpp::export(rng = false)]]
Rcpp::List heterogeneity_scan_cpp(const Rcpp::List& plan,
    const Rcpp::NumericMatrix& X, const std::string method,
    const int min_tips, const bool keep_cores, SEXP delta_control = R_NilValue,
    const bool signal_test = false) {
  if (method != "K" && method != "lambda" && method != "D" && method != "Delta") Rcpp::stop("unknown clade signal method");
  const Rcpp::List trees = plan["trees"], tips = plan["tips"];
  const int n_clade = trees.size(), p = X.ncol();
  if (tips.size() != trees.size() || min_tips < 3) Rcpp::stop("invalid clade scan dimensions");
  if (p < 1 || static_cast<double>(n_clade) * p > std::numeric_limits<int>::max())
    Rcpp::stop("too many clade/trait output rows");
  const int nrow = n_clade * p;
  Rcpp::NumericVector estimate(nrow, NA_REAL);
  Rcpp::NumericVector signal_p(nrow, NA_REAL);
  Rcpp::IntegerVector n_obs(nrow);
  Rcpp::CharacterVector status(nrow, "insufficient_data"), message(nrow,
    "fewer than " + std::to_string(min_tips) + " observed species");
  Rcpp::LogicalVector blocked(n_clade, false);
  Rcpp::List cores(n_clade);
  Rcpp::List data(n_clade);
  Rcpp::List delta_details(method == "Delta" ? n_clade : 0);
  Rcpp::List delta_contexts(method == "Delta" && keep_cores ? n_clade : 0);
  // Keep the namespace unloaded on complete-data scans. RObject protects the
  // callback once the first missing-data group actually needs ape's pruning.
  Rcpp::RObject keep_tip;
  Rcpp::RObject canonicalize;
  struct Group { std::vector<int> keep, columns; };
  for (int i = 0; i < n_clade; ++i) {
    Rcpp::checkUserInterrupt();
    Rcpp::List tree = trees[i];
    if (method == "Delta") {
      if (canonicalize.isNULL()) canonicalize =
        Rcpp::Environment::namespace_env("fastphylosig")[".safe_canonicalize_core"];
      tree = Rcpp::Function(canonicalize)(tree);
    }
    const Rcpp::IntegerVector tip_ids = tips[i];
    const int nt = tip_ids.size();
    const Rcpp::CharacterVector labels = tree["tip.label"];
    if (labels.size() != nt) Rcpp::stop("inconsistent clade labels/tips");
    for (int tip : tip_ids)
      if (tip < 1 || tip > X.nrow()) Rcpp::stop("clade tip outside data matrix");
    std::vector<Group> groups;
    std::unordered_map<std::string, int> masks;
    std::vector<int> eligible;
    for (int j = 0; j < p; ++j) {
      std::string mask((nt + 7) / 8, '\0');
      std::vector<int> keep;
      bool invalid = false;
      for (int tip = 0; tip < nt; ++tip) {
        const double value = X(tip_ids[tip] - 1, j);
        if (!std::isnan(value)) { keep.push_back(tip); mask[tip / 8] |= static_cast<char>(1U << (tip % 8)); }
        if (!std::isfinite(value) && !R_IsNA(value)) invalid = true;
      }
      const int row = i * p + j; n_obs[row] = keep.size();
      if (n_obs[row] < min_tips) continue;
      eligible.push_back(j);
      if (invalid) { status[row] = "invalid_trait"; message[row] = "trait contains non-finite values other than NA"; continue; }
      status[row] = "ok"; message[row] = NA_STRING;
      const auto found = masks.find(mask);
      if (found == masks.end()) {
        masks.emplace(mask, groups.size()); groups.push_back(Group{std::move(keep), std::vector<int>(1, j)});
      } else groups[found->second].columns.push_back(j);
    }
    if (eligible.empty()) continue;
    const Rcpp::IntegerMatrix edge = tree["edge"];
    const Rcpp::NumericVector length = tree["edge.length"];
    if (edge.ncol() != 2 || length.size() != edge.nrow()) Rcpp::stop("inconsistent crown edges/lengths");
    for (int e = 0; e < edge.nrow(); ++e)
      if (length[e] == 0.0 && (method != "lambda" || edge(e, 1) <= nt)) blocked[i] = true;
    if (blocked[i]) {
      for (int j : eligible) { status[i * p + j] = "invalid_subtree"; message[i * p + j] = NA_STRING; }
      continue;
    }
    const bool retain_core = keep_cores && method != "Delta";
    Rcpp::List local_cores(retain_core ? p : 0);
    Rcpp::List local_data(retain_core ? p : 0);
    Rcpp::List local_delta(method == "Delta" ? p : 0);
    Rcpp::List delta_trees(method == "Delta" ? groups.size() : 0);
    Rcpp::List delta_data(method == "Delta" ? groups.size() : 0);
    Rcpp::List delta_columns(method == "Delta" ? groups.size() : 0);
    Rcpp::List delta_rows(method == "Delta" && keep_cores ? groups.size() : 0);
    int delta_group_index = 0;
    try {
      for (const Group& g : groups) {
        Rcpp::List local_tree = tree;
        if (static_cast<int>(g.keep.size()) != nt) {
          if (keep_tip.isNULL()) keep_tip = Rcpp::Environment::namespace_env("ape")["keep.tip"];
          Rcpp::CharacterVector keep_labels(g.keep.size());
          for (std::size_t k = 0; k < g.keep.size(); ++k) keep_labels[k] = labels[g.keep[k]];
          local_tree = Rcpp::Function(keep_tip)(tree, keep_labels);
        }
        const Rcpp::CharacterVector retained = local_tree["tip.label"];
        Rcpp::NumericMatrix local(retained.size(), g.columns.size());
        // keep.tip retains original tip order; verify rather than assume it.
        if (retained.size() != static_cast<int>(g.keep.size())) Rcpp::stop("unexpected retained clade tips");
        for (int tip = 0; tip < retained.size(); ++tip) {
          if (retained[tip] != labels[g.keep[tip]]) Rcpp::stop("unexpected retained clade tip order");
        }
        for (std::size_t k = 0; k < g.columns.size(); ++k)
          for (int tip = 0; tip < retained.size(); ++tip)
            local(tip, k) = X(tip_ids[g.keep[tip]] - 1, g.columns[k]);
        if (method == "Delta") {
          delta_trees[delta_group_index] = local_tree;
          delta_data[delta_group_index] = local;
          delta_columns[delta_group_index] = Rcpp::wrap(g.columns);
          if (keep_cores) {
            Rcpp::IntegerVector rows(retained.size());
            for (int tip = 0; tip < retained.size(); ++tip) rows[tip] = tip_ids[g.keep[tip]];
            delta_rows[delta_group_index] = rows;
          }
          ++delta_group_index;
          continue;
        }
        const Rcpp::List compiled = compile_tree_cpp(local_tree["edge"], local_tree["edge.length"], Rcpp::wrap(retained.size()));
        Rcpp::List fit;
        if (method == "K") fit = fast_k_tree_batch_cpp(compiled, local, 64);
        else if (method == "lambda") {
          Rcpp::NumericVector bound(1, lambda_bound(compiled));
          fit = fast_lambda_tree_optimize_cpp(compiled, local, bound, signal_test, false, 101, 64, 1, 0.0001220703125);
        }
        const Rcpp::NumericVector values = method == "D" ? Rcpp::NumericVector(g.columns.size(), NA_REAL) : Rcpp::as<Rcpp::NumericVector>(fit[method == "K" ? "K" : "lambda"]);
        Rcpp::LogicalVector valid;
        Rcpp::CharacterVector notes;
        if (method == "lambda") { valid = fit["valid"]; notes = fit["status"]; }
        const Rcpp::NumericVector probabilities = method == "lambda" && signal_test ?
          Rcpp::as<Rcpp::NumericVector>(fit["P"]) : Rcpp::NumericVector(0);
        for (std::size_t k = 0; k < g.columns.size(); ++k) {
          const int j = g.columns[k], row = i * p + j;
          estimate[row] = values[k];
          if (method == "lambda" && signal_test && valid[k]) signal_p[row] = probabilities[k];
          if (keep_cores) {
            local_cores[j] = compiled;
            if (method != "K" || signal_test) local_data[j] = Rcpp::NumericVector(local.column(k));
          }
          if (method == "D") {
            int ones = 0;
            for (int tip = 0; tip < local.nrow(); ++tip) {
              if (local(tip, k) != 0.0 && local(tip, k) != 1.0) Rcpp::stop("D requires binary states coded as 0/1");
              ones += local(tip, k) == 1.0;
            }
            const Rcpp::IntegerVector ptr = compiled["child_ptr"];
            bool binary = true;
            for (int node = local.nrow(); node + 1 < ptr.size(); ++node)
              binary = binary && ptr[node + 1] - ptr[node] == 2;
            if (!binary) { status[row] = "invalid_subtree"; message[row] = "D heterogeneity requires bifurcating crowns; resolve polytomies explicitly"; }
            else if (ones == 0 || ones == local.nrow()) { status[row] = "undefined"; message[row] = "D is undefined for a single-state crown"; }
          }
          if ((method == "K" && !std::isfinite(values[k])) || (method == "lambda" && !valid[k])) {
            status[row] = "undefined";
            message[row] = method == "K" ? Rcpp::String("K is undefined for a constant or numerically degenerate trait") : Rcpp::String(notes[k]);
          }
        }
      }
      if (method == "Delta") {
        Rcpp::Function batch = Rcpp::Environment::namespace_env("fastphylosig")[".heterogeneity_delta_batch"];
        const Rcpp::List fit = batch(delta_trees, delta_data, delta_columns, delta_control, p,
          keep_cores ? Rcpp::RObject(delta_rows) : Rcpp::RObject(R_NilValue));
        const Rcpp::NumericVector values = fit["estimate"];
        const Rcpp::CharacterVector statuses = fit["status"], messages = fit["message"];
        const Rcpp::List details = fit["details"];
        for (const Group& g : groups) for (int j : g.columns) {
          const int row = i*p+j;
          estimate[row] = values[j]; status[row] = statuses[j];
          message[row] = messages[j]; local_delta[j] = details[j];
        }
        if (keep_cores) delta_contexts[i] = fit["contexts"];
      }
    } catch (const Rcpp::internal::InterruptedException&) {
      throw;
    } catch (const std::exception& e) {
      for (int j : eligible) { const int row = i * p + j; estimate[row] = signal_p[row] = NA_REAL; status[row] = "invalid_subtree"; message[row] = e.what(); }
    }
    if (retain_core) { cores[i] = local_cores; data[i] = local_data; }
    if (method == "Delta") delta_details[i] = local_delta;
  }
  return Rcpp::List::create(Rcpp::_["estimate"] = estimate, Rcpp::_["n_obs"] = n_obs,
    Rcpp::_["status"] = status, Rcpp::_["message"] = message,
    Rcpp::_["blocked"] = blocked, Rcpp::_["cores"] = cores, Rcpp::_["data"] = data,
    Rcpp::_["delta_details"] = delta_details, Rcpp::_["delta_contexts"] = delta_contexts,
    Rcpp::_["P_signal"] = signal_p);
}

// Local K tests reuse the scan's exact NA groups and the production bounded-
// memory permutation kernel. No R callback/preflight/NA pruning is repeated.
// This separate call preserves the established overall/heterogeneity/local
// RNG ordering; performing tests in the initial scan would change that stream.
// [[Rcpp::export]]
Rcpp::List heterogeneity_k_signal_cpp(const Rcpp::List& cores,
    const Rcpp::List& data, const Rcpp::LogicalVector& eligible,
    const int nsim, const int simulation_chunk) {
  const int nc = cores.size();
  if (nc < 1 || data.size() != nc || eligible.size() % nc != 0 || nsim < 1)
    Rcpp::stop("invalid local K signal dimensions");
  const int p = eligible.size() / nc;
  Rcpp::NumericVector probability(eligible.size(), NA_REAL);
  Rcpp::IntegerVector successful(eligible.size(), 0);
  for (int i = 0; i < nc; ++i) {
    Rcpp::checkUserInterrupt();
    if (Rf_isNull(cores[i])) continue;
    const Rcpp::List local_cores = cores[i], local_data = data[i];
    if (local_cores.size() != p || local_data.size() != p)
      Rcpp::stop("inconsistent local K scan cache");
    std::unordered_map<SEXP, int> by_core;
    std::vector<std::vector<int> > columns;
    for (int j = 0; j < p; ++j) if (eligible[i * p + j] == TRUE) {
      SEXP key = local_cores[j];
      if (key == R_NilValue) Rcpp::stop("missing eligible local K core");
      const auto found = by_core.find(key);
      if (found == by_core.end()) {
        by_core.emplace(key, columns.size()); columns.push_back(std::vector<int>(1, j));
      } else columns[found->second].push_back(j);
    }
    for (const std::vector<int>& ids : columns) {
      const Rcpp::List core = local_cores[ids[0]];
      const int nt = Rcpp::as<int>(core["n_tip"]);
      Rcpp::NumericMatrix X(nt, ids.size());
      for (std::size_t k = 0; k < ids.size(); ++k) {
        const Rcpp::NumericVector values = local_data[ids[k]];
        if (values.size() != nt) Rcpp::stop("inconsistent local K trait cache");
        std::copy(values.begin(), values.end(), X.begin() + k * nt);
      }
      const Rcpp::List fit = fast_k_tree_permutation_cpp(core, X, nsim,
        R_NilValue, 64, false, true, 1, simulation_chunk);
      const Rcpp::NumericVector P = fit["P"];
      const Rcpp::IntegerVector n = fit["nsim_successful"];
      for (std::size_t k = 0; k < ids.size(); ++k) {
        probability[i * p + ids[k]] = P[k]; successful[i * p + ids[k]] = n[k];
      }
    }
  }
  return Rcpp::List::create(Rcpp::_["P"] = probability, Rcpp::_["nsim_successful"] = successful);
}

// All signal-reference bootstraps share the same plus-one tail and MCSE.
// [[Rcpp::export(rng = false)]]
Rcpp::List heterogeneity_bootstrap_summary_cpp(const std::vector<double>& draws,
    const double observed, const bool keep_null) {
  return fastphylosig::hetero_bootstrap_summary(draws, observed, keep_null);
}

// Unknown tail indicators remain in [0,1]. These algebraic bounds retain the
// requested denominator; they are not confidence intervals or imputed Delta.
// [[Rcpp::export(rng = false)]]
Rcpp::List heterogeneity_bootstrap_bounds_cpp(const std::vector<double>& draws,
                                             const double observed) {
  int finite = 0, exceed = 0;
  for (double value : draws) if (std::isfinite(value)) {
    ++finite;
    if (value >= observed) ++exceed;
  }
  const int n = draws.size();
  const bool valid = n >= 2 && std::isfinite(observed);
  return Rcpp::List::create(
    Rcpp::_["finite_exceedance_count"] = valid ? exceed : NA_INTEGER,
    Rcpp::_["P_heterogeneity_lower"] = valid ? (1.0 + exceed) / (n + 1.0) : NA_REAL,
    Rcpp::_["P_heterogeneity_upper"] = valid ? (1.0 + exceed + n - finite) / (n + 1.0) : NA_REAL);
}

// One pass replaces per-category R vapply closures on every reference fit.
// Keep isTRUE semantics, including absent/NA flags, and never enter RNG scope.
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector heterogeneity_delta_flags_cpp(const Rcpp::List& full,
                                                 const Rcpp::List& local) {
  Rcpp::IntegerVector flags(8);
  flags.names() = Rcpp::CharacterVector::create("nsim_MCMC_warning", "nsim_optimizer_warning",
    "nsim_numerical_recovery", "nsim_entropy_recovery", "nsim_other_warning",
    "nsim_calculation_failure", "nsim_unavailable", "nsim_single_state");
  const char* fields[] = {"MCMC_warning", "optimizer_warning", "numerical_recovery",
    "entropy_recovery", "other_warning", "calculation_failed"};
  auto is_true = [](SEXP x) {
    return TYPEOF(x) == LGLSXP && Rf_xlength(x) == 1 && LOGICAL(x)[0] == TRUE;
  };
  int retried = 0, recovered = 0;
  auto accumulate = [&](const Rcpp::List& fit) {
    const Rcpp::List details = fit["details"];
    for (int i = 0; i < 6; ++i)
      if (details.containsElementNamed(fields[i]) && is_true(details[fields[i]])) flags[i] = 1;
    if (!std::isfinite(Rcpp::as<double>(fit["estimate"]))) flags[6] = 1;
    if (fit.containsElementNamed("single_state") && is_true(fit["single_state"])) flags[7] = 1;
    if (details.containsElementNamed("optimizer_retry") && is_true(details["optimizer_retry"])) retried = 1;
    if (details.containsElementNamed("optimizer_recovered") && is_true(details["optimizer_recovered"])) recovered = 1;
  };
  accumulate(full);
  for (R_xlen_t i = 0; i < local.size(); ++i) accumulate(Rcpp::as<Rcpp::List>(local[i]));
  if (retried) flags.attr("optimizer_retry") = Rcpp::IntegerVector::create(
    Rcpp::_["nsim_optimizer_retry"]=retried,Rcpp::_["nsim_optimizer_recovered"]=recovered);
  return flags;
}

// BM generation, prepared K evaluation and online dispersion accumulation
// share one C++ call. Trait/clade draw ordering remains the R implementation's.
// [[Rcpp::export]]
Rcpp::List heterogeneity_calibrate_cpp(const Rcpp::List& compiled,
    const double observed, const int nsim, const int simulation_chunk,
    const bool keep_null) {
  const int n_clade = compiled.size();
  if (n_clade < 2 || nsim < 2 || !std::isfinite(observed)) Rcpp::stop("invalid BM calibration dimensions");
  std::vector<double> mean(nsim, 0.0), m2(nsim, 0.0);
  Rcpp::NumericVector means(n_clade, NA_REAL), lower(n_clade, NA_REAL), upper(n_clade, NA_REAL);
  for (int i = 0; i < n_clade; ++i) {
    Rcpp::checkUserInterrupt();
    const Rcpp::NumericVector draw = heterogeneity_bm_k_cpp(compiled[i], nsim, simulation_chunk);
    std::vector<double> values(draw.begin(), draw.end());
    bool all_good = true;
    for (int b = 0; b < nsim; ++b) {
      const bool good = std::isfinite(values[b]) && values[b] > 0.0;
      all_good = all_good && good;
      const double log_k = good ? std::log(values[b]) : NA_REAL;
      const double delta = log_k - mean[b];
      mean[b] += delta / (i + 1.0);
      m2[b] += delta * (log_k - mean[b]);
    }
    if (all_good) {
      means[i] = finite_mean(values);
      std::sort(values.begin(), values.end());
      lower[i] = fastphylosig::quantile_type7_sorted(values, 0.025);
      upper[i] = fastphylosig::quantile_type7_sorted(values, 0.975);
    }
  }
  Rcpp::NumericVector null(nsim);
  std::vector<double> sorted;
  sorted.reserve(nsim);
  int successful = 0, count = 0;
  for (int b = 0; b < nsim; ++b) {
    null[b] = std::isfinite(m2[b]) ? std::sqrt(std::max(0.0, m2[b] / (n_clade - 1.0))) : NA_REAL;
    if (std::isfinite(null[b])) { ++successful; if (null[b] >= observed) ++count; sorted.push_back(null[b]); }
  }
  double p = NA_REAL, mcse = NA_REAL, bm_mean = NA_REAL, bm_upper = NA_REAL;
  if (successful == nsim) {
    p = (1.0 + count) / (nsim + 1.0);
    const double q = static_cast<double>(count) / nsim;
    mcse = std::sqrt(nsim * q * (1.0 - q)) / (nsim + 1.0);
    bm_mean = finite_mean(sorted);
    std::sort(sorted.begin(), sorted.end());
    bm_upper = fastphylosig::quantile_type7_sorted(sorted, 0.95);
  }
  return Rcpp::List::create(Rcpp::_["clade_mean"] = means, Rcpp::_["clade_lower"] = lower,
    Rcpp::_["clade_upper"] = upper, Rcpp::_["nsim_successful"] = successful,
    Rcpp::_["exceedance_count"] = successful == nsim ? count : NA_INTEGER,
    Rcpp::_["P_BM"] = p, Rcpp::_["MCSE_P_BM"] = mcse,
    Rcpp::_["BM_mean"] = bm_mean, Rcpp::_["BM_upper"] = bm_upper,
    Rcpp::_["null"] = keep_null ? Rcpp::RObject(null) : Rcpp::RObject(R_NilValue));
}

