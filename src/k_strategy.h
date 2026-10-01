#ifndef FASTPHYLOSIG_K_STRATEGY_H
#define FASTPHYLOSIG_K_STRATEGY_H

// Internal strategy selector for the K permutation kernel.
//
// Stage 4 moved the Stage 3C prototype evaluator into src/.  The prototype
// exposed `mode` (production / stage1 / fusion / combined) and an explicit
// `block` as test knobs.  Production must not expose those, so this header
// replaces them with a pure function of the workload:
//
//     select(n_tip, n_trait, nsim, memory_budget_bytes) -> (hybrid, fuse, block)
//
// Every constant below is traceable to a Stage 3C measurement; the provenance
// is recorded next to it.  Nothing here touches the public API, the statistical
// definition, or the RNG stream: the selector only chooses WHICH evaluator
// implementation runs, and all implementations are contract-equivalent (see
// the numerical contract in k_permutation.cpp).
//
// This header deliberately depends on nothing but the C++ standard library, so
// the policy can be unit-tested without R, Rcpp, or a compiled package.

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace kstrategy {

// ---------------------------------------------------------------------------
// Workspace model.
//
// For one trait chunk of width `chunk` and a block of `block` permutations the
// evaluator holds:
//
//     message   n_total * block * chunk doubles
//     state     n_total * block * chunk doubles
//     baseline  block * chunk doubles
//     delta     block * chunk doubles
//     out       nperm * chunk doubles   (nperm <= block)
//
// Only the first two are large.  The Stage 1 hybrid does NOT reduce this: the
// residual passes still need message and state, so the combined kernel's
// workspace equals the fusion-only kernel's.  One rule therefore serves both
// fused modes.
// ---------------------------------------------------------------------------
inline double workspace_bytes(const int n_total, const int block,
                              const int chunk) {
  const double cells = static_cast<double>(n_total) *
    static_cast<double>(block) * static_cast<double>(chunk);
  return 2.0 * cells * 8.0 +
    2.0 * static_cast<double>(block) * static_cast<double>(chunk) * 8.0 +
    static_cast<double>(block) * static_cast<double>(chunk) * 8.0;
}

// Block ladder.  Descending so the first fitting rung wins.
inline const int* block_ladder() {
  static const int ladder[] = {64, 32, 16, 8, 4, 1};
  return ladder;
}
inline const int block_ladder_size() { return 6; }

// Largest rung that fits `budget_bytes`.  Returns 1 when even the smallest
// fused block does not fit, which makes the caller fall back to no fusion.
inline int adaptive_block(const int n_total, const int chunk,
                          const double budget_bytes) {
  const int* ladder = block_ladder();
  for (int k = 0; k < block_ladder_size(); ++k) {
    const int b = ladder[k];
    if (workspace_bytes(n_total, b, chunk) <= budget_bytes) return b;
  }
  return 1;
}

// ---------------------------------------------------------------------------
// Gate constants.
// ---------------------------------------------------------------------------

// Default memory budget for the block workspace: 64 MiB.  This is the value the
// Stage 3C "budget" grid used, so the adaptive rule is reproduced exactly.
inline double default_budget_bytes() { return 64.0 * 1048576.0; }

// Package default for `trait_chunk` (see fast_k()'s signature).
inline int default_trait_chunk() { return 64; }

// Legacy diagnostic retained for compatibility with the internal strategy
// report.  Fusion policy is now based on the total trait count and tree size;
// this historical chunk threshold is not consulted by the selector.
inline int fuse_max_chunk() { return 10; }

// Hybrid has no selector-side nsim threshold.  The production entry point
// validates nsim, while direct selector inspection still reports hybrid=true
// for every input value.
inline int hybrid_min_nsim() { return 0; }

// Conservative v0.3.0 fusion window.  This gate uses total n_trait, never the
// per-pass trait chunk, so callers cannot bypass the safety window by lowering
// trait_chunk.
inline int fusion_max_traits(const int n_tip) {
  if (n_tip <= 500) return 10;
  if (n_tip <= 5000) return 5;
  if (n_tip <= 10000) return 5;
  return 0;
}

// ---------------------------------------------------------------------------
// Strategy.
// ---------------------------------------------------------------------------
struct Strategy {
  // --- the contract required by the Stage 4 brief ---
  bool hybrid = false;   // use the Stage 1 dot-product GLS offset
  bool fuse = false;     // evaluate a block of permutations per traversal
  int block = 1;         // permutations per block (1 when !fuse)

  // --- audit trail, not part of the contract ---
  int chunk = 1;             // min(trait_chunk, n_trait)
  int n_total = 0;           // tree nodes actually used by the rule
  int n_total_estimated = 0; // true when the caller did not supply one
  int rule_block = 1;        // rung chosen by the memory rule, before capping
  bool block_capped_to_nsim = false;
  bool block_capped_to_traits = false;
  bool budget_invalid = false;
  bool gate_chunk = false;   // legacy field; no chunk-width gate is applied
  bool gate_traits = false;  // total n_trait exceeds the tree-size window
  bool gate_tree_size = false; // n_tip > 10000
  bool gate_workspace = false; // even a one-permutation block exceeds budget
  bool gate_nsim = false;    // fusion refused because nsim < 2
  bool hybrid_gate_nsim = false;
  double budget_bytes = 0.0;
  double workspace_bytes = 0.0;
};

// ---------------------------------------------------------------------------
// The selector.
//
// The Stage 4 brief fixes the input vector at (n_tip, n_trait, nsim, memory
// budget).  Two optional refinements are offered so the production call site
// can be exact rather than estimated; both default to "not supplied" so the
// four-argument form is the literal contract:
//
//   n_total_actual  the real node count from the compiled tree.  When omitted,
//                   2*n_tip - 1 is used, which is exact for a rooted
//                   bifurcating tree and an over-estimate for a tree with
//                   multifurcations (a larger n_total means a smaller block,
//                   so the estimate errs towards less memory).
//   trait_chunk     the caller's chunk width.  When omitted, the package
//                   default (64) is used for workspace estimation only; the
//                   fusion benefit gate always uses total n_trait.
// ---------------------------------------------------------------------------
inline Strategy select(const int n_tip, const int n_trait, const int nsim,
                       const double memory_budget_bytes,
                       const int n_total_actual = 0,
                       const int trait_chunk = 0) {
  Strategy s;

  // --- normalise the inputs -------------------------------------------------
  const int p = std::max(1, n_trait);
  const int requested_chunk = trait_chunk > 0 ? trait_chunk
                                              : default_trait_chunk();
  s.chunk = std::max(1, std::min(requested_chunk, p));

  if (n_total_actual > 0) {
    s.n_total = n_total_actual;
    s.n_total_estimated = 0;
  } else {
    const int tips = std::max(2, n_tip);
    s.n_total = std::max(tips + 1, 2 * tips - 1);
    s.n_total_estimated = 1;
  }

  s.budget_bytes = memory_budget_bytes;
  if (!(memory_budget_bytes > 0.0) || !std::isfinite(memory_budget_bytes)) {
    s.budget_bytes = default_budget_bytes();
    s.budget_invalid = true;
  }

  const int steps = std::max(0, nsim);

  // --- hybrid ---------------------------------------------------------------
  // Hybrid is unconditional at the strategy level.  nsim validity is checked
  // by the production entry point, not used to choose the evaluator.
  s.hybrid = true;
  s.hybrid_gate_nsim = false;

  // --- fuse -----------------------------------------------------------------
  s.gate_nsim = steps < 2;
  s.gate_chunk = false;
  s.gate_traits = p > fusion_max_traits(n_tip);
  s.gate_tree_size = n_tip > 10000;
  s.gate_workspace = workspace_bytes(s.n_total, 1, s.chunk) > s.budget_bytes;
  s.fuse = n_tip >= 2 && !s.gate_traits && !s.gate_tree_size &&
    !s.gate_workspace && !s.gate_nsim;

  // --- block ----------------------------------------------------------------
  if (s.fuse) {
    s.rule_block = adaptive_block(s.n_total, s.chunk, s.budget_bytes);
    const int within_nsim = std::min(s.rule_block, steps);
    s.block_capped_to_nsim = s.rule_block > steps;
    const int trait_block_cap = p <= 2 ? 32 : (p <= 10 ? 4 : 1);
    s.block = std::min(within_nsim, trait_block_cap);
    s.block_capped_to_traits = within_nsim > trait_block_cap;
    if (s.block < 2) s.fuse = false;  // no block to fuse
  }
  if (!s.fuse) {
    s.rule_block = 1;
    s.block = 1;
  }
  s.workspace_bytes = workspace_bytes(s.n_total, s.block, s.chunk);
  return s;
}

// ---------------------------------------------------------------------------
// Thread policy.
//
// Fusion is a single-thread memory-locality transformation: it makes one
// traversal serve `block` permutations.  Under OpenMP the production kernel
// already hides that latency by running permutations concurrently, and a fused
// block of 64 permutations cannot be split across 8 threads without turning the
// block into a scheduling unit and paying a reduction per block.  Stage 3C
// measured the fused kernel single-threaded only, so enabling it under
// OpenMP would be an unmeasured change.  Fusion is therefore disabled when more
// than one thread is in use; the hybrid evaluator, which is a per-permutation
// substitution with the identical accumulation order, stays enabled.
// ---------------------------------------------------------------------------
inline Strategy apply_thread_policy(const Strategy& in, const int n_threads) {
  Strategy s = in;
  if (n_threads > 1 && s.fuse) {
    s.fuse = false;
    s.rule_block = 1;
    s.block = 1;
    s.workspace_bytes = workspace_bytes(s.n_total, 1, s.chunk);
  }
  return s;
}

}  // namespace kstrategy

#endif
