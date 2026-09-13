# fastphylosig 0.2.0

`fastphylosig` 0.2.0 provides a faster and more robust workflow for
phylogenetic signal analysis while preserving the established statistical
definitions and public interface.

## Performance

- Substantially improved preparation for large phylogenetic trees.
- Faster repeated analyses when using prepared tree contexts.
- Efficient multi-trait and batch analysis workflows.

## Correctness and Robustness

- Canonicalization Contract V2 provides representation-invariant tree
  handling for supported equivalent tree representations.
- Exact mutation detection protects computationally relevant tree fields.
- A stale-cache correctness fix prevents reuse of invalid retained-subtree
  inspection results.

## Persistence and Reproducibility

- V2 prepared contexts can be persisted with `saveRDS()` and restored with
  `readRDS()`.
- Legacy and unknown prepared-context schemas are rejected explicitly.
- Existing random-number-generation and statistical contracts are preserved.

## API

- The public API remains stable.
- Deprecated `match_phylo_data()` is retained for compatibility.
- Estimator definitions are unchanged.

## Validation Scope

Local qualification was completed on Windows 11 with R 4.6.1. Observed-D
parity was confirmed on representative fixture(s). This release note does not
claim qualification on all platforms, CRAN readiness, universal speedups over
reference implementations, or complete `caper::phylo.d()` output parity.
