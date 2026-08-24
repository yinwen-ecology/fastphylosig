# Post-0.1.0 technical debt

This internal file records work that is not part of the 0.1.0 public contract.
It is excluded from source-package builds by `.Rbuildignore`.

## Prepared-tree integrity

The current prepared context is an R snapshot. Its reuse check fingerprints
tip labels, the edge matrix, branch lengths, and `Nnode` of `ctx$tree`; it neither tracks
later mutation of the original `phylo` object nor covers node labels or other
attributes. A future integrity contract should define the complete set of
calculation-relevant fields, use a compact versioned digest, test every covered
mutation, and state whether unsupported attribute mutation is ignored or
rejected.

## Input identity and ordering

Named inputs are matched by unique species identifiers and reordered to tree-tip
order. Unnamed tables can currently inherit tree-tip order when dimensions
match. Shared table validation now covers duplicate names, default row names,
and wrong-length unnamed input before method-specific conversion. Consider
requiring explicit names for all non-trivial workflows and extend cross-method
tests for scrambled names and controlled permutation/state matrix column order.

## Insufficient data and topology

Delta currently stops when fewer than two species are shared globally, but
returns per-trait insufficient-data rows after NA pruning. D accepts rooted
polytomies through a compatibility path and rejects unary nodes, including
those exposed by pruning. Consolidate these rules in one public validation
schema and add end-to-end tests for global and per-trait failures.

## Monte Carlo results and retention

K, D, and Delta do not expose one uniform retention API. D requires both
`return_sim` and `keep_null`; Delta has only `return_sim`. `MCSE_P` names and
availability also vary, with D exposing two null-specific values. A future API
may normalize these without removing legacy fields. Acceptance tests should
cover disabled tests, partial failures, zero successful draws, boundary P
values, compact results, and plotting preconditions.

## Progress and warnings

`progress = FALSE` suppresses progress and matching messages but deliberately
does not suppress warnings or errors. The naming can invite a stronger
interpretation. Consider a separate condition-policy control only if it can
preserve standard R warning semantics; meanwhile retain tests that warnings are
not hidden by the progress switch.

## Parallel reproducibility

Parallel replay is conditional on a fixed worker count. Delta's PSOCK worker
streams are repartitioned when `ncores` changes, so draw-for-draw null and
P/MCSE equality is not promised across worker counts. Future work should either
assign deterministic RNG substreams per permutation index or make the current
worker-count dependence machine-readable. Validate serial/parallel behavior on
all supported platforms before strengthening the promise.

## Delta controls

The 0.1.0 defaults are `test = FALSE`, `nsim = 1000`, `return_sim = test`,
`mcmc_sim = 10000`, `thin = 10`, `burn = 100`, and `model = "ARD"`. They are
computational defaults, not guarantees of adequate ESS, split R-hat, or Monte
Carlo precision. Revisit defaults only with retained calibration evidence and a
versioned behavior change; do not silently tune runs based on diagnostics.
