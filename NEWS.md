# fastphylosig 0.1.0

- Consolidated repeated trait-table validation and C++ quantile code, removed
  obsolete ACE development timing hooks, and reorganized the English and
  Chinese guides around a one-command start and an explicit advanced workflow.
- Declared compatibility with R 4.1.0 and newer. Windows source builds use the
  matching Rtools toolchain; OpenMP remains optional with a serial fallback.
- Added the method-explicit `fast_signal(tree, data, method, ...)` dispatcher
  for K, lambda, D, and Delta, with public `fast_k()` and `fast_lambda()`.
- Kept `fast_d()`, `fast_delta()`, and `fast_ace()` as specialist entry points;
  methods are never inferred from trait values and no public numerical-engine
  selector was added.
- Made `match_tree_data()` the documented matching API and retained
  `match_phylo_data()` only as a deprecated compatibility alias.
- Added shared tree/data preparation, compact workflow metadata, packed NA-mask
  grouping, retained-subtree validation, and bounded prepared-tree caches while
  preserving specialist result classes and statistical contracts.
- Restricted automatic tree normalization to representation-only changes. It
  never invents roots or branch lengths, jitters zero branches, resolves a
  biological polytomy, ultrametricizes a tree, or silently removes tips.
- Added consistent per-trait `status`, `message`, and `note` output for
  non-finite, insufficient, single-state, constant, and undefined traits.
- Added a shared tree-requirement registry, actionable method-specific errors,
  direct ACE topology/branch guards, root-sensitivity fixtures, and explicit
  zero-terminal/zero-internal branch contracts.
- Added controlled K permutation and D random/Brownian null fixtures, lambda
  fixed/optimized pruning tests, unified-versus-specialist equivalence tests,
  large-offset stability tests, and deterministic raw/prepared comparisons.
- Clarified that `fast_delta()` implements a single-tree statistic with local
  ACE and RNG paths. Cross-language parity and the multiple-tree extension are
  outside the package contract.
- Replaced the large-tree NA-subset cache key with the packed mask key, so
  3000/5000-tip analyses no longer hit R's 10,000-byte environment-name limit.
- Set Delta defaults to `mcmc_sim = 10000`, `thin = 10`, and `burn = 100`.
  Added one aggregated retained-sample warning and one aggregated ESS/R-hat/
  MCSE warning per batch; warnings do not rerun or lengthen chains.
- Reused permutation-invariant ACE workspaces inside Delta and skipped unused
  post-fit Hessian/SE work only in the private Delta ancestral-state helper.
  Public `fast_ace()` still computes and returns standard errors.
- Removed generated build artifacts and obsolete development-only benchmark
  records from the working tree; retained only concise release and reference
  notes needed to explain the current package contract.
- Final release-contract patch: ambiguous simultaneous `x`/`X` inputs now fail
  for D and Delta, ACE rejects invalid `method` values clearly, and unused
  Delta `nsim` input is ignored when `test = FALSE`. Prepared-tree fingerprint
  validation now also covers `Nnode`; targeted regressions and the full local
  check were rerun without changing numerical kernels.

# fastphylosig 0.0.4

- Added `check_tree()` for read-only, method-specific readiness checks for K,
  lambda, D, and Delta. Diagnostics distinguish global tree defects from the
  topology and branch-length contracts of individual statistics.
- Added conservative `resolve_tree()`. It standardizes edge order and internal
  node numbering, rechecks the result, records provenance, and refuses unsafe
  repairs such as random polytomy resolution, branch jitter, or tip deletion.
- Added method-level `progress` messages and stable elapsed-time metadata to
  `fast_signal()`, `fast_d()`, `fast_delta()`, and `fast_ace()`.
- Added explicit R-boundary topology guards for D and ACE. D now rejects unary
  nodes or noncanonical roots instead of silently discarding a branch; ACE now
  rejects mixed unary/polytomous topologies that previously could pass a node-
  count-only check.
- Added boundary and numerical-contract tests for ACE transition evaluation;
  production tolerances and public results remain unchanged.

# fastphylosig 0.0.3

- Made the validated tree-pruning implementations the automatic production
  paths for Blomberg's K and Pagel's lambda. Public callers no longer select
  dense, tree, spectral, or validation engines.
- Added streaming, chunked K randomization with Monte Carlo uncertainty and
  optional null retention, without allocating a full tip-by-simulation matrix.
- Added a tree-linear fixed-lambda likelihood and a batched C++ Brent optimizer
  returning lambda, log-likelihood, GLS mean, sigma2, LR, P, and optional
  likelihood profiles without production VCV/eigendecomposition work.
- Added explicit Monte Carlo uncertainty and successful/failed simulation
  counts for both Fritz-Purvis D null hypotheses.
- Retained Delta ESS, split R-hat, MCSE, iteration, and permutation diagnostics.
- Profiled ACE and Delta before optimization. ACE branch-transition arithmetic
  and repeated Delta permutation calculations are the demonstrated hotspots;
  small PSOCK workloads can be slower than serial execution.
- Added production and component benchmarks plus expanded deterministic,
  stochastic, boundary, offset, topology, and cache regression tests.

# fastphylosig 0.0.2

## Tree preparation and cache

- `prepare_tree()` now compiles O(n) topology and traversal data eagerly while
  creating VCV, Cholesky, and lambda spectral payloads only on demand.
- Added separate structural and numerical caches, a byte-budget numerical LRU,
  packed C++ grouping of trait-wise missing-data masks, and `cache_info()`.
- D, Delta, ACE, and tree-engine K workflows no longer request dense VCV data.

## Continuous signal

- Kept dense/Cholesky K as the default compatibility engine.
- Added opt-in `engine = "tree"` for O(n p) observed Blomberg K and
  `engine = "validate"` for explicit tree-versus-dense comparisons.
- Lambda remains on the validated spectral engine with a dense numerical
  fallback; no unvalidated tree-lambda implementation is exposed.

## D and Delta

- Added bounded, fused D null simulation and traversal with configurable
  `chunk_size` and `keep_null`.
- Added Delta chain means/SDs, ESS, split R-hat, covariance-aware MCSE for
  Delta, permutation P MCSE, successful-count reporting, and safe RNG metadata.

## Validation and documentation

- Added deterministic golden contracts, compiled-tree, packed-mask, streaming,
  tree-K, cache, and Delta-diagnostic tests.
- Added a unified environment-controlled benchmark matrix and rewrote the user
  and expert documentation around explicit compatibility boundaries.
