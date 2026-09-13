# fastphylosig 0.2.0 Final Public Benchmark Protocol

Status: audit-only, single-machine, single-platform qualification.

This directory contains the final user-facing performance and memory evidence
for the exact source commit recorded in the result files.  It is deliberately
separate from the numerical-kernel candidate benchmarks.  No cross-platform
qualification is performed here.

## Scope and invariants

- The benchmark is run from one clean ASCII staging copy of the current source
  tree and one private R library built from that same source.
- Each fixture is generated once from a fixed seed and reused in every paired
  call.  Tree tip order, branch lengths, traits, and matrix columns are fixed.
- Calls are serialized, warmed up, and alternated within each pair.  Formal
  completed cells use at least 10 repeats.  No result is copied from Stage 1,
  Stage 2, or any earlier benchmark.
- A correctness guard runs before timings for each fixture.  A failed guard
  invalidates that cell; it is recorded as `CORRECTNESS_FAIL`, never as a
  speed result.
- Timing is wall-clock elapsed time from `system.time()` and includes the
  public call path specified by the workload.  It is not peak RSS.
- The benchmark does not modify package source, tests, package metadata, or
  public API.

## Workloads

### A. Single-trait raw

`fast_k(tree, x, test = FALSE)` and `fast_lambda(tree, x)` on a raw `phylo`.
The principal sizes are `n = 500, 2000, 5000`; `n = 100` is used for a
reference guard where the reference implementation is feasible.

### B. Prepared repeated analysis

One `ctx <- prepare_tree(tree)` is created outside the timed call.  The same
prepared context is reused for repeated `fast_k(ctx, x, test = FALSE)` and
`fast_lambda(ctx, x)` calls.  Preparation time is reported separately and is
not included in the prepared-call time.

### C. Batch/matrix traits

`fast_k(ctx, X, test = FALSE)` and the equivalent `fast_signal(ctx, data = X,
method = "K", test = FALSE)` are measured for `n = 500` and 8, 32, and 100
continuous traits.  Matrix correctness requires one-column parity with the
corresponding single-trait call and row names exactly equal to tree tips.

### D. Permutation-heavy K and D

Prepared `fast_k(ctx, x, test = TRUE, nsim = 999)` and prepared
`fast_d(ctx, binary, test = TRUE, nsim = 999, return_sim = FALSE)` are measured
for `n = 500, 2000, 5000`.  `ncores = 1` is authoritative.  A separate
`ncores = 2` observation is collected only when it fits the per-cell budget;
it does not create a new reproducibility promise.

### E. Large-tree scaling and memory

Prepared contexts are characterized at `n = 500, 1000, 5000, 10000, 20000`.
For each size the evidence records `object.size(ctx)`, serialized RDS byte
size, save time, and a fresh-process `readRDS()` plus first validation time.
These are measured object/file quantities only and are never described as
peak RSS.

Small reference comparisons use `phytools::phylosig()` for K/lambda and
`caper::phylo.d()` for D only when the package interface accepts the fixed
fixture within the declared budget.  An unavailable or over-budget reference
is recorded as `UNAVAILABLE` or `RESOURCE_LIMIT`; no old timing or extrapolated
speedup is substituted.

## Formal timing protocol

- Warmup: two untimed calls per route and fixture.
- Formal repeats: 10 paired repeats per completed cell.
- Order: route A then B on odd repetitions and B then A on even repetitions;
  a fresh fixed seed is restored before every stochastic call.
- Runtime budget: 120 seconds per light cell, 600 seconds per heavy
  permutation cell, and 300 seconds per memory size.  A timeout or process
  failure produces a row with `RESOURCE_LIMIT` and no timing summary.
- A formal run must finish all declared cells to be called complete.  A
  deliberately narrowed run is useful for smoke testing only and is marked
  `INCOMPLETE` in `benchmark_status.csv`.
- The benchmark writes raw repeat rows as it progresses so an interrupted run
  cannot be mistaken for a complete result.

## Correctness guard

Before a timing cell is accepted, the runner checks raw/prepared parity for
the relevant method, finite/status-valid output, and deterministic replay for
K and D under the same seed and `ncores = 1`.  Matrix cells additionally check
one-column equality.  Reference rows are checked against the frozen
method-specific estimate on the small reference fixture; any mismatch is
reported for investigation and is not used to claim performance.

## Reporting

`final_benchmark.csv` contains one row per call and repeat, including method,
workload, shape, n, trait count, nsim, route, elapsed seconds, status, and
reference metadata.  `memory_characterization.csv` contains only measured
memory lifecycle quantities.  `benchmark_status.csv` is the fail-closed
summary of guards, completed cells, resource limits, and provenance.

The runner records the exact Git commit, package version, R version, compiler
toolchain, BLAS/LAPACK, OpenMP setting, locale, operating system, and hardware
in the CSV metadata and in `MEMORY_CHARACTERIZATION.md`.

