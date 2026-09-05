# fastphylosig 0.2.0 Stage 2B1 Expert Review

## Decision

`PREPARATION_ENGINE_STAGE1 = ACCEPTED_FROZEN`

`RAW_ROUTE_SINGLE_PREPARATION_BOUNDARY = ACCEPTED_FROZEN`

Stage 2B1 meets every correctness, package-check, call-count, and performance
gate. The candidate removes repeated raw-route preparation without weakening
tree validation or changing any estimator, numerical formula, tolerance, RNG,
thread policy, public argument, or public result structure.

Stage 2B2 was not started.

## Provenance

- Frozen before commit: `017f39df1aa520f71fc2a0c5647bd598ff23e824`
- Stage 2B1 candidate commit: `e265b47cd7abbb3c727ec466998c67e8eecdbc98`
- Package version: `0.2.0.9000`
- Runtime: R 4.6.1 UCRT, `x86_64-w64-mingw32`, Windows 11
- Compiler: GCC 14.3.0, C++17
- Candidate installation used `-fopenmp`; no thread or RNG policy changed
- Fixed fixture: `ape::rtree(seed = 20260904 + n)`, postorder
- Fixed trait: normal vector with seed `97531 + n`
- Formal timing: warmup, alternating before/after order, 20 pairs for
  `n < 10000`, 10 pairs for `n >= 10000`, median and IQR
- Child startup and serialization were outside the timer

The installed before and after libraries were loaded in separate child R
sessions. The formal run was serialized; no competing R benchmark was run.

## Production Change

The public `prepare_tree()` wrapper is unchanged. A private
`.prepare_tree_core()` now accepts package-owned inspection and
canonicalization evidence created earlier in the same raw public call.

The raw route is now:

```text
raw phylo
  -> representation-safe canonicalization once
  -> complete structural/method inspection once
  -> private context compilation once
  -> internal validated capability
  -> common prepared analysis route
```

An external prepared call remains:

```text
prepared context supplied by user
  -> one complete integrity fingerprint
  -> internal validated capability
  -> common prepared analysis route
```

The capability is private and cannot be requested through a public argument.
Each independent public prepared call still revalidates its context. Protected
mutations therefore remain detectable on the next call.

## Helper Call Audit

Counts are identical at n = 500 and n = 5000.

| Helper | Raw before | Raw after | Prepared before | Prepared after |
|---|---:|---:|---:|---:|
| `.inspect_tree_core()` | 3 | 1 | 0 | 0 |
| `.safe_canonicalize_core()` | 2 | 1 | 0 | 0 |
| `.canonical_tree_signature()` | 4 | 2 | 0 | 0 |
| `.tree_fingerprint()` | 2 | 1 | 1 | 1 |
| `.prepare_tree_subset()` | 1 | 1 | 0 | 0 |
| `.prepared_tree_subset()` | 1 | 1 | 1 | 1 |

Every remaining scan is required:

- The single inspection establishes structural validity and readiness for K,
  lambda, D, and Delta.
- The single canonicalization is the existing representation-safe raw-tree
  normalization.
- Its two canonical signatures are the unchanged before/after equivalence
  proof; the signature definition was not modified.
- The single fingerprint protects the newly compiled context and its retained
  subtree cache.
- The subset preparation seeds the full structural cache once.
- The retained-subtree lookup enters the unchanged common analysis path once.
- Downstream `.validate_prepared_context()` entries see the private validated
  capability and perform only shallow capability checks, not new fingerprints
  or tree scans.

## Performance

Times are formal medians in milliseconds. IQR values are in parentheses.

| Tips | Raw before | Raw after | Speedup | Prepare after | Prepared after | Raw / (prepare + prepared) | Duplicate residual / raw |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 500 | 186 (2.5) | 86 (1.75) | 2.16x | 85.75 (4.13) | 7 (0.5) | 0.927 | -7.85% |
| 1,000 | 481 (3.25) | 216 (5.25) | 2.23x | 208 (4.25) | 13 (1) | 0.977 | -2.31% |
| 2,000 | 1,412 (18) | 593 (11) | 2.38x | 590 (15.5) | 22 (2) | 0.969 | -3.20% |
| 5,000 | 7,160 (1,191.25) | 2,827.5 (278.75) | 2.53x | 2,797.5 (41.25) | 45 (<0.01) | 0.995 | -0.53% |
| 10,000 | 25,535 (210) | 10,015 (115) | 2.55x | 9,885 (110) | 95 (10) | 1.004 | 0.35% |
| 20,000 | 103,655 (6,617.5) | 37,260 (2,487.5) | 2.78x | 37,230 (935) | 180 (10) | 0.996 | -0.40% |

Negative residuals occur because the equation combines independent route
medians at the Windows timer resolution. They are small and indicate no
measurable duplicate residual, not negative computation.

All mandatory performance gates pass:

- n >= 5000 raw speedup is at least 1.7x.
- n = 20000 raw speedup exceeds the preferred 2.0x target.
- Absolute positive duplicate residual is below 15% of raw time.
- `raw / (prepare + prepared)` is below 1.20 at every n.
- K statistics are exactly equal between raw and prepared routes in every
  formal observation (`max absolute difference = 0`).

## Prepared-Path Stability

| Tips | Prepared before (ms) | Prepared after (ms) | Before / after |
|---:|---:|---:|---:|
| 500 | 7 | 7 | 1.000 |
| 1,000 | 12 | 13 | 0.923 |
| 2,000 | 22 | 22 | 1.000 |
| 5,000 | 45 | 45 | 1.000 |
| 10,000 | 90 | 95 | 0.947 |
| 20,000 | 190 | 180 | 1.056 |

There is no systematic prepared-path regression. Small differences are within
timer resolution and run variability; the prepared integrity fingerprint
remains exactly one call.

## Correctness and Failure Gates

Targeted tests confirm:

- raw/prepared parity for K, lambda, D, and Delta;
- tree and trait input immutability;
- root sensitivity and method-specific polytomy behavior;
- unary-node rejection;
- negative and zero-terminal branch failure behavior;
- malformed edge and invalid `Nnode` rejection;
- NA masks, species matching, duplicate/missing names, and retained-subtree
  behavior through the existing full suite;
- prepared mutation detection and one-fingerprint public boundary;
- warning classes and failure/status semantics through the existing suite.

D depends on random and Brownian null expectations even when significance
testing is disabled. Its raw/prepared parity gate therefore uses identical,
explicit null-state matrices rather than comparing two different RNG streams.
The numerical tolerance remains `1e-12`; no warning is suppressed to conceal a
production failure.

## Verification

- Targeted Stage 2B1 gate: `39 PASS, 0 FAIL, 0 WARN, 0 SKIP`
- Full testthat: `6902 PASS, 0 FAIL, 0 WARN, 0 SKIP`
- `R CMD check --no-manual --timings`: `Status: OK`
- Compiler build and package load: PASS
- Public API additions or changes: none
- Scientific, numerical, tolerance, RNG, or thread-policy changes: none

The first check attempt inherited `LC_ALL=C.UTF-8`, which Windows R cannot
initialize, and stopped during DESCRIPTION metadata checking. The same tarball
was rerun with the Windows-supported `LC_ALL=C`; it completed with `Status: OK`.
This was an execution-environment correction, not a source change.

## Evidence Index

Formal evidence is stored under `benchmarks/stage2b1/results/`:

- `stage2b1_before_after_summary.csv`
- `stage2b1_raw_preparation_summary.csv`
- `stage2b1_raw_preparation_timings.csv`
- `stage2b1_raw_helper_call_counts.csv`
- `stage2b1_helper_call_contract.csv`
- `stage2b1_provenance.csv`
- `stage2b1_call_audit_provenance.csv`
- `stage2b1_00check.log`
- `stage2b1_testthat.Rout`

The accepted production boundary is frozen. No structural traversal reuse,
fingerprint optimization, numerical optimization, scheduler redesign, or
Stage 2B2 work was performed.
