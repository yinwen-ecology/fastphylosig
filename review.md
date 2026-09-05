# fastphylosig 0.2.0 Stage 2B2B Candidate 1 Expert Review

## Decision

`PREPARATION_ENGINE_STAGE1 = ACCEPTED_FROZEN`

`RAW_ROUTE_SINGLE_PREPARATION_BOUNDARY = ACCEPTED_FROZEN`

`STAGE2B2A_COMPLEXITY_ROOT_CAUSE = CONFIRMED`

`INSPECTION_INTEGER_ADJACENCY = ACCEPTED_FROZEN`

`CANDIDATE_2 = NOT_STARTED`

Candidate 1 passed exact semantic equality, regression, package-check, and
performance gates. It is accepted and frozen. This stage stops here.

## Scope And Provenance

- Old-oracle baseline: `4dc4fa3b4eb6a0815789bdbd5276ec3063dcbda4`.
- Candidate source commit: `1685d8f58dea6dfe0a3e3d519d58581b3bbae5b9`.
- Package version: `0.2.0.9000`.
- Runtime: R 4.6.1 UCRT, `x86_64-w64-mingw32`, Windows 11 build 26200.
- Compiler evidence: GCC/G++ 14.3.0; package compilation used C++17.
- Production scope: `R/check_tree.R` only.
- Test-only oracle: the complete baseline `.inspect_tree_core()` copied from
  the baseline commit into a private test fixture.
- No estimator, numerical formula, tolerance, fingerprint, canonicalization,
  canonical-signature, RNG, thread-policy, result-structure, or public-API
  change was made.

## Implementation Audit

The two structural traversals inside `.inspect_tree_core()` retain the same
validation and visit semantics, but use different private workspaces:

| Concern | Frozen old implementation | Candidate 1 |
|---|---|---|
| Child lookup | named split lists and character keys | integer `child_offsets`, `children`, and edge-row arrays |
| Child order | source edge-row order | source edge-row order |
| Stack | `stack[-length(stack)]` and `c(stack, kid)` | preallocated integer stack plus `top` cursor |
| Connectivity traversal | independent named adjacency | shared integer adjacency |
| Root-distance traversal | second named adjacency | same shared integer adjacency |
| Representation | unchanged | unchanged |

The integer adjacency is filled in original edge-row order. The LIFO traversal
therefore visits children in exactly the same order as the old implementation.
Arbitrary valid internal node numbers, shuffled edge rows, and supported
polytomies remain valid. Cycle detection clears the effective stack and the
`seen` vector exactly as before. The input tree is never modified.

## Exact-Equality Gate

The old and new complete inspection objects were compared with
`expect_identical()`, not with a numerical tolerance or selected booleans.
The gate also compared status, method readiness, root, degree and branch
diagnostics, issue classes and ordering, failure reasons, warning/error classes
and messages, and serialized input trees before and after each call.

Valid fixtures:

- balanced, seeded random, and pectinate trees;
- polytomy, shuffled edge rows, safe internal renumbering, and alternate valid
  root numbering;
- two-tip and 1,000-tip large trees.

Failure fixtures:

- cycle, disconnected tree, unary node, malformed edge matrix;
- invalid `Nnode`, invalid root representation;
- negative branch length and zero terminal branch.

All complete-object and input-immutability comparisons passed for the full
signal set and selected signal subsets.

## Formal Performance Protocol

The formal benchmark used fixed positive-length balanced, random, and
pectinate fixtures at 500, 1,000, 2,000, 5,000, 10,000, and 20,000 tips. One
serialized R process performed warmup, old/new parity checks, fresh fixture
clones, and alternating old/new order. Namespace binding changes, cloning, and
garbage collection were outside the timer.

- 10 paired repeats per cell through 5,000 tips.
- 5 paired repeats per cell at 10,000 and 20,000 tips.
- 36 workload cells and 300 paired timing records.
- Every record: status `ok`, zero warnings, zero errors, zero censoring.
- Results: median and IQR; exact values are retained in
  `benchmarks/stage2b2b1/results/formal/`.

### Direct Inspection

Times are median seconds with IQR in brackets.

| Shape | Tips | Old | New | Speedup |
|---|---:|---:|---:|---:|
| Balanced | 5,000 | 1.255 [0.020] | 0.020 [0.0175] | 62.75x |
| Random | 5,000 | 1.330 [0.0375] | 0.025 [0.010] | 53.20x |
| Pectinate | 5,000 | 1.515 [0.0575] | 0.030 [~0] | 50.50x |
| Balanced | 20,000 | 21.910 [0.350] | 0.100 [0.020] | 219.10x |
| Random | 20,000 | 21.230 [1.420] | 0.110 [0.010] | 193.00x |
| Pectinate | 20,000 | 27.770 [1.330] | 0.170 [0.090] | 163.35x |

The descriptive 10,000-to-20,000 empirical exponents changed from 1.973 to
1.000 for balanced trees, 1.962 to 1.138 for random trees, and 2.077 to 1.766
for pectinate trees. These are empirical diagnostics, not complexity proofs.

### prepare_tree() End To End

Times are median seconds with IQR in brackets.

| Shape | Tips | Old | New | Speedup | Improvement |
|---|---:|---:|---:|---:|---:|
| Balanced | 5,000 | 2.645 [0.0425] | 1.405 [0.0275] | 1.88x | 46.88% |
| Random | 5,000 | 2.800 [0.220] | 1.455 [0.065] | 1.92x | 48.04% |
| Pectinate | 5,000 | 4.740 [0.065] | 3.235 [0.0425] | 1.47x | 31.75% |
| Balanced | 20,000 | 38.920 [1.800] | 17.750 [0.590] | 2.19x | 54.39% |
| Random | 20,000 | 39.160 [2.710] | 18.000 [1.420] | 2.18x | 54.03% |
| Pectinate | 20,000 | 83.910 [2.060] | 55.250 [3.130] | 1.52x | 34.16% |

The remaining pectinate scaling is expected because Candidate 1 did not touch
canonical descendant-key construction. Candidate 2 was not started.

## Small-Tree Resolution Gate

The formal pectinate 500-tip direct cell initially showed 15 ms versus 20 ms,
but both values were at the approximately 10 ms timer resolution and the IQR
was 20 ms. A pre-specified supplemental run timed 50 calls per block over 24
alternating paired blocks, with exact equality outside every timed block.

| Shape | Old, ms/call | New, ms/call | Speedup | Slowdown gate |
|---|---:|---:|---:|---:|
| Balanced | 18.3 [2.05] | 3.3 [0.40] | 5.55x | PASS |
| Random | 18.0 [0.80] | 3.4 [0.20] | 5.29x | PASS |
| Pectinate | 20.6 [0.40] | 3.4 [0.40] | 6.06x | PASS |

No representative workload has a confirmed slowdown greater than 5%.

## Operation Counts

For a binary 20,000-tip tree, both structural traversals together formerly
performed 79,998 character lookups, 79,998 character-key materializations,
79,998 stack pops, and 79,996 stack appends. The old copying-stack workload
depended strongly on tree shape:

| Shape | Old copied stack slots | New stack resize copied slots |
|---|---:|---:|
| Balanced | 1,033,732 | 0 |
| Random | 1,463,388 | 0 |
| Pectinate | 1,599,840,004 | 0 |

The candidate instead performs fixed integer adjacency fills, offset/cursor
reads, and preallocated stack writes. It does not remove either traversal or
any structural validation.

## Regression And Check Gates

- Full `testthat`: **7,412 PASS, 0 FAIL, 0 WARN, 0 SKIP**.
- `R CMD check --no-manual --timings`: **Status: OK**.
- Package check: 0 ERROR, 0 WARNING, 0 NOTE.
- Exact old/new inspection equality: PASS.
- Input immutability: PASS.
- Representative performance slowdown gate: PASS.

Authoritative outputs are retained under
`benchmarks/stage2b2b1/results/`, including the full paired timings,
summaries, operation counts, provenance, `testthat.Rout`, `00install.out`, and
`00check.log`.

## Final Gate

The acceptance requirement was met by both independent routes:

- balanced/random 20,000-tip direct inspection speedup is far above 2x;
- 20,000-tip `prepare_tree()` improvement is 54.39% and 54.03%, with 34.16%
  improvement for pectinate trees;
- no confirmed representative slowdown exceeds 5%.

`INSPECTION_INTEGER_ADJACENCY = ACCEPTED_FROZEN`

Stage 2B2B Candidate 1 ends here. Candidate 2 remains unimplemented.
