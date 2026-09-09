# fastphylosig 0.2.0 Stage 2C Expert Review

Date: 2026-09-09  
Scope: post-V2 wall-time, memory, persistence, and candidate qualification  
Production code changed during Stage 2C: **NO**

## 1. Decision

```text
PREPARATION_ENGINE_STAGE1 = ACCEPTED_FROZEN
RAW_ROUTE_SINGLE_PREPARATION_BOUNDARY = ACCEPTED_FROZEN
INSPECTION_INTEGER_ADJACENCY = ACCEPTED_FROZEN
CANDIDATE_2 = REJECTED
CANONICALIZATION_CONTRACT_V2 = ACCEPTED_FROZEN
PUBLIC_RELEASE_BLOCKER_0_2_0 = CLEARED
COMMON_PREPARATION_ENGINE_0_2_0 = FROZEN

STAGE2C = PASS
PRODUCTION_CODE_CHANGED = NO
NEXT_STAGE = K_PERMUTATION_NULL_ENGINE_DESIGN; D_NULL_ENGINE_ORDER_STATISTIC_DESIGN
```

Stage 2C found no correctness or status regression. It authorizes design work
for two statistic-specific candidates; it does not authorize implementation.

## 2. Provenance

| Item | Evidence |
|---|---|
| Production source | `0636452b931dcfe15ece771b78af3b114c51f200`; production files are unchanged through the Stage 2C commits |
| K/lambda/D run commits | `85f1da4189188b32c306a9b47a5404e77f9c97c2` |
| Delta/context run commit | `6b70b77b97640a0919f248832f57f42bd1b02ce8` |
| Package | fastphylosig `0.2.0.9000` |
| R/platform | R 4.6.1 UCRT, `x86_64-w64-mingw32` |
| OS | `sessionInfo()`: Windows 11 build 26200; `Sys.info()`: Windows 10 x64 |
| Compiler | formal runtime provenance unavailable; the installation console showed GCC 14.3, C++17, `-O2`, but that console log is not an attached artifact |
| OpenMP | formal provenance reported `not independently exposed` / `unavailable`; D ncores=2 was exercised, but no definitive runtime capability flag was recorded |
| BLAS/LAPACK | default matrix products; LAPACK 3.12.1; BLAS identity unavailable at runtime |
| Locale | `LC_ALL=C`, `LC_COLLATE=C`, system code page 65001 |
| Hardware | AMD64 Family 25 Model 80, 8 physical / 16 logical cores |
| Formal execution | serialized; fixed fixtures; warmup; alternating raw/prepared order |
| Evidence directory | `benchmarks/stage2c/results/2026-09-09-post-v2/` |

The estimator groups were run at different **audit-only** commits, but all used
the same installed production source. Speedups below are only calculated
within a method's own paired run; no timing from one commit is divided by a
timing from another commit.

The original `lambda-formal/stage2c_methods_summary.csv` is a superseded empty
artifact caused by an `NA` grouping-key bug in the audit harness. The
authoritative derived summary is
`lambda-formal/stage2c_methods_summary_rebuilt.csv` (60 grouped rows from all
375 raw timing rows). Raw timings were not altered.

## 3. Unified Budget

Times are medians in milliseconds. `Prepared` is compute time; one-shot
prepared end-to-end time is the sum of the independently measured preparation
and prepared-compute medians. It is a derived comparison, not a per-pair
identity with the raw median. Phase values marked `diag` are inclusive probes
and are not additive. Lower bounds come from paired simulation-count
differences.

| Method | Workload | Raw | Prepare | Prepared | Validation | Kernel/optimizer | Simulation | Startup | Packaging | Reference | Speedup | Context cost |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| K | test=FALSE, n=5000 | 445 | 440 | 65 | diag only | diag only | 0 | 0 | diag only | RESOURCE_LIMIT | NA | 13,434,896 B (12.81 MiB) at n=5000 |
| K | test=TRUE, n=10000, nsim=9999 | 7290 | 950 | 6460 | 110 diag | 98.2% combined null engine | 98.2% combined | 0 | 0 diag | not run | NA | 26,834,896 B (25.59 MiB) at n=10000 |
| lambda | random, n=5000 | 410 | 380 | 60 | non-authoritative diag | non-authoritative diag | NA | 0 | non-authoritative diag | RESOURCE_LIMIT | NA | 13,434,896 B (12.81 MiB) at n=5000 |
| D | n=2000, nsim=9999, 1 thread | 4835 | 140 | 4770 | diag only | included | >=4580 | 0 | diag only | not run | NA | not measured at n=2000 |
| Delta | n=500, nsim=199, serial | 6710 | 30 | 6690 | diag only | MCMC/ACE | >=5970 | 0 | diag only | not run | NA | 1,375,896 B (1.31 MiB) at n=500 |
| Delta | n=500, nsim=199, 2 workers | 3890 | 30 | 3850 | diag only | MCMC/ACE | dominant | about 200 | diag only | not run | 1.73x vs serial prepared | 1,375,896 B (1.31 MiB) at n=500 |
| context | balanced, n=20000 | NA | 1750 | NA | read+first validation 590 | NA | NA | NA | saveRDS 230 | NA | NA | object 53,654,896 B (51.17 MiB); RDS 39,878,633 B (38.03 MiB) |

The fused K phase contains serial RNG, permutation generation, trait gather,
null-K evaluation, and exceedance reduction. The individual RNG, null-K, and
copy fractions were not separately identified without intrusive tracing;
their authoritative combined heavy-cell share is 98.2%.

## 4. K Results

### Single-Trait Crossover

| n | phytools | fast_k raw | prepare | prepared compute | prepared end-to-end | fast_signal |
|---:|---:|---:|---:|---:|---:|---:|
| 50 | 0 | 10 | 0 | 0 | 0 | 10 |
| 100 | 0 | 15 | 10 | 0 | 10 | 10 |
| 300 | 5 | 25 | 30 | 10 | 40 | 25 |
| 500 | RESOURCE_LIMIT | 35 | 40 | 10 | 50 | 40 |
| 1000 | RESOURCE_LIMIT | 80 | 75 | 20 | 95 | 75 |
| 2000 | RESOURCE_LIMIT | 160 | 145 | 20 | 165 | 160 |
| 5000 | RESOURCE_LIMIT | 445 | 440 | 65 | 505 | 470 |
| 10000 | RESOURCE_LIMIT | 860 | 950 | 100 | 1050 | 1050 |
| 20000 | RESOURCE_LIMIT | 2240 | 2290 | 240 | 2530 | 2120 |

No crossover was observed in the measurable reference range: phytools was
still faster at n=300. Reference cells at n>=500 were preflighted as
`RESOURCE_LIMIT` because of dense repeated cubic work, so the large-n
crossover is **not established** and no speedup is extrapolated.

### Batch K

At n=5000 and 100 traits, raw/prepared matrix calls took 590/155 ms, while
raw/prepared column loops took 45,810/5,295 ms. Prepared matrix batching was
34.2x faster than the prepared column loop. Equivalent batch advantages were
retained at n=500 and n=2000. Every phytools batch reference cell was
`RESOURCE_LIMIT`; no reference speedup is claimed.

### Permutation K

| n | nsim | raw | prepared | raw/prepared |
|---:|---:|---:|---:|---:|
| 500 | 199 / 999 / 9999 | 45 / 60 / 280 | 20 / 30 / 250 | 2.25x / 2.00x / 1.12x |
| 2000 | 199 / 999 / 9999 | 175 / 250 / 1155 | 50 / 115 / 1015 | 3.50x / 2.17x / 1.14x |
| 5000 | 199 / 999 / 9999 | 485 / 710 / 3455 | 125 / 360 / 3010 | 3.88x / 1.97x / 1.15x |
| 10000 | 199 / 999 / 9999 | 1090 / 1520 / 7290 | 220 / 750 / 6460 | 4.95x / 2.03x / 1.13x |

For heavy prepared K, the fused permutation C++ phase occupied 97.2-100% of
the probe; the n=10000/nsim=9999 share was 98.2%. The all-grid median (40.8%)
in the harness candidate table mixes light and heavy cells and is not the
heavy-K gate specified for Stage 2C. That table remains a light-to-heavy grid
diagnostic; `stage2c_final_candidate_decisions.csv` is the authoritative
machine-readable Stage 2C authorization.

### Prepared V2 Diagnostic

At n=20000, the direct prepared call was 280 ms. The diagnostic probe was
230 ms and attributed 200 ms inclusively to protected context validation/V2
evidence, 60 ms to fingerprinting, 220 ms to `.prepare_analysis()`, and 10 ms
to the K kernel. These nested durations overlap and the probe/direct ratio was
0.82, so they localize the V2 cost but do not form an additive budget. The
evidence supports exact integrity validation as the source of the previously
observed roughly 40 ms regression, but does not support weakening it or an
exact-equivalent optimization candidate yet.

## 5. Lambda Results

Raw/prepared medians in milliseconds:

| Shape | n=100 | n=500 | n=2000 | n=5000 | n=10000 |
|---|---:|---:|---:|---:|---:|
| balanced | 10/0 | 40/10 | 160/30 | 420/70 | 880/140 |
| random | 20/0 | 40/10 | 150/30 | 410/60 | 880/130 |
| pectinate | 10/0 | 40/10 | 170/30 | 440/80 | 1110/140 |

The pectinate n=10000 path completed normally; the former recursive-depth bug
did not recur. The one authoritative phase probe, pectinate raw n=5000,
assigned 400/440 ms (90.9%) to preparation. Unique likelihood-evaluation
counts were not exposed by the production optimizer, and other phase probes
exceeded the 5% tracing limit. Therefore no lambda-specific candidate meets
the evidence gate. References completed only at n=100; larger cells are
`RESOURCE_LIMIT`.

## 6. D Results

At n=2000, prepared one-thread times were 190, 540, and 4770 ms for nsim 199,
999, and 9999. The added 9800 simulations account for at least 4580/4770 =
96.0% of the heavy workload. Two threads took 310, 1030, and 9600 ms, so the
current two-thread route is slower on this fixture. Startup is not a PSOCK
cost: `ncores` controls the OpenMP kernel. All 483 formal timing rows were OK,
with no warning, error, or resource-limit status.

## 7. Delta Results

At n=500, prepared serial time rose from 60 ms (`test=FALSE`) to 720 ms at
nsim=19 and 6690 ms at nsim=199. The nsim increment accounts for at least
89.2% of the nsim=199 wall time. Two workers reduced the prepared nsim=199
median to 3850 ms (1.74x) but did not help nsim=19 because approximately
200 ms of PSOCK startup dominates the small job.

The authoritative prepared serial n=100/nsim=19 probe reported inclusive MCMC,
observed-ACE, and permutation-ACE shares of 61.5%, 33.3%, and 89.7%; these
overlap and must not be summed. Ten raw/prepared warning events at n=500,
nsim=199, serial were paired instances of `singular convergence (7)` from
transition-rate optimization. Estimates and warnings matched between raw and
prepared paths, so this is not a Stage 2C warning/status regression. Candidate
4A remains rejected. The n=100/nsim=19 two-worker **trace probes** returned
errors and were marked non-authoritative; the corresponding formal timing
calls succeeded and remain in the 93/93 OK total.

## 8. Context Memory And Persistence

| n | context bytes | fingerprint | snapshot | cache payload | RDS bytes | save ms | read+validate ms |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 500 | 1,375,896 | 20,674 | 20,624 | 2,413,832 | 1,017,237 | 20 | 10 |
| 1000 | 2,714,896 | 41,175 | 41,125 | 4,788,960 | 2,005,522 | 10 | 20 |
| 5000 | 13,434,896 | 209,175 | 209,125 | 23,798,032 | 9,955,664 | 60 | 180 |
| 10000 | 26,834,896 | 419,176 | 419,126 | 47,559,360 | 19,893,329 | 110 | 310 |
| 20000 | 53,654,896 | 849,176 | 849,126 | 95,102,016 | 39,878,633 | 230 | 590 |

All fresh-process reads and first validations passed. Growth was approximately
linear through n=20000. The 20,000-tip context occupied 53,654,896 bytes
(51.17 MiB), so context memory is not yet a demonstrated practical limit at
that size. `cache payload` is the sum of measured package-owned objects, not
peak RSS, and cannot be extrapolated to larger trees or memory-constrained
systems.

## 9. Correctness Guard

K raw/prepared parity, controlled permutations, edge-order invariance, input
immutability, lambda parity, D controlled-null accounting, Delta fixed-seed
parity, and raw/prepared warning/status parity all passed. K had 7/7 guard
checks; lambda, D, and Delta each passed their method guard. Formal status was
K PASS, lambda 285 OK plus 90 explicit reference `RESOURCE_LIMIT` rows, D
483/483 OK, Delta 93/93 OK, and context 5/5 fresh-process validations OK.

## 10. Candidate Ranking

| Candidate | Hotspot share | Ideal ceiling | Practical saving estimate | Correctness risk | Engineering risk | Decision |
|---|---:|---:|---:|---|---|---|
| K permutation/null engine: eliminate per-replicate permutation allocation and two full-index copies while preserving serial Fisher-Yates draws and static result order | 98.2% heavy prepared cell | 98.2% | 10-30% | Medium: RNG/thread invariance and exact exceedance counts | Medium | GO for candidate design; HIGH_PRIORITY |
| D random/Brownian null engine: replace each full Brownian sort with exact two-order-statistic selection and reuse bounded buffers | >=96.0% heavy prepared cell | >=96.0% | 15-35% | Medium: exact type-7 threshold, ties, RNG draw order | Medium | GO for candidate design |

The K design must preserve the exact permutation stream, observed-first rule,
thread invariance, tail comparison, P/MCSE, and returned simulation order. The
D design must return the same two type-7 order statistics and threshold, keep
all R RNG calls in the same order, and preserve random/Brownian accounting.

Delta permutation ACE ranks third: its heavy share qualifies, but only three
formal repeats were used, the path produced paired optimizer diagnostics, and
an exact batching design must preserve PSOCK RNG/task order. With the two-item
cap, it is NO-GO for 0.2.0. Prepared validation and lambda are also NO-GO:
neither has both authoritative phase evidence and a proven exact-equivalent
engineering design.

## 11. Required Answers

1. **K crossover:** not observed through n=300; n>=500 reference cells hit the declared resource limit, so the crossover is unresolved.
2. **K test=FALSE:** no statistic-specific numerical optimization is justified; preparation/integrity work dominates.
3. **K test=TRUE hotspot:** the fused permutation/null C++ engine, 98.2% in the representative heavy prepared cell.
4. **Prepared V2 regression:** exact protected-state validation and V2 evidence/fingerprint comparison; diagnostic phases overlap, so no additive 40 ms claim is made.
5. **Lambda bottleneck:** raw-tree preparation (90.9% in the authoritative pectinate n=5000 probe), not a verified likelihood/optimizer hotspot.
6. **D bottleneck:** random/Brownian null generation and evaluation; at least 96.0% in the heavy one-thread workload.
7. **Delta bottleneck:** MCMC plus permutation ACE; at least 89.2% by paired nsim differences. PSOCK startup matters mainly for small jobs.
8. **Context memory:** not a demonstrated limit through n=20000; peak RSS was not measured.
9. **Final candidates:** K permutation/null engine and D random/Brownian null engine.
10. **Candidate risks and ceilings:** recorded in the ranking table; both are design-only GO decisions.
11. **Stage status:** `STAGE2C = PASS`.
12. **Next stage:** design and equality-gate the two named candidates, one at a time; do not implement automatically from this report.
