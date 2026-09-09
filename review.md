# fastphylosig 0.2.0 Stage 2D1 Expert Review

## 1. Decision

```text
K_PERMUTATION_PROTOTYPE = FAIL
CANDIDATE_A = REJECT
CANDIDATE_B = REJECT
CANDIDATE_C = ALREADY_PRESENT / NOT_WORTH_COMPLEXITY
PRODUCTION_CANDIDATE_SET = NONE
PRODUCTION_INTEGRATION = NOT_RUN
PRODUCTION_CODE_CHANGED = NO
```

Candidates A and B are scientifically exact private prototypes, but neither
produced the required end-to-end saving. They must not be integrated merely
because they reduce allocation or copy counts. Candidate C has no remaining
target: production already reads trait values through permutation indices and
does not construct a complete permuted-trait vector.

## 2. Frozen Scope And Provenance

| Item | Evidence |
|---|---|
| Production source | `baea4cf3306b2e9c314d446782dfbe0b98263db6` |
| Package | fastphylosig `0.2.0.9000` |
| Runtime | R 4.6.1 UCRT, Windows x86-64, 16 logical cores |
| Compiler | GCC 14.3.0 through the installed Rtools toolchain |
| Formal timing | fixed fixtures, warmup, serialized cells, alternating order, 10 paired repeats |
| Prepared grid | `n=500/2000/5000/10000`; `nsim=199/999/9999`; 1 and 2 threads |
| Production files changed | none in `R/`, `src/`, `tests/`, `DESCRIPTION`, `NAMESPACE`, or public API |

The prototype lives only under `benchmarks/stage2d1/prototype/` and is loaded
with `Rcpp::sourceCpp()`. It is not registered in package bindings. The local
R startup emitted host `C.UTF-8` locale warnings; the package calls and audit
captures were warning-free.

## 3. Frozen RNG Contract

The default internal K permutation path uses R's RNG through `R::runif()` and
a custom Fisher-Yates shuffle. The identity first replicate consumes no random
draws. Every later replicate consumes exactly `n - 1` uniform draws in the
same order. Permutations are generated serially before OpenMP evaluation.

Stage 2D1 preserves the existing promise: identical seed, retained-species
order, arguments, chunk setting, and worker configuration replay exactly.
Existing tests additionally demonstrate current chunk and thread invariance;
this audit preserves that behavior but does not create a broader public RNG
promise across platforms, RNG kinds, or arbitrary future schedulers.

## 4. Production Hot Loop

For each internally generated replicate, production performs:

1. one `n`-integer permutation allocation;
2. one identity reset and `n - 1` Fisher-Yates steps;
3. one `n`-integer generation-to-chunk copy;
4. one `n`-integer chunk-to-worker copy;
5. direct trait access through permutation indices, with no full trait gather;
6. four numeric K work-vector allocations per `compute_one()` call;
7. one exceedance update per trait and optional null-result storage.

At `n=5000`, `nsim=999`, one trait, the measured audit counts were:

| Mode | Index allocations | Generation to chunk | Chunk to worker | Trait gather allocations | Peak workspace proxy |
|---|---:|---:|---:|---:|---:|
| Exact private oracle | 1007 | 4,995,000 | 4,995,000 | 0 | 2.76 MB |
| A | 9 | 4,995,000 | 4,995,000 | 0 | 2.76 MB |
| B | 1 | 4,995,000 | 0 | 0 | 2.74 MB |

All modes retained 4,995,000 identity-initialization elements, 4,989,002
shuffle swaps, 25,001,000 indexed trait reads, and 4,000 numeric workspace
allocations. Candidate A therefore removes repeated index allocation but not
the dominant K work. Candidate B also removes the second index copy, reducing
the peak proxy by only 20 KB (`0.72%`). Memory remains bounded by the chunk;
no `n_tip x nsim` permutation matrix is created.

## 5. Exact Scientific Gates

Final Candidate B validation produced 1,169 gate records and 321,700
replicate-level ordered-null comparisons. Results were:

| Gate | Result |
|---|---|
| Controlled-permutation null K, replicate by replicate | exact PASS |
| Ordered internal-RNG null K | 321,700/321,700 exact |
| Maximum absolute null-K difference | 0 |
| Inclusive-tail exceedance mismatches | 0 |
| Observed K, P, MCSE_P | exact PASS |
| Requested/successful/failed accounting | exact PASS |
| Same-seed, same-ncores replay | exact PASS |
| Chunk boundaries and current chunk policy | exact PASS |
| 1- and 2-thread configurations | exact PASS |
| NA masks and retained analysis sets | PASS |
| Large offset, constant, near-constant, two-tip | PASS |
| Malformed controlled permutations and failure capture | PASS |
| Trait/tree input immutability | PASS |

The full ordered-null CSV had SHA-256
`972f5e1bf7ee245aa3a20ce49da12fca3366de1ae06918a332944f076b8a7442`.
Its compact manifest is retained in the repository rather than committing a
35 MB replicate dump.

## 6. Timing Calibration

The generic prototype entry point adds a `PermView` dispatch that production
does not need. A production-versus-private-oracle calibration on four heavy
cells showed the private oracle was 7.1-9.7% slower. Counting and ordered-null
storage were disabled for timing, but this structural prototype overhead
remained.

Candidate percentages therefore use the exact private oracle compiled in the
same translation unit. This symmetric comparison isolates the allocation and
copy changes. The production calibration is retained separately and no
prototype timing is presented as an installed-package speedup.

## 7. Heavy Prepared Timing

Representative formal medians are seconds:

| Candidate | n | nsim | Threads | Oracle | Candidate | Speedup | Reduction |
|---|---:|---:|---:|---:|---:|---:|---:|
| A | 2,000 | 9,999 | 1 | 1.050 | 1.050 | 1.000x | 0.0% |
| A | 5,000 | 9,999 | 1 | 2.600 | 2.575 | 1.010x | 1.0% |
| A | 10,000 | 9,999 | 1 | 5.220 | 5.230 | 0.998x | -0.2% |
| A | 10,000 | 9,999 | 2 | 3.220 | 3.190 | 1.009x | 0.9% |
| B | 2,000 | 9,999 | 1 | 1.030 | 1.030 | 1.000x | 0.0% |
| B | 5,000 | 9,999 | 1 | 2.595 | 2.555 | 1.016x | 1.5% |
| B | 10,000 | 9,999 | 1 | 5.180 | 5.165 | 1.003x | 0.3% |
| B | 10,000 | 9,999 | 2 | 4.250 | 4.220 | 1.007x | 0.7% |

After excluding timer-resolution zero/NA cells, the largest finite B speedup
among heavy cells was 1.075x in a short `n=500`, `nsim=9999`, two-thread cell.
No heavy cell reached 1.20x and no two heavy cells reached 15% reduction. Both
A and B therefore fail the predefined performance gate by a wide margin.

## 8. Light, Raw, And Batch Guards

The unchanged raw public route completed every grid cell without captured
package warnings or errors. At `n=10000`, one-thread raw medians were 1.00,
1.42, and 5.75 seconds for `nsim=199`, `999`, and `9999` respectively. These
are guard timings, not prototype speedups.

The unchanged public `test=FALSE` matrix route passed for 1/8/32/100 traits;
medians were 0.03/0.05/0.04/0.05 seconds. The private `test=TRUE` batch guard
also had exact oracle parity for every trait count. At 32 and 100 traits,
Candidate B measured 1.17x and 1.00x. The 1- and 8-trait cells were 0-20 ms
and produced direction-changing timer-resolution noise, so they do not
support a confirmed regression or a performance claim.

Several short prepared `nsim=199` cells also crossed the nominal 5% slowdown
line at 10 ms timer resolution. The machine-readable light gate therefore
remains conservatively `FAIL`; it is not waived. This has no effect on the
decision because the independent heavy-workload gate already failed.

## 9. Candidate Decisions

| Candidate | Exactness | Engineering effect | Performance | Decision |
|---|---|---|---|---|
| A: reusable generation buffer | PASS | 1007 to 9 index allocations | no material heavy saving | REJECT |
| B: A plus direct chunk-row evaluation | PASS | removes 4,995,000 index copies in audit cell | about 0-1.6% on representative heavy cells | REJECT |
| C: fused trait gather/evaluation | not applicable | production already has direct indexed reads | no remaining gather hotspot | ALREADY_PRESENT / NOT_WORTH_COMPLEXITY |

## 10. Final Recommendation

The smallest evidence-supported production candidate set is **none**. The
index allocation and copy work is real, but it is not a meaningful fraction
of the present heavy K permutation wall time. Integrating A or B would add
maintenance surface without meeting the agreed benefit threshold.

Stop Stage 2D1 here. Do not modify production bindings, redesign RNG, loosen
scientific tolerances, or start D optimization as part of this stage.

## Evidence Index

- [RNG contract](benchmarks/stage2d1/K_RNG_CONTRACT_STAGE2D1.md)
- [source hot-loop audit](benchmarks/stage2d1/SOURCE_HOT_LOOP_AUDIT.md)
- [protocol](benchmarks/stage2d1/PROTOCOL.md)
- [machine-readable final decision](benchmarks/stage2d1/results/FINAL_DECISION.csv)
- [Candidate A summary](benchmarks/stage2d1/results/candidate-a/stage2d1_prepared_summary.csv)
- [Candidate B summary](benchmarks/stage2d1/results/candidate-b/stage2d1_prepared_summary.csv)
- [final correctness status](benchmarks/stage2d1/results/candidate-b/stage2d1_correctness_status.csv)
- [batch exactness and timing](benchmarks/stage2d1/results/candidate-b/stage2d1_batch_test_true_summary.csv)
- [operation counters](benchmarks/stage2d1/results/operation-audit/operation_counters.csv)
- [ordered-null manifest](benchmarks/stage2d1/results/ORDERED_NULLS_MANIFEST.csv)
