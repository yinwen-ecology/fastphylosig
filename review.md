# fastphylosig 0.2.0 Stage 2D2B Expert Review

## 1. Final Decision

```text
TRUE_BLOCKED_IMPLEMENTATION = VALID
CONTROLLED_PERMUTATION_EXACTNESS = PASS
IDENTITY_FIRST_ACCOUNTING = PASS
RNG_REPLAY = PASS
EVALUATION_ONLY_GATE = FAIL
FULL_PIPELINE = NOT_RUN_BY_GATE
LIGHT_BATCH_PARALLEL_PERFORMANCE_GUARDS = NOT_RUN_BY_GATE

K_BLOCKED_EVALUATOR_PROTOTYPE = FAIL
BLOCKED_EVALUATION = REJECTED_NOT_ENOUGH_EVALUATOR_GAIN
K_OPTIMIZATION_0_2_0 = CLOSED_FINAL
NEXT_STAGE = D_NULL_ENGINE
```

The remediated Candidate B is a genuine blocked evaluator and is scientifically
exact over the bounded audit. It nevertheless fails the predefined performance
authorization gate: only one representative heavy cell reached 1.35x, while
at least two were required. It must not enter production.

## 2. Why The Previous Candidate Was Invalid

The Stage 2D2 Candidate B did not test blocked K evaluation. Its attempted
block passes used different arithmetic, discarded their result, and then
called frozen `compute_one()` for every returned null K. Its identity-first
P/MCSE accounting also differed from production. That result proved only
`IMPLEMENTATION_INVALID`; it could not establish a numerical no-go.

Stage 2D2B removes those defects. The candidate's observed and null K values
come directly from its own `compute_block()`. No candidate result is copied
from the oracle, production entry point, or frozen `compute_one()`.

## 3. True Block Architecture

The prototype uses node-major, replicate-minor bounded storage:

```text
message[N x B]
state[N x B]
baseline[B]
delta[B]
permutations[n x B]
```

For each traversal node, the evaluator advances the independent replicates in
the block. Within each replicate, it preserves the frozen postorder/preorder,
child order, long-double accumulation order, division form, and final K
operation order. Supported block sizes are 2, 4, 8, and 16. The candidate does
not allocate `N x nsim` numerical state.

Permutation generation remains serial descending Fisher-Yates. The identity
first replicate consumes no RNG draws. Batching changes only evaluator layout,
not draw order or RNG partition.

## 4. Provenance And Hard Counters

| Item | Evidence |
|---|---|
| Frozen production source | `57d7cccb3889fed603b8454c6b8963d28bc33705` |
| Package | fastphylosig `0.2.0.9000` |
| Prototype SHA-256 | `9e6d6f1a507bd4c03b0fe15058e749fede5cd555910fa5a07041e633e364a114` |
| Runtime | R 4.6.1 UCRT, Windows x86-64 |
| Compiler | GCC 14.3.0 through the installed Rtools toolchain |
| Oracle/candidate build | same `Rcpp::sourceCpp()` translation unit |
| Production files changed | none in `R/`, `src/`, tests, DESCRIPTION, NAMESPACE, or public API |

Every successful candidate call was required to report:

```text
candidate_old_compute_one_calls = 0
oracle_call_count = 0
candidate_compute_count > 0
```

All 388 formal counter records passed. The timing runner independently rejected
missing or non-zero counters before accepting a pair.

## 5. Controlled Exactness

Formal coverage included balanced, random, and pectinate trees; `n=50`, 500,
2,000, and 5,000; blocks 2/4/8/16; and identity, single-swap, reverse, cyclic,
fixed-random, and adversarial permutation layouts.

| Gate | Result |
|---|---|
| Controlled cases | 288/288 PASS |
| Controlled failures | 0 |
| Ordered replicate comparisons | 8,283 |
| Ordered null K | bitwise exact |
| Maximum absolute null-K difference | 0 |
| Exceedance, P, MCSE | exact PASS |
| Candidate result provenance counters | PASS |
| Input immutability | PASS |
| Failure-semantics cases | 5/5 PASS |
| Batch correctness guards | 1/8/32/100 traits PASS |

The detailed ordered-result file has SHA-256
`0e705dc8eefbe026cf8b0fe52d10a77dc299878da958dda5d04888b464c40b6e`.
The compact correctness status has SHA-256
`5105ebe383da0c8b73b7f7a658ab104ebf11526a1bc1949e61a8c7918f6f880c`.

## 6. Identity-First And RNG Integration

The candidate reproduces the frozen accounting contract:

- replicate 1 is identity when `include_observed=TRUE` and permutations are
  internally generated;
- replicate 1 consumes zero RNG draws;
- later replicates use the same descending Fisher-Yates draws and order;
- requested, successful, failed, returned order, inclusive exceedance, P, and
  MCSE match the oracle;
- constant/non-finite observed K produces `successful=0` and `failed=nsim`;
- 96 same-seed RNG replay cases across current one- and two-thread test
  configurations passed exactly.

No tolerance was widened, and no tail or Monte Carlo definition changed.

## 7. Evaluation-Only Performance

Evaluation timing used pre-generated controlled permutations, excluding RNG
generation. Oracle and candidate were compiled together. Every cell used a
fixed fixture, warmup, serialized execution, alternating call order, and 10
paired repeats. The formal grid contained 72 cells across three tree shapes,
`n=2000/5000/10000`, `nsim=999/9999`, and four block sizes.

Across the representative heavy cells (`n>=5000`, `nsim=9999`), median
performance aggregated by block was:

| Block | Oracle median | Candidate median | Median speedup | Median reduction |
|---:|---:|---:|---:|---:|
| 2 | 3.648 s | 3.815 s | 1.003x | 0.3% |
| 4 | 3.718 s | 3.435 s | 1.179x | 15.2% |
| 8 | 3.680 s | 3.148 s | 1.280x | 21.9% |
| 16 | 3.630 s | 3.208 s | 1.272x | 21.3% |

Blocks 8 and 16 were within 3%, so the smaller `B=8` was selected. It was also
slightly faster in the aggregate and uses half the bounded workspace.

Selected-block heavy results were:

| Shape | n | Oracle | Candidate B=8 | Speedup | Reduction |
|---|---:|---:|---:|---:|---:|
| balanced | 5,000 | 2.260 s | 2.050 s | 1.102x | 9.3% |
| balanced | 10,000 | 4.470 s | 4.310 s | 1.037x | 3.6% |
| random | 5,000 | 2.890 s | 2.265 s | 1.276x | 21.6% |
| random | 10,000 | 6.080 s | 4.735 s | 1.284x | 22.1% |
| pectinate | 5,000 | 2.660 s | 1.965 s | 1.354x | 26.1% |
| pectinate | 10,000 | 5.295 s | 4.030 s | 1.314x | 23.9% |

Only pectinate `n=5000` reached the required evaluator speedup of 1.35x.
The required count was at least two representative heavy cells. Therefore:

```text
evaluation_heavy_pass_cells = 1
EVALUATION_ONLY_GATE = FAIL
```

## 8. Full Pipeline And Regression Gates

The protocol authorizes full-pipeline timing only after the evaluation gate
passes. It did not pass, so the full runner stopped before candidate
compilation or timing and recorded:

```text
FULL_PIPELINE = NOT_RUN_BY_GATE
reason = evaluation gate is not PASS (status=FAIL)
```

The light, `test=FALSE`, performance batch, and 1/2/4/8-thread scaling guards
depend on a passing full-pipeline result. Their runner likewise stopped before
timing and recorded:

```text
LIGHT_BATCH_PARALLEL_PERFORMANCE_GUARDS = NOT_RUN_BY_GATE
reason = full-pipeline gate is not PASS (status=NOT_RUN)
```

This is not a PASS claim for unrun work. Exact batch correctness and one-/two-
thread RNG replay were covered by the scientific gate; performance guards and
four-/eight-thread interaction were not run because the candidate had already
failed authorization.

## 9. Memory

At `n=10000`, the candidate reported the following upper working payloads:

| Block | Permutation workspace | Numeric workspace | Total candidate workspace |
|---:|---:|---:|---:|
| 2 | 0.08 MB | 0.64 MB | 0.72 MB |
| 4 | 0.16 MB | 1.28 MB | 1.44 MB |
| 8 | 0.32 MB | 2.56 MB | 2.88 MB |
| 16 | 0.64 MB | 5.12 MB | 5.76 MB |

The memory gate passed. Workspace grows linearly with `N * B`, remains bounded
independently of `nsim`, and completed the `n=10000`, `nsim=9999` workload.
The externally supplied controlled-permutation matrix is benchmark input and
is excluded from candidate working-memory claims.

## 10. Expert Conclusion

Stage 2D2B corrects the scientific and provenance defects of the first blocked
attempt. It demonstrates that block-interleaved traversal can reduce evaluator
time for some large random and pectinate trees. The gain is not sufficiently
general: balanced trees improve only 3.6-9.3%, and only one heavy cell meets
the explicit 1.35x authorization threshold.

The candidate is therefore rejected for insufficient evaluator gain, not for
incorrectness. No Stage 2D3 integration is authorized. This is the final K
optimization stop for 0.2.0; the next permitted investigation is the D null
engine.

## Evidence Index

- [exactness protocol](benchmarks/stage2d2b/EXACTNESS_PROTOCOL.md)
- [performance protocol](benchmarks/stage2d2b/PERFORMANCE_PROTOCOL.md)
- [correctness status](benchmarks/stage2d2b/results/correctness-formal/stage2d2b_correctness_status.csv)
- [correctness summary](benchmarks/stage2d2b/results/correctness-formal/stage2d2b_correctness_summary.csv)
- [counter evidence](benchmarks/stage2d2b/results/correctness-formal/stage2d2b_correctness_counters.csv)
- [evaluation status](benchmarks/stage2d2b/results/evaluation-formal/stage2d2b_benchmark_status.csv)
- [evaluation summary](benchmarks/stage2d2b/results/evaluation-formal/stage2d2b_benchmark_summary.csv)
- [memory evidence](benchmarks/stage2d2b/results/evaluation-formal/stage2d2b_benchmark_memory.csv)
- [full-pipeline gate](benchmarks/stage2d2b/results/full-gated/stage2d2b_benchmark_status.csv)
- [guard gate](benchmarks/stage2d2b/results/guards-gated/stage2d2b_benchmark_status.csv)
- [machine-readable final decision](benchmarks/stage2d2b/results/FINAL_DECISION.csv)
