# fastphylosig 0.2.0 Stage 2D2 Expert Review

## 1. Final Decision

```text
STAGE2D2_PHASE_DECOMPOSITION = PASS
CANDIDATE_A = REJECT
CANDIDATE_B = REJECT
CANDIDATE_C = REJECT
PRODUCTION_CANDIDATE_SET = NONE
PRODUCTION_CODE_CHANGED = NO

K_NUMERICAL_NULL_PROTOTYPE = FAIL
K_OPTIMIZATION_0_2_0 = CLOSED_FINAL
NEXT_STAGE = D_NULL_ENGINE
```

The evaluation phase is a genuine heavy-workload bottleneck, but none of the
three exact bounded candidates satisfies the predefined performance and
scientific gates. No prototype is suitable for production integration.

## 2. Scope And Provenance

| Item | Evidence |
|---|---|
| Audited production snapshot | `57d7cccb3889fed603b8454c6b8963d28bc33705` |
| Branch | `codex/fastphylosig-0.2.0-dev` |
| Package | fastphylosig `0.2.0.9000` |
| Runtime | R 4.6.1 UCRT, Windows x86-64 |
| Compiler | GCC 14.3.0 through the installed Rtools toolchain |
| Formal timing | warmup, serialized cells, alternating order, 10 paired repeats |
| Production files changed | none in `R/`, `src/`, `tests/`, `DESCRIPTION`, or `NAMESPACE` |

All candidates are private `Rcpp::sourceCpp()` prototypes under
`benchmarks/stage2d2/prototype/`. Each oracle and candidate was compiled in
the same translation unit. The ASCII staging directory was intentionally not
a Git checkout, so raw runner provenance records `source_commit=UNAVAILABLE`;
the production source was copied from the Git snapshot above. Prototype
SHA-256 values are recorded by the machine-readable provenance files.

## 3. Phase Decomposition

Generation-only performs identity initialization and the frozen descending
Fisher-Yates stream. Evaluation-only consumes pre-generated controlled
permutations and calls the production-equivalent K evaluator. Full-pipeline
timing includes both. These are independent timings and are not assumed to be
strictly additive.

Heavy-cell medians are seconds; parentheses contain IQR.

| n | nsim | Threads | Generation | Evaluation | Full | Generation | Evaluation | Integration residual |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5,000 | 9,999 | 1 | 0.596 (0.010) | 2.119 (0.024) | 2.782 (0.033) | 21.4% | 76.1% | 2.4% |
| 5,000 | 9,999 | 2 | 0.584 (0.021) | 1.100 (0.014) | 1.882 (0.021) | 31.1% | 58.4% | 10.5% |
| 10,000 | 9,999 | 1 | 1.169 (0.017) | 4.655 (0.019) | 5.941 (0.022) | 19.7% | 78.4% | 2.0% |
| 10,000 | 9,999 | 2 | 1.186 (0.023) | 2.387 (0.037) | 3.753 (0.018) | 31.6% | 63.6% | 4.8% |

The evaluation fraction exceeds the 50% authorization threshold in every
representative heavy cell. OpenMP one- and two-thread execution and RNG replay
both passed. P/MCSE locals in the phase harness are not part of the returned
checksum and may be compiler-elided, so the small integration residual is not
treated as separately authoritative. This limitation cannot overturn the
evaluation gate because evaluation alone remains 58.4-78.4% of full time.

## 4. `compute_one()` Numerical Audit

Each replicate creates four local numeric workspaces:

| Workspace | Length per active trait chunk | Purpose |
|---|---:|---|
| `message` | `N * c` doubles | upward Gaussian messages |
| `state` | `N * c` doubles | downward conditional states |
| `baseline` | `c` doubles | baseline-relative precision control |
| `delta` | `c` doubles | phylogenetic GLS mean offset |

For each trait and replicate, the kernel performs `5*n + 1` indexed trait
reads and `2*N + 3*e` structural node-loop visits. Its source-level arithmetic
proxy is `12*n + 14*e + 2*i + 3`, excluding addressing, finite checks, casts,
and hardware fusion.

The tree cache already hoists `parent`, CSR children, traversal orders,
branches, `aggregate`, `outgoing`, `sum_inv`, and K normalization. The GLS
mean, messages, states, numerator, and residual phylogenetic energy are
permutation-dependent and cannot be replaced with set-level trait moments.
The only clear residual tree-only arithmetic is:

- `alpha = branch[node] * outgoing[node]`, twice per internal non-root node
  and trait chunk;
- repeated branch divisions whose replacement with reciprocal multiplication
  would alter floating operation form and therefore was not attempted.

## 5. Candidate A: Reusable Numeric Workspaces

Candidate A allocates worker-local `message`, `state`, `baseline`, and `delta`
buffers once, then resets the fields required by the next replicate. Debug
fill and canary checks were used to detect stale state and bounds errors.

Exactness result:

| Gate | Result |
|---|---|
| Bounded exactness checks | 1,410/1,410 PASS |
| Ordered replicate rows | 103,275 exact |
| Null K, observed K, P, MCSE, exceedance | maximum absolute difference 0 |
| Requested/successful/failed accounting | exact PASS |
| RNG replay, thread/chunk settings | exact PASS |
| Input immutability | PASS |
| Malformed controlled permutations | 3/3 PASS |
| Workspace safety cells | 54/54 PASS |
| Stale-state/canary failures | 0/0 |

Formal one-thread heavy medians are seconds:

| Shape | n | Oracle | Candidate A | Speedup | Reduction |
|---|---:|---:|---:|---:|---:|
| balanced | 5,000 | 2.375 | 2.350 | 1.011x | 1.1% |
| balanced | 10,000 | 4.970 | 4.880 | 1.018x | 1.8% |
| random | 5,000 | 2.895 | 2.850 | 1.016x | 1.6% |
| random | 10,000 | 6.185 | 6.165 | 1.003x | 0.3% |
| pectinate | 5,000 | 3.160 | 2.990 | 1.057x | 5.4% |
| pectinate | 10,000 | 5.960 | 5.835 | 1.021x | 2.1% |

No cell approaches the required 20% median reduction; the maximum is 5.4%.
Candidate A is exact but fails the performance gate and is rejected.

## 6. Candidate B: Blocked Evaluation

The phase gate authorized a bounded blocked prototype, but the implementation
failed before performance testing. Its attempted block passes are not a valid
replacement for production `compute_one()`: they use different traversal and
arithmetic, discard their own numerical result, and call the original
`compute_one()` separately for every replicate. Its P/MCSE/accounting logic
also does not reproduce the frozen identity-first contract.

The bounded smoke gate recorded 450 checks and 201 failures. Ordered null K
appeared unchanged only because the prototype fell back to production
`compute_one()` for the returned K values; P differed by as much as 0.1667,
MCSE differed by as much as 0.0314 in the displayed failing cells, and failure
semantics did not match. Per the exact-first rule, no performance benchmark
was run. Candidate B is rejected for scientific-contract failure, not merely
for insufficient speed.

## 7. Candidate C: Exact Invariant Hoisting

Candidate C narrowly precomputes the existing double-precision expression
`branch[node] * outgoing[node]` once per engine call and reuses it in the two
downward passes. It does not replace division with multiplication, rearrange
reductions, alter traversal order, or hoist permutation-dependent quantities.

Its exactness gate passed all 1,410 checks and all 103,275 ordered replicate
rows with zero difference in K, null K, P, MCSE, exceedance, or accounting.
Formal one-thread heavy medians were:

| Shape | n | Oracle | Candidate C | Speedup | Reduction |
|---|---:|---:|---:|---:|---:|
| balanced | 5,000 | 2.325 | 2.315 | 1.004x | 0.4% |
| balanced | 10,000 | 4.780 | 4.770 | 1.002x | 0.2% |
| random | 5,000 | 2.795 | 2.470 | 1.132x | 11.6% |
| random | 10,000 | 5.980 | 5.280 | 1.133x | 11.7% |
| pectinate | 5,000 | 2.875 | 2.890 | 0.995x | -0.5% |
| pectinate | 10,000 | 5.890 | 5.805 | 1.015x | 1.4% |

Across both one- and two-thread runs, the best reduction was 11.7% and the
best speedup was 1.133x. No cell reached 20% or 1.25x, and the effect is not
topology-general. Candidate C is therefore rejected despite exactness.

## 8. Scientific And API Guards

Candidates A and C preserved the ordered null sequence, inclusive
`sim_K >= observed K` tail, P, MCSE, identity-first replicate, requested and
successful counts, return ordering, RNG replay, thread/chunk behavior, trait
masks, retained sets, large-offset traits, near-constant traits, constant
failure, two-tip cases, and batch traits. The `test=FALSE` public matrix guard
covered 1/8/32/100 traits. The formal runner reported:

```text
light_nsim199_guard = PASS
test_false_batch_guard = PASS
parity_gate = PASS
warnings_or_errors = NONE
```

No tolerance, estimator, RNG stream, thread policy, public result, or API was
changed. Candidate B failed these gates and is excluded.

## 9. Memory And Parallel Scaling

Candidate A uses workspace bounded by worker count, total nodes, and active
trait chunk. Reported one-trait workspace payload was about 320 KB at
`n=5,000` and 640 KB at `n=10,000`, with eight workspace allocations and
10,000 replicate reuses per heavy call. It does not allocate an
`n_tip * nsim` state matrix.

Candidate A was the best fully safe reusable-workspace prototype, so it was
used for the required parallel diagnostic. Candidate-versus-oracle speedups
remained negligible:

| n | Threads | Oracle | Candidate A | Speedup |
|---:|---:|---:|---:|---:|
| 5,000 | 2 | 1.910 | 1.890 | 1.011x |
| 5,000 | 4 | 1.370 | 1.375 | 0.996x |
| 5,000 | 8 | 1.055 | 1.025 | 1.029x |
| 10,000 | 2 | 3.985 | 3.980 | 1.001x |
| 10,000 | 4 | 2.765 | 2.770 | 0.998x |
| 10,000 | 8 | 2.215 | 2.110 | 1.050x |

The unchanged oracle itself scales with supported OpenMP threads, but
workspace reuse does not improve that scaling enough to justify production
complexity. No scheduler, thread policy, or RNG partition was changed.

## 10. Candidate Matrix

| Candidate | Exact gate | Performance gate | Decision |
|---|---|---|---|
| A: thread-local reusable workspaces | PASS | FAIL, maximum reduction 5.4% | REJECT |
| B: blocked/batched evaluator | FAIL | NOT RUN by rule | REJECT |
| C: exact tree-invariant `alpha` hoist | PASS | FAIL, maximum reduction 11.7% | REJECT |

There is no exact candidate meeting the formal performance threshold. The
smallest production proposal is therefore the empty set.

This is the required hard stop for K optimization in 0.2.0. Do not continue
with GPU, SIMD rewrites, approximate K, alternative RNG, unsafe validation
paths, or another K candidate. The next bounded investigation is the D null
engine.

## Evidence Index

- [phase protocol](benchmarks/stage2d2/PHASE_DECOMPOSITION_PROTOCOL.md)
- [`compute_one()` audit](benchmarks/stage2d2/COMPUTE_ONE_AUDIT.md)
- [phase budget](benchmarks/stage2d2/results/phase-audit/stage2d2_phase_budget.csv)
- [phase status](benchmarks/stage2d2/results/phase-audit/stage2d2_phase_status.csv)
- [Candidate A exactness](benchmarks/stage2d2/results/correctness-a/stage2d2_correctness_status.csv)
- [Candidate A performance](benchmarks/stage2d2/results/candidate-a-benchmark/stage2d2_candidate_benchmark_summary.csv)
- [workspace safety](benchmarks/stage2d2/results/workspace-safety/stage2d2_workspace_safety_status.csv)
- [Candidate A parallel scaling](benchmarks/stage2d2/results/candidate-a-parallel/stage2d2_candidate_benchmark_summary.csv)
- [Candidate B failed gate](benchmarks/stage2d2/results/correctness-b-smoke/stage2d2_correctness_status.csv)
- [Candidate C exactness](benchmarks/stage2d2/results/correctness-c/stage2d2_correctness_status.csv)
- [Candidate C performance](benchmarks/stage2d2/results/candidate-c-benchmark/stage2d2_candidate_benchmark_summary.csv)
- [ordered-null evidence manifest](benchmarks/stage2d2/results/ORDERED_NULLS_MANIFEST.csv)
- [machine-readable final decision](benchmarks/stage2d2/results/FINAL_DECISION.csv)
