# fastphylosig 0.2.0 Stage 2E1A Expert Review

## 1. Final decision

```text
COMMON_PREPARATION_ENGINE_0_2_0 = FROZEN
CANONICALIZATION_CONTRACT_V2 = ACCEPTED_FROZEN
K_OPTIMIZATION_0_2_0 = CLOSED_FINAL
lambda optimization = CLOSED_FOR_0_2_0
Delta new optimization = CLOSED_FOR_0_2_0

D_EXACT_ORDER_STATISTIC = NO_GO
D_OPTIMIZATION_0_2_0 = CLOSED
NEXT_STAGE = FINAL_0_2_0_VALIDATION
PRODUCTION_CODE_CHANGED = NO
```

The exact order-statistic candidate is not authorized. In all 18 representative
heavy cells, sort plus type-7 thresholding consumed less than 20% of complete
prepared D wall time. The observed range was 15.88-18.73%; no cell reached the
20% borderline band and none reached the 30% GO gate. This closes the last
statistic-specific performance direction for 0.2.0.

## 2. Evidence and provenance

| Item | Result |
|---|---|
| Audited Git source | `03d72813f9da08b27f924aa37437ab2715cf8dfe` |
| Package | `fastphylosig 0.2.0.9000` |
| Runtime | R 4.6.1 UCRT, Windows x86-64 |
| Formal design | fixed fixtures, warmup, serialized, alternating order |
| Formal grid | 3 shapes x 3 n x 3 prevalence x 2 nsim = 54 cells |
| Repeats | 10 per cell |
| Timing rows | 540/540 PASS, 0 warnings/errors |
| Representative heavy cells | 18 (`n >= 2000`, `nsim = 9999`) |
| Correctness harness | 1299/1299 PASS |
| Authoritative route | prepared D, `ncores = 1` |
| Additional observation | prepared D, `ncores = 2` |
| Production files changed | none |

The private C++ harness includes the frozen production `src/d_stream.cpp` in
the same translation unit. Static source evidence was read from an ASCII copy
of the same source snapshot because the local Windows R locale could not read
the Chinese repository path reliably. SHA-256 values for `d_stream.cpp`,
`fast_signal.cpp`, and `numeric_utils.h` are recorded in the source audit CSV;
the Git provenance is recorded separately from the real repository.

## 3. Frozen D contract

The following contract was extracted from current production source, tests,
and documentation and was not changed.

- Observed D uses the frozen inverse-branch-length contrast traversal. Final
  D is `(observed - mean_brownian) / (mean_random - mean_brownian)`.
- The random null permutes the retained, already recoded binary state vector.
  Its strict tail is `random_null < observed`.
- The Brownian null draws continuous tip values with R's RNG, obtains the
  current type-7 threshold, and maps `value < threshold` to state 1. Its strict
  tail is `brownian_null > observed`.
- Binary recoding retains the current factor-level order and prevalence
  definition. NA handling, species matching, retained-tip masks, and
  fewer-than-two-species failure semantics are unchanged.
- `P_random`, `P_Brownian`, their aliases, MCSE, and requested/successful/failed
  accounting retain their current formulas and result fields.
- RNG draws remain serial. Contrast columns may use the frozen static OpenMP
  schedule; no RNG, scheduler, or thread-policy change was made.
- Rooted polytomies retain the compatibility path; unary and invalid trees
  retain their existing rejection behavior.
- Literal strict tails remain sensitive to a one-ULP observed/null boundary.
  Equality is not counted, and no epsilon was introduced.
- Tree and trait inputs remain immutable, and prepared-context mutation checks
  remain owned by the frozen common preparation engine.

The formal harness directly checked Brownian states, type-7 thresholds, binary
states, null D values, strict exceedance directions, P, MCSE, accounting,
same-seed replay, one-/two-thread replay, public prepared `fast_d()` parity,
rooted-polytomy compatibility, and one-ULP strict-tail behavior.

## 4. Type-7 requirement

For ascending values of length `n`, production computes:

```text
h     = 1 + (n - 1) * prevalence
lo    = floor(h)
gamma = h - lo
```

It uses the minimum at the lower endpoint, the maximum at the upper endpoint,
and otherwise only `values[lo - 1]` and `values[lo]`. Thus type-7 needs at most
two adjacent order statistics. The complete sorted vector has no later use in
observed D, null D, P, MCSE, accounting, RNG, or result packaging.

An exact non-probabilistic selection design is technically feasible, but this
fact is only operation evidence. It does not override the wall-time gate and
does not authorize Stage 2E1B.

## 5. Representative heavy cells

Times are median complete prepared D seconds with IQR in parentheses. Fractions
are independent phase medians divided by the authoritative `ncores=1` median.
They are descriptive and are not forced to sum to 100%. `Sort/Brownian` is the
full-sort share of the complete Brownian path; `Sort+Q/total` is the qualifying
sort-plus-type-7 share of the complete D pipeline.

| Shape | n | Prev. | Total s (IQR) | Random | Brownian gen. | Sort/Brownian | Sort+Q/total | Binary+D |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| balanced | 2000 | 0.10 | 5.760 (0.020) | 35.3% | 36.9% | 33.8% | 18.1% | 20.4% |
| random | 2000 | 0.10 | 5.915 (0.227) | 29.3% | 36.9% | 33.3% | 18.3% | 19.9% |
| pectinate | 2000 | 0.10 | 5.485 (0.018) | 29.4% | 38.8% | 30.3% | 16.0% | 20.2% |
| balanced | 5000 | 0.10 | 16.010 (0.260) | 25.6% | 34.1% | 34.9% | 18.4% | 21.5% |
| random | 5000 | 0.10 | 16.095 (0.345) | 25.6% | 34.2% | 35.0% | 18.0% | 20.8% |
| pectinate | 5000 | 0.10 | 15.080 (0.225) | 27.2% | 36.5% | 32.4% | 15.9% | 22.1% |
| balanced | 2000 | 0.50 | 5.860 (0.025) | 35.5% | 37.1% | 34.2% | 18.1% | 18.9% |
| random | 2000 | 0.50 | 5.965 (0.107) | 33.9% | 36.8% | 33.6% | 18.1% | 19.4% |
| pectinate | 2000 | 0.50 | 5.510 (0.055) | 35.6% | 38.9% | 31.7% | 17.6% | 21.2% |
| balanced | 5000 | 0.50 | 16.395 (0.200) | 24.3% | 33.9% | 35.8% | 18.0% | 21.0% |
| random | 5000 | 0.50 | 16.645 (0.122) | 24.7% | 33.6% | 35.3% | 17.7% | 20.5% |
| pectinate | 5000 | 0.50 | 14.995 (0.205) | 27.4% | 36.6% | 32.7% | 16.8% | 21.6% |
| balanced | 2000 | 0.90 | 5.830 (0.040) | 34.7% | 37.4% | 34.1% | 18.4% | 20.7% |
| random | 2000 | 0.90 | 5.845 (0.073) | 34.8% | 37.7% | 33.0% | 18.0% | 20.0% |
| pectinate | 2000 | 0.90 | 5.525 (0.018) | 35.3% | 39.3% | 31.8% | 18.7% | 21.7% |
| balanced | 5000 | 0.90 | 16.310 (0.282) | 25.0% | 33.9% | 36.3% | 17.9% | 21.2% |
| random | 5000 | 0.90 | 16.590 (0.180) | 24.6% | 33.0% | 35.1% | 17.5% | 20.8% |
| pectinate | 5000 | 0.90 | 14.965 (0.193) | 27.5% | 36.3% | 32.7% | 16.2% | 21.9% |

Across these cells, the median sort-plus-threshold fraction was 17.99%. The
paired repeat-level estimate was nearly identical: 15.95-18.73%, median
17.93%. This agreement, together with one full sort per Brownian replicate and
no downstream use of the sorted vector, provides independent timing, paired,
and operation-count support for the conclusion.

## 6. Speedup ceilings

| Assumption | Heavy-cell speedup range | Median | Best-case wall-time saving |
|---|---:|---:|---:|
| Full sort is free | 1.204-1.240x | 1.229x | 19.4% |
| Selection costs 20% of old sort | 1.157-1.183x | 1.175x | 15.5% |
| Selection costs 30% of old sort | 1.135-1.157x | 1.150x | 13.5% |
| Selection costs 50% of old sort | 1.093-1.107x | 1.103x | 9.7% |

These are Amdahl ceilings, not predicted production gains. They assume zero
new allocation, validation, copying, compiler, and platform cost. The 50%
residual-cost scenario remains below the required 15% practical saving, and
the qualification share already fails the independent `<20%` stop rule.

## 7. Parallel observation

`ncores=2` was observation-only. In the 18 heavy cells it took
`1.665-2.105x` the `ncores=1` wall time (median `1.779x`), so two threads were
slower in this current workload. This does not authorize a scheduler or
OpenMP redesign, and the authoritative candidate decision remains based on
`ncores=1`.

## 8. Expert conclusion

Brownian full sorting is visible within the Brownian subpath, accounting for
30.34-36.28% of that subpath. It is not large enough in complete prepared D:
sort plus type-7 thresholding is only 15.88-18.73% of total D time across all
representative heavy cells. Brownian generation, random-null work, and
binary-state D traversal each remain substantial independent costs.

The exact order-statistic replacement therefore fails the predeclared wall-
time gate. Random-null cost is also not a candidate: this audit identified no
exact algorithmic redundancy and did not establish a realistic total saving
of at least 15%. No D prototype, buffer reuse candidate, or parallel redesign
is authorized.

Stage 2E1A ends with `NO_GO`. All statistic-specific performance optimization
for fastphylosig 0.2.0 is closed. The only permitted next step is final 0.2.0
release validation.

## Evidence index

- [frozen D contract](benchmarks/stage2e1a/D_CONTRACT_AUDIT.md)
- [type-7 requirement](benchmarks/stage2e1a/TYPE7_REQUIREMENT_AUDIT.md)
- [harness protocol](benchmarks/stage2e1a/HARNESS_PROTOCOL.md)
- [performance protocol](benchmarks/stage2e1a/PERFORMANCE_PROTOCOL.md)
- [authoritative correctness status](benchmarks/stage2e1a/results/correctness-authoritative/stage2e1a_correctness_status.csv)
- [authoritative correctness checks](benchmarks/stage2e1a/results/correctness-authoritative/stage2e1a_correctness.csv)
- [formal benchmark status](benchmarks/stage2e1a/results/formal-authoritative/stage2e1a_benchmark_status.csv)
- [formal timing summary](benchmarks/stage2e1a/results/formal-authoritative/stage2e1a_phase_summary.csv)
- [heavy qualification](benchmarks/stage2e1a/results/formal-authoritative/stage2e1a_heavy_qualification.csv)
- [heavy aggregate](benchmarks/stage2e1a/results/formal-authoritative/stage2e1a_heavy_aggregate.csv)
- [source type-7 audit](benchmarks/stage2e1a/results/formal-authoritative/stage2e1a_source_type7_audit.csv)
- [machine-readable final decision](benchmarks/stage2e1a/results/formal-authoritative/FINAL_DECISION.csv)
