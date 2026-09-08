# fastphylosig 0.2.0 Stage 2C Expert Review Template

This is a report skeleton for the post-V2 full performance reprofile. Replace
every bracketed value with evidence from one Stage 2C run. Do not mark a gate
PASS from a theoretical estimate or from timings copied from another commit.

## 0. Frozen state

```text
PREPARATION_ENGINE_STAGE1 = ACCEPTED_FROZEN
RAW_ROUTE_SINGLE_PREPARATION_BOUNDARY = ACCEPTED_FROZEN
INSPECTION_INTEGER_ADJACENCY = ACCEPTED_FROZEN
CANDIDATE_2 = REJECTED
CANONICALIZATION_CONTRACT_V2 = ACCEPTED_FROZEN
PUBLIC_RELEASE_BLOCKER_0_2_0 = CLEARED
COMMON_PREPARATION_ENGINE_0_2_0 = FROZEN
```

Stage 2C production change: `NONE`.

## 1. Provenance

| Item | Value |
|---|---|
| Full commit | [full 40-character hash] |
| Package version | [value] |
| R version/platform | [value] |
| OS and machine | [value] |
| Compiler | [value] |
| OpenMP status | [not exercised / exact status] |
| BLAS/LAPACK | [value] |
| LC_ALL / LC_COLLATE | [value] |
| Hardware | [CPU/core information] |
| Formal output directory | [path] |

## 2. Unified Stage 2C table

All estimator speedups must be derived from paired/reference timings in the
same Stage 2C run. `RESOURCE_LIMIT` is a valid result and is not a speedup.
For context-only rows, use `memory/context cost` and do not invent estimator
fractions.

| method | workload | raw total | prepare | prepared total | validation | kernel/optimizer | simulation | startup | packaging | reference time | speedup | memory/context cost |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| K | [n/traits/test/nsim] | [ms] | [ms/%] | [ms] | [ms/%] | [ms/%] | [ms/%] | [ms/%] | [ms/%] | [ms or RESOURCE_LIMIT] | [x or NA] | [bytes] |
| lambda | [shape/n] | [ms] | [ms/%] | [ms] | [ms/%] | [ms/%] | [NA] | [NA] | [ms/%] | [ms or RESOURCE_LIMIT] | [x or NA] | [bytes] |
| D | [n/nsim/workers] | [ms] | [ms/%] | [ms] | [ms/%] | [ms/%] | [ms/%] | [ms/%] | [ms/%] | [ms or RESOURCE_LIMIT] | [x or NA] | [bytes] |
| Delta | [n/nsim/workers] | [ms] | [ms/%] | [ms] | [ms/%] | [ms/%] | [ms/%] | [ms/%] | [ms/%] | [ms or RESOURCE_LIMIT] | [x or NA] | [bytes] |
| context | [balanced n] | [NA] | [ms] | [NA] | [read + validation ms] | [NA] | [NA] | [NA] | [NA] | [NA] | [NA] | [object/cache/RDS bytes] |

## 3. K single-trait crossover

Report the same fixed fixtures for `phytools::phylosig(method="K",
test=FALSE)`, `fast_k(raw)`, `prepare_tree() + fast_k(prepared)`, and
`fast_signal(raw, method="K")` at n = 50, 100, 300, 500, 1000, 2000, 5000,
10000, and 20000. State separately:

```text
prepared compute = [value]
prepared end-to-end = prepare + compute = [value]
```

Do not label prepared compute as one-shot user performance.

## 4. K batch and permutation budget

For batch workloads, report raw/prepared matrix and column-loop timings for 1,
8, 32, and 100 traits at n = 500, 2000, and 5000. For permutation workloads,
report raw/prepared at nsim = 199, 999, and 9999 where resources allow, with
the following fractions:

```text
simulation_fraction = [value]
null_K_fraction = [value]
RNG_fraction = [value]
memory_copy_fraction = [value]
```

`K_PERMUTATION_OPTIMIZATION` is `AUTHORIZED_FOR_CANDIDATE_DESIGN` only when
permutation/null evaluation is at least 50% of heavy K wall time; at least 70%
is `HIGH_PRIORITY`. This stage profiles/ranks only and implements nothing.

## 5. Prepared-K V2 regression audit

Decompose the observed prepared-path change into external protected-state
validation, canonical fingerprint/snapshot comparison, cache lookup, trait
handling, retained-subtree lookup, K kernel, and result packaging. Record the
confirmed regression threshold and evidence count. Integrity coverage and
mutation detection must not be weakened to improve this number.

## 6. Lambda budget

Report raw, prepared, and reference paths for balanced, random, and pectinate
trees at n = 100, 500, 2000, 5000, and 10000. Include preparation, context
validation, unique likelihood evaluation count, fixed-lambda kernel, optimizer,
and result construction. Confirm that the pectinate recursion failure is
absent. Any failure is reported with its exact status and call context.

## 7. D budget

Report raw/prepared at n = 100, 500, and 2000 with nsim = 199, 999, and 9999
where resources allow, including observed D, random null, Brownian null, RNG,
context validation, parallel startup, and packaging. Include one worker and
the current public parallel policy where supported. Do not change RNG streams.

## 8. Delta budget

Report `test=FALSE` and permutation-heavy nsim = 19 and 199 at n = 50, 100,
and 500 where resources allow, using the current MCMC settings. Decompose ACE
observed fit, MCMC, diagnostics, permutation ACE, worker startup, context
validation, and packaging for serial and the current two-worker configuration.
Do not use short runs as biological inference.

## 9. Context memory and persistence

Attach `stage2c_context_memory.csv`, `stage2c_provenance.csv`, and
`stage2c_sessionInfo.txt`. Report the measured values for n = 500, 1000, 5000,
10000, and 20000:

| n | context object bytes | fingerprint bytes | protected snapshot bytes | cache payload bytes | saveRDS bytes | saveRDS ms | fresh readRDS + first validation ms |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 500 | [value] | [value] | [value] | [value] | [value] | [value] | [value] |
| 1000 | [value] | [value] | [value] | [value] | [value] | [value] | [value] |
| 5000 | [value] | [value] | [value] | [value] | [value] | [value] | [value] |
| 10000 | [value] | [value] | [value] | [value] | [value] | [value] | [value] |
| 20000 | [value] | [value] | [value] | [value] | [value] | [value] | [value] |

Use **measured object/cache bytes**. Do not call these values peak RSS or
total process memory unless a separate peak-RSS measurement exists.

## 10. Correctness guard

Before interpreting performance, record K parity, lambda parity, D controlled-
null parity, Delta fixed-seed parity, raw/prepared parity, warning/status
parity, and any anomaly. An unexplained correctness anomaly changes the report
to `INCOMPLETE` and stops candidate interpretation.

## 11. Candidate ranking (maximum two)

List no more than two statistic-specific candidates. Each candidate must have:

| Candidate | Measured hotspot share | Ideal end-to-end ceiling | Estimated practical saving | Correctness risk | Engineering risk | GO/NO-GO |
|---|---:|---:|---:|---|---|---|
| [name] | [>=50% or NO-GO] | [>=20% or NO-GO] | [value] | [value] | [value] | [GO/NO-GO] |
| [name] | [>=50% or NO-GO] | [>=20% or NO-GO] | [value] | [value] | [value] | [GO/NO-GO] |

Candidate authorization requires measured hotspot >= 50% of wall time, ideal
end-to-end saving of approximately 20% or more, a clear route that preserves the statistical
contract, no larger common bottleneck, and reasonable complexity/maintenance
cost. Allowed categories are K permutation/null engine, prepared exact-
validation efficiency, lambda likelihood/optimizer path, D random/Brownian
null engine, Delta permutation ACE engine, and Delta worker startup/scheduler.
This section authorizes design at most; it does not authorize implementation.

## 12. Final questions and gate

Answer these twelve questions explicitly:

1. Where is the 0.2.0 K crossover?
2. Is K `test=FALSE` still worth optimizing?
3. What is the real K `test=TRUE` hotspot?
4. What explains the prepared V2 regression?
5. Where is lambda mainly slow now?
6. Where is D mainly slow now?
7. Where is Delta mainly slow now?
8. Is context memory a real limitation?
9. Which at most two candidates are worth final development?
10. For each candidate, what are hotspot share, ideal ceiling, practical saving, correctness risk, engineering risk, and GO/NO-GO?
11. Is `STAGE2C = PASS` or `INCOMPLETE`?
12. What is `NEXT_STAGE` (at most two candidates)?

Final status:

```text
STAGE2C = [PASS / INCOMPLETE]
NEXT_STAGE = [at most two named candidates / STOP]
PRODUCTION_CODE_CHANGED = NO
```

Do not implement any candidate automatically after completing this review.
