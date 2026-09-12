# Stage 2E1A performance and decision protocol

Stage 2E1A is the final D-specific profiling qualification before release
validation. It is an audit-only stage. The runner does not change `R/`, `src/`,
tests, package metadata, public APIs, RNG code, OpenMP policy, or production
estimators. No D candidate is implemented by this stage.

## Frozen boundary

The runner assumes the following production state is frozen:

```text
COMMON_PREPARATION_ENGINE_0_2_0 = FROZEN
CANONICALIZATION_CONTRACT_V2 = ACCEPTED_FROZEN
K_OPTIMIZATION_0_2_0 = CLOSED_FINAL
K blocked evaluator = REJECTED_NOT_ENOUGH_GAIN
lambda optimization = CLOSED_FOR_0_2_0
Delta optimization = CLOSED_FOR_0_2_0
```

Only D null timing is examined. K, lambda, and Delta are not reopened.

## Invocation

The formal run is serialized in one R process:

```text
Rscript benchmarks/stage2e1a/run_stage2e1a_benchmark.R --formal <repo> <out>
```

The small protocol check is:

```text
Rscript benchmarks/stage2e1a/run_stage2e1a_benchmark.R --smoke <repo> <out>
```

Formal defaults are:

```text
tree shapes: balanced, random, pectinate
n:            500, 2000, 5000
optional n:  10000 when FASTPHYLOSIG_STAGE2E1A_INCLUDE_10000=true
nsim:         999, 9999
prevalence:   0.10, 0.50, 0.90
repeats:      10 paired repetitions per cell
authoritative ncores: 1
observation ncores:   2
chunk size:   128 columns for phase probes
```

Environment variables may narrow or extend a local run. A run that narrows
the formal grid is not evidence that the complete formal grid passed; the
provenance records the actual grid.

## Correctness gate

No timing is interpreted before the correctness files report `status=PASS`.
The gate uses a fixed 13-tip balanced tree and controlled permutations. It
checks all of the following against the current production functions:

1. Brownian continuous-state generation follows the same R RNG stream and
   traversal order.
2. Production `brownian_tree_threshold_cpp()` and the audit translation unit
   produce identical type-7 thresholded binary states.
3. The audit sort-only, sort-plus-threshold, and thresholded binary probes have
   the expected dimensions and exact binary output.
4. Production random null states, Brownian states, observed D, means, D, both
   tail P values, both MCSE values, successful/failed counts, and retained null
   vectors are identical for a controlled fixed input.
5. The public prepared `fast_d()` path replays exactly under the same seed.
6. The strict tail definition is checked directly (`random < observed`,
   `Brownian > observed`).

Any mismatch or warning/error fails closed and no formal timing is authorized.

## Independent phase probes

The runner reports independent measurements. Their elapsed times are not forced
to add to the public total because the probes use separate calls and bounded
128-column chunks to avoid retaining an `n x nsim` continuous-state matrix.
This is intentional: phase rows are diagnostic measurements, while complete
prepared `fast_d()` timing is the wall-time reference.

The probes are:

| Label | Measurement |
| --- | --- |
| `random_null_s` | Production `.phylo_d_random_states()` on bounded chunks, without replacement-weight changes. |
| `brownian_generation_s` | Audit-only production-equivalent Brownian continuous-state traversal, before sorting. |
| `brownian_sort_s` | Full ascending sort of pre-generated continuous columns only. |
| `brownian_sort_threshold_s` | Full sort plus the current type-7 threshold arithmetic, without binary conversion. |
| `binary_D_s` | Thresholded binary conversion plus production `phylo_d_sums_cpp()` for Brownian columns. |
| `complete_brownian_s` | Production `.phylo_d_brownian_states()` including generation, full sort, type-7 threshold, and binary conversion. |
| `complete_D_ncores1_s` | Prepared public `fast_d(test=TRUE, nsim=..., ncores=1)`; authoritative total. |
| `complete_D_ncores2_s` | The same prepared public call with `ncores=2`; observation only. |

Every cell is warmed once, run serially, and measured with fixed fixtures. The
phase order alternates on successive repetitions to reduce drift. No worker
processes are created by the benchmark.

## Reported fractions

Fractions use medians from the same cell and authoritative `ncores=1` public
total:

```text
sort_fraction_of_complete_Brownian = median(brownian_sort_s) /
                                     median(complete_brownian_s)

sort_threshold_fraction_of_total_D  = median(brownian_sort_threshold_s) /
                                     median(complete_D_ncores1_s)

random_fraction_of_total_D          = median(random_null_s) /
                                     median(complete_D_ncores1_s)

Brownian_generation_fraction        = median(brownian_generation_s) /
                                     median(complete_D_ncores1_s)

binary_D_fraction                   = median(binary_D_s) /
                                     median(complete_D_ncores1_s)
```

The independent component sum is recorded for diagnosis only. It is not an
authoritative decomposition identity.

## Type-7 and exact-selection contract

The current production Brownian helper performs one ascending full sort and
then calls the existing type-7 helper. For a probability `p`, type-7 uses

```text
h = 1 + (n - 1) * p
lo = floor(h)
gamma = h - lo
values[lo - 1] and values[lo]
```

The runner performs a static source audit for this contract and also verifies
the resulting binary matrix exactly. It does not change threshold arithmetic.
The fact that only two adjacent order statistics are consumed after sorting
does not, by itself, authorize a replacement. An exact, non-probabilistic
selection design is not implemented in Stage 2E1A and therefore the runner
cannot return `GO` merely because sorting is expensive.

## Decision rules

The formal qualification uses representative heavy cells with `n >= 2000`
and `nsim = 9999`.

`D_EXACT_ORDER_STATISTIC = GO` would require all of:

1. At least two valid heavy cells with
   `sort_threshold_fraction_of_total_D >= 0.30`.
2. The source audit confirms that the current complete sort only feeds the
   exact type-7 order statistics used by production.
3. A separately reviewed exact non-probabilistic selection design exists.

Stage 2E1A does not implement condition 3, so it never silently authorizes a
candidate. A result in the 20%-30% range, or a result above 30% without the
design review, is `NEED_MORE_EVIDENCE` and goes to a design review only.

If the share is below 20% and no timing/equivalence failure exists:

```text
D_EXACT_ORDER_STATISTIC = NO_GO
D_OPTIMIZATION_0_2_0 = CLOSED
NEXT_STAGE = FINAL_0_2_0_VALIDATION
```

If the result is `NEED_MORE_EVIDENCE`, the runner sets
`D_OPTIMIZATION_0_2_0 = CONTINUE` only for the bounded design-review decision;
it does not compile or run a candidate.

## Speedup ceilings

For each cell the runner reports a theoretical free-sort ceiling:

```text
free_sort_speedup = 1 / (1 - sort_fraction_of_total_D)
```

It also reports conservative ceilings when an exact selection still pays
20%, 30%, or 50% of the old sort cost:

```text
speedup(q) = 1 / (1 - sort_fraction_of_total_D * (1 - q))
```

The 50% residual-sort value is the conservative practical summary. These are
ceilings, not promises: they assume no new allocation, conversion, validation,
cache, compiler, or platform cost.

The random-null fraction is reported but does not authorize a random-state
buffer candidate. Such a candidate would require independent proof of exact
algorithmic redundancy and at least 15% realistic end-to-end saving.

## Outputs

The output directory contains:

```text
correctness/stage2e1a_correctness_checks.csv
correctness/stage2e1a_correctness_status.csv
stage2e1a_source_type7_audit.csv
stage2e1a_phase_timings.csv
stage2e1a_phase_summary.csv
stage2e1a_decision.csv
stage2e1a_provenance.csv
stage2e1a_benchmark_status.csv
```

No output is a production artifact. `ncores=2` is observation only; OpenMP,
scheduler, and RNG policy are not modified. If `NO_GO` is reached, stop D
statistic-specific optimization and proceed directly to final 0.2.0
validation. Do not reopen K, lambda, Delta, or implement a D candidate from
this profiling run.
