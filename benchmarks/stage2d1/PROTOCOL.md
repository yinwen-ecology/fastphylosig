# fastphylosig 0.2.0 Stage 2D1 Protocol

## Scope and freeze

Stage 2D1 is a private design and prototype audit for the K permutation/null
engine. It does not modify `R/`, `src/`, `tests/`, the package namespace, the
public API, the statistical estimator, tolerance, RNG policy, thread policy,
or the common preparation engine.

The frozen inputs are:

```text
PREPARATION_ENGINE_STAGE1 = ACCEPTED_FROZEN
RAW_ROUTE_SINGLE_PREPARATION_BOUNDARY = ACCEPTED_FROZEN
INSPECTION_INTEGER_ADJACENCY = ACCEPTED_FROZEN
CANONICALIZATION_CONTRACT_V2 = ACCEPTED_FROZEN
COMMON_PREPARATION_ENGINE_0_2_0 = FROZEN
STAGE2C = PASS
K_PERMUTATION_OPTIMIZATION = AUTHORIZED_HIGH_PRIORITY
```

The prototype is test-only. No result from this directory authorizes
production integration automatically. Candidate A, B, and C must each be
accepted, rejected, or marked not worth complexity on their own evidence.

## Files and ownership

The audit files are confined to this directory:

- `stage2d1_common.R` contains fixtures, oracle adapters, result comparison,
  provenance, and workspace proxies.
- `run_stage2d1_correctness.R` runs controlled-permutation, RNG replay,
  edge-case, failure, NA-mask, input-immutability, and batch guards.
- `run_stage2d1_benchmark.R` runs timing only after a correctness PASS file is
  present and only when explicitly invoked with `--formal`.
- `K_RNG_CONTRACT_STAGE2D1.md` records the extracted RNG contract used by the
  scripts.

These scripts are audit harnesses. They must not be copied into package
`R/`, `src/`, or `tests/` as part of this stage.

## Prototype hook contract

The scripts can load a prototype from either:

```text
FASTPHYLOSIG_STAGE2D1_HOOK_FILE=<R file>
FASTPHYLOSIG_STAGE2D1_PROTO_CPP=<sourceCpp file>
```

When the environment variable is absent, the harness auto-discovers the
checked-in audit prototype at
`benchmarks/stage2d1/prototype/k_permutation_prototype.cpp` when that file is
present.

The file is sourced into a private environment. The default function names
are `stage2d1_proto_controlled`, `stage2d1_proto_rng`, and the optional
`stage2d1_proto_batch`. Function names can be changed with:

```text
FASTPHYLOSIG_STAGE2D1_PROTO_CONTROLLED
FASTPHYLOSIG_STAGE2D1_PROTO_RNG
FASTPHYLOSIG_STAGE2D1_PROTO_BATCH
```

The controlled and RNG hooks should accept the same named arguments as the
production C++ binding:

```r
function(compiled_tree, X, nsim, permutations, trait_chunk,
         return_sim, include_observed, n_threads, simulation_chunk, ...)
```

The RNG hook may additionally accept a seed argument. Its name defaults to
`seed` and can be changed with `FASTPHYLOSIG_STAGE2D1_SEED_ARG`. If the hook
does not declare that argument, the harness uses `set.seed()` around the call
and restores the caller's RNG state.

The prototype should return a list containing the following canonical fields,
or aliases recognized by `stage2d1_common.R`:

```text
K, P, MCSE_P, exceedance_count, nsim_requested,
nsim_successful, nsim_failed, success, sim_K
```

`sim_K` is expected to be `nsim x traits`; a transposed `traits x nsim`
matrix is accepted and normalized only for comparison. No n-tip by nsim
matrix is allowed inside the candidate. An optional
`stage2d1_workspace_bytes` field or attribute can be added for a measured
candidate workspace report; otherwise the harness records an analytic proxy.

## Frozen RNG contract

The contract used by this audit is deliberately narrower than a new RNG
promise:

- Same seed, same retained species order, same arguments, same
  `simulation_chunk`, and same `ncores` are required for old-versus-prototype
  replay.
- A supplied permutation matrix is a controlled-input path. It has `nsim`
  rows and one 1-based permutation of `1:n` per row.
- Internal RNG mode uses the current C++ Fisher-Yates implementation backed by
  R's RNG (`R::runif`). With the existing `include_observed` behavior, the
  first row is identity and subsequent rows are random permutations.
- The current source generates internal-RNG permutations serially before
  OpenMP evaluation. The source comment says this makes the draw stream
  independent of thread count and `simulation_chunk`; the harness records
  this separately and does not broaden the public contract beyond existing
  documentation and tests.
- The README requires seed, explicit stochastic inputs when supplied, and
  `ncores` to be held fixed for stochastic replay. A changed `ncores` value is
  not used to invent a new worker-count-invariant contract.

See `K_RNG_CONTRACT_STAGE2D1.md` for the source anchors and the exact wording
used by the harness.

## Current production hot-loop audit

The following is a source-level audit of the current `src/k_permutation.cpp`
path. Counts are per replicate unless stated otherwise; they describe the
old oracle and are not claims about a prototype implementation.

### Internal RNG mode

For each replicate in a generated block:

1. A temporary `std::vector<int> perm(n)` is constructed in the serial
   generation loop.
2. The identity vector is assigned into `perm` before Fisher-Yates, so one
   full index-length copy occurs. Fisher-Yates then performs the same draw
   order in place.
3. The completed permutation is copied into the bounded
   `chunk_perms` buffer, adding one full index-length copy per replicate.
4. In the OpenMP branch, each worker owns a `perm(n)` buffer. Each replicate
   is copied from `chunk_perms` into that worker buffer before evaluation.
   In the serial branch, one `perm(n)` buffer is reused but receives the same
   per-replicate row copy.
5. `compute_one()` receives the permutation by const reference and gathers
   trait values through the index vector. Its message, state, baseline, delta,
   and output workspaces are sized for each trait chunk on each call.
6. Exceedance is updated once per trait after the null K value is computed.
   The optional `sim_K` matrix is allocated once for the whole call and is
   written by replicate/trait; it is not allocated per replicate.

The bounded generated-index payload is approximately
`4 * n * min(simulation_chunk, nsim)` bytes, plus one active `4 * n` worker or
serial permutation buffer. This is an analytic proxy, not peak RSS.

### Controlled permutation mode

The caller supplies the permutation matrix. The C++ binding validates every
row, then the serial path allocates one reusable `perm(n)` and copies one row
per replicate. The OpenMP path allocates one `perm(n)` per worker and copies
one row per replicate. There is no R-level Fisher-Yates draw in this branch.
The supplied matrix itself is caller-owned input and is reported separately
from candidate workspace.

### Candidate interpretation

- Candidate A may remove repeated generated-index allocation/copy only if it
  preserves identity reset, Fisher-Yates draw count/order, and the current
  chunk/RNG behavior.
- Candidate B may remove only a proven redundant full-index copy. An alias is
  not acceptable if a later Fisher-Yates or gather operation can mutate the
  values used for returned `sim_K` or failure handling.
- Candidate C may fuse index-based trait gather with K evaluation only if the
  floating-point operation order and all exact/tolerance gates remain valid.
  A mismatch in any ordered null K rejects C; no tolerance is widened.

## Correctness invocation

Correctness must be run before any formal benchmark. It is explicit to avoid
accidental heavy execution:

```text
Rscript run_stage2d1_correctness.R <repo> <output> --formal
```

A bounded smoke run is available for hook wiring:

```text
Rscript run_stage2d1_correctness.R <repo> <output> --smoke
```

The formal grid covers `n = 50, 500, 2000, 5000`, `nsim = 1, 2, 199, 999`,
multiple seeds, `ncores = 1` and `2` when available, and chunk boundaries
including `127, 128, 129, 255, 256, 257` plus `nsim - 1`, `nsim`, and
`nsim + 1`. It includes controlled null rows, ordered exceedance booleans,
P, MCSE, requested/successful/failed accounting, malformed permutations,
large offsets, constant and near-constant traits, two-tip trees, four NA mask
families, input immutability, and the test=FALSE matrix sizes 1/8/32/100.

Required output files include:

```text
stage2d1_correctness_status.csv
stage2d1_correctness_rows.csv
stage2d1_ordered_nulls.csv
stage2d1_failure_accounting.csv
stage2d1_na_mask_guard.csv
stage2d1_batch_guard.csv
stage2d1_correctness_provenance.csv
```

The status is `PASS` only when the prototype hooks are present and all core
comparison, replay, failure, NA, and batch guards pass. Missing hooks are
written as `NOT_RUN` and terminate the script; they are never interpreted as
success.

## Formal timing invocation

Only after the same output directory contains
`stage2d1_correctness_status.csv` with `status = PASS`:

```text
Rscript run_stage2d1_benchmark.R <repo> <output> --formal
```

The benchmark is serialized. It uses fixed prepared fixtures, one warmup per
route/cell, alternating old/prototype order, the same seed for both calls,
median and IQR summaries, and the prepared grid `n = 500, 2000, 5000, 10000`
with `nsim = 199, 999, 9999`. Formal repeats default to 10 and can be changed
only for an explicitly documented run with
`FASTPHYLOSIG_STAGE2D1_BENCH_REPEATS`.

Correctness always compares the installed production kernel with the selected
prototype mode and retains the ordered null sequence. Timing disables audit
counters and ordered-null storage. A separate production-versus-private-oracle
calibration found 7.1-9.7% generic prototype dispatch overhead in representative
heavy cells. Candidate A/B timing therefore sets
`FASTPHYLOSIG_STAGE2D1_BENCH_BASELINE=private_oracle`: both sides then use the
same compiled entry point and differ only in the candidate allocation/copy
policy. The calibration remains separate, and private-oracle ratios must not
be presented as installed-package speedups.

The raw public `fast_k()` route is retained as an end-to-end guard and is not
called a candidate speedup. The test=FALSE matrix path at n=500 with
1/8/32/100 traits is recorded as an unchanged production protection guard;
Stage 2D1 does not replace that path.

Formal outputs include prepared old/prototype timings, raw route guard
timings, test=FALSE batch guard timings, memory proxies, paired medians/IQR,
and the acceptance status. A missing comparable candidate row is a failed
performance gate, not a zero or an extrapolated speedup.

## Acceptance gates

The candidate performance gate requires all of the following:

- controlled and RNG correctness status is PASS;
- at least two heavy prepared cells with `nsim >= 999` show median reduction
  of at least 15%;
- at least one heavy cell has speedup at least 1.20x;
- no light prepared cell shows confirmed slowdown greater than 5%;
- the test=FALSE/batch protection guard has no warning or runtime failure;
- candidate workspace remains bounded and is reported separately from the
  caller-supplied controlled matrix and optional returned null matrix.

If any gate fails, `K_PERMUTATION_PROTOTYPE` is `FAIL` and no production
integration is permitted. If the prototype is not yet wired, the status is
`NOT_RUN` at the correctness stage and formal timing is blocked.

For the completed run, short prepared cells near the Windows timer resolution
triggered the conservative light-cell flag. This does not rescue a candidate:
the heavy performance gate is independently required and failed for A and B.

## Stage 2D1 stop condition

After correctness and, if explicitly requested, formal timing, stop. Do not
modify production bindings, integrate the prototype, start D optimization,
redesign RNG, alter public APIs, or begin another candidate stage.

At completion, the expert review should report the frozen RNG contract,
allocation/copy audit, Candidate A/B/C decisions, ordered-null and accounting
parity, memory proxies, prepared/raw timings, light and batch regression
guards, and the minimum production candidate set. A `PASS` in this directory
is not a release or production-integration commit.

## Private batch `test=TRUE` timing guard

The focused batch guard is separate from the broad timing script and is
explicitly opt-in:

```text
Rscript run_stage2d1_batch_benchmark.R <repo> <output> --formal
```

It requires the same output directory to contain a correctness status of
`PASS`. The formal workload is fixed to a balanced tree with `n = 500`,
`nsim = 199`, `ncores = 1`, `trait_chunk = 64`, `simulation_chunk = 128`,
and 1/8/32/100 trait columns. The checked-in prototype is called in private
`oracle`, `A`, and `B` modes. Oracle-versus-A and oracle-versus-B exact parity
is checked with the ordered null retained before timing. Each candidate has
one warmup and ten same-seed paired repeats with alternating call order.

The script writes machine-readable exact-parity, raw timing, paired
median/IQR, summary, memory-proxy, status, and provenance CSV files named with
the `stage2d1_batch_` prefix. It does not alter production code or authorize
candidate integration. `--smoke` is a reduced wiring check and is not formal
evidence.
