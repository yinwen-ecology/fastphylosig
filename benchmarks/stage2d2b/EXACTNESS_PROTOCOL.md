# Stage 2D2B Exactness Protocol

This directory contains audit-only infrastructure for the one authorized
Candidate B remediation.  The correctness runner does not modify `R/`, `src/`,
package tests, the public API, or frozen Stage 2D2 evidence.  It is fail-closed:
the oracle and candidate must be sourceCpp-bound exports from one translation
unit.

## Required environment

Set these variables before invoking R:

```text
FASTPHYLOSIG_STAGE2D2B_PROTO_CPP=<one-TU oracle and Candidate B source file>
FASTPHYLOSIG_STAGE2D2B_GATE=AUTHORIZED
FASTPHYLOSIG_STAGE2D2B_LIBRARY=<audited fastphylosig library, optional>
```

The translation unit must expose either a generic export accepting
`mode="oracle"` and `mode="B"`, or separate sourceCpp exports named
`stage2d2b_k_permutation_oracle` (or the historical Stage 2D2 oracle name)
and `stage2d2_k_permutation_candidate_b` (or an equivalent Candidate B name).

The candidate result must report these hard provenance counters, at top level
or under `candidate_metadata`/`counters`:

* `candidate_old_compute_one_calls` exactly `0`;
* `oracle_call_count` exactly `0`; and
* `candidate_compute_count` strictly positive for a non-empty successful call.

Missing counters are a failed provenance gate, not evidence of correctness.

## Invocation

Run from an ASCII staging tree when the repository path contains non-ASCII
characters:

```powershell
Rscript benchmarks/stage2d2b/run_stage2d2b_correctness.R --smoke <repo> <out>
Rscript benchmarks/stage2d2b/run_stage2d2b_correctness.R --correctness <repo> <out>
```

Smoke uses n=50, three tree shapes, three controlled patterns, blocks 2/4,
and nsim=7.  Formal mode uses balanced/random/pectinate trees, n=50/500/2000/
5000, all six controlled patterns, blocks 2/4/8/16, nsim=19, and replay at
native thread counts 1/2.  Environment variables `FASTPHYLOSIG_STAGE2D2B_N`,
`_SHAPES`, `_PATTERNS`, `_BLOCKS`, `_TRAITS`, `_NSIM`, `_REPLAY_NSIM`, and
`_REPLAY_THREADS` may narrow a local diagnostic run; they must not be used to
represent the formal evidence unless the provenance records the narrowed grid.

## Gate order

1. Controlled permutations run first.  Each identity, single-swap, reverse,
   cyclic-shift, fixed-random, and adversarial row is compared in order against
   the same-TU frozen oracle.  K, ordered null K, status, failure reason,
   exceedance count, P, MCSE, requested/successful/failed counts, and all
   identity metadata are compared bitwise where numeric representation allows.
2. Identity-first internal-RNG calls are compared against the oracle and then
   replayed under the same restored seed.  The first row must equal observed K;
   the accounting denominator must distinguish `nsim` from randomizations.
3. Malformed permutation dimensions, duplicate/out-of-range indices, nonfinite
   traits, and invalid block size must retain the expected error/status and
   warning semantics.
4. Candidate inputs are snapshotted before evaluation.  Tree structure,
   branch lengths, labels, `Nnode`, trait matrix, and supplied permutations must
   remain identical after evaluation.
5. A controlled batch guard covers trait widths 1/8/32/100 without adding RNG
   or timing claims.

The candidate ordered null matrix must be `nsim x traits`; no candidate result
may be copied from the oracle.  Any exceedance, P, or MCSE mismatch fails the
run even when the K vector happens to match.

## Evidence files

The runner writes only to `<out>`:

* `stage2d2b_correctness_status.csv`;
* `stage2d2b_correctness_summary.csv` and `_replicates.csv`;
* `stage2d2b_correctness_counters.csv`;
* `stage2d2b_correctness_failures.csv`;
* `stage2d2b_correctness_immutability.csv`;
* `stage2d2b_batch_guard.csv`;
* `stage2d2b_rng_replay.csv`; and
* `stage2d2b_correctness_provenance.csv`.

Formal acceptance requires status `PASS`, zero controlled failures, bitwise
replicate parity, passing accounting and replay gates, passing failure and
immutability gates, all hard counters valid, and `production_code_changed=NO`.
Performance is a separate protocol and cannot rescue a correctness failure.
