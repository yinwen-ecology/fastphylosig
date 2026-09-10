# Stage 2D2B performance protocol

This directory contains audit-only timing infrastructure for the one final
Candidate B remediation. It does not change `R/`, `src/`, package tests, or
the public API. The runner is deliberately fail-closed: it will not call an
oracle or candidate until the separately produced Stage 2D2B exactness status
file contains `status=PASS`.

## Required candidate contract

Set `FASTPHYLOSIG_STAGE2D2B_PROTO_CPP` to one sourceCpp translation unit that
contains both a frozen production-equivalent oracle and the true blocked
candidate. The translation unit must expose either:

* one generic export, such as `stage2d2b_k_null_prototype(...)`, accepting
  `mode="oracle"` and `mode="B"`; or
* two sourceCpp exports, one named `stage2d2b_k_null_oracle` (or
  `stage2d2b_k_permutation_oracle`) and one named
  `stage2d2b_k_null_candidate_b` (or
  `stage2d2b_k_permutation_candidate_b`).

The runner rejects pairs that are not sourceCpp-bound from the same
translation unit. Candidate results must report the audit counters in the
returned list, preferably under `candidate_metadata`:

* `candidate_old_compute_one_calls` (must be exactly zero);
* `oracle_call_count` (must be exactly zero on the candidate path); and
* `candidate_compute_count` (positive for a successful non-empty evaluation).

`blocks_processed`, `workspace_peak_bytes`, and equivalent explicit memory
metadata should also be returned. Missing hard-counter metadata is a failed
candidate provenance gate, not evidence of correctness.

## Exactness prerequisite

Run the independent bounded exactness harness first and set
`FASTPHYLOSIG_STAGE2D2B_EXACTNESS_STATUS` to its status CSV. The default is
`benchmarks/stage2d2b/results/exactness/stage2d2b_exactness_status.csv`.
The status must be `PASS`; a missing file, a non-PASS status, a non-zero old
`compute_one` counter, or a source hash mismatch produces `NOT_RUN` and no
timings. The exactness harness is responsible for controlled permutation,
identity-first accounting, RNG replay, and failure-semantics fixtures.

## Invocation order

All formal timing is serialized in one R process. Each cell uses one fixed
tree, trait matrix, and pre-generated controlled permutation matrix.

1. Evaluation-only gate:

   `Rscript run_stage2d2b_benchmark.R --formal <repo> <out>`

   This evaluates `n=2000,5000,10000`, `nsim=999,9999`, and
   `block_size=2,4,8,16` at one native thread. The default shape set is
   balanced/random/pectinate. It performs one warmup per side, then at least
   ten alternating oracle/candidate pairs. Generation is outside the timed
   interval. At least two heavy cells (`n >= 5000`, `nsim=9999`) must have
   speedup `>= 1.35` after selecting a block.

2. Full-pipeline gate, only after the evaluation status is `PASS`:

   `Rscript run_stage2d2b_benchmark.R --full <repo> <out>`

   This uses the frozen internal Fisher-Yates path with `n=2000,5000,10000`,
   `nsim=9999`, and the selected block. At least two representative heavy
   cells must improve by `>=20%`, with at least one speedup `>=1.25`.

3. Regression and parallel guards, only after the full status is `PASS`:

   `Rscript run_stage2d2b_benchmark.R --guards <repo> <out>`

   This records `nsim=199`, `test=FALSE`, batch widths 1/8/32/100, NA and
   numerical edge traits, constant/two-tip/retained-subset cases, and native
   thread counts 1/2/4/8. Guard results are not used to manufacture a heavy
   speedup. A confirmed representative slowdown over 5%, warning, error, or
   failed parity is a guard failure.

`--smoke` is a small, serialized protocol check with the same prerequisite
gate. Environment variables listed below can narrow a formal run for local
debugging; a formal run refuses fewer than ten paired repeats.

## Timing and memory rules

The runner uses warmup, alternating order, fixed seeds, and no R worker
processes. Controlled evaluation passes a pre-generated integer permutation
matrix and sets `return_sim=FALSE`; the matrix generation time is recorded
separately. Full-pipeline mode passes no matrix and restores the same R seed
around both calls. The candidate and oracle calls are never nested.

The candidate workspace must be bounded by `n * block_size` (or an explicitly
equivalent bounded layout), never `n * nsim`. The memory CSV records the
reported candidate metadata and the independent index upper bound
`n * block_size * sizeof(int)`. The supplied controlled permutation matrix and
returned result object are not counted as candidate working memory.

Block selection minimizes the median candidate time over heavy evaluation
cells. If blocks 8 and 16 differ by less than 3%, block 8 is selected even if
16 is marginally faster. Selection is made only after valid timing and hard
counter checks.

## Outputs

Each invocation writes only audit evidence under its requested output
directory:

* `stage2d2b_benchmark_status.csv`;
* `stage2d2b_benchmark_raw.csv`;
* `stage2d2b_benchmark_summary.csv`;
* `stage2d2b_benchmark_memory.csv`; and
* `stage2d2b_benchmark_provenance.csv`.

The runner does not write package artifacts and does not delete historical
evidence. A performance failure is recorded as `FAIL`/`REJECTED`; full and
guard modes cannot run from a failed prerequisite.

The formal grid is intentionally not run by merely sourcing this file. It is
opt-in through the command line and should be run only after the exactness
status is final.
