# fastphylosig 0.2.0 Stage 2B2B Candidate 1

This directory contains a read-only paired benchmark for the Candidate 1
inspection rewrite. It compares the frozen old inspector copied from baseline
commit `4dc4fa3` with the current `.inspect_tree_core()` implementation. The
public `prepare_tree()` route is also compared by swapping the namespace
binding before the timer starts and restoring it after the timer ends.

## Protocol

The default grid is balanced, random, and pectinate fixtures at `n = 500,
1000, 2000, 5000, 10000, 20000`. Every fixture has deterministic tip labels,
unit positive branch lengths, and a fixed edge construction. Random fixtures
use seed `2072000 + n`. Each cell performs one unmeasured old/new warmup and
parity check, then serial paired timings on fresh serialized tree copies. Pair
order alternates `old > new`, `new > old`. The default repeat minimum is 10
paired repeats through `n = 5000` and 5 at `n >= 10000`.

Run from an ASCII staging copy of the repository and an ASCII output path:

```text
Rscript benchmarks/stage2b2b1/run_paired.R <ascii-repo> <ascii-output>
```

The script never edits package source or tests. It is intentionally serial;
run independent shape or size partitions as separate R processes when needed.

## Bounded partitions

Environment variables use the `FASTPHYLOSIG_STAGE2B2B1_` prefix:

- `SHAPES=balanced,random,pectinate`
- `N_GRID=500,1000,5000`
- `REPEATS=10`
- `REPEATS_BY_N=500:10,10000:5,20000:5`
- `PARTITION=pectinate-large` (a provenance label)
- `LIB=<installed library>` (avoids source loading when needed)
- `TIMING_TIMEOUT=<seconds>` (records censored calls instead of hanging)
- `QUICK=1` (small smoke defaults: `n = 20, 50, 100`, two pairs; it relaxes
  the formal repeat minimum for smoke validation)

`SHAPES`, `N_GRID`, and `REPEATS_BY_N` are intended for serialized partitions.
Formal mode rejects repeat counts below 10 for `n <= 5000` or below 5 for
larger configured sizes.

## Outputs

- `paired_timings.csv`: one row per shape, size, workload, and paired repeat,
  including order, status, elapsed time, warning/error fields, and new/old
  slowdown.
- `paired_summary.csv`: median and IQR for old and new timings, speedup,
  slowdown, and successful/failed pair counts.
- `empirical_p_before_after.csv`: adjacent-size empirical exponents for the
  old (`before`) and new (`after`) paths.
- `representative_slowdown_gate.csv`: cell and all-shape-median checks at
  representative `n = 5000` (or the largest configured size at or below 5000)
  against a strict 5 percent slowdown threshold.
- `operation_counts.csv`: deterministic source-aligned operation counts. Old
  character-list lookup and key-materialization counts are separated from
  copied existing stack slots for pop and append operations. New integer
  adjacency fill, cursor, offset, child/edge-row reads, and preallocated-stack
  writes are recorded; new stack resize copied slots are zero.
- `fixture_manifest.csv`: shape, size, root, edge count, and compact checksum.
- `provenance.csv`: exact candidate/base commit IDs, source and oracle MD5s,
  R/compiler/platform details, package versions, paths, protocol, and grid.

The accepted formal evidence is retained in `results/formal/`. The complete
test, install, and package-check transcripts are `results/testthat.Rout`,
`results/00install.out`, and `results/00check.log`.

## Small-cell timer resolution

`run_small_resolution.R` resolves sub-20-ms cells by timing 50 inspection
calls per block over 24 alternating old/new paired blocks. Binding changes,
fixture cloning, and exact-equality checks remain outside the timed block. Its
accepted output is retained in `results/small_resolution/`.
