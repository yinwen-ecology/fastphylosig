# Preparation Engine Optimization, Stage 1

This directory preserves the fixed before evidence and the paired benchmark
protocol for the fastphylosig 0.2.0 preparation-engine work. The package build
excludes `benchmarks/`; these files are development evidence, not installed
package data.

## Baseline provenance

- Source commit: `278c01a7d2e588ff4f82581c0efd0b4a16465c73`
- Source version at measurement: `0.1.0`
- R: 4.6.1 on Windows
- Compiler: GCC 14.3.0 with OpenMP
- Protocol: warmup followed by 20 formal repeats for every K grid and batch
  cell; medians are stored in the CSV files under `before/`.

The 0.2.0 Stage 1 candidate changes must preserve fixed fixtures and use
alternating before/after order with at least 10 paired repeats (20 preferred).

## Before evidence

- `k_grid_summary.csv`: K raw/prepared/unified paths for n=50--5000.
- `k_batch_summary.csv`: matrix and column-loop paths at n=500.
- `k_exclusive_phase_budget.csv`: low-overhead K component budget.
- `k_prepared_duplicate_call_counts.csv`: validation/helper call counts.
- `lambda_profile_summary.csv` and `lambda_profile_detail.csv`.
- `d_profile_summary.csv` and `d_instrumented_split.csv`.
- `delta_profile_summary_complete.csv`.

No statistical estimator, numerical formula, tolerance, RNG policy, thread
policy, or public API was changed while collecting this evidence.

## Candidate gates

- Candidate 1, iterative canonical signature: 25 targeted assertions and
  6,825 full-suite assertions passed. Recursive and iterative signatures were
  identical for balanced, pectinate, random, polytomous, edge-shuffled, and
  safely renumbered fixtures. Pectinate raw lambda at n=1,000 changed from an
  expression-depth error to a successful fit. The K n=5,000 paired result was
  0.98x, so this candidate is retained only as the required P0 correctness fix.
- Candidate 2, single-boundary validation: 9 targeted assertions and 6,834
  full-suite assertions passed. A prepared public call performs one complete
  fingerprint; each independent call validates again. All protected-field
  mutation tests passed.
- Candidate 3, retained-subtree evidence: 25 targeted assertions and 6,859
  full-suite assertions passed. Inspection evidence is created once per cached
  retained mask, cache identity failures rebuild the entry, and the structural
  cache remains bounded.

## Authoritative Stage 1 paired result

`stage1_final_paired_timings.csv` and `stage1_final_paired_summary.csv` compare
the audited 0.1.0 source with the final Stage 1 source. Every workload used 10
alternating before/after pairs after symmetric warmup. Short workloads used
10--20 calls inside each timed repeat and report the per-call time; n=5,000
workloads used one call per repeat.

| Workload | Before median (ms) | After median (ms) | Speedup |
|---|---:|---:|---:|
| K prepared n=500 | 51.5 | 7.0 | 7.36x |
| K prepared n=5000 | 1780 | 50 | 35.60x |
| K raw n=5000 | 7820 | 6480 | 1.21x |
| K prepared matrix n=500, 100 traits | 60 | 27 | 2.22x |
| lambda prepared n=500 | 54.5 | 14.5 | 3.76x |

The maximum absolute statistic difference was zero in every paired workload.
The prepared n=500 helper audit changed `.tree_fingerprint` from 9 to 1 call,
`.inspect_tree_core()` from 2 to 0 calls, and `.prepared_tree_subset()` from 6
to 1 call per public analysis. The before and after counts are preserved in
`before/k_prepared_duplicate_call_counts.csv` and
`stage1_final_call_counts.csv`.
