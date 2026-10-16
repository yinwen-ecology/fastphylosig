# fast_hetero D/lambda performance candidate

Version remains 0.3.3.9021. Candidate SHA256 0A97F8379DED8F6C194BFEC3790019E01833731F9BDAD6746BB43B0EC0E8E016; baseline SHA256 69FD3343F36254F5BFDE4C3E9A70514F221E70CF6B33C475E149D715640392D9. Exactly two production files differ; K code, R APIs/defaults, node-bound and ARD correctness fixes, tests and fixtures are unchanged.

## Changes and scientific constraints

D: reuse long-double child transition mixtures in count_weights, retaining the original coefficient i/j convolution/accumulation order. The 65-point q grid/refinement, truncation, likelihood scaling, conditional null model and all draws are unchanged. Two reusable O(m) long-double scratch arrays add bounded workspace; no memory-saving claim is made.

Lambda: the internal Brent implementation optionally returns fx paired with its returned best x. Bootstrap validity uses that fx instead of re-evaluating the same objective. The shared four-argument C++ entry point remains available. This is NOT reuse of the last-evaluated point. Search iterations, endpoints, tolerances, objective/REML definitions and failure handling are unchanged.

## Controlled measurements

Same Windows machine, R 4.6.1, GCC 14.3.0, existing dependencies, single OpenMP thread. 1024 tips, 1 or 6 traits, automatic clades 24–128 tips, shared NA patterns, K/lambda/D test FALSE/TRUE; nsim=19 and d_nsim=49 in BOTH versions. No default or simulation count was reduced for speed.

Data read/tree creation/prepare_tree, seed setting and (in the later lambda/combined comparisons) GC are outside timing. Three rounds alternate baseline/candidate process order; each case/path records a first invocation separately, then three reused invocations, giving 9 paired reused measurements per case/path. First invocation means first call of that scenario within its batch, not always a completely cold R/DLL process. Both raw and prepared paths are measured. Baseline CPU profiling additionally used a fresh R process per case.

No full package check/compilation ran concurrently with the CPU timing comparisons. Unrelated existing user analyses were not stopped or altered. Short millisecond-level cases are near timer resolution and are not treated as meaningful speedups. Results are bounded workload evidence, not a universal speed claim.

D-only: 108 complete statistical result objects and final RNG states bitwise identical. Lambda-only versus D candidate: 108 identical. Original baseline versus combined candidate: **324 identical**, including every clade/heterogeneity/overall/selection/null result, diagnostics/status fields and final .Random.seed. Only runtime timing and call/provenance fields were removed from comparison; no statistical field or numeric tolerance was relaxed. Inputs are retained in inputs.rds.

| case | path | baseline median s | candidate median s | paired speedup median | paired IQR |
|---|---|---:|---:|---:|---|
| D-TRUE-multi | prepared | 0.450 | 0.340 | 1.324 | 1.294–1.333 |
| D-TRUE-multi | raw | 0.430 | 0.330 | 1.333 | 1.303–1.344 |
| D-TRUE-NA | prepared | 0.410 | 0.310 | 1.323 | 1.323–1.419 |
| D-TRUE-NA | raw | 0.410 | 0.300 | 1.345 | 1.300–1.367 |
| lambda-TRUE-multi | prepared | 0.670 | 0.640 | 1.031 | 1.031–1.062 |
| lambda-TRUE-multi | raw | 0.660 | 0.630 | 1.048 | 1.016–1.048 |
| lambda-TRUE-NA | prepared | 0.620 | 0.590 | 1.052 | 1.051–1.068 |
| lambda-TRUE-NA | raw | 0.620 | 0.590 | 1.051 | 1.034–1.069 |

The independent first D comparison predated GC outside every timed invocation; its wider timing spread is retained and distinguished from the final combined comparison. Descriptive/single-trait paths do not show stable general gains. Read all per-case distributions in the CSV files, including controls.

## CPU, allocations and memory

Baseline sampled CPU: lambda multi-trait calibration about 93%; D calibration about 92%. K NA scanning includes about 45% in ape pruning. These are R-visible native boundaries, not claimed line-level native samples. Full-function before/after results corroborate the selected arithmetic reuse.

Rprofmem reports R allocations, not std::vector allocations. The separate memory comparison uses the same fixed number of calls in fresh R processes for both versions; actual R PIDs (not the launcher) are sampled via Windows GetProcessMemoryInfo. Peak working set/commit includes setup, namespace and profiler overhead and remains observational under normal system pressure. MEMORY_COMPARISON.csv retains all values; added D scratch memory is explicit.

## K decision

K factors were already outside the block loop. At identical 129 simulations, chunk sizes 8/32/129 gave 0.11/0.11/0.10 seconds in a small diagnostic and identical scientific results/RNG. This is not a statistically established block-reuse benefit and does not isolate zeroing/allocator CPU. NA pruning and R result assembly are visible costs; changing the biological pruning path is higher risk. **K remains byte-identical**, with no K speedup claimed.

## Reproduction and validation status

Scripts paired_batch.R, compare_pairs.R and run_pairs.py are supplied; set R_SCRIPT to a platform-appropriate Rscript executable. For Linux unpack this evidence at the project root with docs/hetero-performance-20261016/, compile baseline-source.tar.gz into a separate library, apply D-only.patch / lambda-only.patch separately to reproduce stages, then compile the frozen candidate archive into another library. Use the same inputs.rds and pair counts. profile_runner.py contains Windows-specific process-memory instrumentation; profile_case.R itself provides portable CPU/R allocation profiling.

The full Windows candidate check has been launched with --no-manual --run-donttest (full tests enabled) and records its independent completion status under final-windows-full-check/. At this evidence handoff it may still be running; no completed candidate check is claimed until that state/log is inspected. PDF manual is not checked. No LaTeX/system software installation, default-library update, formal lambda benchmark, publication or CRAN action occurred.

## Windows completion update
The frozen candidate full Windows check completed with Status OK (0 ERROR, 0 WARNING, 0 NOTE). Testthat: 13679 PASS, 0 FAIL, 0 WARN, 1 SKIP (cross-install RNG reference libraries absent). The standalone external bridge separately skipped its build-excluded harness. PDF remains unchecked. See final-windows-full-check/WINDOWS_FINAL_REPORT.md and raw logs.
