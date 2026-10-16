# Development validation: 2026-10-13

## Subsequent scientific revision (0.3.2.9000)

The [four-task revision](four-goals-20261013/README_zh.md) corrects the two
Delta issues identified below. Nine independent posterior-integral cases and
the full corrected-package test suite pass. All four tasks are now complete:
7000 calibration datasets, exact corrected-target acceleration parity,
cutting guidance and executable tree plots. The final 0.3.2.9000 source check
reports Status: OK, with 10025 assertions passing and one existing skip.
The workspace library was updated with its previous version backed up.
Numerical/RNG equivalence records
below describe the historical 0.3.1.9000 target; they do not apply across the
intentional scientific correction. Heterogeneity calibration is recorded
separately in the new revision directory.

## Previous optimization record (0.3.1.9000)

Package version: **0.3.1.9000**. This is a local development addition; no
GitHub or CRAN release was published. Existing release archives were preserved.
The verified Windows/R 4.6.1 installation was also copied into the workspace's
`.rlib_heterogeneity` library and loaded successfully there, including a BM
calibration smoke call. The user's default R library was not changed.

## Scientific readiness audit: important Delta qualifications

The [2026-10-13 scientific audit](scientific-readiness-20261013/README_zh.md)
identifies two pre-existing Delta issues that acceleration parity does not
resolve: the package reports mean(beta)/mean(alpha), whereas the author's R
implementation reports mean(beta/alpha); and its log-normal MH proposal lacks
the proposal-density ratio required for the documented exponential-prior
posterior. Same-chain comparisons and an independent analytic detailed-balance
check are archived with that audit. No production estimator was changed in
this assessment. Existing numerical equivalence records must not be interpreted
as validation of the claimed posterior or comprehensive heterogeneity-test
calibration. These method issues should be resolved before such claims.

## Current mixed-length ARD performance audit

The [mixed-length ARD audit](ard-mixed-performance-20261013/README_zh.md)
records **10,005 passed assertions**, zero failures/errors/test warnings and
one existing external-version skip. **125 strict paired comparisons** retain
results, full optimization trajectories, Delta chains, RNG state and following
draws. Three mixed-length default fits cover **6,757 objective calls**. Checks
include 2/4/8/20 states, random/balanced/ladder trees, NA tips, branch scales,
six RNG combinations, reference checks and two-process PSOCK execution.

The existing exact exponent cache now combines with a compact immutable
worklist of uncached branches. Large mixed workloads batch only that work;
small/fully cached workloads remain serial. Odd batch widths avoid starting
another parallel region for one leftover branch. The **2 MiB** batch budget
and automatic maximum of **two threads** are retained. Pruning, normalization,
matrix products, optimizer calls and RNG order remain unchanged.

Three sequential medians give **6.37 -> 4.70 s** (1.36x) for a 4,096-tip,
four-state ARD fit with 1% rounded branches, **3.71 -> 3.00 s** (1.24x) with
50% rounded branches, and **15.44 -> 13.75 s** (1.12x) for a 50% mixed Delta
two-trait scan plus overall estimates. Rounded fraction is not the exact cache
fraction. Fully distinct/fully repeated cases are essentially unchanged;
CPU and wall time are reported separately because parallelism costs CPU.

A 32,768-tip isolated native experiment gives **96 runs / 2,880 strictly
identical objective values and RNG states** across batch budgets and expression
prototypes. Fourfold larger batches yield only about 2%-6% wall-time gains;
direct writes show small/inconsistent gains. Both remain diagnostic experiments.

Final clean source compilation, linking and examples return **Status: OK**;
tests ran separately against that checked installation in `.rlib_heterogeneity`.
**133 files** match workspace, staging and archive contents, allowing
build-inserted DESCRIPTION metadata and CRLF-to-LF normalization. Current source
archive SHA-256:
`5f68cb851a89a4558911abe292ff279e98aa7704b674c7149c672d3825752648`.

## Previous distinct-length ARD performance audit

The [distinct-length ARD audit](ard-unique-performance-20261013/README_zh.md)
records **9,883 passed assertions**, zero failures/errors/test warnings and one
existing external-version skip. **119 strict paired comparisons** retain
results, full optimizer trajectories, MCMC chains, RNG states and following
draws. The optimizer comparisons cover **13,370 objective calls** across
two/four/eight states. Additional checks cover large exponent batches,
random/balanced/ladder trees, 2/4/8/20 states, branch scales, missing tips,
six RNG combinations and two-process Delta PSOCK execution.

The existing ACE complex exponent expression now runs in bounded **2 MiB**
OpenMP batches for large distinct-length likelihood-only fits. Serial pruning,
normalization, matrix products and optimizer order are retained. The default
uses at most **two threads**; `options(fastphylosig.ace_threads = 1L)` selects
serial fits. Small crowns and Delta PSOCK inner fits remain serial.

Three sequential repetitions give median **5.14 -> 3.46 s** (1.49x) for a
2,048-tip/four-state ARD fit and **12.22 -> 9.68 s** (1.26x) for Delta local
scans plus overall estimates on two traits. Parallel execution can increase
total CPU use; the audit reports CPU and wall time separately. Small
256-tip/eight-state fits give **3.18 -> 3.18 s**. No optimizer or MCMC settings
were shortened. Exactness is verified on Windows/R 4.6.1, rather than assuming
cross-platform libm bit equality.

Clean source compilation, linking and examples return **Status: OK**; tests
ran separately against the workspace installation. **132 source files** match
workspace, staging and archive contents, allowing build-inserted DESCRIPTION
metadata and CRLF-to-LF normalization. Source archive SHA-256 at that audit:
`e0bff86fea1a65448f921bfc29b3edf412d61516a15af28a24ebed96993f100c`.
The checked installation is copied to `.rlib_heterogeneity`; final paired
optimizer and Delta/RNG checks are recorded in the audit directory.

## Previous Delta performance audit

The [Delta performance audit](delta-performance-20261013/README_zh.md) records
**9,759 passed assertions**, zero failures/errors/test warnings and one existing
cross-install environment skip. **160 strict paired comparisons** preserve
complete V2 graphs, ER/ARD fits, Delta signal/reference tests, ancestral
probabilities, MCMC traces, RNG state and subsequent random draws. Comparisons
include three normal RNG algorithms and two uniform RNG kinds.

Native canonicalization gives median **0.60 -> 0.03 s** for 8,192 tips. ARD fits
on equal/rounded branch lengths reuse exact exponent vectors and give
**2.64 -> 0.76 s** and **4.57 -> 1.47 s**, respectively, for 2,048 tips/four
states. Random distinct lengths give a smaller fitting gain. Delta local scans
plus overall estimates on 2,048 tips/two traits give **1.24 -> 0.80 s** for ER
and **14.51 -> 11.97 s** for ARD. Each is a three-repetition median. Defaults,
optimizer settings, floating-point expression order and draw order are retained.
Native fit scratch is explicitly released at exit.
Two-chain MCMC with five million iterations per chain gives **7.39 -> 5.86 s**,
with exact summaries, RNG state and following draws.

Compilation, linking and examples return **Status: OK**; a final archive check
against that installed production code also returns **Status: OK**. Tests ran
separately against the workspace installation. The source manifest verifies
**131 files** against staging and archive contents, allowing R's documented
DESCRIPTION metadata insertion and CRLF-to-LF build normalization. Current
archive SHA-256 at that audit:
`10ae00f78698c6393db87cfabe3ea238ad8419d93ca2bb730f9c7849fa89d727`.

## Previous deep performance audit

The [deep performance audit](deep-performance-20261013/README_zh.md) records
**9,522 passed assertions**, zero failures/errors/test warnings and one existing
environment skip. All **167 paired comparisons are strictly identical**,
including four-method fits/probabilities/RNG and 144 ER/ARD ancestral probability
scenarios with varied state counts, branch scales, shapes and missing tips.
K local tests reuse compiled scan groups and the specialist C++ permutation
kernel. ACE shares native topology validation and defers independent row
normalizations while preserving floating-point fixed points/two-cycles.
Automatic selection skips whole selected/undersized traversal intervals.

Three sequential repetitions give median **7.08 -> 0.36 s** for 32,768 tips,
16 traits and K scan + 19 local tests; **9.45 -> 5.57 s** for 8,192 tips, two
categorical traits and ER Delta scan + overall estimate (default MCMC, no tests).
Fixed-rate ACE ancestral reconstruction on a 32,768-tip random tree gives
**5.89 -> 0.06 s**; this is a kernel gain, not an end-to-end Delta gain.
The directory includes cutting coverage comparisons and recommended bounds.
Source compilation/linking/examples return **Status: OK** (tests run separately).
All 130 workspace source/documentation/test files match the packaged staging
source by SHA-256. That workspace install produced the paired results.
The archive SHA-256 at that audit was:
`102111f5534040d95c6d11d9d56f1a8e801fde57ec4fe71df8050a91e007656e`.

## Earlier correctness/performance audit

The [earlier audit](audit-20261013/README.md) records **9,468 passed assertions**,
zero failures/errors/test warnings and one existing skip. It fixes local
K/lambda eligibility and unifies Delta reference MCSE with the shared native
plus-one summary. Lambda local LRTs reuse the first fit; Delta tests reuse scan
contexts and crown-sized permutation data; Mk simulation pools exact repeated
branch lengths. Valid estimates/probabilities, retained draws and final RNG
states agree with the prior development implementation.

Three paired sequential repetitions gave median 3.25 -> 0.81 s for a 12K-tip,
four-trait lambda scan with local LRTs, and 27.17 -> 23.56 s for an 8K-tip,
three-trait ER Delta scan with five local permutations (MCMC 1200/5/100).
Repeated-length Mk preparation benefits strongly; random lengths do not.
The scripts, workload limits, parity and timing tables are in the audit directory.
Clean source compilation, symbols and examples return **Status: OK**. The
updated workspace library reproduces the refreshed source demonstration in a
fresh process. Default-MCMC (10,000/10/100) old/new installed comparisons also
retain estimates, test probabilities, null draws and RNG, excluding the
intentional reference MCSE correction.

## Earlier Delta integration

The current source and workspace library support K/lambda/D/Delta scans.
Local/overall estimates, categorical composition, descriptive contrasts and
tree plots are primary; reference-model checks are optional supplements.

- Complete freshly compiled source suite: **9,428 passed assertions**, zero
  failures/errors/test warnings, one existing external cross-version RNG skip.
- Final Delta-only regression after original-state reference labels were added:
  **102 passed assertions**, zero failures/errors/test warnings.
- Native ER simulation matches analytic R transitions for 3/4/6 states; ARD
  matches an independent uniformization oracle. Whole-tree Delta bootstrap
  draws match R simulations followed by public specialist refits. Random-tree
  canonicalization and interleaved NA groups retain estimator/final-RNG parity.
- Clean final source check (compilation, symbols and examples, tests separately
  as above): **Status: OK**. A fresh workspace-library process reproduces the
  source demonstration's clade/overall tables, reference summaries and all
  retained null draws exactly; plotting leaves its RNG unchanged.
- The new optional Delta local/overall reference check is approximate under a
  common fitted ER/ARD Mk model, not a rate-shift test or a comprehensively
  calibrated universal equal-signal test. Broader false-positive-rate and power
  validation remains outstanding and is explicitly documented in the API.

See [implementation and records](delta/README.md), `delta/full-tests.log`,
`delta/full-test-results.csv`, `delta/final-delta-tests.log`,
`delta/hetero-delta-test-results.csv`, `delta/source-check.log`,
`delta/installed-smoke.log`, and the [simulated tree illustration](delta/tree-delta.png).

## Earlier lambda/D numerical and regression checks

- Full test suite against the current lambda/D extended, freshly compiled installation:
  **9,323 passed assertions, 0 failures, 0 errors, 0 test warnings, 1 skipped test**.
  The skip is the existing cross-version 0.2.0 RNG check, which needs two
  externally supplied package libraries.
- The new tests compare subtree K with an independent dense covariance
  calculation, lambda with the existing specialist, BM simulations with
  path-summed branch increments and covariance, and simulated dispersion
  probabilities with a dense Monte Carlo oracle.
- Trait-specific missingness, constant/non-finite traits, original node IDs,
  explicit nested nodes, work limits, prepared contexts, scale invariance,
  print/plot behavior, and RNG consistency across simulation chunks are covered.
- An equal-sized, balanced-tree BM calibration used 100 independently
  generated trait datasets and 199 reference replicates per dataset.
  **5/100** rejected at P <= 0.05. This is a bounded calibration example,
  not qualification across every tree shape or evolutionary process.

See `test-totals.csv`, `test-results.csv`, `bm-calibration.csv`, and
`calibration.log` for the initial implementation results. The first acceleration
suite is recorded in `acceleration/test-totals.csv`, `acceleration/test-results.csv`
and `acceleration/full-tests.log`. The latest C++ implementation is recorded in
`cpp-acceleration/test-totals.csv`, `cpp-acceleration/test-results.csv` and
`cpp-acceleration/full-tests.log`. The current short API/tree-view build is
recorded in `tree-view/test-totals.csv`, `tree-view/test-results.csv` and
`tree-view/full-tests.log`.

## Initial implementation smoke timings (before acceleration)

Windows, R 4.6.1; one named continuous trait; automatic disjoint selection with
min_tips = 20 and max_tips = 200. These are single-run end-to-end timings,
including preparation, with no universal performance guarantee.

| Tree tips | Clades | Tree coverage | BM reference | Elapsed |
|---:|---:|---:|---:|---:|
| 2,000 | 17 | 98.05% | No | 0.61 s |
| 2,000 | 17 | 98.05% | 199 replicates | 0.84 s |
| 50,000 | 449 | 98.82% | No | 24.85 s |

This initial 50K workload validated the descriptive K scan only. Lambda and BM
calibration at 50K were not included in that run. See `smoke-timings.csv` and
`sessionInfo.txt`.

## Acceleration: matched end-to-end benchmarks

Each scenario used the same input, seeds, clade bounds (20–200 tips), and
three measured runs in the same Windows/R 4.6.1 environment. Timings include
tree validation, clade selection, data matching, estimation and calibration.
The table reports medians; speedups depend on tree shape, trait count,
missingness, number of replicates and the machine. These are serial runs.

| Scenario | Before | After | Speedup |
|---|---:|---:|---:|
| K, 50,000 tips, one trait, 449 clades | 14.64 s | 1.56 s | 9.38× |
| K + 199 BM replicates, same 50K input | 18.82 s | 2.83 s | 6.65× |
| K + 199 BM replicates, 2,000 tips, eight traits, three missingness patterns | 1.72 s | 0.51 s | 3.37× |
| Lambda, 2,000 tips, four traits | 0.62 s | 0.12 s | 5.17× |

All **12 matched runs** produced identical clade tables, heterogeneity
summaries, selections, retained BM dispersion values and final `.Random.seed`.
Timing and call fields were excluded from this comparison. Estimated K/lambda
and BM-null maximum absolute differences were zero on these workloads.
Additional tests compare the direct clade path with the public specialists
for missingness, large trait offsets, constant/Inf/NaN traits, zero branches
and protected prepared-context mutation. Invalid node metadata is checked
before compilation; method-specific branch restrictions remain local to crowns.

The acceleration removes unnecessary whole-tree canonical contexts and
repeated public-wrapper preparation per crown. It uses integer tip indexing,
shares compiled subtrees among traits with identical NA masks, and reuses them
for calibration. BM draws and branch scaling remain in R; C++ performs the
forward traversal with one scratch vector rather than a node-by-replicate
value matrix. The K estimator, lambda optimizer, sample dispersion, BM reference
and plus-one probability are unchanged. The public interface is unchanged.

Evidence: `acceleration/profile-summary.txt`, `acceleration/*-timings.csv`,
`acceleration/speedups.csv`, `acceleration/parity.csv` and
`acceleration/source-check.log` (fresh source check: **Status: OK**).

Reproduce from the workspace, using the preserved baseline installation and
the accelerated installation (neither changes the default user library):

```powershell
Rscript.exe docs/signal-heterogeneity/acceleration/benchmark.R .rlib_heterogeneity_baseline baseline
Rscript.exe docs/signal-heterogeneity/acceleration/benchmark.R .rlib_heterogeneity accelerated
Rscript.exe docs/signal-heterogeneity/acceleration/compare.R
Rscript.exe docs/signal-heterogeneity/acceleration/test.R
```

Load the accelerated version in a fresh R session:

```r
library(fastphylosig, lib.loc = "C:/Users/wenyi/Documents/ChatGPT/系统发育信号R包/.rlib_heterogeneity")
```

## C++ pipeline and experience from the five specialist functions

The scan now uses C++ for descendant indexing, disjoint selection, crown
extraction, packed exact NA-mask grouping, per-crown/trait batch orchestration,
BM generation, prepared K evaluation, online dispersion, quantiles and tail
counts. R retains the input/output interface, species matching, protected
prepared-context validation, observed per-trait SD and actionable diagnostics.
NA pruning uses ape's existing behavior through a C++ callback, preserving its
tip order, root collapse and unary-node handling. No new estimator was introduced.

The implementation applies concrete lessons from the existing sources:

| Source | Applied here |
|---|---|
| `src/k_tree.cpp`, `src/k_strategy.h` (`fast_k`) | The same stable two-pass K evaluator; tree factors computed once per crown simulation, reusable message/state buffers, and a 64 MiB block-buffer budget. |
| `src/lambda_tree.cpp` (`fast_lambda`) | The existing likelihood/Brent kernel, exact terminal-feasible lambda bound and status reporting; no optimizer replacement. |
| `src/d_stream.cpp` (`fast_d`) | Serial R RNG draws consumed in bounded native batches, without a whole-tree-by-simulation matrix. |
| `src/fast_signal.cpp` (`fast_ace`) | Preallocated workspaces reused across repeated traversals; repeated allocation and setup removed. |
| `src/delta_diagnostics.cpp` (`fast_delta`) | Compiled summaries and optional retained outputs; collecting summaries does not request extra random draws. |

The shared K refactor preserves its operation ordering and translation by a
baseline tip before both residual solves. All existing specialist tests were
rerun. The BM tail still uses strict numerical `H_sim >= H` with the same
plus-one formula; the separate permutation-tail tie tolerance was not imported.
No new parallel execution was added.

Matched qualification used two rounds with opposite arm order, three measured
runs per arm per round, and the same inputs and seeds. End-to-end medians compare
the previously accelerated implementation (`.rlib_heterogeneity_stage1`) with
the C++ pipeline:

| Scenario | Previous acceleration | C++ pipeline | Further speedup |
|---|---:|---:|---:|
| K, 50K tips, one trait | 1.350 s | 0.075 s | 18.00× |
| K, 50K tips, 32 traits | 1.575 s | 0.190 s | 8.29× |
| K + BM199, 50K tips | 2.715 s | 1.060 s | 2.56× |
| K + BM999, 50K tips | 7.875 s | 5.030 s | 1.57× |
| K + BM199, 2K tips, eight traits with missingness | 0.510 s | 0.330 s | 1.55× |
| Lambda, 2K tips, four traits | 0.110 s | 0.050 s | 2.20× |

A 5K-tip comb tree was also checked, but its C++ timings reach the timer's
resolution; no reliable speedup claim is made for that short scenario. Coverage
and actual selected work must be considered for imbalanced trees. These results
are specific to the recorded Windows/R environment and workloads.

All **42 paired runs** produced identical result tables, summaries, selections,
retained null values and final RNG states (excluding timing/call fields). All
estimate, BM-envelope and null-dispersion absolute differences were zero.

Final cold-start verification found and corrected a native RNGScope boundary:
loading an R dependency inside a non-random native call could overwrite the
dependency's restored seed when the outer scope exited. Deterministic native
entries now use `rng = false`; the two simulation entries retain RNGScope.
Four cold/warm pairs in eight fresh R processes confirm identical point/BM
results and RNG states under both Inversion and Box-Muller normal generation.
The final build was additionally rerun on all seven benchmark scenarios,
three times each: **21/21** results remain identical to the paired baseline.
Its medians are 0.07 s for K/50K, 1.14 s for BM199/50K and 5.11 s for BM999/50K.
See `cpp-acceleration/final-parity.log` and `cold-start-parity.log`.

The native tests also cover random/balanced/star/deep trees, shuffled edge rows
and internal IDs, exact mask sharing, large trait offsets, two uniform RNG kinds,
two normal RNG kinds, chunk invariance and a request exceeding the buffer budget.

In a bounded memory scenario (8,192 tips, two 4,096-tip crowns, BM999 and an
oversized `simulation_chunk`), independent R processes had observed peak working
sets of approximately **401 MiB before / 183 MiB after**. The PowerShell driver
reads the OS peak counter every 40 ms. This is a one-run-per-arm diagnostic,
including R and package/input memory, rather than a universal bound. The 64 MiB
policy limits dominant active block buffers; fixed tree arrays, inputs, outputs,
and R's allocation/GC behavior require additional memory. At least one replicate
is processed when that replicate alone exceeds the budget. Simulated K values,
null values, probabilities and final seeds matched in this scenario.

Evidence and reproduction scripts are in `cpp-acceleration/`: `benchmark.ps1`,
`benchmark.R`, `compare.R`, `speedups.csv`, `parity.csv`, `memory.ps1`,
`memory.csv`, `cold-start.ps1`, `test.R` and `source-check.log`. The latest source
check returned **Status: OK**; the full suite was run separately. The public
function signature and result layout are unchanged, and the workspace library
and source archive have been updated.

## Code cleanup and further performance audit (2026-10-13)

Removed four obsolete R helpers and two obsolete native entries from production,
retaining independent extraction/fit oracles in test helpers. BM reference
propagation is now entirely in R, independent of the production simulator.
Removed unused K arrays, duplicate matrix checks and unused BM diagnostics;
complete-data scans defer loading ape. The public interface is unchanged.

An additional precision/interpolation cache was measured and then withdrawn:
it offered only small BM changes and slightly slower multi-trait throughput.
Its separate records are in `cleanup/cache-experiment/`; they do not describe
the final source. The retained cleanup, verification, final benchmarks and
prioritized performance directions are documented in
[`cleanup/AUDIT.md`](cleanup/AUDIT.md). Final check/test logs are in `cleanup/`.

## Short API, tree annotation and overall signal (2026-10-13)

The public entry point is now `fast_hetero()`; the earlier development name
is absent from exports, the namespace and help aliases. Tests and active
documentation use the short name. Historical benchmark libraries retain their
original interfaces; a benchmark-only lookup supports both historical and
current installations without adding an alias to the package.

`plot(h, tree = tree)` colours crown-clade K/lambda estimates on the original
tree. The native mapping handles shuffled edges/internal IDs, descendant
identity checks, canonical plotting IDs, stem ownership and nested crowns
(including unavailable inner crowns). Default rightwards/fan plots reserve
annotation space; node labels refer to the original supplied tree.

Optional `overall`/`overall_test` compute signal on the complete supplied tree,
including tips outside selected crowns, through `fast_k()`/`fast_lambda()`.
Tests compare estimates, probabilities, lambda LR, trait-specific observed
counts and final seeds with these existing production APIs. Overall K
randomizations run after the BM dispersion calculation, preserving its results.
Plot annotations distinguish the complete-tree signal test from dispersion
testing, explicitly identify unperformed/unavailable tests, and consume no RNG.
Lambda supports local estimates, SD(lambda), overall lambda and its LRT against
zero; a lambda dispersion test is still unavailable.

The 128-tip simulated example has overall K = 0.6030019, randomization
P = 0.005025126 (199 draws); SD(log K) = 1.576193, BM P = 0.005
(199 simulations). Overall lambda = 0.8295435, LRT P = 1.184975e-13;
SD(lambda) = 0.6817395. These are simulated illustration results.
All 21 repeated benchmark cases have identical clade tables, dispersion
summaries, selections, retained null values and RNG states compared with the
pre-rename cleanup build (`tree-view/parity.csv`). Whole-tree estimation/testing
remains optional and was disabled for this comparison.
The final full suite passed 9,284 assertions with no failures/errors/warnings
and one existing skip; the fresh source check returned **Status: OK**.
On the recorded 50K-tip/449-crown input, 30 native annotation mappings took
0.10 s in isolation (approximately 3.3 ms per mapping). This excludes tree
compilation, species matching and graphics rendering, and is not end-to-end
plot time or a universal performance bound.
The checked installation and DLL were copied to `.rlib_heterogeneity` and
`src/fastphylosig.dll`. Separate fresh R processes successfully loaded the
workspace installation and the source via `pkgload::load_all(compile = FALSE)`,
ran K/lambda overall tests and rendered annotated trees. Both confirm that the
old function is absent. See `tree-view/workspace-installed.log` and
`workspace-source.log`.
See `tree-view/README.md`, the example PNG/PDF, CSVs and validation logs.
An initial check caught an incorrect namespace qualifier on the P-value
formatter; it was corrected to `base::format.pval()` before final validation.

## Package checks and archive

- Clean source compilation and installation succeeded (`install.log`).
- Source and installed-package checks with
  `R CMD check --no-install --no-tests --no-manual --library=<clean library>`
  returned **Status: OK** (`check-installed.log`). Examples ran; tests were
  executed separately in the full suite above.
- A separate source-archive check, including fresh compilation and installation,
  with `R CMD check --no-tests --no-manual` returned **Status: OK**
  (initial: `source-check.log`; final accelerated archive:
  `cpp-acceleration/source-check.log`; latest cleanup:
  `cleanup/source-check.log`; current API/tree-view build:
  `tree-view/source-check.log`). Tests were executed separately as described above.
- An earlier source-check reinstall failed in its temporary build directory
  with a missing make target (`check.log`); the retry in a fresh shorter
  temporary directory completed successfully. Both logs are retained.
- The final source archive is `fastphylosig_0.3.1.9000.tar.gz`. It contains
  the package sources and tests, without unrelated research environments or
  existing release/build directories.
- Earlier archive SHA-256 (correctness/performance audit):
  `b162172da11c20069a916a6edc82d14bcea68707957a4d467262e14ee0d58548`.

The illustrative plot (`demo.png`) uses **simulated data**: four 32-tip clades
with deep-split structure, BM, independent values, and alternating terminal
values. Their observed SD(log K) is 1.576 and their BM reference P is 0.005
(199 replicates, no simulated dispersion at least as large as observed).
Grey bars are BM reference envelopes, not confidence intervals for K.

## Shared-lambda and count-conditioned binary D extension (2026-10-13)

`fast_hetero(method = "lambda"/"D", test = TRUE)` now returns
`P_heterogeneity` and Monte Carlo uncertainty. Lambda fits a common covariance
parameter with separate crown means/variances, then uses parametric bootstrap
SD(lambda). D fits a shared symmetric ER rate, conditioning on each observed
state count in both likelihood and simulation; polynomial pruning and downward
sampling retain those counts without rejection. D normalization and count
probabilities are reused in C++. These are approximate model-reference tests,
not a universal equality-of-signal test or an original Fritz–Purvis test.

Independent validation covers dense Gaussian likelihood/bootstrap, R D
contrasts and normalization, exact ancestral-state enumeration and conditional
SD(D) moments, labels, missingness, RNG, failures and annotated tree plots.
The full suite passed **9,323 assertions**, 0 failures/errors/test warnings,
with one existing external bridge skip. A fresh source-archive check including
compilation returned **Status: OK** (`lambda-d/source-check.log`). Workspace
installation and source-DLL smoke tests passed in fresh R processes.

Bounded calibration on four 32-tip crowns used 100 datasets per scenario and
199 bootstrap draws. At alpha .05, null rejection rates were 4–7% for lambda
and 0–7% for D (0 failures). D was conservative at low common rates; the
selected contrast gave 87% rejection for lambda and 29% for D. These limited
simulations do not establish general nominal calibration or power. Data,
binomial intervals, scripts, illustrative plots and logs are in `lambda-d/`.
D count conditioning increases work with crown size; testing is restricted
to 1000 observed tips per crown, with default automatic crowns <=200 tips.

Isolated full scans on the recorded 50K-tip/440-crown random tree took 6.54 s
for lambda and 8.64 s for D, using **19 heterogeneity replicates** (D
normalization: 199). On the 5K-tip/42-crown tree they took 0.64 s and 0.62 s.
These measurements include matching, tree compilation, selection and native
calibration, and exclude package compilation, overall signal and graphics.
They are absolute times for these inputs and replicate settings, not timings
for the default 999 draws or a universal speedup claim. All cases returned
`status = "ok"`; see `lambda-d/benchmark.csv` and session metadata.

## Shared-kernel reuse audit (2026-10-13)

The D public and heterogeneity paths now share prepared contrast factors;
lambda fitting reuses workspaces, and binary conditional simulations reuse
count-probability buffers. Five before/after cases, including results, overall
tests, retained null draws and final RNG state, are identical. The fresh full
suite passed **9,329 assertions**, with no failures/errors/test warnings and
one existing skip. Clean source compilation and examples returned **Status:
OK**; installed and source-DLL smoke tests passed in fresh R processes.
See `kernel-reuse/README.md`, `parity.csv`, `tests.log`, `source-check.log`,
`workspace-smoke.log` and `source-smoke.log`.

Sequential timings show modest differences and overlap, including drift in
the unchanged K control. They do not establish a general speedup. The documented
benefits are shared code, fewer repeated factors and smaller D batch buffers.
This earlier audit covered K/lambda/D. Delta has since been integrated; see the
current validation above. Superseded duplicate integration code was removed
from `kernel-reuse/delta-draft`, retaining design and literature-review notes.

## Delta scan timing

Earlier alternating measurements on the integration source used 50,000 tips and
449 non-overlapping 20--200-tip crowns, one four-state trait, ER and specialist
default MCMC controls (10,000 iterations, thinning 10, burn 100). Three measured
repetitions per engine gave median **18.98 s** for `fast_hetero` versus **36.84 s**
for direct specialist calls per extracted subtree (**1.94x**). Native extraction
is included; specialist extraction is excluded. Estimates match within 1e-10
and final RNG states are identical; all 449 crowns are estimable. Whole-tree
estimates and all tests were excluded. The 5K case has no clear speedup
(2.35 s versus 2.41 s), and this does not establish the same gain for ARD,
overall Delta, reference simulations or randomizations. See `delta/benchmark.csv`,
`delta/benchmark.log`, session/source hashes and the reproducible script.
# Follow-up audit (0.3.3 development)

The five-goal audit is in `five-goals-20261013/`. Modern saved-trace diagnostics
are independently validated against pinned posterior reference functions.
Lambda reference fitting now defaults to pooled REML; explicit ML retains the
historical reference. Earlier calibration numbers below refer to their original
versions and controls. Historical `nsim_MCMC_warning` included ACE warnings;
the follow-up separates actual MCMC flags and `nsim_fit_warning`. Do not read
old aggregate warning percentages as pure MCMC failure rates.

