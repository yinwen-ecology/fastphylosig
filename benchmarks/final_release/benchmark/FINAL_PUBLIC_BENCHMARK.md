# fastphylosig 0.2.0 Final Public Benchmark

Status: **PASS for the single-machine qualification scope**.

This report records one fresh, serialized run from the frozen source commit
`084e979af137355b4f976faff6e0079237f23f72`. No historical timings or
cross-platform results are included.

## Provenance

- Package at timing run: `fastphylosig 0.2.0.9000` (development metadata)
- Commit: `084e979af137355b4f976faff6e0079237f23f72`
- Branch recorded by the harness: `codex/fastphylosig-0.2.0-dev`
- R: `R 4.6.1 (2026-06-24 ucrt)`
- Platform reported by R: Windows `10` compatibility identifier, `x86-64`;
  host OS was Windows 11 x64 (build 26200)
- C++17 compiler recorded by R: `g++` (installation output identified GCC 14.3.0)
- OpenMP capability: `capabilities("openmp")` was
  `disabled_or_unavailable` in the R harness; the final install transcript
  contains `-fopenmp` for C++ compilation and linking. Timing used serial
  public settings (`ncores=1`) and makes no cross-core claim.
- BLAS/LAPACK: the benchmark harness did not expose named BLAS/LAPACK
  entries; the final R session recorded default matrix products and LAPACK
  `3.12.1`
- Locale: `C`
- Logical cores: 16; formal calls were serialized and D used `ncores=1`
- Formal repeats: 10 per completed route/cell, with warmup and alternating
  route order

The formal grid used balanced, random, and pectinate fixed fixtures; sizes
500, 2000, and 5000; matrix workloads at n=500 with 8, 32, and 100 traits;
and K/D permutation-heavy workloads with `nsim=999`. A separate fixed n=100
external-reference guard was run where optional packages were available.

## Completeness And Guards

- Formal timing rows: 1260 non-reference rows, all `PASS`
- Completed timing cells: 126, each with exactly 10 observations
- Timing warnings: 0
- Guard rows: 9/9 `PASS`
- Memory rows: 5/5 `PASS`
- Formal completeness: `PASS`
- `phytools` reference guards: 6/6 K/lambda rows `PASS`
- `caper` reference guards: 3 rows marked `UNAVAILABLE`; the installed
  `phylo.d` call returned `Names column '1L' not found in data frame 'data'`.
  These optional rows were not used for timing or speed claims.

The guards covered raw/prepared parity, deterministic replay, matrix parity,
and input immutability before the corresponding timing cells were retained.

## Public Workloads

The following values are medians in seconds with IQR in parentheses. Each
entry aggregates the three tree shapes and 10 repeats at the stated size;
raw and prepared routes were paired and alternated.

### A/B: Single Trait

| n | K raw | K prepared | lambda raw | lambda prepared | fast_signal K raw | fast_signal K prepared |
|---:|---:|---:|---:|---:|---:|---:|
| 500 | 0.060 (0.010) | 0.010 (0.020) | 0.080 (0.010) | 0.030 (0.000) | 0.055 (0.020) | 0.010 (0.010) |
| 2000 | 0.230 (0.040) | 0.030 (0.010) | 0.315 (0.060) | 0.115 (0.020) | 0.245 (0.080) | 0.030 (0.010) |
| 5000 | 0.580 (0.140) | 0.070 (0.020) | 0.795 (0.080) | 0.285 (0.030) | 0.560 (0.060) | 0.065 (0.010) |

### C: Batch Matrix, n=500

| traits | K raw | K prepared | fast_signal K raw | fast_signal K prepared |
|---:|---:|---:|---:|---:|
| 8 | 0.050 (0.010) | 0.010 (0.020) | 0.050 (0.010) | 0.010 (0.020) |
| 32 | 0.060 (0.010) | 0.010 (0.010) | 0.055 (0.010) | 0.020 (0.010) |
| 100 | 0.060 (0.000) | 0.020 (0.010) | 0.060 (0.010) | 0.020 (0.020) |

### D: Permutation-Heavy, nsim=999

| n | K raw | K prepared | D raw | D prepared |
|---:|---:|---:|---:|---:|
| 500 | 0.090 (0.020) | 0.050 (0.010) | 0.220 (0.000) | 0.170 (0.010) |
| 2000 | 0.350 (0.070) | 0.160 (0.010) | 0.905 (0.060) | 0.720 (0.070) |
| 5000 | 0.955 (0.190) | 0.440 (0.040) | 2.555 (0.200) | 1.970 (0.160) |

These are public end-to-end timings. They include the route's documented
validation and result construction; they are not numerical-kernel timings.

## Memory Characterization

The prepared context was measured with R `object.size`, saved as an
uncompressed RDS, and read and validated in a fresh R process. These are
object and serialization measurements, not peak RSS measurements.

| n | object.size (bytes) | RDS (bytes) | save (s) | fresh read + first validation (s) |
|---:|---:|---:|---:|---:|
| 500 | 1375896 | 1017237 | 0.00 | 0.01 |
| 1000 | 2714896 | 2005522 | 0.03 | 0.04 |
| 5000 | 13434896 | 9955664 | 0.08 | 0.23 |
| 10000 | 26834896 | 19893329 | 0.16 | 0.37 |
| 20000 | 53654896 | 39878633 | 0.28 | 0.75 |

## Evidence Files

- `final_benchmark.csv` (1269 rows including reference rows), SHA-256:
  `9a547fa1990adfe60d416f17535cf2bd4f5470fa23b4e69cd0a441fb5f3c8f41`
- `memory_characterization.csv`, SHA-256:
  `1ae6810a3b2da3274fb4e719f012336343fa7e62b04f6815c3097b1ef668e36a`
- `benchmark_status.csv`, SHA-256:
  `8da422a10fb016ecd39c7a2159dea6d94bd68cae7cb7fb997f277e6dc5fc011f`
- `benchmark_manifest.txt`
- `MEMORY_CHARACTERIZATION.md`
- `FINAL_BENCHMARK_PROTOCOL.md`

The timing run predates the final one-line `DESCRIPTION` version transition
from `0.2.0.9000` to `0.2.0`; `git diff 084e979..c084ced` is empty for
production R/C++ source, tests, and estimator documentation, so these are
code-equivalent release timings. The raw CSV retains its original package
version field rather than being relabeled.

Cross-platform qualification was not run. Linux, macOS, R-devel, alternate
BLAS, and other toolchain combinations remain deferred validation items.
