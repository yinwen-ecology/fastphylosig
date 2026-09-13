# fastphylosig 0.2.0 Final Release Provenance

Date: 2026-09-13

The source/tarball provenance below remains valid, but the current overall
release decision is superseded by
`FINAL_LOCAL_RELEASE_EVIDENCE_AUDIT.md`, which records a documentation blocker
and environment-blocked fresh R gates. This does not alter the source commit or
the existing artifact hash.

## Release source

- Release source commit: `c084ced0c9e44966369f6b0106d185d40909237b`
- Branch: `codex/fastphylosig-0.2.0-dev`
- Package version: `0.2.0`
- Release source snapshot: `fastphylosig-source-c084ced.tar`
- Release source snapshot SHA-256: `a20e5d2065565b4eae6c1f312f8c249df034b69bd72a95701b957f53d95980db`

The source snapshot is a deterministic `git archive` of the release source
commit. The final package tarball was built from this snapshot. Evidence files
under `benchmarks/final_release/` are repository-only and are excluded from
the package archive.

## Local build environment

- R: `R version 4.6.1 (2026-06-24 ucrt)`
- R platform: `x86_64-w64-mingw32/x64`
- Host OS: Windows 11 Home China, build `26200`, 64-bit
- CPU: AMD Ryzen 7 5800H with Radeon Graphics, 8 physical cores / 16 logical
  processors, reported maximum 3201 MHz
- Installed memory: 14,888,460,288 bytes (about 13.9 GiB)
- Toolchain: Rtools45; GCC/G++/GNU Fortran `14.3.0`; GNU Make `4.4.1`
- C++ standard: C++17 (`-std=gnu++17`)
- OpenMP build: enabled in the final install transcript (`-fopenmp` in C++
  compile and link commands). The benchmark harness reported
  `capabilities("openmp")` as unavailable; no cross-platform or no-OpenMP
  qualification is claimed here.
- BLAS/LAPACK: R default matrix products; LAPACK `3.12.1`; no separately
  named BLAS implementation was reported by the captured `extSoftVersion()`.
- Locale: `LC_COLLATE=C`, `LC_CTYPE=C`, system code page `65001`, timezone
  `Asia/Shanghai`. The host startup environment attempted `C.UTF-8`, which
  Windows R could not select; formal runs explicitly used the supported `C`
  locale.
- Runtime dependencies in the final smoke session: `ape 5.8.1`, `Rcpp 1.1.2`,
  `RcppArmadillo` as linked by the package, and base `parallel`/`stats`.

## Production-change audit

The V2 acceptance baseline is commit `0636452`. A source-level diff from that
baseline through the release source contains only:

- the release metadata version change in `DESCRIPTION`;
- a roxygen clarification in `R/prepare_tree.R` and its generated Rd wording,
  changing the ACE example to `fast_ace(x, prepared = ctx)`.

The executable bodies in `R/`, all `src/` C++ files, tests, estimator
definitions, numerical formulas, tolerances, RNG paths, thread policy, and
public exports are unchanged after V2 acceptance. The development freeze point
before the version transition was `084e979af137355b4f976faff6e0079237f23f72`;
`git diff 084e979..c084ced -- R src tests man README.md USAGE_zh.md NEWS.md`
is empty, proving that the final transition changed only `DESCRIPTION`.

Optimization states remain closed for 0.2.0:

```text
COMMON_PREPARATION_ENGINE_0_2_0 = FROZEN
CANONICALIZATION_CONTRACT_V2 = ACCEPTED_FROZEN
K_OPTIMIZATION_0_2_0 = CLOSED_FINAL
D_OPTIMIZATION_0_2_0 = CLOSED_FINAL
LAMBDA_OPTIMIZATION_0_2_0 = CLOSED
DELTA_OPTIMIZATION_0_2_0 = CLOSED
STATISTIC_SPECIFIC_OPTIMIZATION_0_2_0 = CLOSED
```

## Cross-platform boundary

Cross-platform qualification was intentionally not run in this release task.
Linux, macOS, R-devel, alternate BLAS, and no-OpenMP fallback remain
`NOT_RUN/DEFERRED`; they are not represented as PASS and are not release
blockers for this explicitly local release scope.
