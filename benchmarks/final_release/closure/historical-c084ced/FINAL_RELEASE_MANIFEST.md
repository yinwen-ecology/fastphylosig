# fastphylosig 0.2.0 Final Release Manifest

Release decision date: 2026-09-13

## Superseding audit status

The gate table below preserves the preceding local release-hardening evidence.
The later `FINAL_LOCAL_RELEASE_EVIDENCE_AUDIT.md` supersedes its overall
decision: current status is
`FASTPHYLOSIG_0_2_0 = NOT_READY_PENDING_FINAL_LOCAL_AUDIT`.
The authoritative commit and tarball identity remain unchanged. No production
code was modified.

## Authoritative identity

- `PRIOR_LOCAL_RELEASE_EVIDENCE = RELEASE_READY (SUPERSEDED)` (local Windows/R 4.6.1 scope)
- Authoritative release source commit:
  `c084ced0c9e44966369f6b0106d185d40909237b`
- Branch: `codex/fastphylosig-0.2.0-dev`
- Package version: `0.2.0`
- Authoritative tarball: `fastphylosig_0.2.0.tar.gz`
- Tarball SHA-256: `c3f05d59055c72eec6a9765cdd9a68803ab325583d72fc7d221ec25f028bba63`
- Source snapshot SHA-256:
  `a20e5d2065565b4eae6c1f312f8c249df034b69bd72a95701b957f53d95980db`

## Verification hashes

- Final testthat log `testthat.Rout`: `f4bfd17577aef8a98e6f05271a53c4cbbf708498f676fc305d4d7c024f4ef73f`
- Final check log `00check.log`: `ec9ccdbfcdf106de3f40d074861f879869669d1ab8f17bca97313e2e28650318`
- Scientific validation log: `607ca7f0576eca98667791ac4b8a20ad76d84801a54e3f2de6c7ebc96e657289`
- Scientific status CSV: `193d0d3501c32701e7af1bf660b536a29a42b679a2e8db88ab59f7a97b794b2d`
- Scientific summary CSV: `1cb7cc5c95323f542558580b2b4564e040ded40cea9cc234cef2b3835b65a2b1`
- Final benchmark timing CSV: `9a547fa1990adfe60d416f17535cf2bd4f5470fa23b4e69cd0a441fb5f3c8f41`
- Final benchmark memory CSV: `1ae6810a3b2da3274fb4e719f012336343fa7e62b04f6815c3097b1ef668e36a`
- Final benchmark status CSV: `8da422a10fb016ecd39c7a2159dea6d94bd68cae7cb7fb997f277e6dc5fc011f`
- Final benchmark summary `FINAL_PUBLIC_BENCHMARK.md`:
  `7914d3171d16feedcef640ccfaf0379cd17befbcf4feabc28100dc609f82ef7e`
- Memory summary `MEMORY_CHARACTERIZATION.md`:
  `94a95d8a1efb125e4ea7f073cc3c3e75ccdb8bb96ea867c697d9c8cb1bd30379`

## Gate summary

| Gate | Result | Evidence |
|---|---|---|
| Full testthat | PASS: 8336 pass, 0 fail, 0 warning, 0 skip | `testthat.Rout` |
| `R CMD check --no-manual --timings` | PASS: Status OK, 0 ERROR/WARNING/NOTE | `00check.log` |
| Clean install smoke | PASS: raw and prepared public entry points | `clean_install_smoke.log` |
| Scientific parity and contract suite | PASS: 22/22, 0 unexpected warnings | `scientific/summary.csv` |
| API and contract documentation audit | PASS for API/V2 contracts; release-version wording is BLOCKER in current audit | `audits/audit_status.csv`, `FINAL_DOCUMENT_CONSISTENCY_AUDIT.md` |
| Final public benchmark correctness | PASS: 9/9 guards, 126/126 cells, 10 repeats | `benchmark/benchmark_status.csv` |
| Memory characterization | PASS: n=500, 1000, 5000, 10000, 20000 | `benchmark/memory_characterization.csv` |
| Cross-platform qualification | NOT RUN / DEFERRED by request | no claim made |
| Current final-local evidence audit | BLOCKED: fresh R gates unavailable; stale release wording requires correction | `FINAL_LOCAL_RELEASE_EVIDENCE_AUDIT.md` |

The public benchmark was run at the pre-version freeze commit
`084e979af137355b4f976faff6e0079237f23f72` with package version `0.2.0.9000`.
The final release commit differs from that commit only in `DESCRIPTION`'s
version field; the production R/C++ source, tests, and estimator behavior are
identical. This is recorded explicitly rather than relabeling the raw CSV.

## Historical artifacts

The repository-root `fastphylosig_0.1.0.tar.gz` and its earlier check records
remain historical evidence and are superseded. They are not used as the
0.2.0 artifact. No historical artifact is deleted.

## Deferred qualification

Linux/macOS, R-devel, alternate BLAS/LAPACK, no-OpenMP fallback, and
cross-core deterministic equivalence remain deferred. Non-blocking technical
debt is recorded in `audits/POST_0.2.0_TECH_DEBT.md`.
