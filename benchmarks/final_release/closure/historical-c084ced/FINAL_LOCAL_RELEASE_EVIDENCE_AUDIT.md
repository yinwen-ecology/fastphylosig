# fastphylosig 0.2.0 Final Local Release Evidence Audit

Audit date: 2026-09-13

## Decision

```text
FINAL_TARBALL_IDENTITY = PASS
EXACT_TARBALL_SCIENTIFIC_VALIDATION = ENVIRONMENT_BLOCKED
FULL_TEST_SUITE_ON_EXACT_INSTALLED_ARTIFACT = ENVIRONMENT_BLOCKED
R_CMD_CHECK_NO_MANUAL_TIMINGS = PASS (existing evidence)
R_CMD_CHECK_AS_CRAN = NOT_RUN (environment blocked)
ALL_PUBLIC_EXPORTS_SMOKE_ON_FINAL_AUDIT_LIB = ENVIRONMENT_BLOCKED
HIGH_RISK_REGRESSION_ON_EXACT_INSTALLED_ARTIFACT = ENVIRONMENT_BLOCKED
API_NAMESPACE_AUDIT = PASS (static)
DEPENDENCY_AUDIT = PASS (declarations); RcppArmadillo version not captured
DOCUMENT_CONSISTENCY = BLOCKER (stale 0.2.0.9000 release wording)
TARBALL_HYGIENE = PASS
FINAL_HASH_MANIFEST = PASS for existing artifacts
FASTPHYLOSIG_0_2_0 = NOT_READY
```

No production R/C++ code, estimator, RNG policy, thread policy, or public API
was changed. No cross-platform qualification or new performance work was run.

## 1. Authoritative identity

| Item | Result |
|---|---|
| Release commit | `c084ced0c9e44966369f6b0106d185d40909237b` |
| Branch | `codex/fastphylosig-0.2.0-dev` |
| Package version | `0.2.0` |
| Tarball | `fastphylosig_0.2.0.tar.gz` |
| Tarball SHA-256 | `c3f05d59055c72eec6a9765cdd9a68803ab325583d72fc7d221ec25f028bba63` |
| Source archive SHA-256 | `a20e5d2065565b4eae6c1f312f8c249df034b69bd72a95701b957f53d95980db` |

`git diff` shows no uncommitted change in `DESCRIPTION`, `NAMESPACE`, `R/`,
`src/`, `tests/`, `man/`, `README.md`, `USAGE_zh.md`, or `NEWS.md`. The current
working tree contains the expected repository-only review/evidence edits and
unrelated historical directories; these are not production source changes.

## 2. Exact-tarball dynamic gates

The previously recorded final smoke and scientific runs used the final commit
and final package version, but this audit could not create the newly required
`FINAL_AUDIT_LIB`. Invoking the installed R 4.6.1 executable was denied by the
default sandbox; the required escalation was rejected because the tool usage
limit had been reached. No workaround or indirect execution was attempted.

Consequently the following are retained as historical evidence, not falsely
reclassified as a fresh exact-tarball run:

- scientific suite: `22/22 PASS`, 0 fail/skip/unexpected warnings;
- full source-tree test: `8336 PASS`, 0 fail/warn/skip;
- existing final clean-install smoke: PASS for the major raw/prepared paths;
- existing `R CMD check --no-manual --timings`: Status OK, 0 ERROR/WARNING/NOTE.

The fresh-library scientific rerun, installed-artifact full test, all-export
smoke, `R CMD check --as-cran`, and ordinary manual check remain
`ENVIRONMENT_BLOCKED` in this audit. The `caper` reference rows in the existing
benchmark are `UNAVAILABLE` because the installed call rejected the fixture;
they were not used as a false PASS.

## 3. Static package and source audit

- `NAMESPACE` contains exactly 13 intended top-level exports and registered
  S3 methods only; forbidden public selectors (`fast_sig`, `fast_delta2`,
  `fast_ace2`, `backend`, `engine`, `unsafe`, `trust`, `skip_check`) are absent.
- `DESCRIPTION` declares runtime `ape`, `parallel`, `Rcpp`, and `stats`, links
  `Rcpp` and `RcppArmadillo`, and puts `caper`, `phytools`, and `testthat` in
  `Suggests`. C++17 and optional OpenMP are declared.
- Default worker settings are serial (`ncores = 1`) for K, D, Delta, and the
  dispatcher; documentation records that small parallel jobs may be slower.
- The memory protocol uses `saveRDS(..., compress = FALSE)`, so the existing
  “uncompressed RDS size” wording is accurate.
- The extracted tarball contains 106 entries and is 231,978 bytes. It contains
  no benchmark directory, manuscript, Rcheck/Rout payload, local path, or
  temporary audit library. The `stage2b2b` files are test-only oracle fixtures
  and are classified `DEVELOPMENT_ONLY`; the largest payload is a retained
  test fixture, not a debug or benchmark artifact.

The full source/debug/path classification is in
`FINAL_DOCUMENT_CONSISTENCY_AUDIT.md` and `audits/PACKAGE_HYGIENE_AUDIT.md`.

## 4. External reference evidence

Existing local evidence records R 4.6.1, `ape` 5.8.1, and `Rcpp` 1.1.2.
`phytools` K/lambda guards passed in the public benchmark. `caper::phylo.d()`
was attempted but its current call returned a fixture/data-column error, so
`CAPER_REFERENCE_PARITY = NOT_RUN`; no caper PASS is claimed. The installed
`RcppArmadillo` version was not captured in the existing session record.

## 5. Tag status

`v0.2.0` does not currently exist. Because the strict dynamic gates are not
complete and the tree is not clean of repository-only evidence edits, no tag
was created. Once the documentation is corrected, the tarball rebuilt, and all
gates pass, the exact command is:

```powershell
git tag -a v0.2.0 c084ced0c9e44966369f6b0106d185d40909237b -m "fastphylosig 0.2.0"
```

## 6. Required next action

The release is not marked locally ready by this audit. First resolve the
documentation consistency blocker, then rerun the dynamic gates from a fresh
temporary library after rebuilding the tarball. Do not start 0.2.1/0.3.0 or
any estimator/performance work as part of this audit.
