# fastphylosig 0.2.0 Final Scientific Validation

## Decision

The focused scientific suite completed 22 case assertions:

- PASS: 22
- FAIL: 0
- SKIP: 0
- Unexpected warning cases: 0
- Fail-closed scientific validation: **PASS**
- Release blocker: **NO**

The suite includes a stable ACE parity fixture for the release gate. A
separate deliberately ill-conditioned fixture is retained as an explicit
expected-warning case; its warning is checked as a non-convergence state and
is not suppressed.

## Provenance

- Source commit supplied to the audit: `c084ced0c9e44966369f6b0106d185d40909237b`
- Package version: `0.2.0`
- R: 4.6.1 (2026-06-24 ucrt)
- Platform: Windows 11 x64, `x86_64-w64-mingw32`
- `ape`: 5.8.1
- Matrix/LAPACK: default / LAPACK 3.12.1
- Loaded package: final release smoke library from the preceding evidence run
  (the fresh `FINAL_AUDIT_LIB` was not created in the current audit; see the
  environment-blocked disposition below)
- Locale used by the audit: `LC_COLLATE=C`, `LC_CTYPE=C`
- Cross-platform qualification: **NOT RUN**, as requested

The machine emitted R startup messages that `C.UTF-8` could not be selected;
the audit then ran under the recorded `C` locale. These are environment
startup messages, not package warnings captured in `status.csv`.

## Coverage

The reproducible script covers:

- K estimate, inclusive upper-tail P, MCSE, controlled permutations, thread
  invariance, large-offset values, NA accounting, and raw/prepared parity.
- lambda estimate, log likelihood, LR/P fields, NA matching, and
  raw/prepared parity.
- D estimate, fixed random and Brownian nulls, MCSE and strict-tail
  accounting, input immutability, raw/prepared parity, and a rooted
  polytomy.
- Delta estimate, saved-chain diagnostics, P/MCSE/status, fixed permutation
  paths, NA pruning, raw/prepared parity, and same-seed replay.
- ACE log likelihood, rates, SE, ancestral likelihoods, stable raw/prepared
  parity, and the explicit `CI=FALSE` likelihood-only path.
- V2 internal-node renumbering, edge-row shuffling, delimiter/Unicode labels,
  heterogeneous branch lengths, branch association, fingerprints, and
  K/lambda/D/Delta/ACE propagation.
- Structural-field mutation rejection, root sensitivity, input immutability,
  V1/unknown schema rejection, and saveRDS/readRDS fresh-session gates.

## ACE Warning Disposition

The stable ACE fixture is the same seed-31 ER parity scenario used by
`tests/testthat/test-fast-signal.R`. Both `fast_ace(CI=TRUE)` and
`fast_ace(CI=FALSE)` returned `convergence = 0` with no warning. The
`ape::ace` comparison also emitted no warning. The fitted log likelihood was
`-42.4173076049793`; the rate difference from `ape::ace` was
`3.34e-08`.

The audit also retains the deterministic cyclic-state/high-heterogeneity
fixture from the earlier audit. Its branch lengths span `0.01` to `10`, the
ER rate is approximately `218.63`, and the likelihood is nearly flat at the
reported solution. `fast_ace(CI=FALSE)` returned finite results with
`convergence = 1` and message `singular convergence (7)`. This warning is
recorded in `status.csv` with `warning_expected=TRUE` and is validated as an
explicit non-convergence state. Independent diagnostics show that different
starting values reach rates from about 213 to 256 while changing the deviance
only at roughly the `1e-8` scale; an independent R eigen implementation agrees
with the C++ objective to about `1e-10` or better. This is an ill-conditioned
fixture/optimizer termination case, not evidence of a changed ACE estimator.

The current `0.2.0` `CI=FALSE` path skips the Hessian/SE calculation and
the rerun captured no `NaNs produced` warning. The explicit non-convergence
warning remains visible; no warning suppression, tolerance relaxation, or
production algorithm change was made.

The warning provenance, including convergence code and fitted values, is in
`ace_warning_provenance.csv`. The full diagnosis is in
`ACE_RELEASE_BLOCKER_DIAGNOSIS.md`.

## Evidence Files

- `run_scientific_validation.R`: reproducible focused suite, including the
  fresh-process `readRDS()` child mode.
- `status.csv`: case-level status, expected-warning classification, and
  captured warnings.
- `summary.csv`: machine-readable fail-closed summary (`PASS`).
- `ace_warning_provenance.csv`: stable ACE, pathological ACE, and `ape::ace`
  warning/convergence provenance.
- `ACE_RELEASE_BLOCKER_DIAGNOSIS.md`: ACE root-cause and release disposition.
- `scientific_validation.log`: complete captured run log.
- `environment.txt`: R/package/platform/session provenance.
- `scientific_context_v2.rds`, `fresh_session_result.rds`, and
  `fresh_session_result.rds.log`: persistence and fresh-session evidence.

## Reproduction

Run from the current final-audit staging directory. Install the authoritative
tarball into a new private library before invoking the validation script; do
not reuse a development library:

```powershell
Set-Location 'C:\Users\wenyi\.codex\visualizations\2026\05\26\019e6293-544a-7c11-9cb6-8324d897ca01\final-local-audit-20260913'
$env:FINAL_AUDIT_LIB = 'C:\Users\wenyi\.codex\visualizations\2026\05\26\019e6293-544a-7c11-9cb6-8324d897ca01\final-local-audit-20260913\lib'
$env:FASTPHYLOSIG_LIBRARY = $env:FINAL_AUDIT_LIB
& 'C:\Users\wenyi\AppData\Local\R\R-4.6.1\bin\R.exe' CMD INSTALL `
  --library=$env:FINAL_AUDIT_LIB `
  'C:\Users\wenyi\.codex\visualizations\2026\05\26\019e6293-544a-7c11-9cb6-8324d897ca01\final-local-audit-20260913\fastphylosig_0.2.0.tar.gz'
& 'C:\Users\wenyi\AppData\Local\R\R-4.6.1\bin\Rscript.exe' --vanilla `
  'run_scientific_validation.R' `
  --repo . --out 'scientific' `
  --source-commit c084ced0c9e44966369f6b0106d185d40909237b
```

Expected outcome for the current package state is process exit status 0, with
22 PASS assertions, no unexpected warnings, and
`scientific_validation=PASS`.

This command template was not executable in the current audit session because
the sandbox denied launching the local R executable after the usage-limit
escalation rejection; the fresh-install run is therefore recorded as
`ENVIRONMENT_BLOCKED` in the final-local audit.

## Scope Boundary

This release-validation update changed only evidence and audit files under
`benchmarks/final_release/scientific/`. It did not modify production R/C++
code, existing package tests, DESCRIPTION, documentation, or the public API.
No cross-platform qualification or performance benchmark was performed.
