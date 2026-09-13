# fastphylosig 0.2.0 Final Local Release Evidence Audit

Audit date: 2026-09-13

## Decision

```text
FINAL_TARBALL_IDENTITY = PASS
EXACT_TARBALL_SCIENTIFIC_VALIDATION = PASS (22/22)
FULL_TEST_SUITE = PASS (8336/0/0/0)
R_CMD_CHECK_NO_MANUAL_TIMINGS = PASS (Status: OK)
ALL_PUBLIC_EXPORTS_SMOKE = PASS (13/13)
CAPER_OBSERVED_D_PARITY = PASS
API_NAMESPACE_AUDIT = PASS
DEPENDENCY_AUDIT = PASS
DOCUMENT_CONSISTENCY = PASS
TARBALL_HYGIENE = PASS
FASTPHYLOSIG_0_2_0 = PUBLIC_RELEASE_READY_LOCAL_VALIDATION
```

## Identity

| Item | Value |
|---|---|
| Source commit | `3f84775bc96e00820aa39522245cb63175f79a5d` |
| Package version | `0.2.0` |
| Tarball SHA-256 | `d0178b870c6ab77158a434a9fd03c871814c49e360a3d288b4dda218385aeb67` |
| Source archive SHA-256 | `162603bd8aec228b692cd14181a0c942cd8736de2d167797449e13924cd1effd` |
| Local runtime | Windows 11, R 4.6.1, Rtools45/GCC 14.3 |

## Dynamic gates

The tarball built and installed successfully into a fresh private library.
The focused scientific suite reported 22 PASS, no failure or skip, and no
unexpected warning. The explicit ill-conditioned ACE non-convergence warning
remains expected and visible. Full testthat reported 8336 PASS with no
failure, warning, or skip. The authoritative no-manual check returned
`Status: OK`.

The all-exports smoke invoked all 13 namespace exports. The deprecated
`match_phylo_data` alias emitted exactly its expected warning; no other export
warning or error occurred. caper 1.0.4 and fast_d returned the same observed D
contrast sum, 2, with difference 0 and no warning. This caper result is scoped
to observed D only.

## Executed environment-limited checks

The `--as-cran` run completed package installation, dependencies, code,
examples, and tests, and no longer reports an unstated `mvtnorm` dependency.
It ended with one ERROR, one WARNING, and three NOTEs because `pdflatex` is not
installed and because of CRAN/environment checks. The ordinary manual check
ended with the same one ERROR/one WARNING from unavailable `pdflatex`.

These are retained as exact evidence and are not called successful checks.
They do not block this explicitly local, non-CRAN release decision. CRAN and
cross-platform readiness are not claimed.

## Static gates

- `DESCRIPTION` declares R >= 4.1.0 and the complete runtime/test dependency
  set, including `mvtnorm` in `Suggests`.
- `NAMESPACE` exports exactly the intended 13 functions.
- Current public documentation contains no stale 0.2.0.9000 development label.
- The tarball contains 106 entries and no benchmark, paper, docs, Rcheck,
  Rout, compiled-object, local-path, or temporary-library payload.
- No executable package file changed after the prior release candidate.

## Superseded evidence

The `c084ced`/`c3f05d...` and `146c03c`/`366b57...` candidates are historical
only. Their original tarballs and evidence are retained under the closure
archive; neither is authoritative for 0.2.0.
