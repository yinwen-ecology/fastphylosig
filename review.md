# fastphylosig 0.2.0 Final Release Expert Review

Review date: 2026-09-13

## Decision

```text
FASTPHYLOSIG_0_2_0 = PUBLIC_RELEASE_READY_LOCAL_VALIDATION
AUTHORITATIVE_RELEASE_SOURCE = 3f84775bc96e00820aa39522245cb63175f79a5d
AUTHORITATIVE_RELEASE_ARTIFACT_SHA256 = d0178b870c6ab77158a434a9fd03c871814c49e360a3d288b4dda218385aeb67
RELEASE_TAG = v0.2.0 -> 3f84775bc96e00820aa39522245cb63175f79a5d
CROSS_PLATFORM_QUALIFICATION = NOT_RUN / DEFERRED
CRAN_READINESS = NOT_CLAIMED
```

This decision is limited to the completed Windows 11 / R 4.6.1 local release
qualification. It does not claim Linux, macOS, R-devel, alternate BLAS/LAPACK,
or no-OpenMP validation.

## Frozen source

The authoritative package source is commit
`3f84775bc96e00820aa39522245cb63175f79a5d` on
`codex/fastphylosig-0.2.0-dev`. Relative to the preceding `c084ced` candidate,
only `DESCRIPTION`, `README.md`, `NEWS.md`, and `docs/index.md` changed.
`R/`, `src/`, tests, man pages, `NAMESPACE`, statistical definitions,
tolerances, RNG, thread policy, and public result structures are unchanged.

The final metadata repair declares `mvtnorm` in `Suggests`, matching the
existing test-only `mvtnorm::rmvnorm()` call. Public release wording now says
0.2.0 rather than 0.2.0.9000.

## Validation gates

| Gate | Result |
|---|---|
| Build and fresh install | PASS under R 4.6.1, Rtools45/GCC 14.3, C++17, OpenMP |
| Focused scientific suite | 22 PASS, 0 fail, 0 skip, 0 unexpected warning |
| Full source testthat | 8336 PASS, 0 fail, 0 warning, 0 skip |
| `R CMD check --no-manual --timings` | PASS, `Status: OK` |
| All 13 public exports smoke | PASS from the fresh final library |
| Deprecated alias behavior | PASS; one expected `match_phylo_data` warning |
| `caper::phylo.d()` observed-D parity | PASS; 2 vs 2, difference 0, no warning |
| Tarball hygiene | PASS; 106 entries, no development evidence leakage |

`R CMD check --as-cran` completed package installation, code checks,
dependency checks, examples, and tests successfully. Its final nonzero status
is entirely attributable to unavailable `pdflatex` (one PDF-manual ERROR and
one WARNING) plus three CRAN/environment notes. The ordinary manual check
reproduced only the same unavailable-`pdflatex` ERROR/WARNING. These logs are
retained exactly and are not mislabeled as `Status: OK`. They are not a blocker
for this explicitly local, non-CRAN release scope.

## Scientific disposition

K, lambda, D, Delta, and ACE retained their frozen estimators and contracts.
The focused suite covers controlled nulls, P/MCSE accounting, raw/prepared
parity, NA masks, large offsets, root sensitivity, V2 canonicalization,
mutation rejection, persistence, and same-seed replay. The deliberately
ill-conditioned ACE fixture continues to emit its explicit expected
non-convergence warning; no warning was suppressed or reclassified.

The optional caper check is intentionally narrow: it validates only the
observed D contrast sum. It does not claim parity for P values, simulated null
distributions, or performance.

## Authoritative hashes

| Artifact | SHA-256 |
|---|---|
| `fastphylosig_0.2.0.tar.gz` | `d0178b870c6ab77158a434a9fd03c871814c49e360a3d288b4dda218385aeb67` |
| exact source archive | `162603bd8aec228b692cd14181a0c942cd8736de2d167797449e13924cd1effd` |
| no-manual `00check.log` | `842e2ab998f918b30cf5c937baf4e96b3f3b833a5150c3fedeee4690a293ee0a` |
| full `testthat.Rout` | `402e4d2106ad4b540e0b74be6364e527dddcd4fe84ece71cab7feff1828ffa44` |
| scientific validation log | `4ca0fb5dec1acdcb994f6ce128ddc134bba0daf5069834d5656e2b1654bfbbf6` |
| all-exports smoke log | `631400f6a52f96c25c53db61cef018500fa3996f82d3b96c30f2f9f594fa62d0` |
| caper observed-D log | `8d2d34c69a7ad8a14e74599b77017abe8acdba3eae26803778d6b401090f0313` |

The older `c084ced` / `c3f05d...` package and the intermediate `146c03c` /
`366b57...` package are `SUPERSEDED/HISTORICAL`. They remain preserved under
the final-release closure evidence and are not presented as final artifacts.

The annotated `v0.2.0` tag was verified to reference the authoritative source
commit above. Repository-only evidence is recorded in a later commit and does
not alter the tagged package source.

## Expert conclusion

Within the tested local scope, the release evidence is coherent and supports
publishing 0.2.0. The only incomplete qualifications are explicitly deferred
platform/toolchain checks and local PDF-manual generation without a TeX
installation. No 0.2.1/0.3.0 work or new optimization is authorized by this
review.
