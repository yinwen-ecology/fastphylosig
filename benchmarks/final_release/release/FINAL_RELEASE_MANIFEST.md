# fastphylosig 0.2.0 Final Release Manifest

Decision date: 2026-09-13

## Authoritative identity

- Release source commit: `3f84775bc96e00820aa39522245cb63175f79a5d`
- Branch: `codex/fastphylosig-0.2.0-dev`
- Annotated tag: `v0.2.0` -> `3f84775bc96e00820aa39522245cb63175f79a5d`
- Package version: `0.2.0`
- Artifact: `fastphylosig_0.2.0.tar.gz`
- Artifact SHA-256:
  `d0178b870c6ab77158a434a9fd03c871814c49e360a3d288b4dda218385aeb67`
- Exact source archive SHA-256:
  `162603bd8aec228b692cd14181a0c942cd8736de2d167797449e13924cd1effd`
- Artifact size: 231,949 bytes
- Archive entries: 106

Exactly one artifact is authoritative. The following are retained only as
historical evidence:

- `c084ced0c9e44966369f6b0106d185d40909237b` tarball SHA
  `c3f05d59055c72eec6a9765cdd9a68803ab325583d72fc7d221ec25f028bba63`:
  `SUPERSEDED_PRE_FINAL`.
- `146c03c29dd15b8e367d0a16a4f39771a896f668` tarball SHA
  `366b57ac54d518afcc8f288a2f1cf2855c52aad4ce84e1ae675889363f07562b`:
  `SUPERSEDED_DEPENDENCY_DECLARATION_INCOMPLETE`.

## Gate summary

| Gate | Result | Evidence |
|---|---|---|
| Exact build and fresh install | PASS | `exact_tarball/build.log`, `install.log` |
| Scientific validation | 22/22 PASS; no fail/skip/unexpected warning | `scientific/scientific_exact_tarball/` |
| Full testthat | 8336 PASS; 0 fail/warn/skip | `testthat.Rout`, `testthat_status.csv` |
| No-manual package check | PASS, `Status: OK` | `00check.log` |
| All exported functions | 13/13 PASS | `exact_tarball/all_exports_smoke.log` |
| caper observed D | PASS, exact equality on scoped fixture | `exact_tarball/caper_parity.log` |
| API/document contract | PASS | `audits/audit_status.csv` |
| Tarball hygiene | PASS | `exact_tarball/tar_listing.txt` |
| Cross-platform | NOT RUN / DEFERRED | no claim made |

The `--as-cran` and ordinary manual logs are retained as executed evidence,
not as successful checks. Both reached the manual PDF phase after package
code/dependency/tests passed, then failed because `pdflatex` is unavailable.
The as-cran run also records three CRAN/environment notes. This release is not
being qualified for CRAN in this task.

## Verification hashes

- no-manual `00check.log`:
  `842e2ab998f918b30cf5c937baf4e96b3f3b833a5150c3fedeee4690a293ee0a`
- `testthat.Rout`:
  `402e4d2106ad4b540e0b74be6364e527dddcd4fe84ece71cab7feff1828ffa44`
- `testthat_status.csv`:
  `7410cd924240b80a10a91c70f2e654e35289dcb88b826013d9e312a62b255114`
- scientific log:
  `4ca0fb5dec1acdcb994f6ce128ddc134bba0daf5069834d5656e2b1654bfbbf6`
- all-exports smoke:
  `631400f6a52f96c25c53db61cef018500fa3996f82d3b96c30f2f9f594fa62d0`
- caper observed-D parity:
  `8d2d34c69a7ad8a14e74599b77017abe8acdba3eae26803778d6b401090f0313`

## Decision

```text
AUTHORITATIVE_RELEASE_ARTIFACT = fastphylosig_0.2.0.tar.gz
FASTPHYLOSIG_0_2_0 = PUBLIC_RELEASE_READY_LOCAL_VALIDATION
GITHUB_RELEASE_READY = TRUE
GIT_PROVENANCE = PASS (v0.2.0 verified)
CROSS_PLATFORM_QUALIFICATION = NOT_RUN / DEFERRED
CRAN_READINESS = NOT_CLAIMED
```
