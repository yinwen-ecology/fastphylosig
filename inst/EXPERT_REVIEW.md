# fastphylosig 0.1.0 expert review

**Review date:** 2026-08-24  
**Scope:** final local release-freeze review. This evidence file is excluded
from built package archives.

## Freeze decision

`FASTPHYLOSIG_0.1.0 = RELEASE_FROZEN`

`GITHUB_RELEASE_READY = TRUE`

The production estimators, defaults, numerical tolerances, traversal,
thread/RNG policies, and public function set are frozen. No new numerical
candidate, algorithm optimization, parser/traversal refactor, benchmark
expansion, or performance claim is authorized by this review.

## Authoritative release evidence

Exactly one artifact is authoritative:

| Field | Authoritative result |
|---|---|
| Artifact | `AUTHORITATIVE_RELEASE_ARTIFACT` |
| Package | `fastphylosig 0.1.0` |
| Tarball SHA-256 | `0A992F524C3CC5F5199423EFFC19E9399D97D761D857913BC2E38DA120AC4BD5` |
| Source snapshot SHA-256 | `AAE43848AFFE734BDA8644E1758F47526CD07109F9B9D5922F785F8E3B337E52` |
| Check-log SHA-256 | `0FDD674951C42907D6289AAFFA9D482C1F6BBB338F84C63F5BF23ABCB01FE295` |
| Command | `R CMD check --no-manual --timings fastphylosig_0.1.0.tar.gz` |
| Check result | `Status: OK`; 0 ERROR, 0 WARNING, 0 NOTE |
| testthat | 6,791 pass; 0 fail; 0 warn; 0 skip |
| Clean-install smoke | PASS |

The source snapshot digest is calculated over the 88 files extracted from the
authoritative tarball. Each record is the forward-slash relative path, a tab,
the lowercase per-file SHA-256, and a newline. Records are sorted by relative
path, encoded as UTF-8 without BOM, concatenated, and hashed with SHA-256.

The successful Windows check used `LC_ALL=C`, `LC_CTYPE=C`, and `LANG=C`
because the managed parent process supplied the unsupported locale
`C.UTF-8`. It used `R_PARALLEL_PORT=11673` after an independent two-worker
PSOCK smoke established that a non-default port was available. These are
check-environment settings; package code and runtime defaults were unchanged.

## Maintenance refactor audit

`MAINTENANCE_REFACTOR = SEMANTICALLY_EQUIVALENT`

- Trait-table conversion and naming checks were consolidated. Valid inputs
  retain their values, order, NA behavior, and downstream statistical path.
  Ambiguous or invalid species names now fail earlier. The historical named-
  vector column name `x` returned by `match_tree_data()` is retained.
- The shared type-7 quantile helper retains exactly
  `h = 1 + (n - 1) * p`, `floor(h)`, the same linear interpolation, endpoints,
  sorted indices, and tie behavior. It was moved, not redefined.
- ACE timing-hook removal deleted development instrumentation only. It did not
  change likelihood, pruning, RNG, tolerance, ordering, or thread behavior.
- The deleted tree-normalization wrapper contained only a direct call to
  `.safe_canonicalize_core()`; both call sites now invoke the same function
  directly.

No numerical semantic change was found. `PRODUCTION_FREEZE_BROKEN` is false.

## Release-contract audit

- `fast_d()` and `fast_delta()` reject simultaneous `x` and `X`.
- `fast_ace()` rejects invalid `method` values with a clear message.
- `fast_delta(test = FALSE)` does not validate or warn about unused `nsim`.
- Prepared-tree fingerprints include computationally relevant `Nnode`.
- Node labels and other non-computational attributes are not claimed as part
  of the fingerprint contract.
- `fast_d()` help now states that safe root-numbering and edge-order
  representation differences are normalized automatically.
- README, `USAGE_zh.md`, package help, and the implementation agree on the
  general and advanced workflows.

`RELEASE_CONTRACT_FIXES = PASS`

## Clean-install smoke

The authoritative tarball was installed into an empty temporary R library.
A separate fresh R session loaded `fastphylosig` and successfully called:

`fast_k()`, `fast_lambda()`, `fast_d()`, `fast_delta()`, `fast_ace()`,
`check_tree()`, and `fast_signal()`.

`CLEAN_INSTALL_SMOKE = PASS`

## Superseded evidence

The following records are retained only for provenance. None is a final or
authoritative release artifact:

| Status | testthat | Tarball SHA-256 | Check-log SHA-256 |
|---|---:|---|---|
| `SUPERSEDED` post-cleanup artifact | 6,790 | `2C1D569781436394E604BB9CF7291364EEAFDDD81A7B93D9A3BA5A28C35C1DA5` | `7F683F78333D2A5CF7CC63A47118C8F24959698316B2B4F51958A557473DE6A0` |
| `SUPERSEDED` minimal release patch | 6,787 | `D92D790851657B36D9EA5A9B5AA0DCD817A702F65CAA960A0D12EE0D0FDCA2B5` | `DB1AD001BDB047755845CBF1778D7758C16151BF98A0E57FA2E26C68DE19B335` |

Earlier pre-patch and compatibility artifacts remain explicitly historical in
`RC_READINESS.md`.

## Git and deferred validation

The release source is committed with message `fastphylosig 0.1.0 release`
and identified by tag `v0.1.0`. The commit/tag themselves are the Git
provenance; this committed file intentionally does not self-reference a commit
hash.

CRAN, Linux, macOS, R-devel, alternate BLAS, and no-OpenMP validation remain
deferred. No PASS is claimed for them.
