# fastphylosig 0.1.0 release readiness

**Status date:** 2026-08-24

## Release state

| Gate | State |
|---|---|
| `FASTPHYLOSIG_0.1.0` | `RELEASE_FROZEN` |
| `GITHUB_RELEASE_READY` | `TRUE` |
| `CURRENT_R_VALIDATION` | `PASS` |
| `MAINTENANCE_REFACTOR` | `SEMANTICALLY_EQUIVALENT` |
| `RELEASE_CONTRACT_FIXES` | `PASS` |
| `CLEAN_INSTALL_SMOKE` | `PASS` |
| `CRAN_VALIDATION` | `DEFERRED` |
| `LINUX_VALIDATION` | `DEFERRED` |
| `MACOS_VALIDATION` | `DEFERRED` |
| `R_DEVEL_VALIDATION` | `NOT_RUN_ENVIRONMENT` |
| `SECOND_BLAS_VALIDATION` | `NOT_RUN_ENVIRONMENT` |
| `NO_OPENMP_VALIDATION` | `DEFERRED` |

This is a local GitHub-release freeze, not a CRAN or cross-platform claim.

## Authoritative artifact

Only the following evidence is authoritative:

| Field | Result |
|---|---|
| Artifact status | `AUTHORITATIVE_RELEASE_ARTIFACT` |
| Environment | Windows x86-64, R 4.6.1, Rtools 4.5, GCC/G++ 14.3.0 |
| Build | `R CMD build` from fresh ASCII staging |
| Check | `R CMD check --no-manual --timings`; `Status: OK` |
| Check totals | 0 ERROR, 0 WARNING, 0 NOTE |
| testthat | 6,791 pass; 0 fail; 0 warn; 0 skip |
| Tarball SHA-256 | `0A992F524C3CC5F5199423EFFC19E9399D97D761D857913BC2E38DA120AC4BD5` |
| Source snapshot SHA-256 | `AAE43848AFFE734BDA8644E1758F47526CD07109F9B9D5922F785F8E3B337E52` |
| Check-log SHA-256 | `0FDD674951C42907D6289AAFFA9D482C1F6BBB338F84C63F5BF23ABCB01FE295` |
| Clean-install smoke | PASS for package load and all seven required entry points |

The source snapshot algorithm is defined in `docs/validation_manifest.md`.
The successful check used the standard `C` locale and fixed an available
managed-environment PSOCK port. Neither setting changes package behavior.

## Scientific and API freeze

No statistic definition, numerical formula, type-7 interpolation, NA
behavior, ordering, boundary classification, tolerance, RNG stream policy,
thread policy, estimator, or public function was changed by the maintenance
refactor.

The trait helper consolidation rejects ambiguous species identifiers earlier
while preserving valid-input calculations. The type-7 helper is a literal
relocation of the existing formula. ACE timing hooks were development-only.
The tree-normalization wrapper was a one-line internal alias.

The final release-contract verification covers simultaneous D/Delta `x` and
`X`, ACE `method` validation, unused Delta `nsim` when testing is disabled,
the prepared-tree `Nnode` fingerprint, the explicitly limited node-label
contract, and user documentation. The historical `match_tree_data()` named-
vector column `x` is preserved.

## Clean-install smoke

The authoritative tarball was installed into a new empty R library. A second,
fresh R session loaded the package and ran the minimum requested workflow:

- `fast_k()`
- `fast_lambda()`
- `fast_d()`
- `fast_delta()`
- `fast_ace()`
- `check_tree()`
- `fast_signal()`

Result: `CLEAN_INSTALL_SMOKE = PASS`.

## Superseded and historical evidence

Every prior tarball/check record is non-authoritative and retained only as
history:

| Status | Context | testthat | Tarball SHA-256 | Check evidence |
|---|---|---:|---|---|
| `SUPERSEDED` | Post-cleanup, before final contract correction | 6,790 | `2C1D569781436394E604BB9CF7291364EEAFDDD81A7B93D9A3BA5A28C35C1DA5` | `7F683F78333D2A5CF7CC63A47118C8F24959698316B2B4F51958A557473DE6A0` |
| `SUPERSEDED` | Minimal release patch | 6,787 | `D92D790851657B36D9EA5A9B5AA0DCD817A702F65CAA960A0D12EE0D0FDCA2B5` | `DB1AD001BDB047755845CBF1778D7758C16151BF98A0E57FA2E26C68DE19B335` |
| `SUPERSEDED_HISTORICAL` | Current-R pre-patch | 6,777 | `3C6A50FE5319206B4389E85167871569A9265E6C255BCF6B7DB08F5C0FB0BC1F` | historical local check |
| `SUPERSEDED_HISTORICAL` | R 4.1.0 compatibility reference | 6,777 | `B52FB5E391DF03695D72B07C0B5496CBF2266AE407C66DCB8CFEBDDEDE8BF722` | historical local check |
| `SUPERSEDED_HISTORICAL` | Earlier current-R artifact | not retained | `8BFFF70A59A1C6A59E4705AB1940DBC31930B9734BBF2A795BB4331718693BF5` | historical only |
| `SUPERSEDED_HISTORICAL` | Documentation refresh artifact | not retained | `057D024531C535F367F878DEFD7C5116F85E0E46F98B7807E4CEA45A9AE4EF4C` | historical only |
| `SUPERSEDED_UNREFERENCED` | Existing temporary tarball, mtime 2026-08-24 03:08 | not established | `E43D683C981B93277F04A062B70B51A3582DFA1172C65B87893F64FA457D60D2` | no matching authoritative log |

Historical pre-patch check-log SHA-256:
`22FF9CEAE0F7203B17C4C835CA7122538E16D41B2173DE2F980FD7FBF5DEC53A`.
Historical pre-patch test-output SHA-256:
`DEF40CDEC72AD744CCA37F415F12B371693ACBD39F0666AD339412F5D8A85AD6`.
The former source snapshot
`EF89EDAE96B32A0DCEB25260877727F9E93C9BC75CF9215313DFD4E4920BE191`
is `SUPERSEDED_HISTORICAL` and is not the current source snapshot.

No historical evidence was deleted. Historical PASS results are not promoted
to current platform claims.

## Rejected and deferred work

Candidate 4A remains rejected. Its tolerance must not be relaxed and its code
must not be integrated into production.

The following remain deferred and are not PASS:

- CRAN submission checks
- Linux and macOS CI
- R-devel
- alternate BLAS
- no-OpenMP build

No Candidate 5, ACE/Delta optimization, parser refactor, D traversal refactor,
many-NA grid expansion, or large performance benchmark is part of this freeze.

## Git provenance

The frozen source is committed with message `fastphylosig 0.1.0 release` and
tagged `v0.1.0`. The Git commit and tag are the provenance record. This file
does not embed its own commit hash.
