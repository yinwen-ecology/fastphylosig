# fastphylosig 0.2.0 Document Consistency Audit

Audit date: 2026-09-13

## Scope

This audit compares the authoritative source commit
`3f84775bc96e00820aa39522245cb63175f79a5d`, the rebuilt tarball, public guides,
Rd files, `DESCRIPTION`, and `NAMESPACE`. It does not qualify other platforms.

## Results

| Area | Result | Evidence |
|---|---|---|
| Package version | PASS | Source and tarball report 0.2.0 |
| R version contract | PASS | `R (>= 4.1.0)` |
| Public API | PASS | Exactly 13 exports; current names agree across guides and Rd |
| Method signatures/defaults | PASS | R formals and Rd usage agree |
| Release wording | PASS | README, NEWS, docs index, and USAGE_zh contain no stale 0.2.0.9000 development wording |
| Dependencies | PASS | `mvtnorm` is declared in Suggests for its test-only namespace call |
| V2/mutation/persistence contract | PASS | Claims remain scoped to the implemented protected-field and schema behavior |
| RNG/thread wording | PASS | No cross-core universal determinism or always-faster parallel claim |
| Benchmark provenance | PASS | 084e979/0.2.0.9000 is explicitly historical, not relabeled |
| Platform boundary | PASS | Cross-platform and CRAN qualification are not claimed |

The lambda summary tables describe the available likelihood/LR capability;
the formal signatures and examples correctly retain `test = FALSE` as the
default. This is not treated as a claim that LR fields are always computed.

## Tarball check

The rebuilt package embeds the corrected README, NEWS, DESCRIPTION, and
dependency declaration. A content scan found no user-machine path, `.codex`
path, `sourceCpp`, debugger/browser hook, working-directory mutation, locale
mutation, or audit-library variable.

```text
DOCUMENT_CONSISTENCY = PASS
DOCUMENTATION_RELEASE_BLOCKER = NO
```
