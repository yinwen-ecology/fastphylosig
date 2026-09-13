# fastphylosig 0.2.0 Document Consistency Audit

Audit date: 2026-09-13

## Scope

This is a read-only cross-document audit of the frozen `c084ced` source and
the authoritative `fastphylosig_0.2.0.tar.gz`. It does not perform
cross-platform qualification, performance work, or estimator changes.

## Results

| Area | Result | Evidence / disposition |
|---|---|---|
| Package version | PASS | `DESCRIPTION` and tarball `DESCRIPTION` both report `0.2.0`. |
| Public API names | PASS | README, `USAGE_zh`, website index, Rd, and `NAMESPACE` use the 13 intended top-level exports. |
| Method defaults and boundaries | PASS | K/lambda/D/Delta/ACE defaults and tree/trait contracts agree across Rd and guides. |
| RNG and worker wording | PASS | Guides describe fixed seeds/inputs/worker counts and warn that small parallel jobs may be slower; no universal speed claim. |
| V2 context and RDS contract | PASS | Protected fields, mutation behavior, schema rebuild rule, and `saveRDS()`/`readRDS()` wording agree. |
| Memory terminology | PASS | The benchmark uses `saveRDS(..., compress = FALSE)` and calls the reported value uncompressed RDS size. |
| Benchmark provenance | PASS | Development metadata `0.2.0.9000` at `084e979` is explicitly retained as code-equivalent evidence, not relabeled. |
| Platform boundary | PASS | Local Windows/R 4.6.1 scope is stated; cross-platform qualification is consistently deferred. |
| Release version wording | BLOCKER | `README.md:241`, `NEWS.md:1`, and `docs/index.md:135-137` still call the package `0.2.0.9000`/“development version” while the release artifact is `0.2.0`. Correct these strings, then rebuild and re-check the tarball. |

The version wording finding is documentation-only but affects files included in
the package archive (`README.md` and `NEWS.md`). It therefore cannot be fixed
without invalidating the current tarball and refreshing its build/check/hash
evidence.

## Static source scan

The extracted authoritative tarball has no hard-coded local path, benchmark
output, `.codex` path, `sourceCpp`, browser/debugger call, `setwd()`,
`Sys.setenv()`, `LC_ALL`, or `FASTPHYLOSIG_LIBRARY` dependency. A small
test-only oracle is named `stage2b2b` and is classified `DEVELOPMENT_ONLY`; it
contains no benchmark output or user-machine path. User-facing
`print()`/`cat()` calls and private implementation names such as `engine` and
`unsafe` are expected and are not public selectors or debug hooks.

## Required remediation

1. Change the stale release/development wording in README, NEWS, and the
   website index while preserving the documented API and scope.
2. Rebuild the package from the corrected source and refresh all artifact
   hashes and dynamic release evidence.
