# fastphylosig 0.2.0 Final Release Expert Review

Review date: 2026-09-13

> **Superseding audit note.** The artifact and gate results below were the
> preceding release-hardening record. They are superseded for the current
> local-release decision by
> `benchmarks/final_release/release/FINAL_LOCAL_RELEASE_EVIDENCE_AUDIT.md`:
> the exact-tarball rerun was environment-blocked and stale `0.2.0.9000`
> release wording was found in the public documentation. The current status
> is `NOT_READY_PENDING_FINAL_LOCAL_AUDIT`; no production code
> was changed.

## 1. Final commit and version

`fastphylosig 0.2.0` is frozen at release source commit
`c084ced0c9e44966369f6b0106d185d40909237b` on
`codex/fastphylosig-0.2.0-dev`. The source snapshot used for the package is a
`git archive` of that commit. No executable production change was made after
the V2 acceptance baseline; the final transition changed the package version
only.

## 2. API audit

`PUBLIC_API_AUDIT = PASS`. The documented exported surface is the intended
13-entry surface: `fast_signal`, `fast_k`, `fast_lambda`, `fast_d`,
`fast_delta`, `fast_ace`, `check_tree`, `resolve_tree`, `match_tree_data`,
`prepare_tree`, `cache_info`, `plot_signal`, and deprecated
`match_phylo_data`. No unintended `backend`, `engine`, `unsafe`, `trust`,
`skip_check`, `fast_sig`, `fast_delta2`, or `fast_ace2` public entry was found.

## 3. Scientific parity

`SCIENTIFIC_PARITY = PASS`. The final installed tarball passed 22/22 focused
cases with 0 failure, 0 skip, and 0 unexpected warning. K, lambda, D, Delta,
and ACE estimates and diagnostics were checked against their frozen contracts,
including MCSE/P accounting, log likelihood/LR, ancestral likelihoods, strict
tails, polytomy behavior, and warning/status semantics. The deliberately
ill-conditioned ACE fixture is retained as an expected explicit
non-convergence case; no warning was suppressed.

## 4. Raw/prepared parity

`RAW_PREPARED_PARITY = PASS` for normal, NA, large, and representation-shuffled
fixtures across K, lambda, D, Delta, and ACE. Input objects remain immutable.
Protected prepared-context mutations are rejected before cache or numerical
work.

## 5. V2 integrity and persistence

`CANONICALIZATION_V2 = PASS` and `MUTATION_SAFETY = PASS`. Internal-node
renumbering, edge-row order, alternate root representation, delimiter and
Unicode labels, and heterogeneous branch lengths preserved canonical branch
association and fingerprints in the tested contract. `saveRDS()`/`readRDS()`
reuse passed in a fresh R session; mutation after read was rejected, and V1 or
unknown schemas requested a rebuild.

## 6. Reproducibility and parallel audit

`REPRODUCIBILITY_AUDIT = PASS` for same seed, same arguments, and same worker
configuration on representative K, D, and Delta fixtures. No public promise
of identical trajectories across different `ncores` was added. D documentation
correctly cautions that small parallel jobs may be slower; it does not claim
that more cores are always faster.

## 7. Full tests

The final exact-source run recorded `8336 PASS / 0 FAIL / 0 WARN / 0 SKIP` in
`benchmarks/final_release/release/testthat.Rout`.

## 8. R CMD check

`R CMD check --no-manual --timings fastphylosig_0.2.0.tar.gz` returned
`Status: OK` with 0 ERROR, 0 WARNING, and 0 NOTE. The check used R 4.6.1 and
the final C++17 build. An earlier environment-only run with unsupported
`C.UTF-8` startup variables is retained separately and is not the authoritative
check log.

## 9. Clean-install smoke

`CLEAN_INSTALL_SMOKE = PASS` from a fresh private library and fresh R session.
Both raw and prepared paths exercised `check_tree`, `prepare_tree`,
`fast_k`, `fast_lambda`, `fast_d`, `fast_delta`, `fast_ace`, and
`fast_signal`, with package version `0.2.0`.

## 10. Cross-platform status

`CROSS_PLATFORM_QUALIFICATION = NOT_RUN / DEFERRED` by explicit instruction.
Linux, macOS, R-devel, alternate BLAS/LAPACK, no-OpenMP fallback, and CI
qualification are not represented as PASS. This local release review does not
make an all-platform claim.

## 11. Final public benchmark summary

`FINAL_BENCHMARK_CORRECTNESS = PASS`: 126 workload cells, 10 observations per
cell, 1260 non-reference timing rows, 9/9 guards, and 0 timing warnings. The
serialized benchmark covered raw single-trait, prepared repeated, matrix
batch, permutation-heavy K/D, and large-tree routes at n=500--5000. The
benchmark was run at the pre-version commit `084e979...` with package metadata
`0.2.0.9000`; the final commit differs only in `DESCRIPTION`'s version field,
so the measured production code is identical. `phytools` guards passed; the
optional `caper` call was recorded as unavailable and was not used for speed
claims.

## 12. Memory and package-size summary

Prepared-context memory characterization passed at n=500, 1000, 5000, 10000,
and 20000. Measured `object.size` ranged from 1,375,896 to 53,654,896 bytes;
uncompressed RDS size ranged from 1,017,237 to 39,878,633 bytes. These are
object/serialization measurements, not peak RSS. The authoritative tarball is
231,978 bytes with 106 entries and contains no benchmark, manuscript, docs,
Rcheck, Rout, or temporary-log payload.

## 13. Documentation audit

`DOCUMENTATION_CONTRACT_AUDIT = PASS` for the API and V2 contracts. The
superseding local evidence audit found a separate release-version consistency
blocker: README, NEWS, and the website index still contain `0.2.0.9000`
development wording while the artifact is `0.2.0`. The prepared ACE wording
was corrected to the explicit `fast_ace(x, prepared = ctx)` form.

## 14. Final hashes

| Artifact | SHA-256 |
|---|---|
| `fastphylosig_0.2.0.tar.gz` | `c3f05d59055c72eec6a9765cdd9a68803ab325583d72fc7d221ec25f028bba63` |
| source snapshot archive | `a20e5d2065565b4eae6c1f312f8c249df034b69bd72a95701b957f53d95980db` |
| final `00check.log` | `ec9ccdbfcdf106de3f40d074861f879869669d1ab8f17bca97313e2e28650318` |
| final `testthat.Rout` | `f4bfd17577aef8a98e6f05271a53c4cbbf708498f676fc305d4d7c024f4ef73f` |
| scientific validation log | `607ca7f0576eca98667791ac4b8a20ad76d84801a54e3f2de6c7ebc96e657289` |
| final benchmark timing CSV | `9a547fa1990adfe60d416f17535cf2bd4f5470fa23b4e69cd0a441fb5f3c8f41` |
| final benchmark memory CSV | `1ae6810a3b2da3274fb4e719f012336343fa7e62b04f6815c3097b1ef668e36a` |
| final benchmark status CSV | `8da422a10fb016ecd39c7a2159dea6d94bd68cae7cb7fb997f277e6dc5fc011f` |

## 15. Unresolved deferred qualification

The preceding hardening run reported no release-blocking local issue. The
superseding final-local audit now tracks the release-version wording blocker
and environment-blocked fresh R gates. Other platforms/toolchains, alternate
BLAS/LAPACK, no-OpenMP fallback, and future optimizer robustness remain
deferred. The historical 0.1.0 tarball remains `SUPERSEDED/HISTORICAL` and was
not deleted.

## 16. Release decision

```text
FASTPHYLOSIG_0_2_0 = NOT_READY_PENDING_FINAL_LOCAL_AUDIT
PRIOR_LOCAL_RELEASE_EVIDENCE = RELEASE_READY (SUPERSEDED)
RELEASE_SCOPE = LOCAL_WINDOWS_R4.6.1
COMMON_PREPARATION_ENGINE_0_2_0 = FROZEN
CANONICALIZATION_CONTRACT_V2 = ACCEPTED_FROZEN
K_OPTIMIZATION_0_2_0 = CLOSED_FINAL
D_OPTIMIZATION_0_2_0 = CLOSED_FINAL
LAMBDA_OPTIMIZATION_0_2_0 = CLOSED
DELTA_OPTIMIZATION_0_2_0 = CLOSED
```

This review is the final 0.2.0 release-hardening record. No 0.2.1 or 0.3.0
work is started by this task.
