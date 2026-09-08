# fastphylosig 0.2.0 Canonicalization Contract V2-B Expert Review

Package state: `0.2.0.9000` on `codex/fastphylosig-0.2.0-dev`  
Production implementation commits: `7888f38`, `34ae5de`  
Comparison baseline: `7d0d6e38756f271e48a74a4ac890f64b186674f3`  
Local qualification: R 4.6.1 UCRT, Windows 11 x64, GCC 14.3.0,
Rtools 4.5, OpenMP enabled

## 1. Persistence Decision

`V2_CONTEXT_PERSISTENCE_DECISION = SELF_CONTAINED_EXACT_SNAPSHOT_SUPPORT`.
A valid V2 context may be saved with `saveRDS()` and reused after `readRDS()`
in a fresh R process. Fresh-process K and lambda reuse and mutation-after-load
rejection passed. The decision and pre-integration audit are recorded in
`benchmarks/stage2b2c_v2b/V2_CONTEXT_PERSISTENCE_DECISION.md`.

## 2. Production Canonicalization

V2 is the production canonicalizer. It preserves public tip IDs and the stored
`tip.label`, assigns the root to `n_tip + 1`, orders children by the exact
unsigned UTF-8 byte order of each subtree's minimum descendant tip, numbers
remaining internal nodes in canonical preorder, and emits edges in canonical
postorder. Branch lengths follow the original biological child edge through
endpoint remapping. Original internal-node IDs, original edge-row positions,
locale collation, delimiter strings, and probabilistic hashes are not ordering
truth. Candidate 2 remains rejected and was not revived.

## 3. Fingerprint And Protected Snapshot

The canonical fingerprint and protected source snapshot are separate exact,
versioned structured encodings. Both cover the type, shape, order, and content
of `tip.label`, `Nnode`, `edge`, and `edge.length`; branch values use exact
binary64 bytes. A locked package-owned integrity record holds the authoritative
schema, snapshot, fingerprint, fixed metadata, and cache identities. External
prepared-context reuse validates this evidence before cache lookup or a
numerical kernel. Delimiter collapse and probabilistic hash collisions are not
part of the V2 integrity decision.

## 4. Old Context Handling

Missing, V1, altered, or unknown context schema markers are rejected with a
clear instruction to run `prepare_tree()` again. V2 does not silently migrate
old contexts, accept a legacy fingerprint, reconstruct evidence from mutable
public fields, or reuse an old cache. A raw tree remains the supported route
for rebuilding a context.

## 5. Representation Invariance

Exact gates passed for balanced, random, pectinate, polytomous, two-tip,
delimiter-containing, prefix, punctuation, accented, Unicode, and
heterogeneous-branch fixtures. Safe internal-node renumbering, shuffled edge
rows, alternate valid root numbering, and combined representation changes
produced identical canonical computational trees and fingerprints. Input trees
and traits remained immutable. `node.label` is non-computational metadata and
may be remapped with internal nodes; it is not part of the integrity contract.

## 6. Locale Audit

Canonical trees and fingerprints were identical to the C-locale oracle under
all six available test cases: current C, explicit C, English US CP1252, German
CP1252, Chinese Simplified CP936, and `en_US.UTF-8`. Therefore
`V2_LOCALE_INVARIANCE = PASS` for the tested host and locales. This is not a
claim of qualification for every operating system, R build, or locale.

## 7. Mutation And Cache Safety

Mutations of `tip.label`, `edge`, `edge.length`, or `Nnode` reject before cache
or estimator use. Tests also reject synchronized tampering with public
fingerprint/snapshot fields, protected fixed metadata, cache replacement,
missing integrity records, and changed schema markers. Cache identity and the
current exact snapshot are checked against the locked record on every external
public boundary. Valid V2 serialization recreates a self-contained record;
stale-cache propagation is blocked by validation.

## 8. Estimator Parity

Raw/prepared and representation-equivalent parity passed for K, lambda, D,
Delta, and ACE, including controlled permutation or null inputs where methods
are stochastic. Checked fields include estimates, P values, MCSE, logLik/LR,
rates, ancestral likelihoods, statuses, warnings, failure semantics, and
retained species where applicable. ACE likelihood rows are restored to the
source internal-node order at the public boundary. No estimator definition,
numerical formula, tolerance, RNG stream, thread policy, or public result
structure changed.

## 9. Tests, Check, And Clean Smoke

Final testthat result: `8336 PASS / 0 FAIL / 0 WARN / 0 SKIP` in 59.9 seconds.
`R CMD check --no-manual --timings` on the built
`fastphylosig_0.2.0.9000.tar.gz` ended with `Status: OK`, hence 0 ERROR,
0 WARNING, and 0 NOTE. The authoritative check used `LC_ALL=C`; an earlier
attempt inherited an unsupported Windows `C.UTF-8` startup setting and stopped
at DESCRIPTION metadata, without a package-code failure. A fresh-library,
fresh-session tarball install loaded the package and ran `check_tree()`,
`prepare_tree()`, `fast_k()`, `fast_lambda()`, `fast_d()`, `fast_delta()`,
`fast_ace()`, and `fast_signal()` successfully. Evidence is in
`benchmarks/stage2b2c_v2b/formal/`.

## 10. Performance Evidence

Serialized alternating benchmarks used fixed balanced, random, and pectinate
fixtures; 10 repeats below 10,000 tips and 5 repeats at 10,000/20,000 tips.
At 20,000 tips V2 canonicalization took 1.36 s balanced, 1.47 s pectinate, and
1.72 s random, versus 17.25 s and 25.14 s for the two completed Candidate-1
comparators; the Candidate-1 pectinate comparator failed before a median.
V2 `prepare_tree()` took 1.97/2.92/2.90 s and raw K took 2.03/3.08/2.02 s for
balanced/pectinate/random, representing 8.98x to 10.67x completed preparation
speedups and 8.80x to 9.09x completed raw-K speedups at that size.

The formal five-observation prepared-K signal exceeded 25% on two 20,000-tip
shapes and therefore triggered a targeted confirmation. Twenty paired
observations with five inner iterations measured Candidate 1 versus V2 at
172 versus 209 ms for pectinate (+21.5%) and 168 versus 208 ms for random
(+23.8%). Neither confirmed regression exceeds the 25% release-pause rule.
Correctness was not rolled back. Timing evidence is tied to implementation
commit `34ae5de`; it is not a universal performance guarantee.

## 11. Complexity And Context Size

For a 20,000-tip pectinate tree, V2 recorded 39,998 child-tuple entries,
39,999 descriptor entries, 39,998 edges, and zero materialized descendant
payloads. The exact fingerprint and snapshot were 849,176 and 849,126 bytes;
the full prepared context was 53,654,896 bytes. Persistent descriptor and
adjacency entry counts are linear in tree size, with no n-by-n matrix,
quadratic descendant-label payload, or probabilistic identity table. These are
measured storage facts, not a proof of universal runtime complexity.

## 12. Documentation Audit

`README.md`, `USAGE_zh.md`, `NEWS.md`, `docs/index.md`, relevant Rd pages, and
`inst/POST_0.1.0_TECH_DEBT.md` now describe the V2 normalization, exact
integrity boundary, old-schema rejection, persistence support, branch
association, and tested-locale limitation. Public K/lambda documentation no
longer implies that production paths allocate dense covariance payloads.
Historical 0.1.0 release evidence in `inst/EXPERT_REVIEW.md`,
`inst/RC_READINESS.md`, and `docs/validation_manifest.md` remains historical
and is not presented as V2-B evidence.

## 13. Canonicalization Contract V2 Decision

`CANONICALIZATION_CONTRACT_V2 = ACCEPTED`.

All atomic production-integration gates passed: deterministic canonical form,
exact identity, exact mutation detection, schema rejection, persistence,
representation and locale fixtures, estimator parity, complete tests, package
check, clean install, performance threshold, and linear persistent evidence.

## 14. Public Release Blocker Decision

`PUBLIC_RELEASE_BLOCKER_0_2_0 = CLEARED`.

The delimiter-collision and representation-safety blocker identified in Stage
2B2C is removed by the accepted production V2 contract. This decision clears
that blocker only. Linux, macOS, R-devel, alternate BLAS, no-OpenMP, and other
deferred platform qualifications are not converted to PASS. Stage 2C remains
unauthorized and was not started.
