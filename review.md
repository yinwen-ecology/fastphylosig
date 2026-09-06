# fastphylosig 0.2.0 Stage 2B2C Expert Review

## Decision

This was a read-only correctness-contract audit of production source commit
`654eb8e` with fastphylosig `0.2.0.9000` and R 4.6.1 UCRT on Windows.
Production code, public APIs, estimator definitions, RNG policy, threading,
and existing test expectations were not changed.

`PREPARATION_ENGINE_STAGE1 = ACCEPTED_FROZEN`

`RAW_ROUTE_SINGLE_PREPARATION_BOUNDARY = ACCEPTED_FROZEN`

`INSPECTION_INTEGER_ADJACENCY = ACCEPTED_FROZEN`

`CANDIDATE_2 = REJECTED`

`CANONICALIZATION_CONTRACT_BUG = YES`

`PUBLIC_RELEASE_BLOCKER_0_2_0 = YES`

`NEXT_STEP = CANONICALIZATION_CONTRACT_V2`

## 1. Delimiter-Collision Fixture

The audit used four unique, non-empty labels that are legal under the current
input contract:

```text
clade A: "a",    "b\rc"
clade B: "a\rb", "c"
```

Two rooted binary trees were constructed with identical topology, branch
lengths, tip labels, and tip-to-tip patristic distances, but different internal
node numbering. A third tree retained the same biology and changed only edge
row order while preserving each edge-length association. All fixtures were
accepted as rooted, structurally valid trees.

Under the frozen implementation, both non-root clades produce the identical
descendant key `"a\rb\rc"`. The collision is therefore demonstrated with valid
inputs; it is not based on duplicate or empty labels.

`DELIMITER_COLLISION_FIXTURE = CONFIRMED`

Evidence: `benchmarks/stage2b2c/results/canonical/delimiter_collision_keys.csv`
and `representation_comparisons.csv`.

## 2. Internal Renumbering Invariance

`CANONICAL_RENUMBERING_INVARIANCE = FAIL`

For the biologically equivalent internally renumbered tree, `identical()` was
false for:

- canonical edge matrix;
- canonical edge-length vector and canonical edge-length association;
- before and after canonical signatures;
- before and after tree fingerprints.

Tip labels and `Nnode` remained identical, and both canonicalizations reported
their safety checks as passing. Numeric mapping metadata was identical, but it
did not identify that the tied numeric node IDs referred to different clades.
The original-node-ID tie break is therefore representation dependent when
descendant keys collide.

## 3. Edge-Order Invariance

`CANONICAL_EDGE_ORDER_INVARIANCE = FAIL`

For edge rows reordered together with their branch lengths, biological
equivalence and canonical edge-length association were preserved. However,
the canonical edge matrix, edge-length sequence, after signature, and final
fingerprint were not identical. The current canonicalization is safe in the
narrow sense that it preserves the weighted biological tree, but it does not
produce one edge-row-independent canonical object.

This distinction matters because the public fingerprint encodes edge rows and
edge lengths in their current order.

## 4. Locale Dependence

`CANONICALIZATION_LOCALE_DEPENDENT = TRUE`

The current normal `LC_COLLATE` was `C`. The audit also successfully switched
only `LC_COLLATE` to these supported locales:

- `English_United States.1252`;
- `German_Germany.1252`;
- `Chinese (Simplified)_China.936`;
- `en_US.UTF-8`.

The fixture included lowercase and uppercase ASCII, prefixes, accented labels,
punctuation, numeric-looking labels, and separator-containing labels. Relative
internal ordering changed from:

```text
C:       19,17,13,18,16,14,15
non-C:   19,17,16,13,18,14,15
```

Canonical edges and canonical fingerprints differed from the C-locale result.
A separate deterministic, fully binary fixture (`seed = 20260906`) confirmed
that the canonical fingerprint differs between C and Windows English while
remaining valid for every audited estimator.

Evidence: `locale_comparisons.csv` and
`results/downstream/locale_fixture.csv`.

## 5. Fingerprint Propagation

Independently prepared, representation-equivalent trees had different public
fingerprints even though their readiness status and retained-tip cache key were
the same. Representation identity is therefore fragmented across internal
renumbering and edge-row order.

More importantly, the delimiter encoding creates a direct integrity failure.
Starting from a prepared context, the audit changed the protected tip-label
vector from:

```text
c("a", "b\rc", "a\rb", "c")
```

to another unique, non-empty vector:

```text
c("a\rb", "c", "a", "b\rc")
```

The concatenated fingerprint was identical. Prepared-context validation
returned successfully instead of rejecting the mutation. The stale prepared
cache then produced:

```text
prepared K:     1.41323926109659
raw mutated K:  1.22889253471746
absolute diff:  0.184346726379132
relative diff:  0.130442686849847
```

`FINGERPRINT_MUTATION_COLLISION = TRUE`

`FINGERPRINT_MUTATION_REJECTED = FALSE`

`STALE_CACHE_NUMERICAL_PROPAGATION = TRUE`

Evidence: `fingerprint_mutation_collision.csv`, `prepared_contexts.csv`, and
`fingerprint_cache_identity.csv`.

## 6. K, Lambda, D, Delta, And ACE Propagation

For the unmutated representation-equivalent collision fixtures, 753 comparison
rows covered raw and prepared routes, controlled K permutations, lambda
likelihood/LR, controlled D random and Brownian nulls, fixed-seed serial Delta,
and ACE with CI. All comparisons passed with no comparison errors:

| Method | Compared outputs | Maximum absolute difference | Maximum relative difference |
|---|---|---:|---:|
| K | estimate, P, MCSE, status, retained species | 0 | 0 |
| lambda | estimate, logLik, LR, P, status, retained species | 0 | 0 |
| D | estimate, null P/MCSE/accounting, status, retained species | 0 | 0 |
| Delta | estimate, P/MCSE/diagnostics, status, retained species | 0 | 0 |
| ACE | logLik, aligned ancestral likelihoods, rates, SE, status | 0 | 0 |

The locale-sensitive fully binary fixture was then run under C and Windows
English using identical permutations, Brownian states, Delta seed, and one
worker. `check_tree()`, raw/prepared K, lambda, D, Delta, and ACE all completed;
their reported estimates, likelihood fields, P/MCSE fields, statuses, warning
classes, failure semantics, and retained-species fields were exactly equal.
The canonical fingerprint itself still differed across locales, as recorded by
the fixture provenance.

Thus ordinary independent analyses are numerically invariant on these fixed
fixtures. The user-visible numerical failure occurs when the fingerprint
collision permits a mutated prepared context to reuse stale cached evidence.

Evidence: `results/downstream/propagation_comparisons.csv` and
`locale_propagation.csv`.

## 7. Documentation Contract Audit

`DOCUMENTATION_CONTRACT_MISMATCH = CONFIRMED`

| Claim | Location | Classification | Reason |
|---|---|---|---|
| Normalization changes representation only | `README.md:74-76,134-138`; `USAGE_zh.md:64-68,115-117`; `NEWS.md:18-20` | SUPPORTED, narrow meaning | Weighted topology, tips, root, polytomies, and branch lengths are preserved. |
| Edge order and internal IDs are canonicalized/standardized | `NEWS.md:54-56`; `man/resolve_tree.Rd:22-33`; `man/fast_d.Rd:13-16` | TOO_STRONG if read as a unique canonical form | Equivalent representations can retain different canonical edges and fingerprints. |
| Protected tree-field changes are rejected on context reuse | `README.md:127-130`; `USAGE_zh.md:107-110`; `man/prepare_tree.Rd:33-36` | TOO_STRONG | A valid separator-label mutation preserved the fingerprint and was accepted. |
| Prepared cache is invalidation-safe | `R/prepare_tree.R:5-8` | TOO_STRONG | The demonstrated collision reused stale computational evidence. |
| C locale does not change package behavior | `inst/RC_READINESS.md:42-43` | TOO_STRONG | Canonical internal ordering and canonical fingerprints changed by locale. |
| Cross-platform qualification remains deferred | `README.md:228-234`; `docs/index.md:121-124` | SUPPORTED | The package does not claim completed cross-platform qualification there. |
| Representation-safe canonicalization preserves method contracts | `R/analysis_preparation.R:553-556` | NOT FULLY SUPPORTED | Estimators were stable for independent inputs, but mutation safety failed. |

Documentation narrowing alone would not repair the stale-cache numerical
failure, so this is not only Case C.

## 8. Contract-Bug Decision

`CANONICALIZATION_CONTRACT_BUG = YES`

The findings satisfy Case B. Internal numbering and edge order can produce
different canonical/fingerprint identities under the current representation-
safe wording, locale changes canonical ordering, and the same delimiter design
allows a protected-field mutation to evade integrity validation. The observed
stale-cache numerical propagation is a correctness defect, not merely an
internal performance limitation.

Candidate 2 remains rejected. Its compact descendant-key proposal is neither
revived nor used as a remedy.

## 9. Release-Blocker Decision

`PUBLIC_RELEASE_BLOCKER_0_2_0 = YES`

Version 0.2.0 should not proceed to release or Stage 2C optimization while a
documented prepared-context integrity check can accept a computationally
relevant tree mutation and return a stale numerical result. The accepted Stage
1, Stage 2B1, and Candidate 1 work remains frozen; the blocker is confined to
the canonicalization/fingerprint contract exposed by this audit.

No full test suite, package check, or performance benchmark was run because
this stage changed no production or existing test code. Both audit scripts
completed with R exit status 0 on their final valid fixtures.

## 10. Required Next Step

Proceed with a separately designed `CANONICALIZATION_CONTRACT_V2`, not Stage
2C. Requirements are:

- collision-free ordering truth;
- locale-independent comparison;
- independence from original internal node IDs;
- independence from edge-row order;
- deterministic output;
- no quadratic descendant-label payload;
- no probabilistic hash as ordering truth;
- unambiguous, domain-separated or length-delimited fingerprint encoding;
- complete mutation detection for all computationally relevant protected
  fields under the documented contract;
- input immutability;
- unchanged K, lambda, D, Delta, and ACE estimator definitions;
- exact fixtures for delimiter labels, locale switches, internal renumbering,
  edge-row shuffling, branch-row association, and stale-cache rejection.

Stage 2C is not authorized. Stop after this audit.
