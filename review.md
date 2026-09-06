# fastphylosig 0.2.0 Stage 2B2B Candidate 2 Expert Review

## Decision

`PREPARATION_ENGINE_STAGE1 = ACCEPTED_FROZEN`

`RAW_ROUTE_SINGLE_PREPARATION_BOUNDARY = ACCEPTED_FROZEN`

`STAGE2B2A_COMPLEXITY_ROOT_CAUSE = CONFIRMED`

`INSPECTION_INTEGER_ADJACENCY = ACCEPTED_FROZEN`

`CANDIDATE_2_PROOF_GATE = FAIL`

`CANDIDATE_2 = REJECTED`

No Candidate 2 production code was written. The required exact-ordering proof
failed before prototype integration: the current contract is a locale-sensitive
comparison of delimiter-collapsed strings, not an ordering of descendant sets.
The authorized compact set/rank candidates cannot preserve that contract for
all currently valid labels. Candidate 1 remains the frozen production state.

## 1. Candidate-1 Post-Freeze Baseline

The baseline used Candidate 1 production commit
`1685d8f58dea6dfe0a3e3d519d58581b3bbae5b9`, package version `0.2.0.9000`,
and R 4.6.1 UCRT on Windows. Fixed balanced, seeded-random, and pectinate
fixtures were timed in one serialized R process. Calls used warmup, alternating
serialized fixture copies and helper order, and median/IQR reporting. There
were 144 successful raw timings and no failed cells.

Times below are median seconds with IQR in brackets.

| Tips | Helper | Balanced | Random | Pectinate |
|---:|---|---:|---:|---:|
| 2,000 | descendant keys | 0.170 [0.020] | 0.170 [0.020] | 0.470 [0.010] |
| 2,000 | canonicalization | 0.330 [0.010] | 0.330 [~0] | 0.630 [0.010] |
| 2,000 | `prepare_tree()` | 0.340 [0.010] | 0.370 [0.010] | 0.660 [~0] |
| 5,000 | descendant keys | 0.720 [0.010] | 0.720 [0.020] | 2.370 [0.010] |
| 5,000 | canonicalization | 1.360 [0.040] | 1.390 [0.010] | 2.970 [0.320] |
| 5,000 | `prepare_tree()` | 1.410 [0.030] | 1.450 [0.010] | 3.110 [0.100] |
| 10,000 | descendant keys | 2.330 [0.065] | 2.390 [0.095] | 10.580 [0.015] |
| 10,000 | canonicalization | 4.560 [0.025] | 4.750 [0.200] | 12.820 [0.100] |
| 10,000 | `prepare_tree()` | 4.720 [0.105] | 4.910 [0.055] | 13.190 [0.945] |
| 20,000 | descendant keys | 7.940 [0.080] | 8.430 [0.170] | 42.020 [1.775] |
| 20,000 | canonicalization | 15.950 [0.515] | 18.160 [0.430] | 51.360 [0.355] |
| 20,000 | `prepare_tree()` | 16.660 [0.050] | 18.970 [0.195] | 52.940 [0.250] |

For pectinate trees from 10,000 to 20,000 tips, the descriptive empirical
exponents were 1.990 for descendant keys, 2.002 for canonicalization, and
2.005 for `prepare_tree()`. These are empirical diagnostics, not mathematical
proofs.

## 2. Exact Old Ordering Contract

The contract was extracted from `R/analysis_preparation.R`, not inferred from
function names:

1. `.edge_children_index()` preserves the child order induced by edge rows.
2. A tip contributes its exact `tip.label` string.
3. An empty internal node contributes `!empty:<node>`.
4. Cycle fallback contributes `!cycle:<node>` through the existing iterative
   descendant traversal.
5. Every internal node concatenates all child memo values, applies base R
   `sort(values)` with its current defaults, and stores the complete sorted
   descendant-label vector.
6. `.descendant_keys_iterative()` collapses that vector using the literal
   separator `"\r"`.
7. `.safe_canonicalize_core()` orders non-root internal nodes with
   `order(keys, remaining)`: the collapsed key is primary and the original
   internal node ID is the numeric tie-break.
8. The structural root is always placed first.
9. All children of a polytomy participate. Child order is normally erased by
   the final label sort on valid trees, but remains part of traversal and
   failure behavior.
10. Neither `sort()` nor `order()` fixes a method, encoding, or collation.
    Ordering therefore inherits the active base R `LC_COLLATE` behavior.

The package currently accepts unique, non-empty labels containing punctuation
and control separators. It does not forbid an embedded `"\r"`.

## 3. Exact Representation Assessment

No bounded compact representation met the proof gate:

- Minimum descendant rank is insufficient: sets such as `{a,b}` and `{a,c}`
  share the same minimum but have different old keys.
- Any fixed-length prefix can be defeated by clades sharing that prefix and
  differing later.
- A complete global-rank vector retains the same
  `sum_v descendant_count(v)` payload and therefore fails the complexity gate.
- A DFS interval is invalid because label order and topology order are
  independent; a clade need not be contiguous in globally sorted label order.
- A probabilistic hash cannot be the ordering truth under the explicit task
  contract.
- A persistent rope or sparse set might reduce storage, but no bounded design
  was found that proves exact equivalence to base R's locale-sensitive
  comparison of the complete collapsed strings. It therefore cannot enter an
  exact-equivalence prototype as an asserted solution.

## 4. Collision And Ordering-Drift Proof

The collapse operation is not an injective encoding of descendant labels. In
the C locale, two valid binary clades can contain these four distinct labels:

```text
clade A: "a",    "b\rc"
clade B: "a\rb", "c"
```

Both old keys are exactly `"a\rb\rc"`. The old implementation therefore uses
the original internal node ID to break the tie. A collision-free descendant-
set representation distinguishes the clades and can reverse that tie, causing
canonical mapping, edge order, node labels, and metadata to drift. Treating
the sets as equal without reproducing the exact collapsed string would instead
make the proposed representation non-exact.

This contract was also executed against the frozen source with
`benchmarks/stage2b2b2/verify_ordering_contract.R`. The smoke confirmed a valid
rooted binary tree, identical collapsed keys, numeric node-ID tie-breaking,
safe canonicalization, and serialized input immutability. The machine-readable
result is retained in `benchmarks/stage2b2b2/results/contract/`.

Locale adds a separate obstruction: ordering label tokens independently does
not prove the same result as collating their complete separator-joined string.
Prefix-like, punctuation, accent/case, encoding, and separator-adjacent labels
may be treated differently by the active collation. The old locale is not
fixed or included in a canonical fingerprint.

## 5. Exact-Equivalence Result

`CANDIDATE_2_PROTOTYPE = NOT_INTEGRATED`

`EXACT_ORDERING_PROOF = FAIL`

The proof gate failed before a candidate could legitimately be presented as an
exact implementation. Accordingly, no production binding was replaced and no
old oracle was removed. Equality requirements, label rules, canonical
signature semantics, child ordering, and failure semantics were not relaxed.

The required adversarial fixture plan remains documented for any future,
independently authorized design: balanced/random/pectinate/polytomy/two-tip,
edge shuffling, safe internal renumbering, alternate root numbering, exact
branch-row association, prefix/numeric/punctuation/separator labels, delimiter
collisions, duplicate-name failure, malformed topology, invalid `Nnode`, and
invalid root representation.

## 6. Pectinate Payload Operations

For the current pectinate fixture, `.descendant_keys_iterative()` excludes the
root and stores descendant payload for clades of size 2 through `n - 1`:

```text
old payload = 2 + 3 + ... + (n - 1) = n(n - 1)/2 - 1
```

| Tips | Old descendant-label payload visits | Candidate 2 |
|---:|---:|---:|
| 5,000 | 12,497,499 | not implemented |
| 10,000 | 49,994,999 | not implemented |
| 20,000 | 199,989,999 | not implemented |

No new `n x n` matrix, bitset, probabilistic ordering hash, or quadratic
persistent structure was introduced.

## 7. Canonicalization Before And After

| Shape, 20,000 tips | Before, median [IQR] | After |
|---|---:|---:|
| Balanced | 15.950 [0.515] s | not run: exact gate failed |
| Random | 18.160 [0.430] s | not run: exact gate failed |
| Pectinate | 51.360 [0.355] s | not run: exact gate failed |

The protocol explicitly authorizes formal after timing only after exact
equivalence passes. Reporting an after value here would violate that gate.

## 8. prepare_tree() Before And After

| Shape, 20,000 tips | Before, median [IQR] | After |
|---|---:|---:|
| Balanced | 16.660 [0.050] s | not run: no production candidate |
| Random | 18.970 [0.195] s | not run: no production candidate |
| Pectinate | 52.940 [0.250] s | not run: no production candidate |

The required pectinate 1.5x or 30% production acceptance threshold was not
evaluated because correctness and complexity are prerequisite gates.

## 9. Small-Tree Slowdown Gate

`SMALL_TREE_SLOWDOWN_GATE = NOT_RUN`

There is no Candidate 2 implementation to compare. The Candidate 1 frozen
small-tree evidence remains PASS and was not rerun or relabeled as Candidate 2
evidence.

## 10. Test And Package Check

No production or package test file changed in this stage, so the full suite and
package check were not repeated after the proof-gate rejection. The unchanged
Candidate 1 authoritative evidence remains:

- full `testthat`: 7,412 PASS, 0 FAIL, 0 WARN, 0 SKIP;
- `R CMD check --no-manual --timings`: Status: OK;
- package check: 0 ERROR, 0 WARNING, 0 NOTE.

These results are historical frozen-production evidence, not falsely labeled
as a Candidate 2 integration check.

## 11. Final Status

`CANDIDATE_2 = REJECTED`

Reason: no proposed smaller representation simultaneously preserved the full
delimiter-collapsed, locale-sensitive old ordering contract and removed the
quadratic descendant-label payload. The rejection occurred before production
integration, exactly as required by the correctness-first gate.

Candidate 1 remains unchanged and frozen. Stage 2C was not started.

Formal baseline data and provenance are retained under
`benchmarks/stage2b2b2/results/baseline/`; executable contract evidence is
retained under `benchmarks/stage2b2b2/results/contract/`.
