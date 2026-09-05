# fastphylosig 0.2.0 Stage 2B2A Expert Review

## Decision

`PREPARATION_ENGINE_STAGE1 = ACCEPTED_FROZEN`

`RAW_ROUTE_SINGLE_PREPARATION_BOUNDARY = ACCEPTED_FROZEN`

`STAGE2B2A_PRODUCTION_CODE_CHANGE = NONE`

`STAGE2B2A_COMPLEXITY_ROOT_CAUSE = CONFIRMED`

`STAGE2B2B_IMPLEMENTATION = NOT_STARTED`

The current `prepare_tree()` behavior is near quadratic on this R/runtime, but
not because it performs a per-node full scan of the edge matrix. Two distinct
mechanisms explain the result: named-list child lookup plus repeatedly copied R
stacks in structural inspection, and full descendant-key materialization in
canonicalization. The latter becomes an additional dominant cost on pectinate
trees.

## Provenance And Protocol

- Audit source: `21dd4e522de125d253ed5884b3e1fcc47c693f2a`.
- Installed production source: `e265b47cd7abbb3c727ec466998c67e8eecdbc98`.
- There is no production/test diff between those commits in `R/`, `src/`,
  `DESCRIPTION`, `NAMESPACE`, `man/`, or `tests/`.
- Package version: `0.2.0.9000`.
- Runtime: R 4.6.1 UCRT, `x86_64-w64-mingw32`, Windows.
- Compiler: GCC/G++ 14.3.0, C++17.
- Fixed positive-length balanced, seeded random, and pectinate trees.
- Serialized execution: one non-overlapping R process per shape; no benchmark
  processes ran concurrently.
- Five formal repeats through n=5,000 and three at n=10,000/20,000, after
  warmup and alternating identical fixture copies. Every cell has median/IQR.
- A 120-second per-call guard was active; no formal observation was censored or
  failed.
- Direct helper timings are independent inclusive timings. They must not be
  added as if they were disjoint components.

## Tree-Shape Scaling

Values are median seconds with IQR in brackets.

| Tips | Balanced | Random | Pectinate |
|---:|---:|---:|---:|
| 500 | 0.070 [0.010] | 0.080 [0.020] | 0.110 [0.000] |
| 1,000 | 0.210 [0.020] | 0.200 [0.000] | 0.280 [0.010] |
| 2,000 | 0.600 [0.020] | 0.600 [0.010] | 0.940 [0.000] |
| 5,000 | 2.810 [0.080] | 2.830 [0.070] | 4.820 [0.070] |
| 10,000 | 10.100 [0.135] | 10.070 [0.315] | 20.390 [0.060] |
| 20,000 | 38.010 [1.970] | 37.330 [0.110] | 80.560 [1.490] |

The descriptive public-path exponents are:

| Shape | p, 5k to 10k | p, 10k to 20k |
|---|---:|---:|
| Balanced | 1.846 | 1.912 |
| Random | 1.831 | 1.890 |
| Pectinate | 2.081 | 1.982 |

These exponents diagnose observed scaling; they are not mathematical proofs.

## Major Helpers

At n=20,000:

| Shape | Inspection, s | p 10k-20k | Canonicalization, s | p 10k-20k | Signature, s |
|---|---:|---:|---:|---:|---:|
| Balanced | 20.85 | 1.960 | 16.47 | 1.837 | 3.85 |
| Random | 20.22 | 1.937 | 16.58 | 1.797 | 3.50 |
| Pectinate | 27.02 | 2.014 | 52.45 | 1.998 | 3.53 |

From 5k to 10k, inspection p is 2.000, 2.000, and 2.020 for
balanced, random, and pectinate trees. Canonicalization p is 1.709, 1.748, and
2.135 respectively. Thus inspection is near quadratic for every shape, while
canonicalization has an additional shape-specific quadratic term.

Canonical signature p is 1.879/1.902 (balanced), 1.837/1.807 (random), and
1.570/1.849 (pectinate) over 5k-10k/10k-20k. The signature cost is therefore
also superlinear here, but it is not responsible for the pectinate-specific
gap. Selected n=20,000 direct sub-helper medians are:

| Phase | Balanced, s | Random, s | Pectinate, s |
|---|---:|---:|---:|
| Connectivity traversal | 10.08 | 9.53 | 12.99 |
| Root-distance traversal | 10.56 | 10.11 | 13.59 |
| Descendant-key construction | 8.46 | 7.78 | 44.42 |
| Signature descendant traversal | 3.75 | 3.50 | 3.60 |

Root identification, edge normalization, tip mapping, postorder reordering,
and one parent/child construction were each 0-0.02 s at this scale. They are
not current wall-time hotspots when measured independently.

## Operation Evidence

The public preparation call count is stable at every measured size and shape:
one `.inspect_tree_core()`, one `.safe_canonicalize_core()`, two
`.canonical_tree_signature()` calls, three `.edge_children_index()` calls, one
`.descendant_keys_iterative()`, and five complete parent/child reconstructions
when the two inspection `split()` constructions are included.

For a valid binary n=20,000 fixture (`E = 39,998`):

- per-node full-edge searches in inspection and canonicalization: **0**;
- static full-edge/vector pass inventory: inspection 8, one signature 5,
  canonicalization body 7;
- structural inspection traversals: 2, connectivity and root distance;
- root descendant traversals: 2 inside canonicalization, one per signature;
- canonical descendant-key passes: 1;
- named adjacency lookups across the measured preparation operations: 159,992;
- descendant-key sort calls and string collapses: 19,998 each.

The decisive shape-dependent operation count is descendant-tip payload visits:

| Shape | Signature tip visits | Canonical total | Descendant-key payload only |
|---|---:|---:|---:|
| Balanced | 20,000 | 307,232 | 267,232 |
| Random | 20,000 | 396,178 | 356,178 |
| Pectinate | 20,000 | 200,029,999 | 199,989,999 |

Instrumentation wrappers were used only for call counts, never for the direct
timing tables. At n=20,000 their median calibration overhead was between
-0.281% and +0.284%, with maximum absolute overhead 1.351%. Smaller fast cells
had timer-resolution noise above 5%; their instrumented wall times are not used
as evidence, while their deterministic operation counts remain valid.

## Source-Level Cause

`.inspect_tree_core()` constructs child groups twice and both traversals use
`children[[as.character(node)]]`, `stack[-length(stack)]`, and
`c(stack, kid)` (`R/check_tree.R:307-325` and `R/check_tree.R:405-420`). There is
no `which()`, `match()`, `%in%`, or `edge[edge[, 1] == node, ]` inside these
loops. The stack operations copy active R vectors; named-list lookup may scan a
growing name table. Because balanced and random inspection also have p near 2,
the measured evidence implicates name lookup/copy cost rather than tree depth
alone. Pectinate depth raises the constant through additional stack copying.

`.descendant_keys_iterative()` memoizes graph visits, but each internal node
still combines, sorts, stores, and finally collapses its complete descendant
label vector (`R/analysis_preparation.R:192-264`). Its cumulative payload is
approximately `sum_v s(v)`: near n log n for balanced/random trees and near
n-squared for a pectinate tree. At n=20,000 this one direct helper took 44.42 s,
55.14% of pectinate `prepare_tree()` wall time.

The fixed number of full-edge passes and five adjacency reconstructions increase
constants but cannot explain p near 2 by themselves. The compiled tree path in
`src/tree_core.cpp` already uses indexed arrays and a reserved explicit stack;
no corresponding superlinear scan was found there.

## Evidence Reuse Safety

| Class | Evidence | Reuse decision |
|---|---|---|
| A: invariant | tip labels/counts, topology-valid booleans, root degree, polytomy/unary counts, branch diagnostics | Safe as scalar or label-keyed facts after the existing proof |
| B: remappable | parent/child arrays, adjacency, root ID, root distances, edge-associated lengths, preorder/postorder, method readiness | Safe only with exact node-ID and edge-row remaps |
| C: representation proof/workspace | DFS stacks, before/after signatures, canonical map proof, edge-isomorphism records, canonicalization metadata, fingerprint | Not safe to reuse across canonicalization |

An index is safe when built and consumed within one unchanged representation.
Blind pre-to-post canonicalization reuse is unsafe because canonicalization may
change edge order, internal IDs, root representation, and pruning order. The
before/after signatures and final fingerprint remain required.

## Candidate Ranking

### Candidate 1: integer adjacency plus cursor stacks for inspection

- **Work:** replace character-keyed child lookup and copying R stack mutation in
  the two inspection traversals, without removing any validation.
- **Expected complexity:** the traversal/index component moves from observed
  O(n-squared)-like behavior toward O(n + E); fixed validation passes remain.
- **Measured share:** inspection is 54.85%, 54.17%, and 33.54% of n=20,000
  preparation for balanced, random, and pectinate trees. The narrower sum of the
  two traversal bodies gives theoretical ceilings of 54.30%, 52.61%, and
  32.99%. Actual savings will be lower.
- **Risk:** low to medium; preserve visit order, cycle handling, overflow,
  issue order, and early-failure semantics.
- **Equality fixtures:** balanced/random/pectinate, polytomy, shuffled edges,
  safe renumbering, cycles/disconnection, unary/malformed trees, invalid root,
  zero/negative branches, and exact complete inspection-object equality.
- **Gate:** **GO** for a bounded future candidate.

### Candidate 2: exact compact descendant-key ordering

- **Work:** eliminate repeated materialization/sort/collapse of every ancestor's
  full descendant label payload in `.descendant_keys_iterative()`.
- **Expected complexity:** target an exact collision-free persistent key or
  comparator whose construction is O(n log n)-like rather than proportional to
  the pectinate `sum_v s(v)`. This is a private computation change only;
  `.canonical_tree_signature()` semantics and output must remain byte-identical.
- **Measured share:** descendant-key timing is 22.26%, 20.84%, and 55.14% of
  n=20,000 preparation for balanced, random, and pectinate trees. These are the
  ideal end-to-end ceilings; actual savings will be lower.
- **Risk:** medium to high because exact lexical ordering, arbitrary tip labels,
  polytomies, cycles/fallback behavior, and stable tie-breaking must not drift.
- **Equality fixtures:** exact old/new internal mapping and edge order for all
  shape fixtures, shuffled edges, safe renumberings/root representations,
  polytomies, adversarial labels, before/after signature identity, branch-edge
  association, and canonicalization failure reasons.
- **Gate:** **GO** for a bounded prototype with old-vs-new exact-equality gates.

Cross-canonicalization reuse of a complete structural evidence package remains
`NO-GO / NEED_MORE_EVIDENCE`; its remapping risk is higher and the reusable
fraction has not been isolated from the two confirmed hotspots.

## Final Gate

`STAGE2B2A_AUDIT = PASS`

`CANDIDATE_1 = GO`

`CANDIDATE_2 = GO_FOR_BOUNDED_EXACT_EQUIVALENCE_PROTOTYPE`

`CROSS_REPRESENTATION_EVIDENCE_REUSE = NO_GO_NEED_MORE_EVIDENCE`

No estimator, numerical formula, tolerance, RNG behavior, thread policy,
fingerprint semantics, canonical-signature semantics, public API, production
R/C++ file, or test was changed. Stage 2B2A stops here.
