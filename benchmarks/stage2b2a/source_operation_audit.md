# fastphylosig 0.2.0 Stage 2B2A Source Operation Audit

## Scope

This is a read-only source audit. It does not modify production code, tests,
the public API, statistical definitions, numerical tolerances, RNG behavior, or
thread policy. The structured findings are in
`source_operation_audit.csv`.

The audit covers the preparation path in:

- `R/check_tree.R`
- `R/analysis_preparation.R`
- `R/prepare_tree.R`
- `R/match_phylo_data.R` for downstream comparison
- `src/tree_core.cpp` for the compiled-tree path

Let `n` be the number of tips, `N` the total node count, and `E` the edge
count. For a rooted binary tree, `N = 2n - 1` and `E = 2n - 2`. The notation
`S_stack(T)` below means the sum of the R stack-vector lengths copied by all
pop/append mutations in a traversal. The notation `s(v)` means the number of
descendant tips below internal node `v`.

## Findings

### 1. Inspection has an explicit shape-dependent vector-copy hazard

`.inspect_tree_core()` does not perform `edge[edge[,1] == node, ]`, `which()` or
`match()` inside its graph loops. It builds a `children` index with `split()`
and then visits every edge once for connectivity and once for root distances.
However, both traversals mutate an ordinary R vector on every operation:

```r
stack <- stack[-length(stack)]
stack <- c(stack, kid)
```

The deletion and append allocate/copy the active stack repeatedly. The graph
work is therefore linear in node/edge visits, but the copy work is
`O(S_stack(T))`. A deep or pectinate tree can keep a long list of pending
sibling nodes, giving a quadratic-like `S_stack(T)`. A balanced tree keeps the
frontier much smaller, so the same code is closer to an `n log n`-like copy
profile. Random trees depend on their depth and DFS frontier. Root-distance
traversal repeats this same pattern independently.

This is the strongest source-level explanation for the large inspection
component observed in the Stage 2A data. It is repeated vector copying, not a
missed adjacency lookup.

There is a separate lookup cost that must not be conflated with a full edge
scan. Both inspection traversals and the canonical helpers retrieve children
with `children[[as.character(node)]]` (or `children[[key]]`) from a named list
created by `split()`. The number of such lookups is linear: approximately one
per visited internal node in each traversal, and approximately two per
internal node in the enter/exit descendant-key walk. The package does not own
an integer-indexed child array in these R helpers. Consequently, the constant
time assumption for each name lookup depends on the active R implementation's
name-resolution behavior. If the names table is hashed, the expected cost is
O(1); if a lookup scans the names vector, the per-lookup upper bound is the
number of parent groups and the aggregate can become another O(n^2)-like term.
This source audit therefore records named-list lookup counts separately and
marks their asymptotic cost `NEEDS_MEASUREMENT`; no claim that they are
quadratic is made without the required operation/timing evidence.

### 2. Canonicalization has descendant-payload superlinear work

The iterative traversal itself is depth-safe and memoized. In
`.descendant_keys_iterative()`, each completed internal node is assembled by
combining its child keys, sorting the resulting descendant vector, and later
turning that vector into a string:

```r
values <- unlist(lapply(kids, function(child) get(...)))
assign(key, sort(values), envir = memo)
...
paste(get(as.character(node), memo, inherits = FALSE), collapse = "\\r")
```

Memoization prevents re-traversing the same graph, but it does not prevent
copying a descendant payload once for every ancestor. The aggregate cost is

`sum_v O(s(v) log s(v))` for the per-node sorts, plus character-copy work
proportional to `sum_v character_count(s(v))` for the strings. A pectinate tree
maximizes these sums and can approach `O(n^2 log n)` comparison/copy work; a
balanced tree is substantially smaller (roughly `n log^2 n` for the sort
term). This is a canonicalization-side superlinear candidate distinct from
the inspection stack issue.

`.canonical_tree_signature()` itself has ordinary linear scans and a final
tip-label sort. `.safe_canonicalize_core()` also performs two signature calls,
one before and one after normalization, and sorts edge records for exact
isomorphism checks. These are repeated `O(E)`/`O(E log E)` passes, but they are
not themselves quadratic. The descendant-key aggregation is the only direct
R source pattern in this path that repeatedly materializes all ancestors'
descendant payloads.

### 3. No per-node full-edge search was found

The source search found no occurrence in the audited graph loops of:

- `which()` or `match()` over the edge matrix per node;
- `%in%` over edge rows per node;
- `edge[edge[,1] == node, ]` or equivalent full edge-matrix subsets;
- explicit ancestor search or descendant search by rescanning the edge list.

The production R helpers use `split()`-based child grouping and direct indexed
vectors. The C++ compiler in `src/tree_core.cpp` builds parent/child arrays and
uses a reserved explicit stack; its parse, validation, and traversal loops are
`O(N+E)` with no repeated full-edge scan. Therefore rewriting the compiler in
C++ would not address the observed source-level risk.

The absence of full-edge rescans does not imply that all indexing is free:
`children[[key]]` is a named-list lookup at each internal visit. Its count is
linear, but its implementation-dependent name-resolution cost must be measured
before a candidate is ranked as a confirmed complexity fix.

### 4. Other repeated work is linear or downstream

`.inspect_tree_core()` performs several necessary full-vector checks for labels,
edge ranges, degrees, topology, and branch geometry. The public preparation
core also performs required validation, fingerprint construction, and one
compiled-tree build. These are constant numbers of linear passes.

Trait matching uses `%in%`, `setdiff()`, and `match()` a constant number of
times, while the fallback NA grouping uses `which()` over trait columns. Those
operations scale with data rows/columns and are not the cause of full
`prepare_tree()` scaling. Retained-subtree cache misses similarly operate after
matching and are not part of the full-tree preparation root cause.

## Static operation-count summary

| Area | Full edge scans / indexed passes | Nested or copied work | Shape sensitivity |
|---|---:|---|---|
| Inspection vector validation | At least 10 linear edge/branch passes, plus node scans | None beyond vector allocation | Low |
| Inspection connectivity | One indexed child iteration per edge | `S_stack(T)` from R stack pop/append | High; deepest for pectinate |
| Inspection root distances | One indexed row iteration per edge | A second `S_stack(T)` term | High; deepest for pectinate |
| Named-list child lookup | About `N` per traversal | O(1) expected if hashed; up to O(P) if linear name matching | Needs measurement |
| Canonical signature | Several linear scans; `O(E log E)` ID sort | One root traversal and tip-label sort | Low to moderate |
| Canonical descendant keys | No repeated edge scan | `sum_v s(v) log s(v)` plus repeated strings; about `2M` named-list child lookups | High; pectinate worst case |
| Canonical edge-isomorphism proof | Two full record builds and sorts | `2 O(E log E)` | Low |
| C++ compile path | Indexed edge/node loops | None asymptotically superlinear | Low |

These are static operation-site counts, not instrumented wall-time counts. A
formal operation-counter run should count stack-vector copied elements and
descendant payload lengths separately; counting only function entries would
miss both superlinear mechanisms.

## Evidence and candidate implications

The source supports two bounded production candidates, subject to the required
tree-shape timings, operation counts, and parity fixtures:

1. **Replace the two inspection R stack mutations with an indexed/cursor stack.**
   This targets the confirmed copied-volume term in connectivity and
   root-distance traversal and changes that component toward `O(N+E)`. The
   theoretical end-to-end ceiling is the measured inspection share, not the
   whole preparation time, because other checks remain required. Risk is low
   if traversal order, cycle handling, root-distance overflow behavior, and
   all failure statuses are held exactly constant. Required equality fixtures:
   balanced, random, pectinate, polytomy, unary/malformed trees, zero and
   negative branches, root renumberings, and mutation/failure contracts.

2. **Share or safely remap structural traversal evidence between inspection and
   canonicalization.** This would target the duplicate parent/child and
   traversal evidence, including the descendant-key work where the
   representation is unchanged. The theoretical end-to-end ceiling is bounded
   by the overlapping inspection/canonicalization shares; it cannot be added
   naively. Risk is materially higher because canonicalization can change edge
   order, internal IDs, and root representation. Evidence must be classified
   as representation-invariant, safely remappable, or unsafe before reuse, and
   canonical signature semantics must remain byte-for-byte identical.

No fingerprint rewrite or signature representation rewrite is justified by
this source audit. No production candidate is authorized from source evidence
alone: tree-shape scaling and low-overhead operation counters are still needed
for `GO`; otherwise the decision remains `NO-GO / NEED_MORE_EVIDENCE`.

## Audit status

- `PRODUCTION_CODE_CHANGE`: NONE
- `TEST_CHANGE`: NONE
- `PUBLIC_API_CHANGE`: NONE
- `SOURCE_OPERATION_AUDIT`: COMPLETE
- `CANDIDATE_AUTHORIZATION`: `NO-GO / NEED_MORE_EVIDENCE` pending empirical shape and counter evidence
