# fastphylosig 0.2.0 Stage 2B2A: Evidence Reuse Audit

## Scope and status

This is a read-only source audit for the preparation engine. No production
code, tests, public API, statistical definition, tolerance, RNG policy, or
thread policy was changed. No formal Stage 2B2A benchmark was run by this
subtask. The Stage 2B1 boundary remains frozen:

`PREPARATION_ENGINE_STAGE1 = ACCEPTED_FROZEN`

`RAW_ROUTE_SINGLE_PREPARATION_BOUNDARY = ACCEPTED_FROZEN`

The evidence matrix is in `evidence_reuse_matrix.csv`; the static operation
inventory is in `operation_count_audit.csv`.

## Audit basis

The audit follows the current source at the Stage 2B1 accepted tree. The raw
preparation route canonicalizes first, inspects the canonicalized tree once,
then calls the private preparation core with package-owned evidence. The
public `prepare_tree()` route still performs its own complete inspection and
the existing canonicalization probe. The relevant implementation locations
are:

- `R/check_tree.R:50-423` for structural inspection, connectivity, and root distances.
- `R/analysis_preparation.R:127-419` for adjacency, iterative descendant traversal, canonical signatures, and safe canonicalization.
- `R/prepare_tree.R:17-132,155-200` for preparation and compiled structural entries.

## A/B/C classification

### A. Representation invariant

These values describe the biological tree rather than its storage order or
internal node IDs, provided the existing safety proof has established unchanged
tip labels and branch associations:

- tip labels and `n_tip`;
- `Nnode` and node/edge counts;
- edge validity predicates;
- unique-root and connectivity predicates as booleans;
- root degree, polytomy, and single-child counts;
- finite/nonnegative/zero/near-zero branch diagnostics and scalar ranges;
- descendant tip-label sets, considered as label-defined sets only.

An A item may be reused as a scalar or label-keyed fact. An ID-indexed array
derived from it is not automatically A; it falls under B.

### B. Representation dependent but safely remappable

The following facts remain scientifically meaningful after canonicalization but
their storage coordinates can change:

- `parent[]`, `child[]`, and `edge_for_node[]`;
- children adjacency keyed by parent ID;
- indegree/outdegree arrays;
- the numeric root ID;
- root-distance vectors indexed by node ID;
- edge-associated branch-length vectors;
- preorder, postorder, and pruningwise edge orders;
- method readiness when it includes D/Delta's canonical root-ID requirement;
- the inspection object when it is carried into a different representation.

These can only be reused after applying the exact internal-node and edge-row
mapping. Direct reuse is unsafe. In particular, inspection creates two
different adjacency payloads: child IDs for connectivity and edge-row indices
for root distances (`R/check_tree.R:304-307,405-420`). A single generic
`children` object cannot replace both.

### C. Unsafe across canonicalization

The following are proof or workspace state tied to one exact representation:

- mutable `seen`/stack traversal workspaces;
- canonicalization's `before` and `after` signatures;
- the canonical internal-node mapping itself;
- edge-isomorphism records used to prove the transformation;
- the canonicalization info attribute;
- exact tree fingerprint strings.

These are not interchangeable with inspection evidence and must not be used to
bypass a second representation-specific proof. The exact fingerprint contract
also means a fingerprint must be recomputed after canonicalization.

## Shared work found in source

### Inspection

`.inspect_tree_core()` performs two structural adjacency constructions and two
iterative traversals:

1. `split(child, parent)` followed by a connectivity DFS;
2. `split(seq_along(child), parent)` followed by a root-distance DFS.

Both traversals pop with `stack <- stack[-length(stack)]` and append with
`stack <- c(stack, kid)`. The visit count is linear (`one node/edge visit`),
but the R vector mutation copies the current stack. A deep pectinate tree can
keep O(n) pending nodes, making the aggregate copy work O(n^2)-like. The same
source pattern occurs at `R/check_tree.R:313,325` and `409,420`.

The source contains no per-node `which()`, `match()`, `%in%`, or full edge-matrix
subset inside these DFS loops. Therefore the confirmed source-level concern is
repeated vector copying, not an `edge[edge[,1] == node,]` search.

### Canonicalization

`.safe_canonicalize_core()` invokes `.canonical_tree_signature()` before and
after the mapping. Each signature builds an adjacency index and traverses all
descendants of the root. The canonicalization body builds another adjacency
index and calls `.descendant_keys_iterative()` for every non-root internal.

The important superlinear pattern is in
`.descendant_keys_iterative()` (`R/analysis_preparation.R:192-264`): each
internal node combines and sorts the full descendant-key value for its
children. The total work is proportional to

`sum over internal nodes of descendant-set size`.

For a balanced binary tree this is approximately `O(n log n)` before string
constant factors. For a pectinate tree it is `O(n^2)`-like because almost every
internal node carries a descendant set of order n. Random trees should fall
between these extremes and must be measured in the Stage 2B2A shape profile.

Other canonicalization costs are linear or `O(E log E)` source-level work:
mapping, node-label remapping, branch-length sorting, and two edge-record
sorts. The two signature calls are semantically required by the current
canonicalization safety proof; this audit does not authorize removing either.

## Adjacency/index reuse design audit

A package-owned index could be represented as:

`parent[]`, `child_offsets[]`, `children[]`, `edge_for_node[]`,
`preorder[]`, and `postorder[]`.

On a fixed representation, this index would support O(1)-average child lookup
and array-backed iterative traversals. It could remove the repeated `split()`
calls and, if stack storage is preallocated, avoid the vector-copy behavior.

The canonicalization boundary is the key constraint:

1. canonicalization may change edge-row order;
2. canonicalization may change internal node numbering;
3. canonicalization may call `ape::reorder.phylo(..., "postorder")`;
4. D/Delta readiness observes the canonical root ID.

Therefore an index built before canonicalization is B, not A. It can be
remapped only with an exact node map plus edge-row permutation, and the
edge-associated branch lengths must follow the same edge map. If the
transformation cannot provide both maps, the index is unsafe to reuse. An
index built after canonicalization is valid for the post-canonical inspection
and preparation core, but it cannot be used as evidence about the
pre-canonical representation's safety proof.

The safe design boundary is consequently either:

- build and consume one index entirely within one unchanged representation;
  or
- have canonicalization emit explicit remappable evidence, then rebuild the
  post-canonical method-readiness view from the remapped arrays.

No design here permits skipping mutation detection or the required public
prepared-context fingerprint.

## Candidate ranking for a future production stage

Only two candidates are justified by this source audit. They are proposals,
not implementations or authorization to start Stage 2B2B.

### Candidate 1: array-backed structural traversal

**Precise duplicated/superlinear work.** The inspection connectivity and
root-distance DFS loops each mutate growing/shrinking R vectors. On a
pectinate tree, aggregate stack copying can be O(n^2)-like despite linear edge
visit counts.

**Expected complexity change.** Keep the same traversal semantics but use a
preallocated integer stack or an index cursor; use one adjacency index per
fixed representation. The copy component becomes O(n), leaving the required
edge/node visits linear.

**Affected wall-time share.** Stage 2A direct microprobes at n=20,000 reported
`structural_inspection` about 20.49 s versus a preserved public
`prepare_tree()` total about 38.515 s. This is an inclusive, non-additive
ranking signal, not a component sum. The maximum end-to-end ceiling is
therefore at most about 53.2% of that public preparation total, and only the
copy-related portion is actually recoverable.

**Risk.** Medium: traversal order and cycle/connectedness failure semantics
must remain identical. The stack implementation must preserve the current
edge order and early-failure behavior.

**Required equality fixtures.** Balanced, random, pectinate, polytomy,
shuffled edge order, safe internal-node renumbering, malformed disconnected
graphs, cycles, unary nodes, invalid root numbering, and large pectinate
stress at n=1,000/5,000/20,000. Compare complete inspection objects and issue
tables, not only readiness booleans.

**Authorization.** `GO` for a bounded future candidate only after shape-specific
operation-count and wall-time evidence confirms stack-copy dominance. No source
implementation is included here.

### Candidate 2: one remappable structural-evidence package across the raw boundary

**Precise duplicated work.** Inspection and canonicalization independently
construct parent/child adjacency, root/degree metadata, and traversal support.
The two representations are not identical, so only the explicitly invariant
or remapped portions can be shared.

**Expected complexity change.** Build one package-owned index on the
pre-canonical representation, carry an exact node/edge remap through safe
canonicalization, and consume the remapped index for the post-canonical
inspection. This can reduce repeated `split()`/degree construction and avoid
repeating a full structural traversal where the transformed evidence is
provably equivalent. It must not remove the canonical before/after signature
proof or final fingerprint.

**Affected wall-time share.** Stage 2A direct canonicalization was about
16.995 s at n=20,000 (about 44.1% of the preserved public total), while
structural inspection was about 20.49 s. These figures overlap; 44.1% is only
an inclusive upper bound, not an expected saving. The actual reusable subset
needs operation-count evidence first.

**Risk.** High: edge order, internal numbering, postorder order, root-ID
requirements, branch-length association, and issue ordering all cross the
representation boundary. Incorrect remapping could change D/Delta readiness
or failure semantics.

**Required equality fixtures.** In addition to Candidate 1 fixtures, compare
pre/post canonical mappings, exact fingerprints, canonical signatures,
complete inspection issue order, root distances after remap, branch-length
edge association, raw/prepared parity for K/lambda/D/Delta, and mutation
detection across every protected field.

**Authorization.** `NEED_MORE_EVIDENCE` until a shape-specific profile proves
that the remappable subset is a material wall-time share and an equality audit
shows no representation or diagnostic drift.

## Conclusions for Stage 2B2A

- The source contains no per-node full-edge `which`/`match` search in the
  inspected DFS loops.
- Inspection has two required linear traversals, but both use
  shape-sensitive R vector stack copying; pectinate trees can make this
  O(n^2)-like.
- Canonicalization has a clearer shape-sensitive repeated descendant-key
  operation: `sum |desc(v)|`, O(n log n)-like for balanced trees and O(n^2)-like
  for pectinate trees.
- A shared adjacency/index is feasible only within one representation or with
  an exact node/edge remap. It cannot be blindly shared across canonicalization.
- Candidate 1 is the highest-priority bounded production investigation.
- Candidate 2 is higher risk and remains `NEED_MORE_EVIDENCE`.

Overall gate:

`STAGE2B2A_EVIDENCE_REUSE_AUDIT = COMPLETE`

`PRODUCTION_CODE_CHANGE = NONE`

`STAGE2B2B_AUTHORIZATION = NO-GO_PENDING_SHAPE_PROFILE`

The audit ends here. No Stage 2B2B implementation, benchmark, or production
refactor was started.
