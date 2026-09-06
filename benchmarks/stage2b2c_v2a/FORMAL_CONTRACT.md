# fastphylosig 0.2.0 Canonicalization Contract V2-A

## Status and scope

This document is the normative design for a future V2 canonicalization
contract. It is an audit and design artifact only.

It does not replace the frozen V1 production implementation. It does not
modify production source, tests, public documentation, estimators, RNG
behavior, thread policy, or public APIs. No V2 prototype is integrated by
this document.

Current baseline:

    PREPARATION_ENGINE_STAGE1 = ACCEPTED_FROZEN
    RAW_ROUTE_SINGLE_PREPARATION_BOUNDARY = ACCEPTED_FROZEN
    INSPECTION_INTEGER_ADJACENCY = ACCEPTED_FROZEN
    CANDIDATE_2 = REJECTED
    CANONICALIZATION_CONTRACT_BUG = CONFIRMED
    PUBLIC_RELEASE_BLOCKER_0_2_0 = TRUE
    STAGE_2C = NOT_AUTHORIZED

The contract is based on the V1 source and the Stage 2B2C audit recorded in
commits dd8ec25 and 654eb8e. The audit demonstrated:

- descendant keys can collide for legal unique tip labels because V1 uses
  delimiter concatenation with "\r";
- V1 canonical ordering uses the active LC_COLLATE;
- V1 uses original internal node IDs as a tie break;
- V1 fingerprints can accept a protected tip-label mutation and permit stale
  prepared-cache reuse.

The V2 design separates canonical identity from exact prepared-context
mutation validation. A canonical identity may be representation invariant
while an already prepared context still rejects a change in its source
representation.

## Contract vocabulary

Let T be a finite rooted directed tree with a unique biological root,
unique tip labels, and one biological branch record per parent-child edge.
Input node IDs and edge-row positions are source representation details, not
biological identity.

B(s) denotes the exact byte sequence obtained from the validated UTF-8
encoding of a tip label s. Byte length means the number of bytes, not the
number of characters. A V2 canonical representation is a new copy; the
source phylo object is never modified.

The contract has two distinct byte encodings:

1. C2(T), the canonical identity encoding, is invariant to the permitted
   representation changes.
2. S2(T), the source protected snapshot, records the exact protected fields
   of one prepared context before canonicalization and is used for mutation
   rejection.

C2(T) and S2(T) must never be conflated. A digest of either encoding is only
an acceleration hint unless the complete exact bytes are also checked.

## A. Tip label identity

1. A tip label is identified by the exact validated UTF-8 bytes B(s).
2. No Unicode normalization, case folding, trimming, transliteration, or
   locale collation is applied. Thus precomposed and decomposed Unicode
   sequences are distinct labels unless their byte sequences are identical.
3. NA, missing labels, conversion failures, invalid UTF-8, and empty labels
   are failures before descriptor construction. Non-empty whitespace is label
   content and is not trimmed or silently rejected by this contract.
4. Duplicate detection is byte-exact after UTF-8 conversion. It is not based
   on sort(character), order(character), case-insensitive comparison, or the
   active locale.
5. Separator-containing labels, including labels containing carriage return
   or line feed, are legal when they satisfy the label contract. Their bytes
   are length-delimited in every V2 encoding.
6. The original tree$tip.label values are not rewritten. The comparison key
   is separate from the stored label.

## B. Tip label ordering

The canonical tip rank is a 1-based integer assigned by the following
locale-independent comparator:

1. Compare the UTF-8 byte sequences as unsigned bytes from left to right.
2. At the first differing byte, the smaller unsigned byte sorts first.
3. If one sequence is a strict prefix of the other, the shorter sequence sorts
   first.
4. Equal byte sequences are a duplicate-label failure, never a tie to resolve.

The implementation must use an explicit byte comparator or an equivalent
specified bytewise primitive. It must not use the semantics of R
sort(character) or order(character) as canonical truth. It must not change
the process locale to obtain the result.

Canonical tip rank is an internal comparison key, not a renumbering of public
tip IDs. V2 retains source tip IDs 1, ..., n_tip and the original tip.label
vector exactly, including order and stored string values. This preserves the
index semantics of controlled permutation/null inputs while making internal
ordering independent of locale.

## C. Exact rooted-subtree identity

### C.1 Descriptor definition

For every tip with canonical tip rank r, define the descriptor:

    Tip(r)

For every internal node with children v1, ..., vk, order the children with
the V2 rule in Section D and define:

    Node(k, d1, ..., dk)

where d1, ..., dk are the child descriptors in Section D order.

Descriptor equality is structural equality of the kind tag, arity, and the
ordered child descriptors. It is not equality of a hash, printable string,
source node ID, or edge-row position.

### C.2 Compact storage

The prototype must build an immutable descriptor table. A descriptor record
contains only:

- a kind tag (Tip or Node);
- a tip rank for a tip, or an arity for an internal node;
- integer references to the ordered child descriptors for an internal node.

The descriptor references are handles, not canonical truth. Handle allocation
order must not affect comparison, output, or fingerprint values. Exact
descriptor comparison follows descriptor content recursively through the
table. It must not expand a descriptor into its complete descendant-tip label
sequence.

Descriptor IDs are separate from output node IDs. Tip descriptor IDs are tip
ranks; internal descriptor IDs follow the canonical internal traversal from
Section E. They are therefore derived from label bytes and rooted structure,
never from source handles. Exact descriptor equality recursively follows the
descriptor records, not the numeric handle alone. An optional intern table may
use a hash to find candidates, but exact typed tuple equality must decide
identity.

### C.3 Exactness obligation

For a valid tree with unique non-empty tip labels, descriptor identity is
injective over distinct rooted subtrees:

1. A tip descriptor identifies one unique tip byte sequence through its rank.
2. By induction on subtree height, equal internal descriptors have the same
   arity and equal ordered child descriptors.
3. Therefore equal descriptors have the same rooted labeled subtree.
4. Two distinct subtrees of one valid tree have disjoint non-empty descendant
   tip sets, so they cannot have equal descriptors when tip labels are unique.

If a descriptor tie is nevertheless observed, the prototype must report an
ambiguous canonicalization failure. It must not fall back to the original
internal node ID or input edge-row order.

## D. Child ordering

For each node, define min_tip_rank as its tip rank when it is a tip and as the
minimum min_tip_rank among its children when it is internal. Children of one
internal node are sorted by increasing min_tip_rank.

This is an exact total order for valid unique-tip-label trees: child subtrees
are disjoint and non-empty, so their minimum tip ranks cannot tie. A tie is an
ambiguous-canonicalization failure. The rule is independent of branch length,
source node ID, source edge-row position, descriptor handle allocation, and
locale. It needs one integer per node and never materializes descendant-label
vectors or pairwise descriptor comparisons.

## E. Internal node numbering

Canonical numbering is assigned after child order is fixed:

1. Tips retain IDs 1, ..., n_tip and their original tip.label order.
2. The root receives internal ID n_tip + 1.
3. Starting at the root, visit children in canonical child order.
4. Assign the next internal ID when an internal node is first entered
   (deterministic pre-order).
5. Continue until every internal node has exactly one canonical ID.

The resulting internal IDs are a consequence of the rooted canonical
representation. Original internal IDs are used only as temporary lookup keys
while reading the source edge list. They are not a tie break and never enter
C2(T).

Nnode in the canonical output is the validated number of internal nodes.
The canonical root is always n_tip + 1. A non-contiguous or inconsistent
source numbering is remapped when the graph is otherwise valid; a malformed
graph is rejected.

Node labels are not computational identity. If a future prototype carries
node.label, it may remap them through the internal-node map, but they are not
part of canonical ordering or the V2 fingerprint unless a later schema
explicitly promotes them to protected state.

## F. Edge-row ordering

Canonical edge rows are emitted by an ordered postorder traversal of the
canonical rooted representation:

1. Traverse children in canonical child order.
2. On leaving each non-root node, emit its parent-to-child edge record.
3. Emit no record based on the original row index.

The exact traversal and emission rule is part of the V2 schema. The output
edge matrix therefore depends only on canonical root, canonical node IDs, and
canonical child order. Shuffling source edge rows cannot change it.

The prototype must construct adjacency from the source graph without using
source row order as a tie break. If two source rows claim the same biological
parent-child edge, the graph is malformed and canonicalization fails.

## G. Root placement and rootedness

The source must contain exactly one structural root:

- one internal node has indegree zero;
- every other node has the required parent relationship;
- all nodes are reachable from that root;
- there is no self-loop or directed cycle;
- the root is not a tip.

The root's source ID is ignored after validation. The canonical root is fixed
to n_tip + 1, and root placement is independent of source numbering and
edge-row order.

An unrooted tree or an unsupported root representation is rejected. V2 does
not biologically reroot a tree, select an outgroup, invent a root, or change
the direction of an edge. Rooting remains an explicit user operation outside
this contract.

## H. Polytomy handling

An internal descriptor stores its actual child arity k. Rooted polytomies are
retained without resolution, binary expansion, or arbitrary child insertion.
Their children are ordered by Section D only.

Unary internal nodes are not collapsed or expanded. If current inspection
marks a unary node as method-incompatible, V2 preserves that readiness and
failure behavior. A zero-child internal node, disconnected component, cycle,
or other malformed graph is rejected rather than converted into a legal tree.

Canonicalization never changes whether a method supports a polytomy or unary
node. Method-specific readiness remains the responsibility of the existing
inspection contract.

## I. Branch-length association

Each source edge is treated as one record:

    (source parent node, source child node, exact branch length)

The branch length is remapped together with the biological parent-child edge
when canonical internal-node IDs are assigned. Tip IDs are retained.

The V2 implementation must prove a one-to-one association between source edge
records and canonical edge records. It must not sort edge.length separately,
recompute a length from a path, or match lengths by row position after an
edge reorder.

Heterogeneous branch lengths are required fixtures. Zero lengths retain their
existing semantics. Negative, missing, non-finite, or length-mismatched branch
values fail according to the existing tree contract. No branch jitter,
rounding, or tolerance is introduced.

Branch lengths do not participate in subtree label identity or child ordering.
They are carried by the already identified biological edge and are included
in the canonical fingerprint.

## J. Fingerprint V2 encoding

### J.1 Canonical identity encoding C2

C2(T) is an exact, printable-or-raw serialization of the canonical
representation. A printable form such as hex or base64 is only a transport
encoding of the exact bytes. The binary form is authoritative.

The encoding uses a fixed domain tag and schema version followed by typed
length-delimited fields in this fixed order:

1. domain tag: fastphylosig.canonical.tree;
2. canonical schema version: integer 2;
3. canonical tip count and the canonical tip-label sequence;
4. canonical internal-node count (Nnode);
5. canonical edge count and the canonical edge endpoint matrix;
6. canonical branch-length count and the row-aligned branch-length sequence.

Each field has an explicit field tag and byte length. Each vector has an
explicit element count. Each string has an explicit byte length followed by
the exact UTF-8 bytes. Each node endpoint, count, and Nnode uses a fixed
unsigned big-endian integer representation. Each branch length uses its exact
IEEE-754 binary64 representation in fixed big-endian byte order after finite
non-negative validation. No delimiter is structural syntax.

The edge and branch-length fields are emitted in the same canonical row
order. Thus C2(T) distinguishes a different branch association even when the
multiset of branch lengths is unchanged.

### J.2 Protected source snapshot S2

S2(T) is a separate exact encoding used only for prepared-context integrity.
It includes a source-representation domain tag, schema version, and the
following protected fields in source order:

1. tip-label vector length and each exact UTF-8 label byte sequence;
2. edge matrix dimensions and every endpoint in source row order;
3. branch-length vector length and every exact binary64 value in source row
   order;
4. Nnode;
5. the validated field-shape metadata needed to distinguish vector lengths and
   matrix dimensions.

The snapshot encoding is length-delimited and type/domain separated. A
delimiter-containing label cannot collide with a different label vector.
node.label and other non-computational metadata are outside this snapshot
under the current documented contract. Promoting such a field requires a new
schema and a new snapshot definition.

S2(T) is compared byte-for-byte. A digest may be stored as a lookup hint, but
a digest match alone never authorizes context reuse.

### J.3 Fingerprint claims

The V2 contract may claim that equal C2(T) means equal canonical protected
representation only after exact serialization equality. It may not claim
collision resistance from a digest alone. A cryptographic digest is optional
and non-authoritative.

The canonical identity is independent of active LC_COLLATE, original internal
node IDs, and input edge-row order. The source snapshot deliberately retains
source order so that an already prepared context can reject any protected
field mutation or representation change before reuse.

## K. Exact protected-snapshot mutation validation

### K.1 External boundary

Every external prepared-tree public call must perform this sequence:

1. Verify the private V2 schema marker and package-owned context token.
2. Read the current protected fields from the context's tree.
3. Validate their basic shape and value contract.
4. Encode the current fields as S2_current.
5. Compare S2_current byte-for-byte with the package-owned stored snapshot.
6. Reject before any structural-cache lookup, numerical-cache lookup, or
   estimator call if the comparison fails.
7. Create a private validated token for downstream helpers after the exact
   comparison succeeds.

The downstream private route may trust only that token. The token is not a
public argument and there is no trust=, unsafe=, skip_check=, engine=, or
backend= escape hatch.

### K.2 Snapshot ownership

The stored snapshot must be a deep-copied raw byte object held in a
package-owned private binding or equivalent package-owned immutable record.
The public context list must not be the only copy of the snapshot. A missing,
unrecognized, or altered package-owned token is an integrity failure, not a
request to skip validation.

Structural cache entries must retain the parent context identity and retained
mask contract. A parent snapshot failure invalidates the entire route before
any retained-subtree evidence is reused. Cache hits do not weaken the
external snapshot check.

### K.3 Required rejection behavior

At minimum, these mutations must be rejected on the next prepared public
call:

- any change to tip.label, including a legal unique delimiter collision;
- any endpoint in edge;
- any branch value in edge.length;
- any change to Nnode;
- any protected shape or length change.

The error must clearly request a fresh prepare_tree() call. No estimator
result, warning, status, or stale cache value may be returned after rejection.

## L. Old-context compatibility policy

V2 prepared contexts carry a private schema marker equivalent to:

    context_schema_version = 2
    canonical_contract_version = 2
    protected_snapshot_version = 2

The marker, token, and snapshot are package-owned state. A V1 context that
contains only the old delimiter fingerprint, an absent marker, an unknown
marker, or an incompatible package-owned token must not be silently
interpreted as V2.

The required behavior for an old or incompatible context is a clear
structured rejection that asks the user to call prepare_tree() again on the
underlying raw tree. Automatic migration is prohibited because it could
inherit the V1 collision-prone integrity state. No public argument is added
to select a schema or bypass this rejection.

V1 remains the frozen downstream statistic oracle until a separately
authorized V2 prototype has passed all gates. V2 integration is a separate
decision and is not implied by this document.

## Normative construction procedure

The future private prototype must implement the following conceptual stages
without mutating the source tree:

1. Validate the rooted weighted graph and exact tip-label byte contract.
2. Build canonical tip ranks using the byte comparator.
3. Build an explicit acyclic parent/child representation and an iterative
   postorder, then compute one minimum tip rank per node.
4. Order each child list by minimum tip rank and build one exact descriptor
   record per source node using the resulting child references.
5. Traverse from the canonical root in canonical child order to assign
   internal IDs and emit canonical edge records.
6. Transfer each branch length by its biological source parent-child record.
7. Encode C2(T) from the canonical fields.
8. For a prepared context, separately encode and retain S2(T) in
   package-owned state.

No stage may call sort(character) or order(character) for canonical truth. No
stage may materialize complete descendant-tip label vectors. No stage may use
original internal IDs or source edge-row positions as hidden tie breaks.

## Failure semantics and non-repair guarantees

V2 canonicalization must fail closed for:

- non-phylo or malformed input;
- missing, invalid, empty, or duplicate tip labels;
- invalid UTF-8 conversion;
- malformed edge matrix or endpoint values;
- disconnected, cyclic, self-loop, or multiply-parented graph;
- missing, mismatched, negative, or non-finite branch lengths;
- inconsistent Nnode;
- missing or ambiguous structural root;
- an exact descriptor tie that violates the uniqueness proof.

Failure must not silently:

- reroot a biological tree;
- resolve or remove a polytomy;
- collapse or create a unary node;
- delete a tip;
- alter a label;
- jitter, round, or replace a branch length;
- accept a hash collision as identity;
- use the source node ID or source edge-row position as a tie break.

Existing public failure classes, warnings, statuses, estimator definitions,
RNG streams, thread policy, and numerical tolerances remain unchanged at
integration time. If a representation-only wording must change, that is a
separate documentation patch and not part of V2-A.

## Complexity and storage contract

For V nodes, E edges, and L total UTF-8 label bytes, the prototype must
record:

- one descriptor entry per source node, plus any exact intern entries;
- exactly one child-reference entry per biological edge;
- one canonical edge record per biological edge;
- one byte sequence per stored label, not one copy per ancestor.

The persistent structural representation must be O(V + E + L) in storage up
to constant-size table metadata. It must not contain:

- one complete descendant-label vector per internal node;
- an n_tip by n_tip matrix or bitset;
- a quadratic persistent descriptor key;
- a delimiter-concatenated descendant payload;
- a probabilistic hash used as identity or ordering truth.

Explicit traversal stacks and bounded sort workspaces are permitted. A
prototype must report descriptor entries, child-tuple entries, peak temporary
proxy storage, and any optional digest-hint storage. These counts are
complexity evidence, not a performance claim.

## Exact prototype proof obligations

The prototype may be reported as passing only if all of the following are
demonstrated with exact comparisons:

1. canonical descriptor equality is structural and handle-allocation
   independent;
2. canonical child order is independent of locale, source node IDs, and edge
   rows;
3. canonical numbering and edge rows are deterministic;
4. every canonical edge retains its biological branch association;
5. C2(T) is byte-identical for representation-equivalent valid inputs;
6. S2(T) rejects every protected mutation before cache reuse;
7. no digest collision can authorize identity or mutation acceptance;
8. no stored descendant payload is quadratic;
9. the source input is byte-for-byte unchanged;
10. frozen K, lambda, D, Delta, and ACE behavior remains within the existing
    scientific parity contract.

## Required equivalence fixtures

The audit/prototype fixture matrix must include all of the following:

| ID | Fixture |
| --- | --- |
| A | balanced, random, and pectinate rooted trees |
| B | shuffled source edge rows with branch rows shuffled together |
| C | safe internal-node renumbering |
| D | alternate valid root numbering representation |
| E | rooted polytomy |
| F | two-tip tree |
| G | delimiter labels "a", "b\rc", "a\rb", "c" |
| H | prefix labels |
| I | punctuation labels |
| J | numeric-looking labels |
| K | accented and other valid Unicode labels |
| L | every supported LC_COLLATE available in the test environment |
| M | heterogeneous branch lengths, including distinct equal-prefix values |
| N | malformed and failure fixtures |

The delimiter fixture must verify both clade descriptors and the protected
snapshot encoding. It must include the previous V1 collision and prove that
the V2 length-delimited encoding distinguishes the two legal label vectors.

For every valid representation-equivalent pair, compare with identical():

- canonical edge matrix;
- canonical edge-length vector and row association;
- canonical tip labels;
- Nnode;
- canonical internal mapping;
- canonical edge ordering;
- complete C2(T) bytes.

For invalid fixtures, compare failure class, structured issue code, warning
class, and no-mutation behavior. Do not weaken an existing expectation merely
to accept a V2 representation.

## Invariance gates

The required gate names are:

    V2_RENUMBERING_INVARIANCE
    V2_EDGE_ORDER_INVARIANCE
    V2_LOCALE_INVARIANCE
    V2_DELIMITER_SAFETY
    V2_BRANCH_ASSOCIATION
    V2_MUTATION_REJECTION
    V2_NO_STALE_CACHE_PROPAGATION

All must be PASS before a prototype can be considered. A locale that cannot
be installed must be recorded as NOT_RUN_ENVIRONMENT; it must not be treated
as a pass. The supported locale set must be recorded in fixture provenance.

## Estimator preservation gate

V2 canonical fields are allowed to differ from V1 fields. The biological
input and the statistical definitions are not allowed to differ. Against the
frozen V1 production oracle, the prototype must compare:

- K estimate, controlled permutation values, P, MCSE, and status;
- lambda estimate, log-likelihood, LR, boundaries, and status;
- D observed statistic, random/Brownian accounting, P/MCSE, and status;
- Delta observed fit, MCMC diagnostics, permutation accounting, and status;
- ACE log-likelihood, aligned ancestral likelihoods, rates, SE, and status.

Random comparisons must use identical controlled stochastic inputs or an
explicitly paired fixed RNG configuration. Correlation or visual similarity is
not an equality gate. Existing scientific tolerances must be reused; no
tolerance is widened in V2-A.

The prepared mutation fixture must additionally prove:

    protected mutation -> reject before cache lookup
    stale cache numerical propagation -> impossible

## Documentation impact audit

No formal documentation is changed in V2-A. A later integration/documentation
patch will need to review at least:

- README.md;
- USAGE_zh.md;
- NEWS.md;
- man/prepare_tree.Rd;
- man/resolve_tree.Rd;
- method help files describing representation handling;
- docs/index.md;
- inst/EXPERT_REVIEW.md;
- inst/RC_READINESS.md.

After a passing integration, the package could safely state, if tested:

- canonical identity is independent of the supported locale set;
- canonical identity is independent of internal node numbering and edge-row
  order;
- delimiter-containing labels are safe under length-delimited encoding;
- prepared-context mutation is rejected before cache reuse;
- canonicalization changes representation only and does not perform biological
  repair.

The package must not state without separate evidence:

- cross-platform qualification for every OS, compiler, R version, or BLAS;
- invariance of node.label, which is excluded from the current computational
  contract;
- acceptance or repair of biologically invalid trees;
- collision immunity from a digest alone;
- Unicode normalization or equivalence beyond exact UTF-8 bytes;
- numerical equality for branch-length bit patterns that the contract treats
  as different states.

## V2-A decision

This file completes the formal contract design and its proof obligations. It
does not implement the exact prototype or report any gate as passed.

    V2_FORMAL_CONTRACT = DRAFT_COMPLETE
    V2_EXACT_PROTOTYPE = NOT_IMPLEMENTED
    CANONICALIZATION_CONTRACT_V2_PROTOTYPE = NOT_ASSESSED
    CANDIDATE_2 = REJECTED
    PUBLIC_RELEASE_BLOCKER_0_2_0 = UNCHANGED_TRUE

The next authorized work, if separately approved, is an audit-only private
prototype and fixture implementation against this contract. Production
canonicalization, fingerprinting, prepared validation, and all downstream
estimators remain frozen until every gate passes.
