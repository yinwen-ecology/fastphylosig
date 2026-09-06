# fastphylosig 0.2.0 Canonicalization Contract V2-A Expert Review

Production oracle: `dd8ec252ba0f2742895403f755951397c0554233`  
Environment: fastphylosig 0.2.0.9000, R 4.6.1 UCRT,
`x86_64-w64-mingw32`  
Scope: formal contract and private exact prototype only. Production
canonicalization, estimators, API, RNG, threading, and tolerances were not
changed. Candidate 2 remains rejected and Stage 2C remains unauthorized.

## 1. V2 formal canonicalization definition

The normative contract is in
`benchmarks/stage2b2c_v2a/FORMAL_CONTRACT.md`. V2 validates one finite rooted
weighted labeled tree, retains public tip IDs and stored `tip.label`, derives
locale-independent label ranks, builds exact bottom-up subtree descriptors,
orders children deterministically, numbers internal nodes in canonical
preorder with root `n_tip + 1`, and emits edges in canonical postorder. Source
internal IDs and source edge-row positions are lookup data only. No biological
rerooting, polytomy resolution, branch repair, or tip deletion is permitted.

## 2. V2 label ordering definition

Each valid label is converted with `enc2utf8()` only to form a comparison key.
The exact UTF-8 bytes are compared as unsigned bytes from left to right; the
shorter byte sequence sorts first when it is a strict prefix. No Unicode
normalization, case folding, trimming, `sort(character)`, `order(character)`,
or active collation is canonical truth. Missing, empty, conversion-invalid,
or byte-duplicate labels fail. The input and output `tip.label` vectors remain
`identical()`; carriage returns, line feeds, punctuation, prefixes, accented
text, and Unicode remain legal label content.

## 3. V2 subtree identity algorithm

A tip descriptor is `Tip(byte-rank)`. An internal descriptor is a typed record
containing arity and references to its canonically ordered child descriptors.
Descriptor IDs are separate from public node IDs: tip descriptor IDs are
byte-ranks, and internal descriptor IDs follow canonical internal traversal.
Children are ordered by minimum descendant tip rank. The minima of disjoint,
non-empty child subtrees cannot tie when tip labels are unique; an observed tie
fails closed. Equality follows the typed descriptor table exactly. No digest,
delimiter string, original node ID, or row position determines identity.

## 4. V2 internal numbering and edge ordering

Public tip IDs remain `1:n_tip`, preserving controlled permutation/null index
semantics. The root is always `n_tip + 1`; other internal IDs follow canonical
preorder through children ordered by minimum UTF-8 tip rank. Edge rows follow
canonical postorder. Each edge length is carried by its original biological
child-edge identity before endpoint remapping. Seven representation pairs,
including shuffled rows, internal renumbering, alternate root numbering, and
combined changes, had `identical()` edge matrices, lengths, label vectors,
Nnode, internal mapping, descriptor tables, ordering metadata, roots, and
fingerprints.

## 5. V2 fingerprint encoding

The V2 fingerprint is the complete exact canonical byte stream, rendered as
hex, not a probabilistic digest. It uses domain tag
`fastphylosig.canonical.tree`, schema 2, typed fields, unsigned big-endian
32-bit lengths/counts/endpoints, explicit matrix dimensions, length-delimited
UTF-8 labels, and big-endian IEEE-754 binary64 branch values. It covers
`tip.label`, `Nnode`, `edge`, `edge.length`, and their shape/order metadata.
Delimiter bytes have no syntactic role, so different valid protected states
cannot collide through concatenation ambiguity.

## 6. Exact mutation-safety mechanism

Mutation safety is independent of the canonical fingerprint. Preparation
stores an exact source-order S2 snapshot under a private registry token and
also places a copy in the prototype context. At validation, schema and token
are checked first, then the public snapshot and current protected fields are
compared with the package-owned snapshot before any cache or estimator route.
Mutations of `tip.label`, `edge`, `edge.length`, and `Nnode` were all rejected.
Input immutability and stored-label preservation also passed for all 14 valid
fixtures.

## 7. Context schema policy

The private prototype requires `context_schema_version = 2`,
`canonical_contract_version = 2`, `protected_snapshot_version = 2`, and a
known package-owned token. A missing, V1, unknown, or altered marker/token is
rejected with an instruction to prepare the tree again. Unknown old contexts
are not migrated and cannot be silently interpreted as V2. Production context
schema is unchanged in this stage.

## 8. Invariance fixture results

All exact gates passed: `V2_RENUMBERING_INVARIANCE`,
`V2_EDGE_ORDER_INVARIANCE`, `V2_LOCALE_INVARIANCE`,
`V2_DELIMITER_SAFETY`, and `V2_BRANCH_ASSOCIATION`. Coverage included balanced,
random, pectinate, polytomous, two-tip, heterogeneous-length, delimiter,
prefix, punctuation, numeric-looking, accented, and Unicode fixtures. Exact
locale equality was observed under C, English US CP1252, German CP1252,
Chinese Simplified CP936, and `en_US.UTF-8` collations available on the host.
Eleven malformed/failure fixtures were rejected without repair.

## 9. Estimator parity results

The gate compared frozen raw production results directly with V2-canonical
representations using the same biological trait mapping and controlled
stochastic inputs. All 234 V2 rows passed existing parity rules: K 45/45,
lambda 36/36, D 54/54, Delta 45/45, and ACE 54/54. The maximum absolute
difference among comparable numerical fields was 0. No tolerance, estimator,
RNG stream, thread policy, status contract, or public result structure was
changed.

## 10. Delimiter mutation blocker

The Stage 2B2C protected-label mutation was reproduced exactly:
`c("a", "b\rc", "a\rb", "c")` became
`c("a\rb", "c", "a", "b\rc")`. Under V2, both the source snapshot and
canonical fingerprint changed, and context validation rejected the mutation
before cache reuse. `V2_MUTATION_REJECTION = PASS` and
`STALE_CACHE_PROPAGATION = IMPOSSIBLE_BY_VALIDATION`. The delimiter-based
release blocker is removed in the private V2 design, but remains a production
blocker until separately authorized integration is completed.

## 11. Complexity and stored payload

Counter evidence covered balanced, random, and pectinate trees at 32, 64,
128, 256, 512, and 1024 tips. For a binary tree with 1024 tips, V2 stored 2046
child-tuple entries and 2047 descriptor entries, with peak child-arity proxy 2.
The same linear counts held for all three shapes; `payload_max / n_tip^2`
decreased to 0.001952. The prototype reports zero materialized descendant-label
vectors and contains no n-by-n matrix, bitset, quadratic persistent key, or
hash-based identity. These are storage-complexity findings, not performance
claims.

## 12. V2 prototype decision

`CANONICALIZATION_CONTRACT_V2_PROTOTYPE = PASS`.

All formal consistency, exact representation, locale, delimiter, branch,
mutation, stale-cache, failure, estimator, immutability, stored-label, and
linear-payload gates passed. Evidence is under
`benchmarks/stage2b2c_v2a/results/v2a/`. This decision authorizes neither
production integration nor Stage 2C. A future integration stage must update
README, `USAGE_zh`, NEWS, `man/`, and generated docs to describe only the
tested representation invariance, locale independence, schema rejection, and
exact mutation detection; it must not claim new performance or broader
cross-platform validation.
