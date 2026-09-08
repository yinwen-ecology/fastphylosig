# fastphylosig 0.2.0 V2-B
# Context Persistence Decision

Status: `DECIDED`

Decision: `SELF_CONTAINED_EXACT_SNAPSHOT_SUPPORT`

Scope: read-only audit of the 0.1.0 prepared context before V2 production
integration. No production or test file was changed by this audit.

## 1. Decision

V2 should preserve the currently observed ability to persist a prepared
context with `saveRDS()` and reuse it after `readRDS()` in a fresh R session.
The V2 context must therefore carry a self-contained, exact protected-state
snapshot that survives R serialization. Validation after deserialization must
compare the current protected tree state with that snapshot before any
structural-cache lookup, numerical-cache lookup, or estimator call.

The V2-A prototype's process-local registry token must not be the sole
authority for an accepted serialized context. A token that exists only in the
R process that created the context would make a valid V2 RDS fail after a
fresh-session load. If a package-private token is retained as an accelerator,
it must be optional or rehydratable; exact serialized state remains the final
integrity decision.

V1 contexts are a separate compatibility case. A context with the old
delimiter-based fingerprint and no recognized V2 schema/snapshot must be
rejected explicitly with an instruction to call `prepare_tree()` again on the
raw tree. No silent migration, silent rebuild, fallback fingerprint, or cache
reuse is allowed. This is a structured version rejection, not an implicit
acceptance of the old integrity contract.

## 2. Existing production evidence

The audit used the installed `fastphylosig` 0.1.0 package with R 4.6.1 on
Windows. The package source at the audited HEAD was
`7d0d6e38756f271e48a74a4ac890f64b186674f3`.

### 2.1 Fresh-session round trip

The following experiment was run in two separate `Rscript` processes:

```r
tr <- ape::rtree(5)
ctx <- prepare_tree(tr)
saveRDS(ctx, file = path)

# fresh R process
ctx <- readRDS(path)
fast_k(ctx, setNames(seq_len(5), ctx$tree$tip.label), test = FALSE)
fast_lambda(ctx, setNames(seq_len(5), ctx$tree$tip.label), test = FALSE)
cache_info(ctx)
```

Observed results:

| Check | Result |
|---|---|
| `serialize(ctx, NULL)` | `raw`, 11,940 bytes for the small prepared fixture after a K call |
| class after `readRDS()` | `fastphylosig_tree`, `list` |
| `fast_k()` after fresh-session load | completed with a valid `phylosig` result |
| `fast_lambda()` after fresh-session load | completed with a valid `phylosig` result |
| `cache_info()` after fresh-session load | completed; structural entry and counters readable |

This was an actual cross-process load, not merely an in-process
`unserialize()` check.

### 2.2 Serialized object composition

The prepared context contains `cache`, `structural_cache`, `numerical_cache`,
and `cache_meta` as R environments. They were serialized and recreated by
`readRDS()`; the cached compiled-tree payload was an R `list`, not an
`externalptr`. No process-local C++ pointer was observed in the prepared
context on this fixture.

This supports self-contained persistence of the current context data, but it
does not by itself provide an integrity contract. V1 still stores only the
delimiter-based `fingerprint` and has no schema or protected-snapshot version.

### 2.3 Mutation gate after load

In independent fresh-session loads, each of the following mutations was made
before calling `fast_k()`:

| Mutated field | Observed result |
|---|---|
| `tree$tip.label` | rejected before estimation |
| `tree$edge` | rejected before estimation |
| `tree$edge.length` | rejected before estimation |
| `tree$Nnode` | rejected before estimation |

The observed message was:

```text
the prepared tree was modified after caching; call prepare_tree(tree) again.
```

Thus the current V1 route does not return a stale K result for these ordinary
post-load mutations. A `node.label`-only change remained accepted and the K
result was unchanged, matching the documented non-computational metadata
boundary.

## 3. Source and contract audit

### 3.1 Production implementation

- `R/prepare_tree.R:17-33` accepts an existing prepared context and validates
  it, or constructs a new context.
- `R/prepare_tree.R:88-90` creates the cache environments.
- `R/prepare_tree.R:111-132` stores the current context fields, including the
  old `fingerprint`, but no V2 context, canonical-contract, or snapshot
  version marker.
- `R/prepare_tree.R:155-164` defines the V1 fingerprint by collapsing labels,
  edge endpoints, formatted branch lengths, and `Nnode` with delimiters. This
  is not an unambiguous V2 encoding.
- `R/prepare_tree.R:559-575` validates only the V1 context shape and the
  equality of the current V1 fingerprint to `ctx$fingerprint`.
- `R/check_tree.R:85-92` contains an internal compatibility comment for
  older serialized contexts with absent inspection evidence, but this is not
  a documented cross-session persistence guarantee and does not establish a
  V2 schema policy.

### 3.2 Documentation

`README.md:124-130`, `USAGE_zh.md:107-110`, and `man/prepare_tree.Rd` describe
the prepared object as a snapshot and document the protected fingerprint
fields. They do not promise `saveRDS()`/`readRDS()` persistence, mention
fresh-session reuse, or define old-context migration. `NEWS.md` describes
prepared caches but contains no serialization contract.

The wording is therefore compatible with preserving persistence, but it does
not make persistence an existing explicit public guarantee.

### 3.3 Tests

No test under `tests/testthat/` asserts a `saveRDS()`/`readRDS()` prepared
context contract. Serialization calls found in the repository are fixture
cloning or PSOCK benchmark plumbing, not prepared-context persistence tests.

The absence of a test means the behavior was previously unqualified, not that
it is safe to break silently. The observed public object is serializable and
currently reusable across processes.

## 4. V2 persistence requirements

The production integration should satisfy all of the following:

1. `prepare_tree()` returns a V2 context containing recognized schema markers
   for context, canonical contract, and protected snapshot versions.
2. The protected snapshot is self-contained, exact, length/dimension aware,
   and encodes the protected fields defined by the V2 contract. It must
   survive `saveRDS()`/`readRDS()` without relying on a process-local registry.
3. On every external prepared-tree public call, schema and snapshot validation
   occurs before structural or numerical cache access and before an estimator.
4. A missing, altered, malformed, unknown, or old snapshot/schema produces a
   clear structured error requesting `prepare_tree()`; it never silently
   falls back to V1 fingerprint semantics.
5. After a fresh-session load, `fast_k()` and `fast_lambda()` must work on the
   unchanged context. At least one post-load mutation test must prove
   rejection before estimator/cache use; the full V2 gate should cover
   `tip.label`, `edge`, `edge.length`, `Nnode`, and a combined mutation.
6. The package-owned exact state must not be represented solely by a
   process-local token. If a private token or registry is used, a fresh session
   must either reconstruct it from the serialized exact state or use a
   self-contained validation path.
7. V1/missing/unknown schema contexts must be rejected clearly and must not be
   auto-migrated, because their delimiter fingerprint can contain the already
   confirmed collision defect.

## 5. Compatibility and safety conclusion

`SELF_CONTAINED_EXACT_SNAPSHOT_SUPPORT = REQUIRED_FOR_V2`.

The evidence does not show that cross-session persistence was an explicitly
documented 0.1.0 promise, and there is no existing regression test. However,
the behavior is real, public-context shaped, and successful in two fresh
processes for both K and lambda. V2 can preserve it without preserving the V1
collision-prone integrity definition. A V2 context can be safely accepted
after deserialization only when its exact snapshot/schema contract is present
and validates; otherwise it must fail closed with a re-preparation message.

`STALE_CACHE_REUSE_AFTER_V2_VALIDATION = PROHIBITED`.

The persistence decision does not authorize any estimator, numerical, RNG,
thread, or cache-policy changes beyond the atomic V2 contract integration.

## 6. Decision gate

```text
V2_CONTEXT_PERSISTENCE_DECISION = SELF_CONTAINED_EXACT_SNAPSHOT_SUPPORT
V2_SERIALIZED_CONTEXT_ACCEPTANCE = YES, only for valid V2 schema + exact snapshot
V1_SERIALIZED_CONTEXT_MIGRATION = NO
V1_OR_UNKNOWN_CONTEXT = STRUCTURED_REJECT_AND_REPREPARE
PROCESS_LOCAL_REGISTRY_AS_SOLE_AUTHORITY = NO
SILENT_REBUILD_OR_FALLBACK = NO
STALE_CACHE_REUSE = NO
```

This audit is complete. Production V2 integration may proceed only after
honoring the requirements above.
