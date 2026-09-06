# fastphylosig 0.2.0 Canonicalization Contract V2-A audit

This directory contains an audit-only harness for the V2-A formal contract and
exact private prototype. It does not replace the frozen production
canonicalization, does not resurrect Candidate 2, and does not change tests,
estimators, RNG policy, threading, or the public API.

## Run

The main task should create `v2_prototype.R` in this directory before running
the harness. The required private functions are:

```r
v2_canonicalize(tree)
v2_fingerprint(tree)
v2_protected_snapshot(tree)
v2_prepare_context(tree)
v2_validate_context(context)
```

Run from the repository root:

```text
Rscript benchmarks/stage2b2c_v2a/run_v2a_audit.R \
  . benchmarks/stage2b2c_v2a/results/v2a
```

Optional arguments are an installed package library and an explicit source
commit:

```text
Rscript benchmarks/stage2b2c_v2a/run_v2a_audit.R \
  . benchmarks/stage2b2c_v2a/results/v2a C:/path/to/library <commit>
```

The harness is fail-closed. It stops without evidence if the prototype file or
any required interface is missing. Formal results are intentionally not
created by this handoff; the parent task owns the single controlled run.

## Prototype return contract

`v2_canonicalize()` may return a `phylo` object or a list containing a
canonical `phylo` object under `tree`, `canonical_tree`, `phylo`, or
`canonical`. For exact comparison, the result should expose metadata through
`canonical_fields`, `fields`, `fastphylosig_v2_contract`, or
`fastphylosig_v2_metrics` with these values:

```text
edge
edge.length
tip.label
Nnode
internal_mapping
subtree_descriptors
edge_order
root
```

The first four fields may be read from the returned canonical tree. The
canonical internal mapping must be in canonical coordinates, not an
input-node-ID-to-output-node-ID map. Missing required fields are recorded and
fail the exact invariance gates; they are never inferred from a probabilistic
hash. Complexity evidence should be exposed in the result or its metadata as:

```text
child_rank_tuple_entries
descriptor_entries
peak_internal_object_proxy
```

`v2_fingerprint()` must return one deterministic scalar character/raw value or
a list containing `fingerprint`. `v2_protected_snapshot()` must be an exact,
length-aware protected-field snapshot. `v2_prepare_context()` must return a
context containing the protected tree and a private schema/version marker.
`v2_validate_context()` may return `FALSE` or raise an error for invalid,
mutated, or old-schema contexts; invisible success is accepted for an
unchanged context.

## Coverage

The harness creates fixed, in-memory fixtures for:

- balanced, random, and pectinate trees;
- shuffled edge rows;
- internal renumbering, including a changed root number;
- polytomy and two-tip trees;
- delimiter labels `a`, `b\r c`, `a\r b`, `c` without the display spaces;
- prefixes, punctuation, numeric-looking labels, accented and Unicode labels;
- heterogeneous branch lengths;
- malformed, disconnected, cyclic, duplicate-label, missing-label, and bad
  branch-length inputs.

Only `LC_COLLATE` is varied during the locale audit. Unsupported locales are
written as `NOT_RUN_ENVIRONMENT`, never as PASS. All canonical structural
fields and fingerprints are compared with `identical()`.

The frozen production oracle runs fixed-input K, lambda, D, Delta, and ACE
comparisons. The V2 gate compares each V2-canonical representation directly
with the frozen raw reference for the same biological data. Diagnostic
`frozen_raw` representation rows are retained separately and do not define the
V2 gate. K permutations and D null states are supplied explicitly; Delta uses
a fixed seed and a small audit-only MCMC; no benchmark or speed claim is made.

Mutation checks serialize local copies and exercise `tip.label`, `edge`,
`edge.length`, `Nnode`, and the old schema marker. The source fixtures and
production contexts are not modified.

## Evidence files

The controlled run writes:

```text
v2_gate_summary.csv
v2_canonical_details.csv
v2_representation_invariance.csv
v2_locale_invariance.csv
v2_mutation_safety.csv
v2_failure_semantics.csv
v2_estimator_parity.csv
v2_complexity.csv
v2_provenance.csv
```

`v2_gate_summary.csv` reports `V2_PROTOTYPE = PASS` only when exact
renumbering, edge-order, locale, delimiter, branch-association, mutation,
failure, estimator, and non-quadratic payload gates pass. A missing interface,
unsupported fixture, absent complexity counter, stale-cache uncertainty, or
any exact mismatch remains a failure or explicit not-run condition.
