# fastphylosig 0.2.0 Stage 2B1 benchmark harness

This directory contains audit-only scripts for the **raw single-preparation
boundary** candidate.  They do not edit package source, alter tests, or
change a loaded package outside an isolated child process.  Structural
traversal reuse, estimator changes, and numerical optimization are outside
this stage.

## Provenance boundary

The Stage 2A accepted/frozen package is the before reference.  Set
`FASTPHYLOSIG_BEFORE_COMMIT` to `017f39df1aa520f71fc2a0c5647bd598ff23e824`.
Set `FASTPHYLOSIG_AFTER_COMMIT` to the exact candidate build commit.  The
script records both values; it never guesses a commit from a library path.
The installed libraries must be built from the same fixed source fixtures and
must be loaded in separate child R sessions.

## Formal benchmark

Run from an ASCII staging directory when the Windows R session cannot resolve
the workspace's non-ASCII path:

```text
rtk Rscript --vanilla raw_single_preparation_benchmark.R BEFORE_LIB AFTER_LIB OUTPUT_DIR
```

The optional fourth argument `quick` is only a smoke run over n=500, 1,000,
and 2,000 with two pairs; it is not release evidence.  The formal run uses
the fixed grid n=500, 1,000, 2,000, 5,000, 10,000, and 20,000.  By default it
uses 20 before/after pairs for n<10,000 and 10 pairs for n>=10,000, matching
the Stage 2A protocol.  `FASTPHYLOSIG_STAGE2B1_NS` and
`FASTPHYLOSIG_STAGE2B1_REPEATS` may be used only when explicitly recording a
non-formal protocol variant.

Each n uses `ape::rtree(seed=20260904+n)`, postorder edge order, and a named
normal trait generated with seed `97531+n_tip`.  A serialized-identical copy
is passed to the other version.  Each route is warmed once for each version;
odd pairs run before then after and even pairs run after then before.  Short
calls use an inner batch (20, 10, 5, 2, 1 calls for increasing n) and report
per-call time.  Child process startup, argument serialization, and the
prepared context construction are outside the timed prepared-call interval;
the raw route's preparation is inside its `fast_k()` interval.

Routes and output files:

- `raw_fast_k`: end-to-end raw public call, including raw preparation.
- `prepare_tree`: independent preparation total.
- `prepared_fast_k`: context constructed before the timer, then one public
  prepared call.
- `stage2b1_raw_preparation_timings.csv`: every formal observation.
- `stage2b1_raw_preparation_summary.csv`: median/IQR per version and n.
- `stage2b1_before_after_summary.csv`: before/after speedup and acceptance
  calculations.
- `stage2b1_provenance.csv`: exact library, R, compiler metadata, and timing
  protocol.

For each version and n the summary defines:

```text
duplicate_residual_ms = raw_median_ms - prepare_median_ms - prepared_median_ms
raw_overhead_ratio = raw_median_ms / (prepare_median_ms + prepared_median_ms)
duplicate_residual_pct = 100 * duplicate_residual_ms / raw_median_ms
```

The paired gate is evaluated at n>=5,000 as raw median speedup >=1.7x,
at n=20,000 as >=2.0x, and for every reported n as residual <=15% and
`raw_overhead_ratio <= 1.20`.  The script reports these booleans; it does not
silently convert a failed gate into a pass.

## Helper call audit

Run separately after a smoke or formal timing run:

```text
rtk Rscript --vanilla raw_helper_call_audit.R BEFORE_LIB AFTER_LIB OUTPUT_DIR
```

This counts one raw and one prepared `fast_k(test=FALSE)` call at n=500 and
n=5,000 for both installed versions.  The prepared context is built before
the wrappers are installed, so its one external integrity boundary is not
mistaken for raw-route preparation work.  Results are written to
`stage2b1_raw_helper_call_counts.csv`; the Stage 2A baseline and the intended
single-boundary contract are recorded in
`stage2b1_helper_call_contract.csv`.

The expected raw reduction is structural, not cosmetic:

| Helper | Stage 2A raw | Stage 2B1 target | Why the remaining call is needed |
|---|---:|---:|---|
| `.inspect_tree_core()` | 3 | 1 | one post-canonicalization readiness pass |
| `.safe_canonicalize_core()` | 2 | 1 | one representation-safe normalization |
| `.canonical_tree_signature()` | 4 | 2 | before/after signatures inside that normalization |
| `.tree_fingerprint()` | 2 | 1 | one newly-created context integrity fingerprint |
| `.prepare_tree_core()` | 1 | 1 | one private compilation boundary |
| `.prepare_tree_subset()` | 1 | 1 | one full subset cache seed |
| `.prepared_tree_subset()` | 1 | 1 | one retained-mask lookup in common analysis |
| matching / table / NA / matrix helpers | 1 each | 1 each | shared analysis route, unchanged |

The canonical signature target is two calls because the existing type of
canonicalization computes signatures before and after normalization.  This
stage does not change that definition.  The exact result, including helpers
whose invocation count can be zero on a particular fixture, must be taken
from the CSV rather than inferred from this table.

## Acceptance boundary

The parent release audit must combine these files with targeted correctness
tests covering raw/prepared K, lambda, D, and Delta parity; input immutability;
root sensitivity; polytomy and unary/branch-length failures; NA masks; species
matching; mutation detection; warning/status behavior; and unchanged public
API, estimator, tolerance, RNG, and thread policy.  This harness alone cannot
declare the production candidate accepted.

No formal benchmark is run as part of creating this directory.  The output
files listed above are produced only when the maintainer explicitly invokes
the scripts with two verified installed libraries.
