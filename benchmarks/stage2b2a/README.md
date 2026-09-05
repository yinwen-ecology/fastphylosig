# fastphylosig 0.2.0 Stage 2B2A benchmark

This directory contains the read-only preparation-complexity audit for Stage
2B2A. It is deliberately outside the package build and does not modify
production R/C++ code or tests.

## Run modes

From an ASCII staging copy of the repository, run:

```r
Rscript benchmarks/stage2b2a/prepare_complexity_audit.R <repo> <output>
```

The default protocol uses balanced, random, and pectinate fixtures at
`n = 500, 1000, 2000, 5000, 10000, 20000`, one warmup and ten alternating
formal repeats per shape and size. Set `FASTPHYLOSIG_STAGE2B2A_REPEATS`,
`FASTPHYLOSIG_STAGE2B2A_REPEATS_BY_N`, `FASTPHYLOSIG_STAGE2B2A_SHAPES`, or
`FASTPHYLOSIG_STAGE2B2A_N_GRID` to run bounded serialized partitions. A smoke
run can be requested with
`FASTPHYLOSIG_STAGE2B2A_QUICK=1`; it uses `n = 50, 100, 500, 1000, 2000` and
three repeats.

The authoritative local audit in `results/` ran each shape in a separate,
non-overlapping R process. It used five repeats through `n = 5000` and three
repeats at `n = 10000` and `n = 20000`; the latter two sizes were reduced
after isolated pectinate preflight calls reached 27-56 seconds per major
helper. Every result cell has a median and IQR, and no timed call was censored.
`combine_results.R` combines these serialized partitions and recomputes
shape-and-phase empirical exponents without re-running a timing.

Before timing the largest pectinate fixture, probe one helper in an isolated
serial R process. For example, set:

```text
FASTPHYLOSIG_STAGE2B2A_PROBE_HELPER=descendant_keys_iterative
FASTPHYLOSIG_STAGE2B2A_PROBE_SHAPE=pectinate
FASTPHYLOSIG_STAGE2B2A_PROBE_N=20000
FASTPHYLOSIG_STAGE2B2A_PROBE_TIMEOUT=60
```

The helper names are `prepare_tree`, `inspect_tree_core`,
`safe_canonicalize_core`, `canonical_tree_signature`, `edge_children_index`,
`descendant_labels_iterative`, and `descendant_keys_iterative`. A timeout is
reported as censored evidence; it is not silently converted into a timing.

The script is intentionally serial. It records:

- exact source commit and R/platform provenance;
- total `prepare_tree()` timing by shape and size;
- direct timings for `.inspect_tree_core()`,
  `.safe_canonicalize_core()`, `.canonical_tree_signature()`, and the
  traversal helpers used by those functions;
- median, IQR, adjacent-size empirical exponents, and source-derived
  operation counts;
- one non-timed helper-call count audit and a source operation-pattern audit.

The component timings are independent direct helper timings. They are not
expected to add to `prepare_tree()` wall time because some are nested and the
public preparation path also compiles a tree and initializes caches.

`instrumentation_calibration.csv` reports the overhead of the temporary
namespace wrappers used only for call counting. It is not used as an
authoritative timing table. No production binding is left wrapped after a
counting run.
