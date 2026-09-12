# Stage 2E1A D candidate qualification harness

This directory is an audit-only harness for `fastphylosig` 0.2.0 Stage 2E1A.
It does not modify production R or C++ code.  `prototype/d_phase_harness.cpp`
includes `src/d_stream.cpp` into the same translation unit, so Brownian tip
generation, threshold arithmetic, and binary-tree D contrast traversal use the
production implementation.

## Phases

The private exports separate the following measurements:

* `stage2e1a_random_null()`: random-association generation and D evaluation.
* `stage2e1a_brownian_generation()`: continuous Brownian tip generation only.
* `stage2e1a_brownian_phase(..., phase = "sort")`: full sort/copy only.
* `stage2e1a_brownian_phase(..., phase = "sort_threshold")`: sort plus the
  production type-7 threshold calculation.
* `stage2e1a_brownian_phase(..., phase = "binary_d")`: thresholded binary
  conversion plus production D calculation.  Set `return_payload = TRUE` only
  for small correctness fixtures; formal timing leaves it `FALSE` to avoid
  timing Rcpp result materialisation.
* `stage2e1a_complete_d()`: complete random and Brownian null pipeline with
  production-style chunking, tail directions, P, MCSE, and accounting.

These are independent measurements.  They are intentionally not required to
sum to an end-to-end wall time because some phases use fixed states to isolate
sorting or thresholding, while complete calls include bounded chunk setup.

## Correctness order

The correctness runner checks, in order: Brownian replay and downstream binary
states, thresholds, binary conversion, null D values, random `< observed` and
Brownian `> observed` tails, P/MCSE/count accounting, same-seed replay, and
two-thread deterministic replay.  The public `fast_d()` result is compared on
the same prepared binary fixtures.  A small rooted polytomy guard checks that
the production public path remains available for the frozen compatibility
contract without routing that topology through the binary-only harness.

## Timing protocol

The companion formal runner is expected to use fixed balanced, random, and
pectinate trees, warmups, serialized execution, alternating phase order, and
median/IQR summaries.  The authoritative workload is prepared D with
`n = 500, 2000, 5000` (and `10000` only if local resources permit),
`nsim = 999, 9999`, prevalence approximately 0.10, 0.50, and 0.90, and at
least ten paired repetitions for the representative heavy cells.

`ncores = 1` is the qualification path.  `ncores = 2` is observation only;
the harness does not change OpenMP scheduling, RNG streams, threshold
arithmetic, or public API.
