# Stage 2D2 K null-engine phase decomposition

This directory contains a private, audit-only harness. It does not modify the
package namespace, the production C++ sources, the public API, or the test
suite. The harness is deliberately a same-translation-unit build:
`stage2d2_phase_harness.cpp` includes the checked-in
`src/k_permutation.cpp` and calls the production-included
`kperm::compute_one()` evaluator directly.

## Measurements

`run_stage2d2_phase_audit.R` measures three separate phases for fixed
prepared-tree fixtures:

1. `generation`: identity initialization plus the exact current internal RNG
   permutation generation. No null-K evaluator is called.
2. `evaluation`: a bounded family of fixed, controlled permutations is built
   before the timer and passed directly to the production-equivalent
   `kperm::compute_one()` evaluator. Permutation creation and R object
   construction are excluded. The requested thread count is measured, so an
   OpenMP-enabled build can expose the serial versus representative two-worker
   evaluator budget.
3. `full`: the current prepared internal-RNG route after tree parsing/cache
   construction. It computes the observed identity, generates each bounded
   chunk, copies generated indices into the chunk buffer, copies each row into
   the worker buffer, evaluates K, performs exceedance accounting, and creates
   P/MCSE vectors. It does not materialize `sim_K`, matching the production
   `return_sim = FALSE` workload. Chunk processing and per-chunk OpenMP
   startup follow `src/k_permutation.cpp`.

The parse-tree and cache-build work is deliberately outside each phase timer.
The returned metadata records the source of the evaluator and RNG path so the
phase timings are not confused with `prepare_tree()` timing.

The harness declares both `Rcpp::plugins(cpp11)` and
`Rcpp::plugins(openmp)`. The returned `openmp_compiled`,
`openmp_max_threads`, and `n_threads_effective` fields are authoritative for
each cell. A requested two-thread cell is valid as a two-thread measurement
only when `openmp_compiled` is `TRUE` and `n_threads_effective >= 2`; otherwise
the runner records `UNSUPPORTED` for that OpenMP gate and never relabels the
cell as a two-thread result.

## RNG contract

The generation routine mirrors the production loop exactly:

```text
identity <- 0:(n - 1)
for each replicate i:
    perm <- identity                 # fresh vector, as in production
    if include_observed && i == 1:
        keep identity                # consumes no R RNG draws
    else:
        for pos = n - 1 down to 1:
            j <- floor(R::runif(0, pos + 1))
            swap(perm[pos], perm[j])
```

Thus with `include_observed = TRUE` and `n >= 2`, the expected draw and swap
count is `(nsim - 1) * (n - 1)`. The first three generated permutations are
returned as a small trace. The R runner compares the generation-only trace,
checksum, and draw count with the corresponding full-pipeline generation for
the same seed. The full null matrix is never retained for this replay check.

`R::runif` is called under `Rcpp::RNGScope`; the R runner resets `set.seed()`
for every measured generation/full call. No C++ or OpenMP RNG is used. This
preserves the current serial RNG stream and keeps generation independent of
thread count.

## Controlled evaluator fixture

The evaluator phase uses at most eight deterministic permutations, prepared
outside the timer and repeated cyclically to reach `nsim`. Pattern 1 is the
identity; the remaining patterns are cyclic shifts and reversed cyclic
shifts. Each row is a valid permutation of `0:(n - 1)`. This is a wall-time
fixture, not a statistical null distribution, and avoids allocating an
`n x nsim` R matrix at the largest requested cells.

## Timing protocol

Formal mode is explicit:

```powershell
Rscript benchmarks/stage2d2/run_stage2d2_phase_audit.R --formal <repo> <out>
```

The formal defaults are `n = 500, 2000, 5000, 10000`, `nsim = 999, 9999`,
requested threads `1, 2`, and ten repetitions. Each fixture is created once
per `n`; every phase is warmed once per cell; phase order alternates between
forward and reverse on paired repetitions; cells are serialized. The output
contains raw rows and median/IQR summaries. Environment variables beginning
with `FASTPHYLOSIG_STAGE2D2_PHASE_` can narrow or extend a run without
changing the source script.

The harness is not a production candidate and has no acceptance claim by
itself. The Stage 2D2 report must combine these phase budgets with the
separate correctness and candidate timing gates before authorizing any
production change.
