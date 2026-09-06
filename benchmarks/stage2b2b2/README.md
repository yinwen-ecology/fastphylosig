# fastphylosig 0.2.0 Stage 2B2B Candidate 2 baseline

This directory contains the Candidate-1 post-freeze direct baseline. The
script does not modify package source, tests, or private bindings. It measures
the current frozen implementation of:

- `.descendant_keys_iterative()`
- `.safe_canonicalize_core()`
- `prepare_tree()`

Each cell uses a deterministic balanced, random, or pectinate tree with fixed
unit branch lengths. A warmup is run for both serialized fixture copies. Formal
calls alternate copies and helper order in one R process. Results are written
to a caller-supplied ASCII temporary output directory as raw timings, median /
IQR summaries, adjacent-size empirical exponents, fixture metadata, and
provenance.

The default formal grid is `n = 2000, 5000, 10000, 20000`, with five repeats
through `n = 5000` and three repeats at `n >= 10000` when the schedule below is
used. The script accepts a bounded schedule so smoke tests never masquerade as
formal evidence.

```text
FASTPHYLOSIG_STAGE2B2B2_REPEATS_BY_N=2000:5,5000:5,10000:3,20000:3
```

Example from an ASCII staging copy:

```text
Rscript benchmarks/stage2b2b2/run_baseline.R C:/tmp/fastphylosig C:/tmp/stage2b2b2-baseline
```

Useful bounded controls:

```text
FASTPHYLOSIG_STAGE2B2B2_QUICK=1
FASTPHYLOSIG_STAGE2B2B2_SHAPES=balanced
FASTPHYLOSIG_STAGE2B2B2_N_GRID=50,100,500
FASTPHYLOSIG_STAGE2B2B2_REPEATS=2
FASTPHYLOSIG_STAGE2B2B2_SOURCE_COMMIT=1685d8f58dea6dfe0a3e3d519d58581b3bbae5b9
```

For the formal baseline, keep one process active and do not run another
benchmark concurrently. The baseline must be rerun after any Candidate-2
production change; Stage 2B2A timings are historical context only and are not
substituted for this direct baseline.

The accepted post-Candidate-1 baseline is retained in `results/baseline/`.
Candidate 2 failed its exact-ordering proof gate before prototype or production
integration, so no after timing exists and the frozen Candidate 1 production
implementation remains unchanged.

`verify_ordering_contract.R` is a non-benchmark smoke for the delimiter
collision that blocks compact descendant-set ordering. It verifies the frozen
collapsed keys, node-ID tie-break, safe canonicalization, and input
immutability. Its accepted output is retained in `results/contract/`.
