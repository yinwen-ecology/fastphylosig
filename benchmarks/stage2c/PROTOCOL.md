# Stage 2C Context Memory, Persistence, and Provenance Protocol

## Scope

This protocol defines the prepared-context memory, serialization,
fresh-process revalidation, and provenance portion of Stage 2C. Estimator
timings are collected by the separate K and method harnesses in this directory;
their evidence is stored beside the context evidence. None of these harnesses
changes production code, tests, estimators, numerical formulas, RNG behavior,
thread policy, or the common preparation engine.

The V2-B frozen state is an input to this audit:

```text
PREPARATION_ENGINE_STAGE1 = ACCEPTED_FROZEN
RAW_ROUTE_SINGLE_PREPARATION_BOUNDARY = ACCEPTED_FROZEN
INSPECTION_INTEGER_ADJACENCY = ACCEPTED_FROZEN
CANDIDATE_2 = REJECTED
CANONICALIZATION_CONTRACT_V2 = ACCEPTED_FROZEN
PUBLIC_RELEASE_BLOCKER_0_2_0 = CLEARED
COMMON_PREPARATION_ENGINE_0_2_0 = FROZEN
```

No Stage 2C candidate may be implemented from this harness. A correctness
anomaly stops interpretation of the associated performance evidence.

## Invocation

The default mode is a small smoke run. It uses balanced trees with 8 and 32
tips and verifies the complete save/read/validation route without launching
the formal grid.

```text
Rscript run_stage2c_context.R LIBRARY OUTPUT_DIR --smoke
```

The formal context grid is explicit and serialized. It uses n = 500, 1000,
5000, 10000, and 20000. It must not be started accidentally by a smoke run:

```text
Rscript run_stage2c_context.R LIBRARY OUTPUT_DIR --formal
```

An explicit formal subset can be supplied with
`FASTPHYLOSIG_STAGE2C_N_GRID=500,1000`, but only together with `--formal`.
Set `FASTPHYLOSIG_STAGE2C_KEEP_RDS=1` to retain the RDS files under
`OUTPUT_DIR/context_rds`; otherwise each temporary RDS is removed after the
fresh-process check. The child timeout defaults to 120 seconds and can be
changed with `FASTPHYLOSIG_STAGE2C_CHILD_TIMEOUT_SECONDS`. The package library
must be an installed
`fastphylosig` library, not the source tree.

## Fixture and order

The harness constructs one deterministic balanced rooted `phylo` fixture per
n. Tip labels are `sp1`, ..., `spN`; branch lengths are deterministic and
strictly positive. The same fixture is used for preparation, memory, saveRDS,
and readRDS. No estimator is called. The tree shape is recorded as
`balanced`; this is a context-cost profile, not a tree-shape performance
claim.

The formal context run is serial. It does not run multiple benchmark R
processes concurrently and does not include a worker or scheduler startup.

## Measurements

`stage2c_context_memory.csv` records, for each n:

| Field | Meaning |
|---|---|
| `prepare_tree_ms` | One `prepare_tree()` wall-clock interval, including context construction. |
| `context_object_bytes` | `utils::object.size(ctx)` for the in-memory R object. |
| `fingerprint_bytes` | Exact byte length of `ctx$fingerprint`. |
| `protected_snapshot_bytes` | Exact byte length of `ctx$protected_snapshot`. |
| `serialized_context_bytes` | Length of an uncompressed R serialization probe. |
| `cache_info_bytes_used` | The package-reported numerical-cache budget counter. |
| `cache_reported_*_bytes` | Entry-level bytes reported by `cache_info()`. |
| `*_cache_payload_bytes` | Sum of `object.size()` for package-owned cache entries. |
| `cache_payload_bytes` | Structural plus numerical cache payload bytes. |
| `saveRDS_ms` | Uncompressed `saveRDS(version=3)` interval. |
| `saveRDS_bytes` | Resulting RDS file size. |
| `readRDS_ms` | Fresh child process readRDS interval. |
| `first_validation_ms` | Fresh child process `cache_info(ctx)` validation interval after read. |
| `read_plus_first_validation_ms` | Sum from child-process read start through first validation. |
| `fresh_process_status` | Child result status; `ok` is required for a valid persistence observation. |

The two cache byte families are deliberately retained. `cache_info_bytes_used`
is not silently presented as total context memory; it is the package's
numerical-cache accounting field. Conversely, the payload sum is not a process
allocator measurement.

All timing values are wall-clock observations from `proc.time()[["elapsed"]]`.
The harness makes no claim about peak RSS, resident set size, allocator high
water mark, garbage-collection overhead, or total process memory. A report
must use the phrase **measured object/cache bytes** unless peak RSS was measured
by a separate approved tool.

## Fresh-process persistence check

For every successful save, the parent starts a separate `Rscript --vanilla`
child. The child loads the installed package, calls `readRDS()`, then calls
the public `cache_info(ctx)` boundary. `cache_info()` is used because it
performs the prepared-context schema, protected snapshot, and integrity checks
before returning cache metadata and does not invoke an estimator. The child
returns timing and status through a temporary result RDS, and its PID is
recorded. The child installs a transient elapsed limit before loading or
validating; a limit event is recorded as `timeout`, the child writes its result
status, and then exits synchronously. The parent waits for that exit and does
not create a background worker. This is a real process boundary, not an
in-process unserialize shortcut.

The persistence gate is:

```text
saveRDS succeeds
fresh_process_status = ok
readRDS + first validation completes
```

`error`, `timeout`, `process_error`, and `missing_result` are retained as
explicit failure statuses with an error message. They are not converted into
missing timings. The timeout bounds R-level read/validation work; a hard OS
kill is not claimed for a process blocked inside non-interruptible native code.
The current context-only route uses R serialization and package validation, so
a timed-out child is expected to return and terminate through the status path
above.

This context harness does not test estimator parity or mutation rejection;
those belong to the V2-B persistence and correctness gates. A failure here is
reported as a failure, not converted into a missing timing.

## Provenance

`stage2c_provenance.csv` and `stage2c_sessionInfo.txt` record the package
version, full source commit when supplied or resolvable by `git rev-parse
HEAD`, R version/platform/architecture, OS, machine, logical and physical
core counts, compiler variables, BLAS/LAPACK identifiers, locale variables,
and OpenMP-related environment fields.

The context-only path does not execute an OpenMP estimator or simulation. Its
`OpenMP_status` is therefore `not exercised by context-only path`; capability
and environment fields are recorded for provenance and must not be reported as
an OpenMP PASS. If the full commit cannot be resolved, the value is
`<unavailable>`; no short or invented hash is accepted as a full commit.

## Formal-run discipline

The formal grid is one serialized run on fixed fixtures. Do not combine
measurements from different package commits, different fixture definitions,
or different serialization settings. Record the exact output directory and
provenance with any report. The harness uses `compress=FALSE` and
`version=3` for comparable file size/time; this setting is part of the
measurement contract and is not a recommendation about user-facing RDS
defaults.

This protocol does not define a reference speedup. Stage 2C speedups must be
computed only from paired timings collected in the same Stage 2C run. If a
separate estimator benchmark reaches a resource or time limit, report
`RESOURCE_LIMIT` and do not extrapolate a speedup.

## Current execution status

The serialized formal grid completed on 2026-09-09. Its CSV files,
session/provenance records, rebuilt NA-safe method summaries, route budget,
and hotspot gates are stored under
`results/2026-09-09-post-v2/`. The corresponding expert decision is recorded
in the repository-root `review.md`.

The original lambda `stage2c_methods_summary.csv` is superseded because an
`NA` grouping key produced an empty file. Use
`stage2c_methods_summary_rebuilt.csv`. The final candidate authorization is
`stage2c_final_candidate_decisions.csv`; the older K all-grid candidate table
remains diagnostic and does not override the specified heavy-workload gate.

```text
STAGE2C = PASS
PRODUCTION_CODE_CHANGED = NO
```
