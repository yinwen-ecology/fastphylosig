# Stage 2D1 K Permutation Hot-Loop Audit

Status: read-only source audit. No production code, tests, or public API were
changed. This is a pre-prototype accounting document; it contains no timing
claim and no acceptance of a candidate.

Audit snapshot: `baea4cf3306b2e9c314d446782dfbe0b98263db6`.

## Symbols and accounting convention

Let:

```text
n = number of retained tips passed to the K kernel
i = number of internal nodes
N = n_total = n + i
e = N - 1 = number of non-root edges in a connected rooted tree
p = number of traits in X
q = trait_chunk
b = min(simulation_chunk, nsim)
T = effective OpenMP thread count
```

The common fully bifurcating rooted case has `N = 2*n - 1` and
`e = 2*n - 2`, but the formulas below use `N` and `e` so that supported
polytomies remain covered.

"Payload bytes" below means bytes in the copied or initialized element
payload. A full `int` index copy has `n` elements and `4*n` payload bytes;
the corresponding read-plus-write memory traffic is approximately twice that
number. Numerical tree reads/writes, allocator metadata, cache misses, and
Rcpp object overhead are not included in the lower-bound payload totals.

## Which hot path is used

| Caller state | C++ branch | Random draws | Index storage |
|---|---|---:|---|
| Public K test, `permutations = NULL` | Internal RNG, `src/k_permutation.cpp:631-719` | `(nsim - 1)*(n - 1)` when identity is included | bounded `chunk_perms`, size `b*n` |
| Explicit `permutations` matrix | Controlled serial or OpenMP branch, `:539-630` | none | caller's `nsim*n` matrix plus one worker vector |
| Internal validation with generated controlled matrix | Controlled branch | R-side `sample.int` before C++ | R-side matrix plus one worker vector |

The public R wrapper chooses the internal branch for the normal tree engine
when `permutations` is `NULL`; it creates a controlled matrix only for the
internal validation/dense-oracle route or when the user supplies one
(`R/fast_signal.R:633-653`).

## Call-level work and allocations

These occur once per C++ call or once per simulation block, not once per
replicate:

| Work | Source | Internal RNG | Controlled matrix |
|---|---|---:|---:|
| Parse compiled tree and copy canonical arrays | `:525-534`, `parse_tree()` | once/call | once/call |
| Build aggregate/outgoing cache and normalization | `:535`, `build_cache()` | once/call | once/call |
| Identity vector allocation and fill | `:556-560` | one `n`-int vector/call | one `n`-int vector/call |
| `chunk_perms` allocation | `:637-639` | one `b*n`-int vector/call | none |
| Returned simulation matrix | `:571-576` | `nsim*p` doubles only when requested | same |
| Controlled input validation `seen` and full matrix scan | `:539-553` | not used | one `n`-length seen vector; `nsim*n` index checks |
| OpenMP local exceedance counters | `:581-584` or `:658-661` | `T*p` doubles per internal block | `T*p` doubles once for the controlled call |

The controlled input matrix itself is owned by the caller/R wrapper. The C++
`Rcpp::IntegerMatrix` wrapper does not intentionally create a second full
`nsim*n` matrix, but the worker vector still receives one row at a time.

## Per-replicate index work

### Internal RNG mode

For every replicate, including the identity replicate when present:

| Operation | Source | Element pass | Payload bytes | Heap allocation |
|---|---|---:|---:|---:|
| Reset generated index from identity (`perm = identity`) | `:643-648` | `n` ints | `4*n` | one new `n`-int `perm` vector per replicate |
| Fisher-Yates | `:648`, `fisher_yates()` `:489-495` | `n-1` loop iterations | no retained payload | none; `n-1` `R::runif` draws for shuffled rows |
| Index copy A, generated vector to chunk buffer | `:650-653` | `n` ints | `4*n` | none; `chunk_perms` is reused |
| Index copy B, chunk buffer to worker vector | `:669-674` or `:701-705` | `n` ints | `4*n` | worker `perm` is allocated once per block/thread |

Thus the internal path has three complete index-length passes per replicate
(`12*n` payload bytes), of which copies A and B are the two chunk staging
copies. The identity reset is a required state reset, not a second random
draw. The identity row itself performs the reset and staging copies but skips
Fisher-Yates.

The generated `perm` allocation is the only full index heap allocation that is
currently inside the per-replicate generation loop. The evaluator-side worker
vector is reused across rows in a serial block or by each OpenMP worker in a
block.

### Controlled permutation mode

For every replicate, there is no identity reset and no Fisher-Yates. The
controlled branch performs one full row copy from `perm_matrix` to the
zero-based worker `perm` (`src/k_permutation.cpp:616-621` for serial and
`:588-594` for OpenMP):

```text
one n-int pass = 4*n payload bytes per replicate
one worker-vector allocation per call/thread, not per replicate
zero RNG draws in C++
```

The OpenMP controlled branch currently spells this read as
`perm_matrix(i, r)` inside the loop (`:591-594`). This is an Rcpp
`IntegerMatrix` proxy access despite the nearby comment at `:574-575` saying
the threaded path uses a raw column-major payload. It is a concrete source
level hotspot/safety item for a future prototype: the prototype may evaluate
an immutable raw pointer representation, but must first prove identical row
mapping and thread behavior. This audit does not change it.

## Per-replicate K workspace and trait access

`compute_one()` is called once for the observed identity and once for each
null replicate (`:294`, called at `:560`, `:594`, `:621`, `:675`, and `:707`).
For a trait chunk of width `c`, it creates four local vectors:

```text
message:  N*c doubles
state:    N*c doubles
baseline: c   doubles
delta:    c   doubles
```

At `src/k_permutation.cpp:300-309`, these are empty at function entry and are
assigned on the first chunk. Because the first chunk is the largest, later
chunks normally reuse the vector capacity within that call. Consequently the
usual heap allocation count is four vector allocations per `compute_one`
call, while initialization writes across all chunks total:

```text
2*N*p + 2*p doubles = 8*(2*N*p + 2*p) bytes
```

`out.assign(p, NaN)` at `:300-302` adds `p` reset writes. Its capacity belongs
to the caller's persistent `observed_vec`/`kval` vector and is normally
allocated on the first call for that worker, not on every replicate.

There is no `n`-length temporary permuted trait vector. `tip_value()` at
`:287-292` reads `X[perm[tip] + n*col]` directly. The exact indexed trait
access count per replicate is:

```text
baseline first-tip reads:          p
first upward tip message reads:    n*p
first qlinear reads:               n*p
second upward tip message reads:   n*p
numerator reads:                   n*p
denominator tip reads:             n*p
                                  --------
                                   5*n*p + p double reads
```

The direct read pattern is one reason a future fused-gather candidate is not
automatically a win: it must preserve the current baseline-relative order and
large-offset behavior while avoiding a new full `n*p` buffer.

For each trait, `compute_one()` performs four full tree traversal loops and
three reductions (`:328-366`, `:373-395`, `:400-433`, `:442-475`):

```text
node-loop visits:       2*N + 2*(N-1) + (N-1) = 5*N - 3
tip reductions:         2*n
child-edge iterations:  2*e (two upward passes) + e (denominator) = 3*e
```

These are numerical K work, not additional index copies. They repeat for all
`p` traits, in chunks of width `q`, and are the dominant work in the heavy
prepared K permutation cells identified by Stage 2C.

## Exceedance and result work

After each `compute_one()` call, the hot loop reads each of the `p` K values:

- supplied serial: inclusive-tail check and direct `exceedance[j] += 1` at
  `:622-629`;
- supplied OpenMP: the same check into thread-local `local[tid][j]` at
  `:595-603`, then a `T*p` reduction at `:606-609`;
- internal serial/OpenMP: analogous paths at `:677-685` and `:709-715`.

That is `p` inclusive-tail evaluations and at most `p` scalar count updates
per replicate. The tail helper performs finite checks and an 8-epsilon tie
comparison; no allocation occurs.

When `return_sim` is true, the C++ loop writes `p` doubles directly to the
preallocated column-major matrix through `sim_ptr` (`:597-600`, `:624-626`,
`:679-681`, `:711-713`):

```text
8*p payload bytes written per replicate
one nsim-by-p matrix allocation per call
```

There is no separate chunk-result copy in the C++ kernel. `kval` is read
directly for storage and exceedance. For public vector/matrix results,
`R/fast_signal.R:662-665` later copies each retained simulation column into
the result list; that is O(nsim*p) result packaging after the C++ loop, not a
per-replicate K workspace allocation.

## Lower-bound payload example

The table assumes a fully bifurcating rooted tree (`N = 2*n - 1`), one trait,
internal RNG mode, and `return_sim = FALSE`. It includes the three internal
index passes, the direct indexed trait reads, and message/state/baseline/delta
plus output-vector reset writes. It excludes numerical edge/state traffic,
Rcpp object traffic, allocator metadata, and the optional returned simulation
matrix.

| n | N | index passes | trait reads | workspace + output reset | audited payload total |
|---:|---:|---:|---:|---:|---:|
| 500 | 999 | 6,000 B | 20,008 B | 16,008 B | 42,016 B |
| 2,000 | 3,999 | 24,000 B | 80,008 B | 64,008 B | 168,016 B |
| 5,000 | 9,999 | 60,000 B | 200,008 B | 160,008 B | 420,016 B |
| 10,000 | 19,999 | 120,000 B | 400,008 B | 320,008 B | 840,016 B |
| 20,000 | 39,999 | 240,000 B | 800,008 B | 640,008 B | 1,680,016 B |

For `return_sim = TRUE`, add `8*p` bytes of direct result writes per
replicate, plus one `8*nsim*p` result matrix allocation per C++ call. For a
controlled matrix, replace the internal three-pass `12*n` index payload with
the one-pass `4*n` row copy. The explicit caller-owned matrix itself occupies
approximately `4*nsim*n` bytes before R object overhead.

With the public default `simulation_chunk = 128`, the internal bounded index
buffer is `4*n*min(128, nsim)` bytes:

| n | `chunk_perms` payload |
|---:|---:|
| 500 | 256,000 B |
| 2,000 | 1,024,000 B |
| 5,000 | 2,560,000 B |
| 10,000 | 5,120,000 B |
| 20,000 | 10,240,000 B |

This buffer is bounded by `simulation_chunk`, not `nsim`, in the internal
path. It is released with the C++ call.

## Candidate-relevant findings

### Candidate A: reusable permutation-index buffer

The generation loop allocates an `n`-int vector for every internal-RNG
replicate (`:643`). Reusing that vector is representation-safe only if each
replicate is reset from the identity before the same Fisher-Yates draw loop.
This removes the repeated heap allocation, but it does not remove the
identity reset pass or any random draw. The numerical K tree work and the two
chunk staging copies remain unchanged. Correctness gate: exact ordered
permutation trajectory and all existing K/P/MCSE/accounting tests.

### Candidate B: eliminate duplicate chunk staging copies

The internal path performs two staging copies per replicate:

```text
generated perm -> chunk_perms -> worker perm
```

Both are `n`-int passes. A direct immutable view or a worker-local generation
design could remove one or both, but an OpenMP design must keep row ownership,
simulation order, and no-alias guarantees. The controlled threaded branch also
has the separate Rcpp proxy issue documented above. Correctness gate: exact
controlled rows, same-seed internal rows, `simulation_chunk` boundaries, and
threaded `return_sim` order.

### Candidate C: fused trait gather/null evaluation

The current production kernel already fuses trait access with K evaluation:
there is no full permuted trait allocation or `X[perm, ]` copy. It uses direct
indexed reads in `tip_value()` and stores only tree messages/states. A new C
candidate would therefore have to change the index accessor or numerical loop,
not merely remove an existing gather. It carries higher floating-point and
large-offset risk than A/B and is not justified by this source audit alone.

## Audit decision

The source evidence supports investigating A and B as private exact
prototypes. It does not support treating C as an existing allocation hotspot:
the trait gather has already been fused. No candidate is accepted, no
production binding is changed, and no performance gate is evaluated here.
