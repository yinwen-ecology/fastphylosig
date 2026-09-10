# Stage 2D2 K Null `compute_one()` Audit

Status: read-only source audit. No production code, tests, API, or Stage 2D1
files were changed.

Audit snapshot: `57d7ccc` (`src/k_permutation.cpp`). The accounting below is
source-derived; it is not a timing claim and does not authorize a prototype.

## Scope and notation

Let:

```text
n = tree.n_tip
N = tree.n_total
e = N - 1                 # non-root edges in a connected rooted tree
i = N - n                 # internal nodes, including the root
p = ncol                  # trait columns passed to compute_one()
q = max(1, min(trait_chunk, p))
m = ceiling(p / q)        # trait chunks
```

The compiled-tree contract supplies CSR child arrays, `preorder` and
`postorder`. The formulas use `N` and `e`, so they also cover supported
polytomies.

## Call boundary and local workspaces

`fast_k_tree_permutation_cpp()` parses the compiled tree and builds the K
cache once per C++ call (`src/k_permutation.cpp:525-535`). The observed
identity and every null row then call `compute_one()` (`:560`, `:594`, `:621`,
`:675`, `:707`). The internal RNG path generates the permutation outside
`compute_one()`; it does not change this numerical accounting.

At `compute_one()` entry (`:294-310`), four empty local numeric vectors are
created:

| Object | Logical length per active trait chunk | Use |
|---|---:|---|
| `message` | `N * c` doubles | upward Gaussian messages for the current residual convention |
| `state` | `N * c` doubles | downward conditional states |
| `baseline` | `c` doubles | first permuted tip value for each trait; preserves large-offset precision |
| `delta` | `c` doubles | baseline-relative phylogenetic GLS mean offset |

Here `c` is the current chunk width. On the first chunk, each local vector
normally allocates its payload through `assign()` (`:312-317`); later chunks
reuse capacity but still perform the requested fill. Because these vectors
are local to `compute_one()`, a fresh call normally incurs up to four local
payload allocations, even when the previous replicate used the same sizes.
The caller-owned `out` vector is separate and is not counted as one of the
four local workspaces.

For one call, logical payload sizes are:

```text
message + state + baseline + delta = 2*N*p + 2*p doubles
```

The initial fills therefore write `(2*N + 2) * p` zero values across all
chunks (`:314-317`). `state` receives another `2*p` root-zero writes in the
two downward passes (`:355-358`, `:427-430`). `out.assign()` writes `p`
`NaN` values once at entry (`:300-301`), and the final K values overwrite
`p` positions (`:481-484`). Thus the explicit numeric-vector fill/reset
payload is:

```text
zero/fill writes = (2*N + 4) * p
NaN reset writes = p
final result writes = p
```

The zero fill of `message` is fully overwritten by the first postorder pass;
the zero fill of `state` is fully overwritten by the root reset plus each
preorder pass; `baseline` and `delta` are fully overwritten before use. This
is an observation relevant to a reusable-workspace design, not a production
change. The root resets remain logically required when reusing storage.

## Line-by-line numerical passes

For each active chunk of width `c`, the implementation performs the following
work in this fixed order.

| Source lines | Pass | Exact loop work per chunk | Permutation dependence |
|---|---|---|---|
| `319-322` | baseline | `c` direct `X[perm[0], col]` reads and `c` writes | dependent |
| `324-353` | upward pass 1 | postorder over `N`; `n*c` tip reads/subtractions; `e*c` child-edge weighted additions; `(i)*c` internal message divisions | dependent through messages |
| `355-374` | downward pass 1 | root reset `c`; `e*c` preorder state updates | dependent through messages |
| `376-391` | GLS offset | `n*c` indexed trait reads and terminal-branch operations; `c` `delta` writes | dependent; produces phylogenetic mean offset |
| `396-426` | upward pass 2 | same postorder structure as pass 1, on `(x - baseline) - delta` | dependent |
| `427-445` | downward pass 2 | root reset `c`; `e*c` preorder state updates | dependent through residual messages |
| `447-478` | K reductions | `n*c` numerator reads/squares; `e*c` denominator energy, with additional `n*c` tip reads | dependent |
| `479-484` | result | two scalar divisions/checks and one K write per trait | result depends on permutation |

The two upward passes each visit every child edge once in the CSR child loop.
The two downward passes each visit every non-root node once. The final
denominator visits every non-root node once. Therefore, per trait, the
literal structural counts are:

```text
postorder node visits                 2*N
preorder state-update node visits     2*e
preorder denominator node visits      e
child-edge weighted iterations        2*e
qlinear terminal-tip iterations       n
total node-loop visits                2*N + 3*e = 5*N - 3
edge-oriented iterations, if the
qlinear terminal lookup is included    5*e + n
```

The `2*e` child-edge count is separate from the `2*e` downward state updates;
both are real work. The `n` qlinear iterations use the terminal branch
lengths but are not a CSR child loop.

## Indexed trait reads and writes

`tip_value()` (`:287-292`) directly indexes the column-major `X` through
`perm`; no `n`-length permuted trait vector is allocated. Per trait and per
replicate, the calls occur at:

```text
baseline first tip       1
upward pass 1             n
qlinear / delta           n
upward pass 2             n
numerator                 n
denominator tip branch    n
                           --
                           5*n + 1 indexed X reads
```

Across all `p` traits this is exactly `5*n*p + p` `tip_value()` calls, and
the same number of permutation-index reads. The final `out` write count is
`p`; optional `sim_K` writes and exceedance updates occur in the caller loop,
outside `compute_one()` (`:595-603`, `:622-629`, `:677-685`, `:709-715`).

## Floating-operation proxy

Ignoring casts, finite checks, address arithmetic, and compiler-specific
instruction fusion, one trait has the following operation proxy:

```text
upward tip residuals                         2*n subtractions
upward weighted child reductions             4*e (multiply + add, two passes)
upward internal message divisions             2*i divisions
downward state interpolation                  6*e (sub + multiply + add, two passes)
qlinear / GLS offset                          4*n (two sub + divide + add)
delta normalization                           1 division
numerator residual squares                    4*n (two sub + multiply + add)
denominator residual energy                  4*e (sub + multiply + divide + add)
denominator tip residual construction         2*n subtractions
final K normalization                         2 divisions
                                           -------------------------------
                                             12*n + 14*e + 2*i + 3
```

The proxy counts a first `qlinear +=` as an add and treats each scalar
division as one operation; it is a comparison aid, not a hardware FLOP
measurement. In addition, the two downward loops calculate `alpha` for each
non-root node once per chunk (`:362-364`, `:434-436`). That is
`2*m*(i - 1)` tree-only multiplications for internal non-root nodes; tip
`alpha` is the literal constant `1.0`.

## Candidate C invariant classification

### Already hoisted at the C++ call boundary

The following are computed once by `parse_tree()`/`build_cache()` and then
read by every replicate (`:525-535`, `:210-284`):

| Quantity | Classification | Source evidence |
|---|---|---|
| `parent`, `child_ptr`, `children`, `preorder`, `postorder` | `TREE_ONLY` | parsed once into `Tree` (`:108-207`) |
| `branch` and root distances | `TREE_ONLY` | parsed once; root distances are not used by `compute_one()` |
| `aggregate[node]` | `TREE_ONLY` | postorder Schur aggregate built once (`:213-234`) |
| `outgoing[node]` | `TREE_ONLY` | branch transform built once (`:235-246`) |
| `sum_inv` | `TREE_ONLY` | computed once from all-ones interpolation (`:249-269`) |
| `normalization` | `TREE_ONLY` | computed once from trace and `sum_inv` (`:270-282`) |
| `n`, `N`, `p`, `trait_chunk` | `TREE_ONLY`/call configuration | bound at call/compute entry (`:298-302`) |

No `sum(x)` or `sum(x^2)` is computed in this kernel. The input trait set is
unchanged by permutation, but the implementation intentionally evaluates
indexed residuals in the current permutation order; there is no hidden
trait-moment cache to hoist.

### Repeated tree-only quantities inside `compute_one()`

These are genuine residual invariant work, but exact hoisting would need to
preserve the existing division and multiplication operations:

| Quantity | Current repetition | Classification | Exactness note |
|---|---|---|---|
| `alpha = branch[node] * outgoing[node]` | twice per non-root internal node per chunk | `TREE_ONLY` | storing it would remove `2*m*(i-1)` multiplications; values are already derivable from the frozen cache |
| `1 / branch[tip]` in qlinear | once per tip/trait/replicate (`:384-387`) | `TREE_ONLY` coefficient, current division is repeated | replacing division with a precomputed reciprocal is algebraically equivalent but not automatically bitwise equivalent |
| `1 / branch[node]` in denominator | once per non-root node/trait/replicate (`:475-477`) | `TREE_ONLY` coefficient, current division is repeated | same division-versus-multiply risk |
| `qlinear / cache.sum_inv` | once per trait/replicate (`:389-390`) | `TREE_ONLY` denominator with permutation-dependent numerator | reciprocal hoisting changes the operation form and rounding |
| `(num / den) / cache.normalization` | once per trait/replicate (`:481-484`) | `TREE_ONLY` denominator with permutation-dependent `num`, `den` | reciprocal hoisting changes the operation form and rounding |
| `begin`, `end`, `s`, `base`, `pbase` | repeated for each corresponding pass/chunk | `TREE_ONLY`/index metadata | loads/address arithmetic repeat, but no numerical definition changes |

The fixed coefficients `aggregate`, `outgoing`, `sum_inv`, and
`normalization` themselves are already call-level hoisted. The remaining
`alpha` and branch-division work is the only clear tree-only arithmetic
recomputed in the null evaluator, but removing it while retaining exact
floating behavior requires a controlled prototype and equality gate.

### Trait-set invariant under permutation

| Quantity | Classification | Actual source behavior |
|---|---|---|
| `X` dimensions, finite-value validation | `TRAIT_SET_INVARIANT_UNDER_PERMUTATION` | dimensions and finiteness are checked once before the replicate loop (`:527-534`) |
| trait column count and column-major base offsets | `TRAIT_SET_INVARIANT_UNDER_PERMUTATION` | `p` and `n` are fixed; `tip_value()` uses the same `n*col` layout |
| multiset-level `sum(x)` / `sum(x^2)` | mathematically invariant, but not computed | no source line accumulates them; introducing such reductions would be a formula/order change, not a simple hoist |

### Permutation-dependent quantities

`baseline` (`:319-322`), both message fields (`:330-350`, `:402-423`), both
state fields (`:369-374`, `:441-445`), `delta` (`:376-390`), numerator and
denominator (`:447-478`), and the resulting K (`:481-484`) all depend on the
current `perm` through direct `tip_value()` reads. In particular, the
baseline-relative GLS mean is intentionally recomputed for each permutation;
it cannot be replaced by a set-level arithmetic mean without changing the K
definition. The denominator is a residual phylogenetic energy, not a fixed
expected-Brownian denominator.

## Audit decision

```text
FOUR_NUMERIC_WORKSPACES_CONFIRMED = YES
DIRECT_TRAIT_GATHER_IN_COMPUTE_ONE = NO
TREE_CACHE_ALREADY_HOISTED = YES
REPEATED_TREE_ONLY_ARITHMETIC = alpha and branch divisions
PERMUTATION_DEPENDENT_MEAN/ENERGY = YES
CANDIDATE_C_INVARIANT_HOISTING = PARTIAL_EVIDENCE_ONLY
PRIVATE_C_PROTOTYPE_AUTHORIZED_BY_THIS_AUDIT = NO
PRODUCTION_CODE_CHANGED = NO
```

The source supports a narrowly scoped future exact prototype for reusable
numeric workspaces or preserving the existing divide operations while storing
`alpha`. It does not support hoisting `sum(x)`, `sum(x^2)`, the GLS offset,
or the residual denominator as if they were permutation-invariant. Any C
prototype must first prove ordered-null, P/MCSE, large-offset, NA, and
bitwise/tolerance parity before timing.
