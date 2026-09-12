# fastphylosig 0.2.0 Stage 2E1A

## D frozen contract audit

Status: `SOURCE_CONTRACT_AUDIT = PASS`

Scope: read-only audit of the current D implementation, tests, and user-facing
documentation. No production code, tests, or user documentation were changed.
The source snapshot audited was Git `03d72813f9da08b27f924aa37437ab2715cf8dfe`
(`fastphylosig` `0.2.0.9000`). The direct smoke used the installed package at
R 4.6.1. The existing D test files were inspected as contract evidence. A
`testthat::test_file()` attempt was blocked by the non-ASCII repository path
being unconvertible under the local R locale; this is an environment/path
limitation, not a reported D test failure.

## Frozen production state

The following earlier decisions remain unchanged:

```text
COMMON_PREPARATION_ENGINE_0_2_0 = FROZEN
CANONICALIZATION_CONTRACT_V2 = ACCEPTED_FROZEN
K_OPTIMIZATION_0_2_0 = CLOSED_FINAL
lambda optimization = CLOSED_FOR_0_2_0
Delta new optimization = CLOSED_FOR_0_2_0
```

This audit does not authorize an implementation. It only records the D
contract and the feasibility boundary for a possible exact order-statistic
design.

## Source map

| Contract item | Current source/test evidence |
|---|---|
| Public route, matching, NA grouping, result fields | `R/fast_d.R:1-457` |
| Binary recoding and prevalence | `R/fast_d.R:460-485` |
| Random and Brownian state adapters | `R/fast_d.R:501-550` |
| Null summary, tails, MCSE, accounting | `R/fast_d.R:555-617` |
| Observed and null D contrast sums | `src/fast_signal.cpp:650-930` |
| Brownian sorting and type-7 thresholding | `src/fast_signal.cpp:1150-1263`, `src/d_stream.cpp:320-415` |
| Type-7 arithmetic | `src/numeric_utils.h:23-38` |
| User contract | `README.md`, `USAGE_zh.md`, `man/fast_d.Rd` |
| Regression evidence | `tests/testthat/test-d-*.R`, `test-root-sensitivity-contract.R`, `test-canonicalization-v2-estimator-parity.R` |

## Frozen D contract

### Observed statistic

`fast_d()` first forms a binary state vector after trait-wise NA retention and
matching. `phylo_d_sums_cpp()` computes the observed contrast sum in the first
column. For a binary internal node, the ancestral value is the inverse
branch-length weighted mean of the two child values. The contribution is the
sum of the two absolute deviations from that mean. The effective branch length
is the harmonic combination of the child branch lengths and is added to the
parent edge before the next ancestral step. The resulting value is
`observed` in the result object.

For a rooted polytomy, the compatibility path is used. It makes the existing
reference-column split and inverse-branch-length weighted group calculation
recorded in `src/fast_signal.cpp`; it does not silently resolve the biological
polytomy into a new binary tree. The public R preparation guard rejects unary
nodes, including a unary node exposed by matching or trait-wise NA pruning.

The final D estimate is exactly:

```text
D = (observed - mean_brownian) / (mean_random - mean_brownian)
```

If either null mean is non-finite or the denominator is zero, D is `NA` and a
human-readable note records why. The observed contrast is still retained in
the result when available.

### Random null

There are two production-compatible inputs:

1. `random_states` is an explicit finite `n_species x nsim` matrix. It is
   dimension-checked and passed as numeric state values to the D contrast
   kernel; it is not regenerated or recoded by `.binary_state()`.
2. Without explicit states, `.phylo_d_random_states()` generates one sample
   of the observed, already recoded state vector per simulation. With
   `rnd.bias`, the named weights are reordered to the retained tree tips and
   passed to R `sample()` without replacement (the default size is the state
   vector length). This preserves the current weighted random-association
   contract and its input validation.

The random tail is the literal strict comparison `null < observed`. Equality
is not an exceedance. The vector result calls this value `Pval1`; the additive
unambiguous alias is `P_random`.

### Brownian null

There are likewise two paths:

1. `brownian_states` is an explicit finite `n_species x nsim` state matrix and
   bypasses Brownian generation and thresholding while preserving the supplied
   columns exactly.
2. Otherwise the compatibility path generates Brownian tip values on tree
   edges with `R::rnorm()`, sorts each draw, calculates the threshold, and
   converts the tip values to binary states. Ordinary rooted binary trees with
   no explicit state/permutation/weight controls use `phylo_d_stream_cpp()`,
   which performs the same generation, thresholding, contrast traversal, and
   online accumulation in bounded chunks. Polytomies, controlled inputs, and
   weighted random-null inputs use the compatibility route.

The Brownian tail is the literal strict comparison `null > observed`. Equality
is not an exceedance. The vector result calls this value `Pval0`; the additive
unambiguous alias is `P_Brownian`.

### Binary conversion and prevalence

`.binary_state()` removes missing values before checking the number of states.
It rejects infinite values, more than two distinct states, and a single
remaining state. It applies R factor ordering to the two valid states, stores
the counts in `states_table`, maps the first factor level to numeric `0` and
the second to numeric `1`, and records the first-level proportion as
`prop_state1`. Thus labels such as `10/20` are recoded to `0/1`, but the
factor-level ordering remains the defining mapping for character or numeric
labels.

The Brownian threshold receives that recorded proportion without inversion.
The generated Brownian binary state is `1` when the Brownian tip value is
strictly less than the threshold and `0` otherwise. This is the current source
convention and must not be changed by an optimization candidate.

### Brownian threshold and type-7 arithmetic

The current helper in `src/numeric_utils.h` implements R type-7 arithmetic on
an ascending vector of length `n`:

```text
h     = 1 + (n - 1) * probability
lo    = floor(h)
gamma = h - lo
```

It returns the first value for `lo <= 1`, the last value for `lo >= n`, and
otherwise returns:

```text
(1 - gamma) * values[lo - 1] + gamma * values[lo]
```

The production callers perform `std::sort()` once per Brownian draw and then
call this helper once. The exact same helper is used by the bounded streaming
path and the compatibility threshold path. No alternate interpolation rule,
rounding shortcut, or epsilon is present in the frozen implementation.

### Strict tails, P, MCSE, and accounting

For retained finite null values, `.d_null_summary()` uses:

```text
P_random   = sum(random < observed) / n_successful_random
P_Brownian = sum(brownian > observed) / n_successful_brownian
MCSE       = sqrt(P * (1 - P) / n_successful)
```

The two nulls are summarized independently and are never combined into one P
value. A P of exactly 0 or 1 has MCSE exactly 0. Non-finite values are
excluded from the relevant denominator, reported through `nsim_failed_*`, and
described in `note`/`message`. The streaming summary has the same numerical
formula and reports all requested draws as successful when its validated
finite online summary is valid.

The following fields are part of the current public result contract:

```text
nsim_requested
nsim_successful_random
nsim_successful_brownian
nsim_failed_random
nsim_failed_brownian
Pval1 / P_random
Pval0 / P_Brownian
MCSE_P_random
MCSE_P_Brownian
status, note, message
```

With `test = FALSE`, P and MCSE fields are set to `NA`, but the current code
still runs the null calibration needed for D and retains the successful/failed
counts. Full null vectors are retained only when both `return_sim = TRUE` and
`keep_null = TRUE`.

### RNG and replay

The default streaming route uses R's RNG through Rcpp: `R::runif()` drives
Fisher-Yates random association and `R::rnorm()` drives Brownian edge draws.
The compatibility Brownian generator also uses `R::rnorm()`, while the R
random-state adapter uses R `sample()`. Explicit state/permutation matrices
isolate the corresponding stochastic input and are the preferred controlled
replay mechanism.

RNG draws are serial in `phylo_d_stream_cpp()`; `ncores` is used only by the
independent contrast traversal. The binary batch kernel uses static OpenMP
column scheduling and does not use a shared RNG state. Therefore changing
`ncores` must not be described as changing the RNG stream, although platform
floating-point differences can still affect a statistic at a strict boundary.
Changing a stochastic input, seed, or worker/thread configuration is not an
equivalent replay.

### Rooted polytomies and tree boundaries

The public documentation states that D accepts rooted polytomies through the
compatibility path and rejects unary internal nodes. `phylo_d_has_polytomy()`
routes a parent with more than two children away from the binary streaming
kernel. The R-level preparation layer also requires finite, strictly positive
branch lengths for D, performs matching and canonicalization on a package-owned
copy, and does not guess a biological root. Root placement remains
scientifically meaningful and is covered by `test-root-sensitivity-contract.R`.

The tests cover ordinary binary trees, rooted non-ultrametric polytomies,
noncanonical root representations, malformed/unary topology guards, zero or
negative branch rejection, and raw/prepared parity. The canonicalization V2
tests also verify D outputs and retained species under representation-only
changes.

### NA and retained-species contract

Global species matching is performed before trait-wise NA grouping. Each unique
retained-tip mask is prepared and inspected once by the shared preparation
layer. A retained subset with fewer than two species or an unsafe unary result
does not enter the D kernel; its row receives an explicit note/status. Valid
traits receive `n_species` and `n_removed_na`, while the result metadata keeps
matched and removed species information. Inputs are not modified.

The D-specific tests cover raw/prepared controlled-null parity, retained NA
subsets, fewer-than-two observations, non-finite trait handling, duplicate or
missing species contracts through the shared input tests, and mutation-safe
prepared contexts.

## One-ULP strict-tail limitation

`tests/testthat/fixtures/d_polytomy_tail_boundary_audit.csv` is a fixed audit
of the known IEEE-754 boundary. It uses identical Brownian state bits and an
identical topology. At draw 132, the caper and fast observed sums differ by
one ULP while the null sum is identical:

```text
caper observed = 0x40133f1a53b535dd
fast observed  = 0x40133f1a53b535dc
null           = 0x40133f1a53b535dd
```

The literal `null > observed` result is consequently `FALSE` for caper and
`TRUE` for fastphylosig at that draw. The fixture records 102 versus 103
Brownian exceedances over 199 draws. It explicitly classifies the difference
as IEEE boundary rounding, with no state, topology, or tail-semantics
mismatch. This is a documented machine-precision limitation, not permission
to add an epsilon, alter a tail, or silently coalesce equal values.

## Full-sort usage audit

The production source contains three D Brownian sort sites:

| Site | Sort use | Data used after sort |
|---|---|---|
| `src/fast_signal.cpp:brownian_threshold_cpp()` | sort one supplied Brownian column | only the type-7 threshold; the original column is used for binary conversion |
| `src/fast_signal.cpp:brownian_tree_threshold_cpp()` | sort one generated Brownian tip vector | only the type-7 threshold; `tip_values` is used for binary conversion |
| `src/d_stream.cpp:brownian` loop | sort one generated Brownian tip vector per chunk column | only the type-7 threshold; `brownian_tip_values` is used for binary conversion |

There is no later production use of the sorted vector for D observed sums,
random nulls, P, MCSE, accounting, result packaging, or diagnostics. For
probability `p` and `n` finite values, the type-7 threshold needs at most the
two adjacent order statistics at indices `lo - 1` and `lo` (one statistic at
the endpoint or when `gamma == 0`). Therefore the current full sort is an
implementation route to obtain a small local order-statistic set; it is not a
hidden requirement of any later D calculation.

## Read-only exact-selection design feasibility

An exact non-probabilistic replacement is technically feasible, but was not
implemented or authorized in Stage 2E1A. A bounded design would:

1. preserve the existing `h`, `lo`, `gamma`, endpoint branches, and type-7
   interpolation expression byte-for-byte;
2. for an interior non-integer quantile, select the lower order statistic with
   `std::nth_element` and obtain the upper adjacent statistic by a second
   selection or a suffix minimum;
3. skip the second selection when `gamma == 0` and skip interpolation at the
   endpoints; and
4. immediately threshold the original unsorted Brownian vector, as today.

For finite doubles this preserves the selected order-statistic values and does
not consume RNG, change binary conversion, or alter the D kernel. The expected
per-draw complexity would change from `O(n log n)` to `O(n)` plus a linear scan;
the random null and observed contrast work are unchanged. This is a design
feasibility result only.

Required equality fixtures before any future authorization include all three
tree shapes, binary and polytomy paths, `p = 0` and `p = 1`, integer and
fractional type-7 positions, duplicate and constant values, near-ties, large
offsets, controlled states, exact threshold bits, strict-tail boundary rows,
P/MCSE/accounting, input immutability, and RNG draw-for-draw replay. The
candidate must also specify its handling of signed zero and any non-finite
internal input before claiming bitwise equivalence. `std::nth_element` must
not be treated as a license to change the comparator, tie semantics, or
threshold arithmetic.

The design is therefore:

```text
EXACT_ORDER_STATISTIC_DESIGN = TECHNICALLY_FEASIBLE
IMPLEMENTATION_IN_STAGE_2E1A = NO
AUTHORIZATION_FROM_CONTRACT_AUDIT = NO
```

## Audit conclusion

The current D contract is internally explicit and covered by source and
regression evidence. The only known strict-tail discrepancy is the documented
one-ULP IEEE-754 boundary fixture. Full sorting is confined to Brownian
threshold generation and serves no downstream purpose beyond at most two
adjacent type-7 order statistics. Whether a future exact selection candidate
is worth qualifying remains a performance-gate question; this audit supplies
the correctness requirements but does not make a GO decision.

```text
D_CONTRACT_AUDIT = PASS
FULL_SORT_DOWNSTREAM_USE = NONE
EXACT_SELECTION_READ_ONLY_DESIGN = FEASIBLE
D_CANDIDATE_IMPLEMENTED = NO
PRODUCTION_CODE_CHANGED = NO
```
