# fastphylosig 0.2.0 Stage 2E1A

## Type-7 order-statistic requirement audit

Status: `TYPE7_REQUIREMENT_AUDIT = PASS`

This is a read-only source audit. No production code, tests, or user-facing
documentation were changed. The audited source snapshot is Git
`03d72813f9da08b27f924aa37437ab2715cf8dfe`, package version `0.2.0.9000`.
The companion D contract report is
`benchmarks/stage2e1a/D_CONTRACT_AUDIT.md`.

## Authoritative helper

`src/numeric_utils.h:23-38` defines `quantile_type7_sorted()` for an ascending
`std::vector<double>` of length `n`:

```text
h     = 1 + (n - 1) * probability
lo    = floor(h)
gamma = h - lo
```

The helper returns `values.front()` for `lo <= 1`, `values.back()` for
`lo >= n`, and otherwise computes:

```text
(1 - gamma) * values[lo - 1] + gamma * values[lo]
```

This is the frozen type-7 arithmetic. A future implementation must retain the
same floating-point expression order, endpoint branches, and probability
value. It must not replace this with a different quantile API, a rank
rounding rule, or an epsilon-adjusted interpolation.

## Production call sites

The helper is called at exactly these D threshold sites:

| Function | Input being ordered | Frequency | Post-sort use |
|---|---|---:|---|
| `brownian_threshold_cpp()` in `src/fast_signal.cpp:1150-1170` | each supplied Brownian sample column | once per column | threshold only |
| `brownian_tree_threshold_cpp()` in `src/fast_signal.cpp:1220-1263` | each generated Brownian tip vector | once per simulation | threshold only |
| Brownian loop in `src/d_stream.cpp:371-389` | each generated Brownian tip vector | once per simulation column | threshold only |

In each call site, the vector is sorted, the helper is called, and the
original unsorted sample/tip vector is used for state conversion. The sorted
vector is not used by:

- the observed D contrast;
- the random null;
- the Brownian contrast after thresholding;
- P values or MCSE;
- successful/failed accounting;
- RNG bookkeeping;
- result packaging; or
- diagnostics.

There is therefore no downstream requirement for the complete sorted vector.

## Exact order-statistic count

For an interior probability, type-7 needs only:

```text
lower = v_(lo)
upper = v_(lo + 1)
```

where `v_(k)` is the k-th order statistic in one-based notation. In the source
indices these are `values[lo - 1]` and `values[lo]`. The cases are:

| Condition | Required value(s) |
|---|---|
| `lo <= 1` | minimum / first order statistic |
| `lo >= n` | maximum / last order statistic |
| `1 < lo < n`, `gamma == 0` | lower statistic only; interpolation contributes no second value |
| `1 < lo < n`, `gamma != 0` | two adjacent statistics, lower and upper |

Thus the full sort is not serving a hidden second consumer. It is only a
current implementation method for obtaining one or two adjacent order
statistics for each Brownian replicate.

## State conversion after threshold

After obtaining the scalar threshold, each production path loops over the
original tip order and applies the literal comparison:

```text
binary_state[i] = (original_tip_value[i] < threshold) ? 1 : 0
```

The original order is consequently preserved for the C++ D traversal. A
selection candidate must not threshold the selected or partially rearranged
buffer in place unless it proves that the original tip-order values are still
available unchanged. The threshold comparison is strict; changing it to `<=`
would alter the frozen Brownian tail contract at exact ties.

## Read-only replacement design

The technically feasible design is an exact selection path with no random
choice:

```text
compute h, lo, gamma using the current arithmetic
if lo <= 1: lower = minimum
else if lo >= n: upper = maximum
else:
    select values[lo - 1] exactly
    if gamma != 0: select values[lo] exactly
interpolate with the current type-7 expression
threshold the unchanged original tip vector
```

`std::nth_element` can provide the interior selected values, with a second
selection or a suffix minimum for the adjacent value. For finite values this
is non-probabilistic and preserves the mathematical order-statistic values;
the expected per-draw work is linear rather than `O(n log n)`. The selection
buffer can be a copy, or an explicitly separated scratch buffer, provided the
original tip-order values remain unchanged for thresholding.

This report does not implement or authorize that design. It also does not
claim that a particular standard-library implementation gives a universal
worst-case bound. A future candidate must document its comparator and its
handling of duplicate values, signed zero, and any non-finite input before
claiming exact equivalence.

## Required equality gates for any future candidate

Before authorization, compare old sort and candidate selection at the
threshold-bit and binary-state level for:

1. balanced, random, and pectinate trees;
2. binary and rooted-polytomy compatibility paths;
3. `n = 2`, `20`, and representative production sizes;
4. prevalence values near `0`, `0.5`, and `1`;
5. type-7 endpoints (`p = 0`, `p = 1`);
6. integer type-7 positions (`gamma = 0`);
7. fractional positions (`gamma != 0`);
8. repeated values, constant values, near ties, and large offsets;
9. controlled Brownian state matrices and generated Brownian draws;
10. strict-tail boundary rows, including the fixed one-ULP fixture;
11. observed D, random/Brownian null sums, P, MCSE, and accounting;
12. fixed-seed RNG replay and unchanged draw counts; and
13. input immutability, NA masks, retained species, status, and warnings.

The candidate must retain the current requirement that valid public Brownian
inputs are finite. Direct internal helpers that do not enforce this at their
own boundary must not be used as evidence to weaken the public contract.

## Risk boundary

The selection design changes only how the threshold order statistics are
located. It must not change:

```text
RNG or draw order
prevalence calculation
binary-state mapping
strict < threshold conversion
strict random < observed tail
strict Brownian > observed tail
type-7 h/lo/gamma arithmetic
contrast traversal
P/MCSE/accounting
public API or result structure
```

Possible representation-level risks are non-stable rearrangement, duplicate
ties, signed-zero bit patterns, and accidental loss of the original tip
order. These are correctness gates, not reasons to silently loosen equality
or tail tests.

## Decision

```text
FULL_SORT_REQUIRED_FOR_DOWNSTREAM_D = NO
TYPE7_ORDER_STATISTICS_REQUIRED = AT_MOST_TWO_ADJACENT_VALUES
EXACT_NON_PROBABILISTIC_SELECTION = TECHNICALLY_FEASIBLE
PRODUCTION_IMPLEMENTATION_THIS_STAGE = NO
AUTHORIZATION = PENDING_STAGE_2E1A_WALL_TIME_GATES
```

Only if the formal Stage 2E1A timing gates show that sort plus threshold is at
least 30 percent of complete prepared D time in at least two representative
heavy cells, and the expected conservative end-to-end saving is worthwhile,
may a separate Stage 2E1B prototype be considered. The timing gate must not
be inferred from allocation counts or from this source audit alone.
