# Stage 2D1 K RNG Contract Audit

Status: read-only source, documentation, and test audit. No production code,
tests, or public API were changed for this audit.

Audit snapshot: `baea4cf3306b2e9c314d446782dfbe0b98263db6`.

## Scope and conclusion

This document freezes the contract that a Stage 2D1 private prototype must
preserve. It does not authorize an implementation. The current code has two
different permutation sources:

1. The normal tree production path generates permutations in C++ with R's RNG
   (`R::runif`) and a package-local Fisher-Yates loop.
2. The controlled path consumes a caller-supplied integer permutation matrix.
   The R wrapper validates and coerces that matrix before the C++ call; the
   C++ evaluator does not draw random numbers in this mode.

The ordinary public K path uses the first mode when `test = TRUE` and
`permutations = NULL`. The internal validation/dense-oracle route and an
explicit user `permutations` matrix use the second mode.

## Frozen statistical contract

The following facts are established from `R/fast_signal.R`,
`src/k_permutation.cpp`, and `src/numeric_utils.h`:

- The observed value is computed by `compute_one()` with the identity index
  vector. The production K tree engine and its baseline-relative arithmetic
  are not changed by the permutation loop.
- The null is a uniform permutation of the retained trait rows. The R
  preparation layer has already matched species and applied the current
  trait-specific NA mask; the K kernel only receives the retained `X` matrix.
  No excluded tip is reintroduced by the permutation loop.
- In internal RNG mode, the first replicate is the identity when
  `include_observed = TRUE`. Public `fast_k()` passes this flag when no
  controlled matrix is supplied. The remaining `nsim - 1` replicates are
  shuffled.
- In controlled mode, every supplied row is used exactly as supplied. An
  identity row is special only by its values, not by its position. The C++
  call sets `include_observed = FALSE` for this mode.
- The upper tail is inclusive. `inclusive_upper_tail()` first accepts
  `value >= observed`; for finite values below the observed value it accepts
  only an 8-machine-epsilon tie neighborhood scaled by the larger magnitude.
  Non-finite values do not count.
- For a finite observed K, `P = exceedance_count / nsim`. With an included
  identity, MCSE uses the random-only count and `nsim - 1` denominator for its
  Bernoulli estimate, then scales by `nsim`. Controlled mode uses all `nsim`
  supplied rows. These formulas are in `src/k_permutation.cpp` and are also
  mirrored by `.k_permutation_mcse()`.
- `nsim_requested` is the requested number of rows. For a finite observed K,
  `nsim_successful` is currently set to `nsim`; for a non-finite observed K it
  is set to zero and K/P/MCSE are returned as `NA`. The kernel has no separate
  per-null failed-row vector. A prototype must not silently change this
  accounting contract.
- `return_sim` controls storage of the ordered null K matrix in the C++
  binding. The public wrapper normalizes `return_sim` and `keep_null` into the
  internal storage flag. When retained, the rows remain in simulation order.
- Constant and numerically degenerate observed traits remain undefined. The
  baseline-relative subtraction in `compute_one()` is part of the large-offset
  stability contract and must not be replaced by an operation-order-changing
  gather.

## RNG source and draw order

### Public default tree path

`R/fast_signal.R:633-653` leaves `tree_permutations` as `NULL` for the normal
tree engine when the caller does not supply `permutations`. It then calls
`fast_k_tree_permutation_cpp()` with `include_observed = TRUE`.

In `src/k_permutation.cpp:631-719`, the internal path:

- initializes each generated permutation from the identity;
- keeps the identity as replicate 1 when requested;
- calls `fisher_yates()` for each remaining replicate;
- draws exactly `n - 1` uniforms for each shuffled replicate;
- generates all draws serially before any optional OpenMP evaluation;
- evaluates the generated rows in their original global order.

`fisher_yates()` at `src/k_permutation.cpp:489-495` is not `std::shuffle`.
Each iteration calls `R::runif(0, i + 1)` and uses `floor()` to select the
swap index. The source does not call `GetRNGstate()`/`PutRNGstate()` itself;
the Rcpp/R runtime owns the R RNG state. `set.seed()` therefore controls this
stream when it is set immediately before the call and no unrelated code
consumes the stream in between.

For internal mode, the expected draw count is:

```text
include_observed = TRUE:  (nsim - 1) * (n - 1)
include_observed = FALSE: nsim * (n - 1)
```

The identity replicate consumes no random draw. The generated permutation
buffer is reset from the identity for every replicate; a previous shuffled
state is never used as the next replicate's starting state.

### Controlled permutation path

`R/fast_signal.R:640-643` sends an explicit matrix through
`.permutation_matrix()`. That helper validates integer-valued, finite,
1-based, in-range indices and checks that every row is a permutation. The
matrix is then passed to the C++ controlled branch. No R or C++ RNG is used
after the matrix is supplied. This is the preferred oracle route because the
ordered permutation sequence is an explicit input.

The C++ controlled branch validates the matrix again at
`src/k_permutation.cpp:539-553`, then copies one row into a zero-based worker
index vector for each replicate. The second validation is a safety boundary;
it must not be removed by a prototype unless an equivalent private proof is
introduced and separately accepted.

## Thread and chunk contract

The public documentation is deliberately narrower than a blanket
worker-count-invariance promise:

- `README.md:201-224` asks users to record the seed, explicit stochastic
  inputs, thread count, and chunk sizes. It says stochastic parallel replay
  requires seed, explicit inputs when supplied, and `ncores` held fixed.
- `USAGE_zh.md:159-170` gives the same guidance. It explicitly warns that
  Delta worker-count changes can repartition its RNG streams, but makes no
  analogous promise that all public K calls are invariant to changing
  `ncores`.
- `man/fast_k.Rd:20-27` documents `permutations`, `ncores`,
  `simulation_chunk`, and `keep_null`, but does not promise cross-platform or
  cross-worker RNG identity.

The current tests establish a stronger, binding-level behavior for internal
K RNG mode:

- `tests/testthat/test-k-permutation-kernel.R:152-190` resets the same seed and
  compares `simulation_chunk = 7` with `64`, then compares `n_threads = 1`
  with `2`. It requires the ordered `sim_K` sequence to agree to the existing
  tolerance, and requires exact P/exceedance agreement.
- `tests/testthat/test-golden-contracts.R:138-166` supplies a controlled matrix
  and compares one versus two threads, including ordered simulation values,
  P, and exceedance.
- `tests/testthat/test-k-production-validation.R:227-277` compares every
  controlled null K value against the dense GLS oracle and derives P from the
  same inclusive tail.

Therefore the Stage 2D1 prototype contract is:

```text
same seed + same arguments + same ncores + same R/RNG environment
    -> same ordered internal RNG trajectory and derived P/MCSE;
same supplied permutation matrix + same arguments
    -> same ordered controlled trajectory and exact accounting.
```

The prototype must preserve the tested internal binding behavior across the
current `simulation_chunk` and `n_threads` fixtures. It must not add a new
claim that changing `ncores` is a public, cross-platform stochastic
invariance guarantee. Cross-R version, RNGkind, compiler, BLAS, or floating-
point environment equivalence is not promised by the current documentation.

No `std::shuffle`, thread-local RNG, scheduler-dependent draw, or random draw
occurs in the OpenMP regions of `src/k_permutation.cpp`. Internal RNG draws
are completed before the parallel evaluation. This is why current direct
binding tests can be chunk- and thread-invariant without requiring a new
parallel RNG design.

## What is and is not an oracle

The existing controlled tests provide a suitable test-only numerical oracle:

- compare every ordered `sim_K` row, not only final P;
- compare the observed K;
- compare inclusive-tail booleans or their exact exceedance count;
- compare P, MCSE, requested/successful counts, and returned row order;
- compare raw/prepared and NA-mask paths using the same supplied matrix.

For a future private prototype, a same-seed oracle must additionally capture
the ordered permutation index trajectory, or capture an equivalent exact
representation before K evaluation. Comparing only the final P is
insufficient because different trajectories can produce the same exceedance.

## Contract risks a prototype must guard

1. Reusing a permutation buffer must still copy/reset from identity before
   every replicate. Starting from the previous shuffle would change both the
   permutation distribution's realized trajectory and the number/order of
   swaps.
2. Moving Fisher-Yates into an OpenMP region would introduce a new RNG/thread
   contract and is outside Stage 2D1.
3. Directly fusing trait access must preserve the current `tip_value()` index
   order and the per-trait baseline subtraction. A full gathered trait vector
   is not equivalent by default for large-offset inputs.
4. Removing a staging copy must not let a worker mutate storage that is still
   used by another replicate or by `return_sim` packaging.
5. Changing `include_observed`, `n_randomizations`, or the identity
   contribution changes the P/MCSE contract even if the K estimates remain
   unchanged.
6. A user-supplied permutation matrix is already an external allocation. A
   prototype must not create a second `n x nsim` copy merely to remove a small
   worker-vector copy.

## Audit decision

The frozen contract is sufficiently explicit for a controlled private
prototype. Candidate design may investigate allocation reuse and duplicate
index staging, but it must first preserve the exact internal RNG draw order
and the existing controlled-permutation oracle. No production integration is
authorized by this document.

