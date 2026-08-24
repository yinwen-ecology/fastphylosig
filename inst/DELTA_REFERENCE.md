# Delta reference and production contract

## Methodological target

`fast_delta()` implements the standard, single-tree Delta statistic (`delta_S`)
for categorical traits described by Borges et al. (2019):

- Borges R, Machado JP, Gomes C, Rocha AP, Antunes A. *Measuring
  phylogenetic signal between categorical traits and phylogenies*.
  *Bioinformatics* 35:1862--1869.
  [doi:10.1093/bioinformatics/bty800](https://doi.org/10.1093/bioinformatics/bty800).

Ribeiro et al. (2023) describe the later Python implementation and extend the
method to account for a collection of uncertain trees (`delta_E`):

- Ribeiro D, Borges R, Rocha AP, Antunes A. *Testing phylogenetic signal with
  categorical traits and tree uncertainty*. *Bioinformatics* 39:btad433.
  [doi:10.1093/bioinformatics/btad433](https://doi.org/10.1093/bioinformatics/btad433).

The public function accepts one tree per call. It does **not** implement the
multiple-tree `delta_E` extension, average across a posterior tree sample, or
use PastML. The 2023 paper is therefore methodological and implementation
provenance for the single-tree core, not a claim that every workflow in that
paper is supported.

Borges et al. found poor test sensitivity for small trees and recommended at
least 20 species. `fast_delta()` retains estimates for 2--19 matched,
non-missing species but emits a warning; that warning is not a validation of
small-sample inference.

## Pinned Python artifact

The external source artifact is pinned to:

- repository: [diogo-s-ribeiro/delta-statistic](https://github.com/diogo-s-ribeiro/delta-statistic);
- commit: [`b104a2a7c1b7fed572b782587a811d69719e033a`](https://github.com/diogo-s-ribeiro/delta-statistic/tree/b104a2a7c1b7fed572b782587a811d69719e033a);
- core source: `Delta-Python/delta_functs.py`;
- source-audit date: 2026-08-13.

At that commit, `Delta-webapp/runtime.txt` specifies Python 3.11.2, while
`Delta-Python/requirements.txt` pins Numba 0.57.0, NumPy 1.24.3, SciPy 1.10.1,
PastML 1.9.40, and rpy2 3.5.1. These files are the repository's recorded
environment metadata; they are not a lockfile for this R package and do not
establish that a local end-to-end reference run has completed.

The earlier R implementation linked by Borges et al. is
[mrborges23/delta_statistic](https://github.com/mrborges23/delta_statistic).
No commit or dependency set is pinned here, so no version-specific result is
attributed to that repository.

## Production computation

For each valid trait, `fast_delta()` performs the following computation:

1. It matches species, prunes trait-specific missing values, and requires a
   conventionally rooted, fully binary tree with finite, strictly positive
   branch lengths and at least two observed states.
2. It estimates ancestral-state probabilities by maximum likelihood with the
   package's `fast_ace()` likelihood under the requested `ER` or `ARD` model.
   This is not the pinned repository's PastML ancestral-reconstruction path.
3. It transforms each ancestral probability vector with the selected `LSE`,
   normalized Shannon entropy (`SE`), or normalized Gini impurity (`GINI`).
4. It runs two Metropolis-Hastings chains for the beta-distribution shape
   parameters, using exponential initial values and log-normal proposals, and
   reports `Delta = mean(beta) / mean(alpha)` over the pooled saved draws.
5. When `test = TRUE`, it permutes tip states, repeats ancestral reconstruction
   and MCMC, and reports the inclusive upper-tail proportion
   `mean(null Delta >= observed Delta)` among finite successful fits. There is
   no plus-one correction. Failed permutation fits are excluded from the
   denominator and reported separately. This `P` is a Monte Carlo tail
   proportion (described by Ribeiro et al. as a proxy), not an exact analytic
   P-value.

The public defaults are `mcmc_sim = 10000` iterations per chain, `burn = 100`,
and `thin = 10`, yielding 991 saved draws per chain under the inclusive burn
rule. `mcmc_sim` is the public argument; no separate benchmark-only control is
part of the package API.

ESS, split R-hat, covariance-aware `MCSE_Delta`, permutation `MCSE_P`, chain
accounting, and successful/failed permutation counts are package diagnostics
and extensions. They are not outputs of `Delta-Python/delta_functs.py` and
cannot be required to match fields that the pinned Python core does not
produce. Defaults do not guarantee adequate mixing or permutation precision;
warnings never lengthen or rerun a chain automatically.

## Equivalence boundary

The pinned Python core and the production C++ sampler use the same entropy
families, two-chain beta-shape model, and Delta ratio at the method level.
They are not the same executable implementation. They use different RNGs,
the production sampler evaluates acceptance ratios in log space, and their
iteration/control flow is not stepwise identical. Consequently, equal numeric
seeds do not define matched random streams, and bitwise or draw-by-draw parity
is not part of the production contract.

No Python implementation, benchmark runner, or cross-language result file is
shipped with this package. The pinned commit is retained as methodological
provenance only; completed cross-language core or end-to-end numerical
equivalence is **not claimed**.

A future comparison must distinguish:

- deterministic entropy agreement on fixed probability matrices, with
  boundary perturbations controlled explicitly;
- stochastic agreement of Delta summaries within predeclared Monte Carlo
  uncertainty, not same-seed equality;
- end-to-end sensitivity to ancestral reconstruction (`fast_ace()` versus
  PastML or another declared reference), tree preprocessing, and failed fits;
- permutation-tail agreement using identical explicit permutations while
  separately accounting for sampler randomness; and
- reference-process timing, compilation/warm-up time, ancestral fitting, and
  package timing.

It must record the exact Python executable and installed dependencies and
retain the generated comparison artifact. Matching only an entropy or MCMC
helper is not end-to-end validation.

## Release boundary

The current release keeps the production estimator and its documented
diagnostics stable. Historical optimization experiments and generated audit
tables were removed from the working tree; the current contract is summarized
in `EXPERT_REVIEW.md` and `RC_READINESS.md`.
