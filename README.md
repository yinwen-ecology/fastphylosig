# fastphylosig

[Package website](https://yinwen-ecology.github.io/fastphylosig/) |
[GitHub repository](https://github.com/yinwen-ecology/fastphylosig) |
[Issues](https://github.com/yinwen-ecology/fastphylosig/issues)

`fastphylosig` provides accelerated R and Rcpp/RcppArmadillo implementations
of selected phylogenetic signal statistics. It is designed for repeated trait
calculations on the same tree while preserving the statistical target of the
reference method.

The package currently provides:

- `fast_signal()` for continuous-trait Blomberg's K and Pagel's lambda.
- `fast_d()` for the Fritz and Purvis D statistic on binary traits.
- `fast_delta()` for Delta on categorical traits.
- `fast_ace()` for the discrete ancestral reconstruction used by Delta.
- `match_phylo_data()` for explicit tree/data matching.
- `plot_signal()` for method-specific diagnostic plots.

K and lambda were validated against `phytools::phylosig()` 1.0.3 with
`se = NULL`. D was validated against `caper::phylo.d()` 1.0.1. A newer
reference-package release should be revalidated before compatibility is
claimed for that release.

## Installation

Install the current GitHub version directly in R:

```r
install.packages("remotes")
remotes::install_github("yinwen-ecology/fastphylosig")
```

For local development from the package root:

```r
install.packages(".", repos = NULL, type = "source")
library(fastphylosig)
```

The package needs a working C++ toolchain. OpenMP is optional. Without OpenMP,
all results remain valid but `ncores > 1` does not accelerate the C++ K and D
kernels.

## Continuous traits

Named vector input returns a compact `phylosig`-like object:

```r
set.seed(1)
tree <- ape::rtree(50)
x <- ape::rTraitCont(tree)

fit_k <- fast_signal(
  tree, x, method = "K", test = TRUE, nsim = 1000
)
fit_lambda <- fast_signal(
  tree, x, method = "lambda", test = TRUE
)

fit_k$K
fit_k$P
fit_lambda$lambda
fit_lambda$P
```

Pass a matrix or data.frame to reuse tree work across traits:

```r
X <- matrix(
  rnorm(50 * 100), nrow = 50,
  dimnames = list(tree$tip.label, paste0("trait_", 1:100))
)

K_table <- fast_signal(tree, X, method = "K", test = FALSE)
lambda_table <- fast_signal(tree, X, method = "lambda", test = TRUE)
```

For K, `ape::vcv.phylo()` and the Cholesky factor are computed once for every
distinct retained-species pattern. The C++ kernel evaluates the same GLS mean,
quadratic form, and normalization used by `phytools::phylosig()`. For lambda,
the same transformed covariance, ML variance estimate, likelihood,
optimization interval, and likelihood-ratio test are used. Measurement-error
models (`se != NULL`) are not implemented.

By default, a single lambda fit stores a profile-likelihood grid and an
approximate 95% interval for plotting. Batch lambda fits skip this extra work
unless `lambda_profile = TRUE`.

## Binary traits

```r
x_binary <- as.integer(x > median(x))
names(x_binary) <- tree$tip.label

fit_d <- fast_d(tree, x_binary, nsim = 1000, ncores = 4)
fit_d$DEstimate
fit_d$Pval1  # random-association null
fit_d$Pval0  # Brownian-threshold null
```

`fast_d()` computes the same contrast-sum calibration as `caper::phylo.d()`.
Brownian states are simulated directly along tree edges, avoiding construction
of a dense VCV matrix for the null simulation. Valid binary labels are recoded
internally to equally spaced 0/1 states, so labels such as `0/1`, `1/2`, and
`10/20` give the same D result.

The returned vector object is `phylo.d`-like but intentionally omits caper's
large `NodalVals` and embedded `comparative.data` fields. Matrix input returns
one summary row per trait.

## Categorical traits

```r
x_cat <- sample(letters[1:3], ape::Ntip(tree), replace = TRUE)
names(x_cat) <- tree$tip.label

fit_delta <- fast_delta(
  tree, x_cat, test = TRUE, nsim = 100,
  mcmc_sim = 5000, thin = 10, burn = 100,
  ncores = 4
)
```

`fast_delta()` follows the single-tree Delta workflow: discrete ancestral-state
probabilities, entropy transformation, and the two-chain MCMC calculation. The
default `ace_engine = "fast"` uses `fast_ace()`; use `ace_engine = "ape"` for
reference checks. `model` supports `"ER"` and `"ARD"`.

When `test = TRUE`, `n_failed_sim` reports permutation fits that failed. The P
value uses only finite permutation results and is `NA` if all permutations
fail.

`fast_ace()` is a focused implementation of
`ape::ace(type = "discrete", method = "ML")` for rooted, fully dichotomous
trees and ER/ARD models. It returns optimized rates, numerical standard errors,
and ancestral likelihoods. Names in `x` must match all tree tips exactly.

## Data matching and missing values

Species names are required as vector names or matrix/data.frame row names. If
the number and order of unnamed observations exactly match the tree tips, the
package reports that it is assuming tree-tip order.

```r
matched <- match_phylo_data(tree, X)
matched$report
matched$tree_tips_removed
matched$data_rows_removed
```

Every result table reports matched species, removed tree tips, removed data
rows, and trait-wise NA removals. Missing values are handled trait by trait.
Traits with the same retained species share tree calculations. Raw input
objects are not modified.

Advanced `permutations` matrices must have exactly `nsim` rows and one column
per retained species. Every row must be a complete 1-based permutation. In a
batch with different NA patterns, one matrix cannot describe groups with
different retained-species counts.

## Plotting

```r
plot_signal(fit_k)
plot_signal(fit_lambda)
plot_signal(fit_d)
plot_signal(fit_delta)
```

- K and Delta show the permutation distribution, observed value, and P tail.
- Lambda shows the likelihood profile, lambda estimate, lambda 0 and 1,
  likelihood-ratio summary, and approximate 95% profile interval.
- D overlays the Brownian-threshold and random-association calibration
  distributions on one D axis, with observed D and the D = 0/1 references.

Distribution plots require `test = TRUE` and `return_sim = TRUE`. For D, a zero
extreme count is displayed as `P < 1/nsim`, not `P = 0`.

## Parallel use

- Observed K and lambda for one trait are dominated by dense tree-matrix work;
  package-level multicore execution usually provides little benefit.
- K randomization and D null simulations use OpenMP through `ncores`.
- Delta permutation tests use separate R worker processes through `ncores`.
- For many traits, passing one matrix is more important than wrapping scalar
  calls in an external parallel loop.

Use at most the number of physical/logical cores available to the R session.
On an 8-core machine, `ncores = 8` is the practical upper setting.

## Validation and benchmarks

Run the package checks from the package root:

```r
devtools::test()
devtools::check(args = "--no-manual")
```

Core tests compare K and lambda with `phytools`, D contrast sums with `caper`,
and ER/ARD likelihoods, rates, standard errors, and ancestral probabilities
with `ape::ace()`.

The retained scripts in `benchmarks/` cover distinct questions:

- `benchmark_chol_solver_K_lambda.R`: Cholesky K/lambda numerical kernels.
- `benchmark_multicore_single_trait.R`: one-trait K and D scaling by core count.
- `benchmark_fast_ace_full_workflow.R`: accelerated ACE and full Delta workflow.
- `benchmark_delta_evolution_scenarios.R`: ER/ARD Delta scenarios and the
  original Python implementation.

Benchmark outputs are written under `benchmarks/results/` and are excluded from
the source package.

## Compatibility boundaries

Alternative engines such as castor, phylosignal, or phylolm are not used by
default. They should only be added after separate validation against the same
K, lambda, log-likelihood, likelihood-ratio, and P-value targets.

The compact result objects preserve the main estimates and tests needed for
analysis and `plot_signal()`. They do not reproduce every closure or internal
field stored by `phytools`, `caper`, or `ape`.
