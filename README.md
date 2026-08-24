# fastphylosig

[Package website](https://yinwen-ecology.github.io/fastphylosig/) |
[GitHub](https://github.com/yinwen-ecology/fastphylosig) |
[Issues](https://github.com/yinwen-ecology/fastphylosig/issues) |
[Chinese guide](USAGE_zh.md)

`fastphylosig` calculates phylogenetic signal for continuous, binary, and
categorical traits. Most analyses need one function: `fast_signal()`. The
method-specific functions remain available when you want to inspect and
control each step.

## Install

```r
install.packages("remotes")
remotes::install_github("yinwen-ecology/fastphylosig")
library(fastphylosig)
```

Source installation requires a C++ toolchain. OpenMP is optional.

## Start here: one command

Once you have a `phylo` tree and named trait data, choose a method and run:

```r
set.seed(1)
tree <- ape::rtree(30)
trait <- stats::setNames(stats::rnorm(30), tree$tip.label)

fit <- fast_signal(
  tree, data = trait, method = "K", progress = FALSE
)
summary(fit)
```

That is enough for a standard analysis. `fast_signal()` checks the tree,
matches species, handles missing values trait by trait, and runs the selected
calculation.

Choose `method` from the question and trait type:

| Target | `method` | Input | Returns |
|---|---|---|---|
| Blomberg's K | `"K"` | continuous | K estimate; optional permutation test |
| Pagel's lambda | `"lambda"` | continuous | lambda likelihood and LR test |
| Fritz--Purvis D | `"D"` | binary | random and Brownian null comparisons |
| Delta | `"Delta"` | categorical | Delta estimate, MCMC diagnostics, optional permutation test |

The method is always explicit; it is not guessed from trait values. Use
`print(fit)`, `summary(fit)`, `as.data.frame(fit)`, or `plot_signal(fit)` to
inspect the result. Unknown or method-inapplicable arguments produce an error.

## How the workflow runs

![fastphylosig technical roadmap](docs/technical-roadmap.svg)

The tree and trait data feed either workflow. The general workflow performs
the preparation boxes automatically. The advanced workflow exposes the same
preparation order before a method-specific calculation. The editable
[draw.io source](docs/technical-roadmap.drawio) and
[Mermaid source](docs/technical-roadmap.md) are included for website use.

## Advanced use: prepare, then calculate

Use the explicit path when you need an auditable tree/data record or controls
that belong to one method:

```r
checked <- check_tree(tree, signal = "K")
checked

# Safe representation changes only; an explicit outgroup is required if the
# selected method needs a root and the tree is unrooted.
tree_ready <- resolve_tree(tree, signal = "K")
matched <- match_tree_data(tree_ready, data = trait)
ctx <- prepare_tree(matched$tree)
```

After preparation, call the method-specific function directly:

| Function | Trait | Statistical target |
|---|---|---|
| `fast_k()` | continuous | Blomberg's K and permutation test |
| `fast_lambda()` | continuous | Pagel's lambda likelihood and LR test |
| `fast_d()` | binary | Fritz and Purvis D with random and Brownian nulls |
| `fast_delta()` | categorical | Delta, MCMC diagnostics, and permutation test |

```r
k_fit <- fast_k(ctx, x = matched$data, test = TRUE, nsim = 999,
                return_sim = TRUE, progress = FALSE)
lambda_fit <- fast_lambda(ctx, trait, test = TRUE,
                          lambda_profile = TRUE, progress = FALSE)
binary <- stats::setNames(as.integer(trait > stats::median(trait)),
                          tree$tip.label)
d_fit <- fast_d(ctx, binary, test = TRUE, nsim = 99,
                return_sim = TRUE, keep_null = TRUE, progress = FALSE)
categorical <- stats::setNames(rep(c("a", "b", "c"), length.out = 30),
                               tree$tip.label)
delta_fit <- fast_delta(ctx, categorical, test = TRUE,
                        mcmc_sim = 10000, thin = 10, burn = 100,
                        return_sim = TRUE, progress = FALSE)
```

The method-specific controls are deliberately separate: K uses `permutations` and
`simulation_chunk`; lambda uses `lambda_profile`; D uses separate random and
Brownian state controls plus `keep_null`; Delta uses MCMC and permutation
controls. `fast_ace()` separately exposes the ancestral-state reconstruction
used by Delta. The result
classes and legacy fields remain stable.

## Repeated analyses

Prepare a tree once when it will be used repeatedly:

```r
ctx <- prepare_tree(tree)
X <- cbind(trait_1 = trait, trait_2 = trait + stats::rnorm(length(trait)))
fits <- fast_signal(ctx, data = X, method = "K", verbose = FALSE)
cache_info(ctx)
```

Preparation compiles topology and traversal data. Dense covariance,
Cholesky, and spectral payloads are not created by the production K, lambda,
D, or Delta paths. A prepared context is a snapshot: later changes to the
original `phylo` object do not update it. On reuse, its fingerprint checks the
cached tree's tip labels, edge matrix, branch lengths, and `Nnode` value and
rejects changes to those fields. Node labels and other non-computational
attributes are not part of this fingerprint.

## Tree and data details

`check_tree()` is read only. `resolve_tree()` may reorder edges and renumber
internal nodes on a copy when this changes representation only. These
functions do not guess a biological root, select an outgroup, invent or
jitter branch lengths, resolve biological polytomies, ultrametricize a tree,
or silently delete tips. `match_phylo_data()` remains a deprecated alias for
`match_tree_data()`.

Method boundaries are intentional:

- K and D require finite, strictly positive branches. D accepts rooted
  polytomies through its compatibility path but rejects unary internal nodes,
  including unary nodes exposed after matching or NA pruning.
- Lambda permits a zero internal branch only when its transformed terminal
  variances remain valid; zero terminal branches are rejected.
- Delta and `fast_ace()` require a conventionally rooted, fully binary tree
  with finite, strictly positive branches.
- Root placement is supplied by the user and can affect the statistic.
- Measurement-error input (`se != NULL`) is not implemented.

Species identifiers should be unique vector names or matrix/data-frame row
names. Their input order need not match the tree: matched data are reordered
to tree-tip order. An unnamed input is accepted only when its entries or rows
already correspond to `tree$tip.label` order. Missing values are then handled per
trait. Each unique retained-tip mask is validated and cached once. Inputs are
not modified. A wrong positional order without identifiers cannot be detected
automatically and can assign traits to the wrong species.

## Inference and diagnostics

K randomization supports explicit permutation matrices, bounded chunks, and
optional null retention. D reports separate P values and Bernoulli Monte Carlo
standard errors (`MCSE_P_random` and `MCSE_P_Brownian`) for its random and
Brownian nulls, using the finite successful count for each null. These fields
are `NA` when `test = FALSE`. D retains both null vectors only when both
`return_sim = TRUE` and `keep_null = TRUE`; their defaults are respectively
`test` and `return_sim`.

Lambda can retain an approximate likelihood profile. The nominal 95% interval
uses a regular grid and the cutoff
`max(logL) - 0.5 * qchisq(0.95, df = 1)`; threshold crossings are linearly
interpolated. It is a grid/profile approximation, not an exact confidence
interval.

Delta defaults include `test = FALSE`, `nsim = 1000`, `return_sim = test`,
`mcmc_sim = 10000`, `thin = 10`, `burn = 100`, and `model = "ARD"`; therefore
no permutation P value or null vector is produced unless testing is enabled.
With `test = TRUE`, `MCSE_P` uses the number of finite successful permutations,
and `return_sim = TRUE` retains the Delta null vector. Interpret ESS, split
R-hat, `MCSE_Delta`, and, when present, `MCSE_P` together with the point
estimate. Fewer than two species shared globally by the tree and data is a
call-level error. After matching, a trait with fewer than two non-missing
species receives an explicit per-trait insufficient-data result; a completed
trait with 2--19 retained species contributes to one batch warning. Non-finite
diagnostics, ESS below 20, or split R-hat above 1.1 also produce one batch
warning. Warnings never rerun or extend the chain automatically.

## Reproducibility and performance

Use `set.seed()` or controlled permutations for stochastic comparisons.
Record the package version, tree/data fingerprints, simulation counts,
threads, chunk sizes, and whether preparation time is included.

For throughput:

1. Submit multiple traits as a matrix instead of looping in R.
2. Reuse a `prepare_tree()` context.
3. Keep null distributions only when they are needed.
4. Use `ncores` for K randomization, binary-tree D, or sufficiently large
   Delta permutation jobs; small parallel jobs may be slower than serial work.
5. Set `progress = FALSE` to suppress progress stages and species-matching
   messages. Statistical warnings and errors are still signalled. Results
   retain timing metadata.

For stochastic parallel paths, replay requires the seed, explicit simulation
inputs (when supplied), and `ncores` to be held fixed. In particular, changing
Delta's process count repartitions worker RNG streams and is not guaranteed to
reproduce null draws or derived P/MCSE values draw for draw. Explicit
permutation/state matrices isolate the corresponding randomization step, but
Delta's MCMC fits remain stochastic. Record the worker configuration with the
result. Small jobs may be faster in serial.

## Plotting

```r
plot_signal(k_fit)
plot_signal(lambda_fit)
plot_signal(d_fit)
plot_signal(delta_fit)
```

Null-distribution plots require retained simulations. Use
`return_sim = TRUE` for K and Delta; D requires both `return_sim = TRUE` and
`keep_null = TRUE`.

## Compatibility and release status

The package is version 0.1.0. Production estimators, defaults, numerical
tolerances, and the public API are frozen. Local source-tarball checks pass on
Windows with R 4.6.1/Rtools 4.5; the package declares and retains compatibility
with R 4.1.0 and above. Public-release and cross-platform qualification are
deliberately deferred and are not implied by the local check.
