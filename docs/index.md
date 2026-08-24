---
layout: default
title: fastphylosig
---

# fastphylosig

Phylogenetic signal statistics for continuous, binary, and categorical traits.
Most analyses need one function; advanced users can inspect and control each
preparation step.

[Source](https://github.com/yinwen-ecology/fastphylosig) |
[Issues](https://github.com/yinwen-ecology/fastphylosig/issues) |
[Chinese guide](https://github.com/yinwen-ecology/fastphylosig/blob/main/USAGE_zh.md)

## Install

```r
install.packages("remotes")
remotes::install_github("yinwen-ecology/fastphylosig")
library(fastphylosig)
```

Source installation requires a C++ toolchain. OpenMP is optional.

## Start here: one command

Choose a method and run:

```r
set.seed(1)
tree <- ape::rtree(30)
trait <- setNames(stats::rnorm(30), tree$tip.label)
fit <- fast_signal(tree, data = trait, method = "K", progress = FALSE)
summary(fit)
```

The following traits are used in the advanced examples below:

```r
binary <- setNames(as.integer(trait > stats::median(trait)), tree$tip.label)
categorical <- setNames(rep(c("a", "b", "c"), length.out = 30),
                        tree$tip.label)
```

`fast_signal()` checks the tree, matches species, handles missing values trait
by trait, and runs the selected method. Use `print(fit)`, `summary(fit)`,
`as.data.frame(fit)`, or `plot_signal(fit)` to inspect the result.

## Advanced use: prepare, then calculate

Use the explicit path when tree readiness and species removals need to be
reviewed, or when you need controls belonging to one method.

```r
checked <- check_tree(tree, signal = "K")
checked
tree_ready <- resolve_tree(tree, signal = "K")
matched <- match_tree_data(tree_ready, data = trait)
ctx <- prepare_tree(matched$tree)
```

After preparation, call the method-specific function directly:

| Function | Trait | Purpose |
|---|---|---|
| `fast_k()` | continuous | Blomberg's K |
| `fast_lambda()` | continuous | Pagel's lambda |
| `fast_d()` | binary | Fritz and Purvis D |
| `fast_delta()` | categorical | Delta and MCMC diagnostics |

```r
k <- fast_k(ctx, x = matched$data, test = TRUE, nsim = 999,
            progress = FALSE)
lambda <- fast_lambda(ctx, trait, test = TRUE,
                      lambda_profile = TRUE, progress = FALSE)
d <- fast_d(ctx, binary, test = TRUE, nsim = 99,
            return_sim = TRUE, keep_null = TRUE, progress = FALSE)
delta <- fast_delta(ctx, categorical, test = TRUE,
                    mcmc_sim = 10000, thin = 10, burn = 100,
                    return_sim = TRUE, progress = FALSE)
```

Methods are never inferred from trait values.

## Technical roadmap

<img src="technical-roadmap.svg" alt="fastphylosig technical roadmap" style="max-width: 100%; height: auto;" />

The same tree and named trait data feed both lanes. `fast_signal()` performs
the preparation steps and method dispatch automatically. Advanced users run
the tree and data preparation steps explicitly, then call `fast_k()`,
`fast_lambda()`, `fast_d()`, or `fast_delta()` directly.

The editable [draw.io source](technical-roadmap.drawio) and the
[Mermaid source](technical-roadmap.md) are kept beside this page.

## Preparation boundaries

Representation normalization does not make biological choices. The package
does not guess roots, invent branch lengths, jitter zero branches, resolve
biological polytomies, or silently remove species.

## Delta controls

The public defaults are `mcmc_sim = 10000`, `thin = 10`, and `burn = 100`.
Users should inspect ESS, split R-hat, `MCSE_Delta`, and `MCSE_P`, and increase
the chain or permutation count when precision is inadequate. Small retained
samples and poor chain diagnostics produce aggregated warnings; chains are not
automatically rerun.

For large stochastic jobs, keep the seed, permutation inputs, and worker count
with the saved result. Changing the worker count can repartition parallel
random-number streams, and small jobs may be faster in serial. These are
reproducibility rules; the public estimators and defaults are unchanged.

## 0.1.0 status

Version 0.1.0 declares compatibility with R 4.1.0 and above. The source package
passes local checks on Windows with R 4.6.1/Rtools 4.5. Public-release,
cross-platform, R-devel, and second-BLAS qualification remain outside this
local validation scope.
