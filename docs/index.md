---
layout: default
title: fastphylosig
---

# fastphylosig

Fast phylogenetic signal statistics for R, with repeated numerical work
accelerated using Rcpp and RcppArmadillo.

[View source on GitHub](https://github.com/yinwen-ecology/fastphylosig) ·
[Report an issue](https://github.com/yinwen-ecology/fastphylosig/issues) ·
[中文使用指南](https://github.com/yinwen-ecology/fastphylosig/blob/main/USAGE_zh.md)

## Install

```r
install.packages("remotes")
remotes::install_github("yinwen-ecology/fastphylosig")
library(fastphylosig)
```

A working C++ toolchain is required because the package contains compiled
Rcpp/RcppArmadillo code. OpenMP is optional.

## Available statistics

| Function | Trait type | Statistic |
|---|---|---|
| `fast_signal()` | Continuous | Blomberg's K and Pagel's lambda |
| `fast_d()` | Binary | Fritz and Purvis D |
| `fast_delta()` | Categorical | Delta statistic |
| `fast_ace()` | Categorical | Discrete ancestral-state reconstruction |
| `match_phylo_data()` | Any supported table | Tree/data matching report |
| `plot_signal()` | Fitted result | Method-specific diagnostic plot |

## Continuous traits

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

For many traits on the same tree, provide one matrix or data frame so matching,
VCV construction, and decompositions can be reused:

```r
X <- matrix(
  rnorm(50 * 100), nrow = 50,
  dimnames = list(tree$tip.label, paste0("trait_", 1:100))
)

K_table <- fast_signal(tree, X, method = "K")
lambda_table <- fast_signal(tree, X, method = "lambda", test = TRUE)
```

## Binary and categorical traits

```r
x_binary <- as.integer(x > median(x))
names(x_binary) <- tree$tip.label
fit_d <- fast_d(tree, x_binary, nsim = 1000)

x_cat <- sample(letters[1:3], ape::Ntip(tree), replace = TRUE)
names(x_cat) <- tree$tip.label
fit_delta <- fast_delta(
  tree, x_cat, test = TRUE, nsim = 100,
  mcmc_sim = 5000, thin = 10, burn = 100
)
```

## Plot results

```r
plot_signal(fit_k)
plot_signal(fit_lambda)
plot_signal(fit_d)
plot_signal(fit_delta)
```

K and Delta use permutation-distribution plots, lambda uses a profile
likelihood plot, and D overlays its Brownian and random calibration nulls.

## Statistical references

- K and lambda preserve the `phytools::phylosig()` statistical targets for
  `se = NULL` and were validated against phytools 1.0.3.
- D follows the Fritz and Purvis calibration used by `caper::phylo.d()` and
  was validated against caper 1.0.1.
- Delta uses accelerated discrete ancestral-state likelihoods and the
  two-chain MCMC workflow.

Species matching and trait-wise NA removals are reported in result tables.
Raw input objects are not modified.

## Author

Yin Wen
