# fastphylosig 0.3.1

This patch improves calibration for rare binary traits, preserves species identities, and makes plots and tree preparation more reliable.

## Fixes

- **D calibration:** restore type-7 quantile interpolation at the lower boundary, including Brownian simulations with a rare first-factor state.
- **Species matching:** preserve explicitly supplied numeric species names in vectors, matrices and data frames.
- **D plots:** use the probability for the selected random or Brownian null and retain strict tail counts through scaling and rounding.
- **Lambda estimation:** align approximately ultrametric classification with ape and use the terminal-edge feasibility bound for unequal tip heights.
- **Large-tree preparation:** support long mask identities with exact cache matching and serialization.
- Align tree readiness checks, strengthen dispatcher regressions, and clarify lambda's serial execution.

## Validation

Local Windows checks passed 9,052 assertions and examples (0 failures, 0 warnings, 1 skipped cross-install test); the cross-install K RNG contract was verified separately. Independent checks cover 36 type-7 quantile cases and 27 Stage 4 comparisons. [Windows, Linux and macOS CI](https://github.com/yinwen-ecology/fastphylosig/actions/runs/37261355664) and the [OpenMP-disabled checks](https://github.com/yinwen-ecology/fastphylosig/actions/runs/37261355621) passed for release commit `9a3cf5dc32e7f5485d92dd11c2bca71141d61f68`. These checks cover code, tests and examples; PDF-manual validation is awaiting the Win-builder report.

The 50,000-tip regression verifies preparation, subset caches, serialization, fallback grouping and K estimation without permutations. The published benchmark figures report the recorded v0.3.0 workloads. [Detailed check record](https://github.com/yinwen-ecology/fastphylosig/blob/main/docs/release-0.3.1/RELEASE_CHECK.md).

## Install

```r
install.packages("remotes")
remotes::install_github("yinwen-ecology/fastphylosig@v0.3.1")
```

Or download [fastphylosig_0.3.1.tar.gz](https://github.com/yinwen-ecology/fastphylosig/releases/download/v0.3.1/fastphylosig_0.3.1.tar.gz) and install with `install.packages("fastphylosig_0.3.1.tar.gz", repos = NULL, type = "source")`. Source installation requires a C++ toolchain.

SHA256: `760F68BC4DAD1ADE0767292A92805E50A5D4F550265DF0228CAB7DE5D432B4FE`.
