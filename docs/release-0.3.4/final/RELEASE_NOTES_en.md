# fastphylosig 0.3.4

Version 0.3.4 provides clade-scale analyses with `fast_hetero()`: select monophyletic crowns using prespecified nodes or size bounds, calculate local and whole-tree phylogenetic signal, retain coverage, observed sample sizes, missingness, state composition and diagnostics, and display estimates on the tree. Tree annotations distinguish ordinary whole-tree signal tests from heterogeneity-reference tests.

Differences between local estimates describe phylogenetic conservatism patterns at the analysed scale. Signal values are not PVR explained variance; they do not by themselves identify stabilizing selection, niche-conservatism mechanisms or changes in evolutionary rates.

| Method | Recommended scope | Heterogeneity-reference status |
|---|---|---|
| K | Local/whole-tree description, ordinary signal tests and BM-reference dispersion testing | BM on a fixed tree and prespecified disjoint crowns; not a universal null for all homogeneous processes |
| Lambda | Local/whole-tree description and ordinary signal tests | Approximate fitted common-local-crown lambda Gaussian reference, with pooled REML by default; common whole-tree lambda need not imply equal local parameters |
| D | Binary-trait local/whole-tree description and the original random/Brownian signal tests | The shared symmetric-ER heterogeneity extension conditional on crown state counts remains experimental |
| Delta | Multistate local/whole-tree description, state composition, local/overall contrasts, fit/MCMC diagnostics and ordinary signal tests | All Delta heterogeneity tests remain outside the recommended workflow; explicit research interfaces are retained and marked experimental |

`fast_hetero(test = FALSE)` remains the default. Ordinary local/whole-tree signal probabilities and heterogeneity-reference probabilities are reported separately. Single-state Delta is not assigned zero. Successful numerical recovery, optimizer retries, remaining optimizer warnings, MCMC warnings and actual failures remain distinguishable.

This version also:

- allows explicit `nodes` to override automatic clade-size bounds, while retaining observed-species `min_tips` and node/memory guards;
- reuses D child transition mixtures and the objective value paired with Brent's best lambda point, without reducing simulation counts or changing reference models;
- evaluates independent ARD recovery likelihoods at matching rate vectors, retaining the original tolerance and requiring improvement over the same run's unfinished fit;
- removes two unused private encoding helpers and refreshes help, executable examples, guides and manual line wrapping.

The signatures and defaults of all 14 public functions match the previous frozen candidate. Twelve complete fixed-input result/diagnostic/RNG contracts match exactly. Speed gains are workload dependent; no universal speedup is claimed.

An independent Windows check included fresh native compilation, the full test suite, examples and PDF/HTML manuals: 0 errors, 0 warnings, 1 incoming `New submission` note; 13,679 test assertions passed and one historical cross-install test was skipped. A separate development K bridge passed 27 four-mode comparisons. The maintainer confirms a successful Linux check; its logs were not independently reviewed locally. Other platforms are not claimed as validated in this release.

Source installation requires a C++17 toolchain; OpenMP is optional:

```r
install.packages("fastphylosig_0.3.4.tar.gz", repos = NULL, type = "source")
```

The accompanying guides and validation summary document scope, clade selection, starting parameters and evidence limits. These notes do not announce CRAN acceptance; the earlier 0.3.1 submission record is maintained separately.
