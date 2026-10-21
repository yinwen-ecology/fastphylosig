# fastphylosig 0.3.5

Release date: 2026-10-21. The implementation is frozen from the reviewed
0.3.4.9003 development snapshot. Final release preparation updates documentation
and release metadata, without changing numerical algorithms or public defaults.

## Corrections and their effect on analyses

| Correction | Affected output | Reanalysis guidance |
| --- | --- | --- |
| ARD full-tree ancestral-state downpass | `fast_ace(model = "ARD", marginal = FALSE)` ancestral probabilities and downstream ARD Delta estimates, permutations, local/overall summaries and experimental reference tests | Recalculate these outputs when they were produced with the released 0.3.4 code. Retain the original inputs, controls and seed, and record the corrected package version. |
| Discrete ACE rate-standard-error scale | ER and ARD rate `se` fields | Recalculate uncertainty summaries that used these SEs. Finite SEs increase by a factor of sqrt(2); point estimates of rates and likelihoods are unchanged by this correction. |

For a row-generator transition matrix, `P[i,j]` is the probability of child
state `j` given parent state `i`. The downpass denominator must therefore be
`P %*% L_child`. The correction also permits a zero message for an impossible
ARD parent state when both its posterior and child contribution are zero.
A positive posterior with zero contribution remains a reported failure.

The ACE optimization objective is `-2 * logLik`. Its Hessian divided by two
is the observed information, so rate SEs use
`sqrt(2 * diag(solve(Hessian)))`. Unidentifiable or boundary fits can still
have unavailable SEs and retain their diagnostics. Delta does not estimate
or use these rate SEs.

The ARD probability correction can change entropy, MCMC acceptance and
endpoint jitter. Reusing an old seed therefore does not guarantee the old
ARD Delta value or subsequent RNG state. This is an intentional consequence
of the mathematical correction. The fitted rate/likelihood objective is not
changed by the downpass correction. ER Delta, K, lambda and D retain their
numerical-result and RNG contracts. The standard-error correction applies
to both ER and ARD ACE, independently of that distinction.

## Performance and cache maintenance

The release retains the reviewed exact exponent cache, resets it on
geometry reconfiguration and accounts for hash-node storage. The old
reconfiguration issue was outside the ordinary R path, which configures an
immutable workspace once. Allocation/evaluation failures cannot publish an
incomplete cached column.

The existing Armadillo GEMV dispatch operates directly on pre-sized,
distinct scratch vectors. Cached columns and real outputs are copied into
existing storage. No approximate exponent, new BLAS, optimizer-budget change
or MCMC-kernel change is introduced. Against the immediately preceding
0.3.4.9002 development snapshot, independent complete-call benchmarks measured
4.88% less time for a 512-tip local/overall ARD Delta scan and 3.43% less time
for a 256-tip rare-state pectinate-tree ARD Delta call. These are incremental
Windows measurements, not a speedup estimate against released 0.3.4.

## Statistical scope and interfaces

Public names, arguments, defaults and result classes are retained. The
recommended `fast_hetero()` scope is unchanged:

| Method | Description and reference inference |
| --- | --- |
| K | Local/overall estimates and BM reference dispersion on prespecified disjoint crowns. |
| lambda | Local/overall estimates; approximate common-local-crown Gaussian reference, with pooled REML by default. Crown scale remains part of interpretation. |
| D | Binary local/overall description. Count-conditioned shared-rate symmetric ER heterogeneity inference remains experimental. |
| Delta | Multistate local/overall description, state composition, fit and MCMC diagnostics. All heterogeneity tests remain experimental and outside the recommended workflow. |

Ordinary signal tests remain separate from heterogeneity reference tests.
All failures, numerical recoveries, optimizer warnings and MCMC warnings
remain visible. Historical ARD heterogeneity calibration does not validate
the corrected estimator. Local signal differences can describe conservatism
patterns at the analysed scale; signal alone does not identify an evolutionary
mechanism or rate shift.

## Verification and source identity

See [validation](VALIDATION.md), the frozen [statistical scope](STATISTICAL_SCOPE.json)
and [source manifest](SOURCE_MANIFEST.json). The checked source package,
validation materials and SHA256 checksums accompany the
[0.3.5 release](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.5).
Earlier release tags and assets remain available for reproducibility.
