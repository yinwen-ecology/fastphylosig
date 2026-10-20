Correct ARD full-tree ancestral probabilities by using the child contribution
`P * L_child` under the row-generator convention, including legal zero-message
boundaries. Correct discrete ACE rate SEs for the `-2 logLik` objective. Released
0.3.4 ARD Delta/ancestral-probability analyses and ER/ARD ACE uncertainty summaries
may need recalculation; fitted rate point estimates are unchanged by these
mathematical corrections. K, lambda, D and ER Delta retain their result/RNG
contracts.

Retain the reviewed exact exponent cache and same-kernel GEMV/copy optimizations.
Freeze the candidate as 0.3.5 with bilingual impact/reanalysis notes, source
manifest, checksum and validation evidence. The preparation itself matches
0.3.4.9003 executable R syntax, native DLL code, exports and defaults.

Statistical scope is unchanged: K's BM reference, lambda's approximate
common-local-crown reference, experimental D heterogeneity inference, and
descriptive Delta with all heterogeneity tests outside the recommended workflow.
Historical ARD calibration is not presented as validation of the corrected
estimator. Failures and all existing diagnostics remain visible.

Validation on Windows/R 4.6.1:

- Full source-package `R CMD check --no-manual --timings`: Status OK; 14,438 pass, 0 fail, 0 warning, 1 optional historical-install skip.
- 25 public workflows, 3,956 objective values, 8 posterior matrices, 8 full MCMC replays and 8 mixed-RNG cases: identical to the reviewed corrected snapshot.
- 14 exports/defaults and all 21 executable R syntax trees unchanged; native `.text` SHA256 identical.
- Independent mathematical regressions and standalone guarded cache lifecycle/reconfiguration/accounting checks pass.

See [candidate materials](https://github.com/yinwen-ecology/fastphylosig/blob/codex/prepare-v0.3.5/docs/release-0.3.5/README.md) for the frozen source
identity and evidence. Existing CI provides separate platform checks.

This is a **draft source candidate**. Do not merge, tag or publish a GitHub
Release as part of this preparation; published 0.3.4 assets remain unchanged.
