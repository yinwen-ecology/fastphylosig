## Prepared submission comments for fastphylosig 0.3.4

This file is prepared for a future submission; no submission is performed by this release-preparation step. The earlier 0.3.1 submission record is maintained separately. Its current review status is not inferred here.

## Scope

Clade-scale signal calculations are provided by fast_hetero(). K uses a BM-reference dispersion test. Lambda's fitted common-local-crown reference is approximate. The binary D count-conditioned ER heterogeneity extension remains experimental. Delta is recommended for local/whole-tree description and fit/MCMC diagnostics; all Delta heterogeneity tests remain outside the recommended workflow. Explicit experimental interfaces are retained for reproducibility.

## Check results

Independent Windows environment: Windows 11 x64, R 4.6.1, GCC 14.3.0.

R CMD check --as-cran --run-donttest --no-multiarch

0 errors | 0 warnings | 1 note

The NOTE is the incoming feasibility message "New submission". Installation, native code, examples, full tests and PDF/HTML manual checks passed. Testthat reported 13,679 passing assertions, no failures or warnings, and one skipped historical cross-install RNG test requiring separate reference/candidate installations. The build-excluded development K bridge was separately executed from the checkout: 27 four-mode comparisons passed with no tail-decision differences.

After the full check, only line wrapping in two Rd pages was revised. All R, C++ and test files in the final archive are byte-identical to the full-check archive. Rd checks and a documentation-only closing check passed; the final manual was additionally generated directly from final source Rd and visually checked. The closing check did not rerun installation or full tests.

The maintainer reports a successful Linux check. Linux logs, R/compiler versions and exact check flags were not independently recorded in this local bundle and are not invented here. Other platforms are not claimed as tested in this round.

Twelve complete fixed-input result/diagnostic/RNG regressions agree with the previous frozen candidate. All 14 exported function signatures and defaults are unchanged from that candidate. These checks do not assert universal error-rate control beyond the documented statistical models.
