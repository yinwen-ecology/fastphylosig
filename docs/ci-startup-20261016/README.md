# CI startup-hook repair

The failure was in the development equivalence test's compiler subprocess.
`R CMD check` sets `R_TESTS=startup.Rs` for its test driver. The installed R
base startup profile sources that path when it is nonempty. The test calls
`Rcpp::sourceCpp()`, which launches `R CMD SHLIB` from a temporary build
directory; the inherited relative path no longer points to the check driver's
startup file. The subprocess exits before the compiler can run.

The development bridge lives under `benchmarks/`, which is excluded from
source archives. An isolated archive check may visibly skip this bridge.
In CI, however, the check directory is inside the repository checkout and
`find_repo_root()` finds the bridge in an ancestor. This exercises the extra
compilation and exposed the environment leak. Crown calculations and the
installed statistical kernels were not the source of the failure.

`compile_stage4_harness()` now temporarily unsets `R_TESTS` for the build and
restores its exact original value on success and error, including the
distinction between an unset variable and an empty string. It preserves the
check driver's startup configuration and leaves four-mode coverage enabled.
The package's R/C++ calculation sources, estimator arguments and RNG calls
are unchanged.

The unpatched test reproduced the missing `startup.Rs` failure. The regression
driver then compared the unpatched clean run with the patched run under four
startup environments: unset, empty, relative, and absolute. For each environment
it compared every native output for three tree shapes, three trait types and
all four evaluator modes, plus the final `.Random.seed`, using `identical()`.
All comparisons passed. Deliberate C++ compilation failures also verified
environment restoration and unchanged RNG state. See
[the scalar results](startup-regression.csv) and
[the regression driver](../../test/test_hybrid_check_startup.R).

The complete local `R CMD check --no-manual` passed with **Status: OK**, no
errors, warnings or notes. Its 13,679 testthat assertions passed, with one
existing cross-install contract skip because the two historical libraries
were not supplied. The check used an independent ASCII checkout containing
the development bridge, so the four-mode test was executed rather than
skipped. See [the final check log](00check.log).
The existing GitHub Actions jobs provide the Linux, macOS, Windows, R-devel,
PDF-manual and OpenMP-disabled checks after the fix is pushed.

The first remote run confirmed the startup repair on Windows, Linux and the
OpenMP-disabled job. It also exposed two independent test assumptions:

- Some ARD fixtures converge on their initial attempt on other BLAS/LAPACK
  configurations. Fresh fits now always receive the independent likelihood
  and probability checks; retry diagnostics are checked when a retry occurs.
  A separate test starts all five captured unfinished fits at locally
  recomputed objectives, guaranteeing recovery-path coverage on every platform.
- R's `log()` and C++ `std::log()` can round Shannon entropy to opposite sides
  of 1 on ARM. The test now checks entropy accuracy at the existing tolerance,
  replays the shared endpoint draws only in an ulp-scale boundary neighbourhood,
  and evaluates the independent R posterior at the exact validated native
  entropy input. It still compares the full chains and final RNG state. The
  arithmetic Gini endpoint tests retain their strict independent expectations.

Both affected test files passed locally after these corrections. These changes
are confined to test assertions and fixtures' documentation; the installed R
and C++ calculations and random-number calls remain unchanged.

To rerun the focused regression from a checkout with a compatible package
installed in the supplied library:

```sh
Rscript --vanilla test/test_hybrid_check_startup.R . /path/to/library
```

An optional third argument supplies the pre-fix standalone script for an exact
before/after comparison. Without it, the driver checks clean and inherited
startup environments against the current test. Test data, seeds and native
calculation calls are identical in both paths. The published release assets
and tag are preserved; this repair affects the repository's test harness.
