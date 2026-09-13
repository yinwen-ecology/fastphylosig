# fastphylosig 0.2.0 Final Release Provenance

Date: 2026-09-13

## Source and artifact chain

The authoritative source is Git commit
`3f84775bc96e00820aa39522245cb63175f79a5d`. A `git archive` of that commit was
extracted into an ASCII-only staging directory and built with `R CMD build`.
The resulting `fastphylosig_0.2.0.tar.gz` has SHA-256
`d0178b870c6ab77158a434a9fd03c871814c49e360a3d288b4dda218385aeb67`.
The source archive SHA-256 is
`162603bd8aec228b692cd14181a0c942cd8736de2d167797449e13924cd1effd`.

The package was installed from that tarball into a new private
`FINAL_AUDIT_LIB`. Scientific validation and all-export smoke loaded
fastphylosig 0.2.0 from that exact library. The full source test used the exact
extracted Git archive; package checks independently installed and tested the
tarball.

## Final source changes

Compared with `c084ced`, the final source changes are limited to:

- public 0.2.0 release wording in `README.md`, `NEWS.md`, and `docs/index.md`;
- adding `mvtnorm` to `DESCRIPTION` `Suggests`, matching the existing test-only
  namespace call.

`git diff c084ced..3f84775 -- R src tests man NAMESPACE` is empty. There is no
estimator, numerical formula, tolerance, RNG, thread-policy, public-API, or
result-structure change.

## Local environment

- R 4.6.1 (ucrt), x86_64-w64-mingw32
- Windows 11 x64, build 26200
- Rtools45; GCC/G++/GNU Fortran 14.3.0
- C++17; OpenMP compile and link flags enabled
- ape 5.8.1
- Rcpp 1.1.2
- RcppArmadillo 15.4.2.1
- mvtnorm 1.4.2
- phytools 2.5.2
- caper 1.0.4
- testthat 3.3.2
- locale used for formal runs: `LC_COLLATE=C`, `LC_CTYPE=C`

## Check boundary

The authoritative package check is `R CMD check --no-manual --timings`, which
returned `Status: OK`. The separately executed `--as-cran` and ordinary manual
checks reached PDF generation but could not run `pdflatex`; `Rdlatex.log`
states `pdflatex is not available`. These raw logs are preserved and are not
represented as successful checks.

## Historical boundary

The `c084ced` and `146c03c` artifacts are superseded, not deleted. Historical
benchmark rows at commit `084e979` and package metadata 0.2.0.9000 remain valid
pre-version performance provenance because executable production code is
unchanged; they are not relabeled as final artifact runs.

The annotated tag `v0.2.0` was created and verified against the authoritative
source commit above. Its tag object is `a9d8f54eab07b1275ef1af620544b08a22e0cda7`.
Repository-only evidence commits follow without changing the tagged package
source or the artifact identity.

## Deferred qualification

Linux, macOS, R-devel, alternate BLAS/LAPACK, no-OpenMP fallback, and broader
cross-platform deterministic equivalence were not run and remain deferred.
