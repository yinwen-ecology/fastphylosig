# fastphylosig 0.3.1 local patch release check

Status: **NOT_READY for the requested complete release qualification**.
The local patch artifact is built and its functional checks pass. The remaining external check-environment requirement is usable LaTeX. CRAN access was recovered in the authorized network environment; see recovery evidence below. No publication or submission has occurred.

## Artifact and provenance

- Final archive: ../../fastphylosig_0.3.1.tar.gz (project root).
- Version: 0.3.1. SHA256: `760F68BC4DAD1ADE0767292A92805E50A5D4F550265DF0228CAB7DE5D432B4FE`.
- Release notes: ../../RELEASE_NOTES_0.3.1.md; packaged changes are in NEWS.md.
- Source HEAD: 8c600ee097cd7a30e5bf3a58c157822fdedfdef6 plus the local patch.
- Windows 11 x64; R 4.6.1 (ucrt); Rtools 4.5/GCC 14.3.0. Pandoc 3.8.3 was
  found in the existing RStudio installation. No LaTeX executable was available
  to R, and tinytex::tinytex_root() returned an empty root.
- All checks reference this exact final archive. Earlier trial archives and
  hashes were superseded, and their results are not final artifact results.
- An isolated ASCII mirror avoids the local R/CJK-path issue. LC_ALL, LC_CTYPE
  and LANG were set to C only in check subprocesses; no persistent system
  configuration or new software was installed.
- All 109 packaged files correspond to the workspace. R CMD build generates
  Packaged/Author/Maintainer and formats DESCRIPTION; declared DCF fields were
  independently matched. RcppExports.cpp receives standard CRLF-to-LF
  normalization. No other source differences were found (source-manifest.json).

## Fixes and bounded additional review

The five fixes are D type-7 lower interpolation, explicit numeric species
identity, D plot null/P/strict-tail correspondence, near-ultrametric lambda
feasibility, and exact indexing of long mask identities.

Additional confirmed issues were repaired: explicit D p_col selection was
being overridden when aliases differed or a custom table column was selected;
check_tree marked factor tip labels and missing/zero Nnode ready before V2
necessarily rejected them; scalar K dispatcher tests were comparing absent
fields; lambda ncores documentation incorrectly promised OpenMP execution.
Valid representations, strict statistical tails and RNG draw order were retained.
Existing invalid-Nnode diagnostic text was retained to preserve the frozen
inspection contract; 510 targeted legacy-inspection assertions and 33 new
boundary assertions passed. No new numerical/statistical algorithm was added.

The 50K regression covers preparation, subsets, persistence, exact mask
identity, C++/R grouping paths and K without permutations. It does not qualify
all methods or 50K simulation workloads. Historical 0.3.0 performance data,
figures and publication records retain their version meaning. No performance
benchmark or Delta algorithm study was added, and the bounded review is not an
exhaustive audit of every possible input or platform.

## Final check results

| Check | Result | Limits / explanation |
|---|---|---|
| Build and fresh isolated install from exact archive | PASS | clean-install.log; clean-install-smoke.log confirms version and isolated path |
| Independent type-7 header vs stats::quantile(type=7) | PASS: 36 checks | Header compiled directly from the extracted final archive; q(0:3,.25)=.75 |
| Full package tests | PASS 9052; FAIL 0; WARN 0; SKIP 1 | package-tests.Rout |
| Examples | PASS | Both final full initial check and recheck ran examples |
| PDF manual | BLOCKED / FAIL | R explicitly reports pdflatex unavailable; indexed attempt WARNING, non-indexed retry ERROR |
| Standard full initial check | 2 ERRORs, 1 WARNING | Four PSOCK tests could not bind port 11169; the other error/warning is the unavailable LaTeX manual |
| Standard recheck after PSOCK recovery | 1 ERROR, 1 WARNING, 1 NOTE | Full tests/examples pass; manual remains blocked; see compiled-code NOTE below |
| Initial online R CMD check --as-cran | INCOMPLETE, process exit 1 | Initial sandbox network failure; superseded by the complete recovery check below |
| Final complete online R CMD check --as-cran | 1 ERROR, 1 WARNING, 3 NOTEs; exit 1 | Tests PASS 9052 / FAIL 0 / WARN 0 / SKIP 1; examples OK; missing pdflatex causes PDF ERROR/WARNING; NOTEs: New submission, reused-install object metadata, leftover manual .tex |
| Limited offline as-cran static diagnostic | 0 ERROR, 0 WARNING, 1 NOTE | Remote incoming checks disabled; installation/tests/examples/manual intentionally omitted; not a complete as-cran pass |
| Independent Stage 4 bridge | PASS: 27 topology/trait/mode comparisons | All four modes; maximum numerical difference 6.0285110237146e-14; zero tail decision differences |
| Cross-version 0.2.0 vs 0.3.1 K RNG contract | PASS | Existing 0.2.0 archive installed in isolated reference library; guard/labels in a temporary copy of the 0.3.0 harness were updated only to expect 0.3.1 |

The original PSOCK port error was environmental. A process-scoped
R_PARALLEL_PORT chosen with an available local listener resolved all 122
assertions in the affected files. The standard recheck then passed the full
suite using port 13138. No package worker behavior was changed.

The standard recheck reuses the successful installation log and exact isolated
installation to avoid recompiling unchanged code. Its compiled-code NOTE says
object-file symbol provenance is unavailable, and lists linked runtime symbols
_assert/_exit/abort/exit that may come from libraries. The full initial check of
this same artifact had object information and reported compiled code OK. This
NOTE is recorded rather than presented as an algorithm failure or silently
removed. The offline diagnostic NOTE is solely inability to remotely verify
current time. The earlier missing-Pandoc NOTE was resolved by pointing the
subprocess PATH at the already installed Pandoc. No New submission NOTE was produced in the initial incomplete incoming attempt. The final recovery attempt reached incoming feasibility and reported New submission.

The package-internal cross-install test skips because its external harness is
excluded from source tarballs. The standalone packaged bridge also reports
SKIPPED because benchmarks/stage4 is excluded. Both were additionally executed
against the final artifact outside the source tarball; the original project
harness and historical evidence were left untouched. The RNG contract checks
internal permutation streams indirectly via full K-null output and final
.Random.seed; hidden internal permutation rows are not exposed by the API.

## Commands and scope

The build and initial complete standard/online attempts used:

```text
R CMD build source
R CMD INSTALL --no-multiarch --library=<fresh-ship-lib> fastphylosig_0.3.1.tar.gz
R CMD check fastphylosig_0.3.1.tar.gz
R CMD check --as-cran fastphylosig_0.3.1.tar.gz
```

After resolving only the process-scoped PSOCK port, the standard recheck used
successful install-log reuse (all tests, examples and manual attempts enabled):

```text
R CMD check --install=check:<ship-install.log> --library=<ship-lib> fastphylosig_0.3.1.tar.gz
```

The explicitly limited offline diagnostic used
_R_CHECK_CRAN_INCOMING_REMOTE_=false and:

```text
R CMD check --as-cran --no-install --no-tests --no-examples --no-manual fastphylosig_0.3.1.tar.gz
```

To complete the requested qualification, make an existing LaTeX toolchain
available (or separately authorize installation), and
rerun full checks on this unchanged SHA256 archive. Windows qualification does
not establish Linux/macOS, R-devel, another BLAS, or every supported R version.
No commit, push, tag, GitHub release, CRAN submission or paper/figure change was
made. All pre-existing unrelated/untracked work was preserved.


## Bounded recovery evidence (same final archive)

The default sandbox's read-only HTTPS probe failed before TLS with curl error
7. The authorized network probe to the same official CRAN endpoint returned
HTTP 200 without certificate errors. No VPN, certificate validation bypass,
security-setting change or persistent network setting was used. A first
recovery check with a process timeout of 30 seconds timed out fetching
Meta/archive.rds. One bounded R probe at 60 seconds succeeded in 29.06 seconds
(list length 28,086). This justified a final complete online attempt with a
90-second process download timeout; incoming feasibility completed in 34
seconds and reported the ordinary New submission NOTE.

No existing TeX toolchain was found on PATH, common TinyTeX/TeX Live/MiKTeX
locations or the Windows uninstall registry. The installed tinytex R package
reports an empty root. The minimal official installation option is
`tinytex::install_tinytex()` (https://yihui.org/tinytex/). This downloads and
installs new software, and later compilation can download missing TeX
packages; it was not run without separate authorization.

The compiled-code NOTE was independently reproduced on a temporary copy of
the clean installation. The actual tools namespace check_compiled_code
checker returned the linked-runtime-symbol NOTE before adding object metadata.
After copying saved symbols.rds from the full initial check of this same source
archive into only that temporary copy, the same checker returned NULL.
The tested DLL stayed byte-identical to the clean installed DLL throughout,
SHA256 7F7F0B0AEB4540518B77286B3C9DF22A18A92BC1E34201A55E36492A430A693E.
Nine original .o files and their symbol metadata remain available in the
initial full-check tree. This confirms the NOTE is missing object provenance
when reusing a clean installation, rather than a newly introduced forbidden
call. The source archive, actual isolated installation and checker were not
modified. See symbol-recovery.log and recovery-R-network.log.


Final complete recovery command (tests/examples/manual all enabled):

```text
R CMD check --as-cran --install=check:<ship-install.log> --library=<ship-lib> fastphylosig_0.3.1.tar.gz
```

The full online check completed with exit 1: 1 ERROR, 1 WARNING, 3 NOTEs.
Tests: PASS 9052, FAIL 0, WARN 0, SKIP 1; examples and HTML manual OK.
ERROR and WARNING both arise from pdflatex being unavailable. The three NOTEs
are New submission, missing object metadata in the reused clean installation,
and the manual .tex left behind by failed PDF creation. No unrelated error or
warning appeared. See recovery-as-cran-final.log and
recovery-as-cran-final-00check.log. Qualification remains NOT_READY until a
usable LaTeX toolchain permits the PDF checks to complete. No archive rebuild,
code change, software installation or persistent environment change occurred
in this recovery.

## Library delivery status

The ordered three-file Library save was attempted with the current official
Library skill workflow. The default sandbox failed at tool discovery with a
network error. An authorized-network retry reached discovery but reported that
the prepare-upload capability was unavailable in the helper environment. The
failure occurred before any upload sessions or Library writes; no confirmed
Library file IDs exist for any of the three requested files. The native app
schemas are advertised in the agent environment, but the helper environment
cannot resolve the required capability. The current skill explicitly forbids
switching a started helper write to direct actions, so no alternative create
was attempted. Files remain intact locally; web-downloadable Library delivery
is blocked by this environment capability mismatch.

## GitHub patch publication scope

After these checks, the user authorized a GitHub-only v0.3.1 commit, tag and
Release. The historical NOT_READY status above concerns complete release
qualification and remains unchanged. No CRAN submission or LaTeX installation
is authorized. The checked source archive is reused with the same SHA256.
Only GitHub release wording in README and excluded publication documentation
was updated after validation; implementation and regression sources remain
identical. Public repository evidence includes this summary and SHA256SUMS;
detailed raw machine-local logs are retained locally and not published.
