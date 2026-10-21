# Final local validation: 0.3.5 release

Date: 2026-10-21. Platform: Windows, R 4.6.1. This record applies to the
final release archive, not to the previously published 0.3.4 archive.
GitHub CI reports its platform results on the release pull request and main branch.

## Package check

`R CMD build --no-manual` and `R CMD check --no-manual --timings` completed
successfully. Final status: **OK**, with no errors, warnings or notes.
The package tests report **14,438 PASS, 0 FAIL, 0 WARN, 1 SKIP**.
The skipped cross-install K test requires a separately configured historical
0.2.0 library. Existing development-only benchmark bridges also retain their
documented source-package exclusions; they are not counted as newly executed.

Evidence: [00check.log](evidence/00check.log),
[installation log](evidence/00install.out),
[test output](evidence/testthat.Rout),
[machine-readable validation](evidence/VALIDATION.json).

The Git copies of session/test output trim trailing presentation whitespace;
original raw logs are retained locally. Test content and counts are unchanged.

## Implementation and interface freeze

- All 31 C++ source/header/build files match the reviewed 0.3.4.9003 snapshot.
- Executable syntax in all 21 R files is identical; the roxygen version label is documentation only.
- All 14 exports, argument names and defaults match the reviewed installation.
- The checked Windows DLL's executable `.text` section matches the reviewed DLL exactly: `6f82d96c45938afd81375380c52ba49c393dda5cea71a53e820395c33f374f78`.
- The statistical scope matches the frozen 0.3.4 scope; only version/status metadata were added.

Evidence: [R syntax](evidence/r-syntax.log),
[public interface](evidence/candidate-interface.log),
[source inputs](SOURCE_MANIFEST.json), [statistical scope](STATISTICAL_SCOPE.json).

## Results and random streams

Comparisons use the corrected 0.3.4.9003 implementation. They do not assert
that the corrected ARD probabilities or ACE SEs reproduce released 0.3.4.

| Check | Result |
| --- | --- |
| 25 public workflows: numeric output, diagnostics/warning text and final seed | All identical |
| 3,956 actual optimizer objective values | All identical |
| 8 reconstructed probability matrices | All identical |
| 8 complete MCMC outputs/traces and final seeds | All identical |
| 8 mixed-RNG cases: MT/L'Ecuyer, Inversion/Box-Muller, Delta/hetero | All identical, including future random draws |

Evidence: [public comparisons](evidence/candidate-contracts.csv),
[trajectory replay](evidence/candidate-replay.csv),
[mixed RNG](evidence/mixed-rng.csv).

Independent mathematical tests are included in the package suite. They use
analytic two-state CTMC probabilities, exhaustive latent-state enumeration
and independently calculated likelihood curvature. Legacy reference-package
ARD probabilities and deviance-Hessian SEs are not used as truth for the
corrected quantities.

## Cache checks

The standalone C++ test, compiled with guarded allocations and debug
iterators, passed lifecycle/failure recovery, geometry reconfiguration,
zero-configuration, eviction-allocation and memory-accounting checks.
After warm-up, replacement allocation count is zero; reported and requested
cache storage both equal 27,688 bytes in this local test. Allocator metadata
and total R process memory are outside that count.

Evidence: [lifecycle](evidence/native-lifecycle.log),
[reconfiguration](evidence/native-reconfigure.log),
[replacement allocations](evidence/native-no-allocations.log),
[accounting](evidence/native-accounting.log).

## Archive identity and publication

The single checked source archive is `fastphylosig_0.3.5.tar.gz`:

```text
eb27317dc527eaa5e0bc6d280e309570ba7c8764f5712746f08009a70ec3e4d1
```

It is retained locally with preparation logs and distributed as a release asset.
`verify_source.py` checks the
194 frozen package inputs, scope and optional archive hash. Public release
notes identify affected ARD outputs and both ER/ARD ACE SE corrections.
No additional statistical method or calibration claim is introduced by this
version-preparation step. The release assets include the source package,
validation-materials ZIP and SHA256 checksums; earlier release assets remain available.
