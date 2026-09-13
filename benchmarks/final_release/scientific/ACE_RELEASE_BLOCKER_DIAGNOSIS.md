# ACE Release-Blocker Diagnosis

## Disposition

- `ACE_RELEASE_BLOCKER = NO`
- `PRODUCTION_CHANGE = NONE`
- `ESTIMATOR/API/RNG/THREAD_POLICY_CHANGE = NONE`
- `CROSS_PLATFORM_QUALIFICATION = NOT_RUN` (per user instruction)

The earlier release-gate warning came from an ill-conditioned audit fixture,
not from a mismatch in the ACE likelihood implementation. The fixture is
kept as an expected-warning regression case so that its explicit
non-convergence remains visible.

## Reproduction

All values below were obtained with the final production code at source commit
`c084ced0c9e44966369f6b0106d185d40909237b`, package `0.2.0`, R 4.6.1,
`ape` 5.8.1, and the current Windows x64 build.

| path | warning | convergence | message | logLik | rate |
|---|---|---:|---|---:|---:|
| stable `fast_ace`, `CI=TRUE` | none | 0 | both X-convergence and relative convergence (5) | -42.4173076049793 | 3.82310196497805 |
| stable `fast_ace`, `CI=FALSE` | none | 0 | both X-convergence and relative convergence (5) | -42.4173076049793 | 3.82310196497805 |
| stable `ape::ace` | none | n/a | n/a | -42.4173076049793 | 3.82310193161385 |
| pathological `fast_ace`, `CI=FALSE` | `transition-rate optimization did not converge: singular convergence (7)` | 1 | singular convergence (7) | -42.8458792736776 | 218.629509083683 |

The stable fixture is the existing seed-31 ER parity scenario from
`tests/testthat/test-fast-signal.R`. It has no warning in either CI mode and
preserves the existing `ape::ace` parity checks. The pathological fixture is
the deterministic 40-tip random tree with branch lengths from `0.01` to `10`
and the cyclic three-state trait used by the earlier final audit.

## Root Cause

The pathological fixture drives the one-parameter ER rate to a very large
value, where the likelihood becomes effectively flat. With the production
default `ip = 0.1`, `nlminb()` returned the explicit PORT termination message
`singular convergence (7)` and `convergence = 1`, while the returned rate,
log likelihood, and ancestral probabilities remained finite.

The observed behavior is conditioning- and starting-value-sensitive:

| start | convergence | fitted rate | deviance |
|---:|---:|---:|---:|
| 1e-6 | 0 | 256.455854627845 | 85.6917585175062 |
| 0.001 | 0 | 252.059942389413 | 85.6917585181196 |
| 0.01 | 1 | 228.823509550455 | 85.6917585295947 |
| 0.1 | 1 | 218.629509083683 | 85.6917585473553 |
| 1 | 0 | 230.938582654148 | 85.6917585274572 |
| 10 | 1 | 212.930941844255 | 85.6917585660465 |
| 100 | 0 | 100.000048465728 | 85.6923127373037 |
| 1000 | 0 | 1000 | 85.6917585160332 |

Thus, different optimizer terminations move the rate substantially while
changing the deviance only at about the `1e-8` scale. Local finite differences
around the default fit changed the objective by only about `2e-11` for a
relative step of `1e-7`, with rapidly diminishing curvature at larger steps.
This is the expected numerical signature of an ill-conditioned, nearly flat
likelihood surface.

The corresponding `ape::ace` fit was also in the same flat region:

- `ape` rate: `218.629767681748`
- fast rate difference: `-2.586e-04`
- log likelihood difference: `1.86e-12`
- `ape` warning: none

An independent R eigen-algebra objective agreed with the C++ objective within
approximately `7.1e-11` over the diagnostic rate grid, including the fitted
region. This does not support a C++ transition/pruning semantic defect.

## CI and NaN Check

The current `fast_ace()` implementation passes `estimate_se = isTRUE(CI)` to
the private fit helper. Consequently, `CI=FALSE` does not execute the
Hessian/SE calculation. The current rerun captured no `NaNs produced` warning
in either CI mode. The only warning in the retained pathological case is the
explicit optimizer non-convergence warning, which remains observable.

No `suppressWarnings()`, test-threshold reduction, optimizer rewrite, or
warning reclassification was used to obtain the passing release summary.

## Release Decision

The pathological fixture is scientifically valid as input but unsuitable as
a zero-warning fixture because it intentionally combines an extreme branch
scale with a cyclic state pattern that produces a nearly unidentifiable ER
rate. It is now marked `warning_expected=TRUE` and must continue to assert
finite output plus the explicit non-convergence message.

The stable seed-31 fixture is the release parity gate. The focused suite now
reports `22 PASS`, `0 FAIL`, `0 SKIP`, and `0 unexpected warning cases`.
Therefore this ACE finding is not a `0.2.0` release blocker.

Any future work on robust starting values or optimizer termination criteria is
post-release numerical work and is outside this final validation task.

## Files and Validation

Updated evidence files:

- `benchmarks/final_release/scientific/run_scientific_validation.R`
- `benchmarks/final_release/scientific/status.csv`
- `benchmarks/final_release/scientific/summary.csv`
- `benchmarks/final_release/scientific/ace_warning_provenance.csv`
- `benchmarks/final_release/scientific/FINAL_SCIENTIFIC_VALIDATION.md`
- this diagnosis file

The scientific validation was rerun against the installed authoritative
`0.2.0` tarball. It completed with exit status 0 and
`scientific_validation=PASS`.
