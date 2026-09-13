# fastphylosig 0.2.0 Final Scientific Validation

## Decision

The exact-tarball focused suite completed 22 cases:

- PASS: 22
- FAIL: 0
- SKIP: 0
- Unexpected warning cases: 0
- Scientific release blocker: NO

## Provenance

- Source commit: `3f84775bc96e00820aa39522245cb63175f79a5d`
- Package version: 0.2.0
- Loaded package: fresh `FINAL_AUDIT_LIB` installed from the authoritative
  tarball
- R: 4.6.1 on Windows 11 x64
- ape: 5.8.1
- Locale: `LC_COLLATE=C`, `LC_CTYPE=C`
- Scientific log SHA-256:
  `4ca0fb5dec1acdcb994f6ce128ddc134bba0daf5069834d5656e2b1654bfbbf6`

## Coverage

- K estimate, inclusive upper-tail P, MCSE, controlled permutations,
  thread-invariant controlled null, large offsets, NA accounting, and
  raw/prepared parity.
- Lambda estimate, log likelihood, LR/P fields, NA matching, input
  immutability, and raw/prepared parity.
- D observed value, fixed random/Brownian nulls, strict-tail P/MCSE accounting,
  rooted polytomy behavior, and raw/prepared parity.
- Delta estimate, MCMC diagnostics, permutation accounting, NA pruning,
  status/warning behavior, and same-seed replay.
- ACE likelihood, rates, SE, ancestral likelihoods, raw/prepared parity, and
  the `CI = FALSE` likelihood-only path.
- Canonicalization V2 representation invariance, branch association,
  fingerprint propagation, mutation rejection, root sensitivity, legacy
  schema rejection, saveRDS/readRDS, and fresh-session reuse.

## Warning disposition

The stable ACE fixture converged without warning for both CI modes. The
deliberately ill-conditioned ACE fixture returned finite results with an
explicit non-convergence warning (`singular convergence (7)`). The suite
requires that warning and records it as expected; it is not suppressed.

Expected Delta low-information/MCMC-diagnostic warnings are also classified
case by case in `status.csv`. The fail-closed summary counts zero unexpected
warning cases.

## Evidence

The authoritative copy is under `scientific/scientific_exact_tarball/` and is
also mirrored in the scientific directory for convenient review:

- `summary.csv`
- `status.csv`
- `scientific_validation.log`
- `ace_warning_provenance.csv`
- `environment.txt`
- fresh-session RDS evidence

Cross-platform scientific qualification was not run and is not implied by
this local result.
