> **Retired V1 audit drivers:** The obsolete runners depended on removed V1
> canonicalization internals and were retired during reference/old-version test
> cleanup. Historical evidence and results are retained; retired source is
> recoverable from pre-cleanup Git history. Current package tests assert V2
> and current inspection contracts directly.
>
> The retired driver is audit_canonical_contract.R. The independent public-workflow audit_downstream_propagation.R is retained.
>
> The original README below is historical context. Commands for retired drivers
> are no longer runnable as written.

# Stage 2B2C Canonicalization Contract Safety Audit

This directory contains audit-only scripts and evidence. It does not contain a
production candidate.

- `audit_canonical_contract.R` tests delimiter collisions, representation
  equivalence, canonical signatures, fingerprints, locale ordering, mutation
  detection, stale-cache propagation, and input immutability.
- `audit_downstream_propagation.R` tests raw/prepared K, lambda, D, Delta, ACE,
  status/warning/failure contracts, fixed stochastic inputs, and a deterministic
  locale-sensitive fully binary fixture.
- `results/canonical/` contains canonicalization and fingerprint evidence.
- `results/downstream/` contains public workflow propagation evidence.

The scripts are intended to run from an ASCII staging directory because the
Windows R locale used for the audit cannot safely normalize every non-ASCII
source path.

Final decision: `CANONICALIZATION_CONTRACT_BUG = YES` and
`PUBLIC_RELEASE_BLOCKER_0_2_0 = YES`. Candidate 2 remains rejected, and Stage
2C was not started.
