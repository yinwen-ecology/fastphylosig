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
