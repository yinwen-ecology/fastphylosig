# Stage 2D1 Evidence

This directory contains the compact, reviewable evidence for the private K
permutation prototype. `candidate-a/` and `candidate-b/` contain exactness,
timing, provenance, raw-route, and batch guards. `calibration/` records the
production-versus-private-oracle timing calibration. `operation-audit/`
contains the independently collected allocation, copy, and workspace counts.

The final Candidate B correctness run produced a 35 MB ordered-null CSV with
321,700 replicate rows. The repository retains its aggregate exactness result
and SHA-256 in `ORDERED_NULLS_MANIFEST.csv`; the raw CSV was intentionally not
committed. All 321,700 `null_exact` values were true, all exceedance decisions
matched, and the maximum absolute null-K difference was zero.

Formal source commit:
`baea4cf3306b2e9c314d446782dfbe0b98263db6`.

No evidence in this directory is a production binding or authorizes automatic
integration.
