# 0.3.5 release materials

This directory records the 0.3.5 source release, mathematical corrections,
reanalysis guidance and validation. Statistical scope and public interfaces
are unchanged from the reviewed patch candidate.

- [English patch notes and reanalysis guidance](RELEASE_NOTES_en.md)
- [中文补丁说明与重算建议](RELEASE_NOTES_zh.md)
- [Final local validation](VALIDATION.md) and [machine-readable results](evidence/VALIDATION.json)
- [Frozen statistical scope](STATISTICAL_SCOPE.json)
- [Frozen package-input manifest](SOURCE_MANIFEST.json)
- [Source archive checksum](SHA256SUMS.txt)

The checked `fastphylosig_0.3.5.tar.gz` archive and validation-materials ZIP are
available from the [0.3.5 release](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.5).
Local copies and preparation logs are retained under `local/` and excluded
from Git. The manifest verifies the package
inputs in this commit. To check them from the repository root, run:

```sh
python docs/release-0.3.5/verify_source.py
```

To verify the exact retained archive, supply its path with `--archive`.
Rebuilding with a different R version or build timestamp can change archive
bytes; the SHA256 identifies the single checked archive, while the manifest
identifies its frozen input sources.

Independent mathematical regressions are part of the ordinary package test
suite, including `test-ace-mathematics.R` and its analytic/enumeration helper.
The standalone cache test in `tests/native/test-ace-exponent-cache.cpp` checks
guarded allocations, failure recovery, reconfiguration, eviction and accounting.
The local Windows check and source assertions are recorded here; GitHub CI
reports its platform results on the release pull request and main branch.

The retained corrected trajectory in `evidence/baseline-captured.rds` is an
input for exact replay, not a new statistical calibration. With an installed
candidate library, replay can be run into a scratch directory containing that
file:

```sh
Rscript --vanilla docs/release-0.3.5/verify_replay.R SCRATCH_DIR CANDIDATE_LIBRARY candidate
```

`verify_contracts.R` compares public results, warning text and final RNG against
a supplied reviewed 0.3.4.9003 installation. `verify_interface.R` compares all
exports and argument defaults between two supplied installations.
`verify_mixed_rng.R` is the Windows dual-DLL check for both RNG kinds and both
normal generators. Its local execution is recorded in the evidence tables.
Historical incremental performance measurements are retained separately in
`evidence/performance-summary.json`; this preparation introduces no numerical
optimization beyond the reviewed development snapshot.
