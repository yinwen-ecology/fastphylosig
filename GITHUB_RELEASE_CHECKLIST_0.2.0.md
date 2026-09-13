# fastphylosig 0.2.0 GitHub Release Checklist

Release target: `v0.2.0`

## Frozen Identity

- [x] Tag name is `v0.2.0`.
- [x] The tag resolves to the authoritative source commit
  `3f84775bc96e00820aa39522245cb63175f79a5d`.
- [x] Artifact filename is `fastphylosig_0.2.0.tar.gz`.
- [x] Artifact SHA-256 is
  `d0178b870c6ab77158a434a9fd03c871814c49e360a3d288b4dda218385aeb67`.

## Release Assets

- [ ] Attach `RELEASE_NOTES_0.2.0.md` as the release notes.
- [ ] Attach the authoritative `fastphylosig_0.2.0.tar.gz` artifact.
- [ ] Verify that no historical `fastphylosig_0.1.0.tar.gz` artifact is
  attached to this release.
- [x] Retain `benchmarks/final_release/release/FINAL_RELEASE_MANIFEST.md`.
- [x] Retain `benchmarks/final_release/release/FINAL_RELEASE_PROVENANCE.md`.

## Publication Gate

Before publishing, re-check the tag target, artifact filename, and the
artifact SHA-256 above against the local files. Mark the GitHub release
non-draft only after all three checks pass and the authoritative artifact is
the asset being attached.

- [ ] Final tag, artifact, and SHA-256 recheck completed immediately before
  publication.
- [ ] Release is marked non-draft only after the recheck passes.

The package source and frozen release evidence must remain unchanged during
publication. This checklist does not claim cross-platform qualification or
CRAN readiness.
