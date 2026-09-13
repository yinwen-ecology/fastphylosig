# fastphylosig 0.2.0 Package Hygiene Audit

Audit date: 2026-09-13

Audit source: `3f84775bc96e00820aa39522245cb63175f79a5d`

## Decision

```text
PACKAGE_HYGIENE_AUDIT = PASS
AUTHORITATIVE_TARBALL_ACCEPTABLE_FOR_0_2_0 = YES
RELEASE_BLOCKER = NO
```

## Authoritative archive

| Field | Value |
|---|---|
| File | `fastphylosig_0.2.0.tar.gz` |
| SHA-256 | `d0178b870c6ab77158a434a9fd03c871814c49e360a3d288b4dda218385aeb67` |
| Size | 231,949 bytes |
| Entries | 106 |
| Embedded version | 0.2.0 |
| R dependency | R >= 4.1.0 |
| Test-only mvtnorm declaration | Present in Suggests |

The tarball contains the intended package payload and no benchmark directory,
paper, website build, check directory, Rout/log evidence, compiled object,
user-machine path, or temporary audit library. `inst/DELTA_REFERENCE.md` is an
intentional installed methodological reference. Stage-labelled files under
`tests/testthat` are regression tests/oracles rather than benchmark output.

`.Rbuildignore` excludes repository evidence, build products, docs/paper,
compiled artifacts, and root/internal review files from the package payload.

## Historical artifacts

The old 0.2.0 tarball at `c084ced` (SHA `c3f05d...`) and intermediate
`146c03c` tarball (SHA `366b57...`) are preserved as superseded evidence. The
0.1.0 tarball also remains historical. None is used for the final 0.2.0 check,
smoke, or scientific result.

Package hygiene is therefore complete for the local release scope.
