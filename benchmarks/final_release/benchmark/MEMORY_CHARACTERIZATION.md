# Memory Characterization

Status: PASS

Single-machine measurements from the exact source commit below. `object.size` is an R object size and is not peak RSS. RDS files were written with `compress = FALSE`; temporary RDS and fresh-process CSV files were removed after validation.

The memory run used the pre-version freeze metadata `0.2.0.9000` at commit
`084e979af137355b4f976faff6e0079237f23f72`. The final `0.2.0` release commit
changes only `DESCRIPTION`'s version field, so the measured production source
is code-equivalent; raw CSV metadata is intentionally unchanged.

- Commit: `084e979af137355b4f976faff6e0079237f23f72`
- Package version: `0.2.0.9000`
- R: `R version 4.6.1 (2026-06-24 ucrt)`
- OS: `Windows 10 x64`
- Compiler: `g++`
- Locale: `C`

| shape | n | status | object.size (bytes) | RDS (bytes) | save (s) | fresh read + first validation (s) |
|---|---:|---|---:|---:|---:|---:|
| balanced | 500 | PASS | 1375896 | 1017237 | 0 |    0.01 |
| balanced | 1000 | PASS | 2714896 | 2005522 |    0.03 |    0.04 |
| balanced | 5000 | PASS | 13434896 | 9955664 |    0.08 |    0.23 |
| balanced | 10000 | PASS | 26834896 | 19893329 |    0.16 |    0.37 |
| balanced | 20000 | PASS | 53654896 | 39878633 |    0.28 |    0.75 |
