# Rank diagnostic oracle fixture

`delta-rank-oracle.rds` contains artificial input alpha/beta traces and 450
expected numeric values computed by upstream `posterior` diagnostic functions.
It does not contain upstream executable code. Generation and base-R dependency
utility substitutions are recorded in
`docs/signal-heterogeneity/five-goals-20261013/posterior-diagnostic-oracle.R`.
The archived reference source and BSD 3-Clause license are excluded from the
package source archive; their checksums are recorded in the follow-up audit.
The production C++ FFT/rank implementation is independently implemented.
