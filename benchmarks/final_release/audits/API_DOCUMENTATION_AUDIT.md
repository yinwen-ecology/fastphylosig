# fastphylosig 0.2.0 API and Documentation Audit

Audit date: 2026-09-13

Audit source: `3f84775bc96e00820aa39522245cb63175f79a5d`

## Result

`API_DOCUMENTATION_AUDIT = PASS`

- `NAMESPACE` contains exactly 13 intended exports: `fast_signal`, `fast_k`,
  `fast_lambda`, `fast_d`, `fast_delta`, `fast_ace`, `prepare_tree`,
  `cache_info`, `match_tree_data`, deprecated `match_phylo_data`,
  `plot_signal`, `check_tree`, and `resolve_tree`.
- No current public selector or export named `fast_sig`, `fast_delta2`,
  `fast_ace2`, `backend`, `engine`, `unsafe`, `trust`, or `skip_check` exists.
- R formals and Rd usage agree for all production functions.
- `match_phylo_data` remains a forwarding compatibility alias and emits the
  documented deprecation warning.
- README, NEWS, docs index, USAGE_zh, and package metadata consistently name
  release 0.2.0. Historical NEWS and benchmark provenance may retain old API
  or 0.2.0.9000 labels only inside clearly historical records.
- The R >= 4.1.0 compatibility contract remains explicit.
- Cross-platform and CRAN qualification are not overstated.

The lambda overview tables name LR testing as a supported capability. The
formal default remains `test = FALSE`, which is correctly shown in code and
help; the tables do not assert that LR fields are unconditional.

`API_CONTRACT_RELEASE_BLOCKER = NO`
