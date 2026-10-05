# fastphylosig 0.3.1 patch notes

GitHub patch release. This is not a CRAN submission and does not claim a clean complete CRAN qualification.

Fixes:

- D: restore type-7 interpolation between the first two order statistics.
  Rare first-factor states no longer produce an all-zero Brownian null.
- Preserve explicitly supplied numeric species names in vectors, matrices and
  data frames. Only genuinely missing or automatic names use positional matching.
- D plots: bind the selected random/Brownian null to its fitted probability;
  preserve strict original-contrast tail counts through scaling and rounding.
- Lambda: align R/C++ approximately ultrametric classification with ape's
  scaled-range criterion. Use the exact terminal-edge feasibility bound;
  phytools' max-height ratio can be slightly too high for unequal tip heights.
- Mask caches: index long identities by exact equality before using environment
  names, retaining short-key compatibility and V2 serialization support.
- Bounded release review: honor explicitly selected D probability columns,
  align check_tree readiness with V2's representation requirements, strengthen
  scalar K dispatcher assertions, and describe lambda's serial execution.

Validation scope:

- Independent regression groups use R quantiles, explicit species
  identity, original contrast tails, dense Pagel covariance likelihoods and
  exact mask identity rather than reproducing the defective branches.
- The 50K test covers preparation, subsets, persistence, fallback grouping and
  K without permutations. It does not qualify all methods or 50K simulations.
- Final artifact checks and their exact statuses are recorded in
  docs/release-0.3.1/RELEASE_CHECK.md, including any environment blockers.
- Existing 0.3.0 figures and performance measurements retain their original
  version meaning. No new performance benchmark was run.
- Local Windows checks do not qualify Linux/macOS, R-devel or another BLAS.

Local Windows validation: 9052 PASS, 0 FAIL, 0 WARN, 1 SKIP; examples and isolated-install smoke pass. Independent type-7 oracle: 36 cases; Stage 4 bridge: 27 comparisons; cross-version K RNG contract: pass.

Complete online --as-cran: 1 ERROR, 1 WARNING, 3 NOTEs. PDF manual checks could not complete because pdflatex was unavailable. Notes concern New submission, missing object provenance when reusing the clean installation (independently diagnosed), and the leftover manual TeX file. This is not an error-free complete as-cran pass. No LaTeX was installed.

Archive SHA256: 760F68BC4DAD1ADE0767292A92805E50A5D4F550265DF0228CAB7DE5D432B4FE. The validated archive is reused unchanged. GitHub release wording in README was updated after validation; package implementation and regression files remain identical to the archive. CI checks use --no-manual, so they do not establish PDF-manual completion.
