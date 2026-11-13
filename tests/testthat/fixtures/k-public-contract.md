# Public K frozen compatibility fixture

k-public-contract.rds was captured from the independently installed repaired
fastphylosig 0.3.5 library final-verification/library (2026-11-11), using
R 4.6.1. Provenance is also embedded in the fixture. Tests never regenerate it.
Inputs are the four-tip golden tree and two fixed numeric traits, rather than
ape::rtree output. RNG kind is fixed explicitly and restored after capture.

Coverage: vector and batch public fast_k results, explicit four-row permutations
and 31 internal permutations, complete null draws, exact initial/final RNG states,
explicit-permutation RNG non-consumption, all result fields and attributes,
classes, field names, matching/workflow metadata and timing schema. Numerical
results use 1e-12 tolerance; RNG state, top-level names and classes are exact.
Only measured timing durations are zeroed (attributes and scalar timing field).
No diagnostic, status, or scientific output is silently dropped.

This replaces the deleted optional historical cross-install RNG gate with a
mandatory self-contained release test. No old package installation or version-
specific environment variable is required. Existing golden algebra and strategy
selector tests independently cover numerical derivations and algorithm branches.
The four-tip fixture does not replace those larger branch tests. Hidden internal
permutation rows remain unobservable: full null draws plus final RNG state are
stream evidence, not direct inspection of permutation rows.

Intentional removals: obsolete version-pair requirements and historical whole-
object byte identity. Known incorrect historical outputs are not frozen here;
the repaired baseline is used instead. This is K coverage, not a claim of complete
object parity for D, lambda, Delta, or heterogeneity. Before replacing this fixture,
review the proposed field-level changes and independently validate corrected
statistics; do not regenerate it automatically from the candidate under test.
