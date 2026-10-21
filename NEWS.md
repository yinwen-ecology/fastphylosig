# fastphylosig 0.3.5

- Correct the ARD ancestral-state downpass to use the child contribution
  `P * L_child` under the row-generator convention. Impossible parent states
  with zero posterior and zero contribution send zero messages. Shared ACE
  fitting, transition matrices, normalization and Delta MCMC kernels are reused.
- ARD ancestral probabilities and downstream Delta estimates, local/overall
  summaries, tests and subsequent RNG may change. The ER downpass retains its
  previous arithmetic; K, lambda, D and ER Delta keep their result/RNG contracts.
- Scale discrete ACE rate standard errors by the information matrix for the
  negative twice log likelihood. This corrects both ER and ARD standard errors;
  fitted rates are unchanged and Delta does not use these standard errors.
- Add independent latent-state enumeration, analytic two-state transitions and
  information-curvature regressions. Reference-package checks remain for valid
  point-likelihood comparisons, not for the corrected ARD probabilities or SEs.
- Cache exact eigenvalue exponent columns for eligible small three-state ARD
  fits. Reset entries on reconfiguration, account for hash-node storage and
  publish only complete columns. The normal R path configures its immutable
  workspace once; its earlier reconfiguration boundary issue was not triggered.
- Reuse Armadillo's existing GEMV dispatch in pre-sized scratch and remove
  repeated vector preparation/copies. Exact exponent evaluation, optimizer
  steps, MCMC, corrected results and RNG contracts are retained. Cache lifecycle,
  memory bounds, FP controls and full-call equivalence have regression checks.
- Retain public interfaces, defaults and statistical scope: K's BM reference,
  lambda's approximate common-local-crown reference, experimental D inference,
  and descriptive Delta with experimental heterogeneity tests outside the
  recommended workflow. Historical ARD calibration does not validate the
  corrected estimator. See docs/release-0.3.5 for impact and reanalysis guidance.

# fastphylosig 0.3.4

- Consolidate `fast_hetero()` clade selection, coverage, trait-wise missing
  values, local/overall summaries and separate signal/reference tests.
  K uses a BM reference; lambda uses an approximate common-local-lambda
  Gaussian reference. D's count-conditioned ER test remains experimental.
  Delta is recommended for local/overall description and MCMC diagnostics;
  its explicit experimental heterogeneity options remain available.
- Explicit nodes override automatic clade-size bounds while `min_tips`
  continues to enforce the observed-species threshold. Automatic selection
  retains ordered size bounds and all node/memory safety checks.
- ARD recovery regressions compare independent likelihoods at matching rate
  vectors, retaining the original tolerance and requiring improvement over
  the same run's unfinished fit.
- Reuse D transition mixtures and Brent's matching best objective value in
  lambda calibration. Fixed-input Windows comparisons preserve statistical
  results and RNG; speed gains depend on workload. K kernels are unchanged.
- Remove two unused private scalar encoding helpers. Refresh documentation,
  source/Rd agreement and executable short examples. No public API is removed.

# fastphylosig 0.3.3.9021 (development)

- Delta heterogeneity tests leave the recommended workflow. Local and overall
  Delta, their diagnostics and ordinary signal tests remain supported. Explicit
  experimental heterogeneity calls retain their existing results and RNG.
- Define K's BM null, lambda's common-local-crown covariance null, and binary
  D's experimental count-conditioned symmetric ER null. Guides distinguish
  observed conservatism patterns from model-reference evidence and mechanisms.
- Document and independently verify lambda's change of scale after crown
  extraction, including unequal crown ages, non-ultrametric geometry and root
  contraction after missing-tip pruning. Estimation and C++ kernels are unchanged.
- Print/tree-plot annotations identify experimental Delta/D inference and the
  local scale of lambda. Recommended examples use descriptive Delta and D;
  all default arguments, numerical results and simulation ordering are retained.

# fastphylosig 0.3.3.9020 (development)

- Explicit nuisance `statistic = "between_clades"` tests SD(log(local Delta)),
  reusing the existing native dispersion kernel. Whole-tree Delta, local/global
  descriptions and all estimator/MCMC calculations are still retained.
- Explicit `nuisance_control$proposal_counts = TRUE` mixes the old proposal
  bank with a fixed count-universal sequential proposal at equal total weights.
  Simulation and outward density evaluation share their actual CDF kernel;
  no observed counts, fitted proposal or lgamma density approximation is used.
  Existing Mk, mixture, likelihood and continuous-domain bound kernels are reused.
- Both additions remain experimental. Count coverage has a polynomial regret
  bound in the independent-composition submodel; it does not cover arbitrary
  slow-mixing ARD correlations. Neither option changes default results or RNG.

# fastphylosig 0.3.3.9019 (development)

- Nuisance diagnostics distinguish a certified final bound of one, all-tail
  banks, unresolved envelopes and gamma clipping. Raw bounds, final width,
  stop reason, state support and duplicate-phenotype tail counts are retained.
- Native support checks select deterministic closed-pair candidates for ARD.
  Existing point likelihoods must verify confidence-set membership and rank;
  candidate values never replace the continuous-domain upper bound. Search
  stops early only when a verified lower bound plus gamma already certifies
  a final probability of one. Reference draws and RNG are unchanged.
- Proposal-density evaluation now uses effective sequential CDF thresholds,
  including rounded overshoots. Legal monotone thresholds retain the previous
  arithmetic; the sampler and default plugin/Delta/ACE/MCMC paths are unchanged.

# fastphylosig 0.3.3.9018 (development)

- Experimental ARD transient envelopes share four slow-edge coordinates with
  the existing two compact base coordinates. A positive first-jump expansion
  retains their noncommuting products and bounds only two-or-more extra jumps
  by a Poisson tail. Existing C++ stochastic-matrix and pruning kernels are reused.
- Subtree-repeat indices are prepared once per nuisance workspace and reused
  across rectangles and boxes. Numeric messages are recalculated at each
  parameter rectangle, using the original multiply/add and scaling order.
- Auto/transient can intersect a complete-domain mass certificate based on
  phenotype normalization, the observation confidence constraint and a
  three-state return-probability floor. Its benefit depends on tail labels;
  it can remain one, and adding gamma may erase a small raw improvement.
  Default plugin estimates and RNG are unchanged. New mass-envelope and
  repeated-message diagnostics are separate from optimizer/MCMC warnings.

# fastphylosig 0.3.3.9017 (development)

- Experimental three-state ARD nuisance bounds retain a common full transient
  near all six multiclass boundary families. Two compact coordinates share
  total base rate and branching ratio across branch lengths; a coupling error
  covers the other four rates without assuming a positive spectral gap or
  bounding the base rates. Existing C++ Taylor pruning and rank certificates
  are reused with size/error safeguards and separate transient diagnostics.
- Fixed, membership-verified ER-direction witnesses can improve the ARD
  auto/transient maximization lower bound, including the complete unbounded
  domain. They do not certify its upper bound. Unsupported boxes retain the
  complete previous enclosure; default plugin estimates and RNG are unchanged.
- Shared-model pruning explicitly retains nonnegative probability messages at
  zero-rate boundaries, avoiding negative subnormal lower endpoints before log.

# fastphylosig 0.3.3.9016 (development)

- Experimental ARD nuisance bounds can share the maximum-forest stationary
  projection across terminal branches in wide and unbounded rate boxes.
  A Dobrushin/reset coupling bound encloses finite-branch transient remainders;
  shared Bernstein polynomials certify likelihood and weighted-tail bounds.
- Count-lumped positive uniformization reuses the interval kernel and caches
  fixed reference coefficients. Default plugin estimates and RNG are unchanged;
  nonmixing boxes keep the previous complete-domain enclosure.

# fastphylosig 0.3.3.9015 (development)

- Experimental ARD nuisance bounds can retain six shared rate coordinates
  through all branch transitions, pruning messages and reference likelihoods.
  First-order Taylor models include rigorous nonlinear/rounding remainders;
  signed tail/other contrasts retain numerator/denominator correlations.
- Finite boxes use an automatic cost gate or explicit Taylor mode. Unbounded
  tails retain the existing forest/Euler enclosure. Uniformization powers are
  shared across branch lengths. Default plugin estimates and RNG are unchanged.

# fastphylosig 0.3.3.9014 (development)

- Experimental Delta nuisance tests now cover unbounded three-state ARD boxes
  through nonnegative forest/resolvent intervals and a uniform Euler error
  enclosure. Continuous search retains all infinite tails without rate caps.
- Stochastic-row constraints tighten pruning bounds. ARD proposals can add
  six fixed, trait-independent rare-state/fast-pair components; the original
  ER-only proposal remains selectable. These changes affect opt-in nuisance
  tests only; default plugin estimates and RNG ordering are retained.
- Boundary coverage does not guarantee useful power or a small numerical gap.

# fastphylosig 0.3.3.9013 (development)

- Delta heterogeneity gains an experimental, opt-in nuisance-domain upper bound
  for explicit ER (2--8 states) and ARD3 state spaces. Corrected importance ranks
  reuse one fixed proposal bank; conservative continuous-rate likelihood boxes
  and an optional finite-sample confidence set handle the common Q uncertainty.
  Search budgets retain an upper bound, never an unverified grid maximum.
- ER covers its infinite-rate endpoint. Unresolved unbounded ARD returns one
  without simulating an unused bank; finite ARD domains must be specified before
  looking at traits. This does not establish useful power for the new method.
- Existing plugin tests, Delta estimates and default RNG ordering are retained.

# fastphylosig 0.3.3.9012 (development)

- ACE's existing real-exponent shortcut also covers positive normal arguments.
  This avoids complex exponential dispatch for tiny positive real eigenvalues
  produced by roundoff. Complex, overflowing, underflowing and nonfinite inputs
  retain their previous evaluation paths; signed zeros are retained. ARD solver
  budgets, recovery guards, candidate selection and RNG boundaries are unchanged.

# fastphylosig 0.3.3.9011 (development)

- Bounded Delta ARD recovery uses the existing stable matrix exponential for
  positive exit rates separated beyond sqrt(machine epsilon). This prevents
  loss of relative spectral accuracy near an absorbing state. Ordinary fits,
  common reference fitting and MCMC retain their paths and RNG boundaries.
- Recovery rechecks its initial objective with the same guarded evaluator;
  original optimizer returns and any initial numerical correction are retained.
  The guard is used consistently for objectives and returned ancestral rows.

# fastphylosig 0.3.3.9010 (development)

- Unfinished default ARD fits in Delta workspaces receive one bounded,
  parameter-scaled continuation, followed only when needed by one bounded
  L-BFGS-B attempt. Both reuse the existing C++ objective and native scratch.
  Successful fits, common multistart, other models and MCMC retain their paths.
- Actual optimizer returns determine candidate selection and convergence.
  Original warnings and retry budgets are retained in optional diagnostics;
  retry/recovery counts are separate from numerical recovery and remaining
  optimizer warnings. Healthy output structures and RNG are retained.

# fastphylosig 0.3.3.9009 (development)

- ARD likelihoods use the existing Armadillo matrix-exponential kernel when
  an eigensystem is defective or ill-conditioned. Ordinary spectral arithmetic
  and MCMC are retained; recovery is deterministic and reported separately
  from optimizer warnings. Default local fitting strategy is retained.
- Multistart prefers an actually converged candidate when its deviance is
  numerically equivalent to an unfinished minimum, with an explicit tolerance
  and recorded selection recovery. Materially better unfinished candidates
  remain unfinished. Successful ordinary selections are retained.
- ACE checks every returned ancestral probability row, repairs only tiny
  boundary roundoff, and rejects material invalidity. Delta validates all rows
  before entropy endpoint jitter, including rows with interior entropy.

# fastphylosig 0.3.3.9008 (development)

- ACE pruning reads child likelihood rows directly into the existing complex
  workspace, avoiding an intermediate transposed real vector for every edge.
  Real exponents in a guarded normal range reuse scalar exp; complex,
  underflowing, overflowing and nonfinite cases retain the original path.
- MCMC acceptance draws reuse R's uniform generator for fixed open-unit bounds,
  retaining endpoint rejection and conditional draw order. Proposals, lgamma,
  accumulations, saved traces, estimators and diagnostic thresholds are retained.
- A predefined independent screen compares generating Q, default nlminb and
  multistart on 400 new inputs, then three cut scales on 300 new inputs. Larger
  crowns improve shuffled-crown detection in the tested 512-tip balanced trees,
  while fitted-Q bias, missing states and diagnostic warnings remain material.
  Defaults, estimators and scientific tests are unchanged; this is finite
  screening evidence, not a claim of general calibration.

# fastphylosig 0.3.3.9007 (development)

- Optional `delta_control$reference_fit="multistart"` improves the search for
  common ARD Q with deterministic starts and continuation of nlminb, reusing one
  prepared ACE objective. Attempts and deviance gains remain visible; the best
  nonconverged fit is never relabelled as converged. Observed Delta is unchanged.
- Reference tail bounds retain all requested draws. Optional
  `reference_unavailable="conservative"` counts single-state reference
  discrepancies in the upper tail without assigning Delta values, dropping or
  replacing draws. Actual computation errors still make P unavailable. Default
  nlminb and strict rules remain; neither option corrects plug-in fitting bias.
- Replicate diagnostic flags use one native pass in place of repeated R closures,
  preserving default estimators, diagnostics and random streams.

# fastphylosig 0.3.3.9006 (development)

- Delta scans and complete-tree estimates separate optimizer warnings, analytic
  ACE recovery, entropy roundoff recovery, other warnings and calculation errors.
  Reference counts flag replicates, with overlapping categories; unavailable
  single-state crowns remain distinct from calculation failures. The retained
  `nsim_fit_warning` column now mirrors `nsim_optimizer_warning`.
- Reference fits reuse a scalar C++ warning evaluator instead of constructing
  per-fit R diagnostic tables, and accumulate replicate counts outside the result
  data frame. Thresholds, reason strings, estimators and random streams are retained.
  Custom R numeric formatting retains the existing warning formatter.
- A predefined independent Delta screen uses 1,600 outer datasets, 199 reference
  simulations and production MCMC settings, with unique branch lengths,
  imbalanced crowns and rare-state ARD. Conditional null rejection is 2.8--6.0%,
  while ARD null unavailability is 5.5--13.0%; shuffled-crown power is 16.2--21.0%.
  These finite, fixed-tree results do not establish general calibration or power.

# fastphylosig 0.3.3.9004 (development)

- Delta rank diagnostics reuse the sorted split traces for full-trace type-7
  quantiles, merging omitted middle observations for odd lengths. Folded ranks
  reuse this order through a monotone merge instead of another full sort.
  Midranks, normal scores, ESS and R-hat arithmetic remain unchanged.
- The existing MCMC kernel counts down to saved iterations instead of using
  per-step modulo. Burn/thin semantics, proposals, acceptance, trace values,
  floating-point sums and random-number call order remain unchanged.

# fastphylosig 0.3.3.9003 (development)

- Delta repairs machine-scale entropy excursions beyond [0,1] only when the
  corresponding probability row is finite, bounded and normalized within a
  small tolerance. It then reuses the existing endpoint jitter and MCMC, warns
  once and reports recovery counts. Healthy entropy arithmetic, estimates and
  RNG remain unchanged; formerly failing calculations can continue. Material
  excursions are errors. The tolerance is capped at 1e-12.
- Entropy calculations cache log(k) once per call while retaining division and
  fuse range checking with endpoint detection, preserving the complete sampler
  outputs and RNG of the initial roundoff-repair implementation.

# fastphylosig 0.3.3.9002 (development)

- ACE recovers failed ER eigensystems with the existing analytic ER transition
  kernel, reusing pruning, ancestral updates and prepared caches. Successful
  eigensystem calculations keep their original arithmetic and RNG. Recovered
  fits warn once and report evaluation counts; previously failed Delta analyses
  can now continue, changing their subsequent RNG and reference probabilities.
  ARD behavior and optimizers are unchanged.

# fastphylosig 0.3.3.9001 (development)

- An explicit experimental Delta common-ER reference fit is available through
  `delta_control = list(model="ER", reference_fit="profile")`. It reuses ACE
  prepared scratch, analytic ER transitions and the lambda Brent implementation,
  compares zero/independent endpoints, and reports fit diagnostics. Default
  nlminb fits, observed Delta estimators and their RNG contract remain unchanged.
  The alternative changes reference Q and probabilities; calibration remains
  exploratory. Infinite ER rates use exact independent-state transitions.

# fastphylosig 0.3.3.9000 (development)

- Delta adds rank-normalized/folded R-hat, bulk/tail ESS and relative MCSE
  for alpha, beta and their ratio. Ratio diagnostics enter automatic flags.
  Historical diagnostics remain available. Reference MCMC and ACE warnings
  have separate counters; diagnostics never resample or change estimates.
- The common lambda reference in fast_hetero now defaults to pooled REML,
  independently checked against dense GLS and held-out null simulations.
  Local/overall ML estimates remain unchanged. lambda_reference="ML" retains
  the historical reference. The new default intentionally changes reference
  fits, null statistics and probabilities; boundary calibration is approximate.
- Delta/ACE analyses reuse native scratch across fits of the same immutable
  tree, state count and model, with a bounded analysis-owned cache and explicit
  cleanup. PSOCK processes use owned scratch. Final ancestral likelihoods use
  the same prepared kernel; modern diagnostic normal scores are cached locally
  and column FFTs reuse preparation with bounded dominant buffers.
  These performance paths preserve the corrected scientific baseline and RNG.
- Expanded whole-tree calibration reports unavailable results, Monte Carlo
  uncertainty and model sensitivity. High-state Delta and sparse ARD reference
  tests remain exploratory; successful numerical checks do not establish
  general nominal error control. See the installed heterogeneity guide.

# fastphylosig 0.3.2.9000 (development)

- Scientific correction: Delta is now the posterior mean of beta/alpha,
  rather than the ratio of posterior means. Its batch-means MCSE and new
  ESS_Delta/split_Rhat_Delta summarize that same ratio trace.
- Corrected lognormal random-walk MH acceptance for independent exponential
  priors on the raw alpha and beta parameters. Invalid proposals are rejected
  once, and chain initialization has explicit RNG sequencing. Delta values,
  tests and RNG trajectories intentionally change from 0.3.1.9000.
- Repeated-length ACE exponent caching now uses its actual 64 MiB byte budget,
  removing the fixed 4096-group cap. Delta caches logs of accepted raw shapes
  so rejected proposals do not repeat logarithms. These performance changes
  preserve the corrected estimator's numerical results and RNG contract.

# fastphylosig 0.3.1.9000 (development)

- ACE fits now combine the exact repeated-length exponent cache with bounded
  parallel batches for uncached branches. A compact immutable worklist avoids
  computing cached exponents again and enables large mixed-length fits; the
  threshold depends on uncached work rather than total tree size. Batch widths
  may be odd, avoiding a separate parallel region for one leftover branch.
  Fully cached and small workloads keep their serial path. Pruning, optimizer
  and RNG order remain unchanged.

- Large ACE likelihood-only fits with distinct branch lengths now evaluate
  independent complex exponent vectors in bounded 2 MiB OpenMP batches, using
  the existing ACE expression and serial pruning/reduction order. The default
  uses at most two threads; `options(fastphylosig.ace_threads = 1L)` selects
  serial execution. Small fits and builds without OpenMP retain the serial path;
  Delta PSOCK workers use serial inner fits. Full objective trajectories, Delta
  outputs and RNG streams are checked against the preceding development build.

- Delta now uses a native V2 canonical graph/UTF-8 byte ordering path while
  retaining source IDs, descriptors and attributes. Rate fits reuse a private
  ACE workspace, transform repeated tip states once per objective, and cache
  exponent vectors for exactly equal branch lengths with bounded storage.
  MCMC caches accepted-state log-gamma values. The optimizer, floating-point
  expression order and random draw order are unchanged; strict result/trace/RNG
  checks accompany the optimization.

- K crown signal tests now call the production C++ permutation kernel on the
  scan's compiled NA-mask groups, removing repeated tree preparation. Automatic
  cuts jump whole selected/undersized traversal intervals and reserve extraction
  buffers. The shared ACE preflight uses the native topology validator; its
  ancestral downpass lazily evaluates repeated row normalizations, including
  floating-point fixed points/two-cycles, instead of scanning the entire matrix
  after every sibling pair. Paired checks retain estimates, tests and RNG streams.

- Local K/lambda signal tests now respect crown eligibility (`min_tips` and
  estimator status). Lambda returns local LRT probabilities in the first scan,
  avoiding a second fit. Delta tests reuse scan NA/ACE contexts and shuffle
  crown-sized vectors without unused permutation trace diagnostics. Shared
  C++ reference summaries now use the same plus-one Monte Carlo error for all
  methods; Mk simulation pools exactly equal branch lengths to reduce cache
  memory and transition preparation. Estimators and RNG order are retained.

- Added multistate Delta to `fast_hetero()`, reusing ACE likelihood workspaces,
  native entropy/MCMC and diagnostics. Local/overall contrasts, categorical
  state counts/frequencies and optional local signal tests accompany tree maps.
- Added an optional new common-Mk reference check of RMS local/overall log
  Delta deviations, with native cached full-tree simulation, uniform-root
  likelihoods integrating NA tips, complete replicate accounting and fitted
  reference metadata. This tests signal deviations, not rate shifts; broad
  statistical calibration remains future work. Descriptive scans remain the
  default. Single-state crowns and failed replicates produce explicit status.
- Tree plots use log Delta and reverse the default D palette so stronger
  signal has a consistent colour direction. D local signal tests reuse its
  normalization simulations without additional RNG draws.

- Shared prepared D contrast factors across the public scalar/batch engines
  and heterogeneity; removed repeated branch adjustments and duplicate
  traversals. Batch columns use one OpenMP region and smaller scratch buffers.
- Reused lambda likelihood scratch buffers in the specialist optimizer and
  heterogeneity bootstrap. D reuses conditional-likelihood/count workspaces,
  fitted transitions and sampling buffers; both bootstraps reuse branch
  standard deviations while retaining their existing RNG order.

- Added lambda heterogeneity calibration with a fitted shared-lambda Gaussian
  reference and binary D scans/calibration with a conditional shared-rate ER
  reference. Both bootstrap engines run in C++ without dense covariances.
- Added `P_heterogeneity`, Monte Carlo error, reference summaries and fitted
  common parameters. K retains its `P_BM` fields as compatible aliases.
- D preserves each crown's observed state count using polynomial pruning and
  downward conditional sampling, reusing random/Brownian normalization
  (`d_nsim`) and count probabilities across draws. Plots distinguish overall
  D signal probabilities from ER heterogeneity.

- Named the heterogeneity API `fast_hetero()`; the earlier development name
  has been removed.
- Added optional complete-tree K/lambda estimates and signal tests through
  `overall`/`overall_test`, reusing the existing native engines. Plots annotate
  the complete-tree signal test separately from the K BM dispersion test.
- Added `plot(h, tree = tree)` to colour crown-clade K/lambda estimates on
  the phylogeny, with source-node labels, logarithmic K colours, grey unavailable
  regions, and deepest-crown precedence for nested selections. Tree annotation
  uses a linear native traversal and checks selected descendant identities.
- Added `fast_hetero()` for continuous-trait K/lambda scans on
  explicit crown clades or automatic disjoint clades with bounded sizes.
  Reports original node identities, per-trait missingness/failures and coverage.
- Added an optional K dispersion calibration against Brownian motion,
  preserving local topology and missingness, using chunked tree simulations
  and the existing compiled K engine without dense covariance matrices.
  BM reference envelopes and the upper-tail probability have explicit model
  interpretations; nested scans remain descriptive.
- Added print/plot methods, dense numerical and simulation oracles, and
  a Chinese guide separating signal heterogeneity from PVR prediction quality.
- Accelerated clade scans with a lightweight checked tree path, grouped
  missing-data kernels and reuse of compiled crown subtrees. BM branch
  propagation now runs in C++, preserving R's random draws and chunk semantics.
- Moved clade indexing, selection, extraction, exact NA-mask grouping and
  batch signal orchestration into C++. Fused BM generation, prepared K
  evaluation and dispersion calibration retain serial R RNG ordering, reuse
  K factors/workspaces and cap dominant block buffers at 64 MiB. NA pruning
  retains ape's established crown/root/unary behavior.
- Removed obsolete native/R heterogeneity entry points, keeping independent
  reference implementations in test helpers. Removed unused K workspaces and
  duplicate shape checks, and
  omitted unused BM fit diagnostics. Complete-data scans defer loading ape.

# fastphylosig 0.3.1

- Aligned read-only tree readiness with V2 requirements for character tip
  labels and a declared positive Nnode. Invalid representations no longer
  appear ready immediately before preparation rejects them.
- Preserved explicit D plot P-value column selection, including differing
  aliases and custom table columns under the automatic null selection.
- Strengthened scalar K dispatcher tests to compare actual fields and
  simulation values, and corrected lambda's serial-execution documentation.
- Preserved explicit numeric species names through both preparation layers;
  positional matching now applies only to missing names or automatic
  data-frame rownames.
- Bound D plot probabilities to the selected null model and preserved fitted
  probabilities and strict contrast-scale tail counts in calibration plots.
- Aligned R and C++ lambda boundaries with ape's scaled-range ultrametric
  tolerance. Approximately ultrametric trees use the exact terminal-edge
  feasibility limit, which can be slightly below phytools' max-height ratio;
  exactly ultrametric and non-ultrametric domains are unchanged.
- Indexed long packed-mask identities by exact equality before using them as
  environment names. Existing short keys and V2 context persistence remain
  supported. The 50K regression covers preparation, subsets, persistence,
  fallback grouping and K without permutations; it does not qualify every
  method or simulation workload at that size.
- Corrected D's type-7 Brownian threshold interpolation between the first
  two order statistics. A first factor level represented by one tip now
  produces one Brownian state per replicate for distinct simulated values,
  instead of an all-zero state vector. Strict threshold and tail rules,
  RNG draw order, and estimator definitions are unchanged.

# fastphylosig 0.3.0

- Replaced the Blomberg's K permutation evaluator's per-permutation tree
  traversals with an adaptive hybrid fused kernel. Where fusion is selected,
  one traversal serves a block of permutations; hybrid evaluation obtains the
  GLS offset from a tree-only dot product instead of a full upward/downward
  pass pair.
- Added an internal strategy selector. It maps the workload (number of tips,
  number of traits, number of permutations, memory budget) onto the evaluator
  variant and the block size. The selector is not part of the public API and
  adds no arguments to any exported function.
- The permutation stream is unchanged. Blocks are a storage unit, not a drawing
  unit: permutations are still generated serially in global order under R's
  RNG. Fixed-seed validation found identical p-values and exceedance counts;
  K values may differ from 0.2.0 at floating-point rounding scale.
- Public API, estimator definitions, the pointwise p-value rule, and the
  threading policy are unchanged. The 0.2.0 evaluator remains in the source
  as an internal reference; the production selector uses hybrid evaluation.

# fastphylosig 0.2.0

- Integrated Canonicalization Contract V2 for deterministic,
  representation-safe tree preparation. Tip IDs and labels are preserved;
  internal nodes and edge rows are canonicalized without using source
  internal-node IDs, source edge order, locale collation, delimiters, or a
  probabilistic hash as ordering truth.
- Replaced the legacy delimiter-based prepared-tree identity with versioned,
  exact structured encodings of `tip.label`, `Nnode`, `edge`, and
  `edge.length`. A package-owned protected record is checked before cached
  evidence or a numerical kernel can be used.
- V2 prepared contexts support `saveRDS()`/`readRDS()` reuse across R sessions.
  Missing, legacy, or unknown context schemas are rejected with an instruction
  to run `prepare_tree()` again.
- Added representation, heterogeneous-branch, locale, mutation, persistence,
  and K/lambda/D/Delta/ACE parity gates. Estimator definitions, RNG streams,
  threading policy, numerical tolerances, and public APIs are unchanged.

# fastphylosig 0.1.0

- Consolidated repeated trait-table validation and C++ quantile code, removed
  obsolete ACE development timing hooks, and reorganized the English and
  Chinese guides around a one-command start and an explicit advanced workflow.
- Declared compatibility with R 4.1.0 and newer. Windows source builds use the
  matching Rtools toolchain; OpenMP remains optional with a serial fallback.
- Added the method-explicit `fast_signal(tree, data, method, ...)` dispatcher
  for K, lambda, D, and Delta, with public `fast_k()` and `fast_lambda()`.
- Kept `fast_d()`, `fast_delta()`, and `fast_ace()` as specialist entry points;
  methods are never inferred from trait values and no public numerical-engine
  selector was added.
- Made `match_tree_data()` the documented matching API and retained
  `match_phylo_data()` only as a deprecated compatibility alias.
- Added shared tree/data preparation, compact workflow metadata, packed NA-mask
  grouping, retained-subtree validation, and bounded prepared-tree caches while
  preserving specialist result classes and statistical contracts.
- Restricted automatic tree normalization to representation-only changes. It
  never invents roots or branch lengths, jitters zero branches, resolves a
  biological polytomy, ultrametricizes a tree, or silently removes tips.
- Added consistent per-trait `status`, `message`, and `note` output for
  non-finite, insufficient, single-state, constant, and undefined traits.
- Added a shared tree-requirement registry, actionable method-specific errors,
  direct ACE topology/branch guards, root-sensitivity fixtures, and explicit
  zero-terminal/zero-internal branch contracts.
- Added controlled K permutation and D random/Brownian null fixtures, lambda
  fixed/optimized pruning tests, unified-versus-specialist equivalence tests,
  large-offset stability tests, and deterministic raw/prepared comparisons.
- Clarified that `fast_delta()` implements a single-tree statistic with local
  ACE and RNG paths. Cross-language parity and the multiple-tree extension are
  outside the package contract.
- Replaced the large-tree NA-subset cache key with the packed mask key, so
  3000/5000-tip analyses no longer hit R's 10,000-byte environment-name limit.
- Set Delta defaults to `mcmc_sim = 10000`, `thin = 10`, and `burn = 100`.
  Added one aggregated retained-sample warning and one aggregated ESS/R-hat/
  MCSE warning per batch; warnings do not rerun or lengthen chains.
- Reused permutation-invariant ACE workspaces inside Delta and skipped unused
  post-fit Hessian/SE work only in the private Delta ancestral-state helper.
  Public `fast_ace()` still computes and returns standard errors.
- Removed generated build artifacts and obsolete development-only benchmark
  records from the working tree; retained only concise release and reference
  notes needed to explain the current package contract.
- Final release-contract patch: ambiguous simultaneous `x`/`X` inputs now fail
  for D and Delta, ACE rejects invalid `method` values clearly, and unused
  Delta `nsim` input is ignored when `test = FALSE`. Prepared-tree fingerprint
  validation now also covers `Nnode`; targeted regressions and the full local
  check were rerun without changing numerical kernels.

# fastphylosig 0.0.4

- Added `check_tree()` for read-only, method-specific readiness checks for K,
  lambda, D, and Delta. Diagnostics distinguish global tree defects from the
  topology and branch-length contracts of individual statistics.
- Added conservative `resolve_tree()`. It standardizes edge order and internal
  node numbering, rechecks the result, records provenance, and refuses unsafe
  repairs such as random polytomy resolution, branch jitter, or tip deletion.
- Added method-level `progress` messages and stable elapsed-time metadata to
  `fast_signal()`, `fast_d()`, `fast_delta()`, and `fast_ace()`.
- Added explicit R-boundary topology guards for D and ACE. D now rejects unary
  nodes or noncanonical roots instead of silently discarding a branch; ACE now
  rejects mixed unary/polytomous topologies that previously could pass a node-
  count-only check.
- Added boundary and numerical-contract tests for ACE transition evaluation;
  production tolerances and public results remain unchanged.

# fastphylosig 0.0.3

- Made the validated tree-pruning implementations the automatic production
  paths for Blomberg's K and Pagel's lambda. Public callers no longer select
  dense, tree, spectral, or validation engines.
- Added streaming, chunked K randomization with Monte Carlo uncertainty and
  optional null retention, without allocating a full tip-by-simulation matrix.
- Added a tree-linear fixed-lambda likelihood and a batched C++ Brent optimizer
  returning lambda, log-likelihood, GLS mean, sigma2, LR, P, and optional
  likelihood profiles without production VCV/eigendecomposition work.
- Added explicit Monte Carlo uncertainty and successful/failed simulation
  counts for both Fritz-Purvis D null hypotheses.
- Retained Delta ESS, split R-hat, MCSE, iteration, and permutation diagnostics.
- Profiled ACE and Delta before optimization. ACE branch-transition arithmetic
  and repeated Delta permutation calculations are the demonstrated hotspots;
  small PSOCK workloads can be slower than serial execution.
- Added production and component benchmarks plus expanded deterministic,
  stochastic, boundary, offset, topology, and cache regression tests.

# fastphylosig 0.0.2

## Tree preparation and cache

- `prepare_tree()` now compiles O(n) topology and traversal data eagerly while
  creating VCV, Cholesky, and lambda spectral payloads only on demand.
- Added separate structural and numerical caches, a byte-budget numerical LRU,
  packed C++ grouping of trait-wise missing-data masks, and `cache_info()`.
- D, Delta, ACE, and tree-engine K workflows no longer request dense VCV data.

## Continuous signal

- Kept dense/Cholesky K as the default compatibility engine.
- Added opt-in `engine = "tree"` for O(n p) observed Blomberg K and
  `engine = "validate"` for explicit tree-versus-dense comparisons.
- Lambda remains on the validated spectral engine with a dense numerical
  fallback; no unvalidated tree-lambda implementation is exposed.

## D and Delta

- Added bounded, fused D null simulation and traversal with configurable
  `chunk_size` and `keep_null`.
- Added Delta chain means/SDs, ESS, split R-hat, covariance-aware MCSE for
  Delta, permutation P MCSE, successful-count reporting, and safe RNG metadata.

## Validation and documentation

- Added deterministic golden contracts, compiled-tree, packed-mask, streaming,
  tree-K, cache, and Delta-diagnostic tests.
- Added a unified environment-controlled benchmark matrix and rewrote the user
  and expert documentation around explicit compatibility boundaries.
