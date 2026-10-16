# Clade-scale phylogenetic signal ------------------------------------------

#' Scan phylogenetic signal across clades
#'
#' Computes K, lambda, binary D or categorical Delta on crowns selected
#' independently of the traits. Describes local signal strength and, optionally,
#' its deviation from the complete-tree estimate.
#' Automatic selection uses maximal, non-overlapping clades within the size
#' bounds. Explicit nodes also permit nested clades for descriptive exploration.
#' Node numbers always refer to the tree supplied to this function, before any
#' internal canonicalization. For prepared contexts they refer to `tree$tree`.
#' Local differences describe the observed pattern of phylogenetic conservatism
#' at the selected scale; they do not identify the evolutionary mechanism.
#' Delta is recommended for local/overall description only. Its explicitly
#' requested heterogeneity tests are retained for experimental reproducibility.
#'
#' @md
#' @param tree A phylo tree or a prepared tree context.
#' @param data Vector, matrix or data frame, named by species. K/lambda require
#'   numeric traits; D accepts two states and Delta two or more states (including
#'   factors/character labels).
#'   Unnamed
#'   data with one row per tree tip are assumed to be in tip order.
#' @param method `"K"` (default), `"lambda"`, `"D"` or `"Delta"` (\eqn{\Delta} accepted).
#' @param nodes Optional unique internal node numbers in the supplied tree.
#'   Explicit nodes override the automatic size bounds; `min_tips` still sets
#'   the minimum number of observed species per trait and clade.
#' @param min_tips Minimum observed species per clade (default 20, at least 3).
#' @param max_tips Maximum total clade size for automatic selection, a finite
#'   integer at least 3. It must be at least `min_tips` only when `nodes = NULL`.
#' @param test Optional reference test of signal variation (default FALSE).
#'   K provides a BM reference; lambda provides an approximate common-local-lambda
#'   reference; binary D's conditional ER extension is experimental.
#'   Delta heterogeneity testing is experimental and not recommended; use
#'   `test = FALSE, overall = TRUE` for Delta description. K/lambda/D calibrate
#'   between-crown dispersion; Delta calibrates local/overall log deviations and
#'   implies `overall = TRUE`. Requires non-overlapping clades and at least two
#'   estimable crowns. It is not a test of evolutionary-rate heterogeneity.
#' @param nsim Number of heterogeneity replicates when `test = TRUE`, at least 2.
#' @param simulation_chunk Maximum K BM replicates evaluated together.
#'   Lambda, D and Delta use one replicate workspace and ignore this setting.
#' @param keep_null Retain the simulated discrepancy values in the result.
#' @param max_total_tips Limit on the sum of selected clade sizes. This guards
#'   against quadratic work when many nested nodes are requested.
#' @param progress Whether to report calculation stages.
#' @param overall Calculate the signal on the complete supplied tree, using all
#'   observed tips for each trait, including tips outside the selected crowns.
#'   Optional to keep large-tree scans inexpensive.
#' @param overall_test Test the complete-tree signal. Implies `overall = TRUE`.
#'   K uses tip randomizations; lambda uses the existing likelihood-ratio test
#'   against lambda = 0; D uses its random/Brownian signal tests; Delta uses tip
#'   randomizations. These differ
#'   from the heterogeneity test.
#' @param lambda_reference Fit the common lambda reference with pooled REML
#'   (default) or ML. Local and overall lambda estimates continue to use ML.
#'   REML reduces nuisance-fitting bias in tested middle-lambda scenarios but
#'   can be conservative near zero. Both reference tests are approximate.
#' @param overall_nsim Number of K/Delta overall randomizations, at least 1.
#'   Defaults to `nsim`; unused for lambda and D.
#' @param d_nsim Random/Brownian replicates for D normalization (default 1000,
#'   at least 2), including when `test = FALSE`. Also used for overall D tests.
#'   Increase independently of `nsim` to reduce normalization Monte Carlo noise.
#' @param delta_control Named list of specialist Delta controls: `model` (ARD
#'   or ER), `entropy` (LSE, SE or GINI), `mcmc_sim`, `thin`, `burn`, `lambda0`
#'   and `proposal_sd`. Defaults match `fast_delta()`. MCMC is not automatically
#'   shortened. `reference_fit = "multistart"` provides deterministic common ARD
#'   fitting. `reference_unavailable = "conservative"` explicitly treats only
#'   single-state reference discrepancies as upper-tail events; default strict
#'   handling retains NA. Both options affect the reference test, not observed
#'   Delta. The optional `reference_fit = "profile"` uses a deterministic
#'   ER log-rate grid and Brent search for the common reference only; it requires
#'   `model = "ER"`. The default remains `"nlminb"`. This experimental alternative
#'   changes Q and reference probabilities, and supports the exact independent
#'   endpoint. Observed local/overall estimators are unchanged. Ignored for other
#'   methods. `reference_test = "nuisance"` instead enables an experimental
#'   finite-sample nuisance-domain bound and requires `test = TRUE` and an
#'   explicit `state_space` (a vector, or a named list with every trait).
#'   In nuisance mode only, `statistic = "between_clades"` uses
#'   SD(log(local Delta)); the default is `"local_global"`.
#'   `nuisance_control$proposal_counts = TRUE` adds a fixed half-weight
#'   count-universal sequential proposal to the unchanged Mk proposal mixture.
#'   Both options are experimental and disabled by default. The count proposal
#'   uses no observed state counts and does not change the common Mk null.
#'   `nuisance_control` sets `gamma` (0.001), `tolerance` (0.01), `max_boxes`
#'   (128), optional `rate_bounds`, fixed `proposal_intensities`
#'   (c(0.01, 0.1, 1, 10, Inf)), `proposal_family` ("auto"),
#'   `boundary_power` (24, allowed 8--30) and `max_bank_bytes` (64 MiB).
#'   `correlation` ("auto", "none", "taylor", "forest" or "transient") retains shared rate
#'   coordinates through finite three-state ARD Taylor models with bounded
#'   remainders. Auto also shares the maximum-forest stationary projection in
#'   mixing wide/unbounded boxes, with certified finite-branch errors. "taylor"
#'   requests nonpoint finite boxes; "forest" requests the mixing projection
#'   subject to size/error safeguards. "transient" shares the full transient
#'   near three-state multiclass boundaries, including shared slow-edge
#'   first-jump terms and a bounded Poisson remainder.
#'   Auto also considers this route when neither earlier model is used.
#'   Both require a small finite bound on the four remaining rates and a small
#'   bank. Auto/transient also use complete-domain probability-mass constraints
#'   when informative; their upper bound may still equal one. Subtree-repeat
#'   messages share indices across rectangles. Fixed verified witnesses may
#'   improve the maximization lower bound.
#'   Nuisance `reference_model$search` also retains raw bounds, final width,
#'   `stop_reason`, `P_one_reason`, `gamma_clipped` and `bank_diagnostics`.
#'   A certified final lower bound of one permits early termination; a large
#'   unresolved envelope is reported separately and never implies homogeneity.
#'   Other boxes retain the existing complete
#'   boundary envelope. The default plugin path is unchanged.
#'   These controls must be specified before inspecting the traits. ER supports
#'   2--8 declared states and its complete nonnegative rate domain. ARD requires
#'   exactly three declared states; its complete unbounded domain is covered
#'   by forest/resolvent intervals with a uniform transition error enclosure.
#'   `proposal_family = "auto"` adds six fixed asymmetric components for ARD;
#'   `"ER"` retains the original proposal. More components can improve rare-state
#'   and fast-pair weights at a cost to ER efficiency. `boundary_power = p`
#'   gives a transition enclosure error of 32/2^p, not a P-value error bound.
#'   Finite ARD bounds are a 6-by-2 lower/upper matrix in native rate order
#'   (B->A, C->A, A->B, C->B, A->C, B->C for state_space = c("A","B","C")).
#'   Their validity only covers that predeclared domain. `reference_fit` and
#'   `reference_unavailable` apply only to the default `reference_test = "plugin"`.
#' @param clade_test Calculate each crown's ordinary signal test (default FALSE).
#'   These probabilities do not test differences between crowns. D reuses the
#'   normalization draws; K/lambda/Delta use their specialist test conventions.
#'   Only estimable crowns meeting `min_tips` receive probabilities.
#' @param clade_nsim Number of K/Delta local randomizations, at least 1; defaults
#'   to `overall_nsim`. D uses `d_nsim`; lambda needs no simulation count.
#'
#' @details For K, dispersion is the unweighted sample standard deviation of
#'   log(K) across estimable clades. For lambda it is the sample standard
#'   deviation of lambda; D uses SD(D), and Delta uses SD(log(Delta)). These
#'   dispersion summaries are descriptive. Clade size, tree shape and trait-specific
#'   missingness are preserved in all references.
#'   Larger K/lambda/Delta and smaller D usually describe stronger phylogenetic
#'   clustering within a method. Report the observed difference and its scale
#'   even without a reference test; a reference P evaluates its stated model.
#'   Signal strength alone does not establish stabilizing selection, niche
#'   conservatism as a process, or a lower evolutionary rate.
#'
#'   Lambda uses a parametric bootstrap under a fitted shared lambda Gaussian
#'   model, with a separate profiled mean and diffusion variance for each crown.
#'   Location/scale invariance permits unit-variance, zero-mean simulations.
#'   The common lambda domain is the intersection of the crown-specific domains.
#'   The null is common lambda measured from each retained crown root, not
#'   common lambda on the complete tree. On an ultrametric tree with crown
#'   height h and ancestral depth s, whole-tree lambda in `[0,1]` induces the
#'   contrast-equivalent local parameter lambda*h/(h+(1-lambda)*s). Unequal
#'   crown ages can therefore give different local parameters under a common
#'   whole-tree model. Non-ultrametric crowns need not admit a scalar mapping;
#'   pruning that moves the crown root also changes the comparison scale.
#'   Predefine comparable biological scales and interpret local lambda on them.
#'   The common reference uses pooled REML by default; `lambda_reference = "ML"`
#'   retains the historical fit. Local and overall lambda remain ML estimates.
#'   The test is approximate because the common parameter is estimated.
#'   Historical ML calibration found 8% rejection in a middle-lambda scenario.
#'   Independent REML confirmation gave 5.2% there (1000 datasets, 199 replicates),
#'   1.2% at zero and 4% at one (500 datasets each). Near-zero behaviour remains
#'   conservative. Treat this as a model diagnostic;
#'   exact nominal error control is not established.
#'
#'   D's experimental reference uses a shared-rate symmetric two-state ER Mk model with stationary root
#'   probabilities (1/2, 1/2), conditioned on the observed state count in each
#'   estimable crown. The common rate is fitted using the same conditional
#'   likelihood. Polynomial pruning and downward sampling preserve state counts
#'   without rejection. D normalization is shared between observed and simulated
#'   data; the test conditions on these Monte Carlo means. Minority coding and
#'   a fixed orientation preserve binary-label invariance. The positive rate
#'   domain starts at 1e-8 divided by the mean crown branch length; infinity is
#'   the independent-state endpoint. D testing supports up to 1000 observed tips
#'   per crown; its count probabilities require more work than linear pruning.
#'   The ER reference is a new
#'   model-reference heterogeneity extension, not the original D signal test.
#'   D requires strictly positive, bifurcating crowns. Small or strongly
#'   imbalanced binary crowns can give unstable D normalization; the original D
#'   paper evaluated robustness on datasets with more than 50 species.
#'
#'   Delta estimates each crown with the same ACE, entropy and MCMC estimator
#'   as `fast_delta()`. Local state spaces drop unobserved states; state counts
#'   and frequencies retain original labels. Single-state crowns are unavailable,
#'   not assigned zero. With `overall = TRUE`, signed log(local/overall) ratios
#'   and their equal-crown root mean square are reported. The whole-tree estimate
#'   is calculated directly, not averaged from the crowns. Small crowns and
#'   differing state support may affect comparisons; the original Delta study
#'   recommends at least 20 species. Traditional and rank-normalized/folded
#'   R-hat, bulk/tail ESS, ratio MCSE and relative MCSE are retained for local and
#'   overall Delta. Ratio diagnostics enter reference-replicate warning counts.
#'
#'   Delta heterogeneity tests are outside the recommended workflow. Their
#'   implementation is retained unchanged for explicit experimental calls;
#'   the installed HETERO_NUISANCE_zh.md records the research controls.
#'   An explicitly requested plugin reference uses RMS(log(Delta_local/Delta_overall)),
#'   distinct from between-crown dispersion. One common ER/ARD Mk rate matrix is
#'   fitted on the original tree, integrating missing tips, with a uniform root
#'   mixture. Each complete-tree replicate preserves the observed NA masks and
#'   selected crowns, and refits all local and overall Delta estimators. The
#'   fitted parameters and root convention are in `reference_model`. Transition
#'   probabilities are cached by exact branch length in C++ (up to 64 MiB), and ACE workspaces are reused
#'   by retained tree/state count. Rates are reference-model parameters, not
#'   evidence of rate shifts. Analysis-owned ACE scratch is bounded by
#'   `options(fastphylosig.ace_cache_bytes = 64 * 1024^2)`; zero disables
#'   cross-fit caching. Active fits and allocator overhead are outside this bound.
#'   The expanded calibration includes random/missing and unequal crowns,
#'   2--8 states and asymmetric models. High-state and sparse ARD reference
#'   inference remains exploratory; state scarcity can prevent a reference P.
#'   These checks do not establish general nominal error control.
#'   This approximate reference check has numerical
#'   tests and a bounded three-state whole-tree calibration: two ER scenarios
#'   gave 4.5% and 5% rejection (200 datasets each, 199 replicates, default MCMC).
#'   An asymmetric ARD scenario was conservative with 12.5% unavailable results;
#'   rare-state data can make all reference tests unavailable. These cases do
#'   not establish validity across tree shapes, state counts or missingness.
#'   `reference_fit="multistart"` optionally fits common ARD Q using the existing
#'   prepared ACE objective with deterministic scaled starts and continuation.
#'   Attempts and deviance gains are reported; it does not change observed Delta.
#'   This changes the reference Q, null draws, P and subsequent RNG when selected.
#'   `reference_unavailable="conservative"` optionally counts undefined
#'   single-state reference discrepancies in the upper tail, preserving every
#'   requested draw. It returns an explicitly labelled conservative upper tail;
#'   Delta itself remains undefined for those crowns. Actual calculation errors
#'   still make P unavailable. Defaults retain nlminb and strict failure handling.
#'   `P_heterogeneity_lower`/`P_heterogeneity_upper` describe the unknown reference
#'   tail contribution, not a confidence interval; they coincide when complete.
#'   All size-eligible crowns must be estimable; state-based exclusions are not
#'   silently fixed. Any failed replicate makes the reference probability
#'   unavailable. `nsim_MCMC_warning` reports replicates with MCMC flags.
#'   `nsim_optimizer_warning` counts optimizer warnings; `nsim_fit_warning`
#'   is its compatibility column. `nsim_numerical_recovery` counts ACE analytic
#'   recovery and `nsim_entropy_recovery` counts entropy roundoff recovery,
#'   separately from `nsim_other_warning`, `nsim_calculation_failure` (actual
#'   errors/invalid estimates) and `nsim_unavailable` (all unavailable replicates,
#'   including single-state crowns). Counts flag a replicate if any local or
#'   overall fit is affected; categories can overlap and must not be summed.
#'   Local/overall results retain the corresponding flags and error stage.
#'   Unfinished default ARD fits within Delta receive bounded deterministic
#'   optimizer recovery. When used, optional `nsim_optimizer_retry` and
#'   `nsim_optimizer_recovered` count affected reference replicates, separately
#'   from remaining warnings; observed fits retain optional retry details.
#'   Healthy output structures and the MCMC path are retained. Recovery can
#'   change estimates and subsequent random streams.
#'   a finite estimate does not imply reliable MCMC. Inspect these alongside
#'   observed diagnostics before interpreting the reference probability.
#'
#'   Delta's opt-in nuisance test by default retains that discrepancy and the existing
#'   ACE/entropy/MCMC estimator. It declares the complete biological state
#'   space, including possible states absent at the tips, and a uniform root.
#'   One trait-independent mixture proposal bank is reused for every Q.
#'   Harrison's corrected importance rank includes the observation's weight.
#'   A likelihood-ratio confidence set with noncoverage at most `gamma` and
#'   Berger--Boos maximization handle common-Q uncertainty. Conservative C++
#'   likelihood intervals cover continuous parameter boxes, including ER and
#'   three-state ARD infinite-rate boundaries. Row-sum constraints tighten
#'   pruning intervals. Finite ARD Taylor models share six coordinates across
#'   branches, states and reference likelihoods; signed tail/other contrasts
#'   retain numerator/denominator correlations. Nonlinear and rounding errors
#'   remain enclosed. Wide and unbounded mixing boxes can share their maximum-
#'   forest stationary projection across terminal branches. A coupling
#'   error encloses finite-branch transients; common Bernstein polynomials
#'   certify likelihood and weighted-tail bounds. Near multiclass boundaries,
#'   two common compact coordinates retain the whole transient of a two-edge
#'   generator; shared first-jump terms retain the four slow edges and only
#'   higher jump counts enter a Poisson remainder. Neither a spectral gap nor
#'   finite base rates are required. Whole-domain probability-mass constraints
#'   may further tighten the upper bound. Unsupported nonmixing boxes retain the
#'   previous enclosure. Fixed feasible witnesses improve only the supremum's
#'   lower bound. A budget limit retains the complete-domain upper
#'   bound rather than a candidate maximum. Finite-sample error control is
#'   per trait, unconditional over original inputs and Monte Carlo draws, and
#'   assumes the stated homogeneous model, fixed tree/missingness/crowns,
#'   preplanned domain/proposal and identical observed/reference estimator.
#'   It does not imply conditional control among successful runs or adequate
#'   power. Six-dimensional ARD may still have a large gap or return one;
#'   finite ARD domains must have an independent scientific justification.
#'   Single-state reference discrepancies enter the upper tail with Delta
#'   remaining NA; actual calculation failures retain NA and are never redrawn.
#'   `P_heterogeneity` is the conservative upper bound. `P_nuisance_lower`,
#'   `P_nuisance_upper` and `nuisance_gap` describe numerical maximization of
#'   the fixed-bank test, not confidence limits or binomial Monte Carlo error;
#'   `MCSE_P_heterogeneity` is NA. Reference diagnostics, domain, proposal,
#'   remaining boxes, timing and best-point weight/ESS diagnostics are retained
#'   in `reference_model`. `max_bank_bytes` guards the integer pattern bank,
#'   not total process memory. See the installed `HETERO_NUISANCE_zh.md` guide.
#'
#'   BM simulations have unit diffusion variance. K is invariant to trait
#'   location and scale, so neither ancestral means nor diffusion variance need
#'   to be fitted. For disjoint clades, independent simulations on their crown
#'   subtrees give the same joint distribution of K as one BM simulation on
#'   the complete tree: shared ancestral increments only add clade constants.
#'
#'   For the default plugin tests, the upper-tail probability is
#'   `(1 + sum(H_sim >= H)) / (nsim + 1)`.
#'   Its Monte Carlo error is `sqrt(nsim * q * (1 - q)) / (nsim + 1)`, where
#'   `q = sum(H_sim >= H) / nsim`; this is separate from Delta's MCMC error.
#'   It tests excess dispersion relative to the stated reference, not a
#'   universal null of equal signal under all evolutionary models. Other
#'   homogeneous processes can reject these references. BM intervals are reference
#'   envelopes, not confidence intervals for the observed K. Failed simulations
#'   make the probability unavailable rather than silently changing its null.
#'   Nodes must be chosen without looking at trait values for calibration to
#'   apply. Adjust for multiple traits when interpreting multiple probabilities.
#'
#'   Tips outside selected clades, including small remnants of imbalanced trees,
#'   are reported in `selection$uncovered_tips`. No dense whole-tree distance or
#'   covariance matrix is constructed. Work depends on the sum of clade sizes,
#'   not just the size of the complete tree. Simulations run serially under R's
#'   RNG; use `set.seed()` for reproducibility.
#'   A practical initial scan uses size-bounded disjoint crowns, with
#'   `min_tips = 50` and `max_tips = 200` (use more than 50 observed tips for D).
#'   These are computational starting points, not calibrated statistical minima.
#'   Compare topology-selected bounds before inspecting trait estimates; wider
#'   bounds can improve coverage but change the scale of the question. Reuse
#'   `selection$nodes` across methods on the same source tree. Highly imbalanced
#'   trees may admit only one eligible crown; larger bounds cannot guarantee
#'   multiple disjoint crowns. Biologically specified nodes are useful when
#'   clade identity is the question, and must be non-overlapping for testing.
#'   Predefine the bounds and inspect coverage, retained counts and state support.
#'   Do not select a partition by its trait estimates or smallest P value. If
#'   multiple partitions are explored, report all as sensitivity analyses;
#'   their individual probabilities do not correct for that selection.
#'   See the installed `HETERO_GUIDE_zh.md` for a complete workflow.
#'
#'   Clade traversal, extraction, NA grouping, signal batches and calibration
#'   use C++ kernels. NA pruning retains ape's behavior. K tree factors and
#'   workspaces are reused within a simulation. D shares prepared contrast
#'   factors with `fast_d()`; lambda shares the likelihood evaluator and owned
#'   scratch buffers with `fast_lambda()`. D count probabilities and fitted
#'   transitions are reused across replicates; fixed-count D means need not be
#'   resimulated. D conditional likelihood is specific to this reference and
#'   is not interchangeable with the ordinary Mk likelihood of `fast_ace()`.
#'   Delta reuses `fast_ace()` likelihoods and `fast_delta()` entropy/MCMC kernels.
#'   K local tests reuse the scan's compiled NA groups and `fast_k()`'s bounded-
#'   memory C++ permutation path. The shared ACE tree validator runs in C++;
#'   its ancestral downpass defers independent row normalizations until needed,
#'   preserving repeated floating-point normalization and the Delta RNG stream.
#'   Scans and Delta tests share retained NA-mask workspaces. Local Delta
#'   randomizations shuffle only the crown data and skip unused trace diagnostics;
#'   observed estimates and common-Mk replicates retain their diagnostics.
#'   `simulation_chunk` is an upper
#'   bound: native batches are reduced to fit a 64 MiB budget for dominant
#'   block buffers, with at least one replicate per block. This is not a total
#'   memory limit; inputs, outputs and fixed tree arrays require additional RAM.
#'
#' @return A `fastphylosig_heterogeneity` list containing `clades` (one row per
#'   clade and trait), `heterogeneity` (one row per trait), `selection` (node IDs,
#'   descendant labels, coverage and overlap), and optionally `null` (reference
#'   discrepancy vectors). `clades` also reports optional local signal tests,
#'   state composition for D/Delta and local/overall contrasts when requested;
#'   `heterogeneity` separates descriptive `dispersion`, `local_global_deviation`,
#'   and the reference `test_statistic`/`statistic` with `P_heterogeneity`.
#'   D contrasts use local minus overall, without taking logarithms of D.
#'   `plot()` draws the clade estimates for one trait.
#'   `plot(fit, tree = tree)` colours the selected crowns on the supplied tree.
#'   `overall` contains complete-tree estimates and optional signal tests when
#'   requested. Plots annotate these separately from the dispersion statistic
#'   and its heterogeneity test. Plotting never computes tests or advances the RNG.
#' @examples
#' set.seed(4)
#' tree <- ape::rtree(80)
#' x <- setNames(rnorm(80), tree$tip.label)
#' fit <- fast_hetero(tree, x, min_tips = 10, max_tips = 30,
#'                    test = TRUE, nsim = 19, progress = FALSE)
#' fit$clades
#' fit$heterogeneity
#' @export
fast_hetero <- function(
    tree, data, method = c("K", "lambda", "D", "Delta"), nodes = NULL,
    min_tips = 20L, max_tips = 200L, test = FALSE, nsim = 999L,
    simulation_chunk = 32L, keep_null = FALSE,
    max_total_tips = 1000000L, progress = interactive(),
    overall = FALSE, overall_test = FALSE, overall_nsim = nsim,
    d_nsim = 1000L, delta_control = list(),
    clade_test = FALSE, clade_nsim = overall_nsim,
    lambda_reference = c("REML", "ML")) {
  if (identical(method, "\u0394")) method <- "Delta"
  method <- match.arg(method)
  lambda_reference <- match.arg(lambda_reference)
  if (method == "Delta") delta_control <- .heterogeneity_delta_control(delta_control)
  if (method == "Delta") {
    ace_pool <- .fast_ace_native_pool()
    attr(delta_control,"ace_pool") <- ace_pool
    on.exit(.fast_ace_pool_close(ace_pool),add=TRUE)
  }
  min_tips <- .heterogeneity_integer(min_tips, "min_tips", 3L)
  max_tips <- .heterogeneity_integer(max_tips, "max_tips",
                                    if (is.null(nodes)) min_tips else 3L)
  max_total_tips <- .heterogeneity_integer(max_total_tips, "max_total_tips", 3L)
  simulation_chunk <- .heterogeneity_integer(simulation_chunk, "simulation_chunk", 1L)
  for (nm in c("test", "keep_null", "overall", "overall_test", "clade_test")) {
    z <- get(nm)
    if (!is.logical(z) || length(z) != 1L || is.na(z)) {
      stop(nm, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (test) nsim <- .heterogeneity_integer(nsim, "nsim", 2L)
  if (method == "Delta" && identical(delta_control$reference_test, "nuisance") && !test)
    stop("reference_test='nuisance' requires test=TRUE.", call. = FALSE)
  if (method == "D") d_nsim <- .heterogeneity_integer(d_nsim, "d_nsim", 2L)
  overall <- overall || overall_test || (test && method == "Delta")
  if (overall_test && method %in% c("K", "Delta"))
    overall_nsim <- .heterogeneity_integer(overall_nsim, "overall_nsim", 1L)
  if (clade_test && method %in% c("K", "Delta"))
    clade_nsim <- .heterogeneity_integer(clade_nsim, "clade_nsim", 1L)
  rt <- .runtime_begin(paste0(method, " clades"), progress = progress)
  on.exit(.runtime_on_exit(rt), add = TRUE)
  .runtime_stage(rt, "Checking tree and selecting clades...")
  source_tree <- if (inherits(tree, "fastphylosig_tree")) .as_prepared_tree(tree)$tree else tree
  core <- .heterogeneity_compile(source_tree, method)
  explicit <- !is.null(nodes)
  if (explicit) {
    if (!is.numeric(nodes) || !length(nodes) || any(!is.finite(nodes)) ||
        any(nodes != floor(nodes)) || any(nodes <= core$n_tip) ||
        any(nodes > core$n_total) || anyDuplicated(nodes)) {
      stop("nodes must contain unique internal node numbers in the supplied tree.",
           call. = FALSE)
    }
    nodes <- as.integer(nodes)
  }
  plan <- heterogeneity_plan_cpp(core, source_tree$tip.label, nodes,
                                min_tips, max_tips, max_total_tips)
  nodes <- plan$nodes
  total <- if (plan$total <= .Machine$integer.max) as.integer(plan$total) else plan$total
  overlap <- plan$overlap
  if (test && overlap) {
    stop("Heterogeneity testing requires non-overlapping clades; nested nodes are descriptive only.",
         call. = FALSE)
  }
  if (test && length(nodes) < 2L) {
    stop("Heterogeneity testing requires at least two selected clades.", call. = FALSE)
  }
  table <- .as_named_trait_table(data, source_tree, verbose = FALSE,
                                 input_name = "data", vector_name = "x")
  X <- if (method %in% c("D", "Delta"))
    .heterogeneity_categorical_matrix(table, source_tree$tip.label, method) else
    .as_trait_matrix(table, source_tree, verbose = FALSE)
  if (!ncol(X) || anyNA(colnames(X)) || any(!nzchar(colnames(X))) ||
      anyDuplicated(colnames(X))) {
    stop("data must contain uniquely named, non-empty trait columns.", call. = FALSE)
  }
  if (method == "Delta" && identical(delta_control$reference_test, "nuisance")) {
    X <- .heterogeneity_delta_state_space(X, delta_control)
    for (s in attr(X, "state_labels")) .heterogeneity_nuisance_domain(delta_control, s)
    if (4 * core$n_tip * (nsim + 1) > delta_control$nuisance_control$max_bank_bytes)
      stop("Nuisance reference bank exceeds max_bank_bytes.", call. = FALSE)
  }
  aligned <- matrix(NA_real_, nrow = core$n_tip, ncol = ncol(X),
                    dimnames = list(source_tree$tip.label, colnames(X)))
  shared <- match(source_tree$tip.label, rownames(X))
  present <- !is.na(shared)
  aligned[present, ] <- X[shared[present], , drop = FALSE]
  if (method %in% c("D", "Delta")) attr(aligned, "state_labels") <- attr(X, "state_labels")
  .runtime_stage(rt, sprintf("Calculating %s in %d clades...", method, length(nodes)))
  fitted <- heterogeneity_scan_cpp(plan, aligned, method, min_tips,
    test || method == "D" || (method %in% c("K", "Delta") && clade_test),
    if (method == "Delta") delta_control else NULL, clade_test)
  # Retain the shared actionable vocabulary for uncommon blocked crowns.
  for (i in which(fitted$blocked)) {
    checked <- .inspect_tree_core(plan$trees[[i]], signal = method)
    rows <- (i - 1L) * ncol(X) + seq_len(ncol(X))
    eligible <- rows[fitted$status[rows] != "insufficient_data"]
    fitted$message[eligible] <- .format_actionable_condition(checked, method)
  }
  clade_cores <- fitted$cores
  tip_sets <- lapply(plan$trees, `[[`, "tip.label")
  clades <- data.frame(node = rep(nodes, each = ncol(X)),
    trait = rep(colnames(X), times = length(nodes)), method = method,
    n_tip_total = rep(plan$n_tip, each = ncol(X)), n_species = fitted$n_obs,
    estimate = fitted$estimate, status = fitted$status, message = fitted$message,
    stringsAsFactors = FALSE)
  clades$P_signal <- clades$P_random <- clades$P_Brownian <- NA_real_
  if (clade_test && method == "lambda") clades$P_signal <- fitted$P_signal
  clades$signal_test <- clade_test
  clades$signal_test_method <- if (!clade_test) "none" else switch(method,
    K = "tip randomization", lambda = "likelihood ratio against lambda = 0",
    D = "D random/Brownian contrast references", Delta = "Delta tip randomization")
  clades$signal_nsim_requested <- if (!clade_test || method == "lambda") 0L else
    if (method == "D") d_nsim else clade_nsim
  clades$signal_nsim_successful <- 0L
  if (method %in% c("D", "Delta")) {
    counts <- lapply(seq_len(nrow(clades)), function(r) {
      i <- (r - 1L) %/% ncol(X) + 1L; j <- (r - 1L) %% ncol(X) + 1L
      y <- aligned[plan$tips[[i]], j]; y <- y[!is.na(y)]
      labels <- attr(X, "state_labels")[[j]]
      stats::setNames(tabulate(as.integer(y) + if (method == "D") 1L else 0L,
        nbins = length(labels)), labels)
    })
    clades$n_states <- vapply(counts, function(z) sum(z > 0L), integer(1))
    clades$state_counts <- I(counts)
    clades$state_frequencies <- I(lapply(counts, function(z) if (sum(z)) z / sum(z) else z * NA_real_))
  }
  if (method == "Delta") {
    flag_prototype <- .delta_calculation_flags()
    fields <- c("n_saved", "ESS_alpha", "ESS_beta", "split_Rhat_alpha",
      "split_Rhat_beta", "ESS_Delta", "split_Rhat_Delta", "MCSE_Delta", "MCMC_warning",
      .delta_rank_fields(), names(flag_prototype))
    for (nm in fields) clades[[nm]] <- vapply(seq_len(nrow(clades)), function(r) {
      i <- (r - 1L) %/% ncol(X) + 1L; j <- (r - 1L) %% ncol(X) + 1L
      d <- fitted$delta_details[[i]]
      if (is.null(d) || is.null(d[[j]])) return(if (nm %in% names(flag_prototype)) flag_prototype[[nm]] else if (nm == "MCMC_warning") FALSE else NA_real_)
      d[[j]][[nm]]
    }, if (nm %in% names(flag_prototype)) flag_prototype[[nm]] else if (nm == "MCMC_warning") logical(1) else numeric(1))
    clades <- .delta_attach_optimizer_retry(clades,lapply(seq_len(nrow(clades)),function(r) {
      i <- (r-1L)%/%ncol(X)+1L;j <- (r-1L)%%ncol(X)+1L
      fitted$delta_details[[i]][[j]]
    }))
  }
  plan$trees <- NULL
  clades$BM_mean <- clades$BM_lower <- clades$BM_upper <- NA_real_
  if (method == "D") {
    clades$D_contrast <- clades$D_mean_random <- clades$D_mean_Brownian <- NA_real_
  }
  reference <- switch(method, K = "Brownian motion",
    lambda = paste0("shared lambda Gaussian model (",lambda_reference,")"), D = "shared ER rate, conditional on observed state counts",
    Delta = paste0("common ", delta_control$model, " Mk model, uniform root; full-tree simulation"))
  if (method == "Delta" && identical(delta_control$reference_test, "nuisance"))
    reference <- paste0("declared common ", delta_control$model,
      " Mk domain, uniform root; corrected importance nuisance upper bound")
  if (method == "Delta" && identical(delta_control$reference_fit,"profile"))
    reference <- paste0(reference, "; experimental ER profile fit")
  summary_rows <- vector("list", ncol(X))
  null_values <- if (test && keep_null) stats::setNames(vector("list", ncol(X)), colnames(X)) else NULL
  if (test) .runtime_stage(rt, sprintf("Calibrating %s dispersion with %d %s replicates...", method, nsim, reference))
  for (j in seq_len(ncol(X))) {
    row_ids <- (seq_along(nodes) - 1L) * ncol(X) + j
    calibration <- NULL
    if (method == "D") {
      candidates <- which(clades$status[row_ids] == "ok")
      if (length(candidates)) {
        calibration <- heterogeneity_d_calibrate_cpp(
          lapply(candidates, function(i) clade_cores[[i]][[j]]),
          lapply(candidates, function(i) fitted$data[[i]][[j]]),
          d_nsim, if (test) nsim else 0L, test, keep_null)
        ids <- row_ids[candidates]
        clades$estimate[ids] <- calibration$estimate
        clades$status[ids] <- calibration$status
        bad <- ids[calibration$status != "ok"]
        clades$message[bad] <- "D random/Brownian normalization is numerically degenerate"
        clades$D_contrast[ids] <- calibration$contrast
        clades$D_mean_random[ids] <- calibration$D_mean_random
        clades$D_mean_Brownian[ids] <- calibration$D_mean_Brownian
        if (clade_test) {
          clades$P_random[ids] <- calibration$P_random
          clades$P_Brownian[ids] <- calibration$P_Brownian
          clades$P_signal[ids] <- calibration$P_random
          clades$signal_nsim_successful[ids[calibration$status == "ok"]] <- d_nsim
        }
      }
    }
    valid <- is.finite(clades$estimate[row_ids]) & clades$status[row_ids] == "ok"
    if (method %in% c("K", "Delta")) valid <- valid & clades$estimate[row_ids] > 0
    selected <- which(valid)
    estimates <- clades$estimate[row_ids[selected]]
    values <- if (method %in% c("K", "Delta")) log(estimates) else estimates
    H <- if (length(values) >= 2L) stats::sd(values) else NA_real_
    if (method == "D" && !is.null(calibration)) H <- calibration$dispersion
    summary_row <- data.frame(trait = colnames(X)[[j]], method = method,
      n_clades = length(selected), n_clades_excluded = sum(!valid),
      dispersion = H, scale = switch(method, K = "SD(log(K))", lambda = "SD(lambda)", D = "SD(D)", Delta = "SD(log(Delta))"),
      BM_mean = NA_real_, BM_upper = NA_real_, P_BM = NA_real_,
      MCSE_P_BM = NA_real_, P_heterogeneity = NA_real_, MCSE_P_heterogeneity = NA_real_,
      null_mean = NA_real_, null_upper = NA_real_,
      test_reference = reference, lambda_common = NA_real_, q_common = NA_real_,
      lambda_reference = if (method == "lambda") lambda_reference else NA_character_,
      d_nsim = if (method == "D") d_nsim else NA_integer_,
      nsim_requested = if (test) nsim else 0L,
      nsim_successful = 0L, exceedance_count = NA_integer_,
      status = if (length(selected) >= 2L) "descriptive" else "insufficient_clades",
      stringsAsFactors = FALSE)
    if (test && is.finite(H) && method != "Delta") {
      local_cores <- lapply(selected, function(i) clade_cores[[i]][[j]])
      if (method == "K") {
        calibration <- heterogeneity_calibrate_cpp(local_cores, H, nsim,
                                                simulation_chunk, keep_null)
        clades$BM_mean[row_ids[selected]] <- calibration$clade_mean
        clades$BM_lower[row_ids[selected]] <- calibration$clade_lower
        clades$BM_upper[row_ids[selected]] <- calibration$clade_upper
        for (field in c("nsim_successful", "exceedance_count", "P_BM",
                        "MCSE_P_BM", "BM_mean", "BM_upper")) {
          summary_row[[field]] <- calibration[[field]]
        }
        summary_row$P_heterogeneity <- calibration$P_BM
        summary_row$MCSE_P_heterogeneity <- calibration$MCSE_P_BM
        summary_row$null_mean <- calibration$BM_mean
        summary_row$null_upper <- calibration$BM_upper
      } else {
        if (method == "lambda") calibration <- heterogeneity_lambda_calibrate_cpp(
          local_cores, lapply(selected, function(i) fitted$data[[i]][[j]]), H, nsim, keep_null,
          reml = lambda_reference == "REML")
        for (field in c("nsim_successful", "exceedance_count", "P_heterogeneity",
                        "MCSE_P_heterogeneity", "null_mean", "null_upper"))
          summary_row[[field]] <- calibration[[field]]
        parameter <- if (method == "lambda") "lambda_common" else "q_common"
        summary_row[[parameter]] <- calibration[[parameter]]
      }
      summary_row$status <- if (calibration$nsim_successful == nsim) "ok" else "simulation_failure"
      if (keep_null) null_values[[j]] <- calibration$null
    }
    summary_rows[[j]] <- summary_row
  }
  names(tip_sets) <- as.character(nodes)
  n_covered <- sum(plan$covered)
  selection <- list(nodes = nodes, tips = tip_sets, overlap = overlap,
    mode = if (explicit) "explicit" else "automatic_disjoint",
    total_clade_tips = total, n_covered = n_covered,
    coverage = n_covered / core$n_tip,
    uncovered_tips = source_tree$tip.label[!plan$covered],
    data_rows_removed = setdiff(rownames(X), source_tree$tip.label),
    tree_tips_without_data = source_tree$tip.label[!present])
  overall_result <- NULL
  if (overall) {
    .runtime_stage(rt, paste0("Calculating complete-tree ", method,
      if (overall_test) " and its signal test..." else "..."))
    overall_result <- .heterogeneity_overall(source_tree, aligned, method,
      overall_test, if (method == "D") d_nsim else overall_nsim, simulation_chunk, delta_control)
  }
  summary <- do.call(rbind, summary_rows)
  if (method == "Delta" && (test || clade_test)) {
    extra <- .heterogeneity_delta_tests(source_tree, core, fitted$delta_contexts, aligned,
      nodes, clades, overall_result, delta_control, min_tips,
      test, nsim, keep_null, clade_test, clade_nsim)
    clades <- extra$clades
    if (test) {
      for (nm in names(extra$summary)) summary[[nm]] <- extra$summary[[nm]]
      if (keep_null) null_values <- extra$null
    }
  } else if (clade_test && method == "K") {
    signal <- heterogeneity_k_signal_cpp(clade_cores, fitted$data,
      clades$status == "ok", clade_nsim, simulation_chunk)
    clades$P_signal <- signal$P
    clades$signal_nsim_successful <- signal$nsim_successful
  }
  effects <- .heterogeneity_effects(clades, summary, overall_result)
  clades <- effects$clades; summary <- effects$summary
  if (method == "Delta" && any(clades$MCMC_warning)) warning(
    "Delta MCMC diagnostics require inspection; see ratio/parameter rank_Rhat, ESS_bulk, ESS_tail and relative_MCSE_Delta in clades.", call. = FALSE)
  attr(delta_control,"ace_pool") <- NULL
  out <- list(clades = clades, heterogeneity = summary,
              overall = overall_result,
              selection = selection, null = null_values,
              reference = if (test) reference else NULL,
              delta_control = if (method == "Delta") delta_control else NULL,
              call = match.call())
  class(out) <- "fastphylosig_heterogeneity"
  .runtime_attach(out, .runtime_close(rt, success = TRUE))
}

.heterogeneity_overall <- function(tree, X, method, test, nsim, chunk, delta_control = list()) {
  # Reuse the validated public tree engines and their exact test conventions.
  fit <- tryCatch(if (method == "K") {
    fast_k(tree, X, test = test, nsim = nsim, simulation_chunk = chunk,
      verbose = FALSE, progress = FALSE)
  } else if (method == "lambda") {
    fast_lambda(tree, X, test = test, lambda_profile = FALSE,
      verbose = FALSE, progress = FALSE)
  } else if (method == "Delta") {
    do.call(fast_delta, c(list(tree = tree, x = X, test = test, nsim = nsim,
      return_sim = FALSE, verbose = FALSE, progress = FALSE),
      delta_control[setdiff(names(delta_control), c("reference_fit", "reference_unavailable",
        "reference_test", "state_space", "nuisance_control", "statistic"))]))
  } else {
    fast_d(tree, X, test = test, nsim = nsim,
      chunk_size = min(512L, max(32L, chunk)), return_sim = FALSE, keep_null = FALSE,
      verbose = FALSE, progress = FALSE)
  }, error = function(e) e)
  out <- data.frame(trait = colnames(X), method = method,
    n_species = as.integer(colSums(!is.na(X))), estimate = NA_real_,
    P = NA_real_, test = test,
    test_method = if (!test) "none" else if (method == "K")
      "tip randomization" else if (method == "Delta") "Delta tip randomization" else if (method == "lambda") "likelihood ratio against lambda = 0 (chi-square, df = 1)" else
      "D random association and Brownian threshold references",
    P_random = NA_real_, P_Brownian = NA_real_,
    LR = NA_real_, MCSE_P = NA_real_,
    nsim_requested = if (test && method != "lambda") nsim else 0L,
    nsim_successful = 0L, status = "blocked", message = "",
    stringsAsFactors = FALSE)
  if (inherits(fit, "error")) {
    out$message <- conditionMessage(fit)
    return(out)
  }
  rows <- match(out$trait, fit$trait)
  out$n_species <- fit$n_species[rows]
  out$estimate <- fit$estimate[rows]
  out$status <- fit$status[rows]
  out$message <- fit$message[rows]
  if (test) out$P <- if (method == "D") fit$Pval1_fast[rows] else fit$P_fast[rows]
  if (test && method == "K") {
    out$MCSE_P <- fit$MCSE_P[rows]
    out$nsim_successful <- fit$nsim_successful[rows]
  }
  if (test && method == "lambda") out$LR <- fit$LR_fast[rows]
  if (test && method == "D") {
    out$P_random <- fit$Pval1_fast[rows]
    out$P_Brownian <- fit$Pval0_fast[rows]
    out$MCSE_P <- fit$MCSE_P_random[rows]
    out$nsim_successful <- fit$nsim_successful_random[rows]
  }
  if (method == "Delta") {
    for (nm in c("n_states", "ESS_alpha", "ESS_beta", "split_Rhat_alpha",
                 "split_Rhat_beta", "ESS_Delta", "split_Rhat_Delta",
                 "MCSE_Delta", .delta_rank_fields(), names(.delta_calculation_flags()))) out[[nm]] <- fit[[nm]][rows]
    for(nm in intersect(.delta_optimizer_retry_fields(),names(fit)))out[[nm]] <- fit[[nm]][rows]
    out$MCMC_warning <- .delta_mcmc_warning_summary(out$trait, out$estimate,
      out$ESS_alpha, out$ESS_beta, out$split_Rhat_alpha, out$split_Rhat_beta,
      out$MCSE_Delta, ess_delta = out$ESS_Delta, rhat_delta = out$split_Rhat_Delta,
      rank_diagnostics = as.list(out[.delta_rank_fields()]))$by_trait$warning
    if (test) {
      out$MCSE_P <- fit$P_MCSE[rows]
      out$nsim_successful <- fit$successful_simulations[rows]
    }
  }
  out
}

.heterogeneity_effects <- function(clades, summary, overall) {
  clades$overall_estimate <- clades$difference_overall <- clades$log_ratio_overall <- NA_real_
  summary$local_global_deviation <- NA_real_
  summary$local_global_scale <- ifelse(summary$method %in% c("K", "Delta"),
    "RMS(log(local / overall))", "RMS(local - overall)")
  if (!is.null(overall)) {
    global <- overall$estimate[match(clades$trait, overall$trait)]
    valid <- clades$status == "ok" & is.finite(clades$estimate) & is.finite(global)
    clades$overall_estimate <- global
    clades$difference_overall[valid] <- clades$estimate[valid] - global[valid]
    positive <- valid & clades$estimate > 0 & global > 0
    if (clades$method[[1L]] %in% c("K", "Delta"))
      clades$log_ratio_overall[positive] <- log(clades$estimate[positive]) - log(global[positive])
    for (j in seq_len(nrow(summary))) {
      rows <- clades$trait == summary$trait[j]
      effect <- if (summary$method[j] %in% c("K", "Delta"))
        clades$log_ratio_overall[rows] else clades$difference_overall[rows]
      effect <- effect[is.finite(effect)]
      if (length(effect)) summary$local_global_deviation[j] <- sqrt(mean(effect^2))
    }
  }
  if (!"test_statistic" %in% names(summary)) {
    summary$test_statistic <- summary$scale
    summary$statistic <- summary$dispersion
  }
  list(clades = clades, summary = summary)
}

.heterogeneity_integer <- function(x, name, minimum) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop(name, " must be a finite integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}

.heterogeneity_compile <- function(tree, method) {
  # Fast success path: the compiled parser owns topology/branch geometry.
  # Invalid inputs fall back to the established, detailed R diagnostics.
  labels <- if (is.list(tree)) tree$tip.label else NULL
  nnode <- if (is.list(tree)) tree$Nnode else NULL
  metadata_ok <- inherits(tree, "phylo") && is.character(labels) &&
    length(labels) >= 2L && !anyNA(labels) && all(nzchar(trimws(labels))) &&
    !anyDuplicated(labels) && is.numeric(nnode) && length(nnode) == 1L &&
    is.finite(nnode) && nnode >= 1 && nnode == floor(nnode)
  core <- if (metadata_ok) tryCatch(
    compile_tree_cpp(tree$edge, tree$edge.length, length(labels)),
    error = function(e) NULL) else NULL
  if (!is.null(core) && identical(as.double(nnode), as.double(core$n_node))) return(core)
  checked <- .inspect_tree_core(tree, signal = method)
  checked$issues <- checked$issues[checked$issues$signal == "all", , drop = FALSE]
  stop(.format_actionable_condition(checked, signal = NULL,
    prefix = "fast_hetero() cannot compile this tree"), call. = FALSE)
}

#' @export
print.fastphylosig_heterogeneity <- function(x, ...) {
  cat("Clade-scale phylogenetic signal\n")
  cat(sprintf("%d clades; %.1f%% of tree tips covered; %s\n",
              length(x$selection$nodes), 100 * x$selection$coverage,
              if (x$selection$overlap) "nested clades" else "disjoint clades"))
  if (!is.null(x$overall)) {
    cat("Overall signal:\n")
    columns <- c("trait", "method", "n_species", "estimate", "P", "P_Brownian", "status")
    print(x$overall[, columns], row.names = FALSE, ...)
  }
  cat("Clade variation (reference-test P is separate from signal P):\n")
  columns <- c("trait", "method", "n_clades", "dispersion", "scale",
    "local_global_deviation", "P_heterogeneity", "status")
  print(x$heterogeneity[, intersect(columns, names(x$heterogeneity))], row.names = FALSE, ...)
  method <- unique(x$heterogeneity$method)
  note <- switch(method[[1L]],
    K = "Reference: dispersion under BM; this does not test all possible equal-signal processes.",
    lambda = "Reference: common local crown lambda (approximate); whole-tree lambda uses a different scale.",
    D = "Reference: experimental common ER conditional on crown state counts; lower D means stronger clustering.",
    Delta = if (any(x$heterogeneity$nsim_requested > 0L))
      "Delta heterogeneity test: experimental, not recommended; local/overall descriptions are retained." else
      "Delta: local/overall signal description; heterogeneity testing is outside the recommended workflow.")
  cat(note, "\n")
  invisible(x)
}

#' Plot clade-scale signal estimates
#' @param x A result from `fast_hetero()`.
#' @param trait Trait column name. Required for multiple traits.
#' @param ... Further arguments passed to `graphics::plot()` (clade view) or
#'   `ape::plot.phylo()` (tree view). For example `type = "fan"` selects a fan.
#' @param tree The original phylo tree, or prepared context, used in the scan.
#'   Required for the tree view. Selected nodes and their descendant labels are
#'   checked against the result. Supply the original tree, with its branch lengths.
#' @param view Either `"clades"` or `"tree"`. Defaults to `"tree"` when `tree`
#'   is supplied, otherwise `"clades"`.
#' @param color_scale Tree colour scale: `"auto"` uses log(K/Delta) or raw lambda/D;
#'   `"raw"` and `"log"` select an explicit scale. Log colours require positive
#'   estimates. Legend labels always show original estimates.
#' @param palette Optional vector of at least two colours for the tree view.
#' @param na_col Colour for unselected or unestimable crowns.
#' @param show_legend Draw a colour legend in the tree view.
#' @param node_labels Label selected crown roots with original node numbers.
#' @param show_summary Annotate complete-tree signal/test results (when
#'   calculated), and clade dispersion/heterogeneity test results. Tests not performed or
#'   unavailable are explicitly identified. No tests are run while plotting.
#' @return Invisibly, the clade table for the clade view, or a list with the
#'   plotting tree, source-to-plot node mapping, clade table, owner indices and
#'   colours for the tree view. BM bars are reference envelopes, not confidence
#'   intervals. Tree colours represent crown-clade estimates, not branch-specific
#'   estimates. More deeply nested selected crowns take precedence, including
#'   unavailable estimates. Crown stems retain the enclosing crown colour.
#' @export
plot.fastphylosig_heterogeneity <- function(x, trait = NULL, ..., tree = NULL,
    view = if (is.null(tree)) "clades" else "tree",
    color_scale = c("auto", "raw", "log"), palette = NULL, na_col = "grey80",
    show_legend = TRUE, node_labels = FALSE, show_summary = TRUE) {
  view <- match.arg(view, c("clades", "tree"))
  traits <- unique(x$clades$trait)
  if (is.null(trait)) {
    if (length(traits) != 1L) stop("Specify trait for a multi-trait result.", call. = FALSE)
    trait <- traits[[1L]]
  }
  if (!is.character(trait) || length(trait) != 1L || is.na(trait) || !trait %in% traits) {
    stop("trait must name a trait in the result.", call. = FALSE)
  }
  d <- x$clades[x$clades$trait == trait, , drop = FALSE]
  if (!is.logical(show_summary) || length(show_summary) != 1L || is.na(show_summary))
    stop("show_summary must be TRUE or FALSE.", call. = FALSE)
  summary_text <- if (show_summary) .heterogeneity_summary_text(x, trait) else character()
  if (view == "tree") return(invisible(.heterogeneity_plot_tree(
    x, d, tree, match.arg(color_scale), palette, na_col, show_legend, node_labels,
    summary_text, list(...))))
  d <- d[is.finite(d$estimate), , drop = FALSE]
  if (!nrow(d)) stop("No finite clade estimates for this trait.", call. = FALSE)
  defaults <- list(x = d$estimate, y = seq_len(nrow(d)), pch = 19,
    yaxt = "n", xlab = paste0(d$method[[1L]], " (", trait, ")"),
    ylab = "Clade node in supplied tree", main = "Clade-scale phylogenetic signal",
    xlim = range(c(d$estimate, d$BM_lower, d$BM_upper), finite = TRUE))
  args <- utils::modifyList(defaults, list(...))
  if (length(summary_text)) {
    old_mar <- graphics::par("mar")
    on.exit(graphics::par(mar = old_mar), add = TRUE)
    graphics::par(mar = pmax(old_mar, c(4.1, 4.1, length(summary_text) + 3, 2.1)))
    main_args <- args[intersect(names(args), c("main", "cex.main", "font.main", "col.main"))]
    args$main <- NULL
  }
  do.call(graphics::plot, args)
  if (length(summary_text)) do.call(graphics::title, c(main_args, list(line = 1)))
  graphics::axis(2, at = seq_len(nrow(d)), labels = d$node, las = 1)
  good <- is.finite(d$BM_lower) & is.finite(d$BM_upper)
  if (any(good)) {
    graphics::segments(d$BM_lower[good], which(good), d$BM_upper[good], which(good), col = "grey55")
    graphics::points(d$estimate, seq_len(nrow(d)), pch = args$pch)
    graphics::mtext("Grey bars: 95% BM reference envelopes", side = 1, line = 3, cex = 0.75)
  }
  if (length(summary_text)) graphics::mtext(summary_text, side = 3,
    line = rev(seq_along(summary_text)) + 2, cex = 0.7)
  invisible(d)
}

.heterogeneity_summary_text <- function(x, trait) {
  number <- function(z) if (is.finite(z)) format(signif(z, 4), trim = TRUE) else "unavailable"
  probability <- function(z, draws = 0L) {
    if (!is.finite(z)) return("unavailable")
    # A Monte Carlo zero is zero exceedances, not evidence for P < 1e-4.
    if (z == 0 && is.finite(draws) && draws > 0L) return(paste0("0/", draws))
    base::format.pval(z, digits = 3, eps = 1e-4)
  }
  h <- x$heterogeneity[x$heterogeneity$trait == trait, , drop = FALSE]
  text <- character()
  if (!is.null(x$overall)) {
    g <- x$overall[x$overall$trait == trait, , drop = FALSE]
    if (nrow(g) == 1L) {
      text <- paste0("Overall ", g$method, " = ", number(g$estimate), " (n = ", g$n_species, ")")
      text <- c(text, if (!g$test) "Signal test: not performed" else if (g$method == "D")
        paste0("Signal P: random = ", probability(g$P_random, g$nsim_successful),
          "; Brownian = ", probability(g$P_Brownian, g$nsim_successful)) else
        paste0(if (g$method %in% c("K", "Delta")) "Randomization P = " else "LRT P (lambda = 0) = ",
          probability(g$P, g$nsim_successful)))
    }
  }
  if (nrow(h) == 1L) {
    text <- c(text, paste0(h$scale, " = ", number(h$dispersion), " (", h$n_clades, " crowns)"))
    label <- switch(h$method, K = "BM dispersion", lambda = "Shared-lambda heterogeneity",
      D = "Conditional ER heterogeneity", Delta = "Common-Mk signal-deviation")
    if(h$method=="Delta"&&identical(h$test_statistic,"SD(log(Delta_local))"))
      label <- "Common-Mk between-crown signal dispersion"
    if ("P_heterogeneity_kind" %in% names(h) &&
        identical(h$P_heterogeneity_kind, "single_state_upper_tail"))
      label <- paste0(label, " conservative upper-tail")
    nuisance <- "P_heterogeneity_kind" %in% names(h) &&
      identical(h$P_heterogeneity_kind, "nuisance_upper_bound")
    if (nuisance) label <- paste0(label, " nuisance upper bound")
    simulations <- if ("P_heterogeneity_kind" %in% names(h) &&
        identical(h$P_heterogeneity_kind, "single_state_upper_tail"))
      paste0(h$nsim_successful, "/", h$nsim_requested, " finite simulations") else
      paste0(h$nsim_successful, " simulations")
    text <- c(text, if (!h$nsim_requested) paste0(label, " test: not performed") else
      paste0(label, " P = ", probability(h$P_heterogeneity), " (", simulations, ")"))
    if (h$method == "lambda" && h$nsim_requested > 0L)
      text <- c(text, "Reference: common local crown lambda; approximate P")
    if (h$method == "D" && h$nsim_requested > 0L)
      text <- c(text, "Experimental reference; conditional on crown state counts")
    if (h$method == "Delta") text <- c(text, if (h$nsim_requested > 0L)
      "Experimental heterogeneity test; not recommended" else
      "Descriptive Delta comparison (recommended)")
    if ("local_global_deviation" %in% names(h) && is.finite(h$local_global_deviation))
      text <- c(text, paste0(h$local_global_scale, " = ", number(h$local_global_deviation)))
    if (h$method == "Delta" && h$nsim_requested > 0L)
      text <- c(text, if (nuisance) paste0("Declared Mk/state space; nuisance gap = ",
        number(h$nuisance_gap), if (!h$nuisance_converged) "; domain unresolved" else "") else
        "Reference test of signal deviation; fitted Mk model, approximate P")
    if (h$method == "Delta") {
      if ("P_heterogeneity_lower" %in% names(h) &&
          is.finite(h$P_heterogeneity_lower) && is.finite(h$P_heterogeneity_upper) &&
          h$P_heterogeneity_lower < h$P_heterogeneity_upper)
        text <- c(text, paste0("Reference tail bounds [", probability(h$P_heterogeneity_lower),
          ", ", probability(h$P_heterogeneity_upper), "]; single-state replicates = ", h$nsim_single_state))
      observed_warning <- any(x$clades$MCMC_warning[x$clades$trait == trait]) ||
        (!is.null(x$overall) && any(x$overall$MCMC_warning[x$overall$trait == trait]))
      reference_warning <- if ("nsim_MCMC_warning" %in% names(h)) h$nsim_MCMC_warning else 0L
      if (observed_warning || reference_warning > 0L) text <- c(text,
        paste0("Inspect Delta MCMC diagnostics: observed flags = ", observed_warning,
               "; reference replicates flagged = ", reference_warning))
    }
  }
  text
}

.heterogeneity_plot_tree <- function(x, d, tree, color_scale, palette, na_col,
    show_legend, node_labels, summary_text, dots) {
  if (is.null(tree)) stop("Supply the original tree for view = 'tree'.", call. = FALSE)
  for (nm in c("show_legend", "node_labels")) {
    value <- get(nm)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      stop(nm, " must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(palette)) {
    palette <- grDevices::hcl.colors(64L, "Viridis")
    if (d$method[[1L]] == "D") palette <- rev(palette)
  }
  if (!is.character(palette) || length(palette) < 2L || anyNA(palette))
    stop("palette must contain at least two colours.", call. = FALSE)
  if (!is.character(na_col) || length(na_col) != 1L || is.na(na_col))
    stop("na_col must be one colour.", call. = FALSE)
  grDevices::col2rgb(c(palette, na_col)) # Validate before changing the graphics device.
  source <- if (inherits(tree, "fastphylosig_tree")) .as_prepared_tree(tree)$tree else tree
  core <- .heterogeneity_compile(source, d$method[[1L]])
  nodes <- x$selection$nodes
  rows <- match(nodes, d$node)
  if (anyNA(rows) || anyDuplicated(d$node) || length(x$selection$tips) != length(nodes))
    stop("Inconsistent clade table and selection.", call. = FALSE)
  d <- d[rows, , drop = FALSE]
  sizes <- lengths(x$selection$tips)
  tip_ids <- match(unlist(x$selection$tips, use.names = FALSE), source$tip.label)
  mapping <- heterogeneity_tree_map_cpp(core, nodes, tip_ids,
    as.integer(c(0, cumsum(sizes))))
  plot_tree <- list(edge = mapping$edge, edge.length = mapping$edge_length,
    tip.label = source$tip.label, Nnode = core$n_node)
  if (!is.null(source$node.label))
    plot_tree$node.label <- source$node.label[mapping$source_internal - core$n_tip]
  if (!is.null(source$root.edge)) plot_tree$root.edge <- source$root.edge
  class(plot_tree) <- "phylo"
  attr(plot_tree, "order") <- "cladewise"

  scale <- if (color_scale == "auto") {
    if (d$method[[1L]] %in% c("K", "Delta")) "log" else "raw"
  } else color_scale
  valid <- d$status == "ok" & is.finite(d$estimate)
  if (scale == "log") valid <- valid & d$estimate > 0
  values <- if (scale == "log") log(d$estimate[valid]) else d$estimate[valid]
  clade_colors <- rep(na_col, nrow(d))
  legend_values <- numeric()
  legend_colors <- character()
  if (length(values)) {
    limits <- if (d$method[[1L]] == "lambda" && scale == "raw")
      range(c(0, 1, values)) else range(values)
    if (limits[[1L]] == limits[[2L]]) {
      color_ids <- rep(ceiling(length(palette) / 2), length(values))
      legend_values <- limits[[1L]]
      legend_colors <- palette[color_ids[[1L]]]
    } else {
      color_ids <- 1L + floor((values - limits[[1L]]) /
        diff(limits) * (length(palette) - 1L))
      color_ids <- pmax(1L, pmin(length(palette), color_ids))
      legend_values <- seq(limits[[1L]], limits[[2L]], length.out = 5L)
      legend_colors <- palette[round(seq(1, length(palette), length.out = 5L))]
    }
    clade_colors[valid] <- palette[color_ids]
  }
  if (scale == "log") legend_values <- exp(legend_values)
  colors <- c(na_col, clade_colors)
  edge_colors <- colors[mapping$edge_owner + 1L]
  node_colors <- colors[mapping$node_owner + 1L]
  args <- utils::modifyList(list(x = plot_tree, edge.color = edge_colors,
    tip.color = node_colors[seq_len(core$n_tip)], edge.width = 2,
    show.tip.label = core$n_tip <= 64L,
    main = paste0("Clade ", d$method[[1L]], " on the phylogeny")), dots)
  # Reserve a side panel, without changing caller-supplied plot limits.
  annotations <- show_legend || length(summary_text) > 0L
  # Wrap labels to physical device width before choosing plot limits. A long
  # MCMC note otherwise pushes a right-anchored legend back onto the tree.
  panel_fraction <- 0.45
  plot_width <- graphics::par("pin")[[1L]]
  text_budget <- max(0.5, panel_fraction * plot_width - 0.35)
  summary_display <- unlist(lapply(summary_text, function(label) {
    width <- graphics::strwidth(label, units = "inches", cex = 0.7)
    if (width <= text_budget) return(label)
    base::strwrap(label, width = max(12L,
      floor(nchar(label, type = "width") * text_budget / width)))
  }), use.names = FALSE)
  labels <- c(summary_display,
    if (show_legend) c("Unselected / unavailable", paste0("Clade ", d$method[[1L]])))
  if (length(labels)) panel_fraction <- min(0.7, max(panel_fraction,
    (max(graphics::strwidth(labels, units = "inches", cex = 0.7)) + 0.35) / plot_width))
  height <- max(core$root_distance)
  if (identical(args$use.edge.length, FALSE))
    height <- max(ape::node.depth(plot_tree, method = 2)) - 1
  if (isTRUE(args$root.edge) && !is.null(source$root.edge)) height <- height + source$root.edge
  if (annotations && height > 0 && is.null(args$x.lim) &&
      (is.null(args$type) || args$type == "phylogram") &&
      (is.null(args$direction) || args$direction == "rightwards")) {
    args$x.lim <- c(0, height / (1 - panel_fraction))
  }
  if (annotations && height > 0 && identical(args$type, "fan")) {
    if (is.null(args$x.lim)) args$x.lim <- c(-height * 1.08,
      height * (2.16 / (1 - panel_fraction) - 1.08))
    if (is.null(args$y.lim)) args$y.lim <- c(-height * 1.08, height * 1.08)
  }
  do.call(ape::plot.phylo, args)
  plot_nodes <- mapping$source_to_plot[nodes]
  ape::nodelabels(node = plot_nodes, pch = 21, bg = clade_colors,
    col = "grey20", cex = 0.7)
  if (node_labels) ape::nodelabels(text = nodes, node = plot_nodes,
    frame = "none", adj = c(1.2, -0.4), cex = 0.65)
  if (show_legend) graphics::legend("topright",
    legend = c(format(signif(legend_values, 3), trim = TRUE), "Unselected / unavailable"),
    col = c(legend_colors, na_col), lwd = 3, bty = "n", cex = 0.7,
    title = paste0("Clade ", d$method[[1L]], if (scale == "log") " (log scale)" else ""))
  if (length(summary_display)) graphics::legend("bottomright", legend = summary_display,
    bty = "n", cex = 0.7, text.col = "grey20")
  invisible(list(tree = plot_tree, clades = d, edge_colors = edge_colors,
    node_colors = node_colors, edge_owner = mapping$edge_owner,
    node_owner = mapping$node_owner, source_to_plot = mapping$source_to_plot,
    color_scale = scale, legend_values = legend_values, legend_colors = legend_colors,
    summary_text = summary_text))
}
