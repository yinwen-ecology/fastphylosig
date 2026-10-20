\name{fast_delta}
\alias{fast_delta}
\title{Fast Delta statistic for categorical traits}
\usage{
fast_delta(tree, x = NULL, test = FALSE, nsim = 1000, se = NULL,
  mcmc_sim = 10000, thin = 10, burn = 100, lambda0 = 0.1,
  proposal_sd = 0.5, entropy = c("LSE", "SE", "GINI"),
  model = "ARD", permutations = NULL, return_sim = test,
  verbose = TRUE, ncores = 1, X = NULL, prepared = NULL,
  progress = interactive())
}
\arguments{
\item{tree}{A conventionally rooted, fully binary \code{"phylo"} object with
finite, strictly positive branch lengths, or a compatible context returned by
\code{prepare_tree()}. The method does not require ultrametricity.}
\item{x, X}{Named categorical trait vector or matrix/data.frame. Matrix and
data.frame row names must be species names. \code{X} is a compatibility alias
for \code{x}; supply only one. Species identifiers are matched by name and
reordered to tree-tip order.}
\item{test}{Whether to compute a permutation P value.}
\item{nsim}{Number of trait permutations when \code{test = TRUE}.}
\item{se}{Unused; retained so unsupported measurement-error input fails
clearly. Use \code{proposal_sd} for the MCMC proposal.}
\item{mcmc_sim}{Number of MCMC iterations per chain and Delta calculation.}
\item{thin}{MCMC thinning interval.}
\item{burn}{First MCMC iteration eligible to be saved.}
\item{lambda0}{Rate of the exponential priors on alpha and beta.}
\item{proposal_sd}{Standard deviation of the log-normal random-walk proposal.}
\item{entropy}{Entropy transformation: \code{"LSE"}, normalized Shannon
entropy \code{"SE"}, or normalized Gini impurity \code{"GINI"}.}
\item{model}{Discrete transition model used by \code{fast_ace()}: \code{"ER"}
or \code{"ARD"}.}
\item{permutations}{Optional integer permutation matrix with exactly
\code{nsim} rows and one column per retained tip.}
\item{return_sim}{Whether to retain simulated Delta values.}
\item{verbose}{Whether to report species matching.}
\item{progress}{Whether to report method-level progress stages and elapsed
runtime. Defaults to \code{interactive()}. Setting this to \code{FALSE}
suppresses progress and species-matching messages; warnings and errors are
still signalled.}
\item{ncores}{Number of worker processes for the Delta permutation test. Small
jobs can be slower in parallel because process startup is not free. Replaying
stochastic results requires a fixed worker count; changing it repartitions
worker RNG streams and is not draw-for-draw reproducible.}
\item{prepared}{Optional object returned by \code{prepare_tree()}; its V2
integrity is validated before cached evidence is reused.}
}
\description{
Computes the single-tree Delta statistic for categorical traits. Ancestral
state probabilities are estimated with \code{fast_ace()}, node uncertainty is
transformed to entropy, and Delta is the posterior
\code{mean(beta / alpha)} from two chains. The lognormal proposal includes its
Hastings density ratio, targeting the Beta likelihood with independent
exponential priors on raw alpha and beta. Ancestral reconstruction,
entropy, and MCMC calculations use Rcpp/RcppArmadillo. The ape comparison path
is retained internally for regression tests, not exposed as a user engine.

Since version 0.3.5, the shared ACE downpass corrects the ARD
child-to-parent contribution under the row-generator convention. This can
change ARD Delta, local/overall summaries, permutation values and subsequent
RNG relative to 0.3.4, including ordinary descriptive use. ER retains its
previous arithmetic and random stream. The ACE rate-standard-error correction
does not enter Delta, which does not estimate or use those standard errors.

Tiny entropy roundoff excursions are repaired only when the corresponding probability
row is finite, bounded in [0,1] and normalized within a small tolerance. Repaired
endpoints use the existing jitter, emit a warning and are described in the result
note. Healthy calculations retain their original numerical and RNG behavior;
formerly failed calculations can continue with a different subsequent RNG state.
}
\value{
For vector input, a \code{phylo_delta}-compatible object; matrix/data.frame
input returns one row per trait. Results report Delta, alpha/beta chain means
and SDs, ESS, split R-hat, batch-means \code{MCSE_Delta} for the ratio trace, requested and
successful iteration counts, and safe RNG metadata. Tests also report P,
\code{MCSE_P}, failed fits, and successful permutation counts. Inspect the
diagnostics before interpreting short runs. Results expose a
\code{timing} attribute with at least \code{total_elapsed} (seconds); vector
results also provide a \code{timing} list field. Per-trait \code{status} and
\code{message}/\code{note} fields distinguish successful, partial,
invalid-trait, insufficient-data, and failed calculations.

Fewer than two species shared globally by the tree and data is a call-level
error. After matching and trait-wise NA pruning, Delta returns a per-trait
insufficient-data result when a trait retains fewer than two species. A batch with one or more
traits retaining 2--19 species emits one warning and still returns estimates.
The metadata records all per-trait retained/removed counts and sample-size
status (including the fewer-than-two error), plus the affected trait details,
in
\code{attr(result, "analysis_metadata")$delta_retained_species}.

For a finite Delta estimate, non-finite chain diagnostics, ESS below 20, or
traditional split R-hat above 1.1 emit one batch warning. The alpha, beta and
beta/alpha sequences also have rank-normalized/folded \code{rank_Rhat_*},
\code{ESS_bulk_*} and \code{ESS_tail_*}. Screening flags rank R-hat above 1.01,
bulk/tail ESS below 100 per original chain, unavailable diagnostics, or
\code{relative_MCSE_Delta = MCSE_Delta / Delta} above 0.05.
These thresholds are diagnostic heuristics, not a convergence guarantee.
Constant traces have unavailable modern diagnostics. This warning never reruns a
chain or alters requested iterations. Per-trait diagnostic warning details
are retained in
\code{attr(result, "analysis_metadata")$delta_mcmc_diagnostics}.

The public defaults include \code{test = FALSE}, \code{nsim = 1000},
\code{return_sim = test}, \code{mcmc_sim = 10000}, \code{thin = 10},
\code{burn = 100}, and \code{model = "ARD"}. Thus permutation P,
\code{MCSE_P}, and retained null values are absent unless testing is
enabled; null values additionally require \code{return_sim = TRUE}.
\code{MCSE_P} uses the finite successful permutation count. These are
computational defaults rather than a guarantee
of adequate mixing for every data set. Users should increase the chain when
ESS, split R-hat, or Monte Carlo error is inadequate.
}
\details{
Native ACE scratch can be reused across fits within one analysis. The
\code{fastphylosig.ace_cache_bytes} option limits retained storage, with a
default of \code{64 * 1024^2} bytes. Zero disables cross-fit caching. Active fits and allocator overhead
are outside this limit. The cache is released on return/error and is absent
from results. PSOCK workers use owned per-fit scratch.
Version 0.3.2.9000 corrects both posterior aggregation and MH acceptance.
Results and RNG trajectories intentionally differ from earlier versions.
\code{ESS_Delta} and \code{split_Rhat_Delta} assess the beta/alpha sequence.
These are traditional split R-hat and bounded-lag ESS, not rank-normalized
or tail diagnostics. The separate rank-normalized fields use pooled mid-ranks,
split chains and Geyer's positive/monotone ESS sequence with FFT autocovariance.
The sampler still uses two chains; longer runs and independent chain checks
are appropriate for final inference. The warning thresholds screen for poor runs;
they do not establish convergence. For final inference, use longer runs,
check all three sequences, and require MCSE small relative to the differences
being interpreted.

The output records \code{optimizer_warning}, \code{numerical_recovery},
\code{numerical_recovery_evaluations}, \code{entropy_recovery},
\code{entropy_recovery_rows}, \code{other_warning}, \code{calculation_failed}
and \code{failure_stage} separately. Recoveries are distinct from optimizer
warnings and calculation errors. Single-state traits remain unavailable and
are not counted as calculation errors.
Defective or ill-conditioned eigensystems use stable matrix exponentials.
Successful default fits retain their optimization path. An unfinished default
ARD fit in the Delta workspace receives at most one parameter-scaled nlminb
continuation and, only if needed, one bounded L-BFGS-B attempt, each limited
to 150 iterations. They reuse the existing C++ objective and do not draw random
numbers. Original optimizer returns, budgets and candidate selection are
retained. Optional \code{optimizer_retry}, \code{optimizer_recovered} and
\code{optimizer_retry_details} fields appear only when a retry occurred;
remaining optimizer warnings are not hidden or relabelled. Common multistart
and public ACE rate fitting are unchanged. ACE permits only machine-scale
probability boundary recovery; Delta validates every probability row before
entropy or endpoint jitter. These abnormal paths can change estimates and
subsequent random streams. Within the corrected estimator, ordinary paths
retain their arithmetic and draw order. Numerical recovery is reported, not treated as a convergence
failure or evidence of model adequacy.

Delta reuses the likelihood-only native ACE fitting workspace. Large uncached
branch-length workloads may use at most two OpenMP threads by default for
independent exponent calculations, combining exact repeated-length caching
with batches of the remaining branches. Set the
\code{fastphylosig.ace_threads} option to \code{1L} to keep those fits serial; see
\code{fast_ace()} for the thread option and bounded storage. Pruning, optimizer
calls on successful fits and MCMC retain their original order. When \code{ncores > 1}, PSOCK
permutation workers use one inner ACE thread to avoid nested parallelism;
the observed fit still follows the option. Changing this ACE option does not
reallocate Delta's random streams. Changing \code{ncores} does.
}
\seealso{\code{fast_signal()}, \code{fast_ace()}}
\examples{
set.seed(4)
tree <- ape::rtree(25)
x <- sample(letters[1:3], ape::Ntip(tree), replace = TRUE)
names(x) <- tree$tip.label
fast_delta(tree, x, mcmc_sim = 10000, thin = 10, burn = 100,
           verbose = FALSE)
}
