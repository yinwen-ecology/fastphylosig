---
layout: default
title: Fast phylogenetic signal analysis in R
description: fastphylosig 0.3.0 measures phylogenetic signal in continuous, binary, and categorical traits.
---

<nav class="page-jumps" aria-label="Page sections">
  <a href="#speed">Measured speed</a>
  <a href="#install">Install</a>
  <a href="#workflow">Workflow</a>
  <a href="#plots">Plot examples</a>
</nav>

<section class="hero" id="top" markdown="1">
<p class="eyebrow">fastphylosig · R package · v0.3.0</p>

# Fast phylogenetic signal analysis in R

<p class="hero-proof"><strong>At 2,000 tips, observed reference/fastphylosig median-time ratios reached 151× for K, 23,591× for lambda, 224× for D, and 33.6× for Delta.</strong> Delta is a call-time comparison only; <a href="#delta-caveat">diagnostic warnings apply</a>.</p>

<p class="hero-lead">Estimate phylogenetic signal in continuous, binary, and categorical traits. Start with one high-level function, or prepare a tree once and call a method directly.</p>

<div class="speed-cards" aria-label="Observed reference-to-fast runtime ratios">
  <article class="speed-card">
    <p class="speed-card__method">Blomberg's K</p>
    <p class="speed-card__value">135–151×</p>
    <p class="speed-card__note">at 2,000 tips</p>
  </article>
  <article class="speed-card">
    <p class="speed-card__method">Pagel's lambda</p>
    <p class="speed-card__value">20,697–23,591×</p>
    <p class="speed-card__note">at 2,000 tips</p>
  </article>
  <article class="speed-card">
    <p class="speed-card__method">Fritz–Purvis D</p>
    <p class="speed-card__value">218–224×</p>
    <p class="speed-card__note">at 2,000 tips</p>
  </article>
  <article class="speed-card speed-card--delta">
    <p class="speed-card__method">Categorical Delta</p>
    <p class="speed-card__value">23.9–33.6×</p>
    <p class="speed-card__note">at 2,000 tips · call time only · <a href="#delta-caveat">diagnostic warning</a></p>
  </article>
</div>

<p class="hero-footnote">Ratios are reference median time ÷ fastphylosig median time. The 2,000-tip ranges span three designed signal scenarios; they are descriptive ranges, not confidence intervals.</p>
</section>

## Matched runtime comparison {#speed}

Four statistics, 50–2,000 tips, three signal scenarios, and 10 paired calls per setting. Tree and data preparation are excluded.

<figure class="benchmark-figure">
  <a class="figure-open" href="{{ '/release-0.3.0/benchmark_reference_four_methods_speedup.png' | relative_url }}" aria-label="Open the full-size four-statistic benchmark figure">
    <img src="{{ '/release-0.3.0/benchmark_reference_four_methods_speedup.png' | relative_url }}" alt="Four aligned panels compare reference-to-fastphylosig runtime ratios for K, lambda, D, and Delta at 50, 100, 500, and 2,000 tips. A dashed horizontal line marks equal runtime; ratios below it mean the reference was faster.">
  </a>
  <figcaption><strong>Measured runtime by statistic and tree size.</strong> The dashed line marks equal median runtime. Click the figure to open the full-size PNG. Results are from Windows 11 x64 with R 4.6.1; setup and tree/data preparation were excluded.</figcaption>
</figure>

## Install and get started {#install}

<div class="install-stack" markdown="1">
<article class="content-card" markdown="1">
### Install v0.3.0

```r
install.packages("remotes")
remotes::install_github("yinwen-ecology/fastphylosig@v0.3.0")
library(fastphylosig)
```

Source installation requires a C++ toolchain. OpenMP is optional.
</article>

<article class="content-card" markdown="1">
### One-function analysis

```r
set.seed(1)
tree <- ape::rtree(30)
trait <- setNames(stats::rnorm(30), tree$tip.label)
fit <- fast_signal(tree, data = trait, method = "K",
                   progress = FALSE)
summary(fit)
plot_signal(fit)
```

The same entry point supports continuous, binary, and categorical traits.
</article>
</div>

## Choose a method

<div class="method-grid">
  <article class="method-card"><h3>K</h3><p>Blomberg's K for continuous traits; supports permutation testing.</p><code>fast_k()</code></article>
  <article class="method-card"><h3>Lambda</h3><p>Pagel's lambda for continuous traits; inspect the likelihood profile when needed.</p><code>fast_lambda()</code></article>
  <article class="method-card"><h3>D</h3><p>Fritz and Purvis D for binary traits, with random and Brownian calibrations.</p><code>fast_d()</code></article>
  <article class="method-card"><h3>Delta</h3><p>Delta for categorical traits with MCMC and Monte Carlo diagnostics.</p><code>fast_delta()</code></article>
</div>

## A convenient path, with explicit control when needed {#workflow}

`fast_signal()` checks and prepares the tree, matches named trait data to tip labels, handles missing values trait by trait, and dispatches to the selected statistic. For repeated analyses or closer inspection, the same steps are available directly.

<figure class="roadmap-figure">
  <a href="{{ '/technical-roadmap.svg' | relative_url }}" aria-label="Open the full-size fastphylosig workflow diagram">
    <img src="{{ '/technical-roadmap.svg' | relative_url }}" alt="Workflow diagram: the high-level fast_signal route performs checks, tree preparation, data matching, and method dispatch; an advanced route exposes each preparation step before calling a method-specific function.">
  </a>
  <figcaption>One high-level call or an explicit preparation path. Open the SVG for the full-size diagram.</figcaption>
</figure>

<details class="advanced-details">
  <summary>Show the explicit preparation sequence</summary>
  <div markdown="1">

```r
checked <- check_tree(tree, signal = "K")
tree_ready <- resolve_tree(tree, signal = "K")
matched <- match_tree_data(tree_ready, data = trait)
ctx <- prepare_tree(matched$tree)

k <- fast_k(ctx, x = matched$data, test = TRUE, nsim = 999,
            progress = FALSE)
```

Methods are selected explicitly; the package does not infer them from trait values. See the [Chinese usage guide](https://github.com/yinwen-ecology/fastphylosig/blob/main/USAGE_zh.md) and the [complete v0.3.0 examples]({{ '/release-0.3.0/README.html' | relative_url }}).
  </div>
</details>

## See the native `plot_signal()` output {#plots}

These examples use the package's own plotting function on simulated demonstration data. Click an image to inspect the original-size PNG; the linked gallery documents all six styles and their diagnostic limits.

<div class="plot-grid">
  <figure class="plot-card">
    <a href="{{ '/release-0.3.0/plot_signal_K_multi_trait.png' | relative_url }}" aria-label="Open full-size multi-trait K plot">
      <img src="{{ '/release-0.3.0/plot_signal_K_multi_trait.png' | relative_url }}" alt="Native plot_signal output showing K permutation distributions for multiple traits.">
    </a>
    <figcaption><strong>K across traits</strong><span>Permutation distributions with observed-value markers.</span></figcaption>
  </figure>
  <figure class="plot-card">
    <a href="{{ '/release-0.3.0/plot_signal_lambda_profile.png' | relative_url }}" aria-label="Open full-size Pagel lambda profile plot">
      <img src="{{ '/release-0.3.0/plot_signal_lambda_profile.png' | relative_url }}" alt="Native plot_signal output showing a Pagel lambda likelihood profile.">
    </a>
    <figcaption><strong>Pagel's lambda</strong><span>Likelihood profile and estimate.</span></figcaption>
  </figure>
  <figure class="plot-card">
    <a href="{{ '/release-0.3.0/plot_signal_D_single_nulls.png' | relative_url }}" aria-label="Open full-size D null comparison plots">
      <img src="{{ '/release-0.3.0/plot_signal_D_single_nulls.png' | relative_url }}" alt="Two native plot_signal panels showing D calibrations against random and Brownian null distributions.">
    </a>
    <figcaption><strong>Fritz–Purvis D</strong><span>Random and Brownian null views.</span></figcaption>
  </figure>
</div>

<p class="gallery-link"><a class="button-link" href="{{ '/release-0.3.0/plot_signal_gallery.html' | relative_url }}">View all six plot styles and the Delta demo</a></p>
<p class="caution-line"><strong>Delta gallery note:</strong> the short-chain Delta panel is for display only. Its diagnostic warning and three-draw permutation P value make it unsuitable for inference.</p>

## Release and validation boundaries

<div class="benchmark-notes">
  <article class="note-card note-card--plain">
    <h3>Small-tree K is slower</h3>
    <p>At 50 and 100 tips, reference K was faster; all observations remain in the figure.</p>
  </article>
  <article class="note-card note-card--caution" id="delta-caveat">
    <h3>Delta diagnostic caveat</h3>
    <p>Reference warnings occurred in 112/120 calls, without final convergence diagnostics; one fast call reached its optimizer limit. Timing ratios do not establish equal accuracy or valid convergence.</p>
  </article>
</div>

<details class="benchmark-details">
  <summary>Benchmark settings and limitations</summary>
  <div markdown="1">

- **K and lambda:** estimates only (`test = FALSE`), without permutation testing.
- **D:** 199 randomizations per call.
- **Delta:** two 10,000-iteration MCMC chains (`thin = 10`, `burn = 100`), without permutation testing.
- Ratios compare median call times; workloads differ across statistics.

See the [full benchmark data, references, settings, and limitations]({{ '/release-0.3.0/README.html' | relative_url }}).
  </div>
</details>

The current release is [fastphylosig v0.3.0](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.0). The 50,000-tip preparation mask-key limitation is deferred to v0.3.x. Benchmark ratios describe the recorded Windows/R environment and these specific workloads; they are not universal speed guarantees or evidence that the four statistics have equal computational work.

<p class="closing-links"><a href="https://github.com/yinwen-ecology/fastphylosig">Source code</a><a href="https://github.com/yinwen-ecology/fastphylosig/issues">Report an issue</a><a href="{{ '/release-0.3.0/README.html' | relative_url }}">Benchmark and usage details</a></p>
