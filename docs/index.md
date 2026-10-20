---
layout: default
title: Fast phylogenetic signal analysis in R
description: fastphylosig 0.3.4 measures phylogenetic signal and compares local and whole-tree patterns with fast_hetero().
lang: en
home_url: /
alternate_url: /index-zh.html
alternate_lang: zh-CN
alternate_label: 中文
---

<nav class="page-jumps" aria-label="Page sections">
  <a href="#install">Install</a>
  <a href="#workflow">Workflow</a>
  <a href="#speed">Measured speed</a>
  <a href="#plots">Plot examples</a>
  <a href="#clade-signal">Clade-scale signal</a>
</nav>

<section class="hero" id="top" markdown="1">
<p class="eyebrow">fastphylosig · R package · v0.3.4</p>

# Fast phylogenetic signal analysis in R

<p class="hero-proof"><strong>At 2,000 tips, fastphylosig computes phylogenetic signal up to 23,591× faster.</strong></p>

<p class="hero-lead">Phylogenetic signal measures how similar a trait is among closely related species. <code>fast_signal()</code> computes K, lambda, D or Delta in one call, and the speed advantage grows with tree size (<a href="#speed">see the benchmark</a>). <code>fast_hetero()</code> compares signal among clades within a large tree (<a href="#clade-signal">see the 512-species case study</a>).</p>

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
    <p class="speed-card__note">at 2,000 tips</p>
  </article>
</div>

<p class="hero-footnote">Ratios are reference median time ÷ fastphylosig median time, so a value above 1 means fastphylosig was faster. The references are phytools 2.5.2 for K and lambda, caper for D, and the published Borges (2019) code for Delta. Each range covers three designed signal scenarios (Brownian/noise mixtures with weights 0, 0.5 and 0.95).</p>
</section>

## Install and get started {#install}

<div class="install-stack" markdown="1">
<article class="content-card" markdown="1">
### Install v0.3.4

The [v0.3.4 release](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.4) is the latest version. It provides the full clade-scale workflow and improves D/lambda performance.

```r
install.packages(
  "https://github.com/yinwen-ecology/fastphylosig/releases/download/v0.3.4/fastphylosig_0.3.4.tar.gz",
  repos = NULL, type = "source"
)
library(fastphylosig)
```

Source installation needs a C++17 toolchain and, on Windows, Rtools. It also requires ape, Rcpp and RcppArmadillo.
</article>

<article class="content-card" markdown="1">
### One call for K, lambda, D or Delta

<p class="code-label">Simulated example data</p>

```r
set.seed(1)
tree <- ape::rtree(30)
trait <- setNames(stats::rnorm(30), tree$tip.label)
binary <- setNames(as.integer(trait > stats::median(trait)), tree$tip.label)
state <- setNames(sample(letters[1:3], 30, replace = TRUE), tree$tip.label)
```

<div class="two-column" markdown="1">
<div markdown="1">

<p class="code-label">K</p>

```r
k <- fast_signal(tree, data = trait, method = "K",
                 test = TRUE, nsim = 999,
                 return_sim = TRUE, progress = FALSE)
```

</div>
<div markdown="1">

<p class="code-label">lambda</p>

```r
lambda <- fast_signal(tree, data = trait, method = "lambda",
                      test = TRUE, progress = FALSE)
```

</div>
<div markdown="1">

<p class="code-label">D</p>

```r
d <- fast_signal(tree, data = binary, method = "D",
                 test = TRUE, nsim = 999,
                 return_sim = TRUE, progress = FALSE)
```

</div>
<div markdown="1">

<p class="code-label">Delta</p>

```r
delta <- fast_signal(tree, data = state, method = "Delta",
                     test = TRUE, nsim = 199, progress = FALSE)
```

</div>
</div>
</article>
</div>

## Technical workflow {#workflow}

<code>fast_signal()</code> combines tree checking, tree–trait matching, missing-value handling and signal calculation. If needed, the four steps can also be run separately.

<figure class="roadmap-figure">
  <a href="{{ '/technical-roadmap.svg' | relative_url }}" aria-label="Open the full-size fastphylosig workflow diagram">
    <img src="{{ '/technical-roadmap.svg' | relative_url }}" alt="Workflow diagram: the high-level fast_signal route performs checks, tree preparation, data matching, and method dispatch; an advanced route exposes each preparation step before calling a method-specific function.">
  </a>
  <figcaption>Open the SVG for the full-size diagram.</figcaption>
</figure>

<details class="advanced-details">
  <summary>Run the four steps separately</summary>
  <div markdown="1">

```r
checked <- check_tree(tree, signal = "K")
tree_ready <- resolve_tree(tree, signal = "K")
matched <- match_tree_data(tree_ready, data = trait)
ctx <- prepare_tree(matched$tree)

k <- fast_k(ctx, x = matched$data, test = TRUE, nsim = 999,
            progress = FALSE)
```

See the [Chinese usage guide](https://github.com/yinwen-ecology/fastphylosig/blob/main/USAGE_zh.md) and the [complete v0.3.0 examples]({{ '/release-0.3.0/README.html' | relative_url }}).
  </div>
</details>

## Runtime against the standard implementations {#speed}

K and lambda are compared with phytools 2.5.2 `phylosig()`, D with caper 1.0.4 `phylo.d()`, and Delta with the [published Borges (2019) code](https://github.com/mrborges23/delta_statistic). Trees of 50 to 2,000 tips, three signal scenarios, and 10 paired replicates per setting (20 timed calls). Timing starts after tree and data preparation.

<figure class="benchmark-figure">
  <a class="figure-open" href="{{ '/release-0.3.0/benchmark_reference_four_methods_speedup.png' | relative_url }}" aria-label="Open the full-size four-statistic benchmark figure">
    <img src="{{ '/release-0.3.0/benchmark_reference_four_methods_speedup.png' | relative_url }}" alt="Four aligned panels compare reference-to-fastphylosig runtime ratios for K, lambda, D, and Delta at 50, 100, 500, and 2,000 tips. A dashed horizontal line marks equal runtime; ratios below it mean the standard implementation was faster.">
  </a>
  <figcaption><strong>Measured runtime by statistic and tree size.</strong> The dashed line marks equal median runtime; ratios below it mean the standard implementation was faster. Click the figure to open the full-size PNG. Results are from Windows 11 x64 with R 4.6.1; setup and tree/data preparation were excluded.</figcaption>
</figure>

## `plot_signal` {#plots}

Plots from simulated data.

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
  <figure class="plot-card">
    <a href="{{ '/release-0.3.0/plot_signal_Delta_categorical_demo.png' | relative_url }}" aria-label="Open full-size Delta permutation demo plot">
      <img src="{{ '/release-0.3.0/plot_signal_Delta_categorical_demo.png' | relative_url }}" alt="Delta permutation demo showing the permuted distribution and the observed Delta marked by a vertical line.">
    </a>
    <figcaption><strong>Delta</strong><span>Permutation distribution and observed value.</span></figcaption>
  </figure>
</div>

<p class="gallery-link"><a class="button-link" href="{{ '/release-0.3.0/plot_signal_gallery.html' | relative_url }}">View all six plot styles and the Delta demo</a></p>

## Compare clades within a large tree with `fast_hetero()` {#clade-signal}

In a large tree (say more than 1,000 tips), a single signal value cannot describe how lineages differ: some clades may be more conserved, others more labile. I developed `fast_hetero()` to split the tree into crown clades (a node and all its descendants), estimate signal for each, and set the local estimates beside the whole-tree value. You can name the clades yourself, or let the function choose the largest non-overlapping clades within a size range. Results keep the real sample size of each clade, its missing observations, its state composition and its fit and MCMC diagnostics.

<figure class="roadmap-figure hetero-roadmap">
  <a href="{{ '/fast-hetero-roadmap.svg' | relative_url }}" aria-label="Open the fast_hetero technical workflow"><img src="{{ '/fast-hetero-roadmap.svg' | relative_url }}" alt="fast_hetero technical workflow, from tree compilation and crown selection through shared C++ kernels to local and whole-tree results, optional reference tests and native plots."></a>
  <figcaption><code>fast_hetero()</code> workflow: select crowns, reuse the method kernels, then combine local estimates, optional whole-tree analysis and diagnostics. <a href="{{ '/fast-hetero-case.html' | relative_url }}#workflow">Follow the workflow in the case study</a> · <a href="{{ '/fast-hetero-roadmap-zh.svg' | relative_url }}">中文路线图</a></figcaption>
</figure>

<div class="hetero-feature" markdown="1">

### Demonstration case

<figure>
  <a href="{{ '/fast-hetero-case.html' | relative_url }}" aria-label="Read the clade-scale signal case study">
    <img src="{{ '/assets/fast-hetero-case/K-tree.png' | relative_url }}" alt="Native fast_hetero tree plot for eight selected crowns in a simulated 512-species tree, with whole-tree K, its ordinary signal P and a separate BM dispersion P.">
  </a>
  <figcaption>Built-in <code>plot()</code> output from the case study; numbers, diagnostics, code and data are on the case page.</figcaption>
</figure>

[Read the full case study]({{ '/fast-hetero-case.html' | relative_url }}){: .button-link}

</div>

## Other notes

<div class="benchmark-notes benchmark-notes--single">
  <article class="note-card note-card--caution" id="delta-caveat">
    <h3>Delta diagnostics</h3>
    <p>Delta ratios compare returned calls, not accuracy. The reference emitted numerical warnings in 112 of 120 calls and exposes no final convergence check, and one fast call reached its optimizer limit. No warned call was removed or re-run, so these ratios do not show that the two implementations are equally accurate or equally converged. Warnings are retained in the <a href="{{ '/release-0.3.0/benchmark_reference_d_delta_runs.csv' | relative_url }}">raw timed-call records</a>, with counts by type in the <a href="{{ '/release-0.3.0/benchmark_reference_four_methods_qa.html' | relative_url }}">benchmark QA report</a>.</p>
  </article>
</div>

The current release is [fastphylosig v0.3.4](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.4). The runtime figures above come from v0.3.0 and were not re-run for v0.3.4, because the core computational code did not change. See the [release notes](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.4) and the [fast_hetero case study]({{ '/fast-hetero-case.html' | relative_url }}).

<p class="closing-links"><a href="https://github.com/yinwen-ecology/fastphylosig">Source code</a><a href="https://github.com/yinwen-ecology/fastphylosig/issues">Report an issue</a><a href="{{ '/release-0.3.0/README.html' | relative_url }}">Benchmark and usage details</a></p>
