---
layout: default
title: R 系统发育信号快速分析
description: fastphylosig 0.3.5 计算系统发育信号，并通过 fast_hetero() 比较支系与整树格局。
lang: zh-CN
home_url: /index-zh.html
alternate_url: /
alternate_lang: en
alternate_label: English
---

<nav class="page-jumps" aria-label="页面目录">
  <a href="#install">安装</a><a href="#workflow">技术路线</a><a href="#speed">实测性能</a><a href="#plots">绘图案例</a><a href="#clade-signal">支系信号</a>
</nav>

<section class="hero" id="top" markdown="1">
<p class="eyebrow">fastphylosig · R 包 · v0.3.5</p>

# R 系统发育信号快速分析

<p class="hero-proof"><strong>在 2,000 个物种/末端节点规模上，fastphylosig 计算系统发育信号的最高可提速达 23,591 倍。</strong></p>

<p class="hero-lead">系统发育信号衡量性状在近缘物种间相似的程度。<code>fast_signal()</code> 用一行代码就能算出 K、lambda、D 或 Delta，树越大，速度优势越明显（<a href="#speed">查看耗时基准</a>）。<code>fast_hetero()</code> 用于比较大型系统发育树内部各支系的信号差异（<a href="#clade-signal">查看 512 物种案例</a>）。</p>

<div class="speed-cards" aria-label="参考实现与 fastphylosig 的实测耗时比">
  <article class="speed-card"><p class="speed-card__method">Blomberg's K</p><p class="speed-card__value">135–151×</p><p class="speed-card__note">2,000 个末端</p></article>
  <article class="speed-card"><p class="speed-card__method">Pagel's lambda</p><p class="speed-card__value">20,697–23,591×</p><p class="speed-card__note">2,000 个末端</p></article>
  <article class="speed-card"><p class="speed-card__method">Fritz–Purvis D</p><p class="speed-card__value">218–224×</p><p class="speed-card__note">2,000 个末端</p></article>
  <article class="speed-card speed-card--delta"><p class="speed-card__method">多状态 Delta</p><p class="speed-card__value">23.9–33.6×</p><p class="speed-card__note">2,000 个末端</p></article>
</div>

<p class="hero-footnote">耗时比为“参考实现中位耗时 ÷ fastphylosig 中位耗时”，比值大于 1 表示 fastphylosig 更快。参考实现分别是 phytools 2.5.2（K、lambda）、caper（D）和 Borges (2019) 的公开代码（Delta）。每个区间覆盖三种预设信号场景（Brownian 与噪声混合，权重 0、0.5、0.95）。</p>
</section>

## 安装与开始分析 {#install}

<div class="install-stack" markdown="1">
<article class="content-card" markdown="1">

### 安装 v0.3.5

[v0.3.5 发布版](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.5) 是最新的版本，提供完整的支系分析流程，并优化了 D/lambda 的计算性能。

```r
install.packages(
  "https://github.com/yinwen-ecology/fastphylosig/releases/download/v0.3.5/fastphylosig_0.3.5.tar.gz",
  repos = NULL, type = "source"
)
library(fastphylosig)
```

源码安装需要 C++17 工具链，Windows 上需要安装 Rtools，并依赖 ape、Rcpp 和 RcppArmadillo。
</article>
<article class="content-card" markdown="1">

### 一行代码完成 K/lambda/D/Delta 分析

<p class="code-label">模拟数据示例</p>

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

## 技术路线 {#workflow}

<code>fast_signal()</code> 整合了系统发育树检查、树与性状匹配、缺失值处理、系统发育信号计算功能。如果需要，也可以单独执行四个步骤。

<figure class="roadmap-figure">
  <a href="{{ '/technical-roadmap-zh.svg' | relative_url }}" aria-label="打开 fast_signal 分析流程图"><img src="{{ '/technical-roadmap-zh.svg' | relative_url }}" alt="fast_signal 高层入口自动检查、准备树、匹配数据与选择计算路径；显式路径提供各准备函数，再调用具体指标。"></a>
  <figcaption>点击 SVG 查看大图。</figcaption>
</figure>

<details class="advanced-details">
  <summary>四个步骤单独执行</summary>
  <div markdown="1">

```r
checked <- check_tree(tree, signal = "K")
tree_ready <- resolve_tree(tree, signal = "K")
matched <- match_tree_data(tree_ready, data = trait)
ctx <- prepare_tree(matched$tree)

k <- fast_k(ctx, x = matched$data, test = TRUE, nsim = 999,
            progress = FALSE)
```

更多示例见[中文使用指南](https://github.com/yinwen-ecology/fastphylosig/blob/main/USAGE_zh.md)与 [v0.3.0 完整案例]({{ '/release-0.3.0/README.html' | relative_url }})。
  </div>
</details>

## 与普通方法运行时间比较 {#speed}

K 和 lambda 对比 phytools 2.5.2 的 `phylosig()`，D 对比 caper 1.0.4 的 `phylo.d()`，Delta 对比 Borges (2019) 的[公开代码](https://github.com/mrborges23/delta_statistic)。50 至 2,000 个末端、3 种信号场景，每种设置 10 组配对重复（共 20 次计时调用）。计时从树与数据准备完成后开始。

<figure class="benchmark-figure">
  <a class="figure-open" href="{{ '/release-0.3.0/benchmark_reference_four_methods_speedup.png' | relative_url }}" aria-label="打开四指标性能基准大图"><img src="{{ '/release-0.3.0/benchmark_reference_four_methods_speedup.png' | relative_url }}" alt="四个面板比较 K、lambda、D 和 Delta 在 50、100、500、2,000 个末端下的参考实现与 fastphylosig 耗时比，虚线表示耗时相同。"></a>
  <figcaption><strong>不同指标和树规模下的实测耗时比。</strong>虚线表示中位耗时相等；比值低于虚线时，普通方法实现更快。点击打开 PNG 大图。计算环境为 Windows 11 x64、R 4.6.1，不包含设置及树与数据准备耗时。</figcaption>
</figure>

## `plot_signal` {#plots}

基于模拟数据的绘图

<div class="plot-grid">
  <figure class="plot-card"><a href="{{ '/release-0.3.0/plot_signal_K_multi_trait.png' | relative_url }}"><img src="{{ '/release-0.3.0/plot_signal_K_multi_trait.png' | relative_url }}" alt="多个连续性状的 K 置换分布和观测值标记。"></a><figcaption><strong>多性状 K</strong><span>置换分布与观测值标记。</span></figcaption></figure>
  <figure class="plot-card"><a href="{{ '/release-0.3.0/plot_signal_lambda_profile.png' | relative_url }}"><img src="{{ '/release-0.3.0/plot_signal_lambda_profile.png' | relative_url }}" alt="Pagel lambda 的似然剖面及其估计值。"></a><figcaption><strong>Pagel's lambda</strong><span>似然剖面与估计值。</span></figcaption></figure>
  <figure class="plot-card"><a href="{{ '/release-0.3.0/plot_signal_D_single_nulls.png' | relative_url }}"><img src="{{ '/release-0.3.0/plot_signal_D_single_nulls.png' | relative_url }}" alt="D 的随机与 Brownian 参考分布。"></a><figcaption><strong>Fritz–Purvis D</strong><span>随机与 Brownian 参考。</span></figcaption></figure>
  <figure class="plot-card"><a href="{{ '/release-0.3.0/plot_signal_Delta_categorical_demo.png' | relative_url }}"><img src="{{ '/release-0.3.0/plot_signal_Delta_categorical_demo.png' | relative_url }}" alt="多状态 Delta 的置换分布演示，竖线标出观测 Delta 的位置。"></a><figcaption><strong>Delta</strong><span>置换分布与观测值位置。</span></figcaption></figure>
</div>

<p class="gallery-link"><a class="button-link" href="{{ '/release-0.3.0/plot_signal_gallery.html' | relative_url }}">查看六种绘图样式与 Delta 示例</a></p>

## 用 `fast_hetero()` 比较大树内不同支系的差异 {#clade-signal}

在大树情况下（例如超过 1,000 个末端），仅凭一个系统发育信号值难以描述内部不同分支间的差异，有的分支可能更保守，有的更发散。因此我开发了 <code>fast_hetero()</code>，把大树拆成冠群（一个节点及其全部后代），逐群估计信号，再把局部估计与整树值并列。冠群可以由你指定，也可以在给定大小范围内让函数选择互不重叠的最大冠群。结果保留每个冠群的实际样本量、缺失观测、状态组成，以及拟合与 MCMC 诊断。

<figure class="roadmap-figure hetero-roadmap">
  <a href="{{ '/fast-hetero-roadmap-zh.svg' | relative_url }}" aria-label="打开 fast_hetero 技术路线图"><img src="{{ '/fast-hetero-roadmap-zh.svg' | relative_url }}" alt="fast_hetero 技术路线：编译树、选择冠群、共享 C++ 内核计算、可选整树与参考检验、汇总结果和原生树图。"></a>
  <figcaption><code>fast_hetero()</code> 先选择冠群、复用方法内核，再汇总局部估计、可选整树分析与诊断。<a href="{{ '/fast-hetero-case-zh.html' | relative_url }}#workflow">结合案例阅读技术路线</a> · <a href="{{ '/fast-hetero-roadmap.svg' | relative_url }}">English diagram</a></figcaption>
</figure>

<div class="hetero-feature" markdown="1">

### 示范案例

<figure>
  <a href="{{ '/fast-hetero-case-zh.html' | relative_url }}" aria-label="阅读支系信号案例"><img src="{{ '/assets/fast-hetero-case/K-tree.png' | relative_url }}" alt="fast_hetero 原生树图：512 个模拟物种的八个冠群，按局部 K 着色，同时显示整树 K、常规信号检验的 P 值和单列的 BM 离散参考 P。"></a>
  <figcaption>案例中的原生 <code>plot()</code> 输出；数值、诊断、代码和数据见案例页。</figcaption>
</figure>

[阅读完整案例]({{ '/fast-hetero-case-zh.html' | relative_url }}){: .button-link}

</div>

## 其他说明

<div class="benchmark-notes benchmark-notes--single">
  <article class="note-card note-card--caution" id="delta-caveat"><h3>Delta 诊断</h3><p>Delta 比值比较的是返回调用，不是精度。参考实现在 120 次调用中有 112 次出现数值告警，且不提供最终收敛检查；fast 实现有一次达到优化器限制。带告警的调用没有被剔除或重跑，因此这些比值不能说明两种实现的精度或收敛程度相当。告警记录保留在<a href="{{ '/release-0.3.0/benchmark_reference_d_delta_runs.csv' | relative_url }}">原始计时记录</a>中，分类汇总见<a href="{{ '/release-0.3.0/benchmark_reference_four_methods_qa.html' | relative_url }}">基准 QA 报告</a>。</p></article>
</div>

当前发布版为 [fastphylosig v0.3.5](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.5)。上面的性能数据来自 v0.3.0，由于核心计算代码未改变，未针对 v0.3.4 重跑。详见[发布说明](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.5)与 [fast_hetero 案例]({{ '/fast-hetero-case-zh.html' | relative_url }})。

<p class="closing-links"><a href="https://github.com/yinwen-ecology/fastphylosig">源码</a><a href="https://github.com/yinwen-ecology/fastphylosig/issues">问题反馈</a><a href="{{ '/release-0.3.0/README.html' | relative_url }}">基准与使用细节</a></p>
