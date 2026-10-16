---
layout: default
title: R 系统发育信号快速分析
description: fastphylosig 0.3.4 计算系统发育信号，并通过 fast_hetero() 比较支系与整树格局。
lang: zh-CN
home_url: /index-zh.html
alternate_url: /
alternate_lang: en
alternate_label: English
---

<nav class="page-jumps" aria-label="页面目录">
  <a href="#speed">实测性能</a><a href="#install">安装</a><a href="#clade-signal">支系信号</a><a href="#workflow">分析流程</a><a href="#plots">绘图案例</a>
</nav>

<section class="hero" id="top" markdown="1">
<p class="eyebrow">fastphylosig · R 包 · v0.3.4</p>

# R 系统发育信号快速分析

<p class="hero-proof"><strong>在 2,000 个末端的基准中，参考实现与 fastphylosig 的中位耗时比最高达到：K 151×、lambda 23,591×、D 224×、Delta 33.6×。</strong>Delta 耗时需结合 <a href="#delta-caveat">MCMC 诊断</a>阅读。</p>

<p class="hero-lead">分析连续、二元和多状态性状的系统发育信号。用 <code>fast_signal()</code> 完成整树分析，用 <code>fast_hetero()</code> 比较支系格局，同时保留覆盖度、诊断与原生树图。</p>

<div class="speed-cards" aria-label="参考实现与 fastphylosig 的实测耗时比">
  <article class="speed-card"><p class="speed-card__method">Blomberg's K</p><p class="speed-card__value">135–151×</p><p class="speed-card__note">2,000 个末端</p></article>
  <article class="speed-card"><p class="speed-card__method">Pagel's lambda</p><p class="speed-card__value">20,697–23,591×</p><p class="speed-card__note">2,000 个末端</p></article>
  <article class="speed-card"><p class="speed-card__method">Fritz–Purvis D</p><p class="speed-card__value">218–224×</p><p class="speed-card__note">2,000 个末端</p></article>
  <article class="speed-card speed-card--delta"><p class="speed-card__method">多状态 Delta</p><p class="speed-card__value">23.9–33.6×</p><p class="speed-card__note">2,000 个末端 · 仅调用耗时 · <a href="#delta-caveat">诊断告警说明</a></p></article>
</div>

<p class="hero-footnote">耗时比为“参考实现中位耗时 ÷ fastphylosig 中位耗时”。2,000 末端的区间覆盖三种预设信号场景。</p>
</section>

## 用 `fast_hetero()` 比较不同支系 {#clade-signal}

指定具有生物学含义的单系冠群，或按大小限制自动选择最大不重叠冠群。函数复用已有 fast_* 内核计算局部信号，按需与整树比较，并按局部估计为树上冠群着色。结果保留实际观测数、缺失信息、分类性状组成，以及拟合和 MCMC 诊断。

<div class="hetero-feature" markdown="1">
<div markdown="1">

### 512 个物种的可复现案例

8 个冠群、3 类性状、4 种指标。沿同一组节点比较 K、lambda、二元 D 和多状态 Delta，并查看一个刻意设置的单状态冠群怎样保留在结果中。

K 使用 BM 离散参考；lambda 的共同局部冠群参考是近似方法。D 的新增异质性检验保留实验性标注；Delta 推荐用于局部和总体描述及诊断。按生物学问题的尺度比较系统发育保守性格局。

[阅读 fast_hetero() 案例]({{ '/fast-hetero-case-zh.html' | relative_url }}){: .button-link}

[English case study]({{ '/fast-hetero-case.html' | relative_url }}) · [技术路线说明]({{ '/fast-hetero-case-zh.html' | relative_url }}#workflow)

</div>
<figure>
  <a href="{{ '/fast-hetero-case-zh.html' | relative_url }}" aria-label="阅读支系信号案例"><img src="{{ '/assets/fast-hetero-case/K-tree.png' | relative_url }}" alt="fast_hetero 原生树图：512 个模拟物种的八个冠群，按局部 K 着色，同时显示整树 K、普通信号 P 和单列的 BM 离散参考 P。"></a>
  <figcaption>案例中的原生 <code>plot()</code> 输出；案例页提供数值、诊断、代码和数据。</figcaption>
</figure>
</div>

<figure class="roadmap-figure hetero-roadmap">
  <a href="{{ '/fast-hetero-roadmap-zh.svg' | relative_url }}" aria-label="打开 fast_hetero 技术路线图"><img src="{{ '/fast-hetero-roadmap-zh.svg' | relative_url }}" alt="fast_hetero 技术路线：编译树、选择冠群、共享 C++ 内核计算、可选整树与参考检验、汇总结果和原生树图。"></a>
  <figcaption><code>fast_hetero()</code> 先选择冠群、复用方法内核，再汇总局部估计、可选整树分析与诊断。<a href="{{ '/fast-hetero-case-zh.html' | relative_url }}#workflow">结合案例阅读技术路线</a> · <a href="{{ '/fast-hetero-roadmap.svg' | relative_url }}">English diagram</a></figcaption>
</figure>

## 配对运行时间比较 {#speed}

基准版本为 v0.3.0：4 种指标、50–2,000 个末端、3 种信号场景，每种设置进行 10 次配对调用。耗时覆盖树与数据准备后的方法调用。

<figure class="benchmark-figure">
  <a class="figure-open" href="{{ '/release-0.3.0/benchmark_reference_four_methods_speedup.png' | relative_url }}" aria-label="打开四指标性能基准大图"><img src="{{ '/release-0.3.0/benchmark_reference_four_methods_speedup.png' | relative_url }}" alt="四个面板比较 K、lambda、D 和 Delta 在 50、100、500、2,000 个末端下的参考实现与 fastphylosig 耗时比，虚线表示耗时相同。"></a>
  <figcaption><strong>不同指标和树规模下的实测耗时比。</strong>虚线表示中位耗时相等；比值低于虚线时，参考实现更快。点击打开 PNG 大图。环境为 Windows 11 x64、R 4.6.1，不包含设置及树与数据准备耗时。</figcaption>
</figure>

## 安装与开始分析 {#install}

<div class="install-stack" markdown="1">
<article class="content-card" markdown="1">

### 安装 v0.3.4

[v0.3.4 发布版](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.4) 提供完整的支系分析流程，保留实验性标注，并减少 D/lambda 的重复计算。以下命令安装已冻结并完成检查的源码附件。

```r
install.packages(
  "https://github.com/yinwen-ecology/fastphylosig/releases/download/v0.3.4/fastphylosig_0.3.4.tar.gz",
  repos = NULL, type = "source"
)
library(fastphylosig)
```

源码安装需要 C++ 工具链；OpenMP 为可选项。
</article>
<article class="content-card" markdown="1">

### 一次调用完成分析

```r
set.seed(1)
tree <- ape::rtree(30)
trait <- setNames(stats::rnorm(30), tree$tip.label)
fit <- fast_signal(tree, data = trait, method = "K",
                   progress = FALSE)
summary(fit)
plot_signal(fit)
```

同一入口支持连续、二元与多状态性状。
</article>
</div>

## 选择分析方法

<div class="method-grid">
  <article class="method-card"><h3>K</h3><p>连续性状的 Blomberg's K，支持随机置换检验。</p><code>fast_k()</code></article>
  <article class="method-card"><h3>Lambda</h3><p>连续性状的 Pagel's lambda，可按需查看似然剖面。</p><code>fast_lambda()</code></article>
  <article class="method-card"><h3>D</h3><p>二元性状的 Fritz–Purvis D，提供随机与 Brownian 参考。</p><code>fast_d()</code></article>
  <article class="method-card"><h3>Delta</h3><p>分类性状的 Delta，提供 MCMC 与 Monte Carlo 诊断。</p><code>fast_delta()</code></article>
</div>

## 简便入口与显式控制 {#workflow}

`fast_signal()` 检查并准备树，将命名性状与末端标签匹配，按性状处理缺失值，再调用指定指标。重复分析或需要进一步检查时，也可以直接执行各步骤。

<figure class="roadmap-figure">
  <a href="{{ '/technical-roadmap-zh.svg' | relative_url }}" aria-label="打开 fast_signal 分析流程图"><img src="{{ '/technical-roadmap-zh.svg' | relative_url }}" alt="fast_signal 高层入口自动检查、准备树、匹配数据与选择计算路径；显式路径提供各准备函数，再调用具体指标。"></a>
  <figcaption>可使用高层入口，也可显式执行准备流程。点击 SVG 查看大图。</figcaption>
</figure>

<details class="advanced-details">
  <summary>查看显式准备流程</summary>
  <div markdown="1">

```r
checked <- check_tree(tree, signal = "K")
tree_ready <- resolve_tree(tree, signal = "K")
matched <- match_tree_data(tree_ready, data = trait)
ctx <- prepare_tree(matched$tree)

k <- fast_k(ctx, x = matched$data, test = TRUE, nsim = 999,
            progress = FALSE)
```

分析方法由用户显式指定。更多示例见[中文使用指南](https://github.com/yinwen-ecology/fastphylosig/blob/main/USAGE_zh.md)与 [v0.3.0 完整案例]({{ '/release-0.3.0/README.html' | relative_url }})。
  </div>
</details>

## 原生 `plot_signal()` 绘图示例 {#plots}

这些 v0.3.0 示例使用包内绘图函数和模拟数据。点击图像可打开原尺寸 PNG；完整画廊说明六种样式及其诊断。

<div class="plot-grid">
  <figure class="plot-card"><a href="{{ '/release-0.3.0/plot_signal_K_multi_trait.png' | relative_url }}"><img src="{{ '/release-0.3.0/plot_signal_K_multi_trait.png' | relative_url }}" alt="多个连续性状的 K 置换分布和观测值标记。"></a><figcaption><strong>多性状 K</strong><span>置换分布与观测值标记。</span></figcaption></figure>
  <figure class="plot-card"><a href="{{ '/release-0.3.0/plot_signal_lambda_profile.png' | relative_url }}"><img src="{{ '/release-0.3.0/plot_signal_lambda_profile.png' | relative_url }}" alt="Pagel lambda 的似然剖面及其估计值。"></a><figcaption><strong>Pagel's lambda</strong><span>似然剖面与估计值。</span></figcaption></figure>
  <figure class="plot-card"><a href="{{ '/release-0.3.0/plot_signal_D_single_nulls.png' | relative_url }}"><img src="{{ '/release-0.3.0/plot_signal_D_single_nulls.png' | relative_url }}" alt="D 的随机与 Brownian 参考分布。"></a><figcaption><strong>Fritz–Purvis D</strong><span>随机与 Brownian 参考。</span></figcaption></figure>
</div>

<p class="gallery-link"><a class="button-link" href="{{ '/release-0.3.0/plot_signal_gallery.html' | relative_url }}">查看六种绘图样式与 Delta 示例</a></p>
<p class="caution-line"><strong>Delta 画廊：</strong>该示例使用短链，含诊断告警，置换 P 仅来自三次模拟；画廊记录了具体设置。</p>

## 验证与基准细节

<div class="benchmark-notes">
  <article class="note-card note-card--plain"><h3>小树上的 K 耗时</h3><p>在 50 和 100 个末端时，参考 K 实现更快；图中保留了这些结果。</p></article>
  <article class="note-card note-card--caution" id="delta-caveat"><h3>Delta 诊断</h3><p>Delta 比值比较调用耗时：参考实现 120 次调用中有 112 次告警，且未提供最终收敛诊断；fast 实现有一次达到优化器限制。具体情况见链接中的完整结果。</p></article>
</div>

<details class="benchmark-details">
  <summary>基准设置</summary>
  <div markdown="1">

- **K 和 lambda：**仅估计（`test = FALSE`），不含置换检验。
- **D：**每次调用 199 次随机化。
- **Delta：**两条各 10,000 步的 MCMC 链（`thin = 10, burn = 100`），不含置换检验。
- 比值按中位调用耗时计算，各指标的计算任务不同。

查看[完整基准数据、参考实现和参数]({{ '/release-0.3.0/README.html' | relative_url }})。
  </div>
</details>

当前发布版为 [fastphylosig v0.3.4](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.4)。独立 Windows 检查通过重新编译、示例、13,679 项断言及 PDF/HTML 手册检查，无错误、无警告，有一条入库检查的“New submission”提示。一项历史跨安装测试跳过；开发版 K 四模式桥接另行完成检查。维护者已确认 Linux 检查通过，本轮未声明其他平台已验证。上面的性能图明确来自 v0.3.0，并非 v0.3.4 的新基准。详见[发布说明](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.4)与 [fast_hetero 案例]({{ '/fast-hetero-case-zh.html' | relative_url }})。

<p class="closing-links"><a href="https://github.com/yinwen-ecology/fastphylosig">源码</a><a href="https://github.com/yinwen-ecology/fastphylosig/issues">问题反馈</a><a href="{{ '/release-0.3.0/README.html' | relative_url }}">基准与使用细节</a></p>
