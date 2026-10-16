---
layout: default
title: "fast_hetero()：支系信号异质性案例"
description: "用 512 个模拟物种、8 个预定冠群，演示 fast_hetero 的 K、lambda、D、Delta、树图和诊断。"
lang: zh-CN
page_class: case-page
home_url: /index-zh.html
alternate_url: /fast-hetero-case.html
alternate_lang: en
alternate_label: English
---

<p class="eyebrow">fastphylosig 0.3.4 · 可复现模拟案例</p>

# 比较不同支系的系统发育信号

整棵树的系统发育信号概括了总体格局，但不同支系可能很不一样。`fast_hetero()` 把树划分为单系冠群，复用已有 fast_* 内核计算局部信号，再汇总总体与局部差异。本例用同一组节点比较 K、λ、二元 D 和多状态 Δ，展示怎样从一张树图走向具体的生态与进化解释。

<nav class="page-jumps" aria-label="案例目录">
  <a href="#workflow">技术路线</a><a href="#design">数据与切割</a><a href="#continuous">K 与 λ</a><a href="#categorical">D 与 Δ</a><a href="#diagnostics">诊断与解释</a><a href="#reproduce">代码和下载</a>
</nav>

<div class="case-facts">
  <div><strong>512</strong><span>模拟物种，数据与代码可下载</span></div>
  <div><strong>8 × 64</strong><span>事先按树确定的等龄单系冠群</span></div>
  <div><strong>4</strong><span>指标；普通信号与异质性参考分别解释</span></div>
</div>

## 从整棵树到支系比较 {#workflow}

`fast_hetero()` 先编译树并确定冠群，再按性状处理缺失观测，复用 K、λ、D、ACE 与 Δ 的 C++ 计算路径。总体估计和检验按需计算，最后把局部值、总体值、离散度与诊断汇总到同一对象。

<figure class="roadmap-figure hetero-roadmap">
  <a href="{{ '/fast-hetero-roadmap-zh.svg' | relative_url }}"><img src="{{ '/fast-hetero-roadmap-zh.svg' | relative_url }}" alt="fast_hetero 技术路线：输入与树编译、预定节点或自动选择冠群、按性状筛选观测、共享 C++ 内核扫描、可选总体及参考检验、结果汇总与树图。"></a>
  <figcaption>技术路线：实线为局部描述的主流程，虚线为按需启用的检验与总体分析。<a href="{{ '/fast-hetero-roadmap-zh.svg' | relative_url }}">打开 SVG</a> · <a href="https://github.com/yinwen-ecology/fastphylosig/blob/main/docs/fast-hetero-roadmap-zh.md">可编辑 Mermaid 源码与实现说明</a></figcaption>
</figure>

## 1. 先确定数据设计与冠群 {#design}

树为平衡、超度量二叉树，每条枝长为 1。按根下深度 3 的节点预先选择 8 个不重叠冠群，每群 64 个物种、高度 6，覆盖整棵树。这一步只使用树，不依据性状结果挑选节点。

C1–C4 的连续性状分别按 BM 模拟，分类性状分别按低转移率 ER 模拟；C5–C8 的性状独立随机赋值，让两组支系形成清楚的对照。另加入两个常见数据情形：C1 的连续性状有 3 个缺失值；C8 的二元性状全为 absent、多状态性状全为 A。

```r
tree <- ape::stree(512, "balanced")
tree$edge.length <- rep(1, nrow(tree$edge))
depth <- ape::node.depth.edgelength(tree)
nodes <- which(seq_along(depth) > 512 & depth == 3)
```

真实研究可按科、属或其他有生物学含义的类群指定节点，尽量比较年龄和样本量相近的冠群。没有预定类群时，可用 `min_tips=50, max_tips=80` 自动选择最大不重叠冠群；在本例完整平衡树上，它对应相同的 8 个冠群。选好主切割后，再用少量预定尺度做敏感性比较。

## 2. 连续性状：局部 K 与 λ {#continuous}

`x` 是按物种命名的连续性状向量。下面对相同的 `nodes` 计算局部及总体结果；`overall_test` 是普通整树信号检验，`test` 是冠群间的参考检验。

```r
k <- fast_hetero(tree, x, method = "K", nodes = nodes,
  min_tips = 50, test = TRUE, nsim = 999,
  overall_test = TRUE, progress = FALSE)

lambda <- fast_hetero(tree, x, method = "lambda", nodes = nodes,
  min_tips = 50, test = TRUE, nsim = 999,
  lambda_reference = "REML", overall_test = TRUE, progress = FALSE)

k$clades          # 实际观测数、局部估计、状态与差异
k$overall        # 普通整树信号估计与检验
k$heterogeneity  # 冠群间离散度及其参考检验
plot(k, tree = tree, show.tip.label = FALSE, node_labels = TRUE)
```

本次实际运行的汇总如下。“总体信号 P”回答整树是否具有信号；“异质性参考 P”回答冠群间离散是否超出相应模型的预期。

<!-- GENERATED_SUMMARY_START -->

| 指标 | 整树估计 | 普通整树信号 P | 整树 n | 局部可用/选定 | 冠群离散 | 异质性参考 P |
|---|---|---|---|---|---|---|
| K | 0.3029 | 0.001001 | 509 | 8 / 8 | SD(log(K)) = 0.5861 | 0.001 |
| lambda | 0.5634 | 3.745e-17 | 509 | 8 / 8 | SD(lambda) = 0.483 | 0.001 |
| D | 0.1256 | random: 0/1000; BM: 0.196 | 512 | 7 / 8 | SD(D) = 0.9703 | 未计算（描述） |
| Delta | 5.323 | 0.0201 | 512 | 7 / 8 | SD(log(Delta)) = 4.785 | 未计算（描述） |

K/λ 的参考 P = 0.001 是 999 次模拟、零次超越时的 +1 秩结果，也是这次模拟的分辨率下限。Δ 的普通整体 P 为 4/199 = 0.0201，MCSE 约为 0.00995。

<!-- GENERATED_SUMMARY_END -->

K 和 λ 都揭示出整树汇总之外的支系差异。K 的离散与给定树形的 BM 参考比较；λ 使用**共同局部 λ** 的近似 bootstrap。本例冠群同龄，便于直接比较；实际数据中，也要同时看冠群年龄，因为 λ 的取值与冠群的计量尺度有关。

<figure class="case-figure">
  <a href="{{ '/assets/fast-hetero-case/K-tree.png' | relative_url }}"><img src="{{ '/assets/fast-hetero-case/K-tree.png' | relative_url }}" alt="本案例的 K 原生树图：八个冠群按局部 K 着色，旁列整体信号和 BM 离散参考检验。"></a>
  <figcaption>图 1：八个冠群按局部 K 着色，节点编号与下表对应。右侧给出整树 K、普通信号检验、冠群离散与 BM 参考 P。点击打开大图。</figcaption>
</figure>

<figure class="case-figure">
  <a href="{{ '/assets/fast-hetero-case/lambda-tree.png' | relative_url }}"><img src="{{ '/assets/fast-hetero-case/lambda-tree.png' | relative_url }}" alt="本案例的 lambda 原生树图，同时显示整树 lambda、普通 LRT 和共同局部 lambda 的近似异质性参考 P。"></a>
  <figcaption>图 2：λ 树图显示相同节点上的局部格局，并单列整树估计和共同局部 λ 的近似参考检验。</figcaption>
</figure>

## 3. 分类性状：比较聚集格局与状态组成 {#categorical}

二元性状用 D，多状态性状用 Δ。本例两者都设 `test=FALSE`，不计算其新增异质性 P；`overall_test=TRUE` 仍计算原有普通整树信号检验。

```r
d <- fast_hetero(tree, binary, method = "D", nodes = nodes,
  min_tips = 51, test = FALSE, d_nsim = 1000,
  overall_test = TRUE, progress = FALSE)

delta <- fast_hetero(tree, multi, method = "Delta", nodes = nodes,
  min_tips = 50, test = FALSE, overall_test = TRUE,
  overall_nsim = 199, delta_control = list(model = "ER"),
  progress = FALSE)
```

Δ 使用 LSE、两条各 10,000 步的链、`burn=100, thin=10`，拟合模型预先定为 ER。普通整体检验使用 199 次随机化，结果同时给出 P 与 Monte Carlo 误差；正式分析可按所需概率精度增加次数。

<!-- GENERATED_CLADES_START -->

| 冠群 | 原始节点 | 连续性状 n | K | λ | D | Δ | 状态 |
|---|---|---|---|---|---|---|---|
| C1 | 516 | 61 | 1.125 | 0.9636 | -0.8426 | 91.24 | 可估计 |
| C2 | 579 | 64 | 0.7893 | 0.9964 | -1.361 | 193.7 | 可估计 |
| C3 | 643 | 64 | 1.111 | 0.9254 | -0.6767 | 19.81 | 可估计 |
| C4 | 706 | 64 | 0.8178 | 0.8336 | 0.544 | 565.9 | 可估计 |
| C5 | 771 | 64 | 0.3365 | 7.933e-05 | 0.8975 | 0.001573 | 可估计 |
| C6 | 834 | 64 | 0.3488 | 7.933e-05 | 1.118 | 0.08829 | 可估计 |
| C7 | 898 | 64 | 0.2847 | 7.933e-05 | 0.5317 | 0.1788 | 可估计 |
| C8 | 961 | 64 | 0.3411 | 0.1361 | NA | NA | D/Δ 单状态 |

<!-- GENERATED_CLADES_END -->

C8 仍列在表中：它的状态全部相同，因此 D/Δ 的局部估计为 `NA`、状态为 `undefined`。这样的结果提醒我们先回看性状组成。K/λ 的 C1 则保留 61 个实际观测，其他冠群为 64 个；`min_tips` 检查剔除缺失后的观测数。分类性状没有相同的缺失，仍为每冠群 64 个物种。

阅读时注意指标方向：较大的 K、λ、Δ 通常表示更强的近缘相似格局；D 则是较小的值表示更强聚集，而且 D 可以小于 0 或大于 1。按同一指标比较支系，再结合状态组成和诊断理解差异。

<figure class="case-figure">
  <a href="{{ '/assets/fast-hetero-case/D-tree.png' | relative_url }}"><img src="{{ '/assets/fast-hetero-case/D-tree.png' | relative_url }}" alt="D 描述性树图：七个可用冠群和一个单状态不可用冠群，整体 random 与 Brownian P 单列，条件 ER 异质性检验未执行。"></a>
  <figcaption>图 3：D 的局部聚集格局与总体随机/Brownian 信号检验。C8 的单状态冠群显示为灰色；0/1000 表示本次模拟中零次超过观测值。</figcaption>
</figure>

<figure class="case-figure">
  <a href="{{ '/assets/fast-hetero-case/Delta-tree.png' | relative_url }}"><img src="{{ '/assets/fast-hetero-case/Delta-tree.png' | relative_url }}" alt="多状态 Delta 描述性树图：局部值用对数色阶表示，图例显示原始数值，单状态冠群为灰色，异质性检验未执行。"></a>
  <figcaption>图 4：Δ 用对数色阶展示跨度较大的局部估计，图例仍标原始数值。C8 已被选择，但因单状态而不可用；未进入冠群的祖先枝也显示为灰色。</figcaption>
</figure>

## 4. 结合诊断，形成结果解释 {#diagnostics}

<!-- GENERATED_DIAGNOSTICS_START -->

| 冠群 | 状态计数 | rank R-hat(Δ) | bulk ESS | tail ESS | 相对 MCSE |
|---|---|---|---|---|---|
| C1 | A=50;B=2;C=12 | 1.001 | 1712 | 1682 | 0.58% |
| C2 | A=53;B=5;C=6 | 1 | 1918 | 1984 | 0.52% |
| C3 | A=40;B=7;C=17 | 0.9996 | 1812 | 1587 | 0.56% |
| C4 | A=62;B=1;C=1 | 1.001 | 1844 | 1809 | 0.80% |
| C5 | A=20;B=23;C=21 | 1.001 | 1747 | 1783 | 1.39% |
| C6 | A=17;B=27;C=20 | 0.9992 | 1925 | 1856 | 0.41% |
| C7 | A=26;B=24;C=14 | 1 | 1599 | 1865 | 0.37% |
| 整体 | 三状态 | 1.003 | 1642 | 1518 | 0.17% |

<!-- GENERATED_DIAGNOSTICS_END -->

本次已估计的 Δ 冠群与整体结果未触发 MCMC 警告。表中给出 Δ 序列的 rank R-hat、bulk/tail ESS 和相对 MCSE；完整对象还保留 α、β 参数诊断，以及优化器、数值恢复和告警信息。实际分析据这些信息决定是否延长链，或进一步检查稀有状态。

本例可以写成：“在预先确定的冠群尺度上，C1–C4 整体表现出比后半部分支系更强的近缘相似格局。K 的冠群间离散超出 BM 参考，λ 的离散也超出共同局部参数的近似参考。D 和 Δ 补充展示了分类性状的局部聚集差异；C8 的单状态组成单独报告。”

> **适用范围。** K 采用 BM 参考，λ 采用共同局部冠群参数的近似参考；D 的新增异质性检验仍为实验性，Δ 推荐局部/总体描述而不推荐异质性检验。信号差异可描述分析尺度上的系统发育保守性格局；若要解释稳定选择、生态位保守或速率机制，还需相应证据。方法细节见页末完整指南。

## 5. 下载与完整复现 {#reproduce}

<div class="case-downloads">
  <a href="{{ '/assets/fast-hetero-case/run-case.R' | relative_url }}">完整 R 脚本</a>
  <a href="{{ '/assets/fast-hetero-case/traits.csv' | relative_url }}">性状 CSV</a>
  <a href="{{ '/assets/fast-hetero-case/tree.nwk' | relative_url }}">Newick 树</a>
  <a href="{{ '/assets/fast-hetero-case/crowns.csv' | relative_url }}">冠群节点</a>
  <a href="{{ '/assets/fast-hetero-case/case-results.rds' | relative_url }}">完整 RDS 结果</a>
  <a href="{{ '/assets/fast-hetero-case/warnings.json' | relative_url }}">全部告警</a>
  <a href="{{ '/assets/fast-hetero-case/session-info.txt' | relative_url }}">运行环境</a>
</div>

脚本需要 fastphylosig 0.3.4、ape 和 jsonlite；固定数据种子为 20261016，各方法另用固定分析种子。下载后在仓库根目录运行脚本，或设置 `FASTPHYLOSIG_CASE_OUT` 指向输出目录。完整脚本包含数据生成及所有参数，是复现的准确信息源；本页代码块为对应步骤摘录。

各方法 CSV 也可直接下载：[K]({{ '/assets/fast-hetero-case/K-clades.csv' | relative_url }})、[λ]({{ '/assets/fast-hetero-case/lambda-clades.csv' | relative_url }})、[D]({{ '/assets/fast-hetero-case/D-clades.csv' | relative_url }})、[Δ]({{ '/assets/fast-hetero-case/Delta-clades.csv' | relative_url }})。将地址中的 `-clades.csv` 换为 `-overall.csv` 或 `-heterogeneity.csv` 可查看对应完整标量表。RDS 保留状态列表和诊断，CSV 不代替完整对象。

本案例使用已验收的冻结源码包生成。脚本还验证了绘图前后的随机数状态一致；数据、结果及运行环境一并提供，方便复核。

[返回中文首页]({{ '/index-zh.html' | relative_url }}) · [安装 0.3.4](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.4) · [完整统计范围指南](https://github.com/yinwen-ecology/fastphylosig/blob/v0.3.4/inst/HETERO_GUIDE_zh.md)
