# fastphylosig 中文使用指南

`fastphylosig` 是用于快速计算连续、二元和分类性状的系统发育信号的R包。由于使用了优化的C++内核，其速度比标准方法（phytools、caper等）要快最高达 23,591 倍（2000物种计算λ，和phytools比较，详见[中文性能介绍](https://yinwen-ecology.github.io/fastphylosig/index-zh.html#speed)），尤其适用于超过10000个tips的大型系统发育树。一般而言，只需要使用一行代码（统一入口函数 `fast_signal()`）即可进行系统发育信号分析；当你需要逐步检查树结构、调整特定方法参数或进行深度诊断时，再调用对应的专业函数。同时，对于大型系统发育树，本R包还提供了 `fast_hetero()` 函数用于分析内部的系统发育信号异质性（参阅[中文异质性介绍](https://yinwen-ecology.github.io/fastphylosig/index-zh.html#clade-signal)）。

当前发布版本为 **0.3.5**。该版本修正了 ARD 模型的祖先状态概率与离散 ACE 速率的标准误，参阅[补丁说明与重算建议](docs/release-0.3.5/RELEASE_NOTES_zh.md)。

---

## 安装

从源码安装需要 C++ 编译工具链（如 Windows 下的 Rtools），并依赖 `ape`、`Rcpp` 和 `RcppArmadillo` 包。OpenMP 为可选的并行加速项。

```r
# 首次安装时先安装依赖
install.packages(c("ape", "Rcpp", "RcppArmadillo"))

# 安装 0.3.5 源码包
install.packages(
  "https://github.com/yinwen-ecology/fastphylosig/releases/download/v0.3.5/fastphylosig_0.3.5.tar.gz",
  repos = NULL, type = "source"
)
library(fastphylosig)

# 也可以通过 remotes 安装指定版本标签
# remotes::install_github("yinwen-ecology/fastphylosig@v0.3.5")
```

> **环境支持**：声明支持 R 4.1.0 及以上版本。正式版已通过 Linux、Windows、macOS 及禁用 OpenMP 的 CI 检查。

---

## 最快用法：一个命令完成分析

准备好 `phylo` 格式的进化树和带有物种名称的性状数据后，指定方法即可直接运行：

```r
set.seed(1)
tree <- ape::rtree(30)
trait <- setNames(stats::rnorm(30), tree$tip.label)

# 显式选择 Blomberg's K
fit <- fast_signal(tree, data = trait, method = "K", progress = FALSE)
summary(fit)
```

`fast_signal()` 会自动检查树结构、匹配物种名、按性状处理缺失值，并执行所选的信号计算。

### 支持的方法速查

函数不会自动选择计算方法，你需要通过 `method` 显式指定：

| 分析目标 | `method` | 适用性状类型 | 主要输出与检验 |
|---|---|---|---|
| Blomberg's K | `"K"` | 连续性状 | K 估计值、可选的 Tip 随机化检验 |
| Pagel's $\lambda$ | `"lambda"` | 连续性状 | $\lambda$ 极大似然估计、相对于 $\lambda = 0$ 的可选似然比检验（LR test） |
| Fritz–Purvis D | `"D"` | 二元性状 | D 统计量、随机零模型与布朗运动（BM）阈值模型检验 |
| Delta ($\Delta$) | `"Delta"` | 多状态分类性状 | $\Delta$ 估计值、MCMC 诊断指标与可选的随机化检验 |

你可以通过 `print(fit)`、`summary(fit)`、`as.data.frame(fit)` 查看结果，保存模拟值或似然剖面后，也可以用 `plot_signal(fit)` 绘制结果图（见下文绘图说明）。如果传入了当前方法不支持的参数，函数会直接报错。

---

## 分析流程

![fastphylosig 技术路线图](docs/technical-roadmap-zh.svg)

树和性状数据是所有分析的共同起点。

- **常规用法**：由 `fast_signal()` 自动在后台完成检查、匹配与预处理。
- **高级用法**：按相同顺序显式整理树和数据，生成计算上下文，再调用专属的 `fast_*` 函数。

> 中文路线图的 [draw.io 源文件](docs/technical-roadmap-zh.drawio) 和 [Mermaid 源文件](docs/technical-roadmap-zh.md) 存放在本仓库的 `docs/` 目录下。

---

## 高级用法：先整理树和数据，再逐步计算

当你需要保留清晰的数据处理记录，或者使用特定方法的专属参数时，可以分步执行：

```r
# 1. 检查树结构（只读）
checked <- check_tree(tree, signal = "K")

# 2. 规范化树表示层（处理边顺序与内部节点编号）
# 注意：如果方法要求有根树而当前树无根，必须由用户显式指定外群（outgroup）
tree_ready <- resolve_tree(tree, signal = "K")

# 3. 匹配树与性状数据（按树端点顺序自动重排）
matched <- match_tree_data(tree_ready, data = trait)

# 4. 建立预处理上下文快照
ctx <- prepare_tree(matched$tree)
```

整理完成后，直接调用对应的专业函数：

| 专业函数 | 适用性状 | 专属参数与功能 |
|---|---|---|
| `fast_k()` | 连续 | 自定义置换矩阵、分块计算、保存零分布 |
| `fast_lambda()` | 连续 | 似然剖面分析（`lambda_profile = TRUE`）、近似置信区间 |
| `fast_d()` | 二元 | 随机/Brownian 状态控制、零分布保存 |
| `fast_delta()` | 多元 | MCMC 链长与抽样控制（`mcmc_sim`、`thin`、`burn`）、链收敛诊断 |

```r
# K：适用于连续性状
fit_k <- fast_k(ctx, x = matched$data, test = TRUE, nsim = 999,
                return_sim = TRUE, progress = FALSE)

# Lambda：适用于连续性状
fit_lambda <- fast_lambda(ctx, trait, test = TRUE,
                          lambda_profile = TRUE, progress = FALSE)

# D：适用于二元性状
binary <- setNames(as.integer(trait > stats::median(trait)), tree$tip.label)
fit_d <- fast_d(ctx, binary, test = TRUE, nsim = 99,
                return_sim = TRUE, keep_null = TRUE, progress = FALSE)

# Delta：适用于多元、离散性状
categorical <- setNames(rep(c("a", "b", "c"), length.out = 30), tree$tip.label)
fit_delta <- fast_delta(ctx, categorical, test = TRUE, nsim = 99,
                        mcmc_sim = 10000, thin = 10, burn = 100,
                        return_sim = TRUE, progress = FALSE)
```

> **补充说明**：`fast_ace()` 是 Delta 内部使用的离散祖先状态最大似然拟合函数，供需要直接提取或检查祖先状态估计的用户使用。

---

## 速度表现

![四种系统发育信号方法的配对运行时间比](docs/release-0.3.0/benchmark_reference_four_methods_speedup.png)

在 2,000 个物种的配对基准中，参考实现与本包的中位耗时比为：K **135–151 倍**、λ **20,697–23,591 倍**、D **218–224 倍**。每个场景包含 10 对计时，不含树和数据准备时间。

Δ 的返回调用耗时比为 **23.9–33.6 倍**；由于部分拟合存在数值或收敛诊断，这一比值不表示两种实现具有相同精度或收敛程度。更多图表见[中文性能介绍](https://yinwen-ecology.github.io/fastphylosig/index-zh.html#speed)，完整设置与诊断见[英文基准材料](docs/release-0.3.0/README.md)。

---

## 系统发育信号异质性（`fast_hetero()`）

`fast_hetero()` 用于比较大型系统发育树内部不同支系间的信号强弱。它会在自动选取的互不重叠子树（或指定节点的子树）上计算系统发育信号（ K、$\lambda$、D 或 $\Delta$），并报告各支系的样本量、缺失值数量、状态组成以及全树覆盖率。

异质性分析的核心是**描述信号强弱在进化树上的空间分布格局**；异质性参考检验则作为可选的补充，与普通的整树信号检验分开报告。

以下另建一棵 256 物种的演示树，并定义三类性状，便于选出多个子树进行比较。

```r
set.seed(2)
tree <- ape::stree(256, type = "balanced")
tree$edge.length <- rep(1, nrow(tree$edge))
trait <- ape::rTraitCont(tree)
binary_trait <- ape::rTraitDisc(tree, model = "ER", k = 2, rate = 0.15)
categorical_trait <- ape::rTraitDisc(tree, model = "ER", k = 3, rate = 0.15)

# 1. 连续性状（以 K 为例）：自动选取 20~200 个 tips 的互不重叠子树
h <- fast_hetero(
  tree, trait, method = "K", min_tips = 20, max_tips = 200,
  test = TRUE, nsim = 999, overall_test = TRUE, progress = FALSE
)
h$clades                              # 各子树估计、样本量与覆盖状态
h$heterogeneity                       # 支系间离散度与异质性检验结果
h$overall                              # 整棵树的全局信号值及检验
h$selection$coverage                  # 选定子树对全树的覆盖率

# 可视化：在原树上为不同支系着色
plot(h, tree = tree, node_labels = TRUE)
plot(h, tree = tree, type = "fan")      # 扇形树展示

# 2. Pagel's lambda 异质性分析
hl <- fast_hetero(tree, trait, method = "lambda",
                  test = TRUE, nsim = 999, overall_test = TRUE, progress = FALSE)
hl$heterogeneity                       # 包含 SD(lambda)、共享 lambda 及异质性 P 值
plot(hl, tree = tree)

# 3. 二元性状 D 异质性分析（条件 ER 模型检验为实验性功能）
hd <- fast_hetero(tree, binary_trait, method = "D", min_tips = 50,
                  test = TRUE, nsim = 999, d_nsim = 1000,
                  overall_test = TRUE, progress = FALSE)
hd$heterogeneity                       # 包含 SD(D)、共用 ER 速率及异质性 P 值
plot(hd, tree = tree)

# 4. 多状态分类性状 Delta：推荐使用描述性流程（test = FALSE）
hdelta <- fast_hetero(tree, categorical_trait, method = "Delta",
                      overall_test = TRUE, overall_nsim = 999,
                      delta_control = list(model = "ER"), progress = FALSE)
hdelta$clades[, c("node", "estimate", "n_states", "log_ratio_overall")]
hdelta$clades$state_counts              # 保留原始状态的物种计数
hdelta$heterogeneity                   # SD(log Δ) 与局部/总体 RMS 偏差
plot(hdelta, tree = tree)

# 若需要显式调用实验性的 Delta 异质性检验：
hdelta_test <- fast_hetero(tree, categorical_trait, method = "Delta",
                           nodes = hdelta$selection$nodes, test = TRUE, nsim = 999,
                           delta_control = list(model = "ER"), progress = FALSE)
hdelta_test$heterogeneity$P_heterogeneity
```

<details>
<summary><strong>异质性检验的统计原理与解释边界</strong></summary>

1. **`P_BM`（K 的异质性检验）**：检验各子树 K 值的离散程度是否显著大于布朗运动（BM）下通常产生的离散度。该检验严格依赖 BM 前提，不能解释为“在任意演化模型下各支系信号均相同”的检验。绘图时的灰色区间为 BM 参考包络线，而非 K 的置信区间。
2. **$\lambda$ 与 D 的检验**：`P_heterogeneity` 是统一的输出字段。$\lambda$ 以“各冠群共享局部 $\lambda$ 参数的近似高斯模型”为参考；D 则以“共用转移速率的二状态 ER 模型”为参考，并以各有效类群的观测状态数为条件。两者均为参数自助法（parametric bootstrap）近似检验。`d_nsim` 控制 D 的随机/BM 标准化模拟次数，`nsim` 控制异质性模拟次数。
3. **$\Delta$ 的偏差检验**：$\Delta$ 在 `test = TRUE` 时自动计算总体估计，检验统计量为局部与总体的均方根偏差 `RMS(log(Δ_local / Δ_overall))`，而非简单计算 `SD(log Δ)`，也不是演化速率的异质性检验。该实验性检验通过在原树上拟合共同 ER/ARD 模型（根状态采用均匀混合，缺失物种在似然中积分）生成模拟数据，属于近似参考检验。如果任何一次模拟估计失败，P 值会返回 `NA` 并报告 `nsim_successful`、`status` 和 `message`，不会静默剔除失败模拟后继续计算。
4. **总体信号与局部信号的关系**：
   - `overall = TRUE` 仅计算整棵树的全局信号值；`overall_test = TRUE` 会同时运行全局信号检验（K 使用 tip 置换检验，$\lambda$ 使用似然比检验，D 显示随机关联与 BM 阈值模型两个 P 值）。
   - 两者默认关闭，避免在大型树上带来额外开销。
   - 两个支系的信号 P 值一个显著、一个不显著，并不等同于两者的信号存在显著差异。
   - D/Δ 的单一状态子树将被标记为“不可估计”（不会赋值为 0）；比较过小的子树或状态组成差异过大的子树时需要特别谨慎。
5. **绘图说明**：`plot()` 函数仅渲染已计算好的结果，不执行新的统计检验；传入 `show_summary = FALSE` 可隐藏图上的统计文本标注。完整示例见[中文案例](https://yinwen-ecology.github.io/fastphylosig/fast-hetero-case-zh.html)，详细方法见[方法、PVR 局部解释度及限制说明](docs/signal-heterogeneity-zh.md)与[中文完整指南](inst/HETERO_GUIDE_zh.md)。

</details>

---

## 树与数据格式要求

### 树结构规范与 `resolve_tree()`

`check_tree()` 仅做只读检查。`resolve_tree()` 仅在数据副本上处理边顺序和内部节点编号等表示层格式：

- **节点编号规则**：保持末端节点 ID 与 `tip.label` 不变，将根节点编号为 `n_tip + 1`；按子树最小后代末端标签的 UTF-8 字节顺序排列子节点，并按前序规则对其余内部节点编号；边按规范后序排列，枝长严格跟随对应的生物学子边。
- **枝长要求**：
  - **K 与 D**：要求枝长有限且严格为正（$> 0$）。D 支持有根多分叉树，但不支持单子节点（包括数据匹配或剔除 NA 后暴露出单子节点的情况）。
  - **$\lambda$**：仅在变换后的末端方差仍有效时允许内部零枝长，末端零枝长会被直接拒绝。
  - **$\Delta$ 与 `fast_ace()`**：严格要求常规有根、完全二叉树，且枝长必须严格为正。
- **根与外群**：根节点位置由用户输入决定，更改根的位置可能会改变最终的统计推断结果。
- **测量误差**：当前版本不支持测量误差输入，参数 `se` 必须保持为 `NULL`。

### 性状数据对齐与缺失值

- **物种标识**：建议为性状数据提供唯一的向量 `names` 或矩阵/数据框的 `row.names`。输入顺序无需与树完全一致，程序会在匹配后自动重排为树的端点顺序。
- **未命名数据的风险**：只有当未命名的向量或矩阵行顺序与 `tree$tip.label` 完全严格一致时，才可以省略名称；如果顺序不一致且未命名，程序无法自动发现错配，会导致性状被分配给错误的物种。
- **缺失值处理**：缺失值会在物种匹配完成后按各性状独立剔除。

---

## 绘图与可视化

| K 的置换零分布（示例） | λ 的似然剖面（示例） |
|---|---|
| ![K 的置换零分布（示例）](docs/release-0.3.0/plot_signal_K_single_trait.png) | ![λ 的似然剖面（示例）](docs/release-0.3.0/plot_signal_lambda_profile.png) |

```r
# 绘制信号结果与零分布或似然剖面
plot_signal(fit_k)
plot_signal(fit_lambda)
plot_signal(fit_d)
plot_signal(fit_delta)
```

> **注意**：绘制零分布图的前提是计算时已保存模拟值：K 和 Delta 需要在计算时指定 `return_sim = TRUE`；D 则需要同时指定 `return_sim = TRUE` 与 `keep_null = TRUE`。λ 的似然剖面图使用 `lambda_profile = TRUE`。

---

<details>
<summary><strong>推断细节与 MCMC 诊断</strong></summary>

1. **Delta 的后验均值与 MCMC 诊断**：
   - Delta 的点估计值为后验均值 `mean(beta / alpha)`，提议分布包含对数正态的 Hastings 密度比校正。
   - 结果中包含 `ESS_Delta`（有效样本量）、`split_Rhat_Delta`（分裂链 R-hat）和 `MCSE_Delta`（比值序列的均值标准误）。当 `test = TRUE` 时还会报告置换检验的 `MCSE_P`。
   - 默认参数为：`test = FALSE`、`nsim = 1000`、`return_sim = test`、`mcmc_sim = 10000`、`thin = 10`、`burn = 100`、`model = "ARD"`。因此默认不计算置换 P 值，也不保存模拟零分布。
   - 物种匹配并剔除缺失值后，若某个性状仅剩 2~19 个有效物种，程序会发出批量警告；若有效物种少于 2 个，该性状将返回 `insufficient_data` 结果，若全局共享物种少于 2 个则直接整次报错。
2. **D 的零模型与标准误**：
   - D 同时报告随机零模型与布朗运动零模型的 P 值，以及基于有限次成功模拟计算的 Monte Carlo 标准误 `MCSE_P_random` 与 `MCSE_P_Brownian`（`test = FALSE` 时返回 `NA`）。
   - 只有在 `return_sim = TRUE` 且 `keep_null = TRUE` 同时满足时，D 才会保存两个零分布向量。
3. **K 的零分布保存**：
   - K 支持显式置换矩阵与分块计算；可用 `return_sim = TRUE` 或 `keep_null = TRUE` 保存零分布。

</details>

---

<details>
<summary><strong>性能优化与可重复性</strong></summary>

### 性能建议

1. **批量性状输入**：多个性状请合并为矩阵或数据框一次性传入，避免在 R 语言中使用 `for` 循环逐列计算。
2. **复用预处理上下文（`ctx`）**：对同一棵树进行多次分析时，使用 `prepare_tree()` 创建上下文并复用。该上下文是创建时刻的受保护快照，包含树结构、边矩阵、枝长和节点数的严格校验合同。支持通过 `saveRDS()` 和 `readRDS()` 跨 R 会话持久化复用。
3. **控制零分布保存**：零分布主要用于绘图，常规分析无需开启 `return_sim = TRUE`，以降低内存占用。
4. **静默运行**：设置 `progress = FALSE` 可关闭进度条和匹配提示（统计警告和错误依然会正常输出），且计算结果中依然会保留完整的计时元数据（timing）。

<details>
<summary><strong>多线程与随机可重复性</strong></summary>

- **多核加速（`ncores`）**：K 随机化与二叉树 D 的 C++ 内核使用 OpenMP 线程，Delta 置换分析使用 R 工作进程，均通过 `ncores` 控制。
- **结果复现规则**：并行随机计算的完全复现，必须**同时固定随机种子（`set.seed`）、显式模拟输入（如有）以及相同的 `ncores` 数值**。更改 Delta 的进程数会改变并行随机数流的分配，使生成的零分布及 P/MCSE 可能发生变化。
- **OpenMP 内部线程控制**：Delta 与 `fast_hetero(method = "Delta")` 中的大型 ACE 拟合默认最多使用 2 个 OpenMP 线程并行计算未缓存枝长的指数向量。可通过 `options(fastphylosig.ace_threads = 1L)` 强制单线程运行，或设为 `0L` 恢复自动选择。调整此选项不会改变随机数流的分配。

</details>

</details>
