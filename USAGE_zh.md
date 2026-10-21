# fastphylosig 中文使用指南

`fastphylosig` 用于计算连续、二元和分类性状的系统发育信号。多数分析只需
使用 `fast_signal()`；需要逐步检查或调整方法参数时，再使用后面的专业函数。

当前发布版为 **0.3.5**。它修正了 ARD 祖先概率和
离散 ACE 速率标准误。使用过 ARD Δ、ARD 全树祖先概率或 ACE 标准误的分析，
请参阅[补丁说明与重算建议](docs/release-0.3.5/RELEASE_NOTES_zh.md)。

## 系统发育信号异质性

使用 `fast_hetero()` 比较支系信号。
它在自动选取的互不重叠子树
或指定节点的子树上计算 K／lambda／二元 D／多状态 Δ，并报告各类群的样本量、
缺失值、状态组成及全树覆盖率。重点是描述信号强弱在树上的分布。
异质性参考检验作为可选补充，与普通信号检验分别报告。

K 使用 BM 参考；λ 使用共同局部冠群参数的近似参考。D 的条件 ER 异质性检验为实验性；Δ 推荐局部/总体描述，异质性检验仅作显式实验调用。原有普通信号检验不受此范围限制。

```r
h <- fast_hetero(
  tree, trait, method = "K", min_tips = 20, max_tips = 200,
  test = TRUE, nsim = 999, overall_test = TRUE, progress = FALSE
)
h$clades
h$heterogeneity
h$overall                              # 整棵树的信号值和信号检验
h$selection$coverage
plot(h)
plot(h, tree = tree, node_labels = TRUE) # 在原树上着色
plot(h, tree = tree, type = "fan")      # 扇形树

hl <- fast_hetero(tree, trait, method = "lambda",
                 test = TRUE, nsim = 999, overall_test = TRUE, progress = FALSE)
hl$heterogeneity                       # SD(lambda)、共享 lambda、异质性 P 值
plot(hl, tree = tree)

# D 条件 ER 异质性检验：实验性；原有普通 D 信号检验不受此限制。
hd <- fast_hetero(tree, binary_trait, method = "D", min_tips = 50,
                 test = TRUE, nsim = 999, d_nsim = 1000,
                 overall_test = TRUE, progress = FALSE)
hd$heterogeneity                       # SD(D)、共用 ER 速率、异质性 P 值
plot(hd, tree = tree)

# 多状态性状：总体与局部信号、局部/总体偏差、树图。
hdelta <- fast_hetero(tree, categorical_trait, method = "Delta",
                     overall_test = TRUE, overall_nsim = 999,
                     delta_control = list(model = "ER"), progress = FALSE)
hdelta$clades[, c("node", "estimate", "n_states", "log_ratio_overall")]
hdelta$clades$state_counts              # 保留原状态标签
hdelta$heterogeneity                   # SD(log Δ)、局部/总体 RMS 偏差
plot(hdelta, tree = tree)

# 需要各子树的普通信号检验时，加 clade_test = TRUE。
# 可选补充：共同 Mk 模型下的局部/总体 Δ 偏差参考检验。
# 以下为显式实验性 Δ 异质性检验调用。
hdelta_test <- fast_hetero(tree, categorical_trait, method = "Delta",
                     nodes = hdelta$selection$nodes, test = TRUE, nsim = 999,
                     delta_control = list(model = "ER"), progress = FALSE)
hdelta_test$heterogeneity$P_heterogeneity
```

`P_BM` 检验 K 的离散程度是否大于各子树在 BM 下通常产生的离散程度。
它有明确的 BM 前提；不能解释为任意演化模型下的“信号相同”检验。
灰色区间是 BM 参考包络，不是 K 的置信区间。
树图的颜色表示类群内部的估计；K／Δ 默认采用对数色阶，lambda／D 采用原值色阶。
D 默认反转调色板，使较强信号的颜色方向与其他指标一致。
灰色表示未选择或无法估计的位置，嵌套选择显示更内层类群的估计。
`P_heterogeneity` 是统一的异质性检验字段，K 的 `P_BM` 继续保留。
lambda 以共同局部冠群 lambda 参数的近似高斯模型为参考；D 以共用转移速率的二状态 ER
模型为参考，并以每个有效类群的观测状态数量为条件。二者为参数自助法
近似检验；不能解释为所有演化模型下通用的“信号相同”检验。
`d_nsim` 控制 D 的随机／BM 标准化次数，`nsim` 控制异质性模拟次数。
模拟保留各类群的状态比例；对应的 D 标准化均值在观测与模拟之间复用。
图上分别标注整棵树的信号值与检验、类群间离散度与异质性检验。
`overall_test = TRUE` 会计算整棵树的信号及检验：K 用 tip 随机化，lambda
用相对于 lambda = 0 的似然比检验；D 显示随机关联与 BM 阈值模型两个
信号 P 值。它们不等同于异质性检验。
仅需总体信号值时用 `overall = TRUE`；两者默认关闭，避免大型树额外计算。
Δ 的 `test = TRUE` 自动计算总体估计；检验统计量为
`RMS(log(Δ_local / Δ_overall))`，不是 SD(log Δ)，也不是速率异质性检验。
共同 ER／ARD 模型在原树上拟合，根状态采用均匀混合，缺失 tips 在似然中积分。
模拟整棵树后，按原缺失模式重新估计各子树和总体 Δ。它是本包新增的近似参考
检验，尚没有全面的误报率／功效校准。若任何模拟估计失败，P 值为 NA，
并报告 `nsim_successful`、`status`、`message`；不删除失败模拟后继续算 P。
检查观测 MCMC 的 ESS、split R-hat、MCSE_Delta，以及模拟的
`nsim_MCMC_warning`。Δ 的普通总体和局部信号检验仍是 tip 随机化。
单状态子树标为不可估计，不赋值为 0；小子树和不同状态组成的比较需要谨慎。
不同类群的信号 P 值一显著一不显著，不等于两者信号显著不同。
绘图只显示已计算的结果，不运行检验。`show_summary = FALSE` 可隐藏统计标注。
详见[方法、PVR 局部解释度及限制](docs/signal-heterogeneity-zh.md)。

## 安装

```r
# 安装 0.3.5 源码包（需要 C++ 工具链及 ape、Rcpp、RcppArmadillo）
install.packages("https://github.com/yinwen-ecology/fastphylosig/releases/download/v0.3.5/fastphylosig_0.3.5.tar.gz", repos = NULL, type = "source")
library(fastphylosig)
```

也可用 `remotes::install_github("yinwen-ecology/fastphylosig@v0.3.5")` 安装固定标签。
保留的 50K 验证覆盖树准备、子集缓存、持久化、fallback 分组及无置换 K，
不代表所有方法与模拟场景均已在 50K 完成验证。

从源码安装需要 C++ 编译工具，OpenMP 为可选加速项。

## 最快用法：一个命令完成分析

准备好 `phylo` 树和带物种名称的性状数据后，选择方法并运行：

```r
set.seed(1)
tree <- ape::rtree(30)
trait <- setNames(stats::rnorm(30), tree$tip.label)

fit <- fast_signal(
  tree, data = trait, method = "K", progress = FALSE
)
summary(fit)
```

对于一般分析，这几行就够了。`fast_signal()` 会依次检查树、匹配物种、按性状
处理缺失值，然后执行所选的系统发育信号计算。

根据研究问题和性状类型选择 `method`：

| 目标 | `method` | 输入 | 主要结果 |
|---|---|---|---|
| Blomberg's K | `"K"` | 连续性状 | K，可选随机化检验 |
| Pagel's lambda | `"lambda"` | 连续性状 | lambda 似然和 LR 检验 |
| Fritz--Purvis D | `"D"` | 二元性状 | 随机零模型和 Brownian 零模型比较 |
| Delta | `"Delta"` | 分类性状 | Delta、MCMC 诊断和可选随机化检验 |

函数不会根据性状值猜测方法。可用 `print(fit)`、`summary(fit)`、
`as.data.frame(fit)` 或 `plot_signal(fit)` 查看结果；不适用于当前方法的参数会
直接报错。

## 分析流程

![fastphylosig 技术路线图](docs/technical-roadmap-zh.svg)

树和性状数据是两种用法的共同起点。一般用法由 `fast_signal()` 自动完成图中的
准备步骤；高级用法按相同顺序显式整理树和数据，再调用对应的 `fast_*` 函数。
中文路线图的 [draw.io 源文件](docs/technical-roadmap-zh.drawio) 和
[Mermaid 源文件](docs/technical-roadmap-zh.md) 与本指南放在同一仓库中。

## 高级用法：先整理树和数据，再计算

需要保留清晰的整理记录，或使用方法专属参数时，可逐步运行：

```r
checked <- check_tree(tree, signal = "K")
checked

# 只处理安全的表示层问题；如果方法需要有根树而当前树无根，
# 必须由用户明确提供 outgroup。
tree_ready <- resolve_tree(tree, signal = "K")
matched <- match_tree_data(tree_ready, data = trait)
ctx <- prepare_tree(matched$tree)
```

完成整理后，直接调用对应的专业函数：

| 函数 | 适用性状 | 专属能力 |
|---|---|---|
| `fast_k()` | 连续 | 自定义置换、分块计算和零分布保存 |
| `fast_lambda()` | 连续 | 似然剖面和近似区间 |
| `fast_d()` | 二元 | 随机/Brownian 状态控制和零分布保存 |
| `fast_delta()` | 分类 | MCMC、随机化检验和链诊断 |

```r
fit_k <- fast_k(ctx, x = matched$data, test = TRUE, nsim = 999,
                return_sim = TRUE, progress = FALSE)
fit_lambda <- fast_lambda(ctx, trait, test = TRUE,
                          lambda_profile = TRUE, progress = FALSE)
binary <- setNames(as.integer(trait > stats::median(trait)), tree$tip.label)
fit_d <- fast_d(ctx, binary, test = TRUE, nsim = 99,
                return_sim = TRUE, keep_null = TRUE, progress = FALSE)
categorical <- setNames(rep(c("a", "b", "c"), length.out = 30),
                        tree$tip.label)
fit_delta <- fast_delta(ctx, categorical, test = TRUE,
                        mcmc_sim = 10000, thin = 10, burn = 100,
                        return_sim = TRUE, progress = FALSE)
```

`fast_ace()` 是 Delta 使用的离散祖先状态最大似然函数，主要面向需要直接检查
祖先状态估计的用户。

## 重复分析

```r
ctx <- prepare_tree(tree)
X <- cbind(trait_1 = trait, trait_2 = trait + stats::rnorm(length(trait)))
fits <- fast_signal(ctx, data = X, method = "K", verbose = FALSE)
cache_info(ctx)
```

prepared context 是创建时的快照：之后修改原始 `phylo` 对象不会更新它。
每次从外部复用时，V2 schema 和精确的受保护快照都会在使用缓存或数值核心
之前检查 context 内树的物种名称、边矩阵、枝长和 `Nnode`；这些字段被修改后
会拒绝复用。`node.label` 和其他不参与计算的属性不在该完整性合同内；内部
节点重编号时，`node.label` 可能随之重排。V2
context 可以通过 `saveRDS()` 保存，并在另一个 R 会话中用 `readRDS()` 后复用；
旧版或无法识别的 context 必须重新运行 `prepare_tree()`。
生产 K、lambda、D 和 Delta 路径不会仅因准备树而建立稠密协方差矩阵。

## 树与数据细节

`check_tree()` 只读。`resolve_tree()` 只在副本上处理边顺序和内部节点编号等
表示层问题。V2 保持末端节点 ID 和 `tip.label` 不变，把根编号为 `n_tip + 1`，
按各子树最小后代末端标签的 UTF-8 精确字节顺序排列子节点，并以前序规则编号
其余内部节点；边按规范后序排列，枝长始终跟随对应的生物学子边。对于已测试
的合法内部节点重编号、边行顺序和 locale，规范计算表示、指纹和估计结果保持
一致。它不会猜测生物学根、选择外群、补造或抖动枝长、解析生物学
多分叉、超度量化树或静默删除物种。

- K 和 D 要求枝长有限且严格为正。D 可通过兼容路径处理有根多分叉，
  但不支持单子节点，包括匹配或逐性状剔除 NA 后暴露的单子节点。
- lambda 仅在变换后的末端方差仍有效时允许零内部枝，零末端枝会被拒绝。
- Delta 和 `fast_ace()` 要求常规有根、完全二叉且枝长严格为正。
- 根位置由用户提供，改变根可能改变统计结果。
- 当前不支持测量误差输入，即 `se` 必须为 `NULL`。

物种标识应为唯一的向量 names 或矩阵/data.frame row names。输入顺序无需与树
一致，匹配后会重排为树端点顺序。只有当未命名向量的元素或未命名表的行
本来就与 `tree$tip.label` 顺序一致时，才可以省略名称。缺失值在物种匹配后按性状处理。
如果省略名称时数据顺序实际不正确，程序无法自动识别这种位置错配，可能把性状
对应到错误的物种。

## 推断与诊断

0.3.2.9000 将 Delta 修正为后验 `mean(beta/alpha)`，并加入对数正态提议的
Hastings 密度比。估计、检验与 RNG 轨迹会与旧版不同。新增 `ESS_Delta` 和
`split_Rhat_Delta`，`MCSE_Delta` 直接评价比值序列的均值误差。
`fast_hetero()` 的切割、总体检验、异质性检验及其校准边界见
[完整指南](inst/HETERO_GUIDE_zh.md)。

K 随机化检验支持显式置换矩阵、分块计算和可选的零分布保存。D 分别报告
随机零模型和 Brownian 零模型的 P 值与 Bernoulli Monte Carlo 标准误
`MCSE_P_random`/`MCSE_P_Brownian`，分母为各自有限且成功的模拟数；
`test = FALSE` 时这些 MCSE 为 `NA`。D 只在 `return_sim = TRUE` 且
`keep_null = TRUE` 同时成立时保存两个 null 向量；两者的默认值分别为 `test` 和
`return_sim`。

Delta 默认包括 `test = FALSE`、`nsim = 1000`、`return_sim = test`、
`mcmc_sim = 10000`、`thin = 10`、`burn = 100` 和 `model = "ARD"`。因此，默认不计算
permutation P 值，也不保存 null。`test = TRUE` 时，`MCSE_P` 使用有限且成功的
置换次数；Delta 零分布的保存由 `return_sim = TRUE` 单独控制，Delta 没有
`keep_null` 参数。

全局共享物种少于 2 个时，Delta 整次调用报错。物种匹配后，某个性状因 NA 仅保留
少于 2 个物种时，返回该性状的 `insufficient_data` 结果；完成计算但只保留
2--19 个物种的性状会汇总为一次批量警告。解释结果时应联合检查 ESS、split R-hat、
`MCSE_Delta` 以及在存在时的 `MCSE_P`。警告不会自动重跑或延长链。

## 支系间信号与系统发育保守性

`fast_hetero()` 默认 `test = FALSE`。推荐对 Δ 使用
`test = FALSE, overall = TRUE`，描述各支系及整树的信号；`overall_test = TRUE`
仍可计算原有整树信号检验。Δ 异质性检验已移出推荐流程，显式研究调用保留。

局部信号不同可以表述为“在所分析尺度上，不同支系表现出不同程度的近缘相似
格局（系统发育保守性）”，并报告指标、样本量和状态组成。这个格局描述不要求
先获得异质性 P；它自身不识别稳定选择、生态位保守过程或低演化速率。

K 的可选 P 以 BM 为参考；λ 的可选 P 近似检验各冠群共享局部 λ，不能解释为
整树共同 λ 的直接检验。不同冠群年龄及缺失后根节点移动会改变 λ 的尺度。
D 的可选异质性 P 是给定状态计数的共同对称 ER 实验性检验。
推荐切割、零假设、结果表和可执行示例见 [完整指南](inst/HETERO_GUIDE_zh.md)。

## 可重复性与性能

1. 多个性状使用矩阵一次传入，避免在 R 中逐列循环。
2. 重复分析时复用 `prepare_tree()` context。
3. 只在绘图或审计需要时保存 null 分布。
4. K 随机化、二叉树 D 和足够大的 Delta permutation 任务可使用 `ncores`。
5. `progress = FALSE` 只关闭进度阶段和物种匹配消息，统计警告和错误仍会正常发出；
   结果仍保留 timing 元数据。

要复现并行随机计算，需同时固定随机种子、显式模拟输入（如有）和 `ncores`。
改变 Delta 的进程数会重新分配并行随机数流，不能保证零分布及其 P/MCSE 逐次一致。
显式置换/状态矩阵可固定对应的随机化步骤，但 Delta MCMC 本身仍是随机的。

Delta 与 `fast_hetero(method = "Delta")` 的大型 ACE 拟合默认最多用两个 OpenMP
线程批量计算未缓存枝长的指数向量，并与完全相同枝长的精确缓存配合使用。
是否启用并行取决于未缓存的计算量；剪枝累积、优化器和 MCMC 的执行顺序保持不变。
`options(fastphylosig.ace_threads = 1L)` 可强制单线程；`0L` 恢复自动选择，正整数
可请求更多线程，实际受 OpenMP 运行时上限约束。小子树、无 OpenMP 的构建和 PSOCK
置换子进程内的拟合使用单线程。增加线程会占用更多 CPU；这个选项不重新分配随机数流。

## 绘图

```r
plot_signal(fit_k)
plot_signal(fit_lambda)
plot_signal(fit_d)
plot_signal(fit_delta)
```

绘制零分布要求计算时保存模拟值：K 和 Delta 使用 `return_sim = TRUE`；D 需要
`return_sim = TRUE` 且 `keep_null = TRUE`。

## 版本兼容性

包声明支持 R 4.1.0 及以上版本。目前已在 Windows 的 R 4.6.1/Rtools 4.5 环境完成
本地源码包检查；这不代表已经完成全部跨平台验证。
