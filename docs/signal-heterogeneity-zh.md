# 树内系统发育信号异质性

> 当前推荐范围（0.3.4）：K 的 BM 参考、λ 的共同局部参数近似参考，
> 以及二元 D／多状态 Δ 的局部与总体描述。D 条件 ER 检验保持实验性；
> Δ 异质性检验已移出推荐流程。最新工作流见
> [安装指南](../inst/HETERO_GUIDE_zh.md)，数学尺度核查见
> [本轮范围收敛记录](signal-heterogeneity/scope-scale-20261016/README_zh.md)。

系统发育信号可以在树内变化。全树统计量把不同位置、不同尺度的性状结构
汇总成一个数，因此全树有信号不代表每个类群都有同样强的信号；全树没有
显著信号也不能排除局部结构。要分析这种变化，需要先明确比较的尺度。

## 三类互补问题

| 问题 | 方法 | 结果的含义 |
|---|---|---|
| 哪些物种附近出现近缘相似性聚集或局部异常？ | LIPA／局部 Moran's I | 物种与其系统发育邻域的局部关联；依赖邻域和权重定义 |
| 哪些类群内部的信号更强？ | 对预定子树计算 K／lambda／D／Δ | 每个类群内部的结构；局部值与全树值分别估计 |
| 信号集中在深层还是浅层分化？ | 系统发育相关图、PSR 曲线 | 距离或谱尺度分解；尺度差异并不直接定位到某个类群 |

这些方向已有方法基础。`phylosignal::lipaMoran()` 提供 tip 级 LIPA，
`phylosignal::phyloSignalINT()` 提供内部节点的子树扫描。后者逐节点计算
信号和置换概率，但各节点常常互相嵌套，扫描本身不是一个独立的异质性检验。

局部 Moran's I 常写为

\[
I_i=z_i\sum_j w_{ij}z_j,
\]

其中 z 是中心化、标准化的性状，w 描述系统发育邻近程度。正值表示同方向
的局部聚集，负值表示物种与邻域方向相反；高–高、低–低等标签还需要结合
性状与邻域均值的符号。它不能直接当作局部 K 或局部方差解释比例。常规
全距离权重矩阵需要二次内存；大型树的稀疏邻域版本需要明确邻域尺度与
稀疏化带来的方法变化。此次没有实现 LIPA。

## 开发版中的功能

本功能属于开发版；安装当前源码后使用。审计使用的本地库及冻结源包路径见
相应版本记录，不要把历史 `.rlib_heterogeneity` 当作当前版本。

函数入口统一为 `fast_hetero()`，此前开发阶段的长名称已移除。

它的估计流程是“按树选出冠群子树 → 调用已有信号内核 → 汇总局部信号”。
K 和 lambda 直接使用 `fast_k()`、`fast_lambda()` 的 C++ 树内核；D 与
`fast_d()` 共享预先计算的对比权重及遍历。内部省去反复调用公开 R 函数所需的
整套校验、匹配和树准备，并按缺失模式批量计算。D 的随机／BM 标准化保留
自己的少数状态编码和模拟顺序，因此同一随机种子下有限模拟的 D 不必与
单独调用 `fast_d()` 完全相同。两者使用同一 D 定义。

`test = TRUE` 还需要在明确的共同参考模型下拟合、模拟和汇总离散度。
“某个类群信号显著、另一个不显著”不能代替异质性检验。lambda 复用似然工作区；
D 复用条件概率、转移概率和抽样工作区，避免每次模拟重新准备。
普通 `fast_ace()` 的 Mk 似然未对状态数量取条件，不能直接代替 D 的条件模型。
目前本接口支持 K、lambda、二元 D 和多状态 Delta。Delta 复用 `fast_ace()`
的树／模型工作区、似然内核以及 `fast_delta()` 的熵、MCMC 和诊断内核。
扫描描述信号异质性，不需要推断速率异质性。局部与总体的普通信号检验，
以及共同模型下的信号偏差参考检验，是不同结果。

分类性状还报告 `n_states`、带原状态标签的 `state_counts` 和
`state_frequencies`。总体信号直接计算；D 的偏差是 `difference_overall`，
Δ 的偏差是 `log_ratio_overall`，不对可为负值的 D 取对数。
`dispersion` 是子树间 SD(D) 或 SD(log Δ)；`local_global_deviation` 是
相对总体的偏差均方根。它们描述观测异质性，不能单独解释为显著性检验。

Δ 推荐使用 `test = FALSE, overall = TRUE` 保留局部、总体及其对比。
显式 `test = TRUE` 的历史实验模式比较局部与总体 Δ 的对数偏差均方根，
以原树上拟合的共同 ER／ARD Mk 模型、均匀根状态为参考。整树模拟保留
局部／总体依赖，复用 C++ 转移概率后，用相同设置重估总体和所有局部 Δ。
该模式已移出推荐流程；详细模型保存在
`heterogeneity$reference_model`。单状态子树、失败模拟和 MCMC 诊断明确报告，
不自动重采样或删掉失败值。全面的误报率与功效验证仍未完成。

```r
library(fastphylosig)
set.seed(42)

h <- fast_hetero(
  tree, data = traits, method = "K",
  min_tips = 20, max_tips = 200,
  test = TRUE, nsim = 999, simulation_chunk = 32,
  keep_null = FALSE, overall_test = TRUE, progress = FALSE
)

h$clades           # node、trait、estimate、n_species、status，以及 BM 参考包络
h$heterogeneity    # dispersion、P_BM、MCSE、有效类群数、模拟成功数
h$overall          # 整棵树的信号值、信号检验 P 值、样本量及状态
h$selection        # 类群物种名单、覆盖率、未覆盖物种、重叠状态
plot(h, trait = colnames(traits)[1])
plot(h, tree = tree, trait = colnames(traits)[1], node_labels = TRUE)
plot(h, tree = tree, trait = colnames(traits)[1], type = "fan")
```

如果研究设计已经确定了要比较的类群，可传入 `nodes = c(...)`。节点编号
指向传入的原始树；准备过的对象则使用 `prepared$tree` 的节点编号。函数
内部对树进行规范化不会改变结果表中报告的原始节点编号。

自动模式从根向末端寻找大小介于 min_tips 与 max_tips 之间的最大子树，
一旦选中便不再选其后代。选择仅依赖树，而不依赖性状值。极不平衡树可能
留下大量不足 min_tips 的末端残余，所以一定检查覆盖率。这里的 coverage
是树的拓扑覆盖率；各性状实际可用样本量在 n_species 中，不能与其混淆。

每个子树从冠群根开始计算，不包含连接该类群与其祖先的 stem 分支。
不同类群整体均值不同，不能单独说明各类群内部的信号强度不同。

lambda 也可以运行异质性检验：

```r
hl <- fast_hetero(
  tree, traits, method = "lambda", min_tips = 20, max_tips = 200,
  test = TRUE, nsim = 999, overall_test = TRUE, progress = FALSE
)
plot(hl, tree = tree, trait = colnames(traits)[1])
```

设置 `test = TRUE` 后，先拟合所有有效类群共用的 lambda，每个类群单独
剖面化均值和扩散方差，再参数自助模拟 SD(lambda)。`P_heterogeneity` 报告
上尾概率，`lambda_common` 报告共用参数，`null_mean`／`null_upper` 给出
参考分布摘要。它是共享 lambda 模型下的近似检验。
lambda 可能达到边界、在小类群中估计不稳定。需要查看某个类群的参数不确定性时，可单独调用
`fast_lambda(subtree, ..., lambda_profile = TRUE)`。显式选择的嵌套节点允许
描述性扫描，但异质性检验只接受互不重叠的类群。

二元性状优先描述 D；另行事先选择的条件 ER 检验为实验性：

```r
hd <- fast_hetero(tree, binary_trait, method = "D", min_tips = 51,
  max_tips = 200, test = FALSE, d_nsim = 1000,
  overall_test = TRUE, progress = FALSE)
hd$clades
hd$heterogeneity
plot(hd, tree = tree)
```

D 的离散度为 SD(D)，参考为共用速率 `q_common` 的对称二状态 ER Mk 模型，
根状态概率为 1/2、1/2。拟合与模拟均以每个有效类群的观测状态数量为条件，
通过多项式剪枝计算条件概率，再向下抽样；无需反复拒绝稀有状态组合。
单一状态的观测类群没有 D，明确列为无法估计。模拟保留各类群的状态比例，
对应的随机／BM 标准化均值由 `d_nsim` 次模拟估计，并在观测和异质性模拟
之间复用。因此检验也以这组有限模拟均值为条件；应提高 `d_nsim` 并检查
结果对标准化误差的敏感性。
`nsim` 控制异质性模拟次数，`MCSE_P_heterogeneity` 仅表示这一层的 Monte
Carlo 误差，不包括拟合参数和标准化曲线的不确定性。
该检验是本包新增的模型参考扩展，并非 Fritz–Purvis 原始 D 检验。
D 较低通常代表更强的聚集，可能小于 0 或大于 1；采用原值色阶。
需要分叉且枝长严格为正的冠群；小类群、极端状态比例可能造成不稳定。
条件速率下限为 1e-8／有效冠群平均枝长，无穷速率对应独立状态。
状态数量的概率计算比通常的树剪枝更费时，因此 D 检验限制每个有效冠群
最多 1000 个观测物种；默认自动分区上限为 200，不受此限制。

`overall = TRUE` 在整棵传入树上计算总体信号，包含未被选入局部类群的物种，
每个性状分别去除缺失值。它不是局部估计的平均。`overall_test = TRUE` 同时
开启总体信号及其检验：K 的 tip 随机化次数由 `overall_nsim` 控制，默认与
`nsim` 相同；lambda 使用已有的 lambda = 0 似然比检验（卡方、自由度 1）。
这两个选项默认关闭，以保留纯局部扫描的性能。总体模型无法估计时，
`overall$status` 和 `message` 报告原因，不丢弃已有的局部结果。

图上分别显示总体信号及其 P 值、局部离散度及其异质性 P 值。
D 的总体信号检验分别报告随机关联与 BM 阈值模型两个概率，采用 `d_nsim`
次数；它们与条件 ER 异质性检验的假设不同。
未开启或无法完成的检验会明确标注；绘图不重算分析、不产生随机数。
使用 `show_summary = FALSE` 隐藏统计标注。

树图给选中类群的冠群分支和根节点着色，表示整个类群的局部 K／lambda／D。
K 默认用对数色阶，lambda 和 D 用原值色阶；图例始终标示未变换的估计值。
灰色表示未选择或无法估计的区域。嵌套节点以更内层类群为准，
包括内层无法估计时显示灰色；类群 stem 保留外层类群的颜色。
树图不把一个类群估计解释成每条分支各自的信号，也不把颜色当作显著性。
使用 `color_scale = "raw"` 可改变 K 的色阶，`show.tip.label = TRUE` 可显示
物种名。较大树默认隐藏物种名，可指定 `type = "fan"`。
传入计算时的原树或准备对象，绘图会核对节点的后代物种，并通过 C++
迭代遍历生成绘图节点映射，处理原始节点编号和边顺序的变化。

## K 的离散度与 BM 参考检验

函数使用未加权统计量

\[
H=\operatorname{SD}_c\{\log K_c\}.
\]

因此每个类群贡献相等，而不是每个物种贡献相等。只有达到 min_tips 且 K
为有效正值的类群进入 H；排除数会明确报告。至少需要两个有效类群。
样本不足、恒定性状、非有限输入或局部树问题会保留为独立状态。

`test = TRUE` 时，在每个实际用于计算的子树上模拟单位扩散方差的 BM，
保留树形、枝长及该性状的缺失模式，重新计算 K 与 H。K 对性状的平移与
尺度变化不敏感，因此无需估计祖先均值或共同扩散方差。在互不重叠类群上
独立模拟具有正确的联合 K 分布：全树 BM 中共享的祖先增量在一个类群内
只是相同的常数，对该类群的 K 没有影响。嵌套类群不能用这个独立模拟论证。

参考概率为

\[
P_{BM}=\frac{1+\#\{H_b\ge H_{obs}\}}{B+1}.
\]

这个上尾检验回答“观察到的 K 离散程度是否超过 BM 在这些子树上通常产生
的离散程度”。它是此功能的明确操作性定义，不能称为所有模型通用的异质性
检验。全树统一的 OU 或统一测量误差过程，也可能因为在不同树形中的表现
而拒绝 BM 参考；显著结果不能直接证明不同类群存在不同演化机制。

BM_mean、BM_lower／BM_upper 描述参考模拟的均值与 95% 中央包络。
它们不是观察 K 的置信区间。heterogeneity$BM_upper 是 H 的单侧 95%
参考分位数。概率采用加一修正，最低可报告值为 1/(B+1)。MCSE_P_BM 是
二项模拟误差的插入式近似，极端计数时会为零，不能据此认为尾概率没有
不确定性。出现非有限模拟时会标记 simulation_failure 并不报告 P_BM。

不能从“类群 A 的信号检验显著，类群 B 不显著”推出两个类群信号显著不同。
它们的样本量和检验功效也不同。按观察结果挑选最高／最低 K 再进行参考
检验同样存在选择偏差；确认性分析应预先指定类群和尺度。多个性状的 P_BM
可用 `p.adjust(..., method = "BH")` 进行探索性校正；H 随分区尺度变化，
不同 min_tips／max_tips 的结果也需要作为敏感性分析明确报告。

## PVR 中哪个类群“解释得更好”

PVR 将性状回归到选定的系统发育特征向量上。全树 R² 受所选特征向量、
模型复杂度和评估方式影响，不等同于 K，也不能解释为系统发育的因果贡献。
纳入过多特征向量可能把模型拟合到饱和。

对于同一个 PVR 模型，可按类群汇总其预测表现。若 yhat 来自交叉验证，
一个清晰的类群级指标是

\[
R^2_{pred,c}=1-\frac{\sum_{i\in c}(y_i-\hat y_i)^2}
                         {\sum_{i\in c}(y_i-\bar y_c)^2}.
\]

这里采用“类群内均值基准”。结果可为负：负值表示模型预测误差比该基准
更大。类群性状几乎不变时，分母很小，R² 极不稳定，应该同时报告 RMSE、
样本量与实际变异量；恒定性状的 R² 不定义。这个汇总指标也不是各类群
分别重拟合得到的训练 R²。

可以直接使用扫描结果中的类群名单计算：

```r
# observed 和 predicted_cv 必须是按物种命名的数值向量。
# predicted_cv 是交叉验证预测；每轮重新选择特征向量和拟合所有参数。
diagnostics <- do.call(rbind, lapply(seq_along(h$selection$tips), function(i) {
  species <- h$selection$tips[[i]]
  y <- observed[species]
  yhat <- predicted_cv[species]
  ok <- is.finite(y) & is.finite(yhat)
  y <- y[ok]
  yhat <- yhat[ok]
  sst <- sum((y - mean(y))^2)
  sse <- sum((y - yhat)^2)
  data.frame(node = h$selection$nodes[i], n = length(y),
    R2_pred = if (length(y) >= 3 && sst > 0) 1 - sse / sst else NA_real_,
    RMSE = if (length(y)) sqrt(mean((y - yhat)^2)) else NA_real_)
}))

# 预测后仍残留哪些类群内结构？
residual <- observed - predicted_cv[names(observed)]
residual_scan <- fast_hetero(
  tree, residual, nodes = h$selection$nodes, min_tips = 20,
  method = "K", progress = FALSE
)
```

留出同类群内物种评估“有近缘训练物种时能否预测”；整类群留出评估
“能否推广到未见类群”。这两个问题不同，应据研究目的选择，并保持特征
向量选择等使用性状信息的步骤在训练折内。训练残差中的低 K 是模型已经
移除了结构的诊断，不能单独证明预测能力好。传统 PVR 对全树距离矩阵和
特征分解的需求仍存在；本功能不使 PVR 本身变成适用于任意大型树的算法。

若目标是解释形成异质性的机制，可以比较同一过程与多类群过程的模型，
例如允许 lambda、OU 参数或测量误差在预定类群间变化，并用合适的参数
bootstrap 或预测验证评价差异。仅改变 BM 的扩散速率会整体改变类群的
性状尺度，K 对这种尺度变化不敏感，因此 K 扫描不是演化速率异质性检测器。
上述机制模型没有包含在此次功能中。

## 大型树的计算范围

树索引采用迭代遍历；子树直接从预序区间提取，避免每个小类群都扫描全树。
K 和 lambda 使用现有的编译树算法，D 的树因子、标准化和条件概率也在 C++
计算并复用。Delta 的 ACE、熵／MCMC、模拟及诊断主计算使用现有或新建
C++ 内核，R 保留 `nlminb` 优化和结果编排。共同 Mk 转移缓存限制为 64 MiB；
各次模拟复用同一状态数组，局部 ACE 工作区按保留树及状态数缓存。
`simulation_chunk` 不改变 Delta 的串行随机数顺序，也不自动缩短 MCMC。
默认自动模式选出的类群互不重叠，
子树总大小不超过全树；显式嵌套扫描用 max_total_tips 限制累计工作量。
K 模拟按 simulation_chunk 分块；lambda／D 每次复用一个模拟向量，
临时内存随当前类群大小和分块宽度增长，
另保存每个性状的 nsim 个离散度累计量。不会建立全树稠密协方差矩阵。

总耗时仍随有效类群总大小、性状数、模拟数增长，lambda 优化还取决于
迭代次数。开发验证覆盖了数值等价性、缺失数据、RNG 分块一致性、一个
确定性的强弱对照例子及大型树的实际 smoke workload；它不是完整的功效
或第一类错误率研究，不能作为任意树与任意演化过程的性能保证。

## 方法来源

- Keck et al. (2016), *phylosignal: an R package to measure, test, and explore
  the phylogenetic signal*. [论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4799788/)
  及官方函数说明：[lipaMoran](https://search.r-project.org/CRAN/refmans/phylosignal/html/lipaMoran.html)、
  [phyloSignalINT](https://search.r-project.org/CRAN/refmans/phylosignal/html/phyloSignalINT.html)。
- Blomberg, Garland & Ives (2003), *Testing for phylogenetic signal in
  comparative data: behavioral traits are more labile*.
  [DOI: 10.1111/j.0014-3820.2003.tb00285.x](https://doi.org/10.1111/j.0014-3820.2003.tb00285.x)。
- Fritz & Purvis (2010), *Selectivity in mammalian extinction risk and threat
  types: a new measure of phylogenetic signal strength in binary traits*.
  [原始 D 论文](https://pubmed.ncbi.nlm.nih.gov/20184650/)。该论文的模拟评价
  支持大于 50 物种时的稳健性；它不定义本包的 ER 异质性扩展。
- Diniz-Filho, de Sant'Ana & Bini (1998), *An eigenvector method for
  estimating phylogenetic inertia*. [PubMed](https://pubmed.ncbi.nlm.nih.gov/28565378/)。
- Diniz-Filho et al. (2012), *Exploring patterns of interspecific variation
  in quantitative traits using sequential phylogenetic eigenvector regressions*.
  [DOI: 10.1111/j.1558-5646.2011.01499.x](https://doi.org/10.1111/j.1558-5646.2011.01499.x)。

以上文献支持局部信号、K 与 PVR 的方法基础。此处的自动分区及 BM 下
SD(log K) 参考校准是本包实现的组合工作流，不能将其归称为上述论文
已经发表的标准异质性检验。共享 lambda 高斯校准和保留状态数量的 D 条件
ER 校准同样是本包新增的模型参考工作流。
