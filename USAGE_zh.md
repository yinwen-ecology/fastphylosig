# fastphylosig 中文使用指南

`fastphylosig` 用于在系统发育树上快速计算连续、二元和分类性状的系统发育信号。连续性状使用 `fast_signal()`，二元性状使用 `fast_d()`，分类性状使用 `fast_delta()`。

## 安装

直接从 GitHub 安装：

```r
install.packages("remotes")
remotes::install_github("yinwen-ecology/fastphylosig")
library(fastphylosig)
```

也可以在本地 R 包根目录运行：

```r
install.packages(".", repos = NULL, type = "source")
library(fastphylosig)
```

在线帮助页面：<https://yinwen-ecology.github.io/fastphylosig/>

该包需要可用的 C++ 编译工具。OpenMP 是可选的；没有 OpenMP 时结果不变，但 K 和 D 的 C++ 多线程不会生效。

## 1. 连续性状：K 和 lambda

单个性状必须带物种名：

```r
set.seed(1)
tree <- ape::rtree(50)
x <- ape::rTraitCont(tree)

fit_k <- fast_signal(
  tree, x, method = "K", test = TRUE, nsim = 1000
)
fit_lambda <- fast_signal(
  tree, x, method = "lambda", test = TRUE
)

fit_k$K
fit_k$P
fit_lambda$lambda
fit_lambda$logL
fit_lambda$P
```

一次计算多个性状时，应把性状放在矩阵或 data.frame 的列中，行名为物种名：

```r
X <- matrix(
  rnorm(50 * 100), nrow = 50,
  dimnames = list(tree$tip.label, paste0("trait_", 1:100))
)

K_table <- fast_signal(tree, X, method = "K", test = FALSE)
lambda_table <- fast_signal(tree, X, method = "lambda", test = TRUE)
```

矩阵输入比逐列循环更快，因为具有相同 NA 模式的性状可以共用树匹配、VCV 矩阵和矩阵分解。当前 K 和 lambda 只支持 `se = NULL`。

单个 lambda 结果默认保存 likelihood profile，批量结果默认不保存，以避免额外计算。如确实需要批量 profile：

```r
lambda_profile_table <- fast_signal(
  tree, X, method = "lambda", test = TRUE,
  lambda_profile = TRUE, lambda_profile_points = 101
)
```

## 2. 二元性状：Fritz & Purvis D

```r
x_binary <- as.integer(x > median(x))
names(x_binary) <- tree$tip.label

fit_d <- fast_d(
  tree, x_binary, test = TRUE, nsim = 1000,
  return_sim = TRUE, ncores = 4
)

fit_d$DEstimate
fit_d$Pval1  # random association null
fit_d$Pval0  # Brownian threshold null
```

二元标签会在内部统一转换为 0/1，因此 `0/1`、`1/2`、`10/20` 或两个字符标签不会改变 D。D 本身依赖 random 和 Brownian 两组模拟均值，所以即使 `test = FALSE`，仍需要通过 `nsim` 进行校准；此时只是不返回 P 值。

## 3. 分类性状：Delta

```r
x_cat <- sample(letters[1:3], ape::Ntip(tree), replace = TRUE)
names(x_cat) <- tree$tip.label

fit_delta <- fast_delta(
  tree, x_cat,
  test = TRUE, nsim = 100,
  mcmc_sim = 5000, thin = 10, burn = 100,
  model = "ARD", ncores = 4
)

fit_delta$delta
fit_delta$P
```

`fast_delta()` 默认使用 `fast_ace()` 计算离散祖先状态概率。与 `ape::ace()` 对照时可设置 `ace_engine = "ape"`。支持 `model = "ER"` 和 `model = "ARD"`。

Delta 的置换检验比 K 和 D 更耗时。建议先用较小的 `nsim` 和 `mcmc_sim` 做测试，再扩大。批量结果中的 `n_failed_sim` 会报告失败的置换拟合；如果全部失败，P 值为 `NA`。

## 4. 检查树和数据匹配

```r
matched <- match_phylo_data(tree, X)
matched$report
matched$tree_tips_removed
matched$data_rows_removed
```

包不会修改原始数据。结果会报告匹配物种数、删除的树 tips、删除的数据行和每个性状因 NA 删除的物种数。NA 按性状分别处理。

高级参数 `permutations` 必须是 `nsim` 行的整数矩阵，每行都必须完整包含一次 `1:n`。如果不同性状的 NA 模式不同，它们保留的物种数也不同，不能共用同一个 permutation 矩阵。

## 5. 绘图

```r
plot_signal(fit_k)
plot_signal(fit_lambda)
plot_signal(fit_d)
plot_signal(fit_delta)
```

- K：随机置换 K 分布、observed K 和右尾 P 区域。
- lambda：profile likelihood、lambda 估计值、lambda = 0/1、LR 检验和近似 95% CI。
- D：在同一坐标轴叠加 Brownian threshold null 与 random association null，并显示 observed D、D = 0/1、`P_random` 和 `P_Brownian`。
- Delta：置换 Delta 分布、observed Delta 和右尾 P 区域。

K、D 和 Delta 的分布图要求计算时使用 `test = TRUE, return_sim = TRUE`。D 的极端模拟次数为 0 时，图中显示 `P < 1/nsim`，不会显示 `P = 0`。

## 6. 如何计算更快

- 多性状优先一次传入矩阵，不要逐列调用。
- K 只需要 P 值、不需要保存模拟分布时，设置 `return_sim = FALSE`。
- K 随机化和 D 模拟可使用 `ncores`；在 8 核电脑上最多使用 `ncores = 8`。
- Delta 的 `ncores` 用于并行置换；单个 observed Delta 不会因此明显加速。
- 单个 observed K 或 lambda 主要受 VCV 和稠密矩阵分解限制，多核收益通常较小。

## 7. 当前限制

- K 和 lambda 尚不支持测量误差 `se != NULL`。
- `fast_ace()` 只支持离散性状、ML、ER/ARD 和有枝长的有根完全二叉树。
- 返回对象保留主要统计量和绘图所需数据，但不会复制 `phytools`、`caper` 或 `ape` 对象中的所有内部字段。
- 当前验证基于 `phytools` 1.0.3、`caper` 1.0.1 以及本机安装的 `ape`；参考包升级后应重新运行测试。

## 8. 检查包

```r
devtools::test()
devtools::check(args = "--no-manual")
```

核心测试会比较 K/lambda 与 `phytools`、D 与 `caper`、ER/ARD 的 ACE 结果与 `ape::ace()`。
