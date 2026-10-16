# fastphylosig 0.3.4 GitHub 发布记录

已于 2026-10-16 14:44:22 UTC（北京时间 22:44:22）正式发布，非草稿、非预发布。

- Release：[fastphylosig 0.3.4](https://github.com/yinwen-ecology/fastphylosig/releases/tag/v0.3.4)
- 网站：[首页](https://yinwen-ecology.github.io/fastphylosig/)
- 新案例：[比较不同支系的系统发育信号](https://yinwen-ecology.github.io/fastphylosig/fast-hetero-case.html)
- 冻结源码提交/tag：[`3cc0095`](https://github.com/yinwen-ecology/fastphylosig/commit/3cc0095a97f81e89272b182af5d73f13e2881b25)，`v0.3.4` 指向该提交。
- 网站及 GitHub 介绍提交：[`c8d9a4f`](https://github.com/yinwen-ecology/fastphylosig/commit/c8d9a4f7ca2d852af89762ce9e6df46567def77b)，已推送到 main 与 release/v0.3.4。

发布附件均已从公开地址重新下载核对：

| 附件 | SHA256 |
|---|---|
| fastphylosig_0.3.4.tar.gz | `8B1A8940101997492D4F3E2AFE696717BFE63CC9905FE91E225A3534754B997D` |
| fastphylosig-0.3.4-release-materials.zip | `46A4DA6056CF2BB1247794145D4E32658BEA100F7438427D4F04D8AE5978B1C3` |
| SHA256SUMS.txt | `BFB2DAD662C991679ECF0BCC5704DE83F86541C25DBFFA3856154F05B2095E9C` |

v0.3.4 的全部 200 个选定源文件已与冻结快照核对；49 个冻结材料文件在 Git 中保留原始字节，目录属性防止跨平台换行转换破坏哈希。安装包没有重建。随后的网站提交只改变 README 与 docs 文件，R、C++ 和测试代码不变。

公开下载的源码包已在新临时库中完成原生编译和安装。K、λ、D、Δ 均通过加载和基础计算检查，局部与整体估计有限，未触发告警。记录见 publication/PUBLIC_INSTALL_SMOKE.csv、PUBLIC_INSTALL_STATUS.json 及 public-install 日志。

GitHub Pages 已构建并部署网站提交。浏览器检查确认桌面和 390 像素窄屏布局；线上另核验首页介绍、案例页、三个结果表、全部章节锚点及 11 个案例文件。案例采用固定模拟的 512 物种、8 个预定冠群，给出原生树图、完整结果、代码与下载数据；行文以操作和结果解释为主，必要适用范围集中说明。

统计定位保持：K 使用 BM 参考，λ 使用共同局部冠群参数的近似参考；D 的新增异质性扩展保持实验性；Δ 推荐局部与总体描述，异质性检验不在推荐流程。早先 CRAN 提交记录未改动。

final/ 内的状态是发布前的不可改写冻结快照；本文件和 publication/ 中的发布回执记录实际发布状态。自动检查进度可在 [GitHub Actions](https://github.com/yinwen-ecology/fastphylosig/actions) 查看。
