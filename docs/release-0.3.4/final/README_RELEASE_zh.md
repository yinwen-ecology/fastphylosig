# fastphylosig 0.3.4 发布交接

状态：最终本地源码已冻结，发布材料已准备；尚未上传、打 tag 或正式发布。

安装源码包为 [fastphylosig_0.3.4.tar.gz](fastphylosig_0.3.4.tar.gz)，SHA256：

```text
8B1A8940101997492D4F3E2AFE696717BFE63CC9905FE91E225A3534754B997D
```

该归档直接复制自验收后的最终包，没有重建或修改安装包字节。唯一后续材料是说明、证据、清单和可验证的源文件恢复快照。版本、默认参数、计算结果、随机流及实验性标注保持冻结范围。

- GitHub 标题：`fastphylosig 0.3.4`，预定 tag：`v0.3.4`。
- 可直接使用的发布正文：[中文](RELEASE_NOTES_zh.md)、[English](RELEASE_NOTES_en.md)。
- 解释与推荐参数：[统计范围](STATISTICAL_SCOPE_zh.md)、[完整指南](inst/HETERO_GUIDE_zh.md)。
- 证据：[验证摘要](VALIDATION_zh.md)、[性能口径](PERFORMANCE_zh.md)、evidence/ 原始记录。
- 发布前后步骤：[发布操作清单](PUBLISH_zh.md)。`cran-comments.md` 是未来提交材料草稿，不触发提交，也不修改早先 0.3.1 的提交记录。

建议公开附件为安装源码包、`fastphylosig-0.3.4-release-materials.zip` 和 `SHA256SUMS.txt`。材料 ZIP 包含发布说明、范围、指南、验证及性能证据。`provenance/` 的源文件恢复快照和仓库状态用于内部发布交接，不必作为 GitHub 附件。

运行完整冻结校验：

```text
python verify_release.py
```

它验证清单中的全部文件、两个公开附件的哈希、源码包安全路径和版本、选定源文件快照及实验性范围。还可运行 `python verify_release.py --workspace <仓库根目录>` 核对工作区，或 `python verify_release.py --repo <仓库目录> --git-ref <待发布提交>` 核对候选提交。Git 检查仅对 UTF-8 文本规范化 CRLF/LF，二进制逐字节比较。

当前 Git HEAD 属于旧 0.3.1 发布工作，0.3.4 有大量尚未提交的源码与测试。因此不能直接对当前 HEAD 打 v0.3.4 tag，或把 GitHub 自动生成的旧提交 source archive 当成本次安装包。请按 provenance/repository-paths.txt 审阅并纳入目标提交，核验通过后再创建 tag；不要执行不加审阅的 `git add -A`。本轮未提交或移动现有分支。

冻结是内容哈希和清单约束，不是声称操作系统已把文件设为不可修改。后续任何字节改动必须生成新的清单并重新判断受影响的验收；正式发布使用这里记录的确切归档。
