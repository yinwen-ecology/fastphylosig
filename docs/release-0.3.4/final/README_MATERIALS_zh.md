# fastphylosig 0.3.4 配套材料

本 ZIP 保存与冻结 0.3.4 源码包对应的发布说明、统计范围、指南和验收证据。它不是新一轮统计校准，也不表示已完成 CRAN 提交或审核。

安装包为 fastphylosig_0.3.4.tar.gz，SHA256：

```text
8B1A8940101997492D4F3E2AFE696717BFE63CC9905FE91E225A3534754B997D
```

解压后在 R 中运行 `install.packages("fastphylosig_0.3.4.tar.gz", repos=NULL, type="source")`。需要已安装的 R 依赖及 C++17 工具链。

- RELEASE_NOTES_zh.md / RELEASE_NOTES_en.md：可用作发布正文。
- STATISTICAL_SCOPE_zh.md：K、λ、D、Δ 的解释、切割尺度和推荐参数。
- VALIDATION_zh.md：实际 Windows 检查、维护者确认的 Linux 结果与证据限制。
- PERFORMANCE_zh.md：保留原版本标识的受控加速实验和历史规模基准。
- inst/：完整异质性指南、实验性 nuisance 记录和现有可执行示例。
- manifests/：安装包逐文件哈希、14 个公开函数接口和原生内核对应关系。
- evidence/：检查日志、结果/随机流回归以及性能原始配对表。历史数据保留原版本标签，不改称 0.3.4 新计时。

统计范围保持：K 为 BM 参考；λ 为共同局部冠群参数的近似参考；D 异质性扩展是实验性；Δ 推荐局部与总体描述，全部异质性 P 不在推荐流程中。保留全部失败、拟合/MCMC 警告和成功数值恢复。不能把软件检查通过解释为对所有进化模型的误报率保证。
