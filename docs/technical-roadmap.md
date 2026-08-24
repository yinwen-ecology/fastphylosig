# fastphylosig 技术路线图

一般路径和高级路径按相同顺序排列：**输入 → 树与数据准备 → 系统发育信号计算 → 结果**。

- **一般使用**：只运行 `fast_signal()`，准备和计算均由包自动完成。
- **高级使用**：先显式检查和准备树与数据，再调用对应的方法函数。

```mermaid
flowchart LR
    I["共同输入<br/>phylo 树 + 带物种名的性状数据"]

    subgraph GENERAL["一般使用：简单的一条命令"]
        G1["fast_signal(tree, data, method)"]
        subgraph AUTO["自动准备"]
            G21["检查树"] --> G22["匹配物种"] --> G23["处理 NA"] --> G24["准备 context"]
        end
        G3["自动计算<br/>按 method 选择 K / lambda / D / Delta"]
        G1 --> G21
        G24 --> G3
    end

    subgraph ADVANCED["高级使用：显式控制同一流程"]
        subgraph PREP["树与数据整理"]
            A21["check_tree()"] --> A22["resolve_tree()<br/>必要时"] --> A23["match_tree_data()"] --> A24["prepare_tree()"]
        end
        A3["选择方法函数<br/>K → fast_k()<br/>lambda → fast_lambda()<br/>D → fast_d()<br/>Delta → fast_delta()"]
        A24 --> A3
    end

    O["统一结果<br/>方法原生结果类与字段<br/>workflow / timing · P / MCSE · ESS / R-hat（按方法）"]

    I --> G1
    I --> A21
    G3 --> O
    A3 --> O

    classDef api fill:#dae8fc,stroke:#6c8ebf,color:#17365d
    classDef prep fill:#e1d5e7,stroke:#9673a6,color:#3f2a47
    classDef method fill:#d5e8d4,stroke:#82b366,color:#244321
    classDef output fill:#ffe6cc,stroke:#d79b00,color:#5a3a00
    class I,G1 api
    class G21,G22,G23,G24,A21,A22,A23,A24 prep
    class G3,A3 method
    class O output
```

网页若不启用 Mermaid，可直接使用同目录的
[`technical-roadmap.svg`](technical-roadmap.svg)；需要继续编辑时打开
[`technical-roadmap.drawio`](technical-roadmap.drawio)。
