# fastphylosig technical roadmap

Routine and advanced use follow the same sequence: **input -> tree and data
preparation -> phylogenetic signal calculation -> result**.

- **Routine use** calls only `fast_signal()`; preparation and calculation are
  automatic.
- **Advanced use** prepares the tree and data explicitly, then calls the
  method-specific function.

```mermaid
flowchart LR
    I["Common input<br/>phylo tree + named trait data"]

    subgraph GENERAL["Routine use: one command"]
        G1["fast_signal(tree, data, method)"]
        subgraph AUTO["Automatic preparation"]
            G21["Check tree"] --> G22["Match species"] --> G23["Handle NA"] --> G24["Prepare context"]
        end
        G3["Automatic calculation<br/>method: K / lambda / D / Delta"]
        G1 --> G21
        G24 --> G3
    end

    subgraph ADVANCED["Advanced use: explicit control of the same workflow"]
        subgraph PREP["Tree and data preparation"]
            A21["check_tree()"] --> A22["resolve_tree()<br/>when needed"] --> A23["match_tree_data()"] --> A24["prepare_tree()"]
        end
        A3["Method-specific function<br/>K: fast_k()<br/>lambda: fast_lambda()<br/>D: fast_d()<br/>Delta: fast_delta()"]
        A24 --> A3
    end

    O["Unified results<br/>native result type and fields<br/>workflow / timing<br/>P / MCSE / ESS / R-hat when applicable"]

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

The static page uses [`technical-roadmap.svg`](technical-roadmap.svg). The
editable source is [`technical-roadmap.drawio`](technical-roadmap.drawio).
