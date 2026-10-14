# General-k total-mass proof map

A **box** is a pair (u, b): output symbol u and input residue b modulo r. The diagram proves a
bound on the **sum over all ideal outcomes**, not on each outcome separately.

```mermaid
%%{init: {"theme":"base","htmlLabels":false,"themeVariables":{"fontFamily":"Arial","fontSize":"17px","primaryColor":"#eef4fa","primaryTextColor":"#183a5a","primaryBorderColor":"#51799d","lineColor":"#51799d"},"flowchart":{"htmlLabels":false,"curve":"linear","rankSpacing":22,"padding":12,"wrappingWidth":620}}}%%
flowchart TD
    A["k edits: at most k + 1 original displacements"]
    B["Original entries: at most (k + 1)r boxes"]
    C["Add I surviving insertions: at most Q = (k + 1)r + I boxes"]
    D["Cauchy-Schwarz: sum of squared counts >= N^2 / Q"]
    E["Fourier averaging: P(G) = (r / N^2) times squared counts"]
    F["P(G_k) >= P(G) >= r / Q >= 1 / (2k + 1)"]
    A --> B
    B --> C
    C --> D
    D --> E
    E --> F
    style F fill:#eef7ef,stroke:#387449,color:#204e2e
```

Render from the repository root:

```sh
mmdc -i report_figures/general_k_proof.md -o report_figures/general_k_proof.svg
```

Mermaid CLI extracts the diagram and writes `general_k_proof-1.svg`, which
Section 4.1 of `glitchy_period_finding_report.typ` includes directly.
