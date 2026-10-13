# Quantum Algorithm for Glitchy Period Finding Problem (Summer 2019 Research with Shelby Kimmel)

* Analyzed the change in the probability distribution of finding the period in a modified Quantum Period Finding problem where there are random insertions or deletions (glitches) of points in the periodic function.
* Developed a robust and efficient algorithm to deal with the glitches. Proved performance guarantee for the simple case of a few glitches. Worked on finding ways to extend the analysis for arbitrary number of glitches. 
* Wrote a numerical simulation using python to analyze the behaviors corresponding to the change in the places of insertions/deletions, and to test out the predictions of the our algorithm and analysis. 

## Contents

| Path | What it is |
|---|---|
| `glitchy_period_finding_report.pdf` | Report, including the general-k proof in §4.1 (source: `glitchy_period_finding_report.typ`) |
| `formal/` | Lean-checked §4.1 proofs, pinned dependencies, and axiom audit; see [formal/README.md](formal/README.md) |
| `report_figures/` | Editable Mermaid proof map and data-derived counting illustration |
| `poster_figures/` | Descriptively named figures cropped from the 2019 poster |
| `glitchy.py` | The simulator (single file) |
| `Summer_Symposium_Poster.pdf` | 2019 poster |
| `Experiment_plots/` | 2019 logs of the worst-case searches |

## Running the simulator

With [uv](https://docs.astral.sh/uv/) the dependencies (numpy, matplotlib) install automatically:

```sh
uv run glitchy.py check                                  # reproduce every number in glitchy_period_finding_report.pdf
uv run glitchy.py plot  --N 11000 --r 100 --k 10         # plot a random placement of 10 glitches
uv run glitchy.py plot  --N 3000 --r 50 --blocks 6 54    # plot a given placement
uv run glitchy.py worst --N 1500 --r 50 --k 4            # worst case over all placements
uv run glitchy.py worst --N 3000 --r 50 --k 3 --signs 1 -1 1   # mixed insertions/deletions
```

The original 2019 scripts (`Plot.py`, `Find_Worst_Case.py`, `Test.py`) are in git history at commit `a9e171b`.

## General-k guarantee (Section 4.1)

For any `k` literal insertions/deletions, when `r` divides `N`, the total
probability at the exact ideal Fourier outcomes is at least

- `r / ((k+1)r + I) >= 1/(2k+1)`, where `I <= k` inserted entries survive;
- independently, `(1-I/N)^2 / (k+1)`.

The report includes a four-step proof, a [Mermaid proof map](report_figures/general_k_proof.md),
and a worked three-glitch count table. These are **total-mass** guarantees, not
the stronger per-outcome conjecture or a period-recovery theorem.

The [Lean project](formal/README.md) checks the literal edit bookkeeping,
Fourier identity, both bounds and their corollaries, Model I bound and sharpness,
and the worked example. It does not claim to formalize the rest of the report.

```sh
cd formal
lake exe cache get       # first use; requires Lean/elan
lake build
lake env lean Audit.lean
```

## Rebuild the write-up

From the repository root:

```sh
typst compile glitchy_period_finding_report.typ
```

The explanatory SVGs are committed. To regenerate them first:

```sh
mmdc -i report_figures/general_k_proof.md -o report_figures/general_k_proof.svg
python3 report_figures/build_general_k_example.py
```
