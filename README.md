# Quantum Algorithm for Glitchy Period Finding Problem (Summer 2019 Research with Shelby Kimmel)

* Analyzed the change in the probability distribution of finding the period in a modified Quantum Period Finding problem where there are random insertions or deletions (glitches) of points in the periodic function.
* Developed a robust and efficient algorithm to deal with the glitches. Proved performance guarantee for the simple case of a few glitches. Worked on finding ways to extend the analysis for arbitrary number of glitches. 
* Wrote a numerical simulation using python to analyze the behaviors corresponding to the change in the places of insertions/deletions, and to test out the predictions of the our algorithm and analysis. 

## Contents

| Path | What it is |
|---|---|
| `glitchy_period_finding_report.pdf` | Status report: model, one-glitch theorem, numerics, open problems (source: `glitchy_period_finding_report.typ`) |
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

To rebuild the write-up: `typst compile glitchy_period_finding_report.typ`.
