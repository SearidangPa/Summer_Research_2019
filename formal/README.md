# Lean proofs for the general-k bounds (Sections 4.1–4.2)

**The general-k literal-edit theorem is checked end to end**, starting with
actual insertions/deletions of an infinite tagged sequence and ending with the
finite Fourier measurement formula. Section 4.1's Fourier identity, displacement and support bounds, and Cauchy–Schwarz
step are proved, not supplied as assumptions. Section 4.2's nondivisible theorem
is also checked end to end: matching moments, literal pair counts, rounded-window
geometry, the finite Fejér kernel, and the constant `133/358`.

Lean 4 and mathlib are pinned to **v4.19.0**. There are no `sorry` placeholders,
custom axioms, or `native_decide` proofs.

## Run the checks

With [elan](https://github.com/leanprover/elan) installed:

```sh
cd formal
lake exe cache get       # fetch mathlib's precompiled proofs on first use
lake build
lake env lean Audit.lean
```

The checked-in `lean-toolchain`, `lakefile.toml`, and `lake-manifest.json` pin the
dependencies. The first command may download substantial toolchain and cache data.
The generated `.lake/` directory is ignored by git.

On macOS, if the cache executable reports a `__DATA_CONST` / `SG_READ_ONLY`
loader error, run the same cache program through Lean's interpreter:

```sh
lake env lean --run .lake/packages/mathlib/Cache/Main.lean get
lake build
lake env lean Audit.lean
```

This workaround changes only cache retrieval, not the proofs or checking kernel.

## Notation: report and Lean

The Lean names use the report's letters, so the two can be read side by side. Docstrings
cite the report's equations by their source labels, for example "Eq. (matching-bound)".

| Report | Lean | Meaning |
|---|---|---|
| `N`, `r`, `k` | `N`, `r`, `A + D` | prefix length, period, number of edits |
| `A`, `D` | `A`, `D` in `Edited A D τ` | insertions, deletions |
| `τ(x)` | `τ x` | tag at position `x`: original index `t ≥ 0`, or `-j` for the `j`-th insertion |
| `w_0, ..., w_(r-1)` | `w : ZMod r → W` | periodic symbols (need not be distinct) |
| `v` | `v : ℤ → W` | inserted symbols, indexed by their negative tags |
| `f_k(x)` | `symbol w v (τ x)`; `literalCycle τ w v` on `ZMod N` | edited prefix |
| `I` | `(inserteds τ N).card` | surviving inserted entries, `I ≤ A` |
| `s = x - t` | `x - τ x` in `Edited.displacement` | displacement of an original entry |
| `u`, `b`, box `(u, b)` | `u : W`, `b : ZMod r`, `box τ w v x` | output symbol, residue, box |
| `n_(u b)`, `m_(u b)` | `prefixCounts τ w v N u b`, `m` | all / retained-original box counts |
| `q`, `Q` (Step 2) | card of the support, `Q` | occupied boxes and their bound |
| `P_f(y)` | `fourierProbability f y`, `cyclicProbability f y` | Eq. (pf), indexed by `ℕ` / `ZMod N` |
| `G`, `G_k` (§4.1) | `goodOutcomes N r`, any `H ⊇ G` | ideal outcomes, a neighbourhood |
| `ℓ_j`, `s_j`, `p_b` (Model I) | `length j`, `shift j`, `grouped shift length b` | blocks `0..k` in Lean, `1..k+1` in the report |
| `h`, `C_h` | `h`, `matchingFraction f h` | lag, circular matching fraction |
| good start `x` | `x ∈ goodPairs τ h L` | `τ x ≥ 0` and `τ (x+h) = τ x + h` |
| `g_c`, `G_k` (§4.2) | `g c` with `NearestPeaks g`, `roundedWindows g k` | rounded peaks, their radius-`k` windows |
| `d_N(y, g) = min_m abs(y - g - m N)` | `∃ m : ℤ, abs(y - g - m*N) ≤ k` | circular distance |
| `t`, `p`, `θ_y` | `t`, `p`, `phaseDistance (r*y)` | scale `2N/(r(2k+1))`, window mass, phase distance |
| `T`, `F_T`, `E_T` | `T`, `fejerSquare T / T = F_T`, `E` | kernel size, Fejér kernel, expectation `E_T/T` |
| `q` (Step 3) | `q` | bound on `F_T/T` outside `G_k` |
| `i` (imaginary unit) | `Complex.I` | written in full; `I` always means surviving insertions |

## What is proved

All theorem names below are in the `Glitchy` namespace.

| Part of the argument | Source | Main declarations |
|---|---|---|
| Literal insertion/deletion histories; displacement in `[-D,A]`; distinct tags; at most `A` surviving insertions | `Glitchy/Edits.lean` | `Edited`, `Edited.displacement`, `Edited.inserted_card_le` |
| At most `(A+D+1)r+I` occupied value–residue boxes; actual fiber counts | `Glitchy/Boxes.lean` | `Edited.original_boxes`, `Edited.all_boxes`, `sum_fiberCounts` |
| Cauchy–Schwarz; both finite-size bounds; the universal constant; few-edit corollary; asymptotic factor | `Glitchy/Counting.lean` | `mass_from_boxes`, `mass_from_retained`, `universal_bound`, `few_edits`, `retained_factor_tendsto` |
| Character orthogonality and exact squared-count identity | `Glitchy/Fourier.lean` | `character_orthogonality`, `residue_parseval`, `residue_mass` |
| Bridge from the actual edited prefix to its complex-exponential Fourier formula; distinct ideal outcomes; assembled bounds | `Glitchy/Literal.lean` | `idealProbability_eq_fourier`, `general_k_glitches`, `general_k_sharp_bounds`, `general_k_neighbourhood` |
| No-surviving-insertions and few-edits corollaries | `Glitchy/Literal.lean` | `general_k_no_surviving_insertions`, `general_k_few_edits` |
| Model I's grouped-length identity, inverse-linear bound, equality case, and absence of a positive uniform exact-set bound | `Glitchy/ModelI.lean` | `modelMass_eq`, `modelI_bound`, `modelI_sharp`, `no_uniform_exact_bound` |
| The illustration's actual three-edit history, count table, surviving insertions, and exact Fourier mass `1/3` | `Glitchy/Examples.lean` | `example_edited`, `example_counts_correct`, `example_fourier_mass` |
| Circular Fourier matching moments, normalization, and the complex-exponential bridge | `Glitchy/Matching.lean` | `cyclic_matching_moment`, `cyclicProbability_sum`, `cyclicProbability_eq_fourier` |
| Matching-pair counts by induction on actual edit histories | `Glitchy/PairCounting.lean` | `Edited.goodPairs_bound`, `Edited.matching_bound` |
| Finite expectation inequality, normalization, and both trigonometric regimes | `Glitchy/Fejer.lean` | `mass_from_test_function`, `short_sine_bound`, `short_fejer_mass`, `long_fejer_lower`, `long_fejer_outside` |
| Squared geometric sum, triangular expansion, matching expectation, and concentration | `Glitchy/FejerKernel.lean` | `fejerSquare_expansion`, `fejerSquare_geometric`, `fejer_expectation_lower`, `fejer_mass_of_matching` |
| Arbitrary nearest-peak choices, circular windows, coverage, and separation | `Glitchy/PeakGeometry.lean` | `canonicalPeaks_nearest`, `roundedWindows_full_of_lt`, `roundedWindows_separation` |
| Complete constant-mass theorem and original inverse-linear guarantee, without divisibility | `Glitchy/Nondividing.lean` | `general_k_nondividing`, `general_k_nondividing_universal`, `general_k_nondividing_canonical` |
| Rounding ties, mixed edits, deleted insertions, and analytic boundary regressions | `Glitchy/NondividingExamples.lean` | `nondividing_tie_regression`, `nondividing_mixed_regression`, `nondividing_long_regression` |

### Reading the main theorem

`general_k_sharp_bounds` uses `A` insertions and `D` deletions rather than a
separate `k`; set `k = A + D`. For any history `h : Edited A D τ`, any symbol
assignments, positive `N`, and positive `r` dividing `N`, it proves both

- total Fourier mass on the **set** of exact ideal outcomes at least
  `r / ((A+D+1)r + I)`;
- the same mass at least `(N-I)^2 / ((A+D+1)N^2)`.

Here `I` is **computed** from surviving inserted tags in the first `N` entries.
The support and displacement conditions are consequences of `h`, not extra
hypotheses. `good_index_injective` prevents double-counting ideal outcomes.
`general_k_neighbourhood` supplies the universal `1/(2k+1)` bound for any finite
superset of the exact ideal outcomes; overlapping neighbourhoods are harmless.

### Reading the nondivisible theorem

`general_k_nondividing` takes an actual history `hedit : Edited A D τ`, arbitrary
periodic and inserted symbol assignments, `N > 0`, `r >= 2`, and `A+D >= 1`.
For **any** choice `g` of nearest integers to `cN/r`, it proves

```text
133/358 <= sum_{y in roundedWindows g (A+D)} fourierProbability f y.val
```

`NearestPeaks g` means exactly `abs(g c - c.val*N/r) <= 1/2`. It is a geometric
specification, not an assumed concentration bound. `canonicalPeaks_nearest`
constructs a valid choice using `round`, and `general_k_nondividing_canonical`
needs no rounding hypothesis. The theorem covers either choice at every tie.
`roundedWindows` uses integer translates by `N` to express circular distance,
and is a finite **set** of outcomes, so overlaps are not double-counted.

There is no divisibility, sample-size, small-edit, or distinct-symbol hypothesis.
As in the report, `r > N` is allowed; the windows then cover everything.
`general_k_nondividing_universal` recovers the original `1/(2k+1)` bound.

The proof follows the report step by step. `Edited.goodPairs_bound` is the
good-pairs induction of Step 1, giving `C_h >= 1-((k+1)h+A)/N`.
`roundedWindows_full_of_lt` is the coverage case `N < (2k+1)r` of Step 2; otherwise
`t >= 2`. `fejer_mass_of_matching` uses the short kernel `T = 2` for `2 <= t <= 6`
and the long kernel `T = ceil(3t/4)` for `t > 6`.

No matching estimate, Fourier identity, Fejér expansion, kernel envelope, or
expectation bound is an extra hypothesis of the final theorem.

### How literal edits are represented

An unedited sequence has tags `τ x = x`. Original tags are nonnegative.
Insertion number `A` gets a fresh tag `-(A+1)`; deleting it simply removes that
tag. `insertAt` and `deleteAt` shift the rest of an infinite sequence exactly
as their names suggest.

The final symbol attached to an original tag is its periodic value.
Inserted tags receive arbitrary values from a separate function. Thus two
insertions may have the same output, or may match an original output, while
their provenance remains distinguishable. This avoids silently assuming that
inserted symbols are new or distinct.

The formal alphabet `W` is finite. This is no substantive restriction for a
fixed finite edit history: the original period and the finitely many inserted
symbols use a finite alphabet, even when the ambient output space is infinite.
The formal lower bound does not even need the original period symbols to be
distinct.

### Exact scope and trust boundary

- The Fourier measurement formula is the mathematical definition of
  `fourierProbability`. Its equivalence to the count-table formula is proved;
  implementation of a quantum circuit is not formalized.
- `modelMass` is the exact-outcome formula of Model I. Normalized block lengths,
  their grouping, and the sharpness calculation are checked. This does not
  verify a physical oracle interpretation of overlapping Model I position sets.
- The asymptotic assertion uses the checked finite-size inequality and
  `retained_factor_tendsto`: when `k/N` tends to zero, `(1-k/N)^2` tends to one.
- This formalization covers the **new general-k arguments in Sections 4.1–4.2**.
  It does not claim to formalize the report's other continuum theorems,
  transfer sketch, simulations, or period-recovery proposal. The radius-k
  per-outcome conjecture remains open.

`Audit.lean` reports transitive axiom dependencies of the main results.
Only Lean's standard `propext`, `Classical.choice`, and `Quot.sound` are used;
there is no `sorryAx` or added mathematical axiom. The small example is evaluated
by kernel-checked `decide`, not an unverified numerical computation.
