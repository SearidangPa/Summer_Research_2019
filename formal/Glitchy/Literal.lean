import Glitchy.Boxes
import Glitchy.Fourier

-- Putting it all together for an actual edited sequence:
--   1. `literal_total`, `literal_retained`, `literal_universal`: lower bounds on the total
--      ideal-outcome probability in its count-table form (Boxes + Counting + Fourier).
--   2. `idealProbability_eq_fourier`: the count-table formula equals the true QFT
--      measurement probability `fourierProbability`.
--   3. `general_k_*`: the final theorems, stated in terms of `fourierProbability`.
--
-- Lean ideas introduced here:
--   * using a theorem by applying it to arguments, like a function
--   * term-mode proofs that chain inequalities with `.trans`
--   * `let` inside a theorem statement, unfolded with `dsimp`

namespace Glitchy
open Finset Complex
-- W must now be finite (`Fintype W`), because we sum over every symbol.
variable {W : Type*} [Fintype W] [DecidableEq W] {r : ℕ} [NeZero r]

-- `(box s w v : Fin N → W × ZMod r)` is a type ascription. It tells Lean which N to use.
/-- Actual value-residue counts of the first N entries after the edits:
how many of the first `N` positions show symbol `a` and are `≡ b (mod r)`. -/
noncomputable def prefixCounts (s : ℕ → ℤ) (w : ZMod r → W) (v : ℤ → W)
    (N : ℕ) (a : W) (b : ZMod r) : ℝ :=
  fiberCounts (box s w v : Fin N → W × ZMod r) univ (a, b)

-- `((A + D : ℕ) + 1 : ℝ)`: compute A + D in ℕ, cast it into ℝ, and add 1.
/-- The box bound connected to Fourier probabilities, not assumed as a hypothesis:
for `h : Edited A D s` and `N > 0`, `r / ((A+D+1)·r + I) ≤ ∑_c P(c)`, where `I` is the
number of inserted entries among the first `N`. -/
theorem literal_total {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) :
    (r : ℝ) / (((A + D : ℕ) + 1 : ℝ) * r + (inserteds s N).card) ≤
      ∑ c : ZMod r, idealProbability N (prefixCounts s w v N) c := by
  -- Replace the sum of probabilities with (r/N²)·∑ (count)² (Fourier.lean).
  rw [residue_mass]
  -- n is the count table, viewed as a function on boxes (pairs (a, b)).
  let n := fiberCounts (box s w v : Fin N → W × ZMod r) univ
  -- The counts add up to N.
  have hsum : ∑ b, n b = (N : ℝ) := by
    simpa [n] using sum_fiberCounts (box s w v : Fin N → W × ZMod r) univ
  -- At most (A+D+1)r + I boxes are occupied.
  have hq : (univ.filter (fun b => n b ≠ 0)).card ≤
      (A + D + 1) * r + (inserteds s N).card := by
    -- The occupied boxes are the image of the box map ...
    rw [support_fiberCounts]
    -- ... and `Edited.all_boxes` bounds the size of that image.
    exact h.all_boxes w v N
  -- Apply the inequality from Counting.lean. Each argument fills in one of its inputs, in
  -- order. The `(by ...)` argument proves Q > 0 (`NeZero.pos r` gives 0 < r).
  have hb := mass_from_boxes n N r ((A + D + 1) * r + (inserteds s N).card)
    hN (by have := NeZero.pos r; positivity) hsum hq
  -- hb sums over pairs (a, b), while the goal has a double sum ∑_a ∑_b.
  -- `Fintype.sum_prod_type` converts between the two, and the `Nat.cast` lemmas fix casts.
  simpa [n, prefixCounts, Fintype.sum_prod_type, Nat.cast_add, Nat.cast_mul] using hb

-- This uses only the original entries: there are N - I of them, and they fill at most
-- (A+D+1)r boxes.
/-- The second bound: `(N - I)² / ((A+D+1)·N²) ≤ ∑_c P(c)`. -/
theorem literal_retained {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) :
    ((N : ℝ) - (inserteds s N).card)^2 /
      (((A + D : ℕ) + 1 : ℝ) * (N : ℝ)^2) ≤
      ∑ c : ZMod r, idealProbability N (prefixCounts s w v N) c := by
  rw [residue_mass]
  -- n counts every entry and m counts only the original entries.
  let n := fiberCounts (box s w v : Fin N → W × ZMod r) univ
  let m := fiberCounts (box s w v : Fin N → W × ZMod r) (originals s N)
  -- There are N - I original entries.
  have hsum : ∑ b, m b = (N : ℝ) - (inserteds s N).card := by
    rw [sum_fiberCounts]
    have hc := original_card s N
    have hc' : ((originals s N).card : ℝ) + (inserteds s N).card = N := by
      exact_mod_cast hc
    -- `linarith`: linear arithmetic over ℝ. Here it just rearranges hc'.
    linarith
  -- The original entries fill at most (A+D+1)r boxes.
  have hq : (univ.filter (fun b => m b ≠ 0)).card ≤ (A + D + 1) * r := by
    rw [support_fiberCounts]
    exact h.original_boxes w v N
  -- `fiberCounts_mono` with originals ⊆ univ supplies 0 ≤ m b ≤ n b for every b.
  have hb := mass_from_retained n m N r (A + D) (inserteds s N).card hN (NeZero.pos r)
    (fun b => fiberCounts_mono _ (subset_univ _) b) hsum hq
  simpa [n, prefixCounts, Fintype.sum_prod_type] using hb

-- This is a term-mode proof (no `by`). It chains two inequalities with `.trans`:
--   1/(2k+1) ≤ r/((k+1)r + I)   by `universal_bound`, which needs I ≤ k = A + D;
--            ≤ ∑_c P(c)          by `literal_total`.
-- I ≤ A + D comes from I ≤ A (`inserted_card_le`) and A ≤ A + D (`Nat.le_add_right`).
/-- The bound that depends only on `k = A + D`: `1/(2(A+D)+1) ≤ ∑_c P(c)`. -/
theorem literal_universal {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) :
    (1 : ℝ) / (2 * (A + D : ℕ) + 1) ≤
      ∑ c : ZMod r, idealProbability N (prefixCounts s w v N) c :=
  (universal_bound r (A + D) (inserteds s N).card (NeZero.pos r)
    ((h.inserted_card_le N).trans (Nat.le_add_right A D))).trans
    (literal_total h w v N hN)

/-- Grouping a finite sum by the residue boxes, before taking absolute squares:
`∑_b count(b) · F(b) = ∑_{x ∈ T} F(f x)`. -/
theorem weighted_fibers {X B : Type*} [Fintype B] [DecidableEq B]
    (f : X → B) (T : Finset X) (F : B → ℂ) :
    (∑ b, (fiberCounts f T b : ℂ) * F b) = ∑ x ∈ T, F (f x) := by
  -- For each box b: count(b) · F(b) = the sum of F(f x) over the x in box b.
  have h (b : B) : (fiberCounts f T b : ℂ) * F b =
      ∑ x ∈ T.filter (fun x => f x = b), F (f x) := by
    unfold fiberCounts
    push_cast
    calc
      -- count · constant = the sum of that constant over the box (`simp` knows this).
      _ = ∑ _x ∈ T.filter (fun x => f x = b), F b := by simp
      -- Inside the box, f x = b, so F b = F (f x).
      _ = _ := by
        apply sum_congr rfl
        intro x hx
        rw [(mem_filter.mp hx).2]
  simp_rw [h]
  -- `sum_fiberwise`: summing over each box and then over all boxes = summing over T.
  exact sum_fiberwise T f (fun x => F (f x))

-- In the statement, `if symbol w v (s x) = a then ... else 0` keeps only the positions
-- that show symbol a.
/-- The actual (unnormalized) Fourier amplitudes at ideal outcomes agree with
the amplitudes computed from the count table:
`∑_{x < N, symbol at x = a} e(c·x) = ∑_b count(a, b) · e(c·b)`. -/
theorem amplitude_by_residue (s : ℕ → ℤ) (w : ZMod r → W) (v : ℤ → W)
    (N : ℕ) (a : W) (c : ZMod r) :
    (∑ x : Fin N, if symbol w v (s x) = a
      then ZMod.stdAddChar (c * (x.val : ZMod r)) else 0) =
    ∑ b : ZMod r, (prefixCounts s w v N a b : ℂ) * ZMod.stdAddChar (c * b) := by
  -- Apply `weighted_fibers` with F(a', b) = (if a' = a then e(c·b) else 0).
  have h := weighted_fibers (box s w v : Fin N → W × ZMod r) univ
    (fun p => if p.1 = a then ZMod.stdAddChar (c * p.2) else 0)
  -- Turn the sum over pairs (a', b) into a double sum ∑_{a'} ∑_b.
  rw [Fintype.sum_prod_type] at h
  -- Simplify the `if`s. Only a' = a survives in the outer sum.
  simp only [mul_ite, mul_zero, ← sum_ite_irrel, sum_ite_eq', mem_univ, if_true] at h
  -- Unfold `box` and `prefixCounts`. h, flipped with `.symm`, is now the goal.
  simpa [box, prefixCounts] using h.symm

-- This is the probability of outcome y when you prepare (1/√N) ∑_{x<N} |x⟩|f(x)⟩ and apply
-- the quantum Fourier transform to the first register:
--   P(y) = (1/N²) · ∑_a |∑_{x : f(x) = a} exp(2πi·x·y/N)|².
-- `I` is the imaginary unit `Complex.I`; `open Complex` lets us write it without `Complex.`.
-- `(x.val : ℂ)` and `(y : ℂ)` cast natural numbers into ℂ.
/-- Fourier measurement for a literal prefix, using the report's positive exponent.
y is a natural number; only its residue modulo N affects the result. -/
noncomputable def fourierProbability {N : ℕ} (f : Fin N → W) (y : ℕ) : ℝ :=
  (1 / (N : ℝ)^2) * ∑ a : W,
    normSq (∑ x : Fin N, if f x = a
      then Complex.exp (2 * Real.pi * I * (x.val : ℂ) * (y : ℂ) / (N : ℂ)) else 0)

-- Why it holds: y/N = c/r when y = c·(N/r).
-- `r ∣ N` means "r divides N". The symbol is typed `\mid`; it is not the keyboard `|`.
-- `c.val` is c as a natural number in {0, ..., r-1}, and `N / r` is division in ℕ.
/-- At the ideal outcome `y = c·(N/r)`, the QFT phase `exp(2πi·x·y/N)` equals
`e(c·x) = exp(2πi·c·x/r)`. -/
theorem ideal_phase (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) (x : Fin N) (c : ZMod r) :
    ZMod.stdAddChar (c * (x.val : ZMod r)) =
      Complex.exp (2 * Real.pi * I * (x.val : ℂ) * (c.val * (N / r) : ℕ) / (N : ℂ)) := by
  -- Write c·x (in ZMod r) as the cast of the natural number c.val·x.val.
  have hc : c * (x.val : ZMod r) = ((c.val * x.val : ℕ) : ZMod r) := by simp
  -- `ZMod.stdAddChar_coe`: e(j mod r) = exp(2πi·j/r) for an integer j. To use it, first view
  -- the cast from ℕ as a cast from ℤ (`← Int.cast_natCast`).
  rw [hc, ← Int.cast_natCast, ZMod.stdAddChar_coe]
  push_cast
  -- Both sides are now exp(...). `congr 1` reduces the goal to equality of the exponents.
  congr 1
  -- The denominators are nonzero, and (N/r)·r = N because r divides N.
  have hn : (N : ℂ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hN)
  have hr : (r : ℂ) ≠ 0 := by exact_mod_cast (NeZero.ne r)
  have hd : ((N / r : ℕ) : ℂ) * r = N := by exact_mod_cast (Nat.div_mul_cancel hrN)
  -- Clear denominators, replace N by (N/r)·r, and finish with algebra.
  field_simp
  rw [← hd]
  ring

-- This bridge turns bounds about count tables into bounds about the actual measurement.
/-- Bridge from the literal sequence's QFT formula to its residue-count formula:
the count-table probability `P(c)` equals the QFT probability at `y = c·(N/r)`. -/
theorem idealProbability_eq_fourier (s : ℕ → ℤ) (w : ZMod r → W) (v : ℤ → W)
    (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) (c : ZMod r) :
    idealProbability N (prefixCounts s w v N) c =
      fourierProbability (fun x : Fin N => symbol w v (s x)) (c.val * (N / r)) := by
  -- Both definitions have the form (1/N²) · ∑_a |amplitude|².
  unfold idealProbability fourierProbability
  -- The 1/N² factors agree, so `congr 1` leaves the sums to compare.
  congr 1
  -- Compare the sums term by term, for each symbol a.
  apply sum_congr rfl
  intro a ha
  -- Replace the count-table amplitude with the position-by-position one.
  rw [← amplitude_by_residue]
  -- Compare what is inside |·|², then term by term for each position x.
  congr 1
  apply sum_congr rfl
  intro x hx
  -- The phases agree by `ideal_phase`.
  rw [ideal_phase N hN hrN x c]

-- MAIN THEOREM, k-only form. Take any sequence obtained by A insertions and D deletions,
-- periodic symbols w (period r), inserted symbols v, and N > 0 with r ∣ N. Then
--   1/(2(A+D)+1) ≤ ∑_{c ∈ ZMod r} P(c·N/r),
-- i.e. the ideal outcomes keep total probability at least 1/(2k+1), with k = A + D.
/-- The complete general-k theorem for literal edits and the QFT formula. -/
theorem general_k_glitches {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) :
    (1 : ℝ) / (2 * (A + D : ℕ) + 1) ≤
      ∑ c : ZMod r,
        fourierProbability (fun x : Fin N => symbol w v (s x)) (c.val * (N / r)) := by
  -- `literal_universal` proves this for `idealProbability`. Rewriting each term with
  -- `idealProbability_eq_fourier` gives exactly this statement.
  simpa only [idealProbability_eq_fourier s w v N hN hrN] using
    literal_universal h w v N hN

/-- Actual ideal outcomes, each counted once:
the set `G = {c·(N/r) : c = 0, ..., r-1}` of natural numbers. -/
noncomputable def goodOutcomes (N r : ℕ) [NeZero r] : Finset ℕ :=
  univ.image (fun c : ZMod r => c.val * (N / r))

/-- Distinct `c` give distinct outcomes `c·(N/r)`, so no outcome is counted twice. -/
theorem good_index_injective (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) :
    Function.Injective (fun c : ZMod r => c.val * (N / r)) := by
  -- N/r > 0, because r ≤ N (r divides N > 0) and r > 0.
  have hm : 0 < N / r := Nat.div_pos (Nat.le_of_dvd hN hrN) (NeZero.pos r)
  -- Suppose a·(N/r) = b·(N/r); that is `hab`.
  intro a b hab
  -- It is enough to show a.val = b.val (`ZMod.val_injective`) ...
  apply ZMod.val_injective
  -- ... which follows by cancelling the nonzero factor N/r.
  exact mul_right_cancel₀ (Nat.ne_of_gt hm) hab

/-- Every ideal outcome is a valid outcome: `c·(N/r) < N`. -/
theorem good_index_lt (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) (c : ZMod r) :
    c.val * (N / r) < N := by
  have hm : 0 < N / r := Nat.div_pos (Nat.le_of_dvd hN hrN) (NeZero.pos r)
  calc
    -- c.val < r (`ZMod.val_lt`), so c.val·(N/r) < r·(N/r).
    _ < r * (N / r) := Nat.mul_lt_mul_of_pos_right (ZMod.val_lt c) hm
    _ = N := by rw [Nat.mul_comm, Nat.div_mul_cancel hrN]

-- `sum_image`: if f is injective, ∑_{y ∈ f(S)} P y = ∑_{c ∈ S} P (f c). Its argument is the
-- injectivity proof, written as a function of a, a proof that a ∈ S, b, a proof that b ∈ S,
-- and `hab : f a = f b`.
/-- Summing over the set `G` is the same as summing over `c ∈ ZMod r`. -/
theorem sum_goodOutcomes (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) (P : ℕ → ℝ) :
    (∑ y ∈ goodOutcomes N r, P y) = ∑ c : ZMod r, P (c.val * (N / r)) := by
  unfold goodOutcomes
  exact sum_image (fun a _ b _ hab => good_index_injective N hN hrN hab)

-- MAIN THEOREM, sharp form. The setting is the same, but the sum is over the set G and both
-- bounds are given:
--   r / ((A+D+1)·r + I) ≤ ∑_{y ∈ G} P(y)   and   (N - I)² / ((A+D+1)·N²) ≤ ∑_{y ∈ G} P(y).
-- `let P := ...` and `let I := ...` inside the statement are local abbreviations.
-- This `I` (the number of surviving insertions) shadows the imaginary unit `I` from
-- `open Complex`: inside this statement the local name wins.
/-- Both sharper bounds, for the actual set G rather than an indexed multiset. -/
theorem general_k_sharp_bounds {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) :
    let P := fourierProbability (fun x : Fin N => symbol w v (s x))
    let I := (inserteds s N).card
    (r : ℝ) / (((A + D : ℕ) + 1 : ℝ) * r + I) ≤ ∑ y ∈ goodOutcomes N r, P y
    ∧ ((N : ℝ) - I)^2 / (((A + D : ℕ) + 1 : ℝ) * (N : ℝ)^2)
      ≤ ∑ y ∈ goodOutcomes N r, P y := by
  -- `dsimp` unfolds the `let`s in the goal.
  dsimp
  -- Turn the sum over the set G into a sum over c ∈ ZMod r.
  rw [sum_goodOutcomes N hN hrN]
  -- One goal for each side of the ∧.
  constructor
  · simpa only [idealProbability_eq_fourier s w v N hN hrN] using
      literal_total h w v N hN
  · simpa only [idealProbability_eq_fourier s w v N hN hrN] using
      literal_retained h w v N hN

-- `fun _ _ => ...` ignores both of its arguments (the symbol and a membership proof).
/-- QFT probabilities are nonnegative. -/
theorem fourierProbability_nonneg {N : ℕ} (f : Fin N → W) (y : ℕ) :
    0 ≤ fourierProbability f y := by
  unfold fourierProbability
  exact mul_nonneg (by positivity) (sum_nonneg fun _ _ => normSq_nonneg _)

-- For example, H could be a union of windows around the ideal outcomes. The windows may
-- overlap, because H is a set and each outcome in it is counted once.
/-- Any neighbourhood containing G inherits the bound, even when its windows overlap. -/
theorem general_k_neighbourhood {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) (hrN : r ∣ N)
    (H : Finset ℕ) (hH : goodOutcomes N r ⊆ H) :
    (1 : ℝ) / (2 * (A + D : ℕ) + 1) ≤
      ∑ y ∈ H, fourierProbability (fun x : Fin N => symbol w v (s x)) y := by
  -- Start from the main theorem ...
  have hb := general_k_glitches h w v N hN hrN
  -- ... rewrite its sum over c as a sum over the set G ...
  rw [← sum_goodOutcomes N hN hrN] at hb
  -- ... and enlarge G to H. Probabilities are ≥ 0, so the sum can only grow.
  exact hb.trans (mass_mono _ (fourierProbability_nonneg _) hH)

-- With I = 0 (no inserted entry survives in the first N positions, e.g. deletions only), the
-- bound improves to 1/(A+D+1).
/-- Deletions only, or any history whose inserted entries all disappear from the prefix. -/
theorem general_k_no_surviving_insertions {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) (hrN : r ∣ N)
    (hI : (inserteds s N).card = 0) :
    (1 : ℝ) / ((A + D : ℕ) + 1) ≤
      ∑ y ∈ goodOutcomes N r, fourierProbability (fun x : Fin N => symbol w v (s x)) y := by
  -- `.1` takes the first half of the ∧ in `general_k_sharp_bounds`.
  have hb := (general_k_sharp_bounds h w v N hN hrN).1
  -- Substitute I = 0 and simplify r/((k+1)r + 0) to 1/(k+1) (`no_insertions`).
  simpa only [hI, Nat.cast_zero, no_insertions r (A + D) (NeZero.pos r)] using hb

-- If A + D ≤ N, the ideal outcomes have total probability ≥ (1 - k/N)²/(k+1), where
-- k = A + D. As k/N → 0 this tends to 1/(k+1) (see `retained_factor_tendsto`).
/-- The finite-N version of the "k small compared with N" corollary. -/
theorem general_k_few_edits {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) (hrN : r ∣ N)
    (hkN : A + D ≤ N) :
    (1 - ((A + D : ℕ) : ℝ) / N)^2 / ((A + D : ℕ) + 1) ≤
      ∑ y ∈ goodOutcomes N r, fourierProbability (fun x : Fin N => symbol w v (s x)) y := by
  -- I ≤ A ≤ A + D.
  have hI := (h.inserted_card_le N).trans (Nat.le_add_right A D)
  -- (1 - k/N)²/(k+1) ≤ (N - I)²/((k+1)N²)   by `few_edits`
  --                  ≤ ∑_{y ∈ G} P(y)        by the second half (`.2`) of the sharp bounds.
  exact (few_edits N (A + D) (inserteds s N).card hN hI hkN).trans
    (general_k_sharp_bounds h w v N hN hrN).2

end Glitchy
