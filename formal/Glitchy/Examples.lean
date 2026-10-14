import Glitchy.ModelI

-- The worked example from the report's figure: a period-3 sequence, with N = 15 positions,
-- after two insertions and one deletion (so k = 3). The numbers are checked by computation.
--
-- Lean ideas introduced here:
--   * `decide`: Lean evaluates a finite, computable statement, and the kernel re-checks it
--   * `norm_num`: arithmetic with explicit numbers
--   * `![...]`: vector and matrix literals
--   * a proof built by applying constructors directly (`example_edited`)

namespace Glitchy
open Finset

-- `![a, b, c]` is a vector, i.e. a function Fin 3 → ℕ. Nesting gives a 5 × 3 table.
-- Rows are symbols: 0, 1, 2 are the periodic symbols and 3, 4 are the two inserted symbols.
-- Columns are positions mod 3. Each entry counts the positions in that box.
/-- The value-residue count table in the three-glitch figure. -/
def exampleCounts : Fin 5 → Fin 3 → ℕ :=
  ![![1, 2, 1], ![1, 1, 3], ![2, 1, 1], ![1, 0, 0], ![0, 1, 0]]

-- `by decide` evaluates both sides and checks that they agree.
/-- The table counts all 15 positions. -/
theorem example_count_total : (∑ w, ∑ b, exampleCounts w b) = 15 := by decide
/-- The sum of the squared entries is 25. -/
theorem example_count_squares : (∑ w, ∑ b, (exampleCounts w b)^2) = 25 := by decide

/-- The total ideal-outcome probability: `(r/N²) · ∑ (count)² = 3/15² · 25 = 1/3`. -/
theorem example_mass :
    (3 : ℝ) / 15^2 * ∑ w, ∑ b, ((exampleCounts w b : ℝ))^2 = 1 / 3 := by
  -- The previous theorem works in ℕ; `exact_mod_cast` moves it into ℝ.
  have hc : (∑ w, ∑ b, ((exampleCounts w b : ℝ))^2) = (25 : ℝ) := by
    exact_mod_cast example_count_squares
  rw [hc]
  -- 3/225 · 25 = 1/3.
  norm_num

-- With k = 3, r = 3, I = 2: the k-only bound is 1/(2k+1) = 1/7, the sharp bound is
-- r/((k+1)r + I) = 3/14, and the true value is 1/3.
/-- The bounds are consistent with the example: `1/7 ≤ 3/14 ≤ 1/3`. -/
theorem example_bounds : (1 : ℝ) / 7 ≤ 3 / 14 ∧ (3 : ℝ) / 14 ≤ 1 / 3 := by norm_num

/-- The old two-glitch universal constant is a corollary of the new one: `k = 2` gives `1/5`. -/
theorem two_glitch_constant : (1 : ℝ) / (2 * 2 + 1) = 1 / 5 := by norm_num

-- Two equal contributions with opposite phase cancel. A single ideal outcome can therefore
-- have probability 0 even though the total over all ideal outcomes is positive.
/-- Cancellation at one exact outcome does not contradict positive total mass. -/
theorem one_outcome_cancels : (1 / 2 : ℂ) + (1 / 2 : ℂ) * (-1) = 0 := by norm_num

-- Read the nested calls from the inside out. Start from the unedited sequence x ↦ x, insert
-- at position 3 (first insertion, A = 0, tag -1), insert at position 7 (A = 1, tag -2),
-- and finally delete position 11.
/-- The literal history depicted in the figure, with unique insertion tags -1 and -2. -/
def exampleSequence : ℕ → ℤ :=
  deleteAt (insertAt (insertAt (fun x => (x : ℤ)) 3 0) 7 1) 11

-- The proof is literally the edit history: the constructors `original`, `insert`, `insert`,
-- `delete`, applied in order. Lean checks that the result has type `Edited 2 1 ...`:
-- A goes 0 → 1 → 2 and D goes 0 → 1.
/-- The example sequence is the result of 2 insertions and 1 deletion. -/
theorem example_edited : Edited 2 1 exampleSequence :=
  Edited.delete (Edited.insert (Edited.insert Edited.original 3) 7) 11

-- An element of `Fin 5` is a pair ⟨value, proof that value < 5⟩. The proof chains
-- b.val < 3 (`ZMod.val_lt`) and 3 < 5 (`by decide`) with `Nat.lt_trans`.
/-- Periodic symbols: residue `b` shows symbol `b ∈ {0, 1, 2}`. -/
def exampleAlphabet (b : ZMod 3) : Fin 5 :=
  ⟨b.val, Nat.lt_trans (ZMod.val_lt b) (by decide)⟩

/-- Inserted symbols: tag `-1` shows symbol `3`; any other tag (here `-2`) shows `4`. -/
def exampleInserted (t : ℤ) : Fin 5 := if t = -1 then 3 else 4

-- `revert a b` turns the goal back into `∀ a b, ...`. That is a statement about finitely many
-- cases (5 × 3), so `decide` can check all of them by computation.
/-- Kernel evaluation checks that the actual edited prefix gives exactly the pictured table. -/
theorem example_counts_correct (a : Fin 5) (b : ZMod 3) :
    ((univ : Finset (Fin 15)).filter
      (fun x => box exampleSequence exampleAlphabet exampleInserted x = (a,b))).card =
      exampleCounts a ⟨b.val, ZMod.val_lt b⟩ := by
  revert a b
  decide

/-- Exactly two inserted entries survive in the first 15 positions (`I = 2`). -/
theorem example_inserted_count : (inserteds exampleSequence 15).card = 2 := by decide

/-- The same table, stated with `prefixCounts` (real numbers). -/
theorem example_prefix_counts (a : Fin 5) (b : ZMod 3) :
    prefixCounts exampleSequence exampleAlphabet exampleInserted 15 a b =
      (exampleCounts a ⟨b.val, ZMod.val_lt b⟩ : ℝ) := by
  unfold prefixCounts fiberCounts
  exact_mod_cast example_counts_correct a b

-- The three ideal outcomes are y = 0, 5, 10 (that is, c·5 for c = 0, 1, 2).
/-- Exact Fourier mass for the actual three-edit sequence in the illustration. -/
theorem example_fourier_mass :
    (∑ c : ZMod 3, fourierProbability
      (fun x : Fin 15 => symbol exampleAlphabet exampleInserted (exampleSequence x))
      (c.val * 5)) = (1 : ℝ) / 3 := by
  -- The bridge theorem for this sequence with N = 15 and r = 3. The two `by decide`s check
  -- 0 < 15 and 3 ∣ 15.
  have he := idealProbability_eq_fourier exampleSequence exampleAlphabet exampleInserted
    15 (by decide) (by decide)
  -- Flip it (`.symm`) so it rewrites `fourierProbability` into `idealProbability`.
  -- Lean accepts `c.val * 5` for `c.val * (15 / 3)` because the two are equal by computation.
  have he' (c : ZMod 3) :
      fourierProbability
        (fun x : Fin 15 => symbol exampleAlphabet exampleInserted (exampleSequence x))
        (c.val * 5) =
      idealProbability 15 (prefixCounts exampleSequence exampleAlphabet exampleInserted 15) c :=
    (he c).symm
  simp_rw [he']
  -- The total ideal mass is (3/15²) · ∑ (count)².
  rw [residue_mass]
  -- Replace the counts with the table entries.
  simp_rw [example_prefix_counts]
  -- `change` restates the goal in a form that is equal by definition. (`ZMod 3` is defined to
  -- be `Fin 3`.) The new form matches `example_mass` exactly.
  change (3 : ℝ) / 15^2 * ∑ a : Fin 5, ∑ b : Fin 3, ((exampleCounts a b : ℝ))^2 = _
  exact example_mass

end Glitchy
