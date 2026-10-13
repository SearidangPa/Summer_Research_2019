import Glitchy.ModelI

namespace Glitchy
open Finset

/-- The value-residue count table in the three-glitch figure. -/
def exampleCounts : Fin 5 → Fin 3 → ℕ :=
  ![![1, 2, 1], ![1, 1, 3], ![2, 1, 1], ![1, 0, 0], ![0, 1, 0]]

theorem example_count_total : (∑ w, ∑ b, exampleCounts w b) = 15 := by decide
theorem example_count_squares : (∑ w, ∑ b, (exampleCounts w b)^2) = 25 := by decide

theorem example_mass :
    (3 : ℝ) / 15^2 * ∑ w, ∑ b, ((exampleCounts w b : ℝ))^2 = 1 / 3 := by
  have hc : (∑ w, ∑ b, ((exampleCounts w b : ℝ))^2) = (25 : ℝ) := by
    exact_mod_cast example_count_squares
  rw [hc]
  norm_num

theorem example_bounds : (1 : ℝ) / 7 ≤ 3 / 14 ∧ (3 : ℝ) / 14 ≤ 1 / 3 := by norm_num

/-- The old two-glitch universal constant is a corollary of the new one. -/
theorem two_glitch_constant : (1 : ℝ) / (2 * 2 + 1) = 1 / 5 := by norm_num

/-- Cancellation at one exact outcome does not contradict positive total mass. -/
theorem one_outcome_cancels : (1 / 2 : ℂ) + (1 / 2 : ℂ) * (-1) = 0 := by norm_num

/-- The literal history depicted in the figure, with unique insertion tags -1 and -2. -/
def exampleSequence : ℕ → ℤ :=
  deleteAt (insertAt (insertAt (fun x => (x : ℤ)) 3 0) 7 1) 11

theorem example_edited : Edited 2 1 exampleSequence :=
  Edited.delete (Edited.insert (Edited.insert Edited.original 3) 7) 11

def exampleAlphabet (b : ZMod 3) : Fin 5 :=
  ⟨b.val, Nat.lt_trans (ZMod.val_lt b) (by decide)⟩

def exampleInserted (t : ℤ) : Fin 5 := if t = -1 then 3 else 4

/-- Kernel evaluation checks that the actual edited prefix gives exactly the pictured table. -/
theorem example_counts_correct (a : Fin 5) (b : ZMod 3) :
    ((univ : Finset (Fin 15)).filter
      (fun x => box exampleSequence exampleAlphabet exampleInserted x = (a,b))).card =
      exampleCounts a ⟨b.val, ZMod.val_lt b⟩ := by
  revert a b
  decide

theorem example_inserted_count : (inserteds exampleSequence 15).card = 2 := by decide

theorem example_prefix_counts (a : Fin 5) (b : ZMod 3) :
    prefixCounts exampleSequence exampleAlphabet exampleInserted 15 a b =
      (exampleCounts a ⟨b.val, ZMod.val_lt b⟩ : ℝ) := by
  unfold prefixCounts fiberCounts
  exact_mod_cast example_counts_correct a b

/-- Exact Fourier mass for the actual three-edit sequence in the illustration. -/
theorem example_fourier_mass :
    (∑ c : ZMod 3, fourierProbability
      (fun x : Fin 15 => symbol exampleAlphabet exampleInserted (exampleSequence x))
      (c.val * 5)) = (1 : ℝ) / 3 := by
  have he := idealProbability_eq_fourier exampleSequence exampleAlphabet exampleInserted
    15 (by decide) (by decide)
  have he' (c : ZMod 3) :
      fourierProbability
        (fun x : Fin 15 => symbol exampleAlphabet exampleInserted (exampleSequence x))
        (c.val * 5) =
      idealProbability 15 (prefixCounts exampleSequence exampleAlphabet exampleInserted 15) c :=
    (he c).symm
  simp_rw [he']
  rw [residue_mass]
  simp_rw [example_prefix_counts]
  change (3 : ℝ) / 15^2 * ∑ a : Fin 5, ∑ b : Fin 3, ((exampleCounts a b : ℝ))^2 = _
  exact example_mass

end Glitchy
