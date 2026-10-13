import Mathlib

/-! Exact Fourier averaging on the r ideal outcomes, with no asymptotic limit. -/
namespace Glitchy
open Finset Complex
variable {r : ℕ} [NeZero r]

/-- Orthogonality of the standard characters, proved rather than assumed. -/
theorem character_orthogonality (b d : ZMod r) :
    (∑ c : ZMod r, (starRingEnd ℂ) (ZMod.stdAddChar (c * b)) * ZMod.stdAddChar (c * d)) =
      if b = d then (r : ℂ) else 0 := by
  simp_rw [← AddChar.map_neg_eq_conj, ← AddChar.map_add_eq_mul]
  have he (c : ZMod r) : -(c * b) + c * d = (d - b) * c := by ring
  simp_rw [he]
  by_cases h : b = d
  · simp [h, ZMod.card]
  · simp only [h, if_false]
    exact AddChar.sum_eq_zero_of_ne_one
      (ZMod.isPrimitive_stdAddChar r (sub_ne_zero.mpr (Ne.symm h)))

/-- Parseval for one row of real-valued residue counts. -/
theorem residue_parseval (n : ZMod r → ℝ) :
    (∑ c : ZMod r, normSq (∑ b : ZMod r, (n b : ℂ) * ZMod.stdAddChar (c * b))) =
      (r : ℝ) * ∑ b : ZMod r, (n b)^2 := by
  apply Complex.ofReal_injective
  push_cast
  simp_rw [normSq_eq_conj_mul_self, map_sum, map_mul, conj_ofReal,
    sum_mul, mul_sum]
  have rearrange (c b d : ZMod r) :
      (n b : ℂ) * (starRingEnd ℂ) (ZMod.stdAddChar (c * b)) *
        ((n d : ℂ) * ZMod.stdAddChar (c * d)) =
      (n b : ℂ) * n d *
        ((starRingEnd ℂ) (ZMod.stdAddChar (c * b)) * ZMod.stdAddChar (c * d)) := by ring
  simp_rw [rearrange]
  rw [sum_comm]
  apply sum_congr rfl
  intro b hb
  rw [sum_comm]
  simp_rw [← mul_sum, character_orthogonality]
  simp [mul_comm, pow_two]

/-- Probability at the ideal outcome indexed by c, expressed in residue counts. -/
noncomputable def idealProbability {W : Type*} [Fintype W]
    (N : ℕ) (n : W → ZMod r → ℝ) (c : ZMod r) : ℝ :=
  (1 / (N : ℝ)^2) *
    ∑ w, normSq (∑ b : ZMod r, (n w b : ℂ) * ZMod.stdAddChar (c * b))

/-- The exact sum-of-squared-counts identity used in Step 1. -/
theorem residue_mass {W : Type*} [Fintype W] (N : ℕ) (n : W → ZMod r → ℝ) :
    (∑ c : ZMod r, idealProbability N n c) =
      (r : ℝ) / (N : ℝ)^2 * ∑ w, ∑ b, (n w b)^2 := by
  simp only [idealProbability, ← mul_sum]
  rw [sum_comm]
  simp_rw [residue_parseval, ← mul_sum]
  ring

theorem idealProbability_nonneg {W : Type*} [Fintype W]
    (N : ℕ) (n : W → ZMod r → ℝ) (c : ZMod r) :
    0 ≤ idealProbability N n c := by
  unfold idealProbability
  apply mul_nonneg (by positivity)
  exact sum_nonneg (fun w _ => normSq_nonneg _)

end Glitchy
