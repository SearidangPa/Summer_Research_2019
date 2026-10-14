import Glitchy.Literal

/-!
Fourier matching moments on an arbitrary finite cyclic domain. Unlike the
residue-count identity, these results do not require a period dividing N.
-/

namespace Glitchy
open Finset Complex

variable {N : ℕ} [NeZero N] {W : Type*} [Fintype W] [DecidableEq W]

/-- The QFT probability, with both positions and outcomes indexed by `ZMod N`. -/
noncomputable def cyclicProbability (f : ZMod N → W) (y : ZMod N) : ℝ :=
  (1 / (N : ℝ)^2) * ∑ a : W,
    normSq (∑ x : ZMod N, if f x = a then ZMod.stdAddChar (y * x) else 0)

/-- Fraction of positions whose symbol agrees with its circular translate. -/
noncomputable def matchingFraction (f : ZMod N → W) (h : ZMod N) : ℝ :=
  (1 / (N : ℝ)) * ∑ x : ZMod N, if f x = f (x + h) then (1 : ℝ) else 0

/-- Character orthogonality in a form convenient for correlation identities. -/
theorem sum_stdAddChar_mul (a : ZMod N) :
    (∑ y : ZMod N, ZMod.stdAddChar (y * a)) = if a = 0 then (N : ℂ) else 0 := by
  simpa [eq_comm] using (character_orthogonality (0 : ZMod N) a)

/-- Expand a squared amplitude into equal-symbol pairs. -/
theorem cyclicProbability_expand (f : ZMod N → W) (y : ZMod N) :
    (cyclicProbability f y : ℂ) = (1 / (N : ℂ)^2) *
      ∑ x : ZMod N, ∑ z : ZMod N,
        if f x = f z then ZMod.stdAddChar (y * (z - x)) else 0 := by
  unfold cyclicProbability
  push_cast
  congr 1
  simp_rw [normSq_eq_conj_mul_self, map_sum, sum_mul, mul_sum]
  rw [sum_comm]
  apply sum_congr rfl
  intro x _
  rw [sum_comm]
  apply sum_congr rfl
  intro z _
  simp only [apply_ite (starRingEnd ℂ), map_zero, ite_mul, mul_ite, zero_mul, mul_zero]
  by_cases hf : f x = f z
  · simp [hf, ← AddChar.map_neg_eq_conj, ← AddChar.map_add_eq_mul,
      show -(y * x) + y * z = y * (z - x) by ring]
  · simp [hf]

/-- Complex Fourier moment equals the real matching fraction. -/
theorem cyclic_matching_moment (f : ZMod N → W) (h : ZMod N) :
    (∑ y : ZMod N, (cyclicProbability f y : ℂ) * ZMod.stdAddChar (y * h)) =
      (matchingFraction f h : ℂ) := by
  simp_rw [cyclicProbability_expand, mul_assoc, sum_mul]
  rw [← mul_sum, sum_comm]
  have inner (x : ZMod N) :
      (∑ y : ZMod N, ∑ z : ZMod N,
        (if f x = f z then ZMod.stdAddChar (y * (z - x)) else 0) *
          ZMod.stdAddChar (y * h)) =
        if f x = f (x - h) then (N : ℂ) else 0 := by
    rw [sum_comm]
    simp_rw [ite_mul, zero_mul, ← AddChar.map_add_eq_mul,
      show ∀ y z : ZMod N, y * (z - x) + y * h = y * (z - x + h) by
        intros; ring, sum_ite_irrel, sum_stdAddChar_mul]
    have hz (z : ZMod N) : z - x + h = 0 ↔ z = x - h := by constructor <;> intro hz <;> linear_combination hz
    simp_rw [hz]
    simp only [sum_const_zero]
    have swap (z : ZMod N) :
        (if f x = f z then if z = x - h then (N : ℂ) else 0 else 0) =
        (if z = x - h then if f x = f z then (N : ℂ) else 0 else 0) := by
      split_ifs <;> rfl
    rw [show (∑ z : ZMod N, if f x = f z then if z = x - h then (N : ℂ) else 0 else 0) =
      ∑ z : ZMod N, if z = x - h then if f x = f z then (N : ℂ) else 0 else 0 from
      sum_congr rfl (fun z _ => swap z)]
    simp
  simp_rw [inner]
  have shift : (∑ x : ZMod N, if f x = f (x - h) then (N : ℂ) else 0) =
      ∑ x : ZMod N, if f x = f (x + h) then (N : ℂ) else 0 := by
    simpa [eq_comm] using
      (Equiv.sum_comp (Equiv.addRight h)
        (fun x : ZMod N => if f x = f (x - h) then (N : ℂ) else 0)).symm
  rw [shift]
  unfold matchingFraction
  push_cast
  simp only [apply_ite Complex.ofReal, Complex.ofReal_one, Complex.ofReal_zero]
  simp_rw [show ∀ p : Prop, [Decidable p] → (if p then (N : ℂ) else 0) =
      (if p then (1 : ℂ) else 0) * N by intros; split_ifs <;> simp]
  rw [← sum_mul]
  have hn : (N : ℂ) ≠ 0 := by exact_mod_cast NeZero.ne N
  field_simp
  ring

/-- Probabilities sum to one even when a proposed period does not divide N. -/
theorem cyclicProbability_sum (f : ZMod N → W) :
    ∑ y : ZMod N, cyclicProbability f y = 1 := by
  have hm := cyclic_matching_moment f 0
  have hn : (N : ℝ) ≠ 0 := by exact_mod_cast NeZero.ne N
  simpa [matchingFraction, ZMod.card, hn] using congrArg Complex.re hm

theorem cyclicProbability_nonneg (f : ZMod N → W) (y : ZMod N) :
    0 ≤ cyclicProbability f y := by
  unfold cyclicProbability
  exact mul_nonneg (by positivity) (sum_nonneg fun _ _ => normSq_nonneg _)

/-- The cyclic indexing is only a reindexing of the literal prefix. -/
noncomputable def cyclePrefix (f : Fin N → W) : ZMod N → W :=
  fun x => f ((ZMod.finEquiv N).symm x)

theorem finEquiv_eq_natCast (x : Fin N) :
    ZMod.finEquiv N x = (x.val : ZMod N) := by
  cases N with
  | zero => exact (NeZero.ne 0 rfl).elim
  | succ n =>
    change x = (x.val : Fin (n + 1))
    apply Fin.ext
    simp

/-- Bridge to the existing complex-exponential definition, for every outcome. -/
theorem cyclicProbability_eq_fourier (f : Fin N → W) (y : ZMod N) :
    cyclicProbability (cyclePrefix f) y = fourierProbability f y.val := by
  unfold cyclicProbability fourierProbability
  congr 1
  apply sum_congr rfl
  intro a _
  congr 1
  rw [← Equiv.sum_comp (ZMod.finEquiv N).toEquiv]
  apply sum_congr rfl
  intro x _
  change (if f ((ZMod.finEquiv N).symm (ZMod.finEquiv N x)) = a
    then ZMod.stdAddChar (y * (ZMod.finEquiv N x)) else 0) = _
  rw [RingEquiv.symm_apply_apply, finEquiv_eq_natCast]
  have hp := ideal_phase N (NeZero.pos N) (dvd_refl N) x y
  simpa using congrArg (fun z : ℂ => if f x = a then z else 0) hp

/-- Real-part form of the matching identity, used by the Fejér argument. -/
theorem cyclic_matching_moment_re (f : ZMod N → W) (h : ZMod N) :
    (∑ y : ZMod N, cyclicProbability f y * (ZMod.stdAddChar (y * h)).re) =
      matchingFraction f h := by
  simpa using congrArg Complex.re (cyclic_matching_moment f h)

end Glitchy
