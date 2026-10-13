import Glitchy.Literal

namespace Glitchy
open Finset Complex
variable {r k : ℕ} [NeZero r]

noncomputable def grouped (shift : Fin (k + 1) → ZMod r)
    (length : Fin (k + 1) → ℝ) (b : ZMod r) : ℝ :=
  ∑ j ∈ univ.filter (fun j => shift j = b), length j

noncomputable def modelMass (shift : Fin (k + 1) → ZMod r)
    (length : Fin (k + 1) → ℝ) : ℝ :=
  (1 / (r : ℝ)) * ∑ c : ZMod r,
    normSq (∑ j, (length j : ℂ) * ZMod.stdAddChar (c * shift j))

theorem grouped_sum (shift : Fin (k + 1) → ZMod r) (length : Fin (k + 1) → ℝ) :
    (∑ b, grouped shift length b) = ∑ j, length j :=
  sum_fiberwise univ shift length

omit [NeZero r] in
theorem grouped_zero (shift : Fin (k + 1) → ZMod r) (length : Fin (k + 1) → ℝ)
    (b : ZMod r) (hb : b ∉ univ.image shift) : grouped shift length b = 0 := by
  have he : univ.filter (fun j => shift j = b) = ∅ := by
    apply eq_empty_iff_forall_not_mem.mpr
    intro j hj
    exact hb (mem_image.mpr ⟨j, mem_univ _, (mem_filter.mp hj).2⟩)
  simp [grouped, he]

theorem grouped_support (shift : Fin (k + 1) → ZMod r) (length : Fin (k + 1) → ℝ) :
    (univ.filter (fun b => grouped shift length b ≠ 0)).card ≤ k + 1 := by
  have hs : univ.filter (fun b => grouped shift length b ≠ 0) ⊆ univ.image shift := by
    intro b hb
    by_contra hn
    exact (mem_filter.mp hb).2 (grouped_zero shift length b hn)
  exact (card_le_card hs).trans (by simpa using (card_image_le (s := univ) (f := shift)))

theorem grouped_amplitude (shift : Fin (k + 1) → ZMod r) (length : Fin (k + 1) → ℝ)
    (c : ZMod r) :
    (∑ b, (grouped shift length b : ℂ) * ZMod.stdAddChar (c * b)) =
      ∑ j, (length j : ℂ) * ZMod.stdAddChar (c * shift j) := by
  unfold grouped
  push_cast
  simp_rw [sum_mul]
  have he (b : ZMod r) :
      (∑ j ∈ univ.filter (fun j => shift j = b),
        (length j : ℂ) * ZMod.stdAddChar (c * b)) =
      ∑ j ∈ univ.filter (fun j => shift j = b),
        (length j : ℂ) * ZMod.stdAddChar (c * shift j) := by
    apply sum_congr rfl
    intro j hj
    rw [(mem_filter.mp hj).2]
  simp_rw [he]
  exact sum_fiberwise univ shift _

/-- Exact Model I total mass as the squared grouped lengths. -/
theorem modelMass_eq (shift : Fin (k + 1) → ZMod r) (length : Fin (k + 1) → ℝ) :
    modelMass shift length = ∑ b, (grouped shift length b)^2 := by
  unfold modelMass
  simp_rw [← grouped_amplitude, residue_parseval]
  have hr : (r : ℝ) ≠ 0 := by exact_mod_cast (NeZero.ne r)
  field_simp

theorem modelI_bound (shift : Fin (k + 1) → ZMod r) (length : Fin (k + 1) → ℝ)
    (hlength : ∑ j, length j = 1) :
    (1 : ℝ) / (k + 1) ≤ modelMass shift length := by
  rw [modelMass_eq]
  exact grouped_lengths _ k ((grouped_sum shift length).trans hlength)
    (grouped_support shift length)

/-- Equal blocks with distinct shift residues attain the exact-set bound. -/
theorem modelI_sharp (shift : Fin (k + 1) → ZMod r) (hi : Function.Injective shift) :
    modelMass shift (fun _ => (1 : ℝ) / (k + 1)) = (1 : ℝ) / (k + 1) := by
  rw [modelMass_eq]
  have hp (j : Fin (k + 1)) :
      grouped shift (fun _ => (1 : ℝ) / (k + 1)) (shift j) = (1 : ℝ) / (k + 1) := by
    simp [grouped, hi.eq_iff, Finset.filter_eq']
  have hs :
      (∑ b ∈ univ.image shift, (grouped shift (fun _ => (1 : ℝ) / (k + 1)) b)^2) =
      ∑ b, (grouped shift (fun _ => (1 : ℝ) / (k + 1)) b)^2 := by
    apply sum_subset (subset_univ _)
    intro b hb hn
    simp [grouped_zero shift _ b hn]
  rw [← hs, sum_image (fun a _ b _ hab => hi hab)]
  simp_rw [hp]
  exact equal_grouped_lengths k

/-- Equal-sign shifts 0,...,k are distinct modulo r=k+1. -/
def canonicalShift (k : ℕ) (j : Fin (k + 1)) : ZMod (k + 1) := (j.val : ZMod (k + 1))

theorem canonicalShift_injective (k : ℕ) : Function.Injective (canonicalShift k) := by
  intro a b hab
  apply Fin.ext
  have hv := congrArg ZMod.val hab
  simpa [canonicalShift, ZMod.val_natCast, Nat.mod_eq_of_lt a.isLt,
    Nat.mod_eq_of_lt b.isLt] using hv

/-- Consequently there is no positive k-independent bound on the exact set G in Model I. -/
theorem no_uniform_exact_bound (eps : ℝ) (heps : 0 < eps) :
    ∃ k : ℕ, modelMass (canonicalShift k) (fun _ => (1 : ℝ) / (k + 1)) < eps := by
  obtain ⟨k, hk⟩ := exists_nat_gt (1 / eps)
  refine ⟨k, ?_⟩
  rw [modelI_sharp _ (canonicalShift_injective k)]
  have hkp : (0 : ℝ) < (k : ℝ) + 1 := by positivity
  have hk' : (1 : ℝ) < (k : ℝ) * eps := (div_lt_iff₀ heps).mp hk
  apply (div_lt_iff₀ hkp).mpr
  nlinarith

end Glitchy
