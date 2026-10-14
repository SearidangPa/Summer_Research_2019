import Glitchy.Matching

/-!
Literal-edit bookkeeping for matching moments. The estimate deliberately uses
A (all insertions), rather than the sharper surviving-prefix count I. Since
A ≤ A+D, it gives exactly the normalized estimate used in Section 4.2.
-/
namespace Glitchy
open Finset

/-- Starts of pairs whose tags are original and differ by their position gap. -/
def goodPairs (s : ℕ → ℤ) (h L : ℕ) : Finset ℕ :=
  (range L).filter fun x => 0 ≤ s x ∧ s (x+h) = s x + (h : ℤ)

private theorem discard_crossings (S : Finset ℕ) (p h : ℕ) :
    S.card ≤ (S.filter (fun x => x+h < p ∨ p ≤ x)).card + h := by
  have hb : (S.filter (fun x => ¬(x+h < p ∨ p ≤ x))) ⊆ Ico (p-h) p := by
    intro x hx
    simp only [mem_filter, not_or, not_lt, not_le, mem_Ico] at hx ⊢
    omega
  have hc := card_le_card hb
  rw [Nat.card_Ico] at hc
  have hs := filter_card_add_filter_neg_card_eq_card (s := S) (fun x => x+h < p ∨ p ≤ x)
  omega

private theorem discard_deleted (S : Finset ℕ) (p h : ℕ) :
    S.card ≤ (S.filter (fun x => x+h < p ∨ p < x)).card + (h+1) := by
  have hb : (S.filter (fun x => ¬(x+h < p ∨ p < x))) ⊆ Icc (p-h) p := by
    intro x hx
    simp only [mem_filter, not_or, not_lt, mem_Icc] at hx ⊢
    omega
  have hc := card_le_card hb
  rw [Nat.card_Icc] at hc
  have hs := filter_card_add_filter_neg_card_eq_card (s := S) (fun x => x+h < p ∨ p < x)
  omega

/-- Inserting a symbol loses at most h old matching pairs; the prefix grows by one. -/
theorem goodPairs_insert (s : ℕ → ℤ) (p A h L : ℕ) :
    (goodPairs s h L).card ≤ (goodPairs (insertAt s p A) h (L+1)).card + h := by
  let S := (goodPairs s h L).filter (fun x => x+h < p ∨ p ≤ x)
  let shift := fun x : ℕ => if x < p then x else x+1
  have hinj : Function.Injective shift := by
    intro x y he
    dsimp [shift] at he
    split_ifs at he <;> omega
  have hsub : S.image shift ⊆ goodPairs (insertAt s p A) h (L+1) := by
    intro y hy
    obtain ⟨x, hx, rfl⟩ := mem_image.mp hy
    have hx' := mem_filter.mp hx
    have hg := mem_filter.mp hx'.1
    have hxL := mem_range.mp hg.1
    rcases hx'.2 with hb | ha
    · have hxp : x < p := by omega
      simp only [goodPairs, mem_filter, mem_range, shift, if_pos hxp]
      refine ⟨by omega, ?_, ?_⟩
      · simpa [insertAt, hxp] using hg.2.1
      · simpa [insertAt, hxp, hb] using hg.2.2
    · have hxp : ¬x < p := by omega
      have h1 : ¬x+1 < p := by omega
      have h2 : x+1 ≠ p := by omega
      have h3 : ¬x+1+h < p := by omega
      have h4 : x+1+h ≠ p := by omega
      simp only [goodPairs, mem_filter, mem_range, shift, if_neg hxp]
      refine ⟨by omega, ?_, ?_⟩
      · simpa [insertAt, h1, h2] using hg.2.1
      · simpa [insertAt, h1, h2, h3, h4, show x+1+h-1=x+h by omega] using hg.2.2
  have hc := card_le_card hsub
  rw [card_image_of_injective _ hinj] at hc
  have hd := discard_crossings (goodPairs s h L) p h
  change (goodPairs s h L).card ≤ S.card + h at hd
  omega

/-- Deleting a symbol loses at most h+1 old pairs; the prefix shrinks by one. -/
theorem goodPairs_delete (s : ℕ → ℤ) (p h L : ℕ) :
    (goodPairs s h (L+1)).card ≤ (goodPairs (deleteAt s p) h L).card + (h+1) := by
  let S := (goodPairs s h (L+1)).filter (fun x => x+h < min p L ∨ min p L < x)
  let shift := fun x : ℕ => if x < p then x else x-1
  have hinj : Set.InjOn shift S := by
    intro x hx y hy he
    have hx' := (mem_filter.mp hx).2
    have hy' := (mem_filter.mp hy).2
    have hxL := mem_range.mp (mem_filter.mp (mem_filter.mp hx).1).1
    have hyL := mem_range.mp (mem_filter.mp (mem_filter.mp hy).1).1
    dsimp [shift] at he
    split_ifs at he <;> omega
  have hsub : S.image shift ⊆ goodPairs (deleteAt s p) h L := by
    intro y hy
    obtain ⟨x, hx, rfl⟩ := mem_image.mp hy
    have hx' := mem_filter.mp hx
    have hg := mem_filter.mp hx'.1
    have hxL := mem_range.mp hg.1
    rcases hx'.2 with hb | ha
    · have hxp : x < p := by omega
      have hbp : x+h < p := by omega
      simp only [goodPairs, mem_filter, mem_range, shift, if_pos hxp]
      refine ⟨?_, ?_, ?_⟩
      · omega
      · simpa [deleteAt, hxp] using hg.2.1
      · simpa [deleteAt, hxp, hbp] using hg.2.2
    · have hxp : ¬x < p := by omega
      have h1 : ¬x-1 < p := by omega
      have h2 : ¬x-1+h < p := by omega
      simp only [goodPairs, mem_filter, mem_range, shift, if_neg hxp]
      refine ⟨by omega, ?_, ?_⟩
      · simpa [deleteAt, h1, show x-1+1=x by omega] using hg.2.1
      · simpa [deleteAt, h1, h2, show x-1+1=x by omega,
          show x-1+h+1=x+h by omega] using hg.2.2
  have hc := card_le_card hsub
  rw [card_image_of_injOn hinj] at hc
  have hd := discard_deleted (goodPairs s h (L+1)) (min p L) h
  change (goodPairs s h (L+1)).card ≤ S.card + (h+1) at hd
  omega

/-- The matching-pair count follows from the edit history, not from a
piecewise-periodicity assumption. It holds for every prefix and every lag. -/
theorem Edited.goodPairs_bound {A D s} (hedit : Edited A D s) (h L : ℕ) :
    L ≤ (goodPairs s h L).card + (A+D)*h + A := by
  induction hedit generalizing L with
  | original => simp [goodPairs]
  | @insert A D s hedit p ih =>
    cases L with
    | zero => omega
    | succ L =>
      have hcount := goodPairs_insert s p A h L
      have hprev := ih L
      nlinarith
  | @delete A D s hedit p ih =>
    have hcount := goodPairs_delete s p h L
    have hprev := ih (L+1)
    nlinarith

/-- Discarding all wrapped pairs costs at most h additional positions. -/
theorem Edited.prefix_goodPairs_bound {A D s} (hedit : Edited A D s) (h N : ℕ) :
    N ≤ (goodPairs s h (N-h)).card + (A+D+1)*h + A := by
  have := hedit.goodPairs_bound h (N-h)
  have : N ≤ (N-h)+h := by omega
  nlinarith [hedit.goodPairs_bound h (N-h)]

section Symbols
variable {W : Type*} [Fintype W] [DecidableEq W]
    {r N : ℕ} [NeZero r] [NeZero N]

/-- The literal prefix viewed as a cyclic function; wrapping is only used when
forming correlations, not when performing edits. -/
noncomputable def literalCycle (s : ℕ → ℤ) (w : ZMod r → W) (v : ℤ → W) : ZMod N → W :=
  fun x => symbol w v (s x.val)

omit [Fintype W] [DecidableEq W] [NeZero r] in
/-- Safe original pairs at a multiple-of-period lag have equal output symbols. -/
theorem goodPairs_symbols (s : ℕ → ℤ) (w : ZMod r → W) (v : ℤ → W)
    (h L x : ℕ) (hrh : r ∣ h) (hx : x ∈ goodPairs s h L) :
    symbol w v (s x) = symbol w v (s (x+h)) := by
  have hg := (mem_filter.mp hx).2
  have hs : 0 ≤ s (x+h) := by rw [hg.2]; exact add_nonneg hg.1 (Int.ofNat_nonneg h)
  have hz : ((h : ℤ) : ZMod r) = 0 := by
    simpa using (ZMod.natCast_zmod_eq_zero_iff_dvd h r).mpr hrh
  simp only [symbol, if_pos hg.1, if_pos hs]
  rw [hg.2, Int.cast_add, hz, add_zero]

omit [Fintype W] [NeZero r] in
/-- The nondivisible matching estimate, derived from actual literal edits.
Using A instead of I is weaker but suffices for the same normalized bound. -/
theorem Edited.matching_bound {A D s} (hedit : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (h : ℕ) (hrh : r ∣ h) :
    1 - (((A+D+1 : ℕ) : ℝ) * h + A) / N ≤
      matchingFraction (literalCycle (N := N) s w v) (h : ZMod N) := by
  classical
  let S := goodPairs s h (N-h)
  let C := (univ : Finset (ZMod N)).filter fun x =>
    literalCycle s w v x = literalCycle s w v (x + (h : ZMod N))
  have hxlt (x : ℕ) (hx : x ∈ S) : x+h < N := by
    have := mem_range.mp (mem_filter.mp hx).1
    omega
  have hsub : S.image (fun x : ℕ => (x : ZMod N)) ⊆ C := by
    intro z hz
    obtain ⟨x, hx, rfl⟩ := mem_image.mp hz
    have hxN : x < N := by have := hxlt x hx; omega
    have hxh := hxlt x hx
    simp only [C, mem_filter, mem_univ, true_and, literalCycle]
    rw [ZMod.val_natCast_of_lt hxN, ← Nat.cast_add, ZMod.val_natCast_of_lt hxh]
    exact goodPairs_symbols s w v h (N-h) x hrh hx
  have hinj : Set.InjOn (fun x : ℕ => (x : ZMod N)) S := by
    intro x hx y hy he
    have hxN : x < N := by have := hxlt x hx; omega
    have hyN : y < N := by have := hxlt y hy; omega
    have := congrArg ZMod.val he
    simpa [ZMod.val_natCast_of_lt hxN, ZMod.val_natCast_of_lt hyN] using this
  have hc := card_le_card hsub
  rw [card_image_of_injOn hinj] at hc
  have hn := hedit.prefix_goodPairs_bound h N
  change N ≤ S.card + (A+D+1)*h + A at hn
  have hnat : N ≤ C.card + (A+D+1)*h + A := by omega
  have hreal : (N : ℝ) ≤ C.card + ((A+D+1 : ℕ) : ℝ)*h + A := by exact_mod_cast hnat
  have hN : (0 : ℝ) < N := by exact_mod_cast NeZero.pos N
  have hmatch : matchingFraction (literalCycle (N := N) s w v) (h : ZMod N) = C.card / N := by
    simp [matchingFraction, C, div_eq_mul_inv, mul_comm]
  rw [hmatch]
  apply (le_div_iff₀ hN).mpr
  have heq : (1 - (((A+D+1 : ℕ) : ℝ)*h+A)/N) * N = N - (((A+D+1 : ℕ) : ℝ)*h+A) := by
    field_simp
  rw [heq]
  linarith
end Symbols

end Glitchy
