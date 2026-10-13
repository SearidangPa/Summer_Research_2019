import Glitchy.Edits
import Glitchy.Counting

namespace Glitchy
open Finset

section Fibers
variable {X B : Type*} [DecidableEq B]

/-- Number of entries of T in one box, as a real number. -/
noncomputable def fiberCounts (f : X → B) (T : Finset X) (b : B) : ℝ :=
  ((T.filter (fun x => f x = b)).card : ℝ)

theorem sum_fiberCounts [Fintype B] (f : X → B) (T : Finset X) :
    (∑ b, fiberCounts f T b) = T.card := by
  unfold fiberCounts
  norm_cast
  exact (card_eq_sum_card_fiberwise (fun x _ => mem_univ (f x))).symm

theorem support_fiberCounts [Fintype B] (f : X → B) (T : Finset X) :
    univ.filter (fun b => fiberCounts f T b ≠ 0) = T.image f := by
  ext b
  simp only [mem_filter, mem_univ, true_and, fiberCounts, Nat.cast_ne_zero,
    card_ne_zero, filter_nonempty_iff, mem_image]

theorem fiberCounts_mono (f : X → B) {S T : Finset X} (hST : S ⊆ T) (b : B) :
    0 ≤ fiberCounts f S b ∧ fiberCounts f S b ≤ fiberCounts f T b := by
  constructor
  · unfold fiberCounts; positivity
  · unfold fiberCounts
    exact_mod_cast card_le_card (filter_subset_filter _ hST)

end Fibers

variable {W : Type*} [DecidableEq W] {r : ℕ} [NeZero r]

/-- Tags select original periodic symbols or arbitrary inserted symbols.
Distinctness of original symbols is not needed for this lower bound. -/
def symbol (w : ZMod r → W) (v : ℤ → W) (t : ℤ) : W :=
  if 0 ≤ t then w (t : ZMod r) else v t

def box (s : ℕ → ℤ) (w : ZMod r → W) (v : ℤ → W) {N : ℕ}
    (x : Fin N) : W × ZMod r :=
  (symbol w v (s x), (x.val : ZMod r))

def originals (s : ℕ → ℤ) (N : ℕ) : Finset (Fin N) :=
  univ.filter (fun x : Fin N => 0 ≤ s x)

def inserteds (s : ℕ → ℤ) (N : ℕ) : Finset (Fin N) :=
  univ.filter (fun x : Fin N => s x < 0)

theorem original_card (s : ℕ → ℤ) (N : ℕ) :
    (originals s N).card + (inserteds s N).card = N := by
  have := filter_card_add_filter_neg_card_eq_card
    (s := (univ : Finset (Fin N))) (fun x : Fin N => 0 ≤ s x)
  simpa [originals, inserteds, not_le] using this

/-- The retained originals use at most (A+D+1)r value-residue pairs. -/
theorem Edited.original_boxes {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) :
    ((originals s N).image (box s w v)).card ≤ (A + D + 1) * r := by
  let choices : Finset (ZMod r × ℤ) := univ ×ˢ Icc (-(D : ℤ)) (A : ℤ)
  let place : ZMod r × ℤ → W × ZMod r := fun p => (w p.1, p.1 + (p.2 : ZMod r))
  have hs : (originals s N).image (box s w v) ⊆ choices.image place := by
    intro b hb
    obtain ⟨x, hx, rfl⟩ := mem_image.mp hb
    have hx' : 0 ≤ s x := (mem_filter.mp hx).2
    refine mem_image.mpr ⟨((s x : ZMod r), (x.val : ℤ) - s x), ?_, ?_⟩
    · exact mem_product.mpr ⟨mem_univ _, mem_Icc.mpr (h.displacement x hx')⟩
    · simp [place, box, symbol, hx', Int.cast_sub]
  have hc : choices.card = r * (A + D + 1) := by
    simp only [choices, card_product, card_univ, ZMod.card, Int.card_Icc]
    congr 1
    omega
  calc
    _ ≤ (choices.image place).card := card_le_card hs
    _ ≤ choices.card := card_image_le
    _ = r * (A + D + 1) := hc
    _ = (A + D + 1) * r := Nat.mul_comm _ _

/-- Inserts can add at most one box apiece, including repeated symbols. -/
theorem Edited.all_boxes {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) :
    ((univ : Finset (Fin N)).image (box s w v)).card ≤
      (A + D + 1) * r + (inserteds s N).card := by
  have hs : (univ : Finset (Fin N)).image (box s w v) ⊆
      (originals s N).image (box s w v) ∪ (inserteds s N).image (box s w v) := by
    intro b hb
    obtain ⟨x, hx, rfl⟩ := mem_image.mp hb
    by_cases h0 : 0 ≤ s x
    · apply mem_union_left
      exact mem_image.mpr ⟨x, mem_filter.mpr ⟨mem_univ _, h0⟩, rfl⟩
    · apply mem_union_right
      exact mem_image.mpr ⟨x, mem_filter.mpr ⟨mem_univ _, by omega⟩, rfl⟩
  exact (card_le_card hs).trans ((card_union_le _ _).trans
    (Nat.add_le_add (h.original_boxes w v N) card_image_le))

end Glitchy
