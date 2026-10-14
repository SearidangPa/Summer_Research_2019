import Glitchy.Edits
import Glitchy.Counting

-- Counting occupied boxes. A box is a pair (u, b) = (symbol, position mod r), and each of the
-- first N positions lands in exactly one box. Counting.lean needs an upper bound on the number
-- of occupied boxes; this file proves it is at most (A+D+1)·r + I, the report's Q in Eq.
-- (box-count) (Section 4.1, Step 2).
--
-- Lean ideas introduced here:
--   * `section` and `variable`: declare arguments once and reuse them in many theorems
--   * `noncomputable def`, for definitions that produce real numbers
--   * `ZMod r` (integers mod r) and the `[NeZero r]` assumption
--   * `ext`, `constructor`, `refine` with `?_` holes, `by_cases`

namespace Glitchy
open Finset

-- A `section` limits how far the `variable` line below applies.
section Fibers
-- `variable` declares arguments that Lean adds automatically to every theorem in the section
-- that mentions them, so `{X B : Type*} [DecidableEq B]` need not be repeated.
-- `[DecidableEq B]`: equality on B can be tested by an algorithm. `filter (fun x => f x = b)`
-- needs that.
variable {X B : Type*} [DecidableEq B]

-- `noncomputable`: Lean cannot run code that produces an exact real number, so it asks us to
-- mark such definitions. Mathematically this is an ordinary definition.
-- Later, `f` is the map "position ↦ its box", and `fiberCounts f T b` is the count in box b.
/-- Number of entries of T in one box, as a real number:
the size of `{x ∈ T | f x = b}` (the "fiber" of `f` over `b`). -/
noncomputable def fiberCounts (f : X → B) (T : Finset X) (b : B) : ℝ :=
  ((T.filter (fun x => f x = b)).card : ℝ)

/-- The counts over all boxes add up to `|T|`, because every element is in exactly one box. -/
theorem sum_fiberCounts [Fintype B] (f : X → B) (T : Finset X) :
    (∑ b, fiberCounts f T b) = T.card := by
  -- Replace `fiberCounts` by its definition.
  unfold fiberCounts
  -- The goal is an equation between casts ℕ → ℝ. `norm_cast` turns it into the same
  -- equation in ℕ.
  norm_cast
  -- Mathlib's `card_eq_sum_card_fiberwise` says |T| = ∑_b |{x ∈ T | f x = b}|, provided
  -- every f x lies in the set we sum over (here `univ`, so `mem_univ` proves it).
  -- `.symm` flips an equation a = b into b = a.
  exact (card_eq_sum_card_fiberwise (fun x _ => mem_univ (f x))).symm

/-- The occupied boxes (count ≠ 0) are exactly the image `f(T)`. -/
theorem support_fiberCounts [Fintype B] (f : X → B) (T : Finset X) :
    univ.filter (fun b => fiberCounts f T b ≠ 0) = T.image f := by
  -- `ext b`: two finite sets are equal if they have the same members b.
  ext b
  -- `simp only` rewrites with exactly these lemmas. The chain is: count ≠ 0 ↔ fiber is
  -- nonempty ↔ ∃ x ∈ T, f x = b ↔ b ∈ T.image f.
  simp only [mem_filter, mem_univ, true_and, fiberCounts, Nat.cast_ne_zero,
    card_ne_zero, filter_nonempty_iff, mem_image]

/-- Counts are nonnegative, and a subset has a smaller count in every box. -/
theorem fiberCounts_mono (f : X → B) {S T : Finset X} (hST : S ⊆ T) (b : B) :
    0 ≤ fiberCounts f S b ∧ fiberCounts f S b ≤ fiberCounts f T b := by
  -- `constructor` splits the goal `P ∧ Q` into a goal P and a goal Q.
  constructor
  -- A natural number cast into ℝ is ≥ 0; `positivity` knows that.
  · unfold fiberCounts; positivity
  -- `filter_subset_filter`: S ⊆ T gives S.filter p ⊆ T.filter p. `card_le_card`: a subset
  -- is no larger. `exact_mod_cast` turns the resulting ℕ inequality into one in ℝ.
  · unfold fiberCounts
    exact_mod_cast card_le_card (filter_subset_filter _ hST)

end Fibers

-- From here on, W is the alphabet of output symbols and r is the period.
-- `ZMod r` is the integers mod r: {0, 1, ..., r-1} with arithmetic mod r.
-- `[NeZero r]` assumes r ≠ 0. (Mathlib defines `ZMod 0` as ℤ itself, which is not finite.)
variable {W : Type*} [DecidableEq W] {r : ℕ} [NeZero r]

-- `w : ZMod r → W` is the periodic pattern (the report's w_0, ..., w_{r-1}) and `v : ℤ → W`
-- gives inserted symbols (the report's v).
-- `(t : ZMod r)` reduces the integer t mod r.
/-- Tags select original periodic symbols or arbitrary inserted symbols.
Distinctness of original symbols is not needed for this lower bound.
An original tag `t ≥ 0` shows `w (t mod r)`, so the unedited sequence has period `r`; an
inserted tag `t < 0` shows `v t`, which may repeat or match an original symbol.
The edited sequence is `f_k(x) = symbol w v (τ x)`, as in the report's Step 2. -/
def symbol (w : ZMod r → W) (v : ℤ → W) (t : ℤ) : W :=
  if 0 ≤ t then w (t : ZMod r) else v t

-- `{N : ℕ}` is implicit because Lean can read N off the type of `x : Fin N`.
-- `W × ZMod r` is the type of pairs, and `(a, b)` builds a pair.
/-- The box of position `x`: its symbol together with `x mod r`. -/
def box (τ : ℕ → ℤ) (w : ZMod r → W) (v : ℤ → W) {N : ℕ}
    (x : Fin N) : W × ZMod r :=
  (symbol w v (τ x), (x.val : ZMod r))

/-- The positions among the first `N` that hold an original tag (`≥ 0`). -/
def originals (τ : ℕ → ℤ) (N : ℕ) : Finset (Fin N) :=
  univ.filter (fun x : Fin N => 0 ≤ τ x)

/-- The positions among the first `N` that hold an inserted tag (`< 0`). Its size is `I`. -/
def inserteds (τ : ℕ → ℤ) (N : ℕ) : Finset (Fin N) :=
  univ.filter (fun x : Fin N => τ x < 0)

/-- Every position is original or inserted: `#originals + #inserteds = N`. -/
theorem original_card (τ : ℕ → ℤ) (N : ℕ) :
    (originals τ N).card + (inserteds τ N).card = N := by
  -- Mathlib: |{x | P x}| + |{x | ¬P x}| = |all|. `have :=` names the fact `this`.
  have := filter_card_add_filter_neg_card_eq_card
    (s := (univ : Finset (Fin N))) (fun x : Fin N => 0 ≤ τ x)
  -- `not_le` rewrites ¬(0 ≤ τ x) as τ x < 0, which matches `inserteds`.
  simpa [originals, inserteds, not_le] using this

-- Idea: an original entry with tag t at position x has box (w (t mod r), x mod r), and
-- x = t + s for a displacement s ∈ [-D, A] (the report's s). So the box is determined by the
-- pair (t mod r, s), and there are r · (A+D+1) such pairs.
/-- The retained originals use at most (A+D+1)r value-residue pairs. -/
theorem Edited.original_boxes {A D τ} (h : Edited A D τ)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) :
    ((originals τ N).image (box τ w v)).card ≤ (A + D + 1) * r := by
  -- All pairs (t mod r, s) with s ∈ [-D, A]. `×ˢ` is the product of two finite sets.
  let choices : Finset (ZMod r × ℤ) := univ ×ˢ Icc (-(D : ℤ)) (A : ℤ)
  -- The box that such a pair determines: (symbol, (t + s) mod r). `p.1` and `p.2` are the
  -- two components of the pair p.
  let place : ZMod r × ℤ → W × ZMod r := fun p => (w p.1, p.1 + (p.2 : ZMod r))
  -- Every box occupied by an original entry is `place` of some choice.
  have hs : (originals τ N).image (box τ w v) ⊆ choices.image place := by
    intro b hb
    -- b is the box of some original position x. `rfl` replaces b by `box τ w v x`.
    obtain ⟨x, hx, rfl⟩ := mem_image.mp hb
    -- x is original, i.e. its tag is ≥ 0.
    have hx' : 0 ≤ τ x := (mem_filter.mp hx).2
    -- `refine` is like `exact`, but each `?_` hole becomes a new goal. We choose the pair
    -- (τ x mod r, x - τ x) and must show (1) it is in `choices` and (2) `place` sends it to
    -- the box of x.
    refine mem_image.mpr ⟨((τ x : ZMod r), (x.val : ℤ) - τ x), ?_, ?_⟩
    -- (1) The displacement x - τ x is in [-D, A] by `Edited.displacement`.
    · exact mem_product.mpr ⟨mem_univ _, mem_Icc.mpr (h.displacement x hx')⟩
    -- (2) Unfold the definitions. The symbol is w (τ x) because τ x ≥ 0, and
    -- τ x + (x - τ x) = x mod r.
    · simp [place, box, symbol, hx', Int.cast_sub]
  -- Count the choices: r residues times (A + D + 1) displacements.
  have hc : choices.card = r * (A + D + 1) := by
    simp only [choices, card_product, card_univ, ZMod.card, Int.card_Icc]
    -- `congr 1` reduces `r * a = r * b` to `a = b`; then omega checks the arithmetic.
    congr 1
    omega
  calc
    _ ≤ (choices.image place).card := card_le_card hs  -- a subset is no larger
    _ ≤ choices.card := card_image_le                  -- an image is no larger than its source
    _ = r * (A + D + 1) := hc
    _ = (A + D + 1) * r := Nat.mul_comm _ _

-- Together with the previous theorem: all of the first N positions occupy at most
-- (A+D+1)·r + I boxes, where I = #inserteds.
/-- Inserts can add at most one box apiece, including repeated symbols. -/
theorem Edited.all_boxes {A D τ} (h : Edited A D τ)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) :
    ((univ : Finset (Fin N)).image (box τ w v)).card ≤
      (A + D + 1) * r + (inserteds τ N).card := by
  -- Every box is occupied by an original entry or by an inserted one.
  have hs : (univ : Finset (Fin N)).image (box τ w v) ⊆
      (originals τ N).image (box τ w v) ∪ (inserteds τ N).image (box τ w v) := by
    intro b hb
    obtain ⟨x, hx, rfl⟩ := mem_image.mp hb
    -- Two cases: x holds an original tag, or it does not.
    by_cases h0 : 0 ≤ τ x
    -- Original: the box is in the left part of the union.
    · apply mem_union_left
      exact mem_image.mpr ⟨x, mem_filter.mpr ⟨mem_univ _, h0⟩, rfl⟩
    -- Inserted: the box is in the right part. `by omega` turns ¬(0 ≤ τ x) into τ x < 0.
    · apply mem_union_right
      exact mem_image.mpr ⟨x, mem_filter.mpr ⟨mem_univ _, by omega⟩, rfl⟩
  -- |all boxes| ≤ |union| ≤ |original boxes| + |inserted boxes| ≤ (A+D+1)r + I.
  -- `card_union_le`: |S ∪ T| ≤ |S| + |T|. `Nat.add_le_add` adds two inequalities.
  exact (card_le_card hs).trans ((card_union_le _ _).trans
    (Nat.add_le_add (h.original_boxes w v N) card_image_le))

end Glitchy
