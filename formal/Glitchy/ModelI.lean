import Glitchy.Literal

-- Model I: a simpler model in which the sequence consists of k+1 blocks. Block j has relative
-- length `length j` (the lengths add up to 1) and is shifted by `shift j` ∈ ZMod r.
-- In the report (Proposition "General glitches in normalized Model I") these are ℓ_j and s_j,
-- and `grouped shift length b` is p_b. The report numbers blocks 1, ..., k+1; Lean uses
-- `Fin (k + 1)` = {0, ..., k}.
-- Results:
--   * the ideal-outcome mass equals ∑_b (total length of the blocks with shift b)²
--     (`modelMass_eq`);
--   * it is ≥ 1/(k+1) (`modelI_bound`);
--   * it equals 1/(k+1) when all shifts differ and all lengths are equal (`modelI_sharp`);
--   * so no positive bound independent of k exists (`no_uniform_exact_bound`).
--
-- Lean ideas introduced here: `omit ... in`, `by_contra`, `Fin.ext`, `congrArg`, and proving
-- an `∃` statement by supplying a witness.

namespace Glitchy
open Finset Complex
-- `Fin (k + 1)` = {0, ..., k}, the block indices.
variable {r k : ℕ} [NeZero r]

/-- The total length of the blocks whose shift is `b`: the report's `p_b`. -/
noncomputable def grouped (shift : Fin (k + 1) → ZMod r)
    (length : Fin (k + 1) → ℝ) (b : ZMod r) : ℝ :=
  ∑ j ∈ univ.filter (fun j => shift j = b), length j

/-- Model I's probability on the exact ideal outcomes, `Pr_{Model I}[Y ∈ G]`:
`(1/r) · ∑_c |∑_j length j · e(c · shift j)|²`. -/
noncomputable def modelMass (shift : Fin (k + 1) → ZMod r)
    (length : Fin (k + 1) → ℝ) : ℝ :=
  (1 / (r : ℝ)) * ∑ c : ZMod r,
    normSq (∑ j, (length j : ℂ) * ZMod.stdAddChar (c * shift j))

-- Term-mode proof: this is exactly Mathlib's `sum_fiberwise`.
/-- The grouped lengths add up to the total length. -/
theorem grouped_sum (shift : Fin (k + 1) → ZMod r) (length : Fin (k + 1) → ℝ) :
    (∑ b, grouped shift length b) = ∑ j, length j :=
  sum_fiberwise univ shift length

-- `omit [NeZero r] in` removes the assumption r ≠ 0 from the `variable` line, for the next
-- theorem only. That theorem does not need it, so it now also covers r = 0. Without `omit`,
-- Lean's linter would warn about an unused assumption.
omit [NeZero r] in
/-- If no block has shift `b`, the grouped length at `b` is `0` (an empty sum). -/
theorem grouped_zero (shift : Fin (k + 1) → ZMod r) (length : Fin (k + 1) → ℝ)
    (b : ZMod r) (hb : b ∉ univ.image shift) : grouped shift length b = 0 := by
  -- The set of blocks with shift b is empty.
  have he : univ.filter (fun j => shift j = b) = ∅ := by
    -- `eq_empty_iff_forall_not_mem`: S = ∅ ↔ no element is in S.
    apply eq_empty_iff_forall_not_mem.mpr
    intro j hj
    -- If j were in it, b = shift j would be in the image of `shift`, contradicting hb.
    -- (`hb : b ∉ ...` is a function from proofs of `b ∈ ...` to `False`.)
    exact hb (mem_image.mpr ⟨j, mem_univ _, (mem_filter.mp hj).2⟩)
  -- Unfold `grouped`. A sum over ∅ is 0.
  simp [grouped, he]

/-- At most `k + 1` grouped lengths are nonzero, because there are only `k + 1` blocks. -/
theorem grouped_support (shift : Fin (k + 1) → ZMod r) (length : Fin (k + 1) → ℝ) :
    (univ.filter (fun b => grouped shift length b ≠ 0)).card ≤ k + 1 := by
  have hs : univ.filter (fun b => grouped shift length b ≠ 0) ⊆ univ.image shift := by
    intro b hb
    -- `by_contra hn`: assume the opposite, `hn : b ∉ univ.image shift`, and derive False.
    by_contra hn
    -- hb says grouped b ≠ 0, but `grouped_zero` says it is 0.
    exact (mem_filter.mp hb).2 (grouped_zero shift length b hn)
  -- |support| ≤ |image of shift| ≤ |Fin (k+1)| = k + 1.
  exact (card_le_card hs).trans (by simpa using (card_image_le (s := univ) (f := shift)))

/-- Grouping the Fourier amplitude by shift value:
`∑_b grouped(b) · e(c·b) = ∑_j length j · e(c · shift j)`. -/
theorem grouped_amplitude (shift : Fin (k + 1) → ZMod r) (length : Fin (k + 1) → ℝ)
    (c : ZMod r) :
    (∑ b, (grouped shift length b : ℂ) * ZMod.stdAddChar (c * b)) =
      ∑ j, (length j : ℂ) * ZMod.stdAddChar (c * shift j) := by
  unfold grouped
  push_cast
  -- Distribute: (∑_j length j) · e(c·b) = ∑_j length j · e(c·b).
  simp_rw [sum_mul]
  -- For the blocks j with shift j = b, replace e(c·b) by e(c · shift j).
  have he (b : ZMod r) :
      (∑ j ∈ univ.filter (fun j => shift j = b),
        (length j : ℂ) * ZMod.stdAddChar (c * b)) =
      ∑ j ∈ univ.filter (fun j => shift j = b),
        (length j : ℂ) * ZMod.stdAddChar (c * shift j) := by
    apply sum_congr rfl
    intro j hj
    rw [(mem_filter.mp hj).2]
  simp_rw [he]
  -- Summing over each shift value and then over all values = summing over all blocks.
  exact sum_fiberwise univ shift _

/-- Exact Model I total mass as the squared grouped lengths. -/
theorem modelMass_eq (shift : Fin (k + 1) → ZMod r) (length : Fin (k + 1) → ℝ) :
    modelMass shift length = ∑ b, (grouped shift length b)^2 := by
  unfold modelMass
  -- Rewrite the amplitude in grouped form (`grouped_amplitude`, right to left), then apply
  -- Parseval (Fourier.lean).
  simp_rw [← grouped_amplitude, residue_parseval]
  -- (1/r) · r = 1, which needs r ≠ 0.
  have hr : (r : ℝ) ≠ 0 := by exact_mod_cast (NeZero.ne r)
  field_simp

/-- Model I's lower bound: if the lengths add up to `1`, the mass is at least `1/(k+1)`. -/
theorem modelI_bound (shift : Fin (k + 1) → ZMod r) (length : Fin (k + 1) → ℝ)
    (hlength : ∑ j, length j = 1) :
    (1 : ℝ) / (k + 1) ≤ modelMass shift length := by
  rw [modelMass_eq]
  -- The grouped lengths add up to 1 and at most k+1 are nonzero, so Cauchy–Schwarz
  -- (`grouped_lengths` in Counting.lean) applies.
  exact grouped_lengths _ k ((grouped_sum shift length).trans hlength)
    (grouped_support shift length)

-- `Function.Injective shift`: all shifts are different.
-- `(fun _ => (1 : ℝ) / (k + 1))` is the constant function; `_` ignores its argument.
/-- Equal blocks with distinct shift residues attain the exact-set bound. -/
theorem modelI_sharp (shift : Fin (k + 1) → ZMod r) (hi : Function.Injective shift) :
    modelMass shift (fun _ => (1 : ℝ) / (k + 1)) = (1 : ℝ) / (k + 1) := by
  rw [modelMass_eq]
  -- At the box b = shift j only block j contributes, so the grouped length is 1/(k+1).
  -- (`hi.eq_iff`: shift i = shift j ↔ i = j. `filter_eq'`: filtering for "= j" leaves {j}.)
  have hp (j : Fin (k + 1)) :
      grouped shift (fun _ => (1 : ℝ) / (k + 1)) (shift j) = (1 : ℝ) / (k + 1) := by
    simp [grouped, hi.eq_iff, Finset.filter_eq']
  -- Only boxes in the image of `shift` contribute; the others are 0 ...
  have hs :
      (∑ b ∈ univ.image shift, (grouped shift (fun _ => (1 : ℝ) / (k + 1)) b)^2) =
      ∑ b, (grouped shift (fun _ => (1 : ℝ) / (k + 1)) b)^2 := by
    apply sum_subset (subset_univ _)
    intro b hb hn
    simp [grouped_zero shift _ b hn]
  -- ... and summing over the image of an injective map = summing over the blocks.
  rw [← hs, sum_image (fun a _ b _ hab => hi hab)]
  simp_rw [hp]
  -- k+1 copies of (1/(k+1))² add up to 1/(k+1).
  exact equal_grouped_lengths k

-- The example for sharpness: r = k+1 and block j is shifted by j.
/-- Equal-sign shifts 0,...,k are distinct modulo r=k+1. -/
def canonicalShift (k : ℕ) (j : Fin (k + 1)) : ZMod (k + 1) := (j.val : ZMod (k + 1))

/-- The shifts `0, ..., k` are different mod `k + 1`. -/
theorem canonicalShift_injective (k : ℕ) : Function.Injective (canonicalShift k) := by
  -- Suppose canonicalShift a = canonicalShift b; that is `hab`.
  intro a b hab
  -- `Fin.ext`: two elements of Fin (k+1) are equal if their values are equal.
  apply Fin.ext
  -- Apply `ZMod.val` to both sides of hab (`congrArg f : x = y → f x = f y`).
  have hv := congrArg ZMod.val hab
  -- val (j mod (k+1)) = j mod (k+1) = j, because j < k+1 (`a.isLt : a.val < k + 1`).
  simpa [canonicalShift, ZMod.val_natCast, Nat.mod_eq_of_lt a.isLt,
    Nat.mod_eq_of_lt b.isLt] using hv

-- To prove `∃ k, P k`, name a specific k and prove P k for it.
/-- Consequently there is no positive k-independent bound on the exact set G in Model I:
for every `ε > 0` some `k` makes the mass smaller than `ε`. -/
theorem no_uniform_exact_bound (eps : ℝ) (heps : 0 < eps) :
    ∃ k : ℕ, modelMass (canonicalShift k) (fun _ => (1 : ℝ) / (k + 1)) < eps := by
  -- The Archimedean property gives a natural number k > 1/ε.
  obtain ⟨k, hk⟩ := exists_nat_gt (1 / eps)
  -- `refine ⟨k, ?_⟩`: use this k as the witness. The goal becomes the statement for k.
  refine ⟨k, ?_⟩
  -- The mass is 1/(k+1) (`modelI_sharp`), so it remains to show 1/(k+1) < ε.
  rw [modelI_sharp _ (canonicalShift_injective k)]
  have hkp : (0 : ℝ) < (k : ℝ) + 1 := by positivity
  -- From 1/ε < k, get 1 < k·ε. (`div_lt_iff₀ heps : a / ε < b ↔ a < b · ε`.)
  have hk' : (1 : ℝ) < (k : ℝ) * eps := (div_lt_iff₀ heps).mp hk
  -- 1/(k+1) < ε ↔ 1 < ε·(k+1), which follows from hk' and ε > 0.
  apply (div_lt_iff₀ hkp).mpr
  nlinarith

end Glitchy
