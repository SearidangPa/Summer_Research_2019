import Mathlib

/-! Exact Fourier averaging on the r ideal outcomes, with no asymptotic limit. -/

-- Discrete Fourier analysis on ZMod r (the integers mod r).
-- `ZMod.stdAddChar` is the function e : ZMod r → ℂ given by e(b) = exp(2πi·b/r).
-- It turns addition into multiplication, e(a + b) = e(a)·e(b), and e(-a) = conj(e(a)).
-- Main results: orthogonality, Parseval's identity, and the identity
--   ∑_c P(c) = (r / N²) · ∑ (count)²
-- that turns "probability on the ideal outcomes" into "sum of squared counts".
--
-- Lean ideas introduced here:
--   * complex numbers: `starRingEnd ℂ` is complex conjugation and `normSq z` is |z|²
--   * `simp_rw` (rewriting inside sums) and `←` (rewriting right to left)
--   * swapping sums with `sum_comm` and comparing sums term by term with `sum_congr`

namespace Glitchy
open Finset Complex
variable {r : ℕ} [NeZero r]

/-- Orthogonality of the standard characters, proved rather than assumed
(the report's identity in Section 4.1, Step 1):
`∑_c conj(e(c·b)) · e(c·b') = r` if `b = b'`, and `0` otherwise. -/
theorem character_orthogonality (b b' : ZMod r) :
    (∑ c : ZMod r, (starRingEnd ℂ) (ZMod.stdAddChar (c * b)) * ZMod.stdAddChar (c * b')) =
      if b = b' then (r : ℂ) else 0 := by
  -- `←` uses each lemma right to left: conj(e(x)) becomes e(-x), then e(x)·e(y) becomes
  -- e(x + y). `simp_rw` (unlike `rw`) can rewrite inside the `∑ c, ...`.
  simp_rw [← AddChar.map_neg_eq_conj, ← AddChar.map_add_eq_mul]
  -- A small local lemma. It takes an argument c and is proved by `ring` (algebra in a
  -- commutative ring).
  have he (c : ZMod r) : -(c * b) + c * b' = (b' - b) * c := by ring
  -- The sum is now ∑_c e((b' - b)·c).
  simp_rw [he]
  by_cases h : b = b'
  -- b = b': every term is e(0) = 1, and there are r terms (`ZMod.card`: |ZMod r| = r).
  · simp [h, ZMod.card]
  -- b ≠ b': `if_false` selects the `else` branch, so the sum must be shown to be 0.
  --   `ZMod.isPrimitive_stdAddChar`: for a ≠ 0, c ↦ e(a·c) is a nontrivial character.
  --   `AddChar.sum_eq_zero_of_ne_one`: a nontrivial character sums to 0 over the group.
  --   `sub_ne_zero.mpr (Ne.symm h)` proves b' - b ≠ 0 from b ≠ b'.
  · simp only [h, if_false]
    exact AddChar.sum_eq_zero_of_ne_one
      (ZMod.isPrimitive_stdAddChar r (sub_ne_zero.mpr (Ne.symm h)))

-- Proof idea: expand |z|² = conj(z)·z into a triple sum over c, b, b', then use orthogonality
-- in c. Every term with b ≠ b' vanishes.
/-- Parseval for one row of real-valued residue counts:
`∑_c |∑_b n b · e(c·b)|² = r · ∑_b (n b)²`. -/
theorem residue_parseval (n : ZMod r → ℝ) :
    (∑ c : ZMod r, normSq (∑ b : ZMod r, (n b : ℂ) * ZMod.stdAddChar (c * b))) =
      (r : ℝ) * ∑ b : ZMod r, (n b)^2 := by
  -- It is easier to work in ℂ. `Complex.ofReal_injective` says the map ℝ → ℂ is injective,
  -- so it is enough to prove the equation after casting both sides into ℂ.
  apply Complex.ofReal_injective
  -- Push the casts inward (the cast of a sum is the sum of the casts, and so on).
  push_cast
  -- |z|² = conj z · z; conj of a sum or product is the sum or product of the conjs; a real
  -- number is its own conj. Then multiply the two sums out into a double sum.
  simp_rw [normSq_eq_conj_mul_self, map_sum, map_mul, conj_ofReal,
    sum_mul, mul_sum]
  -- Regroup each term as (n b · n b') · (conj e(c·b) · e(c·b')).
  have rearrange (c b b' : ZMod r) :
      (n b : ℂ) * (starRingEnd ℂ) (ZMod.stdAddChar (c * b)) *
        ((n b' : ℂ) * ZMod.stdAddChar (c * b')) =
      (n b : ℂ) * n b' *
        ((starRingEnd ℂ) (ZMod.stdAddChar (c * b)) * ZMod.stdAddChar (c * b')) := by ring
  simp_rw [rearrange]
  -- Swap the order of summation (`sum_comm`: ∑_x ∑_y f = ∑_y ∑_x f), putting b outermost.
  rw [sum_comm]
  -- `sum_congr rfl`: two sums over the same set are equal if they agree term by term.
  -- So fix b.
  apply sum_congr rfl
  intro b hb
  -- Swap again so that the sum over c is innermost.
  rw [sum_comm]
  -- Pull the constant n b · n b' out of the sum over c (`← mul_sum`), then apply
  -- orthogonality to that sum.
  simp_rw [← mul_sum, character_orthogonality]
  -- Only b' = b survives, leaving n b · n b · r = r · (n b)². `simp` finishes.
  simp [mul_comm, pow_two]

-- Here `n u b` is a table: how many of the first N positions show symbol u and are ≡ b mod r.
-- `n : W → ZMod r → ℝ` is a function of two arguments. Literal.lean proves this formula
-- equals the actual QFT measurement probability.
/-- Probability at the ideal outcome indexed by c, expressed in residue counts.
This is the report's `P_f(cN/r)` with box counts `n_{ub}` (Section 4.1, Step 1):
`P(c) = (1/N²) · ∑_u |∑_b n u b · e(c·b)|²`. -/
noncomputable def idealProbability {W : Type*} [Fintype W]
    (N : ℕ) (n : W → ZMod r → ℝ) (c : ZMod r) : ℝ :=
  (1 / (N : ℝ)^2) *
    ∑ u, normSq (∑ b : ZMod r, (n u b : ℂ) * ZMod.stdAddChar (c * b))

/-- The exact sum-of-squared-counts identity, the report's Eq. (residue-mass) in Step 1:
`∑_c P(c) = (r / N²) · ∑_u ∑_b (n u b)²` (Parseval applied to each row `u`). -/
theorem residue_mass {W : Type*} [Fintype W] (N : ℕ) (n : W → ZMod r → ℝ) :
    (∑ c : ZMod r, idealProbability N n c) =
      (r : ℝ) / (N : ℝ)^2 * ∑ u, ∑ b, (n u b)^2 := by
  -- Unfold the definition and pull the constant 1/N² out of the sum over c.
  simp only [idealProbability, ← mul_sum]
  -- Swap ∑_c ∑_u into ∑_u ∑_c.
  rw [sum_comm]
  -- Apply Parseval to each row u, then pull the factor r out of the sum over u.
  simp_rw [residue_parseval, ← mul_sum]
  -- (1/N²) · (r · S) = (r / N²) · S.
  ring

/-- Probabilities are nonnegative: `1/N² ≥ 0` and each `|·|² ≥ 0`. -/
theorem idealProbability_nonneg {W : Type*} [Fintype W]
    (N : ℕ) (n : W → ZMod r → ℝ) (c : ZMod r) :
    0 ≤ idealProbability N n c := by
  unfold idealProbability
  -- `mul_nonneg : 0 ≤ a → 0 ≤ b → 0 ≤ a * b`. The first proof is given here
  -- (`positivity`), and `apply` leaves the second as the goal.
  apply mul_nonneg (by positivity)
  exact sum_nonneg (fun u _ => normSq_nonneg _)

end Glitchy
