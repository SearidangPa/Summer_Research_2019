import Mathlib

-- Pure inequalities about finite sums of real numbers. Nothing here knows about sequences,
-- edits, or Fourier transforms; later files plug their numbers into these results.
--
-- Lean ideas introduced here:
--   * `{B : Type*} [Fintype B]`: "for any finite type B" (later, B is the set of boxes)
--   * casts between ℕ and ℝ: `exact_mod_cast`, `push_cast`
--   * inequality tactics: `nlinarith`, `positivity`, and `calc` chains
--   * algebra tactics: `field_simp`, `ring`

namespace Glitchy
open Finset

-- `{B : Type*}`: B is any type. (`*` lets the type live in any universe; ignore that.)
-- `[Fintype B]`: B is finite, so `∑ b, ...` (a sum over all of B) makes sense.
-- `(n : B → ℝ)`: one real number per b. Think of n b as the count in box b.
-- `hq`: the number of b with n b ≠ 0 is at most Q.
-- `(Q : ℝ)`: Q cast from ℕ to ℝ, so that it can multiply a real number.
/-- Cauchy--Schwarz with an upper bound on the number of nonzero boxes:
if at most `Q` of the numbers `n b` are nonzero, then `(∑ n b)² ≤ Q · ∑ (n b)²`. -/
theorem support_energy {B : Type*} [Fintype B] (n : B → ℝ) (Q : ℕ)
    (hq : (univ.filter (fun b => n b ≠ 0)).card ≤ Q) :
    (∑ b, n b)^2 ≤ (Q : ℝ) * ∑ b, (n b)^2 := by
  -- `classical` lets the proof treat every yes/no question, such as whether a real number
  -- is 0, as decidable. `filter` needs that. This is ordinary classical logic.
  classical
  -- S is the support: the b with n b ≠ 0.
  let S := univ.filter (fun b => n b ≠ 0)
  -- Summing over S is the same as summing over all b, since n is 0 outside S.
  have hsum : ∑ b ∈ S, n b = ∑ b, n b := by
    -- `sum_subset h₁ h₂`: if S ⊆ T and f is 0 on T outside S, then ∑_S f = ∑_T f.
    -- `apply` uses it with h₁ supplied; the "0 outside S" condition becomes the new goal.
    apply sum_subset (filter_subset _ _)
    -- Given b with b ∉ S, show n b = 0. Unfolding S in hn gives ¬(n b ≠ 0), which simp
    -- turns into n b = 0.
    intro b hb hn
    simpa [S] using hn
  -- The same for the sum of squares.
  have hsq : ∑ b ∈ S, (n b)^2 = ∑ b, (n b)^2 := by
    apply sum_subset (filter_subset _ _)
    intro b hb hn
    have hz : n b = 0 := by simpa [S] using hn
    -- With n b = 0, the goal (n b)^2 = 0 is a simp computation.
    simp [hz]
  -- Mathlib's Cauchy–Schwarz: (∑_{b∈S} n b)² ≤ |S| · ∑_{b∈S} (n b)².
  -- `(s := S) (f := n)` passes arguments by name.
  have h := sq_sum_le_card_mul_sum_sq (s := S) (f := n)
  -- Turn the sums over S into sums over all b, inside h.
  rw [hsum, hsq] at h
  -- Finish with (∑ n)² ≤ |S| · ∑ n² ≤ Q · ∑ n², using |S| ≤ Q and ∑ n² ≥ 0.
  --   `mul_le_mul_of_nonneg_right : a ≤ b → 0 ≤ c → a * c ≤ b * c`
  --   `by exact_mod_cast hq` restates `hq : S.card ≤ Q` (in ℕ) as the same fact in ℝ.
  --   `sum_nonneg`: a sum of nonnegative terms is nonnegative. `sq_nonneg x : 0 ≤ x²`.
  --   `fun b _ => ...` takes b and an unused membership proof `_`.
  exact h.trans (mul_le_mul_of_nonneg_right (by exact_mod_cast hq)
    (sum_nonneg fun b _ => sq_nonneg (n b)))

-- The right-hand side (r / N²) · ∑ (n b)² is exactly the total probability on the ideal
-- outcomes (see `residue_mass` in Fourier.lean). With Q = (A+D+1)r + I occupied boxes this
-- gives the report's bound r / ((A+D+1)r + I).
/-- The finite-N box bound: no Fourier hypotheses are hidden here.
If the counts `n b` add up to `N > 0` and occupy at most `Q > 0` boxes, then
`r / Q ≤ (r / N²) · ∑ (n b)²`. -/
theorem mass_from_boxes {B : Type*} [Fintype B] (n : B → ℝ)
    (N r Q : ℕ) (hN : 0 < N) (hQ : 0 < Q)
    (hsum : ∑ b, n b = N)
    (hq : (univ.filter (fun b => n b ≠ 0)).card ≤ Q) :
    (r : ℝ) / Q ≤ (r : ℝ) / (N : ℝ)^2 * ∑ b, (n b)^2 := by
  -- Cauchy–Schwarz from above, then substitute ∑ n b = N, giving h : N² ≤ Q · ∑ (n b)².
  have h := support_energy n Q hq
  rw [hsum] at h
  -- N > 0 and Q > 0 as real numbers, from the ℕ facts hN and hQ.
  have hn : (0 : ℝ) < N := by exact_mod_cast hN
  have hq' : (0 : ℝ) < Q := by exact_mod_cast hQ
  -- Divide by Q: N²/Q ≤ ∑ (n b)².
  -- `div_le_iff₀ hq' : a / Q ≤ b ↔ a ≤ b * Q` (valid because Q > 0). `.mpr` goes ←.
  -- `nlinarith [h]` proves N² ≤ (∑ n²) · Q from h.
  have he : (N : ℝ)^2 / Q ≤ ∑ b, (n b)^2 :=
    (div_le_iff₀ hq').mpr (by nlinarith [h])
  -- `calc` proves a chain of steps. On each later line, the leading `_` is the previous
  -- line's right side. A `_` on the far right is filled in from the goal.
  calc
    -- Algebra: r/Q = (r/N²)·(N²/Q). `field_simp` clears the denominators (it uses hn and hq'
    -- to see they are nonzero), and here that alone closes the goal.
    (r : ℝ) / Q = (r : ℝ) / (N : ℝ)^2 * ((N : ℝ)^2 / Q) := by
      field_simp
    -- Multiply `he` by r/N² ≥ 0 (proved by `positivity`).
    -- `mul_le_mul_of_nonneg_left : b ≤ c → 0 ≤ a → a * b ≤ a * c`.
    _ ≤ _ := mul_le_mul_of_nonneg_left he (by positivity)

-- Here m b ≤ n b counts only the original (non-inserted) entries in box b. There are N - I of
-- them, and they occupy at most (k+1)·r boxes.
-- `(hm : ∀ b, 0 ≤ m b ∧ m b ≤ n b)` is a hypothesis that is itself a "for all";
-- `hm b` specializes it to one b.
/-- The retained-entry bound; zero retained entries are allowed.
If `0 ≤ m ≤ n` boxwise, the `m b` add up to `N - I`, and at most `(k+1)r` of them are
nonzero, then `(N - I)² / ((k+1) N²) ≤ (r / N²) · ∑ (n b)²`. -/
theorem mass_from_retained {B : Type*} [Fintype B] (n m : B → ℝ)
    (N r k I : ℕ) (hN : 0 < N) (hr : 0 < r)
    (hm : ∀ b, 0 ≤ m b ∧ m b ≤ n b)
    (hsum : ∑ b, m b = (N : ℝ) - I)
    (hq : (univ.filter (fun b => m b ≠ 0)).card ≤ (k + 1) * r) :
    ((N : ℝ) - I)^2 / (((k : ℝ) + 1) * (N : ℝ)^2)
      ≤ (r : ℝ) / (N : ℝ)^2 * ∑ b, (n b)^2 := by
  -- Cauchy–Schwarz for m: (N - I)² ≤ (k+1)r · ∑ (m b)².
  have h := support_energy m ((k + 1) * r) hq
  rw [hsum] at h
  -- Rewrite the cast ((k + 1) * r : ℕ) as ((k : ℝ) + 1) * (r : ℝ).
  push_cast at h
  have hn : (0 : ℝ) < N := by exact_mod_cast hN
  have hr' : (0 : ℝ) < r := by exact_mod_cast hr
  have hkr : (0 : ℝ) < ((k : ℝ) + 1) * r := by positivity
  -- Divide by (k+1)r.
  have he : ((N : ℝ) - I)^2 / (((k : ℝ) + 1) * r) ≤ ∑ b, (m b)^2 :=
    (div_le_iff₀ hkr).mpr (by nlinarith [h])
  -- ∑ m² ≤ ∑ n², because 0 ≤ m b ≤ n b in each box.
  have hs : ∑ b, (m b)^2 ≤ ∑ b, (n b)^2 := by
    -- `sum_le_sum`: compare the two sums term by term.
    apply sum_le_sum
    intro b hb
    have := hm b
    nlinarith
  calc
    -- The first `_` is the left side of the goal. Rewrite it with algebra.
    _ = (r : ℝ) / (N : ℝ)^2 *
        (((N : ℝ) - I)^2 / (((k : ℝ) + 1) * r)) := by field_simp; ring
    -- Chain he and hs, then multiply by r/N² ≥ 0.
    _ ≤ _ := mul_le_mul_of_nonneg_left (he.trans hs) (by positivity)

-- Why it holds: (k+1)r + I ≤ (k+1)r + kr = (2k+1)r.
/-- The bound that depends only on `k`: if `0 < r` and `I ≤ k`, then
`1 / (2k + 1) ≤ r / ((k+1) r + I)`. -/
theorem universal_bound (r k I : ℕ) (hr : 0 < r) (hI : I ≤ k) :
    (1 : ℝ) / (2 * k + 1) ≤ (r : ℝ) / (((k : ℝ) + 1) * r + I) := by
  -- Move the ℕ hypotheses into ℝ. (For naturals, 0 < r is the same as 1 ≤ r.)
  have hr' : (1 : ℝ) ≤ r := by exact_mod_cast hr
  have hI' : (I : ℝ) ≤ k := by exact_mod_cast hI
  have hk : (0 : ℝ) ≤ k := by positivity
  -- Cross-multiply. `div_le_div_iff₀ : a / b ≤ c / d ↔ a * d ≤ c * b` when b, d > 0.
  -- The two `by positivity` prove that the denominators are positive.
  apply (div_le_div_iff₀ (by positivity) (by positivity)).mpr
  -- Left to prove: 1 · ((k+1)r + I) ≤ r · (2k+1), i.e. I ≤ k·r. That is nonlinear
  -- because of the product k·r, hence `nlinarith`.
  nlinarith

/-- With no surviving insertions (`I = 0`) the bound simplifies:
`r / ((k+1) r + 0) = 1 / (k+1)`. -/
theorem no_insertions (r k : ℕ) (hr : 0 < r) :
    (r : ℝ) / (((k : ℝ) + 1) * r + 0) = 1 / ((k : ℝ) + 1) := by
  -- field_simp needs to know r ≠ 0, so put r > 0 (in ℝ) into context.
  have hr' : (0 : ℝ) < r := by exact_mod_cast hr
  field_simp; ring

-- When k/N → 0, the left side tends to 1/(k+1).
/-- The few-edits form: if `I ≤ k ≤ N` then
`(1 - k/N)² / (k+1) ≤ (N - I)² / ((k+1) N²)`. -/
theorem few_edits (N k I : ℕ) (hN : 0 < N) (hIk : I ≤ k) (hkN : k ≤ N) :
    (1 - (k : ℝ) / N)^2 / ((k : ℝ) + 1) ≤
      ((N : ℝ) - I)^2 / (((k : ℝ) + 1) * (N : ℝ)^2) := by
  have hn : (0 : ℝ) < N := by exact_mod_cast hN
  have hi : (I : ℝ) ≤ k := by exact_mod_cast hIk
  have hk : (k : ℝ) ≤ N := by exact_mod_cast hkN
  -- 0 ≤ N - k ≤ N - I, so (N - k)² ≤ (N - I)².
  have hs : ((N : ℝ) - k)^2 ≤ ((N : ℝ) - I)^2 := by nlinarith
  calc
    -- Rewrite the left side: (1 - k/N)²/(k+1) = (N - k)²/((k+1)N²). `field_simp` clears
    -- denominators, `ring_nf` puts both sides in a normal form, and `simp` finishes.
    _ = ((N : ℝ) - k)^2 / (((k : ℝ) + 1) * (N : ℝ)^2) := by
      field_simp
      ring_nf
      simp
    -- Same denominator and a larger numerator.
    -- `div_le_div_of_nonneg_right : a ≤ b → 0 ≤ c → a / c ≤ b / c`.
    _ ≤ _ := div_le_div_of_nonneg_right hs (by positivity)

/-- Model I's grouped-length bound: if the `p b` add up to 1 and at most `k + 1` are
nonzero, then `1/(k+1) ≤ ∑ (p b)²`. -/
theorem grouped_lengths {B : Type*} [Fintype B] (p : B → ℝ) (k : ℕ)
    (hsum : ∑ b, p b = 1)
    (hq : (univ.filter (fun b => p b ≠ 0)).card ≤ k + 1) :
    (1 : ℝ) / (k + 1) ≤ ∑ b, (p b)^2 := by
  -- Cauchy–Schwarz: 1 = (∑ p)² ≤ (k+1) · ∑ p².
  have h := support_energy p (k + 1) hq
  rw [hsum] at h
  push_cast at h
  -- `(by positivity : (0 : ℝ) < k + 1)` states exactly what `positivity` should prove.
  apply (div_le_iff₀ (by positivity : (0 : ℝ) < k + 1)).mpr
  nlinarith

-- `∑ _ : Fin (k + 1), ...` sums over k+1 indices. The `_` means the index is not used, so
-- every term is the same.
/-- Equal occupied boxes attain the grouped-length bound:
`k+1` copies of `(1/(k+1))²` add up to `1/(k+1)`. -/
theorem equal_grouped_lengths (k : ℕ) :
    (∑ _ : Fin (k + 1), ((1 : ℝ) / (k + 1))^2) = (1 : ℝ) / (k + 1) := by
  -- `simp` turns the sum of a constant into (number of terms) × (constant).
  simp
  field_simp
  ring

-- This is a term-mode proof: no `by`. The proof is one expression.
-- `sum_le_sum_of_subset_of_nonneg : G ⊆ H → (∀ b ∈ H, b ∉ G → 0 ≤ p b) → ∑_G p ≤ ∑_H p`.
-- Its second argument is a function of b and two unused proofs (`_ _`).
/-- Passing from G to any superset uses only nonnegative probabilities. -/
theorem mass_mono {B : Type*} (p : B → ℝ) (hp : ∀ b, 0 ≤ p b)
    {G H : Finset B} (hGH : G ⊆ H) :
    ∑ b ∈ G, p b ≤ ∑ b ∈ H, p b :=
  sum_le_sum_of_subset_of_nonneg hGH (fun b _ _ => hp b)

-- Since I ≤ k, replacing I by k makes the denominator larger and the fraction smaller.
/-- The stronger denominator bound obtained from I <= k before dropping r. -/
theorem insertion_budget_bound (r k I : ℕ) (hr : 0 < r) (hI : I ≤ k) :
    (r : ℝ) / (((k : ℝ) + 1) * r + k) ≤
      (r : ℝ) / (((k : ℝ) + 1) * r + I) := by
  have hr' : (0 : ℝ) < r := by exact_mod_cast hr
  have hi : (I : ℝ) ≤ k := by exact_mod_cast hI
  apply (div_le_div_iff₀ (by positivity) (by positivity)).mpr
  nlinarith

/-- Equivalent forms of the retained-entry bound:
divide the numerator and denominator by `N²`. -/
theorem retained_normalized (N k I : ℕ) (hN : 0 < N) :
    ((N : ℝ) - I)^2 / (((k : ℝ) + 1) * (N : ℝ)^2) =
      (1 - (I : ℝ) / N)^2 / ((k : ℝ) + 1) := by
  have hn : (0 : ℝ) < N := by exact_mod_cast hN
  field_simp
  ring_nf
  simp

-- `Filter.Tendsto f l (nhds a)` is Mathlib's general way to write "f → a". The filter `l`
-- says along what: N → ∞, x → 0, and so on. `nhds a` means "neighbourhoods of a".
-- Proof: the constant 1 tends to 1 (`tendsto_const_nhds`); subtracting δ → 0 gives
-- 1 - δ → 1 - 0 (`.sub h`); squaring gives → (1 - 0)² (`.pow 2`). `simpa` simplifies
-- (1 - 0)² to 1, which matches the goal.
/-- The asymptotic factor used when k/N tends to zero: if `δ → 0` then `(1 - δ)² → 1`. -/
theorem retained_factor_tendsto {X : Type*} {l : Filter X} {delta : X → ℝ}
    (h : Filter.Tendsto delta l (nhds 0)) :
    Filter.Tendsto (fun x => (1 - delta x)^2) l (nhds 1) := by
  simpa using ((tendsto_const_nhds : Filter.Tendsto (fun _ : X => (1 : ℝ)) l (nhds 1)).sub h).pow 2

end Glitchy
