import Mathlib

namespace Glitchy
open Finset

/-- Cauchy--Schwarz with an upper bound on the number of nonzero boxes. -/
theorem support_energy {B : Type*} [Fintype B] (n : B → ℝ) (Q : ℕ)
    (hq : (univ.filter (fun b => n b ≠ 0)).card ≤ Q) :
    (∑ b, n b)^2 ≤ (Q : ℝ) * ∑ b, (n b)^2 := by
  classical
  let S := univ.filter (fun b => n b ≠ 0)
  have hsum : ∑ b ∈ S, n b = ∑ b, n b := by
    apply sum_subset (filter_subset _ _)
    intro b hb hn
    simpa [S] using hn
  have hsq : ∑ b ∈ S, (n b)^2 = ∑ b, (n b)^2 := by
    apply sum_subset (filter_subset _ _)
    intro b hb hn
    have hz : n b = 0 := by simpa [S] using hn
    simp [hz]
  have h := sq_sum_le_card_mul_sum_sq (s := S) (f := n)
  rw [hsum, hsq] at h
  exact h.trans (mul_le_mul_of_nonneg_right (by exact_mod_cast hq)
    (sum_nonneg fun b _ => sq_nonneg (n b)))

/-- The finite-N box bound: no Fourier hypotheses are hidden here. -/
theorem mass_from_boxes {B : Type*} [Fintype B] (n : B → ℝ)
    (N r Q : ℕ) (hN : 0 < N) (hQ : 0 < Q)
    (hsum : ∑ b, n b = N)
    (hq : (univ.filter (fun b => n b ≠ 0)).card ≤ Q) :
    (r : ℝ) / Q ≤ (r : ℝ) / (N : ℝ)^2 * ∑ b, (n b)^2 := by
  have h := support_energy n Q hq
  rw [hsum] at h
  have hn : (0 : ℝ) < N := by exact_mod_cast hN
  have hq' : (0 : ℝ) < Q := by exact_mod_cast hQ
  have he : (N : ℝ)^2 / Q ≤ ∑ b, (n b)^2 :=
    (div_le_iff₀ hq').mpr (by nlinarith [h])
  calc
    (r : ℝ) / Q = (r : ℝ) / (N : ℝ)^2 * ((N : ℝ)^2 / Q) := by
      field_simp
    _ ≤ _ := mul_le_mul_of_nonneg_left he (by positivity)

/-- The retained-entry bound; zero retained entries are allowed. -/
theorem mass_from_retained {B : Type*} [Fintype B] (n m : B → ℝ)
    (N r k I : ℕ) (hN : 0 < N) (hr : 0 < r)
    (hm : ∀ b, 0 ≤ m b ∧ m b ≤ n b)
    (hsum : ∑ b, m b = (N : ℝ) - I)
    (hq : (univ.filter (fun b => m b ≠ 0)).card ≤ (k + 1) * r) :
    ((N : ℝ) - I)^2 / (((k : ℝ) + 1) * (N : ℝ)^2)
      ≤ (r : ℝ) / (N : ℝ)^2 * ∑ b, (n b)^2 := by
  have h := support_energy m ((k + 1) * r) hq
  rw [hsum] at h
  push_cast at h
  have hn : (0 : ℝ) < N := by exact_mod_cast hN
  have hr' : (0 : ℝ) < r := by exact_mod_cast hr
  have hkr : (0 : ℝ) < ((k : ℝ) + 1) * r := by positivity
  have he : ((N : ℝ) - I)^2 / (((k : ℝ) + 1) * r) ≤ ∑ b, (m b)^2 :=
    (div_le_iff₀ hkr).mpr (by nlinarith [h])
  have hs : ∑ b, (m b)^2 ≤ ∑ b, (n b)^2 := by
    apply sum_le_sum
    intro b hb
    have := hm b
    nlinarith
  calc
    _ = (r : ℝ) / (N : ℝ)^2 *
        (((N : ℝ) - I)^2 / (((k : ℝ) + 1) * r)) := by field_simp; ring
    _ ≤ _ := mul_le_mul_of_nonneg_left (he.trans hs) (by positivity)

theorem universal_bound (r k I : ℕ) (hr : 0 < r) (hI : I ≤ k) :
    (1 : ℝ) / (2 * k + 1) ≤ (r : ℝ) / (((k : ℝ) + 1) * r + I) := by
  have hr' : (1 : ℝ) ≤ r := by exact_mod_cast hr
  have hI' : (I : ℝ) ≤ k := by exact_mod_cast hI
  have hk : (0 : ℝ) ≤ k := by positivity
  apply (div_le_div_iff₀ (by positivity) (by positivity)).mpr
  nlinarith

theorem no_insertions (r k : ℕ) (hr : 0 < r) :
    (r : ℝ) / (((k : ℝ) + 1) * r + 0) = 1 / ((k : ℝ) + 1) := by
  have hr' : (0 : ℝ) < r := by exact_mod_cast hr
  field_simp; ring

theorem few_edits (N k I : ℕ) (hN : 0 < N) (hIk : I ≤ k) (hkN : k ≤ N) :
    (1 - (k : ℝ) / N)^2 / ((k : ℝ) + 1) ≤
      ((N : ℝ) - I)^2 / (((k : ℝ) + 1) * (N : ℝ)^2) := by
  have hn : (0 : ℝ) < N := by exact_mod_cast hN
  have hi : (I : ℝ) ≤ k := by exact_mod_cast hIk
  have hk : (k : ℝ) ≤ N := by exact_mod_cast hkN
  have hs : ((N : ℝ) - k)^2 ≤ ((N : ℝ) - I)^2 := by nlinarith
  calc
    _ = ((N : ℝ) - k)^2 / (((k : ℝ) + 1) * (N : ℝ)^2) := by
      field_simp
      ring_nf
      simp
    _ ≤ _ := div_le_div_of_nonneg_right hs (by positivity)

/-- Model I's grouped-length bound. -/
theorem grouped_lengths {B : Type*} [Fintype B] (p : B → ℝ) (k : ℕ)
    (hsum : ∑ b, p b = 1)
    (hq : (univ.filter (fun b => p b ≠ 0)).card ≤ k + 1) :
    (1 : ℝ) / (k + 1) ≤ ∑ b, (p b)^2 := by
  have h := support_energy p (k + 1) hq
  rw [hsum] at h
  push_cast at h
  apply (div_le_iff₀ (by positivity : (0 : ℝ) < k + 1)).mpr
  nlinarith

/-- Equal occupied boxes attain the grouped-length bound. -/
theorem equal_grouped_lengths (k : ℕ) :
    (∑ _ : Fin (k + 1), ((1 : ℝ) / (k + 1))^2) = (1 : ℝ) / (k + 1) := by
  simp
  field_simp
  ring

/-- Passing from G to any superset uses only nonnegative probabilities. -/
theorem mass_mono {B : Type*} (p : B → ℝ) (hp : ∀ b, 0 ≤ p b)
    {G H : Finset B} (hGH : G ⊆ H) :
    ∑ b ∈ G, p b ≤ ∑ b ∈ H, p b :=
  sum_le_sum_of_subset_of_nonneg hGH (fun b _ _ => hp b)

/-- The stronger denominator bound obtained from I <= k before dropping r. -/
theorem insertion_budget_bound (r k I : ℕ) (hr : 0 < r) (hI : I ≤ k) :
    (r : ℝ) / (((k : ℝ) + 1) * r + k) ≤
      (r : ℝ) / (((k : ℝ) + 1) * r + I) := by
  have hr' : (0 : ℝ) < r := by exact_mod_cast hr
  have hi : (I : ℝ) ≤ k := by exact_mod_cast hI
  apply (div_le_div_iff₀ (by positivity) (by positivity)).mpr
  nlinarith

/-- Equivalent forms of the retained-entry bound. -/
theorem retained_normalized (N k I : ℕ) (hN : 0 < N) :
    ((N : ℝ) - I)^2 / (((k : ℝ) + 1) * (N : ℝ)^2) =
      (1 - (I : ℝ) / N)^2 / ((k : ℝ) + 1) := by
  have hn : (0 : ℝ) < N := by exact_mod_cast hN
  field_simp
  ring_nf
  simp

/-- The asymptotic factor used when k/N tends to zero. -/
theorem retained_factor_tendsto {X : Type*} {l : Filter X} {delta : X → ℝ}
    (h : Filter.Tendsto delta l (nhds 0)) :
    Filter.Tendsto (fun x => (1 - delta x)^2) l (nhds 1) := by
  simpa using ((tendsto_const_nhds : Filter.Tendsto (fun _ : X => (1 : ℝ)) l (nhds 1)).sub h).pow 2

end Glitchy
