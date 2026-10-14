import Mathlib

/-!
Finite expectation estimates and the two analytic regimes of the Fejér argument.
These lemmas are independent of the literal-edit bookkeeping.
-/
namespace Glitchy
open Finset

/-- A bounded test function separates a set from its complement. No disjoint
window decomposition is needed: the good outcomes are a finite set. -/
theorem mass_from_test_function {Y : Type*} [Fintype Y] [DecidableEq Y]
    (P K : Y → ℝ) (G : Finset Y) (q L : ℝ)
    (hP : ∀ y, 0 ≤ P y) (hPsum : ∑ y, P y = 1)
    (hK : ∀ y, K y ≤ 1) (hout : ∀ y, y ∉ G → K y ≤ q)
    (hq : q < 1) (hL : L ≤ ∑ y, P y * K y) :
    (L - q) / (1 - q) ≤ ∑ y ∈ G, P y := by
  have hpoint (y : Y) : P y * K y ≤
      (if y ∈ G then P y else 0) + q * (P y - (if y ∈ G then P y else 0)) := by
    by_cases hy : y ∈ G
    · simpa [hy] using mul_le_mul_of_nonneg_left (hK y) (hP y)
    · simpa [hy, mul_comm] using mul_le_mul_of_nonneg_left (hout y hy) (hP y)
  have hsum := sum_le_sum (fun y (_ : y ∈ (univ : Finset Y)) => hpoint y)
  simp only [sum_add_distrib, ← mul_sum, sum_sub_distrib] at hsum
  simp [hPsum] at hsum
  apply (div_le_iff₀ (by linarith : 0 < 1 - q)).mpr
  nlinarith

/-- The lower bound obtained by summing the matching moments. -/
noncomputable def fejerLower (t T : ℝ) : ℝ :=
  1 - 4 * (T^2 - 1) / (9 * t * T) - (T - 1) / (3 * t * T)

/-- Normalization of the edit-pair estimate uses only r ≥ 2 and k ≥ 1. -/
theorem normalised_matching_bound (k r ell t I C : ℝ)
    (hk : 1 ≤ k) (hr : 2 ≤ r) (hell : 1 ≤ ell) (ht : 0 < t) (hI : I ≤ k)
    (hC : 1 - (2 * (k + 1) * ell + 2 * I / r) / ((2 * k + 1) * t) ≤ C) :
    1 - (4 * ell + 1) / (3 * t) ≤ C := by
  have hr0 : 0 < r := by linarith
  have hk0 : 0 ≤ k := by linarith
  have hIr : 2 * I / r ≤ k := (div_le_iff₀ hr0).mpr (by nlinarith)
  have hd : 0 < (2 * k + 1) * t := by positivity
  have hratio : (2 * (k + 1) * ell + 2 * I / r) / ((2 * k + 1) * t) ≤
      (4 * ell + 1) / (3 * t) := by
    apply (div_le_div_iff₀ hd (by positivity : 0 < 3 * t)).mpr
    have hprod := mul_nonneg (by linarith : 0 ≤ k - 1) (by linarith : 0 ≤ 2 * ell - 1)
    have hnum : 3 * (2 * (k + 1) * ell + 2 * I / r) ≤ (4 * ell + 1) * (2 * k + 1) := by
      nlinarith
    nlinarith [mul_nonneg (by linarith : 0 ≤ (4 * ell + 1) * (2 * k + 1) -
      3 * (2 * (k + 1) * ell + 2 * I / r)) ht.le]
  linarith

/-- Chord estimate for sine, directly from its concavity on [0, π]. -/
theorem sin_chord (b u : ℝ) (hb : 0 ≤ b) (hbpi : b ≤ Real.pi)
    (hu : 0 ≤ u) (hu1 : u ≤ 1) :
    u * Real.sin b ≤ Real.sin (u * b) := by
  simpa [smul_eq_mul] using
    strictConcaveOn_sin_Icc.concaveOn.2
      (show (0 : ℝ) ∈ Set.Icc 0 Real.pi from ⟨le_rfl, Real.pi_pos.le⟩)
      (show b ∈ Set.Icc 0 Real.pi from ⟨hb, hbpi⟩)
      (show 0 ≤ 1 - u by linarith) hu (by ring : 1 - u + u = 1)

/-- The elementary sine estimate in the short-kernel regime. -/
theorem short_sine_bound (t : ℝ) (ht : 2 ≤ t) (ht6 : t ≤ 6) :
    (4 : ℝ) / 3 ≤ t * Real.sin (Real.pi / t)^2 := by
  have ht0 : 0 < t := by linarith
  have hx0 : 0 ≤ Real.pi / t := by positivity
  have hx2 : Real.pi / t ≤ Real.pi / 2 :=
    div_le_div_of_nonneg_left Real.pi_pos.le (by norm_num) (by linarith)
  have hs0 := Real.sin_nonneg_of_mem_Icc (show Real.pi / t ∈ Set.Icc 0 Real.pi from
    ⟨hx0, hx2.trans (by linarith [Real.pi_pos])⟩)
  by_cases ht3 : t ≤ 3
  · have hsin : Real.sin (Real.pi / 3) ≤ Real.sin (Real.pi / t) :=
      Real.sin_le_sin_of_le_of_le_pi_div_two (by linarith [Real.pi_pos]) hx2
        (div_le_div_of_nonneg_left Real.pi_pos.le ht0 ht3)
    have hsbase : 0 ≤ Real.sin (Real.pi / 3) := by
      rw [Real.sin_pi_div_three]; positivity
    have hsq := sq_le_sq₀ hsbase hs0 |>.mpr hsin
    rw [Real.sq_sin_pi_div_three] at hsq
    nlinarith [mul_nonneg ht0.le (show 0 ≤ Real.sin (Real.pi / t)^2 - 3/4 by linarith)]
  by_cases ht4 : t ≤ 4
  · have hsin : Real.sin (Real.pi / 4) ≤ Real.sin (Real.pi / t) :=
      Real.sin_le_sin_of_le_of_le_pi_div_two (by linarith [Real.pi_pos]) hx2
        (div_le_div_of_nonneg_left Real.pi_pos.le ht0 ht4)
    rw [Real.sin_pi_div_four] at hsin
    have hsqrt := Real.sq_sqrt (show (0 : ℝ) ≤ 2 by norm_num)
    have hsqrtn := Real.sqrt_nonneg (2 : ℝ)
    have hsq : (1 : ℝ)/2 ≤ Real.sin (Real.pi / t)^2 := by nlinarith
    nlinarith [mul_nonneg ht0.le (show 0 ≤ Real.sin (Real.pi / t)^2 - 1/2 by linarith)]
  · have hc := sin_chord (Real.pi / 4) (4 / t) (by positivity)
        (by linarith [Real.pi_pos]) (by positivity) ((div_le_one ht0).mpr (by linarith))
    have he : 4 / t * (Real.pi / 4) = Real.pi / t := by ring
    rw [he, Real.sin_pi_div_four] at hc
    have hsqrt := Real.sq_sqrt (show (0 : ℝ) ≤ 2 by norm_num)
    have hsqrtn := Real.sqrt_nonneg (2 : ℝ)
    have hc' : 2 * Real.sqrt 2 ≤ t * Real.sin (Real.pi / t) := by
      have := (div_le_iff₀ ht0).mp (show 2 * Real.sqrt 2 / t ≤ Real.sin (Real.pi / t) by
        convert hc using 1 ; ring)
      nlinarith
    have hs : 8 ≤ t^2 * Real.sin (Real.pi / t)^2 := by
      nlinarith [sq_nonneg (t * Real.sin (Real.pi / t) - 2 * Real.sqrt 2)]
    have hnonneg : 0 ≤ t * Real.sin (Real.pi / t)^2 := by positivity
    nlinarith [mul_nonneg (show 0 ≤ 6 - t by linarith) hnonneg]

/-- T=2 gives the stronger bound 3/8. -/
theorem short_fejer_mass (t p : ℝ) (ht : 2 ≤ t) (ht6 : t ≤ 6)
    (hE : fejerLower t 2 ≤ p + (1 - p) * Real.cos (Real.pi / t)^2) :
    (3 : ℝ) / 8 ≤ p := by
  have ht0 : 0 < t := by linarith
  have hs := short_sine_bound t ht ht6
  have htrig := Real.sin_sq_add_cos_sq (Real.pi / t)
  have hlow : fejerLower t 2 = 1 - 5 / (6 * t) := by
    unfold fejerLower
    field_simp
    ring
  rw [hlow] at hE
  have hmul : (1 - p) * (t * Real.sin (Real.pi / t)^2) ≤ 5/6 := by
    have he : 1 - p - (1 - p) * Real.cos (Real.pi / t)^2 ≤ 5 / (6 * t) := by linarith
    have := (le_div_iff₀ (show 0 < 6 * t by positivity)).mp he
    rw [show Real.cos (Real.pi / t)^2 = 1 - Real.sin (Real.pi / t)^2 by linarith] at this
    nlinarith
  by_cases hp : (3 : ℝ) / 8 ≤ p
  · exact hp
  have hplt : p < 3/8 := lt_of_not_ge hp
  have hpn : 0 ≤ 1 - p := by linarith
  nlinarith [mul_le_mul_of_nonneg_left hs hpn]

/-- The ceiling choice T=ceil(3t/4) has a uniform expectation lower bound. -/
theorem long_fejer_lower (t T : ℝ) (ht : 6 ≤ t) (hT : 0 < T)
    (hTupper : T ≤ 3 * t / 4 + 1) : (29 : ℝ) / 54 ≤ fejerLower t T := by
  have ht0 : 0 < t := by linarith
  have he : fejerLower t T = 1 - 4*T/(9*t) - 1/(3*t) + 7/(9*t*T) := by
    unfold fejerLower
    field_simp
    ring
  rw [he]
  have hpos : 0 ≤ 7/(9*t*T) := by positivity
  have h1 : 4*T/(9*t) ≤ 1/3 + 4/(9*t) := by
    apply (div_le_iff₀ (show 0 < 9*t by positivity)).mpr
    have heq : (1/3 + 4/(9*t)) * (9*t) = 3*t+4 := by
      field_simp
      ring
    rw [heq]
    linarith
  have h2 : 7/(9*t) ≤ (7:ℝ)/54 := by
    apply (div_le_iff₀ (show 0 < 9*t by positivity)).mpr
    nlinarith
  have he2 : 4/(9*t) + 1/(3*t) = 7/(9*t) := by ring
  linarith

/-- The geometric-sum outside bound is at most 64/243 in the long regime. -/
theorem long_fejer_outside (t T : ℝ) (ht : 6 ≤ t) (hT : 3*t/4 ≤ T) :
    1 / (T^2 * Real.sin (Real.pi/t)^2) ≤ (64 : ℝ)/243 := by
  have ht0 : 0 < t := by linarith
  have hT0 : 0 < T := by linarith
  have hc := sin_chord (Real.pi/3) (3/t) (by positivity)
    (by linarith [Real.pi_pos]) (by positivity) ((div_le_one ht0).mpr (by linarith))
  have he : 3/t * (Real.pi/3) = Real.pi/t := by ring
  rw [he, Real.sin_pi_div_three] at hc
  have hsqrtn := Real.sqrt_nonneg (3 : ℝ)
  have hsqrt := Real.sq_sqrt (show (0 : ℝ) ≤ 3 by norm_num)
  have hs0 : 0 ≤ Real.sin (Real.pi/t) := le_trans (by positivity) hc
  have hct : 3 * Real.sqrt 3 / 2 ≤ t * Real.sin (Real.pi/t) := by
    have hc' : (3 * Real.sqrt 3 / 2) / t ≤ Real.sin (Real.pi/t) := by
      convert hc using 1 ; ring
    simpa [mul_comm] using (div_le_iff₀ ht0).mp hc'
  have hprod : 9 * Real.sqrt 3 / 8 ≤ T * Real.sin (Real.pi/t) := by
    nlinarith [mul_le_mul_of_nonneg_right hT hs0]
  have hden : (243 : ℝ)/64 ≤ T^2 * Real.sin (Real.pi/t)^2 := by
    nlinarith [sq_nonneg (T * Real.sin (Real.pi/t) - 9 * Real.sqrt 3 / 8)]
  apply (div_le_iff₀ (show 0 < T^2 * Real.sin (Real.pi/t)^2 by linarith)).mpr
  nlinarith

/-- Combine the long-regime expectation and outside bounds. -/
theorem long_fejer_mass (p E q : ℝ) (hE : (29 : ℝ)/54 ≤ E)
    (hq : q ≤ (64 : ℝ)/243) (hp : p ≤ 1) (hupper : E ≤ p + (1-p)*q) :
    (133 : ℝ)/358 ≤ p := by
  have := mul_le_mul_of_nonneg_left hq (show 0 ≤ 1-p by linarith)
  linarith

end Glitchy
