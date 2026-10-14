import Glitchy.Matching
import Glitchy.Fejer

/-! The finite Fejér kernel, defined as a squared character sum. -/
namespace Glitchy
open Finset Complex
variable {N : ℕ} [NeZero N]

/-- Unnormalized squared geometric sum `|∑_{j<T} e(j·a)|²`. With the phase `θ` of `a`, this is
`T · F_T(θ)` for the report's Fejér kernel `F_T`, so the test function `fejerSquare T a / T²`
is the report's `F_T/T`. -/
noncomputable def fejerSquare (T : ℕ) (a : ZMod N) : ℝ :=
  normSq (∑ j ∈ range T, ZMod.stdAddChar ((j : ZMod N) * a))

/-- Triangle inequality bounds the normalized kernel by one (`F_T ≤ T` in the report). -/
theorem fejerSquare_le (T : ℕ) (a : ZMod N) : fejerSquare T a ≤ (T : ℝ)^2 := by
  unfold fejerSquare
  rw [normSq_eq_norm_sq]
  have hb := norm_sum_le (range T) (fun j : ℕ => ZMod.stdAddChar ((j : ZMod N)*a))
  simp only [AddChar.norm_apply, sum_const, card_range, nsmul_eq_mul, mul_one] at hb
  exact pow_le_pow_left₀ (norm_nonneg _) hb 2

/-- Real parts of the character are unchanged by negating their argument. -/
theorem char_re_neg (a : ZMod N) : (ZMod.stdAddChar (-a)).re = (ZMod.stdAddChar a).re := by
  rw [AddChar.map_neg_eq_conj]
  rfl

private theorem fejer_cross (T : ℕ) (a : ZMod N) :
    ((∑ j ∈ range T, ZMod.stdAddChar ((j : ZMod N)*a)) *
      (starRingEnd ℂ) (ZMod.stdAddChar ((T : ZMod N)*a))).re =
      ∑ j ∈ range T, (ZMod.stdAddChar (((j+1 : ℕ) : ZMod N)*a)).re := by
  simp only [sum_mul, re_sum]
  have he (j : ℕ) (hj : j ∈ range T) :
      (ZMod.stdAddChar ((j : ZMod N)*a) *
        (starRingEnd ℂ) (ZMod.stdAddChar ((T : ZMod N)*a))).re =
      (ZMod.stdAddChar (((T-1-j+1 : ℕ) : ZMod N)*a)).re := by
    have hjT := mem_range.mp hj
    rw [← AddChar.map_neg_eq_conj, ← AddChar.map_add_eq_mul]
    have hd : (((T-1-j+1 : ℕ) : ZMod N)*a) = -((j : ZMod N)*a + -((T : ZMod N)*a)) := by
      have hn : T-1-j+1+j = T := by omega
      have hc : ((T-1-j+1 : ℕ) : ZMod N) + j = T := by
        simpa only [Nat.cast_add] using congrArg (fun n : ℕ => (n : ZMod N)) hn
      linear_combination hc * a
    rw [hd, char_re_neg]
  rw [sum_congr rfl he]
  exact sum_range_reflect (fun j : ℕ => (ZMod.stdAddChar (((j+1 : ℕ) : ZMod N)*a)).re) T

/-- Exact triangular expansion of the squared geometric sum: the report's
`F_T(θ) = 1 + 2 ∑_{ℓ=1}^{T-1} (1 - ℓ/T) cos(2πℓθ)`, multiplied by `T` (here `ℓ = j + 1`). -/
theorem fejerSquare_expansion (T : ℕ) (a : ZMod N) :
    fejerSquare T a = (T : ℝ) + 2 * ∑ j ∈ range T,
      ((T : ℝ)-1-j) * (ZMod.stdAddChar (((j+1 : ℕ) : ZMod N)*a)).re := by
  induction T with
  | zero => simp [fejerSquare]
  | succ T ih =>
    have hr : fejerSquare (T+1) a = fejerSquare T a + 1 +
        2 * ∑ j ∈ range T, (ZMod.stdAddChar (((j+1 : ℕ) : ZMod N)*a)).re := by
      unfold fejerSquare
      rw [sum_range_succ, normSq_add, fejer_cross]
      simp [normSq_eq_norm_sq, AddChar.norm_apply]
    rw [hr, ih, sum_range_succ]
    push_cast
    have hs : (∑ j ∈ range T, ((T : ℝ)+1-1-j) *
        (ZMod.stdAddChar (((j+1 : ℕ) : ZMod N)*a)).re) =
        (∑ j ∈ range T, ((T : ℝ)-1-j) *
          (ZMod.stdAddChar (((j+1 : ℕ) : ZMod N)*a)).re) +
        ∑ j ∈ range T, (ZMod.stdAddChar (((j+1 : ℕ) : ZMod N)*a)).re := by
      rw [← sum_add_distrib]
      apply sum_congr rfl
      intros
      ring
    push_cast at hs
    rw [hs]
    ring

/-- The two polynomial sums used by the expectation estimate, proved by induction. -/
theorem fejer_weight_sum (T : ℕ) :
    (∑ j ∈ range T, ((T : ℝ)-1-j)) = (T : ℝ)*(T-1)/2 := by
  induction T with
  | zero => simp
  | succ T ih =>
    rw [sum_range_succ]
    push_cast
    have hs : (∑ j ∈ range T, ((T : ℝ)+1-1-j)) =
        (∑ j ∈ range T, ((T : ℝ)-1-j)) + T := by
      rw [← show (∑ _j ∈ range T, (1 : ℝ)) = T by simp, ← sum_add_distrib]
      apply sum_congr rfl
      intros
      ring
    rw [hs, ih]
    ring

theorem fejer_weighted_lag_sum (T : ℕ) :
    (∑ j ∈ range T, ((T : ℝ)-1-j) * ((j : ℝ)+1)) =
      (T : ℝ)*(T^2-1)/6 := by
  have htri (T : ℕ) : (∑ j ∈ range T, ((j : ℝ)+1)) = (T : ℝ)*(T+1)/2 := by
    induction T with
    | zero => simp
    | succ T ih => rw [sum_range_succ, ih]; push_cast; ring
  induction T with
  | zero => simp
  | succ T ih =>
    rw [sum_range_succ]
    push_cast
    have hs : (∑ j ∈ range T, ((T : ℝ)+1-1-j) * ((j : ℝ)+1)) =
        (∑ j ∈ range T, ((T : ℝ)-1-j) * ((j : ℝ)+1)) +
        ∑ j ∈ range T, ((j : ℝ)+1) := by
      rw [← sum_add_distrib]
      apply sum_congr rfl
      intros
      ring
    rw [hs, ih, htri]
    ring

/-- The geometric sum identity controls the squared kernel away from phase zero. -/
theorem fejerSquare_geometric (T : ℕ) (a : ZMod N) :
    fejerSquare T a * normSq (ZMod.stdAddChar a - 1) ≤ 4 := by
  have hp (j : ℕ) : ZMod.stdAddChar ((j : ZMod N)*a) = ZMod.stdAddChar a ^ j := by
    simpa only [nsmul_eq_mul] using (ZMod.stdAddChar (N := N)).map_nsmul_eq_pow j a
  have hid := congrArg normSq (geom_sum_mul (ZMod.stdAddChar a) T)
  rw [normSq_mul] at hid
  have hnorm : ‖ZMod.stdAddChar a ^ T - 1‖ ≤ 2 := by
    calc
      _ ≤ ‖ZMod.stdAddChar a ^ T‖ + ‖(1 : ℂ)‖ := norm_sub_le _ _
      _ = 2 := by norm_num [norm_pow, AddChar.norm_apply]
  have hb : normSq (ZMod.stdAddChar a ^ T - 1) ≤ 4 := by
    rw [normSq_eq_norm_sq]
    nlinarith [norm_nonneg (ZMod.stdAddChar a ^ T - 1)]
  simpa only [fejerSquare, hp, hid] using hb

/-- Distance of a cyclic phase to the nearest integer, in [0,1/2]. For `a = r·y` this is the
report's `θ_y = min_m |r y/N - m|` (Section 4.2, Step 2). -/
noncomputable def phaseDistance (a : ZMod N) : ℝ :=
  min ((a.val : ℝ)/N) (1-(a.val : ℝ)/N)

theorem phaseDistance_mem (a : ZMod N) : phaseDistance a ∈ Set.Icc (0 : ℝ) (1/2) := by
  have hN : (0 : ℝ) < N := by exact_mod_cast NeZero.pos N
  have hv : (a.val : ℝ) < N := by exact_mod_cast a.val_lt
  have hx : 0 ≤ (a.val : ℝ)/N := by positivity
  have hx1 : (a.val : ℝ)/N ≤ 1 := (div_le_one hN).mpr hv.le
  unfold phaseDistance
  constructor
  · exact le_min hx (by linarith)
  · have := min_le_left ((a.val : ℝ)/N) (1-(a.val : ℝ)/N)
    have := min_le_right ((a.val : ℝ)/N) (1-(a.val : ℝ)/N)
    linarith

theorem stdAddChar_re (a : ZMod N) :
    (ZMod.stdAddChar a).re = Real.cos (2 * Real.pi * ((a.val : ℝ)/N)) := by
  have he := ZMod.stdAddChar_coe (N := N) (a.val : ℤ)
  simp only [Int.cast_natCast, ZMod.natCast_zmod_val] at he
  have hp : 2 * (Real.pi : ℂ) * Complex.I * (a.val : ℂ) / N =
      ((2 * Real.pi * ((a.val : ℝ)/N) : ℝ) : ℂ) * Complex.I := by push_cast; ring
  rw [he, hp, Complex.exp_re]
  simp

/-- Folding the phase into [0,1/2] preserves the cosine. -/
theorem stdAddChar_re_distance (a : ZMod N) :
    (ZMod.stdAddChar a).re = Real.cos (2*Real.pi*phaseDistance a) := by
  rw [stdAddChar_re]
  unfold phaseDistance
  rcases le_total ((a.val : ℝ)/N) (1-(a.val : ℝ)/N) with h | h
  · rw [min_eq_left h]
  · rw [min_eq_right h]
    rw [show 2*Real.pi*(1-(a.val : ℝ)/N) = -(2*Real.pi*((a.val : ℝ)/N)) + 2*Real.pi by ring,
      Real.cos_add_two_pi, Real.cos_neg]

/-- Denominator in the geometric-sum formula. -/
theorem char_sub_one_normSq (a : ZMod N) :
    normSq (ZMod.stdAddChar a - 1) = 4 * Real.sin (Real.pi*phaseDistance a)^2 := by
  rw [normSq_sub]
  simp only [normSq_eq_norm_sq, AddChar.norm_apply, one_pow, norm_one,
    map_one, mul_one, stdAddChar_re_distance]
  rw [show 2*Real.pi*phaseDistance a = 2*(Real.pi*phaseDistance a) by ring, Real.cos_two_mul]
  nlinarith [Real.sin_sq_add_cos_sq (Real.pi*phaseDistance a)]

/-- The short Fejér kernel is exactly a squared cosine. -/
theorem fejerSquare_two (a : ZMod N) :
    fejerSquare 2 a / (2 : ℝ)^2 = Real.cos (Real.pi*phaseDistance a)^2 := by
  rw [fejerSquare_expansion]
  norm_num [sum_range_succ, stdAddChar_re_distance]
  rw [show 2*Real.pi*phaseDistance a = 2*(Real.pi*phaseDistance a) by ring, Real.cos_two_mul]
  ring

/-- Outside a phase window, the short kernel is bounded by cos²(π/t). -/
theorem fejer_short_outside (a : ZMod N) (t : ℝ) (ht : 2 ≤ t)
    (hd : 1/t ≤ phaseDistance a) :
    fejerSquare 2 a / (2 : ℝ)^2 ≤ Real.cos (Real.pi/t)^2 := by
  rw [fejerSquare_two]
  have hmem := phaseDistance_mem a
  have ht0 : 0 < t := by linarith
  have hx : 0 ≤ Real.pi/t := by positivity
  have hy : Real.pi*phaseDistance a ≤ Real.pi/2 := by
    nlinarith [mul_le_mul_of_nonneg_left hmem.2 Real.pi_pos.le]
  have hxy : Real.pi/t ≤ Real.pi*phaseDistance a := by
    simpa [div_eq_mul_inv] using mul_le_mul_of_nonneg_left hd Real.pi_pos.le
  have hc := Real.cos_le_cos_of_nonneg_of_le_pi hx (hy.trans (by linarith [Real.pi_pos])) hxy
  have hc0 := Real.cos_nonneg_of_mem_Icc (show Real.pi*phaseDistance a ∈ Set.Icc (-(Real.pi/2)) (Real.pi/2) from
    ⟨by nlinarith [Real.pi_pos], hy⟩)
  exact (sq_le_sq₀ hc0 (hc0.trans hc)).mpr hc

/-- Outside a phase window, the geometric-sum kernel has the claimed envelope. -/
theorem fejer_long_outside (T : ℕ) (a : ZMod N) (t : ℝ) (hT : 0 < T)
    (ht : 2 < t) (hd : 1/t ≤ phaseDistance a) :
    fejerSquare T a / (T : ℝ)^2 ≤ 1/((T : ℝ)^2 * Real.sin (Real.pi/t)^2) := by
  have ht0 : 0 < t := by linarith
  have hT0 : (0 : ℝ) < T := by exact_mod_cast hT
  have hmem := phaseDistance_mem a
  have hx : 0 < Real.pi/t := by positivity
  have hxpi : Real.pi/t < Real.pi := by
    apply (div_lt_iff₀ ht0).mpr
    nlinarith [Real.pi_pos]
  have hspos := Real.sin_pos_of_pos_of_lt_pi hx hxpi
  have hy : Real.pi*phaseDistance a ≤ Real.pi/2 := by
    nlinarith [mul_le_mul_of_nonneg_left hmem.2 Real.pi_pos.le]
  have hxy : Real.pi/t ≤ Real.pi*phaseDistance a := by
    simpa [div_eq_mul_inv] using mul_le_mul_of_nonneg_left hd Real.pi_pos.le
  have hs := Real.sin_le_sin_of_le_of_le_pi_div_two (by linarith [Real.pi_pos]) hy hxy
  have hsq := (sq_le_sq₀ hspos.le (hspos.le.trans hs)).mpr hs
  have hb := fejerSquare_geometric T a
  rw [char_sub_one_normSq] at hb
  have hnonneg : 0 ≤ fejerSquare T a := normSq_nonneg _
  have hmul := mul_le_mul_of_nonneg_left hsq hnonneg
  have hprod : fejerSquare T a * Real.sin (Real.pi/t)^2 ≤ 1 := by nlinarith
  apply (div_le_div_iff₀ (sq_pos_of_pos hT0) (mul_pos (sq_pos_of_pos hT0) (sq_pos_of_pos hspos))).mpr
  nlinarith [mul_le_mul_of_nonneg_left hprod (sq_nonneg (T : ℝ))]

section Expectation
variable {W : Type*} [Fintype W] [DecidableEq W]

/-- Exact expectation identity: every Fourier coefficient is an actual
matching fraction `C_{ℓr}`, by the already proved correlation identity.
The left side is `T · E_T` in the report's notation. -/
theorem fejerSquare_expectation (f : ZMod N → W) (r : ℕ) (T : ℕ) :
    (∑ y : ZMod N, cyclicProbability f y * fejerSquare T ((r : ZMod N)*y)) =
      (T : ℝ) + 2 * ∑ j ∈ range T,
        ((T : ℝ)-1-j) * matchingFraction f (((j+1 : ℕ)*r : ℕ) : ZMod N) := by
  simp_rw [fejerSquare_expansion, mul_add, mul_sum]
  rw [sum_add_distrib, ← sum_mul, cyclicProbability_sum, one_mul]
  congr 1
  have he (j : ℕ) (y : ZMod N) :
      cyclicProbability f y * (2 * (((T : ℝ)-1-j) *
        (ZMod.stdAddChar (((j+1 : ℕ) : ZMod N)*((r : ZMod N)*y))).re)) =
      2 * (((T : ℝ)-1-j) * (cyclicProbability f y *
        (ZMod.stdAddChar (y * (((j+1 : ℕ)*r : ℕ) : ZMod N))).re)) := by
    have harg : ((j+1 : ℕ) : ZMod N)*((r : ZMod N)*y) =
        y * (((j+1 : ℕ)*r : ℕ) : ZMod N) := by push_cast; ring
    rw [harg]
    ring
  simp_rw [he]
  rw [sum_comm]
  simp_rw [← mul_sum, cyclic_matching_moment_re]

/-- Matching lower bounds imply the report's Eq. (fejer-expectation):
`E_T/T ≥ 1 - 4(T²-1)/(9tT) - (T-1)/(3tT)`. -/
theorem fejer_expectation_lower (f : ZMod N → W) (r T : ℕ) (t : ℝ)
    (hT : 0 < T) (ht : 0 < t)
    (hmatch : ∀ ell : ℕ, 1 ≤ ell →
      1 - (4*(ell : ℝ)+1)/(3*t) ≤ matchingFraction f ((ell*r : ℕ) : ZMod N)) :
    fejerLower t T ≤
      (∑ y : ZMod N, cyclicProbability f y * fejerSquare T ((r : ZMod N)*y)) / (T : ℝ)^2 := by
  have hT0 : (0 : ℝ) < T := by exact_mod_cast hT
  have hb : (∑ j ∈ range T, ((T : ℝ)-1-j) * (1-(4*((j : ℝ)+1)+1)/(3*t))) ≤
      ∑ j ∈ range T, ((T : ℝ)-1-j) * matchingFraction f (((j+1 : ℕ)*r : ℕ) : ZMod N) := by
    apply sum_le_sum
    intro j hj
    have hjT : j+1 ≤ T := by have := mem_range.mp hj; omega
    have hw : 0 ≤ (T : ℝ)-1-j := by
      have hjTr : (j : ℝ)+1 ≤ T := by exact_mod_cast hjT
      linarith
    have hm := hmatch (j+1) (by omega)
    push_cast at hm
    simpa only [Nat.cast_add, Nat.cast_mul, Nat.cast_one] using mul_le_mul_of_nonneg_left hm hw
  have hs : (∑ j ∈ range T, ((T : ℝ)-1-j) * (1-(4*((j : ℝ)+1)+1)/(3*t))) =
      (T : ℝ)*(T-1)/2 - (4*((T : ℝ)*(T^2-1)/6) + (T : ℝ)*(T-1)/2)/(3*t) := by
    have he (j : ℕ) : ((T : ℝ)-1-j) * (1-(4*((j : ℝ)+1)+1)/(3*t)) =
        ((T : ℝ)-1-j) - (4*(((T : ℝ)-1-j)*((j : ℝ)+1)) + ((T : ℝ)-1-j))/(3*t) := by ring
    simp_rw [he]
    rw [sum_sub_distrib, ← sum_div, sum_add_distrib, ← mul_sum,
      fejer_weight_sum, fejer_weighted_lag_sum]
  rw [hs] at hb
  rw [fejerSquare_expectation]
  apply (le_div_iff₀ (sq_pos_of_pos hT0)).mpr
  have he : fejerLower t T * (T : ℝ)^2 =
      (T : ℝ) + 2 * ((T : ℝ)*(T-1)/2 - (4*((T : ℝ)*(T^2-1)/6) + (T : ℝ)*(T-1)/2)/(3*t)) := by
    unfold fejerLower
    field_simp
    ring
  rw [he]
  linarith
/-- Analytic concentration theorem: Steps 3 and 4 of the report's proof in Section 4.2.
Its two hypotheses are Eq. (normalised-matching) and `θ_y ≥ 1/t` outside `G` (Step 2).
The proof takes `T = 2` when `t ≤ 6` and `T = ⌈3t/4⌉` when `t > 6`.
Neither a Fejér identity nor an expectation estimate is assumed. -/
theorem fejer_mass_of_matching (f : ZMod N → W) (r : ℕ) (G : Finset (ZMod N))
    (t : ℝ) (ht : 2 ≤ t)
    (hmatch : ∀ ell : ℕ, 1 ≤ ell →
      1 - (4*(ell : ℝ)+1)/(3*t) ≤ matchingFraction f ((ell*r : ℕ) : ZMod N))
    (hsep : ∀ y, y ∉ G → 1/t ≤ phaseDistance ((r : ZMod N)*y)) :
    (133 : ℝ)/358 ≤ ∑ y ∈ G, cyclicProbability f y := by
  have ht0 : 0 < t := by linarith
  have kernel_le (T : ℕ) (hT : 0 < T) (y : ZMod N) :
      fejerSquare T ((r : ZMod N)*y) / (T : ℝ)^2 ≤ 1 := by
    have hT0 : (0 : ℝ) < T := by exact_mod_cast hT
    exact (div_le_one (sq_pos_of_pos hT0)).mpr (fejerSquare_le T _)
  have expectation_le (T : ℕ) (hT : 0 < T) :
      fejerLower t T ≤ ∑ y : ZMod N, cyclicProbability f y *
        (fejerSquare T ((r : ZMod N)*y) / (T : ℝ)^2) := by
    simpa only [← mul_div_assoc, ← sum_div] using fejer_expectation_lower f r T t hT ht0 hmatch
  by_cases ht6 : t ≤ 6
  · let q := Real.cos (Real.pi/t)^2
    have htrig := Real.sin_sq_add_cos_sq (Real.pi/t)
    have hs := short_sine_bound t ht ht6
    have hq : q < 1 := by dsimp [q]; nlinarith
    have hb := mass_from_test_function (cyclicProbability f)
      (fun y => fejerSquare 2 ((r : ZMod N)*y)/(2 : ℝ)^2) G q (fejerLower t 2)
      (cyclicProbability_nonneg f) (cyclicProbability_sum f)
      (kernel_le 2 (by omega))
      (fun y hy => fejer_short_outside _ t ht (hsep y hy)) hq (expectation_le 2 (by omega))
    have hm := (div_le_iff₀ (show 0 < 1-q by linarith)).mp hb
    have he : fejerLower t 2 ≤ (∑ y ∈ G, cyclicProbability f y) +
        (1-∑ y ∈ G, cyclicProbability f y)*Real.cos (Real.pi/t)^2 := by
      dsimp [q] at hm
      nlinarith
    exact (by norm_num : (133 : ℝ)/358 ≤ 3/8).trans
      (short_fejer_mass t _ ht ht6 he)
  · let T := ⌈3*t/4⌉₊
    have hTlower : 3*t/4 ≤ (T : ℝ) := Nat.le_ceil _
    have hT0 : (0 : ℝ) < T := by linarith
    have hT : 0 < T := by exact_mod_cast hT0
    have hTupper : (T : ℝ) ≤ 3*t/4+1 := by
      exact (Nat.ceil_lt_add_one (by positivity : (0 : ℝ) ≤ 3*t/4)).le
    have hlower := long_fejer_lower t T (by linarith) hT0 hTupper
    have hout (y : ZMod N) (hy : y ∉ G) :
        fejerSquare T ((r : ZMod N)*y)/(T : ℝ)^2 ≤ (64 : ℝ)/243 :=
      (fejer_long_outside T _ t hT (by linarith) (hsep y hy)).trans
        (long_fejer_outside t T (by linarith) hTlower)
    have hb := mass_from_test_function (cyclicProbability f)
      (fun y => fejerSquare T ((r : ZMod N)*y)/(T : ℝ)^2) G (64/243) (29/54)
      (cyclicProbability_nonneg f) (cyclicProbability_sum f)
      (kernel_le T hT) hout (by norm_num) (hlower.trans (expectation_le T hT))
    norm_num at hb ⊢
    exact hb
end Expectation

end Glitchy
