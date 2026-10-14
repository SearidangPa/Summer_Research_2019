import Glitchy.PairCounting
import Glitchy.PeakGeometry

/-!
Section 4.2 end to end: arbitrary literal edit histories, arbitrary inserted
symbols, every nearest-integer tie choice, circular windows counted as a set,
and the existing complex-exponential Fourier measurement formula.
As in the report, `k = A + D` and `t = 2N/(r(2k+1))`.
-/
namespace Glitchy
open Finset
variable {W : Type*} [Fintype W] [DecidableEq W]
    {N r : ℕ} [NeZero N] [NeZero r]

omit [Fintype W] in
/-- The report's Eq. (normalised-matching), `C_{ℓr} ≥ 1 - (4ℓ+1)/(3t)`, follows from the
literal edit history. -/
theorem Edited.normalised_matching {A D τ} (hedit : Edited A D τ)
    (w : ZMod r → W) (v : ℤ → W) (hr : 2 ≤ r) (hk : 1 ≤ A+D)
    (ell : ℕ) (hell : 1 ≤ ell) :
    1 - (4*(ell : ℝ)+1)/(3*(2*(N : ℝ)/(r*(2*(A+D : ℕ)+1)))) ≤
      matchingFraction (literalCycle (N := N) τ w v) ((ell*r : ℕ) : ZMod N) := by
  have hN : (0 : ℝ) < N := by exact_mod_cast NeZero.pos N
  have hr0 : (0 : ℝ) < r := by exact_mod_cast NeZero.pos r
  have hr2 : (2 : ℝ) ≤ r := by exact_mod_cast hr
  have hk1 : (1 : ℝ) ≤ (A+D : ℕ) := by exact_mod_cast hk
  have hk0 : (0 : ℝ) < 2*(A+D : ℕ)+1 := by positivity
  have hb := hedit.matching_bound (N := N) w v (ell*r) (dvd_mul_left r ell)
  have he : (((A+D+1 : ℕ) : ℝ)*(ell*r : ℕ)+A)/N =
      (2*((A+D : ℕ)+1)*ell+2*(A : ℝ)/r)/
        ((2*(A+D : ℕ)+1)*(2*(N : ℝ)/(r*(2*(A+D : ℕ)+1)))) := by
    push_cast
    field_simp
    ring
  rw [he] at hb
  exact normalised_matching_bound (A+D : ℕ) r ell _ A _ hk1 hr2
    (by exact_mod_cast hell) (by positivity) (by exact_mod_cast Nat.le_add_right A D) hb

omit [NeZero r] in
/-- Cyclic reindexing agrees with the already established Fourier definition. -/
theorem literalCycle_probability (τ : ℕ → ℤ) (w : ZMod r → W) (v : ℤ → W)
    (y : ZMod N) :
    cyclicProbability (literalCycle τ w v) y =
      fourierProbability (fun x : Fin N => symbol w v (τ x)) y.val := by
  have hval (x : ZMod N) : ((ZMod.finEquiv N).symm x).val = x.val := by
    cases N with
    | zero => exact (NeZero.ne 0 rfl).elim
    | succ n => rfl
  have hf : literalCycle (N := N) τ w v =
      cyclePrefix (fun x : Fin N => symbol w v (τ x)) := by
    funext x
    simp only [literalCycle, cyclePrefix, hval]
  rw [hf, cyclicProbability_eq_fourier]

/-- The Section 4.2 constant bound, without a divisibility assumption: the report's theorem
"Literal glitches without a divisibility assumption", `Pr[Y ∈ G_k] ≥ 133/358`.

If `N < (2k+1)r` the windows cover everything; otherwise `t ≥ 2` and the Fejér argument applies.
The history supplies the matching estimate. `NearestPeaks` merely specifies
which rounded peaks the windows are centered on; every tie choice works.
There are no kernel, Fourier-identity, or concentration hypotheses. -/
theorem general_k_nondividing {A D τ} (hedit : Edited A D τ)
    (w : ZMod r → W) (v : ℤ → W) (hr : 2 ≤ r) (hk : 1 ≤ A+D)
    (g : ZMod r → ℤ) (hg : NearestPeaks (N := N) g) :
    (133 : ℝ)/358 ≤ ∑ y ∈ roundedWindows (N := N) g (A+D),
      fourierProbability (fun x : Fin N => symbol w v (τ x)) y.val := by
  simp_rw [← literalCycle_probability τ w v]
  by_cases hsize : N < (2*(A+D)+1)*r
  · rw [roundedWindows_full_of_lt g hg (A+D) hsize, cyclicProbability_sum]
    norm_num
  · let t : ℝ := 2*(N : ℝ)/(r*(2*(A+D : ℕ)+1))
    have hN : (0 : ℝ) < N := by exact_mod_cast NeZero.pos N
    have hr0 : (0 : ℝ) < r := by exact_mod_cast NeZero.pos r
    have hk0 : (0 : ℝ) < 2*(A+D : ℕ)+1 := by positivity
    have ht : 2 ≤ t := by
      have hsize' : (2*(A+D : ℕ)+1)*(r : ℝ) ≤ N := by exact_mod_cast (Nat.le_of_not_gt hsize)
      apply (le_div_iff₀ (mul_pos hr0 hk0)).mpr
      nlinarith
    apply fejer_mass_of_matching (literalCycle τ w v) r (roundedWindows g (A+D)) t ht
    · exact fun ell hell => hedit.normalised_matching w v hr hk ell hell
    · intro y hy
      have he : 1/t = ((r : ℝ)*(2*(A+D : ℕ)+1))/(2*N) := by
        dsimp [t]
        field_simp
      rw [he]
      exact roundedWindows_separation g hg (A+D) y hy

/-- The original 1/(2k+1) total-neighbourhood guarantee follows immediately. -/
theorem general_k_nondividing_universal {A D τ} (hedit : Edited A D τ)
    (w : ZMod r → W) (v : ℤ → W) (hr : 2 ≤ r) (hk : 1 ≤ A+D)
    (g : ZMod r → ℤ) (hg : NearestPeaks (N := N) g) :
    (1 : ℝ)/(2*(A+D : ℕ)+1) ≤ ∑ y ∈ roundedWindows (N := N) g (A+D),
      fourierProbability (fun x : Fin N => symbol w v (τ x)) y.val := by
  have hk1 : (1 : ℝ) ≤ (A+D : ℕ) := by exact_mod_cast hk
  have hconst : (1 : ℝ)/(2*(A+D : ℕ)+1) ≤ 133/358 := by
    apply (div_le_iff₀ (by positivity : (0 : ℝ) < 2*(A+D : ℕ)+1)).mpr
    nlinarith
  exact hconst.trans (general_k_nondividing hedit w v hr hk g hg)

/-- A version with no rounding hypothesis: use the concrete canonical peaks. -/
theorem general_k_nondividing_canonical {A D τ} (hedit : Edited A D τ)
    (w : ZMod r → W) (v : ℤ → W) (hr : 2 ≤ r) (hk : 1 ≤ A+D) :
    (133 : ℝ)/358 ≤ ∑ y ∈ roundedWindows (N := N) (canonicalPeaks N r) (A+D),
      fourierProbability (fun x : Fin N => symbol w v (τ x)) y.val :=
  general_k_nondividing hedit w v hr hk (canonicalPeaks N r) canonicalPeaks_nearest

omit [NeZero N] in
/-- When `r ∣ N` the exact ideal outcomes `c·(N/r)` are valid nearest peaks:
they equal the real peaks `cN/r`, so the rounding error is zero. -/
theorem exactPeaks_nearest (hrN : r ∣ N) :
    NearestPeaks (N := N) (fun c : ZMod r => ((c.val * (N / r) : ℕ) : ℤ)) := by
  intro c
  have hr : (r : ℝ) ≠ 0 := by exact_mod_cast NeZero.ne r
  -- `Nat.cast_div` moves the exact division `N / r` from ℕ into ℝ; it needs `r ∣ N`.
  have he : ((c.val * (N / r) : ℕ) : ℝ) = (c.val : ℝ) * N / r := by
    rw [Nat.cast_mul, Nat.cast_div hrN hr, mul_div_assoc]
  simp only [Int.cast_natCast, he, sub_self, abs_zero]
  norm_num

/-- The report's Eq. (nondividing) in the divisible case: when `r ∣ N`, the radius-k
windows `G_k` around the exact ideal outcomes carry mass at least `133/358`.
Eq. (total) gives only `1/(2k+1)` on `G_k`. On the exact set `G` no k-independent bound
holds in Model I (`no_uniform_exact_bound`). -/
theorem general_k_dividing_windows {A D τ} (hedit : Edited A D τ)
    (w : ZMod r → W) (v : ℤ → W) (hr : 2 ≤ r) (hk : 1 ≤ A+D) (hrN : r ∣ N) :
    (133 : ℝ)/358 ≤ ∑ y ∈ roundedWindows (N := N)
        (fun c : ZMod r => ((c.val * (N / r) : ℕ) : ℤ)) (A+D),
      fourierProbability (fun x : Fin N => symbol w v (τ x)) y.val :=
  general_k_nondividing hedit w v hr hk _ (exactPeaks_nearest hrN)

end Glitchy
