import Glitchy.PairCounting
import Glitchy.PeakGeometry

/-!
Section 4.2 end to end: arbitrary literal edit histories, arbitrary inserted
symbols, every nearest-integer tie choice, circular windows counted as a set,
and the existing complex-exponential Fourier measurement formula.
-/
namespace Glitchy
open Finset
variable {W : Type*} [Fintype W] [DecidableEq W]
    {N r : ℕ} [NeZero N] [NeZero r]

omit [Fintype W] in
/-- The normalized matching bound follows from the literal edit history. -/
theorem Edited.normalised_matching {A D s} (hedit : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (hr : 2 ≤ r) (hk : 1 ≤ A+D)
    (ell : ℕ) (hell : 1 ≤ ell) :
    1 - (4*(ell : ℝ)+1)/(3*(2*(N : ℝ)/(r*(2*(A+D : ℕ)+1)))) ≤
      matchingFraction (literalCycle (N := N) s w v) ((ell*r : ℕ) : ZMod N) := by
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
theorem literalCycle_probability (s : ℕ → ℤ) (w : ZMod r → W) (v : ℤ → W)
    (y : ZMod N) :
    cyclicProbability (literalCycle s w v) y =
      fourierProbability (fun x : Fin N => symbol w v (s x)) y.val := by
  have hval (x : ZMod N) : ((ZMod.finEquiv N).symm x).val = x.val := by
    cases N with
    | zero => exact (NeZero.ne 0 rfl).elim
    | succ n => rfl
  have hf : literalCycle (N := N) s w v =
      cyclePrefix (fun x : Fin N => symbol w v (s x)) := by
    funext x
    simp only [literalCycle, cyclePrefix, hval]
  rw [hf, cyclicProbability_eq_fourier]

/-- The Section 4.2 constant bound, without a divisibility assumption.

The history supplies the matching estimate. `NearestPeaks` merely specifies
which rounded peaks the windows are centered on; every tie choice works.
There are no kernel, Fourier-identity, or concentration hypotheses. -/
theorem general_k_nondividing {A D s} (hedit : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (hr : 2 ≤ r) (hk : 1 ≤ A+D)
    (g : ZMod r → ℤ) (hg : NearestPeaks (N := N) g) :
    (133 : ℝ)/358 ≤ ∑ y ∈ roundedWindows (N := N) g (A+D),
      fourierProbability (fun x : Fin N => symbol w v (s x)) y.val := by
  simp_rw [← literalCycle_probability s w v]
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
    apply fejer_mass_of_matching (literalCycle s w v) r (roundedWindows g (A+D)) t ht
    · exact fun ell hell => hedit.normalised_matching w v hr hk ell hell
    · intro y hy
      have he : 1/t = ((r : ℝ)*(2*(A+D : ℕ)+1))/(2*N) := by
        dsimp [t]
        field_simp
      rw [he]
      exact roundedWindows_separation g hg (A+D) y hy

/-- The original 1/(2k+1) total-neighbourhood guarantee follows immediately. -/
theorem general_k_nondividing_universal {A D s} (hedit : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (hr : 2 ≤ r) (hk : 1 ≤ A+D)
    (g : ZMod r → ℤ) (hg : NearestPeaks (N := N) g) :
    (1 : ℝ)/(2*(A+D : ℕ)+1) ≤ ∑ y ∈ roundedWindows (N := N) g (A+D),
      fourierProbability (fun x : Fin N => symbol w v (s x)) y.val := by
  have hk1 : (1 : ℝ) ≤ (A+D : ℕ) := by exact_mod_cast hk
  have hconst : (1 : ℝ)/(2*(A+D : ℕ)+1) ≤ 133/358 := by
    apply (div_le_iff₀ (by positivity : (0 : ℝ) < 2*(A+D : ℕ)+1)).mpr
    nlinarith
  exact hconst.trans (general_k_nondividing hedit w v hr hk g hg)

/-- A version with no rounding hypothesis: use the concrete canonical peaks. -/
theorem general_k_nondividing_canonical {A D s} (hedit : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (hr : 2 ≤ r) (hk : 1 ≤ A+D) :
    (133 : ℝ)/358 ≤ ∑ y ∈ roundedWindows (N := N) (canonicalPeaks N r) (A+D),
      fourierProbability (fun x : Fin N => symbol w v (s x)) y.val :=
  general_k_nondividing hedit w v hr hk (canonicalPeaks N r) canonicalPeaks_nearest

end Glitchy
