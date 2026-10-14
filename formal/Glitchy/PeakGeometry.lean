import Glitchy.FejerKernel

/-! Rounded, circular radius-k windows and their separation from integer phases. -/
namespace Glitchy
open Finset

variable {N r : ℕ} [NeZero N] [NeZero r]

/-- Any choice at half-integer ties is permitted. The integer representatives
need not be reduced modulo N: the neighbourhood definition does that. -/
def NearestPeaks (g : ZMod r → ℤ) : Prop :=
  ∀ c, |(g c : ℝ) - (c.val : ℝ)*N/r| ≤ (1 : ℝ)/2

/-- A concrete nearest-peak choice, using mathlib's rounding convention. -/
noncomputable def canonicalPeaks (N r : ℕ) [NeZero r] : ZMod r → ℤ :=
  fun c => round ((c.val : ℝ)*N/r)

omit [NeZero N] in
theorem canonicalPeaks_nearest : NearestPeaks (N := N) (canonicalPeaks N r) := by
  intro c
  simpa only [canonicalPeaks, abs_sub_comm] using abs_sub_round ((c.val : ℝ)*N/r)

/-- Circular windows, represented as a set, so overlaps are counted once. -/
noncomputable def roundedWindows (g : ZMod r → ℤ) (k : ℕ) : Finset (ZMod N) := by
  classical
  exact univ.filter fun y => ∃ c : ZMod r, ∃ q : ℤ,
    |(y.val : ℤ) - g c - q*(N : ℤ)| ≤ (k : ℤ)

/-- Some integer realizes the folded phase distance. -/
theorem phaseDistance_witness (r : ℕ) (y : ZMod N) :
    ∃ m : ℕ, |(r : ℝ)*y.val/N - m| = phaseDistance ((r : ZMod N)*y) := by
  let b := r*y.val
  let j := b/N
  let a := b%N
  have hN : (0 : ℝ) < N := by exact_mod_cast NeZero.pos N
  have ha : (a : ℝ) < N := by exact_mod_cast Nat.mod_lt b (NeZero.pos N)
  have hb : (b : ℝ) = a + (j : ℝ)*N := by
    dsimp [a, j]
    exact_mod_cast (by simpa [Nat.mul_comm] using (Nat.mod_add_div b N).symm : b = b%N+(b/N)*N)
  have hv : ((r : ZMod N)*y).val = a := by
    rw [← ZMod.natCast_zmod_val y, ← Nat.cast_mul, ZMod.val_natCast]
  have hx : (r : ℝ)*y.val/N = (a : ℝ)/N + j := by
    have hb' : (r : ℝ)*y.val = a + (j : ℝ)*N := by simpa [b] using hb
    apply (div_eq_iff (ne_of_gt hN)).mpr
    field_simp
    simpa using hb'
  unfold phaseDistance
  rw [hv]
  by_cases h : (a : ℝ)/N ≤ 1-(a : ℝ)/N
  · refine ⟨j, ?_⟩
    rw [hx, min_eq_left h]
    simp [abs_of_nonneg (show (0 : ℝ) ≤ a/N by positivity)]
  · refine ⟨j+1, ?_⟩
    rw [hx, min_eq_right (le_of_not_ge h)]
    push_cast
    have ha1 : (a : ℝ)/N ≤ 1 := (div_le_one hN).mpr ha.le
    rw [abs_of_nonpos (by linarith : (a : ℝ)/N + j - ((j : ℝ)+1) ≤ 0)]
    ring

/-- A real peak and its translated rounded representative near any outcome. -/
theorem rounded_peak_witness (g : ZMod r → ℤ) (hg : NearestPeaks (N := N) g)
    (y : ZMod N) :
    ∃ c : ZMod r, ∃ q : ℤ,
      |(y.val : ℝ) - (g c : ℝ) - q*N| ≤
        ((N : ℝ)/r)*phaseDistance ((r : ZMod N)*y) + 1/2 := by
  obtain ⟨m, hm⟩ := phaseDistance_witness r y
  let c : ZMod r := m
  let q : ℤ := ((m/r : ℕ) : ℤ)
  have hN : (N : ℝ) ≠ 0 := by exact_mod_cast NeZero.ne N
  have hr : (0 : ℝ) < r := by exact_mod_cast NeZero.pos r
  have hdecomp : (m : ℝ) = c.val + (q : ℝ)*r := by
    have hm' := Nat.mod_add_div m r
    dsimp [c, q]
    rw [ZMod.val_natCast]
    have hn : m = m%r+(m/r)*r := by simpa [Nat.mul_comm] using hm'.symm
    simpa only [Nat.cast_add, Nat.cast_mul, Int.cast_natCast] using
      congrArg (fun n : ℕ => (n : ℝ)) hn
  have hd : |(y.val : ℝ) - (m : ℝ)*N/r| =
      ((N : ℝ)/r)*phaseDistance ((r : ZMod N)*y) := by
    have he : (y.val : ℝ) - (m : ℝ)*N/r =
        ((N : ℝ)/r)*((r : ℝ)*y.val/N-m) := by field_simp; ring
    rw [he, abs_mul, abs_of_nonneg (by positivity : (0 : ℝ) ≤ N/r), hm]
  refine ⟨c, q, ?_⟩
  have he : (y.val : ℝ) - (g c : ℝ) - (q : ℝ)*N =
      ((y.val : ℝ)-(m : ℝ)*N/r) + ((c.val : ℝ)*N/r - (g c : ℝ)) := by
    rw [hdecomp]
    field_simp
    ring
  rw [he]
  calc
    _ ≤ |(y.val : ℝ)-(m : ℝ)*N/r| + |(c.val : ℝ)*N/r - (g c : ℝ)| := abs_add _ _
    _ ≤ ((N : ℝ)/r)*phaseDistance ((r : ZMod N)*y) + 1/2 := by
      rw [hd, abs_sub_comm]
      exact add_le_add_left (hg c) _

/-- The complement of radius-k rounded windows is separated by k+1/2
from the real peaks. This includes every possible rounding-tie choice. -/
theorem roundedWindows_separation (g : ZMod r → ℤ) (hg : NearestPeaks (N := N) g)
    (k : ℕ) (y : ZMod N) (hy : y ∉ roundedWindows (N := N) g k) :
    ((r : ℝ)*(2*k+1))/(2*N) ≤ phaseDistance ((r : ZMod N)*y) := by
  classical
  obtain ⟨c, q, hb⟩ := rounded_peak_witness g hg y
  have hout : (k : ℤ) < |(y.val : ℤ)-g c-q*(N : ℤ)| := by
    by_contra hh
    apply hy
    simp only [roundedWindows, mem_filter, mem_univ, true_and]
    exact ⟨c, q, by omega⟩
  have hint : (k : ℤ)+1 ≤ |(y.val : ℤ)-g c-q*(N : ℤ)| := by omega
  have hreal : (k : ℝ)+1 ≤ |(y.val : ℝ)-(g c : ℝ)-(q : ℝ)*N| := by exact_mod_cast hint
  have hN : (0 : ℝ) < N := by exact_mod_cast NeZero.pos N
  have hr : (0 : ℝ) < r := by exact_mod_cast NeZero.pos r
  apply (div_le_iff₀ (by positivity : (0 : ℝ) < 2*N)).mpr
  have hdist : (k : ℝ)+1/2 ≤ ((N : ℝ)/r)*phaseDistance ((r : ZMod N)*y) := by linarith
  have hm := mul_le_mul_of_nonneg_right hdist hr.le
  have he : ((N : ℝ)/r)*phaseDistance ((r : ZMod N)*y)*r =
      N*phaseDistance ((r : ZMod N)*y) := by field_simp
  rw [he] at hm
  nlinarith

/-- Strictly overlapping windows cover every outcome. -/
theorem roundedWindows_full_of_lt (g : ZMod r → ℤ) (hg : NearestPeaks (N := N) g)
    (k : ℕ) (hsize : N < (2*k+1)*r) : roundedWindows (N := N) g k = univ := by
  classical
  apply eq_univ_of_forall
  intro y
  by_contra hy
  have hs := roundedWindows_separation g hg k y hy
  have hd := (phaseDistance_mem ((r : ZMod N)*y)).2
  have hN : (0 : ℝ) < N := by exact_mod_cast NeZero.pos N
  have hsize' : (N : ℝ) < (2*k+1)*r := by exact_mod_cast hsize
  have := (div_le_iff₀ (by positivity : (0 : ℝ) < 2*N)).mp (hs.trans hd)
  nlinarith

end Glitchy
