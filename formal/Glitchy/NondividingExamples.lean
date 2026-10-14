import Glitchy.Nondividing

/-! Kernel-checked regression examples for the nondivisible theorem.
These use actual edit constructors, not assumed matching or concentration bounds.
-/
namespace Glitchy
open Finset

private def lowerTiePeaks : ZMod 2 → ℤ := fun c => if c = 0 then 0 else 3
private def upperTiePeaks : ZMod 2 → ℤ := fun c => if c = 0 then 0 else 4

private theorem lowerTie_nearest : NearestPeaks (N := 7) lowerTiePeaks := by
  intro c
  fin_cases c <;> norm_num [lowerTiePeaks, ZMod.val, abs_div]

private theorem upperTie_nearest : NearestPeaks (N := 7) upperTiePeaks := by
  intro c
  fin_cases c <;> norm_num [upperTiePeaks, ZMod.val, abs_div]

/-- Both nearest-integer choices at the peak 7/2 are accepted. The insertion
matches the old symbol 0; no distinctness or novelty assumption is used. -/
theorem nondividing_tie_regression :
    let s := insertAt (fun x => (x : ℤ)) 3 0
    let f := fun x : Fin 7 => symbol (fun c : ZMod 2 => c) (fun _ : ℤ => (0 : ZMod 2)) (s x)
    ((133 : ℝ)/358 ≤ ∑ y ∈ roundedWindows (N := 7) lowerTiePeaks 1, fourierProbability f y.val) ∧
    ((133 : ℝ)/358 ≤ ∑ y ∈ roundedWindows (N := 7) upperTiePeaks 1, fourierProbability f y.val) := by
  dsimp
  constructor
  · exact general_k_nondividing (Edited.original.insert 3) (fun c => c) (fun _ => 0)
      (by norm_num) (by norm_num) lowerTiePeaks lowerTie_nearest
  · exact general_k_nondividing (Edited.original.insert 3) (fun c => c) (fun _ => 0)
      (by norm_num) (by norm_num) upperTiePeaks upperTie_nearest

/-- Three mixed edits with one inserted tag deleted again; N=31, r=2, k=3. -/
theorem nondividing_mixed_regression :
    let s := insertAt (deleteAt (insertAt (fun x => (x : ℤ)) 4 0) 4) 9 1
    (133 : ℝ)/358 ≤ ∑ y ∈ roundedWindows (N := 31) (canonicalPeaks 31 2) 3,
      fourierProbability (fun x : Fin 31 =>
        symbol (fun c : ZMod 2 => c) (fun _ : ℤ => (1 : ZMod 2)) (s x)) y.val := by
  exact general_k_nondividing_canonical
    (((Edited.original.insert 4).delete 4).insert 9)
    (fun c => c) (fun _ => 1) (by norm_num) (by norm_num)

/-- A deletion-only example in the long-kernel regime: t=31/3 > 6. -/
theorem nondividing_long_regression :
    (133 : ℝ)/358 ≤ ∑ y ∈ roundedWindows (N := 31) (canonicalPeaks 31 2) 1,
      fourierProbability (fun x : Fin 31 =>
        symbol (fun c : ZMod 2 => c) (fun _ : ℤ => (0 : ZMod 2))
          (deleteAt (fun n => (n : ℤ)) 12 x)) y.val := by
  exact general_k_nondividing_canonical (Edited.original.delete 12)
    (fun c => c) (fun _ => 0) (by norm_num) (by norm_num)

-- Overlapping windows cover the whole domain, also when k > N.
example : roundedWindows (N := 7) (canonicalPeaks 7 3) 1 = univ :=
  roundedWindows_full_of_lt _ canonicalPeaks_nearest 1 (by norm_num)
example : roundedWindows (N := 7) (canonicalPeaks 7 2) 80 = univ :=
  roundedWindows_full_of_lt _ canonicalPeaks_nearest 80 (by norm_num)

-- The boundary t=2 is handled by the short kernel, not by a strict inequality.
example : (4 : ℝ)/3 ≤ 2*Real.sin (Real.pi/2)^2 :=
  short_sine_bound 2 (by norm_num) (by norm_num)

-- The two regimes agree at t=6; both estimates are checked independently.
example : (4 : ℝ)/3 ≤ 6*Real.sin (Real.pi/6)^2 :=
  short_sine_bound 6 (by norm_num) (by norm_num)
example : (29 : ℝ)/54 ≤ fejerLower 6 5 :=
  long_fejer_lower 6 5 (by norm_num) (by norm_num) (by norm_num)
example : 1/((5 : ℝ)^2*Real.sin (Real.pi/6)^2) ≤ 64/243 :=
  long_fejer_outside 6 5 (by norm_num) (by norm_num)

end Glitchy
