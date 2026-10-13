import Glitchy.Boxes
import Glitchy.Fourier

namespace Glitchy
open Finset Complex
variable {W : Type*} [Fintype W] [DecidableEq W] {r : ℕ} [NeZero r]

/-- Actual value-residue counts of the first N entries after the edits. -/
noncomputable def prefixCounts (s : ℕ → ℤ) (w : ZMod r → W) (v : ℤ → W)
    (N : ℕ) (a : W) (b : ZMod r) : ℝ :=
  fiberCounts (box s w v : Fin N → W × ZMod r) univ (a, b)

/-- The box bound connected to Fourier probabilities, not assumed as a hypothesis. -/
theorem literal_total {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) :
    (r : ℝ) / (((A + D : ℕ) + 1 : ℝ) * r + (inserteds s N).card) ≤
      ∑ c : ZMod r, idealProbability N (prefixCounts s w v N) c := by
  rw [residue_mass]
  let n := fiberCounts (box s w v : Fin N → W × ZMod r) univ
  have hsum : ∑ b, n b = (N : ℝ) := by
    simpa [n] using sum_fiberCounts (box s w v : Fin N → W × ZMod r) univ
  have hq : (univ.filter (fun b => n b ≠ 0)).card ≤
      (A + D + 1) * r + (inserteds s N).card := by
    rw [support_fiberCounts]
    exact h.all_boxes w v N
  have hb := mass_from_boxes n N r ((A + D + 1) * r + (inserteds s N).card)
    hN (by have := NeZero.pos r; positivity) hsum hq
  simpa [n, prefixCounts, Fintype.sum_prod_type, Nat.cast_add, Nat.cast_mul] using hb

theorem literal_retained {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) :
    ((N : ℝ) - (inserteds s N).card)^2 /
      (((A + D : ℕ) + 1 : ℝ) * (N : ℝ)^2) ≤
      ∑ c : ZMod r, idealProbability N (prefixCounts s w v N) c := by
  rw [residue_mass]
  let n := fiberCounts (box s w v : Fin N → W × ZMod r) univ
  let m := fiberCounts (box s w v : Fin N → W × ZMod r) (originals s N)
  have hsum : ∑ b, m b = (N : ℝ) - (inserteds s N).card := by
    rw [sum_fiberCounts]
    have hc := original_card s N
    have hc' : ((originals s N).card : ℝ) + (inserteds s N).card = N := by
      exact_mod_cast hc
    linarith
  have hq : (univ.filter (fun b => m b ≠ 0)).card ≤ (A + D + 1) * r := by
    rw [support_fiberCounts]
    exact h.original_boxes w v N
  have hb := mass_from_retained n m N r (A + D) (inserteds s N).card hN (NeZero.pos r)
    (fun b => fiberCounts_mono _ (subset_univ _) b) hsum hq
  simpa [n, prefixCounts, Fintype.sum_prod_type] using hb

theorem literal_universal {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) :
    (1 : ℝ) / (2 * (A + D : ℕ) + 1) ≤
      ∑ c : ZMod r, idealProbability N (prefixCounts s w v N) c :=
  (universal_bound r (A + D) (inserteds s N).card (NeZero.pos r)
    ((h.inserted_card_le N).trans (Nat.le_add_right A D))).trans
    (literal_total h w v N hN)

/-- Grouping a finite sum by the residue boxes, before taking absolute squares. -/
theorem weighted_fibers {X B : Type*} [Fintype B] [DecidableEq B]
    (f : X → B) (T : Finset X) (F : B → ℂ) :
    (∑ b, (fiberCounts f T b : ℂ) * F b) = ∑ x ∈ T, F (f x) := by
  have h (b : B) : (fiberCounts f T b : ℂ) * F b =
      ∑ x ∈ T.filter (fun x => f x = b), F (f x) := by
    unfold fiberCounts
    push_cast
    calc
      _ = ∑ _x ∈ T.filter (fun x => f x = b), F b := by simp
      _ = _ := by
        apply sum_congr rfl
        intro x hx
        rw [(mem_filter.mp hx).2]
  simp_rw [h]
  exact sum_fiberwise T f (fun x => F (f x))

/-- The actual (unnormalized) Fourier amplitudes at ideal outcomes agree with
the amplitudes computed from the count table. -/
theorem amplitude_by_residue (s : ℕ → ℤ) (w : ZMod r → W) (v : ℤ → W)
    (N : ℕ) (a : W) (c : ZMod r) :
    (∑ x : Fin N, if symbol w v (s x) = a
      then ZMod.stdAddChar (c * (x.val : ZMod r)) else 0) =
    ∑ b : ZMod r, (prefixCounts s w v N a b : ℂ) * ZMod.stdAddChar (c * b) := by
  have h := weighted_fibers (box s w v : Fin N → W × ZMod r) univ
    (fun p => if p.1 = a then ZMod.stdAddChar (c * p.2) else 0)
  rw [Fintype.sum_prod_type] at h
  simp only [mul_ite, mul_zero, ← sum_ite_irrel, sum_ite_eq', mem_univ, if_true] at h
  simpa [box, prefixCounts] using h.symm

/-- Fourier measurement for a literal prefix, using the report's positive exponent.
y is a natural number; only its residue modulo N affects the result. -/
noncomputable def fourierProbability {N : ℕ} (f : Fin N → W) (y : ℕ) : ℝ :=
  (1 / (N : ℝ)^2) * ∑ a : W,
    normSq (∑ x : Fin N, if f x = a
      then Complex.exp (2 * Real.pi * I * (x.val : ℂ) * (y : ℂ) / (N : ℂ)) else 0)

theorem ideal_phase (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) (x : Fin N) (c : ZMod r) :
    ZMod.stdAddChar (c * (x.val : ZMod r)) =
      Complex.exp (2 * Real.pi * I * (x.val : ℂ) * (c.val * (N / r) : ℕ) / (N : ℂ)) := by
  have hc : c * (x.val : ZMod r) = ((c.val * x.val : ℕ) : ZMod r) := by simp
  rw [hc, ← Int.cast_natCast, ZMod.stdAddChar_coe]
  push_cast
  congr 1
  have hn : (N : ℂ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hN)
  have hr : (r : ℂ) ≠ 0 := by exact_mod_cast (NeZero.ne r)
  have hd : ((N / r : ℕ) : ℂ) * r = N := by exact_mod_cast (Nat.div_mul_cancel hrN)
  field_simp
  rw [← hd]
  ring

/-- Bridge from the literal sequence's QFT formula to its residue-count formula. -/
theorem idealProbability_eq_fourier (s : ℕ → ℤ) (w : ZMod r → W) (v : ℤ → W)
    (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) (c : ZMod r) :
    idealProbability N (prefixCounts s w v N) c =
      fourierProbability (fun x : Fin N => symbol w v (s x)) (c.val * (N / r)) := by
  unfold idealProbability fourierProbability
  congr 1
  apply sum_congr rfl
  intro a ha
  rw [← amplitude_by_residue]
  congr 1
  apply sum_congr rfl
  intro x hx
  rw [ideal_phase N hN hrN x c]

/-- The complete general-k theorem for literal edits and the QFT formula. -/
theorem general_k_glitches {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) :
    (1 : ℝ) / (2 * (A + D : ℕ) + 1) ≤
      ∑ c : ZMod r,
        fourierProbability (fun x : Fin N => symbol w v (s x)) (c.val * (N / r)) := by
  simpa only [idealProbability_eq_fourier s w v N hN hrN] using
    literal_universal h w v N hN

/-- Actual ideal outcomes, each counted once. -/
noncomputable def goodOutcomes (N r : ℕ) [NeZero r] : Finset ℕ :=
  univ.image (fun c : ZMod r => c.val * (N / r))

theorem good_index_injective (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) :
    Function.Injective (fun c : ZMod r => c.val * (N / r)) := by
  have hm : 0 < N / r := Nat.div_pos (Nat.le_of_dvd hN hrN) (NeZero.pos r)
  intro a b hab
  apply ZMod.val_injective
  exact mul_right_cancel₀ (Nat.ne_of_gt hm) hab

theorem good_index_lt (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) (c : ZMod r) :
    c.val * (N / r) < N := by
  have hm : 0 < N / r := Nat.div_pos (Nat.le_of_dvd hN hrN) (NeZero.pos r)
  calc
    _ < r * (N / r) := Nat.mul_lt_mul_of_pos_right (ZMod.val_lt c) hm
    _ = N := by rw [Nat.mul_comm, Nat.div_mul_cancel hrN]

theorem sum_goodOutcomes (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) (P : ℕ → ℝ) :
    (∑ y ∈ goodOutcomes N r, P y) = ∑ c : ZMod r, P (c.val * (N / r)) := by
  unfold goodOutcomes
  exact sum_image (fun a _ b _ hab => good_index_injective N hN hrN hab)

/-- Both sharper bounds, for the actual set G rather than an indexed multiset. -/
theorem general_k_sharp_bounds {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) (hrN : r ∣ N) :
    let P := fourierProbability (fun x : Fin N => symbol w v (s x))
    let I := (inserteds s N).card
    (r : ℝ) / (((A + D : ℕ) + 1 : ℝ) * r + I) ≤ ∑ y ∈ goodOutcomes N r, P y
    ∧ ((N : ℝ) - I)^2 / (((A + D : ℕ) + 1 : ℝ) * (N : ℝ)^2)
      ≤ ∑ y ∈ goodOutcomes N r, P y := by
  dsimp
  rw [sum_goodOutcomes N hN hrN]
  constructor
  · simpa only [idealProbability_eq_fourier s w v N hN hrN] using
      literal_total h w v N hN
  · simpa only [idealProbability_eq_fourier s w v N hN hrN] using
      literal_retained h w v N hN

theorem fourierProbability_nonneg {N : ℕ} (f : Fin N → W) (y : ℕ) :
    0 ≤ fourierProbability f y := by
  unfold fourierProbability
  exact mul_nonneg (by positivity) (sum_nonneg fun _ _ => normSq_nonneg _)

/-- Any neighbourhood containing G inherits the bound, even when its windows overlap. -/
theorem general_k_neighbourhood {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) (hrN : r ∣ N)
    (H : Finset ℕ) (hH : goodOutcomes N r ⊆ H) :
    (1 : ℝ) / (2 * (A + D : ℕ) + 1) ≤
      ∑ y ∈ H, fourierProbability (fun x : Fin N => symbol w v (s x)) y := by
  have hb := general_k_glitches h w v N hN hrN
  rw [← sum_goodOutcomes N hN hrN] at hb
  exact hb.trans (mass_mono _ (fourierProbability_nonneg _) hH)

/-- Deletions only, or any history whose inserted entries all disappear from the prefix. -/
theorem general_k_no_surviving_insertions {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) (hrN : r ∣ N)
    (hI : (inserteds s N).card = 0) :
    (1 : ℝ) / ((A + D : ℕ) + 1) ≤
      ∑ y ∈ goodOutcomes N r, fourierProbability (fun x : Fin N => symbol w v (s x)) y := by
  have hb := (general_k_sharp_bounds h w v N hN hrN).1
  simpa only [hI, Nat.cast_zero, no_insertions r (A + D) (NeZero.pos r)] using hb

/-- The finite-N version of the "k small compared with N" corollary. -/
theorem general_k_few_edits {A D s} (h : Edited A D s)
    (w : ZMod r → W) (v : ℤ → W) (N : ℕ) (hN : 0 < N) (hrN : r ∣ N)
    (hkN : A + D ≤ N) :
    (1 - ((A + D : ℕ) : ℝ) / N)^2 / ((A + D : ℕ) + 1) ≤
      ∑ y ∈ goodOutcomes N r, fourierProbability (fun x : Fin N => symbol w v (s x)) y := by
  have hI := (h.inserted_card_le N).trans (Nat.le_add_right A D)
  exact (few_edits N (A + D) (inserteds s N).card hN hI hkN).trans
    (general_k_sharp_bounds h w v N hN hrN).2

end Glitchy
